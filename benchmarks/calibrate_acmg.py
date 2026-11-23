"""Calibrate the model score to ACMG/AMP PP3/BP4 evidence strengths with the ClinGen method.

    pip install -e ".[benchmark]"
    python benchmarks/calibrate_acmg.py [--bootstrap 10000]

Method: a port of the local-posterior calibration of Pejaver et al. (AJHG 2022, PMC9748256), which set
the ClinGen PP3/BP4 thresholds of REVEL and other tools, following their MATLAB code
(github.com/vpejaver/clingen-svi-comp_calibration, local_posterior_probability/):
  * prior probability of pathogenicity 0.0441; evidence strength in Bayes points, a level of n points
    needing a likelihood ratio of 1124^(n/8): Supporting (1 point) 2.406, Moderate (2) 5.790,
    3 points 13.93, Strong (4) 33.53;
  * local posterior probability at each score, from a window around it holding at least 100 ClinVar
    variants and 3% of the rare gnomAD variants (fewer at the ends of the score range), with the
    benign variants weighted to the prior;
  * the threshold of a level is the lowest score above which the posterior never drops below the
    level's target (for BP4, the highest score below which the benign posterior never drops), taken at
    the one-sided 95% bound over bootstrap samples (10,000, as in the paper).
Check of the port: on the paper's own data it reproduces the published thresholds of REVEL, VEST4,
MutPred2 and PolyPhen-2 to within 0.002 (bootstrap noise); REVEL is re-checked on every run.

Data
  * Calibration: the paper's ClinVar 2019 set (Data S1; 11,834 missense variants below 1% in gnomAD,
    training variants of REVEL and other tools removed), the very data REVEL's ClinGen thresholds come
    from. The variants are mapped to GRCh38 through their ClinVar VariationID, annotated by the
    pipeline, and each is scored by a model fitted without its gene (5 folds by gene over the training
    set of train_model.py), because the shipped model was trained on ClinVar variants of these genes.
  * gnomAD set: rare gnomAD v4 missense variants not in ClinVar, sampled at random from 500 genes with
    pathogenic variants (gnomAD GraphQL API). The paper's own gnomAD set is on GRCh37; REVEL calibrated
    with this set instead gives the same thresholds, which shows the substitution is sound.
  * Tests (as the paper's ClinVar 2020 set): the paper's ClinVar 2020 set (Data S3), scored the same
    way; the rare variants of this project's held-out benchmark (genes the shipped model never saw);
    and among those, the ones labelled since 2023.
  * Sensitivity: the same calibration on this project's own training data (ClinVar 2026, >= 2 review
    stars, at most 8 variants per gene and class, rare variants, scored out of fold).

Writes src/missensexplorer/calibration_params.py and benchmarks/results/acmg_calibration.{md,json}.
"""

from __future__ import annotations

import argparse
import datetime as dt
import functools
import gzip
import io
import json
import pprint
import random
import sys
import time
import zipfile
from pathlib import Path

import httpx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import clinvar_benchmark as B  # noqa: E402
import compare_tools as C  # noqa: E402
import train_model as T  # noqa: E402
from missensexplorer import PredictOptions, acmg, model, predict  # noqa: E402
from missensexplorer.cache import default_cache_dir  # noqa: E402
from missensexplorer.model_params import MODEL  # noqa: E402
from missensexplorer.predict import _num  # noqa: E402

PARAMS = ROOT / "src" / "missensexplorer" / "calibration_params.py"
VALIDATION = default_cache_dir() / "validation"
SUPPLEMENT = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9748256/supplementaryFiles"
PEJAVER = {"2019": "mmc2.xlsx", "gnomad": "mmc3.xlsx", "2020": "mmc4.xlsx"}
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_QUERY = """query($gene: String!) { gene(gene_symbol: $gene, reference_genome: GRCh38) {
  variants(dataset: gnomad_r4) { variant_id consequence exome { ac an } genome { ac an } } } }"""
N_GENES, N_GNOMAD = 500, 6000

ALPHA, C_LR = 0.0441, 1124.0        # prior probability of pathogenicity; LR+ of "very strong"
LEVELS = (1, 2, 3, 4)               # Bayes points: Supporting, Moderate, (3), Strong
MIN_POINTS, GNOMAD_FRACTION, STEP = 100, 0.03, 0.001
PUBLISHED = {  # ClinGen thresholds as in acmg.py, by points 1-4 (None = not reached)
    "REVEL": ((0.644, 0.773, 0.879, 0.932), (0.290, 0.183, 0.052, 0.016)),
    "AlphaMissense": ((0.792, 0.906, 0.972, 0.990), (0.169, 0.099, 0.070, None)),
}
SCORES = {"model": "model", "REVEL": "revel", "AlphaMissense": "alphamissense"}


def lr_plus(points: float) -> float:
    return C_LR ** (points / 8)


def target_pathogenic(points: float) -> float:
    lr = lr_plus(points)
    return lr * ALPHA / ((lr - 1) * ALPHA + 1)


def target_benign(points: float) -> float:
    lr = lr_plus(points)
    return lr * (1 - ALPHA) / ((lr - 1) * (1 - ALPHA) + 1)


# --------------------------------------------------------------------------- data

def _download(url: str, f: Path, what: str) -> None:
    f.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {what} ...", flush=True)
    with httpx.stream("GET", url, timeout=600, follow_redirects=True) as r, open(f.with_suffix(".part"), "wb") as fh:
        r.raise_for_status()
        for chunk in r.iter_bytes(1 << 20):
            fh.write(chunk)
    f.with_suffix(".part").rename(f)


def pejaver_table(name: str) -> pd.DataFrame:
    """Data S1-S3 of Pejaver et al. 2022 ("2019", "2020": ClinVar sets; "gnomad"), from Europe PMC."""
    d = VALIDATION / "pejaver2022"
    f = d / f"{name}.csv.gz"
    if f.exists():
        return pd.read_csv(f, low_memory=False)
    if not (d / "supplement.zip").exists():
        _download(SUPPLEMENT, d / "supplement.zip", "the supplement of Pejaver et al. 2022 (80 MB)")
    import openpyxl

    with zipfile.ZipFile(d / "supplement.zip") as z:
        ws = openpyxl.load_workbook(io.BytesIO(z.read(PEJAVER[name])), read_only=True).worksheets[0]
    rows = ws.iter_rows(values_only=True)
    df = pd.DataFrame(rows, columns=next(rows))
    if name == "gnomad":
        df = df[["MAF", "REVEL_score"]]
    df.to_csv(f, index=False)
    return df


def _summary_lines():
    with gzip.open(B.download(), "rt", encoding="utf-8", errors="replace") as fh:
        ix = {k: i for i, k in enumerate(fh.readline().lstrip("#").rstrip("\n").split("\t"))}
        for line in fh:
            yield ix, line.rstrip("\n").split("\t")


def _vcf_key(ix: dict, f: list[str]) -> tuple[str, str, str, str]:
    return (f[ix["Chromosome"]], f[ix["PositionVCF"]], f[ix["ReferenceAlleleVCF"]], f[ix["AlternateAlleleVCF"]])


@functools.lru_cache(maxsize=1)
def clinvar_scan() -> tuple[set[str], set[str]]:
    """Every GRCh38 variant in ClinVar as a gnomAD variant id (chrom-pos-ref-alt), and the genes with
    at least one pathogenic or likely pathogenic variant of any type (the paper's gene filter)."""
    keys, genes = set(), set()
    for ix, f in _summary_lines():
        if f[ix["Assembly"]] != "GRCh38":
            continue
        keys.add("-".join(_vcf_key(ix, f)))
        if f[ix["ClinicalSignificance"]] in B.PATHOGENIC:
            sym = f[ix["GeneSymbol"]]          # "A;B" for variants listed under two genes
            genes.update([sym, *sym.split(";")])
    return keys, genes


def to_grch38(t: pd.DataFrame) -> list[str | None]:
    """Pipeline inputs on GRCh38 for the paper's hg19 variants, through their ClinVar VariationID."""
    keys = [(str(c), str(p), r, a) for c, p, r, a in zip(t["hg19_chr"], t["hg19_pos(1-based)"], t["ref"], t["alt"])]
    wanted, vid, g38 = set(keys), {}, {}
    for ix, f in _summary_lines():
        if f[ix["Assembly"]] == "GRCh37" and (k := _vcf_key(ix, f)) in wanted:
            vid[k] = f[ix["VariationID"]]
    ids = set(vid.values())
    for ix, f in _summary_lines():
        if f[ix["Assembly"]] == "GRCh38" and f[ix["VariationID"]] in ids:
            g38[f[ix["VariationID"]]] = _vcf_key(ix, f)
    out = []
    for k in keys:
        c, p, r, a = g38.get(vid.get(k), (None,) * 4)
        out.append(f"{c}:{p} {r}>{a}" if c and p and p.isdigit() and len(r) == len(a) == 1 else None)
    return out


def annotate(items: list[dict]) -> pd.DataFrame:
    """Predictor scores, gnomAD grpmax and the current PP3/BP4 points of variants given as dicts with
    id, gene, label and input (pipeline input text)."""
    res = predict(text="\n".join(p["input"] for p in items), options=PredictOptions(**T.OPTIONS),
                  progress=T._progress)
    print()
    T.check_sources(res)
    rows = [{"id": p["id"], "gene": r.get("gene") or p["gene"], "label": p["label"], **p.get("extra", {}),
             **{f: _num(r.get(f)) for f in T.FEATURES}, "grpmax": _num(r.get("gnomad_popmax_af")),
             "current_points": acmg.computational(r)[1]}
            for p, r in T._missense_rows(res, items)]
    return pd.DataFrame(rows)


def pejaver_clinvar(year: str) -> pd.DataFrame:
    """The paper's ClinVar 2019 or 2020 set, annotated on GRCh38 (cached)."""
    f = VALIDATION / "pejaver2022" / f"clinvar{year}_annotated.csv.gz"
    if f.exists():
        return pd.read_csv(f)
    t = pejaver_table(year)
    items = [{"id": i, "gene": g, "label": int("athogenic" in s), "input": inp}
             for i, (g, s, inp) in enumerate(zip(t["genename"], t["clnsig"], to_grch38(t))) if inp]
    print(f"ClinVar {year} set: {len(items)} of {len(t)} variants mapped to GRCh38")
    df = annotate(items)
    df.to_csv(f, index=False)
    return df


def held_out_scores(target: pd.DataFrame, folds: int = 5, seed: int = 7) -> np.ndarray:
    """Model score of each variant from a model fitted on the training set without the variant's gene."""
    train = pd.read_csv(T.OUT / "model_training_set.csv.gz")
    genes = sorted(set(train["gene"]) | set(target["gene"]))
    fold = dict(zip(genes, np.random.default_rng(seed).integers(0, folds, len(genes))))
    tf, gf = train["gene"].map(fold).to_numpy(), target["gene"].map(fold).to_numpy()
    out = np.full(len(target), np.nan)
    for k in range(folds):
        tr, sel = train[tf != k], gf == k
        m = T.pipeline().fit(tr[T.FEATURES], tr["label"],
                             logisticregression__sample_weight=T.class_weights(tr["label"].to_numpy()))
        if sel.any():
            out[sel] = m.predict_proba(target.loc[sel, T.FEATURES])[:, 1]
    out[(target["alphamissense"].isna() & target["revel"].isna()).to_numpy()] = np.nan    # as model.score
    return out


def own_set() -> pd.DataFrame:
    """Sensitivity data: this project's rare training variants in genes with a pathogenic variant,
    scored out of fold (5-fold cross-validation split by gene, as in train_model.py)."""
    df = pd.read_csv(T.OUT / "model_training_set.csv.gz", dtype={"clinvar_id": str})
    df["model"] = T.cv_predictions(df, MODEL["training"].get("seed", 7))
    df.loc[df["alphamissense"].isna() & df["revel"].isna(), "model"] = np.nan
    df = df[df["grpmax"].fillna(0) < T.RARE]
    return df[df["gene"].isin(clinvar_scan()[1])].reset_index(drop=True)


def gnomad_gene(http: httpx.Client, gene: str) -> list[list]:
    """[variant_id, allele count, allele number] of the gene's missense variants in gnomAD v4 (cached)."""
    f = VALIDATION / "gnomad_genes" / f"{gene}.json"
    if f.exists():
        return json.loads(f.read_text())
    for attempt in range(6):
        try:
            r = http.post(GNOMAD_API, json={"query": GNOMAD_QUERY, "variables": {"gene": gene}})
            if r.status_code == 429 or r.status_code >= 500:
                raise httpx.HTTPStatusError(str(r.status_code), request=r.request, response=r)
            d = r.json()
            break
        except (httpx.HTTPError, ValueError):
            time.sleep(10 * (attempt + 1))
    else:
        raise SystemExit(f"gnomAD API keeps failing for {gene}; run again later (genes done are cached)")
    out = []
    for v in ((d.get("data") or {}).get("gene") or {}).get("variants") or []:
        if v["consequence"] != "missense_variant":
            continue
        ac = sum((v[k] or {}).get("ac") or 0 for k in ("exome", "genome"))
        an = sum((v[k] or {}).get("an") or 0 for k in ("exome", "genome"))
        if ac > 0 and an > 0:
            out.append([v["variant_id"], ac, an])
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out))
    time.sleep(1.0)  # be gentle with the public API
    return out


def gnomad_set(genes: list[str]) -> pd.DataFrame:
    """Rare gnomAD missense variants not in ClinVar, sampled at random from N_GENES of `genes`."""
    f = VALIDATION / "gnomad_calibration_set.csv.gz"
    if f.exists():
        return pd.read_csv(f)
    rng = random.Random(1)
    single = sorted(g for g in genes if ";" not in g)
    chosen = rng.sample(single, min(N_GENES, len(single)))
    in_clinvar = clinvar_scan()[0]
    pool = []
    with httpx.Client(timeout=120, headers={"User-Agent": "missensexplorer-benchmark"}) as http:
        for i, g in enumerate(chosen, 1):
            print(f"  gnomAD genes: {i}/{len(chosen)}", end="\r", flush=True)
            pool += [(vid, g) for vid, ac, an in gnomad_gene(http, g)
                     if ac / an < T.RARE and vid not in in_clinvar]
    print(f"\n  {len(pool)} rare gnomAD missense variants not in ClinVar in {len(chosen)} genes")
    picked = [{"id": vid, "gene": g, "input": "{}:{} {}>{}".format(*vid.split("-"))}
              for vid, g in rng.sample(pool, min(N_GNOMAD, len(pool)))]
    res = predict(text="\n".join(p["input"] for p in picked), options=PredictOptions(**T.OPTIONS),
                  progress=T._progress)
    print()
    T.check_sources(res)
    rows = [{"variant": p["id"], "gene": p["gene"], **{k: _num(r.get(k)) for k in T.FEATURES},
             "grpmax": _num(r.get("gnomad_popmax_af"))} for p, r in T._missense_rows(res, picked)]
    df = pd.DataFrame(rows)
    df = df[df["grpmax"].fillna(0) < T.RARE].reset_index(drop=True)
    df.to_csv(f, index=False)
    return df


def model_scores(df: pd.DataFrame) -> np.ndarray:
    """Scores of the shipped model."""
    return np.array([np.nan if (ms := model.score(r)) is None else ms.p for r in df.to_dict("records")])


def benchmark_set() -> pd.DataFrame:
    """Rare variants of the held-out benchmark of clinvar_benchmark.py (genes left out of training)."""
    cands = B.load_candidates(B.download())
    old = C.labelled_in_2023()
    items = [{**p, "extra": {"new": p["id"] not in old}} for p in B.sample(cands, 1000, 8, 42)]
    df = annotate(items)
    df["model"] = model_scores(df)
    return df[df["grpmax"].fillna(0) < T.RARE].reset_index(drop=True)


# --------------------------------------------------------------------------- calibration
# Port of local_posterior_probability/*.m; scores are examined on a 0.001 grid (all three lie in 0-1).

GRID = np.round(np.arange(1.0, -STEP / 2, -STEP), 3)        # high to low, as `thrs` in main.m


def _count(a: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    return np.searchsorted(a, hi, "right") - np.searchsorted(a, lo, "left")


def local_posterior(x: np.ndarray, y: np.ndarray, g: np.ndarray, w: float, thrs: np.ndarray = GRID) -> np.ndarray:
    """Posterior probability of pathogenicity at each score in `thrs` (get_both_local_posteriors.m).

    The half-window around a score grows in steps of STEP until it holds at least MIN_POINTS ClinVar
    variants and GNOMAD_FRACTION of the gnomAD variants, both scaled by the share of the window that
    lies within the score range. The condition only gets easier as the window grows, so the smallest
    window is found by bisection instead of step by step."""
    xs, pos, neg, gs = np.sort(x), np.sort(x[y == 1]), np.sort(x[y == 0]), np.sort(g)
    top, bottom = thrs.max(), thrs.min()

    def enough(k: np.ndarray) -> np.ndarray:
        h = k * STEP
        lo, hi = thrs - h, thrs + h
        width = np.where(h > 0, hi - lo, 1.0)
        c = np.where(hi > top, (top - lo) / width, np.where(lo < bottom, (hi - bottom) / width, 1.0))
        return (_count(xs, lo, hi) >= c * MIN_POINTS) & (_count(gs, lo, hi) >= GNOMAD_FRACTION * c * len(gs))

    lo_k = np.zeros(len(thrs), int)
    hi_k = np.full(len(thrs), int(round((top - bottom) / STEP)) + 1)   # a window covering everything
    while (lo_k < hi_k).any():
        mid = (lo_k + hi_k) // 2
        ok = enough(mid)
        hi_k, lo_k = np.where(ok, mid, hi_k), np.where(ok, lo_k, mid + 1)
    h = hi_k * STEP
    p, b = _count(pos, thrs - h, thrs + h), _count(neg, thrs - h, thrs + h)
    return p / (p + w * b)


def crossing(post: np.ndarray, thrs: np.ndarray, target: float) -> np.ndarray:
    """get_all_thresholds.m: going through `thrs` in order, the last score before the posterior first
    drops below `target` (NaN when it starts below the target, or never drops, as in the MATLAB code)."""
    below = post < target
    first = np.argmax(below, axis=-1)
    return np.where(below.any(axis=-1) & (first > 0), thrs[np.maximum(first - 1, 0)], np.nan)


def discounted(t: np.ndarray, kind: str, frac: float = 0.05) -> float | None:
    """get_discounted_thresholds.m: the one-sided 95% bound of the bootstrap thresholds, on the strict
    side (high for PP3, low for BP4); None when more than 5% of the samples do not reach the level."""
    n, bad = len(t), int(np.isnan(t).sum())
    if bad > frac * n:
        return None
    v = np.sort(t[~np.isnan(t)])
    return float((v[::-1] if kind == "pathogenic" else v)[int(np.floor(frac * n)) - bad])


def calibrate(x: np.ndarray, y: np.ndarray, g: np.ndarray, n_boot: int, seed: int = 0) -> dict:
    keep = ~np.isnan(x)
    x, y, g = x[keep], y[keep], g[~np.isnan(g)]
    w = (1 - ALPHA) * (y == 1).sum() / ((y == 0).sum() * ALPHA)     # weight of one benign variant
    post = local_posterior(x, y, g, w)
    rng = np.random.default_rng(seed)
    boots = np.empty((n_boot, len(GRID)))
    for b in range(n_boot):
        i, j = rng.integers(0, len(x), len(x)), rng.integers(0, len(g), len(g))
        boots[b] = local_posterior(x[i], y[i], g[j], w)
    asc = GRID[::-1]
    out = {"n_pathogenic": int((y == 1).sum()), "n_benign": int((y == 0).sum()), "n_gnomad": int(len(g)),
           "pp3": {}, "bp4": {}, "pp3_point_estimate": {}, "bp4_point_estimate": {},
           "posterior": [[float(s), round(float(p), 5)] for s, p in zip(GRID[::10], post[::10])]}
    for n in LEVELS:
        tp, tb = target_pathogenic(n), target_benign(n)
        out["pp3_point_estimate"][n] = _f(crossing(post, GRID, tp))
        out["bp4_point_estimate"][n] = _f(crossing(1 - post[::-1], asc, tb))
        out["pp3"][n] = discounted(crossing(boots, GRID, tp), "pathogenic")
        out["bp4"][n] = discounted(crossing(1 - boots[:, ::-1], asc, tb), "benign")
    return out


def _f(v) -> float | None:
    v = float(v)
    return None if v != v else round(v, 3)


# --------------------------------------------------------------------------- evaluation

def points(s: np.ndarray, pp3: dict, bp4: dict) -> np.ndarray:
    """Bayes points of each score (0 = no evidence or no score), as acmg.computational assigns them."""
    out = np.zeros(len(s))
    for n in sorted(LEVELS):                      # weaker first, so stronger levels overwrite
        if pp3.get(n) is not None:
            out[s >= pp3[n]] = n
        if bp4.get(n) is not None:
            out[s <= bp4[n]] = -n
    out[np.isnan(s)] = 0
    return out


def evidence(pts: np.ndarray, y: np.ndarray) -> dict:
    P, N = pts[y == 1], pts[y == 0]
    r = lambda v: round(float(np.mean(v)), 4)  # noqa: E731
    return {"n_pathogenic": int(len(P)), "n_benign": int(len(N)),
            "pathogenic_pp3": r(P > 0), "pathogenic_pp3_moderate_or_more": r(P >= 2),
            "pathogenic_pp3_strong": r(P >= 4), "pathogenic_bp4": r(P < 0),
            "benign_bp4": r(N < 0), "benign_bp4_moderate_or_more": r(N <= -2),
            "benign_bp4_strong": r(N <= -4), "benign_pp3": r(N > 0),
            "mean_points_pathogenic": round(float(P.mean()), 3), "mean_points_benign": round(float(N.mean()), 3)}


def band_lrs(pts: np.ndarray, y: np.ndarray, rng: np.random.Generator, n_boot: int = 2000) -> dict:
    """Likelihood ratio in each evidence band (for BP4: benign over pathogenic), with a 95% bootstrap
    interval, against the ratio the level requires: the test of the paper's ClinVar 2020 set."""
    def lr(p, yy, k):
        inb = p == k
        a, b = np.mean(inb[yy == 1]), np.mean(inb[yy == 0])
        num, den = (a, b) if k > 0 else (b, a)
        return np.inf if den == 0 else num / den

    idx = [rng.integers(0, len(y), len(y)) for _ in range(n_boot)]
    out = {}
    for k in (4, 3, 2, 1, -1, -2, -3, -4):
        if not (pts == k).any():
            continue
        bs = np.array([lr(pts[i], y[i], k) for i in idx])
        out[k] = {"n_pathogenic": int(((pts == k) & (y == 1)).sum()), "n_benign": int(((pts == k) & (y == 0)).sum()),
                  "lr": _r(lr(pts, y, k)),           # "nearest": no interpolation between infinite ratios
                  "ci": [_r(np.percentile(bs, 2.5, method="nearest")), _r(np.percentile(bs, 97.5, method="nearest"))],
                  "required": round(lr_plus(abs(k)), 2)}
    return out


def _r(v: float) -> float | str:
    return "∞" if np.isinf(v) else round(float(v), 2)




def shares(pts: np.ndarray) -> dict:
    return {int(k): round(float(np.mean(pts == k)), 4) for k in (4, 3, 2, 1, 0, -1, -2, -3, -4)}


def _num_col(d: pd.DataFrame, col: str) -> np.ndarray:
    return pd.to_numeric(d[col], errors="coerce").to_numpy(float)


def port_check(n_boot: int, seed: int) -> dict:
    """REVEL through this port on the paper's own data (Data S1 and S2), for comparison with the
    published thresholds."""
    cv, g = pejaver_table("2019"), pejaver_table("gnomad")
    y = cv["clnsig"].str.contains("athogenic").astype(int).to_numpy()
    return calibrate(_num_col(cv, "REVEL_score"), y, _num_col(g, "REVEL_score"), n_boot, seed)


def evaluate(d: pd.DataFrame, cal: dict, rng: np.random.Generator) -> dict:
    y = d["label"].to_numpy()
    pts = points(_num_col(d, "model"), cal["pp3"], cal["bp4"])
    cur = d["current_points"].fillna(0).to_numpy(float)
    return {"model": evidence(pts, y), "current": evidence(cur, y),
            "bands": band_lrs(pts, y, rng), "current_bands": band_lrs(cur, y, rng)}


def _counts(d: pd.DataFrame) -> dict:
    return {"variants": len(d), "pathogenic": int(d["label"].sum()), "benign": int((d["label"] == 0).sum()),
            "genes": int(d["gene"].nunique())}


# --------------------------------------------------------------------------- main

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bootstrap", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)
    t0 = time.time()

    c19, c20 = pejaver_clinvar("2019"), pejaver_clinvar("2020")
    for d in (c19, c20):
        d["model"] = held_out_scores(d)
    own = own_set()
    gn = gnomad_set(sorted(set(own["gene"])))
    gn["model"] = model_scores(gn)
    bench = benchmark_set()
    print(f"calibration: {_counts(c19)}; gnomAD set: {len(gn)}; tests: 2020 {_counts(c20)}, "
          f"benchmark {_counts(bench)}", flush=True)

    def run(d: pd.DataFrame, col: str) -> dict:
        t1 = time.time()
        out = calibrate(_num_col(d, col), d["label"].to_numpy(), _num_col(gn, col), args.bootstrap, args.seed)
        print(f"  {col}: PP3 {out['pp3']}  BP4 {out['bp4']}  ({time.time() - t1:.0f} s)", flush=True)
        return out

    print("calibration on the ClinVar 2019 set:")
    cal = {name: run(c19, col) for name, col in SCORES.items()}
    print("REVEL on the paper's own data (port check):")
    port = port_check(args.bootstrap, args.seed)
    print(f"  PP3 {port['pp3']}  BP4 {port['bp4']}")
    print("sensitivity: this project's training data:")
    sensitivity = {name: run(own, col) for name, col in (("model", "model"), ("REVEL", "revel"))}

    ours = cal["model"]
    rng = np.random.default_rng(args.seed)
    tests = {"ClinVar 2020 set": c20, "held-out benchmark, rare variants": bench,
             "held-out benchmark, rare variants labelled since 2023": bench[bench["new"]]}
    evaluation = {k: {**_counts(d), **evaluate(d, ours, rng)} for k, d in tests.items()}
    current_gn = np.array([acmg.computational(r)[1] or 0 for r in gn.to_dict("records")], float)
    population = {"model": shares(points(_num_col(gn, "model"), ours["pp3"], ours["bp4"])),
                  "current": shares(current_gn)}

    data = {"calibration": _counts(c19), "gnomad_set": {"variants": len(gn), "genes": int(gn["gene"].nunique())},
            "sensitivity": _counts(own), **{f"test: {k}": _counts(d) for k, d in tests.items()}}
    report = {"model_version": MODEL["version"], "model_trained": MODEL["trained"],
              "calibrated": dt.date.today().isoformat(), "bootstrap": args.bootstrap, "prior": ALPHA,
              "levels": {n: {"likelihood_ratio": round(lr_plus(n), 3),
                             "posterior_pathogenic": round(target_pathogenic(n), 4),
                             "posterior_benign": round(target_benign(n), 4)} for n in LEVELS},
              "data": data, "calibration": cal, "port_check": port, "sensitivity": sensitivity,
              "published": {k: {"pp3": dict(zip(LEVELS, v[0])), "bp4": dict(zip(LEVELS, v[1]))}
                            for k, v in PUBLISHED.items()},
              "evaluation": evaluation, "population": population}
    T.OUT.mkdir(exist_ok=True)
    (T.OUT / "acmg_calibration.json").write_text(json.dumps(report, indent=1))
    (T.OUT / "acmg_calibration.md").write_text(markdown(report))

    params = {"model_version": MODEL["version"], "model_trained": MODEL["trained"],
              "calibrated": report["calibrated"],
              "method": ("Pejaver et al. 2022 (ClinGen) local posterior probability, prior 0.0441, one-sided "
                         f"95% bound over {args.bootstrap} bootstrap samples; ClinVar 2019 set of the paper, "
                         "each variant scored by a model fitted without its gene"),
              "pp3": tuple((t, n) for n, t in sorted(ours["pp3"].items(), reverse=True) if t is not None),
              "bp4": tuple((t, -n) for n, t in sorted(ours["bp4"].items(), reverse=True) if t is not None)}
    PARAMS.write_text('"""ACMG PP3/BP4 thresholds of the model score. Generated by benchmarks/calibrate_acmg.py; '
                      'do not edit by hand."""\n\n'
                      f"CALIBRATION = {pprint.pformat(params, width=110, sort_dicts=False)}\n")
    print(f"wrote {PARAMS.relative_to(ROOT)} and results/acmg_calibration.md in {time.time() - t0:.0f} s")
    return 0


NAMES = {1: "Supporting", 2: "Moderate", 3: "Moderate+ (3 points)", 4: "Strong"}


def _t(v) -> str:
    return "—" if v is None else f"{v:.3f}"


def _evidence_table(ev: dict) -> list[str]:
    rows = [("Pathogenic with PP3", "pathogenic_pp3"), ("… PP3 Moderate or more", "pathogenic_pp3_moderate_or_more"),
            ("… PP3 Strong", "pathogenic_pp3_strong"), ("Pathogenic with BP4 (wrong way)", "pathogenic_bp4"),
            ("Benign with BP4", "benign_bp4"), ("… BP4 Moderate or more", "benign_bp4_moderate_or_more"),
            ("… BP4 Strong", "benign_bp4_strong"), ("Benign with PP3 (wrong way)", "benign_pp3")]
    L = ["| | Model, own calibration | Current (REVEL, else AlphaMissense) |", "|---|---|---|"]
    L += [f"| {label} | {ev['model'][k]:.1%} | {ev['current'][k]:.1%} |" for label, k in rows]
    L.append(f"| Mean points, pathogenic | {ev['model']['mean_points_pathogenic']:+.2f} | "
             f"{ev['current']['mean_points_pathogenic']:+.2f} |")
    L.append(f"| Mean points, benign | {ev['model']['mean_points_benign']:+.2f} | "
             f"{ev['current']['mean_points_benign']:+.2f} |")
    return L


def _band_table(bands: dict, other: dict) -> list[str]:
    L = ["| Band | Pathogenic | Benign | Likelihood ratio (95% CI) | Needed | Current scheme's ratio |",
         "|---|---|---|---|---|---|"]
    for k, b in bands.items():
        code = ("PP3 " if int(k) > 0 else "BP4 ") + NAMES[abs(int(k))]
        o = other.get(k)
        L.append(f"| {code} | {b['n_pathogenic']} | {b['n_benign']} | {b['lr']} ({b['ci'][0]}–{b['ci'][1]}) | "
                 f"{b['required']} | {o['lr'] if o else '—'} |")
    return L


def markdown(rep: dict) -> str:
    cal, d, m = rep["calibration"], rep["data"], rep["calibration"]["model"]
    c = d["calibration"]
    L = [f"# ACMG PP3/BP4 calibration of the MissensExplorer v{rep['model_version']} score", "",
         f"Generated by `benchmarks/calibrate_acmg.py` on {rep['calibrated']}: the method of Pejaver et al. "
         f"(AJHG 2022) that set ClinGen's thresholds for REVEL and other tools, {rep['bootstrap']:,} bootstrap "
         "samples.", "",
         f"- **Calibration data**: the paper's ClinVar 2019 set, {c['variants']:,} variants ({c['pathogenic']:,} "
         f"pathogenic, {c['benign']:,} benign, {c['genes']:,} genes), each scored by a model fitted without "
         "its gene.",
         f"- **gnomAD set**: {d['gnomad_set']['variants']:,} rare gnomAD v4 missense variants not in ClinVar, "
         f"from {d['gnomad_set']['genes']} genes.", "",
         "## Thresholds", "",
         "| Evidence | Likelihood ratio needed | Model score | Point estimate |", "|---|---|---|---|"]
    for n in sorted(LEVELS, reverse=True):
        L.append(f"| PP3 {NAMES[n]} | {rep['levels'][n]['likelihood_ratio']:.2f} | ≥ {_t(m['pp3'][n])} | "
                 f"{_t(m['pp3_point_estimate'][n])} |")
    for n in LEVELS:
        L.append(f"| BP4 {NAMES[n]} | 1/{rep['levels'][n]['likelihood_ratio']:.2f} | ≤ {_t(m['bp4'][n])} | "
                 f"{_t(m['bp4_point_estimate'][n])} |")
    L += ["", "Thresholds are the one-sided 95% bounds over the bootstrap samples (— : level not reached).", "",
          "## Checks of the method", "",
          "| Score and data | PP3 1 / 2 / 3 / 4 points | BP4 1 / 2 / 3 / 4 points |", "|---|---|---|"]

    def row(label, x):
        return (f"| {label} | {' / '.join(_t(x['pp3'][n]) for n in LEVELS)} | "
                f"{' / '.join(_t(x['bp4'][n]) for n in LEVELS)} |")
    pub = rep["published"]
    L += [row("REVEL, published (ClinGen)", pub["REVEL"]),
          row("REVEL, this port on the paper's data", rep["port_check"]),
          row("REVEL, this port, pipeline scores and this gnomAD set", cal["REVEL"]),
          row("AlphaMissense, published (Bergquist et al. 2025)", pub["AlphaMissense"]),
          row("AlphaMissense, this port, same data", cal["AlphaMissense"]),
          row("Model, calibrated on this project's own training data (sensitivity)", rep["sensitivity"]["model"]),
          row("REVEL, on this project's own training data (sensitivity)", rep["sensitivity"]["REVEL"]), ""]
    for name, ev in rep["evaluation"].items():
        L += [f"## Test: {name}", "",
              f"{ev['variants']:,} variants ({ev['pathogenic']:,} pathogenic, {ev['benign']:,} benign) in "
              f"{ev['genes']:,} genes.", "",
              "Likelihood ratio in each evidence band of the model (for BP4: benign over pathogenic):", ""]
        L += _band_table(ev["bands"], ev["current_bands"]) + ["", "Evidence given:", ""]
        L += _evidence_table(ev) + [""]
    L += ["## How often rare variants in the population get each level (gnomAD set)", "",
          "| Points | Model | Current |", "|---|---|---|"]
    for k in rep["population"]["model"]:
        L.append(f"| {int(k):+d} | {rep['population']['model'][k]:.1%} | {rep['population']['current'][k]:.1%} |")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    sys.exit(main())
