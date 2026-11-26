"""Benchmark the MissensExplorer verdict against expert-reviewed ClinVar missense variants.

    python benchmarks/clinvar_benchmark.py [--per-class 1000] [--seed 42]

1. Downloads ClinVar's variant_summary.txt.gz (once, ~450 MB, cached in
   ~/.cache/missensexplorer/benchmark/).
2. Samples germline missense SNVs on GRCh38 with >= 2 review stars (multiple submitters
   without conflicts, expert panel or practice guideline), labelled Pathogenic/Likely
   pathogenic or Benign/Likely benign; at most 8 per gene so a few genes do not dominate.
3. Runs the pipeline on the genomic coordinates WITHOUT ClinVar (it is never an input of the verdict).
4. Reports, for the verdict (model + frequency rule), the model alone, the previous v2.2 model, the v2.1 vote, the
   ClinGen PP3/BP4 call and each predictor: sensitivity (pathogenic called damaging), specificity
   (benign called benign), how many stay uncertain, and the AUC; for all variants, rare variants
   and variants evaluated since 2024.
5. The missense model is trained by train_model.py on other genes only, so this is a held-out test.

Caveat: some predictors used ClinVar-like data in training or calibration (REVEL, and
AlphaMissense for its thresholds), so the figures are optimistic; the "evaluated since
2024" subset is reported separately as a less circular check.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import httpx

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from missensexplorer import PredictOptions, acmg, model, predict
from missensexplorer.cache import default_cache_dir
from missensexplorer.model_params import MODEL
from missensexplorer.predict import WEIGHTS, _modern_present, _num, _votes, classify, prediction_status

URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
GOOD_REVIEW = {"criteria provided, multiple submitters, no conflicts", "reviewed by expert panel", "practice guideline"}
PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}
MISSENSE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)$")
NAME_GENE = re.compile(r"^[A-Z]{2}_\d+(?:\.\d+)?\(([^)]+)\):")
CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
OUT = Path(__file__).resolve().parent / "results"


def download() -> Path:
    d = default_cache_dir() / "benchmark"
    d.mkdir(parents=True, exist_ok=True)
    f = d / "variant_summary.txt.gz"
    if f.exists() and f.stat().st_size > 100_000_000:
        return f
    print("Downloading ClinVar variant_summary.txt.gz ...", flush=True)
    with httpx.stream("GET", URL, timeout=600, follow_redirects=True) as r, open(f.with_suffix(".part"), "wb") as fh:
        for chunk in r.iter_bytes(1 << 20):
            fh.write(chunk)
    f.with_suffix(".part").rename(f)
    return f


def load_candidates(path: Path) -> list[dict]:
    out, seen = [], set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["Assembly"]] != "GRCh38" or f[ix["Type"]] != "single nucleotide variant":
                continue
            sig, rev = f[ix["ClinicalSignificance"]], f[ix["ReviewStatus"]]
            if rev not in GOOD_REVIEW or (sig not in PATHOGENIC and sig not in BENIGN):
                continue
            if "germline" not in f[ix["OriginSimple"]]:
                continue
            m = MISSENSE.search(f[ix["Name"]])
            if not m or m.group(3) in ("Ter", m.group(1)):
                continue
            chrom, pos = f[ix["Chromosome"]], f[ix["PositionVCF"]]
            ref, alt = f[ix["ReferenceAlleleVCF"]], f[ix["AlternateAlleleVCF"]]
            vid = f[ix["VariationID"]]
            if chrom not in CHROMS or not pos.isdigit() or len(ref) != 1 or len(alt) != 1 or vid in seen:
                continue
            seen.add(vid)
            year = re.search(r"(\d{4})", f[ix["LastEvaluated"]] or "")
            out.append(
                {
                    "id": vid,
                    "gene": f[ix["GeneSymbol"]],
                    "label": 1 if sig in PATHOGENIC else 0,
                    "input": f"{chrom}:{pos} {ref}>{alt}",
                    "review": rev,
                    "year": int(year.group(1)) if year else None,
                    "name": f[ix["Name"]],
                }
            )
    return out


def gene_of(name: str, symbol: str) -> str:
    """The variant's gene: the one in its ClinVar name ('NM_000546.6(TP53):c.817C>T ...'), else its first symbol.

    ClinVar's GeneSymbol lists every gene a variant overlaps ('LDLR;MIR6886', 'CCDC40;GAA'), so the
    transcript's gene in the name is the reliable key for grouping variants by gene.
    """
    m = NAME_GENE.match(name or "")
    return m.group(1) if m else (symbol or "").split(";")[0]


def gene_keys(c: dict) -> set[str]:
    """Every gene a candidate can belong to: the gene of its name and each of its ClinVar symbols.

    Used to keep held-out genes out of training: a variant is excluded when ANY of its keys is a
    benchmark gene (an exact match on the GeneSymbol string missed 'LDLR;MIR6886' for LDLR).
    """
    keys = {s.strip() for s in (c.get("gene") or "").split(";") if s.strip() and " " not in s.strip()}
    keys.add(gene_of(c.get("name", ""), c.get("gene", "")))
    return {k for k in keys if k}


def sample(cands: list[dict], per_class: int, per_gene: int, seed: int) -> list[dict]:
    rng = random.Random(seed)
    rng.shuffle(cands)
    picked, n_gene = [], defaultdict(int)
    count = {0: 0, 1: 0}
    for c in cands:
        if count[c["label"]] >= per_class or n_gene[(c["gene"], c["label"])] >= per_gene:
            continue
        picked.append(c)
        n_gene[(c["gene"], c["label"])] += 1
        count[c["label"]] += 1
    return picked


# --------------------------------------------------------------------------- metrics


def auc(scores: list[tuple[float, int]]) -> float | None:
    """Mann-Whitney AUC (ties count half)."""
    pos = [s for s, y in scores if y == 1]
    neg = [s for s, y in scores if y == 0]
    if not pos or not neg:
        return None
    allv = sorted((s, y) for s, y in scores)
    ranks, i = {}, 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        for k in range(i, j + 1):
            ranks.setdefault(allv[k][0], (i + j) / 2 + 1)
        i = j + 1
    r_pos = sum(ranks[s] for s in pos)
    return (r_pos - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def call_metrics(calls: list[tuple[str | None, int]]) -> dict:
    """calls: (call in {'damaging','benign',None/'uncertain'}, label)."""
    P = [c for c, y in calls if y == 1]
    B = [c for c, y in calls if y == 0]
    dec = [(c, y) for c, y in calls if c in ("damaging", "benign")]
    correct = sum((c == "damaging") == (y == 1) for c, y in dec)
    return {
        "n_pathogenic": len(P),
        "n_benign": len(B),
        "sensitivity": sum(c == "damaging" for c in P) / max(1, len(P)),
        "specificity": sum(c == "benign" for c in B) / max(1, len(B)),
        "pathogenic_called_benign": sum(c == "benign" for c in P) / max(1, len(P)),
        "benign_called_damaging": sum(c == "damaging" for c in B) / max(1, len(B)),
        "coverage": len(dec) / max(1, len(calls)),
        "accuracy_when_decisive": correct / max(1, len(dec)),
    }


def verdict_call(label: str | None) -> str | None:
    st = prediction_status(label)
    return {"Damaging": "damaging", "Benign / silent": "benign"}.get(st)


def scheme_call(
    row: dict, weights: dict, hi: float, lo: float, min_weight: float = 0.0, fallback: bool = False
) -> str | None:
    dmg, ben, _ = _votes(row)
    if fallback and not _modern_present(row):  # weak predictors only when nothing better exists
        weights = {**weights, "SIFT": 0.5, "PolyPhen": 0.5}
    dmg = [x for x in dmg if weights.get(x, 0) > 0]
    ben = [x for x in ben if weights.get(x, 0) > 0]
    wd, wb = sum(weights[x] for x in dmg), sum(weights[x] for x in ben)
    if wd + wb == 0 or wd + wb < min_weight:
        return None
    frac = wd / (wd + wb)
    return "damaging" if frac >= hi else "benign" if frac <= lo else None


V21 = {"AlphaMissense": 1.0, "REVEL": 1.0, "EVE": 1.0, "CADD": 0.75, "PolyPhen": 0, "SIFT": 0}

# The v2.2 model (trained on a plain random ClinVar sample), kept to compare with the current one.
V22 = {
    "features": ["alphamissense", "revel", "cadd_phred", "popeve"],
    "fill": [0.3227, 0.3915, 24.5, -3.469],  # v2.2 filled missing scores with the training median
    "mean": [0.4726, 0.455065, 22.068322, -3.591189],
    "scale": [0.385334, 0.341274, 9.409868, 0.764894],
    "coef": [1.150214, 1.752403, 0.837612, -0.827833],
    "intercept": 0.010747,
    "damaging": 0.7,
    "benign": 0.3,
}


def params_call(row: dict, params: dict) -> str | None:
    """Call of a model with the given parameters; the v2.1 vote when the model does not apply."""
    ms = model.score(row, params)
    if ms is None:
        return scheme_call(row, V21, 0.7, 0.3, 0.0, True)
    return "damaging" if ms.p >= params["damaging"] else "benign" if ms.p <= params["benign"] else None


def vote_fraction(row: dict) -> float | None:
    """Weighted share of damaging votes, as the v2.1 verdict computed it (its effect score)."""
    dmg, ben, _ = _votes(row)
    w = V21 if _modern_present(row) else {**V21, "SIFT": 0.5, "PolyPhen": 0.5}
    wd, wb = sum(w[x] for x in dmg), sum(w[x] for x in ben)
    return None if wd + wb == 0 else wd / (wd + wb)


def predictor_calls(row: dict) -> dict[str, str | None]:
    dmg, ben, _ = _votes(row)
    out = {}
    for name in WEIGHTS:
        out[name] = "damaging" if name in dmg else "benign" if name in ben else None
    return out


def fmt(x):
    return "—" if x is None else f"{100 * x:.1f}%" if isinstance(x, float) and x <= 1 else str(x)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-class", type=int, default=1000)
    ap.add_argument("--per-gene", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)

    t0 = time.time()
    cands = load_candidates(download())
    picked = sample(cands, args.per_class, args.per_gene, args.seed)
    print(
        f"{len(cands)} eligible ClinVar missense variants; sampled {len(picked)} "
        f"({sum(p['label'] for p in picked)} pathogenic) from {len({p['gene'] for p in picked})} genes",
        flush=True,
    )

    res = predict(
        text="\n".join(p["input"] for p in picked),
        options=PredictOptions(clinvar=False, protein_context=False, sequences=False, workers=6),
        progress=lambda s, d, n: print(f"  {s}: {d}/{n}", end="\r", flush=True),
    )
    print()
    failed = {k: v for k, v in res.sources.items() if v.startswith("failed")}
    if failed:
        raise SystemExit(f"annotation incomplete, not reporting it: {failed}. Run again (cached answers are kept).")
    by_input = defaultdict(list)
    for r in res.table.to_dict("records"):
        by_input[r["input"]].append(r)
    rows = []
    for p in picked:
        cand = [r for r in by_input.get(p["input"], []) if "missense_variant" in str(r["consequence"])]
        same_gene = [r for r in cand if r["gene"] == p["gene"]]
        r = (same_gene or cand or [None])[0]
        if r is None:
            continue
        rows.append({**p, "row": r})
    print(f"{len(rows)} of {len(picked)} variants annotated as missense by VEP", flush=True)

    def subset(items, cond):
        return [x for x in items if cond(x)]

    for x in rows:
        r = x["row"]
        x["model_label"] = classify(r, {"missense_variant"})[0]  # model without the frequency rule
        x["grpmax"] = _num(r.get("gnomad_popmax_af"))
        x["v21"] = scheme_call(r, V21, 0.7, 0.3, 0.0, True)
        x["v22"] = params_call(r, V22)
        pts = r.get("acmg_points")
        x["clingen"] = (
            None
            if not isinstance(pts, (int, float)) or pts != pts or pts == 0
            else ("damaging" if pts > 0 else "benign")
        )

    def with_freq(call, x, af):  # same rule as the pipeline, ClinGen BA1 exceptions included
        common = x["grpmax"] is not None and x["grpmax"] >= af and not acmg.ba1_exception(x["row"])
        return "benign" if common else call

    ver = f"v{MODEL['version']}"
    methods = {
        f"{ver} default: model + frequency rule (BA1, gnomAD grpmax ≥ 5%)": lambda x: verdict_call(
            x["row"]["prediction"]
        ),
        f"{ver} model alone (no frequency rule)": lambda x: verdict_call(x["model_label"]),
        f"{ver} model + frequency ≥ 1%": lambda x: with_freq(verdict_call(x["model_label"]), x, 0.01),
        "v2.2 model + frequency rule (BA1, 5%)": lambda x: with_freq(x["v22"], x, 0.05),
        "v2.1 vote": lambda x: x["v21"],
        "v2.1 vote + frequency rule (BA1, 5%)": lambda x: with_freq(x["v21"], x, 0.05),
        (
            "PP3/BP4 alone (model score, this project's ClinGen-method calibration)"
            if acmg.MODEL_CALIBRATED
            else "ClinGen PP3/BP4 alone (REVEL, else AlphaMissense)"
        ): lambda x: x["clingen"],
    }
    report: dict = {
        "clinvar_file": URL,
        "seed": args.seed,
        "sampled": len(picked),
        "annotated": len(rows),
        "genes": len({x["gene"] for x in rows}),
        "run": res.run,
        "model": {k: MODEL[k] for k in ("version", "trained", "training", "damaging", "benign")},
    }
    groups = {
        "all": rows,
        "rare (gnomAD grpmax < 0.1% or absent)": subset(rows, lambda x: (x["grpmax"] or 0) < 0.001),
        "evaluated since 2024": subset(rows, lambda x: (x["year"] or 0) >= 2024),
    }
    report["methods"], report["predictors"], report["auc"] = {}, {}, {}
    for gname, items in groups.items():
        report["methods"][gname] = {m: call_metrics([(f(x), x["label"]) for x in items]) for m, f in methods.items()}
        report["predictors"][gname] = {
            name: call_metrics([(predictor_calls(x["row"])[name], x["label"]) for x in items]) for name in WEIGHTS
        }
        aucs = {}
        for name, get in (
            (f"{ver} model score", lambda r: model.score(r).p if model.score(r) else None),
            ("v2.2 model score", lambda r: model.score(r, V22).p if model.score(r, V22) else None),
            ("v2.1 vote fraction", vote_fraction),
            ("AlphaMissense", lambda r: _num(r.get("alphamissense"))),
            ("REVEL", lambda r: _num(r.get("revel"))),
            ("CADD", lambda r: _num(r.get("cadd_phred"))),
            ("EVE", lambda r: _num(r.get("eve"))),
            ("popEVE", lambda r: None if _num(r.get("popeve")) is None else -_num(r.get("popeve"))),
            ("SIFT", lambda r: None if _num(r.get("sift_score")) is None else -_num(r.get("sift_score"))),
            ("PolyPhen-2", lambda r: _num(r.get("polyphen_score"))),
        ):
            sc = [(v, x["label"]) for x in items if (v := get(x["row"])) is not None]
            aucs[name] = {"auc": auc(sc), "coverage": len(sc) / max(1, len(items))}
        report["auc"][gname] = aucs
    report["model_cut_offs"] = {}
    rare_rows = groups["rare (gnomAD grpmax < 0.1% or absent)"]
    for hi, lo in sorted({(0.6, 0.4), (0.7, 0.3), (0.8, 0.2), (model.DAMAGING, model.BENIGN)}):
        name = f"{hi:.2f}/{lo:.2f}" + (" (default)" if (hi, lo) == (model.DAMAGING, model.BENIGN) else "")
        report["model_cut_offs"][name] = {
            g: call_metrics(
                [(params_call(x["row"], {**MODEL, "damaging": hi, "benign": lo}), x["label"]) for x in items]
            )
            for g, items in (("all", rows), ("rare", rare_rows))
        }
    report["seconds"] = round(time.time() - t0, 1)

    OUT.mkdir(exist_ok=True)
    import pandas as pd

    ms_of = {x["id"]: model.score(x["row"]) for x in rows}
    pd.DataFrame(
        [
            {
                "clinvar_id": x["id"],
                "gene": x["gene"],
                "gene_key": gene_of(x.get("name", ""), x["gene"]),
                "label": x["label"],
                "year": x.get("year"),
                "review": x.get("review"),
                "grpmax": x["grpmax"],
                "model_score": ms_of[x["id"]].p if ms_of[x["id"]] else None,
                "prediction": x["row"].get("prediction"),
                "model_label": x["model_label"],
                "v21": x["v21"],
                "v22": x["v22"],
                "acmg_pp3_bp4": x["row"].get("acmg_pp3_bp4"),
                "acmg_points": x["row"].get("acmg_points"),
                **{
                    f: x["row"].get(f)
                    for f in ("alphamissense", "revel", "cadd_phred", "popeve", "eve", "spliceai_max")
                },
            }
            for x in rows
        ]
    ).to_csv(OUT / "clinvar_benchmark_variants.csv.gz", index=False)
    (OUT / "clinvar_benchmark.json").write_text(json.dumps(report, indent=1, default=str))
    (OUT / "clinvar_benchmark.md").write_text(markdown(report))
    print((OUT / "clinvar_benchmark.md").read_text())
    return 0


def _metric_row(name: str, m: dict) -> str:
    return (
        f"| {name} | {fmt(m['sensitivity'])} | {fmt(m['specificity'])} | {fmt(m['pathogenic_called_benign'])} | "
        f"{fmt(m['benign_called_damaging'])} | {fmt(m['coverage'])} | {fmt(m['accuracy_when_decisive'])} |"
    )


HEAD = [
    "| Method | Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | "
    "Accuracy when decisive |",
    "|---|---|---|---|---|---|---|",
]


def markdown(rep: dict) -> str:
    tr = rep["model"]["training"]
    L = [
        "# ClinVar benchmark",
        "",
        f"{rep['annotated']} expert-reviewed (≥ 2 stars) germline missense variants from {rep['genes']} genes, "
        f"GRCh38, sampled with seed {rep['seed']} from ClinVar's variant_summary. ClinVar is not an input of the "
        f"verdict. Run: {rep['run'].get('tool')}, {rep['run'].get('annotation')}, {rep['run'].get('started_utc')}.",
        "",
        f"The missense model (v{rep['model']['version']}, trained {rep['model']['trained']}) was fitted on "
        f"{tr['n']} other ClinVar variants from {tr['genes']} genes; **no gene of this benchmark was used for "
        "training**, so these are held-out results (benchmarks/train_model.py). Its cut-offs "
        f"({rep['model']['damaging']}/{rep['model']['benign']}) were chosen by cross-validation on the training "
        "genes, not on this benchmark.",
        "",
        "Sensitivity = pathogenic variants called *damaging*; specificity = benign variants called *benign*; "
        "coverage = share of variants that get a decisive call (not *uncertain*/*not assessed*).",
        "",
    ]
    for g, ms in rep["methods"].items():
        n = next(iter(ms.values()))
        L += [f"## {g[0].upper() + g[1:]} (n = {n['n_pathogenic']} pathogenic, {n['n_benign']} benign)", "", *HEAD]
        L += [_metric_row(k, m) for k, m in ms.items()]
        a = rep["auc"][g]
        L += [
            "",
            "AUC: "
            + "; ".join(
                f"{k} {v['auc']:.3f}" + (f" (coverage {fmt(v['coverage'])})" if v["coverage"] < 0.995 else "")
                for k, v in a.items()
                if v["auc"] is not None
            )
            + ".",
            "",
        ]
        if g == "all":
            L += [
                "| Predictor alone | Sensitivity | Specificity | Coverage | Accuracy when decisive |",
                "|---|---|---|---|---|",
            ]
            for name, pm in rep["predictors"][g].items():
                L.append(
                    f"| {name} | {fmt(pm['sensitivity'])} | {fmt(pm['specificity'])} | {fmt(pm['coverage'])} | "
                    f"{fmt(pm['accuracy_when_decisive'])} |"
                )
            L.append("")
    for g in ("all", "rare"):
        L += [f"## Model cut-offs ({'all' if g == 'all' else 'rare'} variants, no frequency rule)", "", *HEAD]
        L += [_metric_row(k, m[g]) for k, m in rep["model_cut_offs"].items()] + [""]
    L += [
        "",
        "Caveats: REVEL was trained on disease-variant databases and AlphaMissense's thresholds were calibrated "
        "on ClinVar, so absolute figures are optimistic; the comparison between methods is the robust part. "
        "The frequency rule looks better here than it will in practice because ClinVar curators use the same "
        "gnomAD frequencies to classify variants as benign; it also mislabels pathogenic variants of recessive or "
        "founder conditions (that is why the default is 5% and not 1%).",
        "",
    ]
    return "\n".join(L)


if __name__ == "__main__":
    sys.exit(main())
