"""Compare MissensExplorer with other missense predictors on the held-out ClinVar benchmark.

    python benchmarks/compare_tools.py

1. Takes the benchmark sample of clinvar_benchmark.py (seed 42; genes never used to train the
   MissensExplorer model) and the pipeline's own scores and calls.
2. Fetches the scores of other published predictors for the same variants from dbNSFP, through the
   Ensembl VEP REST API (cached like every VEP answer). When dbNSFP lists one score per transcript,
   the most damaging one is used, as dbNSFP's rank scores do.
3. Reports for every tool the AUC (threshold-free) with a 95% bootstrap interval, on all, rare and
   "new since 2023" variants, both on the variants each tool scores and on the subset that every
   tool scores; and, for the tools ClinGen has calibrated, the calls at those PP3/BP4 thresholds.
4. Flags the tools trained on clinical variant labels (ClinVar, HGMD, UniProt disease variants).
   They are favoured on any ClinVar benchmark: some of the test variants may have been in their
   training data. The "new since 2023" subset limits that advantage: it keeps the variants whose
   pathogenic/benign label did not exist yet in the January 2023 ClinVar release (it does not remove
   variants that were in HGMD before, which REVEL's training used).
5. Also reports the variants absent from gnomAD, those present below 0.01% (grpmax), and both together,
   where allele frequency carries little information: there, tools that use frequency lose most of their
   advantage if it came from frequency (the frequency-only row shows how much is left). Writes the
   per-variant scores to tool_comparison_variants.csv.gz.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import httpx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import clinvar_benchmark as B

from missensexplorer import PredictOptions, acmg, model, predict
from missensexplorer.api import Http, ensembl_release
from missensexplorer.cache import Cache, default_cache_dir
from missensexplorer.inputs import parse_text
from missensexplorer.model_params import MODEL
from missensexplorer.predict import _num
from missensexplorer.vep import Vep

OUT = Path(__file__).resolve().parent / "results"
ARCHIVE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2023/variant_summary_2023-01.txt.gz"
RARE = 0.001

# name -> (source, direction, clinical training labels, uses allele frequency).
# source "row:<column>" = the pipeline's own column; otherwise a dbNSFP field.
# direction +1: higher = more damaging; -1: lower = more damaging.
# Training labels and AF use from the tools' papers (see benchmarks/results/tool_comparison.md).
TOOLS = {
    "AlphaMissense": ("row:alphamissense", 1, "none (population)", False),
    "REVEL": ("row:revel", 1, "HGMD", False),
    "CADD": ("row:cadd_phred", 1, "none (simulated)", False),
    "EVE": ("row:eve", 1, "none", False),
    "popEVE": ("row:popeve", -1, "none (population)", False),
    "SIFT": ("row:sift_score", -1, "none", False),
    "PolyPhen-2": ("row:polyphen_score", 1, "UniProt", False),
    "ClinPred": ("ClinPred_score", 1, "ClinVar", True),
    "MetaRNN": ("MetaRNN_score", 1, "ClinVar", True),
    "BayesDel (no AF)": ("BayesDel_noAF_score", 1, "ClinVar, UniProt", False),
    "BayesDel (with AF)": ("BayesDel_addAF_score", 1, "ClinVar, UniProt", True),
    "VARITY_R": ("VARITY_R_score", 1, "ClinVar, HGMD", False),
    "VARITY_R (leave-one-out)": ("VARITY_R_LOO_score", 1, "ClinVar, HGMD", False),
    "MutPred2": ("MutPred2_score", 1, "HGMD, SwissVar", False),
    "gMVP": ("gMVP_score", 1, "ClinVar, HGMD, UniProt", False),
    "MVP": ("MVP_score", 1, "ClinVar, HGMD, UniProt", False),
    "MetaSVM": ("MetaSVM_score", 1, "UniProt", True),
    "MetaLR": ("MetaLR_score", 1, "UniProt", True),
    "DEOGEN2": ("DEOGEN2_score", 1, "UniProt", False),
    "PHACTboost": ("PHACTboost_score", 1, "ClinVar", False),
    "MutFormer": ("MutFormer_score", 1, "HGMD", False),
    "MutationTaster2021": ("MutationTaster_score", 1, "ClinVar, HGMD (+ look-up)", False),
    "FATHMM-XF": ("fathmm-XF_coding_score", 1, "HGMD", False),
    "MPC": ("MPC_score", 1, "ClinVar (402 variants)", False),
    "PrimateAI": ("PrimateAI_score", 1, "none (population)", False),
    "ESM1b": ("ESM1b_score", -1, "none", False),
    "LIST-S2": ("LIST-S2_score", 1, "few parameters", False),
    "PROVEAN": ("PROVEAN_score", -1, "none", False),
    "MutationAssessor": ("MutationAssessor_score", 1, "none", False),
    "SIFT4G": ("SIFT4G_score", -1, "none", False),
    "DANN": ("DANN_score", 1, "none (simulated)", False),
    "Eigen": ("Eigen-raw_coding", 1, "none", False),
    "phyloP (100 vertebrates)": ("phyloP100way_vertebrate", 1, "none", False),
}

# ClinGen-calibrated thresholds at the supporting level, on each tool's own scale: PP3 from `pp3`
# (score >= pp3, or <= for lower-is-damaging tools) and BP4 from `bp4` (None = no BP4 interval).
# Pejaver et al., AJHG 2022, Table 2; Bergquist et al., Genet Med 2025, Table 1 (AlphaMissense, ESM1b,
# VARITY_R, the latter on leave-one-out scores as in the paper). Calibrated on older versions of some
# tools (CADD before v1.7, ESM1b before dbNSFP's 2024 update), so treat those rows as indicative.
CALIBRATED: dict[str, tuple[float, float | None]] = {
    "REVEL": (0.644, 0.290),
    "AlphaMissense": (0.792, 0.169),
    "BayesDel (no AF)": (0.13, -0.18),
    "MutPred2": (0.737, 0.391),
    "PrimateAI": (0.790, 0.483),
    "VARITY_R (leave-one-out)": (0.675, 0.251),
    "ESM1b": (-10.7, -6.3),
    "CADD": (25.3, 22.7),
    "PolyPhen-2": (0.978, 0.113),
    "SIFT": (0.001, 0.080),
    "MPC": (1.360, None),
    "phyloP (100 vertebrates)": (7.367, 1.879),
}


def dbnsfp_value(tc: dict, field: str, direction: int):
    v = tc.get(field.lower())
    if v is None or v == "invalid_field":
        return None
    vals = [_num(x) for x in str(v).split(",")] if isinstance(v, str) else [_num(v)]
    vals = [x for x in vals if x is not None]
    if not vals:
        return None
    return max(vals) if direction > 0 else min(vals)


def fetch_dbnsfp(rows: list[dict], workers: int = 6) -> dict[str, dict]:
    """DbNSFP scores of every benchmark variant (ClinVar id -> {tool: score})."""
    # MutationTaster_pred is not scored but stays in the request: the cache key depends on the field list
    fields = [t[0] for t in TOOLS.values() if not t[0].startswith("row:")] + ["MutationTaster_pred"]
    params = {"dbNSFP": ",".join(fields), "pick_allele_gene": 1, "mane": 1}
    http, cache = Http(), Cache()
    release = ensembl_release(http, cache, "GRCh38")
    vep = Vep(http, cache, "GRCh38", workers=workers, release=release)
    query = {}
    for x in rows:
        v = parse_text(x["input"]).variants[0]
        query[x["id"]] = v.query
    res = vep.annotate(
        {"region": list(query.values())},
        params,
        lambda s, d, n: print(f"  dbNSFP: {d}/{n}", end="\r", flush=True),
        "dbNSFP",
    )
    print()
    failed = sum("error" in v for v in res.values())
    if failed:
        raise SystemExit(f"{failed} dbNSFP queries failed; run again (answers obtained so far are cached)")
    out = {}
    for x in rows:
        val = res[("region", query[x["id"]])]
        tcs = [
            tc
            for r in val.get("results", [])
            for tc in r.get("transcript_consequences") or []
            if "missense_variant" in (tc.get("consequence_terms") or [])
        ]
        tc = next((t for t in tcs if t.get("gene_symbol") == x["gene"]), tcs[0] if tcs else {})
        scores = {}
        for name, (src, direction, _, _) in TOOLS.items():
            if src.startswith("row:"):
                continue
            scores[name] = dbnsfp_value(tc, src, direction)
        out[x["id"]] = scores
    return out


def labelled_in_2023() -> set[str]:
    """ClinVar variation ids already labelled pathogenic/benign in the January 2023 release."""
    d = default_cache_dir() / "benchmark"
    f = d / "variant_summary_2023-01.txt.gz"
    if not f.exists():
        print("Downloading the January 2023 ClinVar release ...", flush=True)
        with (
            httpx.stream("GET", ARCHIVE, timeout=600, follow_redirects=True) as r,
            open(f.with_suffix(".part"), "wb") as fh,
        ):
            for chunk in r.iter_bytes(1 << 20):
                fh.write(chunk)
        f.with_suffix(".part").rename(f)
    ids = set()
    with gzip.open(f, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(header)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            first = c[ix["ClinicalSignificance"]].split(";")[0].split(",")[0].strip()
            if first in B.PATHOGENIC or first in B.BENIGN:
                ids.add(c[ix["VariationID"]])
    return ids


def frequency_stack():
    """Fit the exploratory model score + allele frequency stack.

    The model score plus log10 gnomAD grpmax as a second feature, fitted on the model's training set
    (benchmarks/results/model_training_set.csv.gz). Shows how much of the gap to the tools that take
    allele frequency as input is the frequency itself. Not used by MissensExplorer.
    """
    import pandas as pd
    from sklearn.linear_model import LogisticRegression

    df = pd.read_csv(OUT / "model_training_set.csv.gz")
    feats = MODEL["features"]
    X = [stack_features({f: (None if pd.isna(r[f]) else r[f]) for f in feats}, r["grpmax"]) for _, r in df.iterrows()]
    return LogisticRegression().fit(np.array(X), df["label"].values)


def stack_features(row: dict, grpmax) -> list[float]:
    ms = model.score(row)
    p = min(max(ms.p if ms else 0.5, 1e-6), 1 - 1e-6)
    af = _num(grpmax) or 0.0
    return [float(np.log(p / (1 - p))), float(np.log10(af + 1e-6))]


def auc(y: np.ndarray, s: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score

    return float(roc_auc_score(y, s)) if len(set(y)) == 2 else float("nan")


def bootstrap(y: np.ndarray, s: np.ndarray, ref: np.ndarray, rng: np.random.Generator, n: int):
    """95% intervals of the AUC and of its difference from the reference score (paired resamples)."""
    a, d = [], []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if 0 < y[i].sum() < len(i):
            a.append(auc(y[i], s[i]))
            d.append(a[-1] - auc(y[i], ref[i]))

    def pct(v: list[float]) -> list[float]:
        if not v:
            return [float("nan"), float("nan")]
        return [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]

    return pct(a), pct(d)


def calls(score: float | None, pp3: float, bp4: float) -> str | None:
    if score is None or np.isnan(score):
        return None
    return "damaging" if score >= pp3 else "benign" if score <= bp4 else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--bootstrap", type=int, default=1000)
    args = ap.parse_args(argv)

    cands = B.load_candidates(B.download())
    picked = B.sample(cands, 1000, 8, args.seed)
    res = predict(
        text="\n".join(p["input"] for p in picked),
        options=PredictOptions(clinvar=False, protein_context=False, sequences=False, workers=6),
        progress=lambda s, d, n: print(f"  {s}: {d}/{n}", end="\r", flush=True),
    )
    print()
    failed = {k: v for k, v in res.sources.items() if v.startswith("failed")}
    if failed:
        raise SystemExit(f"annotation incomplete: {failed}. Run again (cached answers are kept).")
    by_input = defaultdict(list)
    for r in res.table.to_dict("records"):
        by_input[r["input"]].append(r)
    rows = []
    for p in picked:
        cand = [r for r in by_input.get(p["input"], []) if "missense_variant" in str(r["consequence"])]
        r = ([r for r in cand if r["gene"] == p["gene"]] or cand or [None])[0]
        if r is not None:
            rows.append({**p, "row": r})
    ext = fetch_dbnsfp(rows)
    old = labelled_in_2023()

    ours = "MissensExplorer v" + MODEL["version"]
    ours_freq = ours + " + frequency rule (BA1)"
    freq_only = "gnomAD frequency alone (reference)"
    ours_af = ours + " + frequency as a feature (exploratory)"
    names = [ours, ours_freq, *list(TOOLS), freq_only, ours_af]
    labels = {
        ours: "ClinVar (other genes)",
        ours_freq: "ClinVar (other genes)",
        freq_only: "—",
        ours_af: "ClinVar (other genes)",
        **{k: v[2] for k, v in TOOLS.items()},
    }
    uses_af = {ours: False, ours_freq: True, freq_only: True, ours_af: True, **{k: v[3] for k, v in TOOLS.items()}}
    stack = frequency_stack()
    S = np.full((len(rows), len(names)), np.nan)
    grpmax = np.array([_num(x["row"].get("gnomad_popmax_af")) or 0.0 for x in rows])
    for i, x in enumerate(rows):
        ms = model.score(x["row"])
        S[i, 0] = ms.p if ms else np.nan
        common_var = grpmax[i] >= acmg.BA1_AF and not acmg.ba1_exception(x["row"])
        S[i, 1] = np.nan if ms is None else (0.0 if common_var else ms.p)
        for j, (name, (src, direction, _, _)) in enumerate(TOOLS.items(), start=2):
            v = _num(x["row"].get(src[4:])) if src.startswith("row:") else ext[x["id"]].get(name)
            S[i, j] = np.nan if v is None else direction * v
        S[i, -2] = -grpmax[i]
        S[i, -1] = np.nan if ms is None else stack.decision_function([stack_features(x["row"], grpmax[i])])[0]
    y = np.array([x["label"] for x in rows])
    groups = {
        "all": np.ones(len(rows), bool),
        "rare": grpmax < RARE,
        "new since 2023": np.array([x["id"] not in old for x in rows]),
        "absent from gnomAD": grpmax == 0,
        "grpmax below 0.01%": (grpmax > 0) & (grpmax < 1e-4),
        "below 0.01% or absent": grpmax < 1e-4,
    }
    cover = np.mean(~np.isnan(S), axis=0)
    common_cols = [j for j in range(len(names)) if cover[j] >= 0.9]
    common = np.all(~np.isnan(S[:, common_cols]), axis=1)

    rng = np.random.default_rng(args.seed)
    report = {
        "n": len(rows),
        "genes": len({x["gene"] for x in rows}),
        "run": res.run,
        "model": {k: MODEL[k] for k in ("version", "trained", "damaging", "benign")},
        "groups": {g: {"n_pathogenic": int(y[m].sum()), "n_benign": int((y[m] == 0).sum())} for g, m in groups.items()},
        "common_subset": {"tools": [names[j] for j in common_cols], "n": int(common.sum())},
        "tools": {},
    }
    for j, name in enumerate(names):
        entry = {"labels": labels[name], "uses_af": uses_af[name], "coverage": round(float(cover[j]), 4), "auc": {}}
        for g, m in groups.items():
            ok = m & ~np.isnan(S[:, j]) & ~np.isnan(S[:, 0])
            ci, diff = bootstrap(y[ok], S[ok, j], S[ok, 0], rng, args.bootstrap)
            e = {
                "auc": round(auc(y[ok], S[ok, j]), 4),
                "ci": ci,
                "n": int(ok.sum()),
                "ours_same_variants": round(auc(y[ok], S[ok, 0]), 4),
                "difference_ci": diff,
            }
            if j in common_cols:
                cm = m & common
                e["common_subset"] = round(auc(y[cm], S[cm, j]), 4)
            entry["auc"][g] = e
        if name in CALIBRATED:
            pp3, bp4 = CALIBRATED[name]
            direction = TOOLS[name][1]
            bp4 = -np.inf * direction if bp4 is None else bp4
            entry["clingen_calls"] = {
                g: B.call_metrics([(calls(S[i, j], direction * pp3, direction * bp4), y[i]) for i in np.flatnonzero(m)])
                for g, m in groups.items()
            }
        report["tools"][name] = entry
        print(
            f"{name:40s} cov {cover[j]:.3f} "
            + " ".join(f"{g}: {v['auc']:.3f} [{v['ci'][0]:.3f}-{v['ci'][1]:.3f}]" for g, v in entry["auc"].items()),
            flush=True,
        )
    report["ours_calls"] = {
        g: B.call_metrics([(B.verdict_call(rows[i]["row"]["prediction"]), y[i]) for i in np.flatnonzero(m)])
        for g, m in groups.items()
    }

    OUT.mkdir(exist_ok=True)
    table = [
        {
            "clinvar_id": x["id"],
            "gene": x["gene"],
            "label": x["label"],
            "year": x.get("year"),
            "new_since_2023": x["id"] not in old,
            "grpmax": grpmax[i],
            "prediction": x["row"].get("prediction"),
            **{n: (None if np.isnan(S[i, j]) else float(S[i, j])) for j, n in enumerate(names)},
        }
        for i, x in enumerate(rows)
    ]
    import pandas as pd

    pd.DataFrame(table).to_csv(OUT / "tool_comparison_variants.csv.gz", index=False)
    (OUT / "tool_comparison.json").write_text(json.dumps(report, indent=1, default=str))
    (OUT / "tool_comparison.md").write_text(markdown(report))
    print(f"wrote {OUT / 'tool_comparison.md'}")
    return 0


def markdown(rep: dict) -> str:
    g = rep["groups"]
    L = [
        "# MissensExplorer and other missense predictors",
        "",
        f"Same held-out benchmark as [clinvar_benchmark.md](clinvar_benchmark.md): {rep['n']} expert-reviewed "
        f"ClinVar missense variants from {rep['genes']} genes, none of them used to train the MissensExplorer "
        f"model. Other tools' scores come from dbNSFP through the Ensembl VEP REST API "
        f"({rep['run'].get('annotation')}). Produced by `benchmarks/compare_tools.py`.",
        "",
        "- **AUC** (area under the ROC curve) is threshold-free: the chance that a random pathogenic variant "
        "scores higher than a random benign one. 95% bootstrap intervals in brackets.",
        "- **Δ vs ours**: the tool's AUC minus MissensExplorer's on the same variants, with a paired bootstrap "
        "interval. An interval that excludes 0 is a real difference on this benchmark.",
        f"- **Rare**: gnomAD grpmax < 0.1% ({g['rare']['n_pathogenic']} pathogenic, {g['rare']['n_benign']} benign). "
        f"**New since 2023**: variants whose pathogenic/benign label was not yet in the January 2023 ClinVar "
        f"release ({g['new since 2023']['n_pathogenic']} pathogenic, {g['new since 2023']['n_benign']} benign); "
        "tools trained before 2023 cannot have seen these labels.",
        "- **Clinical training labels**: the disease-variant databases a tool was trained on. A tool trained on "
        "ClinVar may have seen some of these very variants (all such tools were trained on ClinVar releases from "
        "2020 or earlier); HGMD and UniProt overlap ClinVar's pathogenic variants. *none* = no clinical labels.",
        "- **Uses allele frequency**: the tool takes gnomAD/1000 Genomes frequencies as input. ClinVar curators "
        "use the same frequencies to call variants benign (ACMG BA1/BS1), so these tools look better on ClinVar "
        "than the protein-level evidence alone would. MissensExplorer applies frequency as a separate, visible "
        "rule instead; the line *+ frequency rule (BA1)* scores the default verdict, and *gnomAD frequency alone* "
        "shows how far frequency goes by itself.",
        "- **MutationTaster2021** calls variants found in ClinVar or HGMD disease-causing by look-up, so on this "
        "benchmark it partly reads the answer.",
        "",
    ]
    for grp in ("all", "rare", "new since 2023", "absent from gnomAD", "grpmax below 0.01%"):
        L += [
            f"## {grp[0].upper() + grp[1:]} variants",
            "",
            "| Tool | Clinical training labels | Uses allele frequency | Coverage | AUC | Δ vs ours |",
            "|---|---|---|---|---|---|",
        ]
        items = sorted(
            [kv for kv in rep["tools"].items() if kv[1]["auc"][grp]["auc"] == kv[1]["auc"][grp]["auc"]],
            key=lambda kv: -kv[1]["auc"][grp]["auc"],
        )
        for name, t in items:
            a = t["auc"][grp]
            d = a["auc"] - a["ours_same_variants"]
            delta = (
                "—"
                if name.startswith("MissensExplorer v") and "frequency" not in name
                else (f"{d:+.3f} [{a['difference_ci'][0]:+.3f}, {a['difference_ci'][1]:+.3f}]")
            )
            bold = "**" if name.startswith("MissensExplorer") else ""
            L.append(
                f"| {bold}{name}{bold} | {t['labels']} | {'yes' if t['uses_af'] else 'no'} | "
                f"{100 * t['coverage']:.0f}% | "
                f"{a['auc']:.3f} [{a['ci'][0]:.3f}–{a['ci'][1]:.3f}] | {delta} |"
            )
        L.append("")
    cs = rep["common_subset"]
    L += [
        f"## Same variants for every tool ({cs['n']} variants scored by all tools with ≥ 90% coverage)",
        "",
        "| Tool | All | Rare | New since 2023 |",
        "|---|---|---|---|",
    ]
    items = [(n, t) for n, t in rep["tools"].items() if "common_subset" in t["auc"]["all"]]
    for name, t in sorted(items, key=lambda kv: -kv[1]["auc"]["all"]["common_subset"]):
        L.append(
            f"| {name} | "
            + " | ".join(f"{t['auc'][grp]['common_subset']:.3f}" for grp in ("all", "rare", "new since 2023"))
            + " |"
        )
    L.append("")
    cal = [(n, t) for n, t in rep["tools"].items() if "clingen_calls" in t]
    if cal:
        L += [
            "## Calls at ClinGen-calibrated thresholds (PP3/BP4, supporting level), rare variants",
            "",
            "| Tool | Pathogenic called damaging | Benign called benign | Pathogenic called benign | "
            "Benign called damaging | Decisive calls | Decisive calls right |",
            "|---|---|---|---|---|---|---|",
        ]
        rows = [(f"**MissensExplorer v{rep['model']['version']} (default verdict)**", rep["ours_calls"]["rare"])]
        rows += [(n, t["clingen_calls"]["rare"]) for n, t in cal]
        for n, m in rows:
            L.append(
                f"| {n} | {100 * m['sensitivity']:.1f}% | {100 * m['specificity']:.1f}% | "
                f"{100 * m['pathogenic_called_benign']:.1f}% | {100 * m['benign_called_damaging']:.1f}% | "
                f"{100 * m['coverage']:.1f}% | {100 * m['accuracy_when_decisive']:.1f}% |"
            )
        L.append("")
    return "\n".join(L)


if __name__ == "__main__":
    sys.exit(main())
