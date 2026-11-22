"""Benchmark the MissensExplorer verdict against expert-reviewed ClinVar missense variants.

    python benchmarks/clinvar_benchmark.py [--per-class 1000] [--seed 42]

1. Downloads ClinVar's variant_summary.txt.gz (once, ~450 MB, cached in
   ~/.cache/missensexplorer/benchmark/).
2. Samples germline missense SNVs on GRCh38 with >= 2 review stars (multiple submitters
   without conflicts, expert panel or practice guideline), labelled Pathogenic/Likely
   pathogenic or Benign/Likely benign; at most 8 per gene so a few genes do not dominate.
3. Runs the pipeline on the genomic coordinates WITHOUT ClinVar (it never votes anyway).
4. Reports, for the verdict and for each predictor: sensitivity (pathogenic called
   damaging), specificity (benign called benign), how many stay uncertain, and the AUC.
5. Compares alternative weights/cut-offs with 2-fold cross-validation.

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

from missensexplorer import PredictOptions, predict  # noqa: E402
from missensexplorer.cache import default_cache_dir  # noqa: E402
from missensexplorer.predict import WEIGHTS, _modern_present, _votes, prediction_status  # noqa: E402

URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
GOOD_REVIEW = {"criteria provided, multiple submitters, no conflicts", "reviewed by expert panel", "practice guideline"}
PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}
MISSENSE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)$")
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
            out.append({"id": vid, "gene": f[ix["GeneSymbol"]], "label": 1 if sig in PATHOGENIC else 0,
                        "input": f"{chrom}:{pos} {ref}>{alt}", "review": rev,
                        "year": int(year.group(1)) if year else None, "name": f[ix["Name"]]})
    return out


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
        "n_pathogenic": len(P), "n_benign": len(B),
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


def scheme_call(row: dict, weights: dict, hi: float, lo: float, min_weight: float = 0.0,
                fallback: bool = False) -> str | None:
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
    print(f"{len(cands)} eligible ClinVar missense variants; sampled {len(picked)} "
          f"({sum(p['label'] for p in picked)} pathogenic) from {len({p['gene'] for p in picked})} genes", flush=True)

    res = predict(text="\n".join(p["input"] for p in picked),
                  options=PredictOptions(clinvar=False, protein_context=False, sequences=False, workers=6),
                  progress=lambda s, d, n: print(f"  {s}: {d}/{n}", end="\r", flush=True))
    print()
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

    report: dict = {"clinvar_file": URL, "seed": args.seed, "sampled": len(picked), "annotated": len(rows),
                    "genes": len({x["gene"] for x in rows}), "run": res.run}
    groups = {"all": rows, "evaluated since 2024": subset(rows, lambda x: (x["year"] or 0) >= 2024)}
    report["verdict"] = {}
    report["predictors"] = {}
    report["auc"] = {}
    for gname, items in groups.items():
        report["verdict"][gname] = call_metrics([(verdict_call(x["row"]["prediction"]), x["label"]) for x in items])
        per = {}
        for name in WEIGHTS:
            per[name] = call_metrics([(predictor_calls(x["row"])[name], x["label"]) for x in items])
        report["predictors"][gname] = per
        aucs = {}
        for name, col, sign in (("Effect score (verdict)", "effect_score", 1), ("AlphaMissense", "alphamissense", 1),
                                ("REVEL", "revel", 1), ("CADD", "cadd_phred", 1), ("EVE", "eve", 1),
                                ("popEVE", "popeve", -1), ("SIFT", "sift_score", -1), ("PolyPhen-2", "polyphen_score", 1)):
            sc = [(sign * float(x["row"][col]), x["label"]) for x in items
                  if x["row"].get(col) is not None and x["row"][col] == x["row"][col]]
            aucs[name] = {"auc": auc(sc), "coverage": len(sc) / max(1, len(items))}
        report["auc"][gname] = aucs

    # ---- alternative schemes, 2-fold cross-validation by gene (no gene in both halves)
    genes = sorted({x["gene"] for x in rows})
    random.Random(args.seed).shuffle(genes)
    half = set(genes[: len(genes) // 2])
    folds = [[x for x in rows if x["gene"] in half], [x for x in rows if x["gene"] not in half]]
    modern = {"AlphaMissense": 1.0, "REVEL": 1.0, "EVE": 1.0, "CADD": 0.75, "PolyPhen": 0, "SIFT": 0}
    v20 = {"AlphaMissense": 1.0, "REVEL": 1.0, "EVE": 1.0, "CADD": 0.75, "PolyPhen": 0.5, "SIFT": 0.5}
    equal = {k: 1.0 for k in WEIGHTS}
    schemes = {}
    for wname, w, fb in (("v2.0 (SIFT/PolyPhen always vote)", v20, False), ("equal weights", equal, False),
                         ("without SIFT/PolyPhen", modern, False),
                         ("SIFT/PolyPhen only as fallback (v2.1 default)", modern, True)):
        for hi, lo in ((0.7, 0.3), (0.6, 0.4), (0.8, 0.2)):
            for mw in (0.0, 1.5):
                key = f"{wname}, cut-offs {hi:.2f}/{lo:.2f}" + (", min. weight 1.5" if mw else "")
                schemes[key] = [call_metrics([(scheme_call(x["row"], w, hi, lo, mw, fb), x["label"]) for x in f])
                                for f in folds]
    report["schemes_cv"] = schemes
    report["seconds"] = round(time.time() - t0, 1)

    OUT.mkdir(exist_ok=True)
    (OUT / "clinvar_benchmark.json").write_text(json.dumps(report, indent=1, default=str))
    (OUT / "clinvar_benchmark.md").write_text(markdown(report))
    print((OUT / "clinvar_benchmark.md").read_text())
    return 0


def markdown(rep: dict) -> str:
    L = ["# ClinVar benchmark", "",
         f"{rep['annotated']} expert-reviewed (≥ 2 stars) germline missense variants from {rep['genes']} genes, "
         f"GRCh38, sampled with seed {rep['seed']} from ClinVar's variant_summary. ClinVar is not used by the "
         f"verdict. Run: {rep['run'].get('tool')}, {rep['run'].get('annotation')}, {rep['run'].get('started_utc')}.",
         "", "Sensitivity = pathogenic variants called *damaging*; specificity = benign variants called *benign*; "
         "coverage = share of variants that get a decisive call (not *uncertain*/*not assessed*).", ""]
    for g, m in rep["verdict"].items():
        L += [f"## Verdict — {g} (n = {m['n_pathogenic']} pathogenic, {m['n_benign']} benign)", "",
              "| Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | Accuracy when decisive |",
              "|---|---|---|---|---|---|",
              f"| {fmt(m['sensitivity'])} | {fmt(m['specificity'])} | {fmt(m['pathogenic_called_benign'])} | "
              f"{fmt(m['benign_called_damaging'])} | {fmt(m['coverage'])} | {fmt(m['accuracy_when_decisive'])} |", ""]
        L += ["| Predictor | Sensitivity | Specificity | Coverage | Accuracy when decisive | AUC (score) |",
              "|---|---|---|---|---|---|"]
        aucs = rep["auc"][g]
        names = {"AlphaMissense": "AlphaMissense", "REVEL": "REVEL", "CADD": "CADD", "EVE": "EVE", "SIFT": "SIFT",
                 "PolyPhen": "PolyPhen-2"}
        for name, pm in rep["predictors"][g].items():
            a = aucs.get(names[name], {}).get("auc")
            L.append(f"| {names[name]} | {fmt(pm['sensitivity'])} | {fmt(pm['specificity'])} | {fmt(pm['coverage'])} | "
                     f"{fmt(pm['accuracy_when_decisive'])} | {'—' if a is None else f'{a:.3f}'} |")
        e = aucs["Effect score (verdict)"]["auc"]
        pe = aucs["popEVE"]["auc"]
        L += ["", f"AUC of the verdict's effect score: {'—' if e is None else f'{e:.3f}'}; "
                  f"popEVE: {'—' if pe is None else f'{pe:.3f}'} (coverage {fmt(aucs['popEVE']['coverage'])}).", ""]
    L += ["## Alternative weights and cut-offs (2-fold cross-validation, split by gene)", "",
          "| Scheme | Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | "
          "Accuracy when decisive |", "|---|---|---|---|---|---|---|"]
    for name, folds in rep["schemes_cv"].items():
        avg = {k: sum(f[k] for f in folds) / len(folds) for k in (
            "sensitivity", "specificity", "pathogenic_called_benign", "benign_called_damaging", "coverage",
            "accuracy_when_decisive")}
        L.append(f"| {name} | {fmt(avg['sensitivity'])} | {fmt(avg['specificity'])} | "
                 f"{fmt(avg['pathogenic_called_benign'])} | {fmt(avg['benign_called_damaging'])} | "
                 f"{fmt(avg['coverage'])} | {fmt(avg['accuracy_when_decisive'])} |")
    L += ["", "Caveat: REVEL was trained on disease-variant databases and AlphaMissense's thresholds were calibrated "
          "on ClinVar, so these figures are optimistic; the *evaluated since 2024* rows are a less circular check.", ""]
    return "\n".join(L)


if __name__ == "__main__":
    sys.exit(main())
