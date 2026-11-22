"""Train the missense model (logistic regression over predictor scores) on ClinVar.

    pip install -e ".[benchmark]"
    python benchmarks/train_model.py [--per-class 3000] [--seed 7]

1. Uses the same ClinVar selection as clinvar_benchmark.py (germline missense SNVs on GRCh38,
   >= 2 review stars, P/LP or B/LB).
2. Leaves out every gene of the benchmark sample, so the benchmark stays a held-out test.
3. Samples a balanced training set (at most 8 variants per gene and class) and annotates it with
   the pipeline (without ClinVar, which is never an input of the verdict).
4. Fits median imputation -> standardisation -> L2 logistic regression on the scores of
   AlphaMissense, REVEL, CADD and popEVE.
   Allele frequency is NOT a feature: it is applied as a separate, visible rule (BA1/BS1).
5. Reports a 5-fold cross-validated AUC split by gene, then writes the coefficients to
   src/missensexplorer/model_params.py (pure Python, no scikit-learn needed at run time).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import pprint
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import clinvar_benchmark as B  # noqa: E402
from missensexplorer import PredictOptions, predict  # noqa: E402

# EVE, SIFT and PolyPhen-2 add nothing on top of these four (cross-validated AUC 0.9766 vs 0.9762 on the
# training set) and, being collinear with them, get coefficients of the "wrong" sign, which would make the
# per-predictor explanation misleading.
FEATURES = ["alphamissense", "revel", "cadd_phred", "popeve"]
PARAMS = ROOT / "src" / "missensexplorer" / "model_params.py"
OUT = Path(__file__).resolve().parent / "results"


def pipeline(c: float = 1.0):
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    # No missing-value indicators: a missing score must never count as evidence either way.
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                         LogisticRegression(C=c, max_iter=5000))


def annotate(picked: list[dict]) -> pd.DataFrame:
    res = predict(text="\n".join(p["input"] for p in picked),
                  options=PredictOptions(clinvar=False, protein_context=False, sequences=False, workers=6),
                  progress=lambda s, d, n: print(f"  {s}: {d}/{n}", end="\r", flush=True))
    print()
    check_sources(res)
    by_input = defaultdict(list)
    for r in res.table.to_dict("records"):
        by_input[r["input"]].append(r)
    rows = []
    for p in picked:
        cand = [r for r in by_input.get(p["input"], []) if "missense_variant" in str(r["consequence"])]
        same_gene = [r for r in cand if r["gene"] == p["gene"]]
        r = (same_gene or cand or [None])[0]
        if r is not None:
            rows.append({"clinvar_id": p["id"], "gene": p["gene"], "label": p["label"], "year": p["year"],
                         **{f: r.get(f) for f in FEATURES}})
    df = pd.DataFrame(rows)
    for f in FEATURES:
        df[f] = pd.to_numeric(df[f], errors="coerce")
    return df, res.run


def check_sources(res) -> None:
    """Never train or benchmark on a partial annotation (e.g. a dropped connection mid-run)."""
    failed = {k: v for k, v in res.sources.items() if v.startswith("failed")}
    if failed:
        raise SystemExit(f"annotation incomplete, not using it: {failed}. Run again: answers already "
                         "obtained are cached, so only the missing ones are fetched.")


def cv_auc(df: pd.DataFrame, seed: int) -> dict:
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold

    X, y, g = df[FEATURES], df["label"].values, df["gene"].values
    p = np.zeros(len(df))
    for tr, te in GroupKFold(n_splits=5, shuffle=True, random_state=seed).split(X, y, g):
        p[te] = pipeline().fit(X.iloc[tr], y[tr]).predict_proba(X.iloc[te])[:, 1]
    hi, lo = p >= 0.7, p <= 0.3
    dec = hi | lo
    return {"auc": round(float(roc_auc_score(y, p)), 4),
            "sensitivity": round(float(np.mean(hi[y == 1])), 4), "specificity": round(float(np.mean(lo[y == 0])), 4),
            "coverage": round(float(np.mean(dec)), 4),
            "accuracy_when_decisive": round(float(np.mean((hi[dec]).astype(int) == y[dec])), 4)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-class", type=int, default=3000)
    ap.add_argument("--per-gene", type=int, default=8)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args(argv)
    t0 = time.time()

    cands = B.load_candidates(B.download())
    bench = B.sample(list(cands), 1000, 8, 42)            # the benchmark sample (same arguments as its defaults)
    held_out = {p["gene"] for p in bench}
    pool = [c for c in cands if c["gene"] not in held_out]
    picked = B.sample(pool, args.per_class, args.per_gene, args.seed)
    print(f"{len(pool)} candidates outside the {len(held_out)} benchmark genes; sampled {len(picked)} "
          f"({sum(p['label'] for p in picked)} pathogenic)", flush=True)

    df, run = annotate(picked)
    print(f"{len(df)} annotated as missense", flush=True)
    coverage = df[FEATURES].notna().mean()
    print("feature coverage:", coverage.round(3).to_dict(), flush=True)
    if coverage["alphamissense"] < 0.9 or coverage["revel"] < 0.9:
        raise SystemExit("AlphaMissense/REVEL coverage below 90%: predictor annotation is incomplete")
    cv = cv_auc(df, args.seed)
    print("5-fold CV by gene:", cv, flush=True)

    model = pipeline().fit(df[FEATURES], df["label"].values)
    imp, sc, lr = model[0], model[1], model[2]
    params = {
        "version": "2.2",
        "trained": dt.date.today().isoformat(),
        "training": {"source": B.URL, "selection": "germline missense SNVs, GRCh38, >= 2 review stars, "
                     "P/LP vs B/LB, at most %d per gene and class, benchmark genes excluded" % args.per_gene,
                     "n": int(len(df)), "pathogenic": int(df["label"].sum()), "genes": int(df["gene"].nunique()),
                     "seed": args.seed, "annotation": run.get("annotation"), "cv_by_gene": cv},
        "features": FEATURES,
        "median": [round(float(x), 6) for x in imp.statistics_],
        "mean": [round(float(x), 6) for x in sc.mean_],
        "scale": [round(float(x), 6) for x in sc.scale_],
        "coef": [round(float(x), 6) for x in lr.coef_[0]],
        "intercept": round(float(lr.intercept_[0]), 6),
        "damaging": 0.7,
        "benign": 0.3,
    }
    PARAMS.write_text(
        '"""Missense model coefficients. Generated by benchmarks/train_model.py; do not edit by hand."""\n\n'
        f"MODEL = {pprint.pformat(params, width=110, sort_dicts=False)}\n")
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "model_training_set.csv.gz", index=False)
    # reference predictions so the pure-Python scorer can be checked against scikit-learn
    ref = df.sample(n=min(40, len(df)), random_state=args.seed)
    ref = ref.assign(expected=model.predict_proba(ref[FEATURES])[:, 1])
    vectors = [{f: (None if pd.isna(r[f]) else float(r[f])) for f in FEATURES} | {"expected": float(r["expected"])}
               for _, r in ref.iterrows()]
    (ROOT / "tests" / "data" / "model_vectors.json").write_text(json.dumps(vectors, indent=1))
    print(f"wrote {PARAMS.relative_to(ROOT)} in {time.time() - t0:.0f} s")
    print(json.dumps({k: params[k] for k in ("features", "coef", "intercept")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
