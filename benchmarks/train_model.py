"""Train the missense model (logistic regression over predictor scores) on ClinVar.

    pip install -e ".[benchmark]"
    python benchmarks/train_model.py [--seed 7]

1. Uses the same ClinVar selection as clinvar_benchmark.py (germline missense SNVs on GRCh38,
   >= 2 review stars, P/LP or B/LB), at most 8 variants per gene and class.
2. Leaves out every gene of the benchmark sample, so the benchmark stays a held-out test. A variant
   is left out when the gene in its ClinVar name (NM_...(GENE):c.) or any of its ClinVar gene
   symbols is a benchmark gene: ClinVar lists overlapping genes together ('LDLR;MIR6886'), and up
   to v2.4 an exact match on that string let 63 variants of 36 benchmark genes into training.
3. Builds a training set that is enriched in rare benign variants:
     * every pathogenic variant;
     * every RARE benign variant (gnomAD grpmax < 0.1% or absent);
     * common benign variants, as many as needed to balance the two classes.
   Rare benign variants are the hard cases, and the ones the model is used for: common variants
   are handled by the frequency rule. In a plain random ClinVar sample only ~12% of the benign
   variants are rare (against ~98% of the pathogenic ones), so a model trained on it learns
   "benign" mostly from common variants. Version 2.2 was trained that way and called more rare
   benign variants damaging than the v2.1 vote did.
4. Annotates the set with the pipeline (without ClinVar, which is never an input of the verdict)
   and fits mean imputation -> standardisation -> L2 logistic regression on the scores of
   AlphaMissense, REVEL, CADD and popEVE, with the classes weighted equally. A missing score is
   filled with the training mean, which is 0 after standardisation: it adds nothing to the score
   (up to v2.4 the median was used, which added a small hidden term).
   Allele frequency is NOT a feature: it is applied as a separate, visible rule (BA1/BS1).
5. Chooses the cut-offs by 5-fold cross-validation split by gene, on rare variants (the training
   genes only; the benchmark plays no part):
     * damaging: the highest cut-off (at least 0.70) at which the model still calls at least as
       many rare pathogenic variants damaging as the v2.1 vote did. Within that limit it calls as
       few rare benign variants damaging as possible; in practice most rare variants of a patient
       are benign, so false alarms are the costlier error. Variants just below the cut-off become
       *uncertain*, never *benign*.
     * benign: the highest cut-off (at most 0.30) at which no more rare pathogenic variants are
       called benign than with the v2.1 vote.
6. Writes the coefficients to src/missensexplorer/model_params.py (pure Python, no scikit-learn
   needed at run time).
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
from missensexplorer.predict import _num  # noqa: E402

# EVE, SIFT and PolyPhen-2 add nothing on top of these four (cross-validated AUC 0.9766 vs 0.9762 on the
# v2.2 training set) and, being collinear with them, get coefficients of the "wrong" sign, which would make
# the per-predictor explanation misleading.
FEATURES = ["alphamissense", "revel", "cadd_phred", "popeve"]
RARE = 0.001  # gnomAD grpmax below 0.1% (or absent)
MIN_DAMAGING, MAX_BENIGN = 0.70, 0.30
PARAMS = ROOT / "src" / "missensexplorer" / "model_params.py"
OUT = Path(__file__).resolve().parent / "results"
OPTIONS = {"clinvar": False, "protein_context": False, "sequences": False, "workers": 6}


def pipeline(c: float = 1.0):
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    # No missing-value indicators: a missing score must never count as evidence either way. The mean
    # (0 after standardisation) makes a missing score contribute exactly nothing to the log-odds.
    return make_pipeline(SimpleImputer(strategy="mean"), StandardScaler(), LogisticRegression(C=c, max_iter=5000))


def _progress(s, d, n):
    print(f"  {s}: {d}/{n}", end="\r", flush=True)


def _missense_rows(res, picked: list[dict]) -> list[tuple[dict, dict]]:
    by_input = defaultdict(list)
    for r in res.table.to_dict("records"):
        by_input[r["input"]].append(r)
    out = []
    for p in picked:
        cand = [r for r in by_input.get(p["input"], []) if "missense_variant" in str(r["consequence"])]
        same_gene = [r for r in cand if r["gene"] == p["gene"]]
        r = (same_gene or cand or [None])[0]
        if r is not None:
            out.append((p, r))
    return out


def frequencies(picked: list[dict]) -> dict[str, float | None]:
    """GnomAD grpmax of each variant from the fast VEP pass (no predictor plugins: quick)."""
    res = predict(
        text="\n".join(p["input"] for p in picked),
        options=PredictOptions(predictors=False, **OPTIONS),
        progress=_progress,
    )
    print()
    check_sources(res)
    return {p["id"]: _num(r.get("gnomad_popmax_af")) for p, r in _missense_rows(res, picked)}


def enrich(picked: list[dict], grpmax: dict) -> list[dict]:
    """Every pathogenic and every rare benign variant, then common benign ones up to balance.

    `picked` is in the seeded random order of B.sample, so the common benign variants taken are a
    random subset.
    """
    known = [p for p in picked if p["id"] in grpmax]
    path = [p for p in known if p["label"] == 1]
    rare_benign = [p for p in known if p["label"] == 0 and (grpmax[p["id"]] or 0) < RARE]
    common_benign = [p for p in known if p["label"] == 0 and (grpmax[p["id"]] or 0) >= RARE]
    return path + rare_benign + common_benign[: max(0, len(path) - len(rare_benign))]


def annotate(picked: list[dict]) -> tuple[pd.DataFrame, dict]:
    res = predict(text="\n".join(p["input"] for p in picked), options=PredictOptions(**OPTIONS), progress=_progress)
    print()
    check_sources(res)
    rows = [
        {
            "clinvar_id": p["id"],
            "gene": p["gene"],
            "gene_key": B.gene_of(p.get("name", ""), p["gene"]),
            "label": p["label"],
            "year": p["year"],
            **{f: r.get(f) for f in FEATURES},
            "grpmax": _num(r.get("gnomad_popmax_af")),
            "v21_call": B.scheme_call(r, B.V21, 0.7, 0.3, 0.0, True),
        }
        for p, r in _missense_rows(res, picked)
    ]
    df = pd.DataFrame(rows)
    for f in [*FEATURES, "grpmax"]:
        df[f] = pd.to_numeric(df[f], errors="coerce")
    return df, res.run


def check_sources(res) -> None:
    """Never train or benchmark on a partial annotation (e.g. a dropped connection mid-run)."""
    failed = {k: v for k, v in res.sources.items() if v.startswith("failed")}
    if failed:
        raise SystemExit(
            f"annotation incomplete, not using it: {failed}. Run again: answers already "
            "obtained are cached, so only the missing ones are fetched."
        )


def class_weights(y: np.ndarray) -> np.ndarray:
    """Equal total weight for the two classes (the score is a ranking from balanced classes)."""
    w = np.empty(len(y))
    for lab in (0, 1):
        w[y == lab] = 0.5 * len(y) / max(1, (y == lab).sum())
    return w


def cv_predictions(df: pd.DataFrame, seed: int) -> np.ndarray:
    """Out-of-fold model scores, 5-fold cross-validation split by gene."""
    from sklearn.model_selection import GroupKFold

    X, y, g = df[FEATURES], df["label"].values, df["gene_key"].values
    p = np.zeros(len(df))
    for tr, te in GroupKFold(n_splits=5, shuffle=True, random_state=seed).split(X, y, g):
        m = pipeline().fit(X.iloc[tr], y[tr], logisticregression__sample_weight=class_weights(y[tr]))
        p[te] = m.predict_proba(X.iloc[te])[:, 1]
    return p


def rates(calls: np.ndarray, y: np.ndarray) -> dict:
    dec = calls != "uncertain"
    return {
        "n_pathogenic": int((y == 1).sum()),
        "n_benign": int((y == 0).sum()),
        "sensitivity": round(float(np.mean(calls[y == 1] == "damaging")), 4),
        "specificity": round(float(np.mean(calls[y == 0] == "benign")), 4),
        "pathogenic_called_benign": round(float(np.mean(calls[y == 1] == "benign")), 4),
        "benign_called_damaging": round(float(np.mean(calls[y == 0] == "damaging")), 4),
        "coverage": round(float(np.mean(dec)), 4),
        "accuracy_when_decisive": round(float(np.mean((calls[dec] == "damaging") == (y[dec] == 1))), 4),
    }


def model_calls(p: np.ndarray, hi: float, lo: float) -> np.ndarray:
    return np.where(p >= hi, "damaging", np.where(p <= lo, "benign", "uncertain"))


def choose_cut_offs(p: np.ndarray, df: pd.DataFrame) -> tuple[float, float]:
    """Cut-offs that keep the v2.1 vote's sensitivity on rare variants and miss no more of them (CV)."""
    rare = (df["grpmax"].fillna(0) < RARE).values
    y, pr = df["label"].values[rare], p[rare]
    v21 = rates(df["v21_call"].fillna("uncertain").values[rare], y)
    keeps_sensitivity = [
        h for h in np.round(np.arange(MIN_DAMAGING, 0.99, 0.01), 2) if np.mean(pr[y == 1] >= h) >= v21["sensitivity"]
    ]
    hi = max(keeps_sensitivity, default=MIN_DAMAGING)
    lo = next(
        lv
        for lv in np.round(np.arange(MAX_BENIGN, 0.0, -0.01), 2)
        if np.mean(pr[y == 1] <= lv) <= v21["pathogenic_called_benign"]
    )
    return float(hi), float(lo)


def cv_report(df: pd.DataFrame, p: np.ndarray, hi: float, lo: float) -> dict:
    from sklearn.metrics import roc_auc_score

    y = df["label"].values
    rare = (df["grpmax"].fillna(0) < RARE).values
    v21 = df["v21_call"].fillna("uncertain").values
    out = {}
    for name, m in (("all", np.ones(len(df), bool)), ("rare", rare)):
        out[name] = {
            "auc": round(float(roc_auc_score(y[m], p[m])), 4),
            "model": rates(model_calls(p[m], hi, lo), y[m]),
            "model_at_0.7/0.3": rates(model_calls(p[m], 0.7, 0.3), y[m]),
            "v2.1 vote": rates(v21[m], y[m]),
        }
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-gene", type=int, default=8)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args(argv)
    t0 = time.time()

    cands = B.load_candidates(B.download())
    bench = B.sample(list(cands), 1000, 8, 42)  # the benchmark sample (same arguments as its defaults)
    held_out = set().union(*(B.gene_keys(p) for p in bench))
    pool = [c for c in cands if not (B.gene_keys(c) & held_out)]
    picked = B.sample(pool, len(pool), args.per_gene, args.seed)
    print(
        f"{len(pool)} candidates outside the {len(held_out)} benchmark genes; {len(picked)} with at most "
        f"{args.per_gene} per gene and class",
        flush=True,
    )
    grpmax = frequencies(picked)
    picked = enrich(picked, grpmax)
    n_rare_benign = sum(p["label"] == 0 and (grpmax[p["id"]] or 0) < RARE for p in picked)
    print(
        f"training set: {sum(p['label'] for p in picked)} pathogenic, {n_rare_benign} rare benign, "
        f"{sum(p['label'] == 0 for p in picked) - n_rare_benign} common benign",
        flush=True,
    )

    df, run = annotate(picked)
    print(f"{len(df)} annotated as missense", flush=True)
    coverage = df[FEATURES].notna().mean()
    print("feature coverage:", coverage.round(3).to_dict(), flush=True)
    if coverage["alphamissense"] < 0.9 or coverage["revel"] < 0.9:
        raise SystemExit("AlphaMissense/REVEL coverage below 90%: predictor annotation is incomplete")
    p = cv_predictions(df, args.seed)
    hi, lo = choose_cut_offs(p, df)
    cv = cv_report(df, p, hi, lo)
    print(f"cut-offs {hi}/{lo}; 5-fold CV by gene:", json.dumps(cv, indent=1), flush=True)

    y = df["label"].values
    model = pipeline().fit(df[FEATURES], y, logisticregression__sample_weight=class_weights(y))
    imp, sc, lr = model[0], model[1], model[2]
    rare = (df["grpmax"].fillna(0) < RARE).values
    params = {
        "version": "2.5",
        "trained": dt.date.today().isoformat(),
        "training": {
            "source": B.URL,
            "selection": "germline missense SNVs, GRCh38, >= 2 review stars, "
            f"P/LP vs B/LB, at most {args.per_gene} per gene and class, benchmark genes excluded (by the gene of the "
            "ClinVar name and every ClinVar gene symbol); every pathogenic "
            "and every rare benign variant (gnomAD grpmax < 0.1%), common benign variants up to "
            "balance; classes weighted equally",
            "n": len(df),
            "pathogenic": int(y.sum()),
            "rare_benign": int(((y == 0) & rare).sum()),
            "genes": int(df["gene_key"].nunique()),
            "seed": args.seed,
            "annotation": run.get("annotation"),
            "cut_off_rule": "rare variants, CV by gene: highest damaging cut-off >= 0.70 that keeps the "
            "v2.1 vote's sensitivity; highest benign cut-off <= 0.30 that calls no more "
            "pathogenic variants benign than the v2.1 vote",
            "cv_by_gene": cv,
        },
        "features": FEATURES,
        "fill": [round(float(x), 6) for x in imp.statistics_],
        "mean": [round(float(x), 6) for x in sc.mean_],
        "scale": [round(float(x), 6) for x in sc.scale_],
        "coef": [round(float(x), 6) for x in lr.coef_[0]],
        "intercept": round(float(lr.intercept_[0]), 6),
        "damaging": hi,
        "benign": lo,
    }
    PARAMS.write_text(
        '"""Missense model coefficients. Generated by benchmarks/train_model.py; do not edit by hand."""\n\n'
        f"MODEL = {pprint.pformat(params, width=110, sort_dicts=False)}\n"
    )
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "model_training_set.csv.gz", index=False)
    # reference predictions so the pure-Python scorer can be checked against scikit-learn
    ref = df.sample(n=min(40, len(df)), random_state=args.seed)
    ref = ref.assign(expected=model.predict_proba(ref[FEATURES])[:, 1])
    vectors = [
        {f: (None if pd.isna(r[f]) else float(r[f])) for f in FEATURES} | {"expected": float(r["expected"])}
        for _, r in ref.iterrows()
    ]
    (ROOT / "tests" / "data" / "model_vectors.json").write_text(json.dumps(vectors, indent=1))
    print(f"wrote {PARAMS.relative_to(ROOT)} in {time.time() - t0:.0f} s")
    print(json.dumps({k: params[k] for k in ("features", "coef", "intercept", "damaging", "benign")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
