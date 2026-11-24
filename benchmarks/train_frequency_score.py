"""Fit the score with frequency: the model score and gnomAD grpmax combined, for ranking variants.

    pip install -e ".[benchmark]"
    python benchmarks/train_frequency_score.py

The specification was fixed earlier, as the exploratory row "frequency as a feature" of
benchmarks/compare_tools.py: a logistic regression (scikit-learn defaults) over
[logit of the model score, log10(grpmax + 1e-6)], fitted on the model's training set
(results/model_training_set.csv.gz) with the shipped model's scores. It is not retuned here, so
the held-out benchmark is not used a second time to choose it.

The score is for ranking a sample's variants. It is not ACMG evidence: frequency has its own
criteria there (BA1, BS1, PM2), and PP3/BP4 must come from a score that leaves frequency out.

Writes src/missensexplorer/frequency_params.py, tests/data/frequency_vectors.json and
benchmarks/results/frequency_score.{md,json}: AUC on the held-out benchmark (all variants, rare
ones, those labelled since 2023) against the model alone, frequency alone and MetaRNN.
"""

from __future__ import annotations

import datetime as dt
import json
import pprint
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import clinvar_benchmark as B  # noqa: E402
import compare_tools as C  # noqa: E402
import train_model as T  # noqa: E402

from missensexplorer import PredictOptions, model, predict  # noqa: E402
from missensexplorer.model_params import MODEL  # noqa: E402
from missensexplorer.predict import _num  # noqa: E402

PARAMS = ROOT / "src" / "missensexplorer" / "frequency_params.py"
VECTORS = ROOT / "tests" / "data" / "frequency_vectors.json"
AF_FLOOR = 1e-6


def fit() -> tuple[dict, pd.DataFrame]:
    from sklearn.linear_model import LogisticRegression

    df = pd.read_csv(T.OUT / "model_training_set.csv.gz")
    rows = df[T.FEATURES].astype(object).where(df[T.FEATURES].notna(), None).to_dict("records")
    df["p"] = [ms.p if (ms := model.score(r)) else 0.5 for r in rows]
    X = np.array([model.frequency_features(p, g) for p, g in zip(df["p"], df["grpmax"], strict=True)])
    lr = LogisticRegression().fit(X, df["label"].to_numpy())
    params = {
        "model_version": MODEL["version"],
        "model_trained": MODEL["trained"],
        "trained": dt.date.today().isoformat(),
        "n": len(df),
        "features": ["logit(model score)", "log10(gnomAD grpmax + 1e-6)"],
        "intercept": round(float(lr.intercept_[0]), 6),
        "coef": [round(float(c), 6) for c in lr.coef_[0]],
        "af_floor": AF_FLOOR,
    }
    sample = df.sample(30, random_state=1)
    vectors = [
        {
            "p": float(p),
            "grpmax": None if pd.isna(g) else float(g),
            "expected": float(lr.predict_proba([model.frequency_features(p, g)])[0, 1]),
        }
        for p, g in zip(sample["p"], sample["grpmax"], strict=True)
    ]
    return params, pd.DataFrame(vectors)


def auc_with_ci(y, s, ref, rng, n_boot=1000) -> dict:
    a, r = C.auc(y, s), C.auc(y, ref)
    d = []
    for _ in range(n_boot):
        i = rng.integers(0, len(y), len(y))
        if 0 < y[i].sum() < len(i):
            d.append(C.auc(y[i], s[i]) - C.auc(y[i], ref[i]))
    return {
        "auc": round(a, 4),
        "vs_model": round(a - r, 4),
        "vs_model_ci": [round(float(np.percentile(d, 2.5)), 4), round(float(np.percentile(d, 97.5)), 4)],
    }


def evaluate(params: dict) -> dict:
    cands = B.load_candidates(B.download())
    picked = B.sample(cands, 1000, 8, 42)
    res = predict(text="\n".join(p["input"] for p in picked), options=PredictOptions(**T.OPTIONS), progress=T._progress)
    print()
    T.check_sources(res)
    old = C.labelled_in_2023()
    rows = []
    for p, r in T._missense_rows(res, picked):
        ms = model.score(r)
        if ms is None:
            continue
        g = _num(r.get("gnomad_popmax_af"))
        rows.append(
            {
                "label": p["label"],
                "new": p["id"] not in old,
                "rare": (g or 0) < T.RARE,
                "model": ms.p,
                "with_frequency": model.with_frequency(ms.p, g, params),
                "frequency_alone": -np.log10((g or 0) + AF_FLOOR),
            }
        )
    df = pd.DataFrame(rows)
    metarnn = json.loads((T.OUT / "tool_comparison.json").read_text())["tools"]["MetaRNN"]["auc"]
    rng = np.random.default_rng(0)
    out = {}
    for name, key, m in (
        ("all variants", "all", np.ones(len(df), bool)),
        ("rare (< 0.1%)", "rare", df["rare"]),
        ("labelled since 2023", "new since 2023", df["new"]),
    ):
        d = df[np.asarray(m)]
        y, ref = d["label"].to_numpy(), d["model"].to_numpy()
        out[name] = {
            "variants": len(d),
            "pathogenic": int(y.sum()),
            "benign": int((y == 0).sum()),
            "model": round(C.auc(y, ref), 4),
            "with_frequency": auc_with_ci(y, d["with_frequency"].to_numpy(), ref, rng),
            "frequency_alone": auc_with_ci(y, d["frequency_alone"].to_numpy(), ref, rng),
            "MetaRNN (tool comparison)": metarnn[key]["auc"],
        }
    return out


def markdown(params: dict, ev: dict) -> str:
    b0, (b1, b2) = params["intercept"], params["coef"]
    L = [
        "# Score with frequency",
        "",
        f"Generated by `benchmarks/train_frequency_score.py` on {params['trained']}, for model "
        f"v{params['model_version']}.",
        "",
        "```",
        f"z = {b0:+.3f} {b1:+.3f}·logit(model score) {b2:+.3f}·log10(gnomAD grpmax + 1e-6)",
        "score with frequency = 1 / (1 + e^(−z))",
        "```",
        "",
        "AUC on the held-out benchmark (genes the model never saw); the difference from the model alone has a "
        "95% paired bootstrap interval.",
        "",
        "| Variants | Model alone | Score with frequency | Frequency alone | MetaRNN |",
        "|---|---|---|---|---|",
    ]
    for name, r in ev.items():
        wf, fa = r["with_frequency"], r["frequency_alone"]
        L.append(
            f"| {name} ({r['variants']:,}) | {r['model']:.3f} | {wf['auc']:.3f} "
            f"({wf['vs_model']:+.3f}; {wf['vs_model_ci'][0]:+.3f} to {wf['vs_model_ci'][1]:+.3f}) | "
            f"{fa['auc']:.3f} | {r['MetaRNN (tool comparison)']:.3f} |"
        )
    return "\n".join(L) + "\n"


def main() -> int:
    params, vectors = fit()
    PARAMS.write_text(
        '"""Score with frequency (model score + gnomAD grpmax). Generated by '
        'benchmarks/train_frequency_score.py; do not edit by hand."""\n\n'
        f"FREQUENCY = {pprint.pformat(params, width=110, sort_dicts=False)}\n"
    )
    VECTORS.write_text(json.dumps(vectors.to_dict("records"), indent=1))
    print("coefficients:", params["intercept"], params["coef"], flush=True)
    ev = evaluate(params)
    (T.OUT / "frequency_score.json").write_text(json.dumps({"params": params, "evaluation": ev}, indent=1))
    (T.OUT / "frequency_score.md").write_text(markdown(params, ev))
    print(markdown(params, ev))
    return 0


if __name__ == "__main__":
    sys.exit(main())
