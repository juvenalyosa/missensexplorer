"""Missense model: logistic regression over the predictor scores.

Trained by benchmarks/train_model.py on expert-reviewed ClinVar variants (benchmark genes left out)
and stored as plain coefficients in model_params.py, so no machine-learning library is needed here.

The score is a 0-1 ranking learned from balanced classes (as for REVEL): it is not the probability
that a variant is pathogenic in a given patient. Each predictor's contribution (log-odds, positive =
towards damaging) is reported so the call can be explained.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .model_params import MODEL

NAMES = {"alphamissense": "AlphaMissense", "revel": "REVEL", "cadd_phred": "CADD", "eve": "EVE",
         "popeve": "popEVE", "sift_score": "SIFT", "polyphen_score": "PolyPhen"}
REQUIRED = ("alphamissense", "revel")  # the model is used only when at least one of these has a score
DAMAGING, BENIGN = MODEL["damaging"], MODEL["benign"]


def _num(x) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


@dataclass
class ModelScore:
    p: float
    contributions: list[tuple[str, float]]  # (predictor, log-odds), largest absolute value first


def score(row: dict) -> ModelScore | None:
    feats = MODEL["features"]
    vals = {f: _num(row.get(f)) for f in feats}
    if all(vals[f] is None for f in REQUIRED):
        return None
    z = MODEL["intercept"]
    contrib = []
    for f, med, mu, sd, c in zip(feats, MODEL["median"], MODEL["mean"], MODEL["scale"], MODEL["coef"]):
        v = vals[f]
        term = c * ((med if v is None else v) - mu) / sd
        z += term
        if v is not None:
            contrib.append((NAMES[f], term))
    contrib.sort(key=lambda x: -abs(x[1]))
    p = 1.0 / (1.0 + math.exp(-z))
    return ModelScore(p, contrib)


def describe(ms: ModelScore, top: int = 3) -> str:
    """'model score 0.93 (REVEL +1.9, CADD +0.8, AlphaMissense +0.7)'."""
    parts = [f"{n} {t:+.1f}" for n, t in ms.contributions[:top] if abs(t) >= 0.05]
    return f"model score {ms.p:.2f}" + (f" ({', '.join(parts)})" if parts else "")
