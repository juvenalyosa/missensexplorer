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

from .frequency_params import FREQUENCY
from .model_params import MODEL

NAMES = {"alphamissense": "AlphaMissense", "revel": "REVEL", "cadd_phred": "CADD", "eve": "EVE",
         "popeve": "popEVE", "sift_score": "SIFT", "polyphen_score": "PolyPhen"}
REQUIRED = ("alphamissense", "revel")  # the model is used only when at least one of these has a score
DAMAGING, BENIGN = MODEL["damaging"], MODEL["benign"]
# The score with frequency is fitted on this model's scores: it applies only to the model it was fitted for.
FREQUENCY_FITTED = (FREQUENCY["model_version"], FREQUENCY["model_trained"]) == (MODEL["version"], MODEL["trained"])


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


def score(row: dict, params: dict | None = None) -> ModelScore | None:
    """Model score of a row; `params` defaults to the shipped model (other values: benchmarks)."""
    m = params or MODEL
    feats = m["features"]
    vals = {f: _num(row.get(f)) for f in feats}
    if all(vals[f] is None for f in REQUIRED):
        return None
    z = m["intercept"]
    contrib = []
    for f, med, mu, sd, c in zip(feats, m["median"], m["mean"], m["scale"], m["coef"]):
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


def frequency_features(p: float, grpmax: float | None) -> list[float]:
    """[logit of the model score, log10 of gnomAD grpmax + floor]; a variant absent from gnomAD's main
    ancestry groups counts as 0 (log10 of the floor, about one allele in a million)."""
    p = min(max(p, 1e-6), 1 - 1e-6)
    af = _num(grpmax) or 0.0
    return [math.log(p / (1 - p)), math.log10(af + FREQUENCY["af_floor"])]


def with_frequency(p: float | None, grpmax: float | None, params: dict | None = None) -> float | None:
    """Model score combined with the gnomAD grpmax frequency (0-1), for ranking variants.

    Not for ACMG PP3/BP4: frequency has its own criteria there (BA1, BS1, PM2), and ClinGen asks for
    computational evidence that does not use it. None when the model gave no score."""
    f = params or FREQUENCY
    if p is None or (params is None and not FREQUENCY_FITTED):
        return None
    z = f["intercept"] + sum(c * x for c, x in zip(f["coef"], frequency_features(p, grpmax)))
    return 1.0 / (1.0 + math.exp(-z))

