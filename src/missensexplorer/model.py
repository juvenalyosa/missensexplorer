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
from typing import Any, TypedDict, cast

from ._types import Row
from .frequency_params import FREQUENCY as _FREQUENCY
from .model_params import MODEL as _MODEL


class ModelParams(TypedDict):
    """What :func:`score` needs: the coefficients and cut-offs of one model (the benchmarks pass others).

    Attributes:
        features: Row columns used as predictors, in the order of the lists below.
        median: Value used for a missing predictor score.
        mean: Mean of each predictor (standardisation).
        scale: Standard deviation of each predictor (standardisation).
        coef: Logistic-regression coefficient of each standardised predictor.
        intercept: Logistic-regression intercept.
        damaging: Score at or above which a variant is called damaging.
        benign: Score at or below which a variant is called benign.
    """

    features: list[str]
    median: list[float]
    mean: list[float]
    scale: list[float]
    coef: list[float]
    intercept: float
    damaging: float
    benign: float


class ShippedModelParams(ModelParams):
    """Shape of ``model_params.MODEL``, written by benchmarks/train_model.py.

    Attributes:
        version: Model version.
        trained: Training date.
        training: Provenance and cross-validated performance of the training run.
    """

    version: str
    trained: str
    training: dict[str, Any]


class FrequencyParams(TypedDict):
    """Shape of ``frequency_params.FREQUENCY``, written by benchmarks/train_frequency_score.py.

    Attributes:
        model_version: Version of the model whose scores it was fitted on.
        model_trained: Training date of that model.
        trained: Date of the fit.
        n: Number of training variants.
        features: Description of the two features (see :func:`frequency_features`).
        intercept: Logistic-regression intercept.
        coef: Logistic-regression coefficients of the two features.
        af_floor: Frequency added to gnomAD grpmax before its log10.
    """

    model_version: str
    model_trained: str
    trained: str
    n: int
    features: list[str]
    intercept: float
    coef: list[float]
    af_floor: float


MODEL = cast(ShippedModelParams, _MODEL)
FREQUENCY = cast(FrequencyParams, _FREQUENCY)

NAMES = {
    "alphamissense": "AlphaMissense",
    "revel": "REVEL",
    "cadd_phred": "CADD",
    "eve": "EVE",
    "popeve": "popEVE",
    "sift_score": "SIFT",
    "polyphen_score": "PolyPhen",
}
REQUIRED = ("alphamissense", "revel")  # the model is used only when at least one of these has a score
DAMAGING, BENIGN = MODEL["damaging"], MODEL["benign"]
# The score with frequency is fitted on this model's scores: it applies only to the model it was fitted for.
FREQUENCY_FITTED = (FREQUENCY["model_version"], FREQUENCY["model_trained"]) == (MODEL["version"], MODEL["trained"])


def _num(x: Any) -> float | None:
    """Return ``x`` as a float, or None when it is missing, not a number or NaN."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


@dataclass
class ModelScore:
    """The model score of one variant, with what each predictor contributed.

    Attributes:
        p: Model score, 0-1 (a ranking, not a probability of pathogenicity; see the module docstring).
        contributions: (predictor, log-odds), largest absolute value first
    """

    p: float
    contributions: list[tuple[str, float]]


def score(row: Row, params: ModelParams | None = None) -> ModelScore | None:
    """Model score of a row; `params` defaults to the shipped model (other values: benchmarks).

    Args:
        row: One results row with the predictor columns (``alphamissense``, ``revel``, ...); a missing
            score is replaced by the training median.
        params: Model coefficients; None for :data:`MODEL`.

    Returns:
        The score, or None when none of :data:`REQUIRED` has a score.
    """
    m = params or MODEL
    feats = m["features"]
    vals = {f: _num(row.get(f)) for f in feats}
    if all(vals[f] is None for f in REQUIRED):
        return None
    z = m["intercept"]
    contrib = []
    for f, med, mu, sd, c in zip(feats, m["median"], m["mean"], m["scale"], m["coef"], strict=True):
        v = vals[f]
        term = c * ((med if v is None else v) - mu) / sd
        z += term
        if v is not None:
            contrib.append((NAMES[f], term))
    contrib.sort(key=lambda x: -abs(x[1]))
    p = 1.0 / (1.0 + math.exp(-z))
    return ModelScore(p, contrib)


def describe(ms: ModelScore, top: int = 3) -> str:
    """Return 'model score 0.93 (REVEL +1.9, CADD +0.8, AlphaMissense +0.7)'.

    Args:
        ms: The model score.
        top: Number of contributions shown at most (those of at least 0.05 log-odds).
    """
    parts = [f"{n} {t:+.1f}" for n, t in ms.contributions[:top] if abs(t) >= 0.05]
    return f"model score {ms.p:.2f}" + (f" ({', '.join(parts)})" if parts else "")


def frequency_features(p: float, grpmax: float | None) -> list[float]:
    """Return the features of the score with frequency.

    [logit of the model score, log10 of gnomAD grpmax + floor]; a variant absent from gnomAD's main
    ancestry groups counts as 0 (log10 of the floor, about one allele in a million).

    Args:
        p: Model score, 0-1 (clipped to [1e-6, 1 - 1e-6]).
        grpmax: gnomAD grpmax allele frequency; None (or NaN) when absent.
    """
    p = min(max(p, 1e-6), 1 - 1e-6)
    af = _num(grpmax) or 0.0
    return [math.log(p / (1 - p)), math.log10(af + FREQUENCY["af_floor"])]


def with_frequency(p: float | None, grpmax: float | None, params: FrequencyParams | None = None) -> float | None:
    """Model score combined with the gnomAD grpmax frequency (0-1), for ranking variants.

    Not for ACMG PP3/BP4: frequency has its own criteria there (BA1, BS1, PM2), and ClinGen asks for
    computational evidence that does not use it. None when the model gave no score.

    Args:
        p: Model score (0-1), or None.
        grpmax: gnomAD grpmax allele frequency; None (or NaN) when absent.
        params: Fitted coefficients; None for :data:`FREQUENCY`, which is used only if it was fitted on the
            shipped model (:data:`FREQUENCY_FITTED`; otherwise the result is None).
    """
    f = params or FREQUENCY
    if p is None or (params is None and not FREQUENCY_FITTED):
        return None
    z = f["intercept"] + sum(c * x for c, x in zip(f["coef"], frequency_features(p, grpmax), strict=True))
    return 1.0 / (1.0 + math.exp(-z))
