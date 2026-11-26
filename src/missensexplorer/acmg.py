"""ACMG/AMP evidence that can be read off the annotation: PP3/BP4 and population frequency.

PP3/BP4 (computational evidence) uses ONE score with calibrated intervals, as ClinGen recommends (no
counting of several tools). Two sources (PredictOptions.acmg_source):
  * "model" (default): the MissensExplorer model score, calibrated with the ClinGen method (Pejaver et
    al., AJHG 2022) on the same ClinVar data as REVEL's ClinGen intervals; see docs/model.md and
    benchmarks/calibrate_acmg.py. The calibration is this project's, not one published by ClinGen.
  * "revel": the intervals ClinGen published: REVEL (Pejaver et al. 2022, as reprinted in Bergquist et
    al. 2025, Table 1); AlphaMissense when REVEL has no score (Bergquist et al., Genet Med 2025).
Points follow the ACMG/AMP point system: Supporting 1, Moderate 2, (3), Strong 4; negative for BP4.
The 3-point level, between Moderate and Strong, is written "Moderate+" (PP3_Moderate+ +3) so that it is
not mistaken for Moderate (2 points).

BP4 means no predicted impact on the protein AND on splicing (Richards et al. 2015). It is therefore
withheld when SpliceAI does not rule out a splicing effect: the ClinGen splicing subgroup takes a
SpliceAI score of 0.1 or less as evidence of no splicing impact (Walker et al., AJHG 2023), so BP4 is
not given above 0.1. PP3 is kept: a damaging protein prediction and a possible splicing effect point
the same way, and PP3 is applied once.

PP3 and the residue-level criteria (PS1/PM5 here; PM1, which a curator adds) can share information,
such as the conservation of the residue; Pejaver et al. (2022) warn that combining a strong PP3 with
them risks double counting, and recommend that PP3 plus PM1 not exceed 4 points. The tool does not cap
the sum (it does not add the points up), but notes the combination (:func:`overlap_note`).

PM2 (absent from controls) is given at Supporting strength (ClinGen SVI, 2020) only when gnomAD itself
confirms the variant is absent and at least half of its exome samples have >= 20 reads at the
position (gnomad.py). Absence meets any disease's maximum credible allele frequency; a variant that
is present but rare may or may not, so then PM2 is left to the user (ClinGen, 2025: thresholds are
disease-specific). PS1/PM5 come from residue.py. PM2 and BS1 may be combined with PP3/BP4 without
limits (Pejaver et al. 2022).

Population frequency uses the highest allele frequency among gnomAD's main genetic ancestry groups
(grpmax: African/African American, Admixed American, East Asian, non-Finnish European, South Asian;
bottlenecked groups excluded, as gnomAD does). BA1 (>= 5%) is stand-alone benign evidence; >= 1% is
flagged as a BS1 candidate, which depends on the disease (too common for a rare dominant disorder,
but frequent variants can be pathogenic in recessive or founder conditions). Variants on ClinGen's
BA1 exception list are never called benign for being common.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, TypeAlias, TypedDict, cast

from . import model
from ._types import Row
from .calibration_params import CALIBRATION as _CALIBRATION
from .model_params import MODEL

Thresholds: TypeAlias = Sequence[tuple[float, int]]
"""Score cut-offs with their points, strongest first: ``((cut-off, points), ...)``."""


class Calibration(TypedDict):
    """Shape of ``calibration_params.CALIBRATION``, written by benchmarks/calibrate_acmg.py.

    Attributes:
        model_version: Version of the model the thresholds were calibrated for.
        model_trained: Training date of that model.
        calibrated: Date of the calibration.
        method: How the thresholds were derived.
        pp3: Model-score cut-offs for PP3 (score >= cut-off), strongest first.
        bp4: Model-score cut-offs for BP4 (score <= cut-off; negative points), strongest first.
    """

    model_version: str
    model_trained: str
    calibrated: str
    method: str
    pp3: Thresholds
    bp4: Thresholds


CALIBRATION = cast(Calibration, _CALIBRATION)

SOURCES = ("model", "revel")
# The thresholds belong to one trained model: after retraining they apply only once recalibrated.
MODEL_CALIBRATED = bool(CALIBRATION["pp3"]) and (CALIBRATION["model_version"], CALIBRATION["model_trained"]) == (
    MODEL["version"],
    MODEL["trained"],
)

# REVEL: Pejaver et al. 2022, as reprinted in Bergquist et al. 2025, Table 1.
REVEL_PP3 = ((0.932, 4), (0.879, 3), (0.773, 2), (0.644, 1))
REVEL_BP4 = ((0.016, -4), (0.052, -3), (0.183, -2), (0.290, -1))
# AlphaMissense: Bergquist et al., Genet Med 2025.
AM_PP3 = ((0.990, 4), (0.972, 3), (0.906, 2), (0.792, 1))
AM_BP4 = ((0.070, -3), (0.099, -2), (0.169, -1))
STRENGTH = {1: "Supporting", 2: "Moderate", 3: "Moderate+", 4: "Strong"}
BP4_SPLICEAI_MAX = 0.1  # SpliceAI above this does not rule out a splicing effect (Walker et al. 2023)

BA1_AF = 0.05
BS1_AF = 0.01
GRPMAX_GROUPS = ("afr", "amr", "eas", "nfe", "sas")

# ClinGen SVI BA1 exception list (Ghosh et al., Hum Mutat 2018, Table 1): common variants classified as
# pathogenic or uncertain. Keyed by dbSNP id and forward-strand alternate allele (None = any allele), so
# the match does not depend on the transcript (BTD p.Asp444His is p.Asp424His on the MANE transcript).
BA1_EXCEPTIONS = {
    ("rs387906242", None): "ACAD9 NM_014049.4:c.-44_-41dupTAAG",
    ("rs72474224", "T"): "GJB2 NM_004004.5:c.109G>A (p.Val37Ile)",
    ("rs1799945", "G"): "HFE NM_000410.3:c.187C>G (p.His63Asp)",
    ("rs1800562", "A"): "HFE NM_000410.3:c.845G>A (p.Cys282Tyr)",
    ("rs11466023", "A"): "MEFV NM_000243.2:c.1105C>T (p.Pro369Ser)",
    ("rs11466024", "T"): "MEFV NM_000243.2:c.1223G>A (p.Arg408Gln)",
    ("rs17089782", "A"): "PIBF1 NM_006346.2:c.1214G>A (p.Arg405Gln)",
    ("rs1800556", "T"): "ACADS NM_000017.3:c.511C>T (p.Arg171Trp)",
    ("rs13078881", "C"): "BTD NM_000060.4:c.1330G>C (p.Asp444His)",
}

# Points of the PS1/PM5 (residue.py) and PM2 codes, which are stored without their points.
CODE_POINTS = {"PS1": 4, "PS1_Moderate": 2, "PM5": 2, "PM5_Supporting": 1, "PM2_Supporting": 1}


def _num(x: Any) -> float | None:
    """Return ``x`` as a float, or None when it is missing, not a number or NaN."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _points(s: float, pp3: Thresholds, bp4: Thresholds) -> int:
    """Return the ACMG/AMP points of score ``s``: positive for PP3, negative for BP4, 0 in between."""
    for lo, p in pp3:
        if s >= lo:
            return p
    for hi, p in bp4:
        if s <= hi:
            return p
    return 0


def _evidence(s: float, pp3: Thresholds, bp4: Thresholds, shown: str) -> tuple[str, int]:
    """Return the PP3/BP4 text and points of score ``s``; ``shown`` names the score in the text."""
    pts = _points(s, pp3, bp4)
    if pts == 0:
        return f"none ({shown}, indeterminate)", 0
    return f"{'PP3_' if pts > 0 else 'BP4_'}{STRENGTH[abs(pts)]} ({shown})", pts


def computational(row: Row, source: str = "model") -> tuple[str | None, int | None]:
    """Return the PP3/BP4 evidence of a row and its points.

    Examples: ('PP3_Strong (model 0.995)', 4), ('none (REVEL 0.500, indeterminate)', 0) or (None, None).

    Args:
        row: One results row; reads the predictor columns (``revel``, ``alphamissense``, ...).
        source: ``"model"`` (the model score with its own calibration, when the shipped model is the
            calibrated one; otherwise as ``"revel"``) or ``"revel"`` (ClinGen's REVEL intervals,
            AlphaMissense when REVEL has no score).

    Returns:
        The evidence text and its points (negative for BP4, 0 when indeterminate), or (None, None)
        when the row has no usable score.

    Raises:
        ValueError: If ``source`` is not one of :data:`SOURCES`.
    """
    if source not in SOURCES:
        raise ValueError(f"unknown PP3/BP4 source {source!r}; use one of {SOURCES}")
    return _splicing_check(row, *_computational(row, source))


def _computational(row: Row, source: str) -> tuple[str | None, int | None]:
    """PP3/BP4 from the protein-level score alone (see :func:`computational`)."""
    if source == "model" and MODEL_CALIBRATED:
        ms = model.score(row)
        if ms is None:
            return None, None
        s = round(ms.p, 3)  # the value shown is the value compared (thresholds lie on a 0.001 grid)
        return _evidence(s, CALIBRATION["pp3"], CALIBRATION["bp4"], f"model {s:.3f}")
    for tool, col, pp3, bp4 in (
        ("REVEL", "revel", REVEL_PP3, REVEL_BP4),
        ("AlphaMissense", "alphamissense", AM_PP3, AM_BP4),
    ):
        value = _num(row.get(col))
        if value is not None:
            return _evidence(value, pp3, bp4, f"{tool} {value:.3f}")
    return None, None


def _splicing_check(row: Row, text: str | None, pts: int | None) -> tuple[str | None, int | None]:
    """Withhold BP4 unless SpliceAI rules out a splicing effect (score at most 0.1; Walker et al. 2023)."""
    sp = _num(row.get("spliceai_max"))
    if pts is not None and pts < 0 and text and (sp is None or sp > BP4_SPLICEAI_MAX):
        code, shown = text.split(" ", 1)[0], text.split("(", 1)[1].rstrip(")")
        why = "no SpliceAI score to rule out" if sp is None else f"SpliceAI {sp:.2f} does not rule out"
        return f"none ({shown}; {code} withheld: {why} a splicing effect)", 0
    return text, pts


def overlap_note(row: Row) -> str | None:
    """Note when a moderate or stronger PP3 meets PS1/PM5, which may share information (Pejaver 2022).

    Args:
        row: One results row with ``acmg_points`` and ``acmg_ps1_pm5`` filled.

    Returns:
        The note, or None.
    """
    pts = _num(row.get("acmg_points"))
    ps = str(row.get("acmg_ps1_pm5") or "")
    if pts is not None and pts >= 2 and ps[:3] in ("PS1", "PM5"):
        return (
            f"PP3 ({pts} points) and {ps.split(' ')[0]} may partly reflect the same residue-level information "
            "(Pejaver et al. 2022): weigh their combination with care; if PM1 is added, PP3 + PM1 should not "
            "exceed 4 points"
        )
    return None


def thresholds_used(evidence: object) -> str | None:
    """Return which thresholds produced a row's PP3/BP4 evidence, for the note shown under it.

    Args:
        evidence: The row's ``acmg_pp3_bp4``, as written by :func:`computational`.

    Returns:
        ``"model"`` (the model score with this project's calibration), ``"clingen"`` (REVEL or
        AlphaMissense at ClinGen's published thresholds), or None when the row has no PP3/BP4, as for
        variants that are not missense.
    """
    text = evidence if isinstance(evidence, str) else ""
    if not text:
        return None
    return "model" if "(model " in text else "clingen"


def ba1_grpmax(row: Row) -> float | None:
    """Return the grpmax frequency that BA1 and the common-variant rule use.

    Ensembl's grpmax, unless gnomAD itself was asked about a BA1 candidate (``_grpmax_2000``, set by
    predict): then the highest frequency among the groups with at least 2,000 alleles, exomes and
    genomes combined, as BA1 requires (Ghosh et al. 2018).
    """
    af = _num(row.get("gnomad_popmax_af"))
    g2 = row.get("_grpmax_2000")
    if af is not None and af > BA1_AF and isinstance(g2, tuple) and g2[0] is not None:
        return float(g2[0])
    return af


def ba1_exception(row: Row) -> str | None:
    """Return the ClinGen BA1 exception this variant matches, if any.

    Args:
        row: One results row; reads ``rsid`` (comma-separated dbSNP ids) and ``alt`` (forward strand).

    Returns:
        The exception's description (gene, transcript change), or None.
    """
    rsids = [x.strip() for x in str(row.get("rsid") or "").split(",") if x.strip()]
    alt = str(row.get("alt") or "").upper()
    for rs in rsids:
        hit = BA1_EXCEPTIONS.get((rs, alt)) or BA1_EXCEPTIONS.get((rs, None))
        if hit:
            return hit
    return None


def pm2(row: Row, info: dict[str, Any] | None) -> tuple[str | None, int | None]:
    """Return PM2 from the direct gnomAD check (`info`, from gnomad.lookup) or Ensembl's frequency.

    Args:
        row: One results row; reads ``gnomad_popmax_af`` and ``gnomad_af``.
        info: The variant's entry from :func:`missensexplorer.gnomad.lookup`, or None when gnomAD was
            not asked.

    Returns:
        The PM2 text and its points (1 for PM2_Supporting, otherwise 0), or (None, None) when no
        frequency information is available.
    """
    from .gnomad import covered

    af = _num(row.get("gnomad_popmax_af")) or _num(row.get("gnomad_af"))
    if af is not None and af >= BS1_AF:
        return "not met (frequent in gnomAD)", 0
    if info and not info.get("found"):
        cov = (info.get("over_20") or {}).get("exome")
        shown = "no exome data" if cov is None else f"{cov:.0%} of exomes"
        if covered(info):
            return f"PM2_Supporting (absent from {info['dataset']}; >= 20 reads in {shown})", 1
        return f"not applicable (absent from {info['dataset']}, but few reads there: >= 20 reads in {shown})", 0
    if af is not None or (info and info.get("found")):
        filt = "" if not info or info.get("passes_filters", True) else ", failing gnomAD quality filters"
        return f"depends on the disease (present in gnomAD{filt}; compare with its maximum credible frequency)", 0
    return None, None


def codes(row: Row) -> str | None:
    """Return the criteria this row meets, with points: 'PP3_Strong +4, PM5 +2, PM2_Supporting +1'.

    Args:
        row: One results row with the ``acmg_pp3_bp4``, ``acmg_points``, ``acmg_ps1_pm5``, ``acmg_pm2``
            and ``population_frequency`` columns already filled.

    Returns:
        The comma-separated codes (BA1 is added as stand-alone benign evidence), or None when no
        criterion is met.
    """
    out = []
    for col, pts_col in (("acmg_pp3_bp4", "acmg_points"), ("acmg_ps1_pm5", None), ("acmg_pm2", None)):
        text = str(row.get(col) or "")
        code = text.split(" ")[0]
        if code[:3] in ("PP3", "BP4", "PS1", "PM5", "PM2"):
            pts = row.get(pts_col) if pts_col else CODE_POINTS.get(code)
            out.append(f"{code} {int(pts):+d}" if pts is not None else code)
    freq = str(row.get("population_frequency") or "")
    if freq.startswith("common: BA1"):
        out.append("BA1 (stand-alone benign)")
    return ", ".join(out) or None


def _pct(af: float) -> str:
    """6.2%, 0.0012%, 0.00009% (two significant digits, never in scientific notation)."""
    v = af * 100
    return f"{v:.2g}%" if v >= 1e-3 else f"{float(f'{v:.2g}'):.10f}".rstrip("0").rstrip(".") + "%"


def frequency(row: Row) -> str | None:
    """Return the population-frequency evidence from the gnomAD grpmax allele frequency.

    Args:
        row: One results row; reads ``chrom``, ``gnomad_popmax_af``, ``gnomad_af``, ``rsid``, ``alt``
            and ``_gnomad`` (the direct gnomAD check, if one was made).

    Returns:
        A text such as 'common: BA1 (grpmax 6.2%)' or 'rare (grpmax 0.002%)', or None for mitochondrial
        variants and for variants with an overall but no grpmax frequency.
    """
    if str(row.get("chrom") or "").upper() in ("MT", "M", "CHRM"):
        return None  # mitochondrial heteroplasmy is not in the gnomAD nuclear frequencies
    af = _num(row.get("gnomad_popmax_af"))
    if af is None:
        if _num(row.get("gnomad_af")) is not None:
            return None
        g = row.get("_gnomad")
        if g and not g.get("found"):
            return f"absent from {g['dataset']}"
        return "not found in gnomAD" if g else "no gnomAD frequency in Ensembl's data (not checked in gnomAD)"
    pct = _pct(af)
    if af >= BS1_AF and ba1_exception(row):
        return f"frequent, but on the ClinGen BA1 exception list (grpmax {pct})"
    if af > BA1_AF and (b := ba1_grpmax(row)) is not None and b <= BA1_AF:
        return (
            f"frequent: BS1 candidate (grpmax {pct}); BA1 not met: {_pct(b)} at most in the gnomAD groups "
            "with >= 2,000 alleles"
        )
    if af > BA1_AF:
        checked = isinstance(row.get("_grpmax_2000"), tuple)
        return f"common: BA1 (grpmax {pct}{'' if checked else '; allele numbers not checked in gnomAD'})"
    if af >= BS1_AF:
        return f"frequent: BS1 candidate (grpmax {pct})"
    return f"rare (grpmax {pct})"
