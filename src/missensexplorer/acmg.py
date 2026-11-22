"""ACMG/AMP evidence that can be read off the annotation: PP3/BP4 and population frequency.

PP3/BP4 (computational evidence) uses ONE predictor with ClinGen-calibrated score intervals, as
ClinGen recommends (no combining of tools):
  * REVEL — Pejaver et al., AJHG 2022 (intervals as reprinted in Bergquist et al. 2025, Table 1);
  * AlphaMissense, when REVEL has no score — Bergquist et al., Genet Med 2025, Table 1.
Points follow the ACMG/AMP point system: Supporting 1, Moderate 2, (3), Strong 4; negative for BP4.
The 3-point intervals count as Moderate under the 2015 guidelines.

Population frequency uses the highest allele frequency among gnomAD's main genetic ancestry groups
(grpmax: African/African American, Admixed American, East Asian, non-Finnish European, South Asian;
bottlenecked groups excluded, as gnomAD does). BA1 (>= 5%) is stand-alone benign evidence; >= 1% is
flagged as a BS1 candidate, which depends on the disease (too common for a rare dominant disorder,
but frequent variants can be pathogenic in recessive or founder conditions). Variants on ClinGen's
BA1 exception list are never called benign for being common.
"""

from __future__ import annotations

REVEL_PP3 = ((0.932, 4), (0.879, 3), (0.773, 2), (0.644, 1))
REVEL_BP4 = ((0.016, -4), (0.052, -3), (0.183, -2), (0.290, -1))
AM_PP3 = ((0.990, 4), (0.972, 3), (0.906, 2), (0.792, 1))
AM_BP4 = ((0.070, -3), (0.099, -2), (0.169, -1))
STRENGTH = {1: "Supporting", 2: "Moderate", 3: "Moderate", 4: "Strong"}

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


def _num(x) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _points(s: float, pp3, bp4) -> int:
    for lo, p in pp3:
        if s >= lo:
            return p
    for hi, p in bp4:
        if s <= hi:
            return p
    return 0


def computational(row: dict) -> tuple[str | None, int | None]:
    """('PP3_Strong (REVEL 0.95)', 4), ('none (REVEL 0.50, indeterminate)', 0) or (None, None)."""
    for tool, col, pp3, bp4, fmt in (("REVEL", "revel", REVEL_PP3, REVEL_BP4, "{:.3f}"),
                                     ("AlphaMissense", "alphamissense", AM_PP3, AM_BP4, "{:.3f}")):
        s = _num(row.get(col))
        if s is None:
            continue
        pts = _points(s, pp3, bp4)
        shown = f"{tool} {fmt.format(s)}"
        if pts == 0:
            return f"none ({shown}, indeterminate)", 0
        code = ("PP3_" if pts > 0 else "BP4_") + STRENGTH[abs(pts)]
        return f"{code} ({shown})", pts
    return None, None


def ba1_exception(row: dict) -> str | None:
    """The ClinGen BA1 exception this variant matches, if any."""
    rsids = [x.strip() for x in str(row.get("rsid") or "").split(",") if x.strip()]
    alt = str(row.get("alt") or "").upper()
    for rs in rsids:
        hit = BA1_EXCEPTIONS.get((rs, alt)) or BA1_EXCEPTIONS.get((rs, None))
        if hit:
            return hit
    return None


def frequency(row: dict) -> str | None:
    """Population-frequency evidence from the gnomAD grpmax allele frequency."""
    if str(row.get("chrom") or "").upper() in ("MT", "M", "CHRM"):
        return None  # mitochondrial heteroplasmy is not in the gnomAD nuclear frequencies
    af = _num(row.get("gnomad_popmax_af"))
    if af is None:
        return "not found in gnomAD" if _num(row.get("gnomad_af")) is None else None
    pct = f"{af * 100:.2g}%"
    if af >= BS1_AF and ba1_exception(row):
        return f"frequent, but on the ClinGen BA1 exception list (grpmax {pct})"
    if af >= BA1_AF:
        return f"common: BA1 (grpmax {pct})"
    if af >= BS1_AF:
        return f"frequent: BS1 candidate (grpmax {pct})"
    return f"rare (grpmax {pct})"
