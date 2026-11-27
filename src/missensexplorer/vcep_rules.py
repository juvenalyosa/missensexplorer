"""The PP3/BP4, PM2 and BA1 rules of ten ClinGen expert panels, encoded by hand from their specifications.

ClinGen's Criteria Specification Registry (CSpec) gives each panel's thresholds as free text. For the ten
panels with most classifications in ClinGen's Evidence Repository (55% of the missense classifications),
the PP3/BP4, PM2 and BA1 rules were read in CSpec (November 2026) and written down here, each next to the
text it comes from. A second encoding made independently from the same texts, during the review of the
MissensExplorer article, gave the same codes for all 3,460 of these panels' classifications.

The app reports what these rules give for a variant (vcep.notes); benchmarks/erepo_vcep_rules.py compares
them with the panels' decisions. Approximations: frequency thresholds that a panel applies to the filtering
allele frequency (FAF) are applied to gnomAD's grpmax allele frequency, which is at least the FAF (so PM2
is given less often, and BA1 more often, than with the FAF); homozygote counts, coverage requirements,
splicing tools other than SpliceAI and exception lists are not applied; the rules are those of the
current specification versions.
"""

# ruff: noqa: E501  (the hand-encoded rules keep one gene per line, next to the text they come from)
from __future__ import annotations

from typing import Any, cast

import numpy as np

# PP3/BP4: list of (threshold, points) on a score, strongest first; ">=" for PP3, "<=" for BP4 unless noted.
REVEL_PEJAVER_PP3 = [(0.932, 4), (0.773, 2), (0.644, 1)]  # as specified (Pejaver et al.), without the +3 level
REVEL_PEJAVER_BP4 = [(0.016, 4), (0.183, 2), (0.290, 1)]

RULES: dict[str, dict] = {
    # Myeloid Malignancy VCEP v3.1 (GN008): PP3 "REVEL score >= 0.88 or SpliceAI >= 0.38"; BP4 "REVEL score
    # < 0.50 AND SpliceAI <= 0.20"; PM2 "Minor allele frequency <= 0.00005 with at least 2000 alleles tested
    # around and 20x coverage at the position" (an absent variant counts only with that coverage);
    # BA1 ">= 0.0015 (0.15%)".
    "RUNX1": {"panel": "Myeloid Malignancy", "spec": "GN008", "score": "revel", "pp3": [(0.88, 1)],
              "pp3_spliceai": 0.38, "bp4": [(0.50, 1)], "bp4_strict": True, "bp4_spliceai_max": 0.20,
              "pm2": ("le", 0.00005), "pm2_coverage": True, "ba1": ("ge", 0.0015)},
    # Monogenic Diabetes VCEP (GN086 GCK v3.1, GN017 HNF1A v3.1, GN085 HNF4A v4.0): PP3 "REVEL score of
    # >= 0.70"; BP4 "REVEL score of <= 0.15"; PM2 "Grpmax FAF <= 1:333,000 (<= 0.000003)"; BA1 "Grpmax FAF
    # >= 1:10,000 (>= 0.0001)".
    **{g: {"panel": "Monogenic Diabetes", "spec": s, "score": "revel", "pp3": [(0.70, 1)], "bp4": [(0.15, 1)],
           "pm2": ("le", 0.000003), "ba1": ("ge", 0.0001)} for g, s in (("GCK", "GN086"), ("HNF1A", "GN017"),
                                                                          ("HNF4A", "GN085"))},
    # Familial Hypercholesterolemia VCEP v1.2 (GN013): PP3 "REVEL score >= 0.75"; BP4 "REVEL score <= 0.5";
    # PM2 (moderate) "PopMax MAF <= 0.0002"; BA1 "PopMax FAF >= 0.005".
    "LDLR": {"panel": "Familial Hypercholesterolemia", "spec": "GN013", "score": "revel", "pp3": [(0.75, 1)],
             "bp4": [(0.5, 1)], "pm2": ("le", 0.0002), "ba1": ("ge", 0.005)},
    # Phenylketonuria VCEP v2.1 (GN006): PP3 and BP4 "as described in Pejaver et al." (REVEL, supporting to
    # strong); PM2 "< 0.0002"; BA1 ">= 0.015".
    "PAH": {"panel": "Phenylketonuria", "spec": "GN006", "score": "revel", "pp3": REVEL_PEJAVER_PP3,
            "bp4": REVEL_PEJAVER_BP4, "pm2": ("lt", 0.0002), "ba1": ("ge", 0.015)},
    # Rett and Angelman-like Disorders VCEP (GN036 MECP2, GN034 CDKL5, GN035 FOXG1, GN032 TCF4 v6.0, GN037
    # UBE3A v7.0, GN033 SLC9A6 v6.0): PP3 "REVEL with a score >= 0.644" (SLC9A6: 0.664); BP4 "REVEL with a
    # score <= 0.290"; PM2 "absent, zero observations"; BA1 ">= 0.000083 in any sub-population".
    **{g: {"panel": "Rett and Angelman-like Disorders", "spec": s, "score": "revel",
           "pp3": [(0.664 if g == "SLC9A6" else 0.644, 1)], "bp4": [(0.290, 1)], "pm2": ("absent", None),
           "ba1": ("ge", 0.000083)}
       for g, s in (("MECP2", "GN036"), ("CDKL5", "GN034"), ("FOXG1", "GN035"), ("TCF4", "GN032"),
                    ("UBE3A", "GN037"), ("SLC9A6", "GN033"))},
    # Severe Combined Immunodeficiency Disease VCEP (v2.2-2.3): PP3 "Do not apply to missense variants"
    # (ADA, DCLRE1C, IL7R, JAK3, RAG1, RAG2, IL2RG); FOXN1 (GN113) PP3 "REVEL >= 0.644" supporting, ">= 0.932"
    # moderate, BP4 "REVEL < 0.290". PM2 and BA1: popmax FAF thresholds per gene.
    **{g: {"panel": "Severe Combined Immunodeficiency", "spec": s, "score": None, "pm2": pm2, "ba1": ("gt", ba1)}
       for g, s, pm2, ba1 in (("ADA", "GN114", ("lt", 0.0001742), 0.00721), ("DCLRE1C", "GN116", ("lt", 0.00003266), 0.00346),
                              ("IL7R", "GN119", ("lt", 0.00004129), 0.00566), ("JAK3", "GN121", ("lt", 0.000115), 0.00447),
                              ("RAG1", "GN123", ("lt", 0.000102), 0.00872), ("RAG2", "GN124", ("lt", 0.0000588), 0.00872),
                              ("IL2RG", "GN129", None, 0.01110))},
    "FOXN1": {"panel": "Severe Combined Immunodeficiency", "spec": "GN113", "score": "revel",
              "pp3": [(0.932, 2), (0.644, 1)], "bp4": [(0.290, 1)], "bp4_strict": True, "pm2": ("le", 0.00002412),
              "ba1": ("gt", 0.00447)},
    # Malignant Hyperthermia Susceptibility VCEP v2.0 (GN012): PP3 moderate "REVEL score of > 0.85"; BP4
    # "REVEL score of < 0.5"; PM2 not specified; BA1 "Popmax allele frequency > 0.0038".
    "RYR1": {"panel": "Malignant Hyperthermia Susceptibility", "spec": "GN012", "score": "revel",
             "pp3": [(0.85, 2)], "pp3_strict": True, "bp4": [(0.5, 1)], "bp4_strict": True, "pm2": None,
             "ba1": ("gt", 0.0038)},
    # RASopathy VCEP (v1.3-2.3): PP3 "REVEL >= 0.7"; BP4 "REVEL <= 0.3"; PM2 "absent from controls"; BA1
    # "filtering allele frequency >= 0.05%".
    **{g: {"panel": "RASopathy", "spec": "RASopathy", "score": "revel", "pp3": [(0.7, 1)], "bp4": [(0.3, 1)],
           "pm2": ("absent", None), "ba1": ("ge", 0.0005)}
       for g in ("SOS1", "PTPN11", "BRAF", "RAF1", "MAP2K2", "MAP2K1", "LZTR1", "NRAS", "SHOC2", "KRAS", "HRAS",
                 "SOS2", "PPP1CB", "MRAS", "RIT1", "RRAS2")},
    # Lysosomal Diseases VCEP: GAA v2.0 (GN010) PP3 "REVEL score > 0.7", BP4 "REVEL score < 0.5", PM2 "< 0.1%
    # (0.001)", BA1 "> 0.01"; IDUA v1.2 (GN091) PP3 "> 0.773" moderate, "0.644 - 0.773" supporting, BP4
    # "< 0.290", PM2 "< 0.025% (0.00025)", BA1 "Grpmax > 0.005".
    "GAA": {"panel": "Lysosomal Diseases", "spec": "GN010", "score": "revel", "pp3": [(0.7, 1)], "pp3_strict": True,
            "bp4": [(0.5, 1)], "bp4_strict": True, "pm2": ("lt", 0.001), "ba1": ("gt", 0.01)},
    "IDUA": {"panel": "Lysosomal Diseases", "spec": "GN091", "score": "revel", "pp3": [(0.773, 2, True), (0.644, 1, False)],
             "bp4": [(0.290, 1)], "bp4_strict": True, "pm2": ("lt", 0.00025), "ba1": ("gt", 0.005)},
    # Limb Girdle Muscular Dystrophy VCEP v2.0 (GN180 DYSF, GN187 CAPN3, GN188 ANO5, GN184-186/189 SGCB/G/D/A):
    # PP3 "REVEL with a score >= 0.7"; BP4 "REVEL with a score <= 0.1 AND SpliceAI with a score <= 0.05"; PM2
    # "< 0.0001" (sarcoglycans "< 0.00009"); BA1 "Grpmax FAF ... > 0.003" (sarcoglycans "> 0.002").
    **{g: {"panel": "Limb Girdle Muscular Dystrophy", "spec": s, "score": "revel", "pp3": [(0.7, 1)],
           "bp4": [(0.1, 1)], "bp4_spliceai_max": 0.05,
           "pm2": ("lt", 0.00009 if g.startswith("SGC") else 0.0001), "ba1": ("gt", 0.002 if g.startswith("SGC") else 0.003)}
       for g, s in (("DYSF", "GN180"), ("CAPN3", "GN187"), ("ANO5", "GN188"), ("SGCA", "GN189"), ("SGCB", "GN184"),
                    ("SGCG", "GN185"), ("SGCD", "GN186"))},
}  # fmt: skip


def _cmp(x: float, op: str, t: float) -> bool:
    return {"ge": x >= t, "gt": x > t, "le": x <= t, "lt": x < t}[op]


def panel_rule_codes(
    row: dict[str, Any], rule: dict[str, Any], af_key: str = "gnomad_popmax_af", pm2_key: str = "acmg_pm2"
) -> dict[str, int]:
    """PP3/BP4 points, PM2 and BA1 (1/0) for one variant under its panel's rule.

    Reads ``revel``, ``spliceai_max``, the gnomAD grpmax allele frequency (``af_key``) and the tool's PM2
    text (``pm2_key``: "absent from gnomAD" in it means absent, whatever the coverage; the rules count it).
    """
    out = {"PP3": 0, "BP4": 0, "PM2": 0, "BA1": 0}
    s = row.get(rule["score"]) if rule.get("score") else None
    sp = row.get("spliceai_max")
    sp = None if sp is None or (isinstance(sp, float) and np.isnan(sp)) else float(sp)
    if s is not None and not (isinstance(s, float) and np.isnan(s)):
        for t, p, *strict in rule.get("pp3", []):  # (threshold, points[, strict]): strict means "above"
            if (s > t) if (strict[0] if strict else rule.get("pp3_strict")) else (s >= t):
                out["PP3"] = p
                break
        if not out["PP3"]:
            for t, p in rule.get("bp4", []):
                ok = (s < t) if rule.get("bp4_strict") else (s <= t)
                if ok and (rule.get("bp4_spliceai_max") is None or (sp is not None and sp <= rule["bp4_spliceai_max"])):
                    out["BP4"] = p
                    break
    if not out["PP3"] and rule.get("pp3_spliceai") is not None and sp is not None and sp >= rule["pp3_spliceai"]:
        out["PP3"] = 1
    af = row.get(af_key)
    af = None if af is None or (isinstance(af, float) and np.isnan(af)) else float(af)
    # absent from gnomAD counts as absent whatever the coverage, as the panels' rules ask, except where the rule
    # itself requires coverage (Myeloid Malignancy: 20x), where the tool's confirmed absence (PM2_Supporting) is used
    pm2_text = str(row.get(pm2_key) or "")
    absent = af is None and (
        pm2_text.startswith("PM2_Supporting") if rule.get("pm2_coverage") else "absent from gnomAD" in pm2_text
    )
    pm2 = rule.get("pm2")
    if pm2:
        out["PM2"] = int(absent if pm2[0] == "absent" else (absent or (af is not None and _cmp(af, *pm2))))
    ba1 = rule.get("ba1")
    if ba1 and af is not None and _cmp(af, *ba1):
        out["BA1"] = 1
    return out


NAMES = {1: "Supporting", 2: "Moderate", 4: "Strong"}


def describe(row: dict[str, Any]) -> str | None:
    """What the panel's encoded rules give for a row ('PP3_Supporting (REVEL 0.812), PM2 not met'), or None."""
    rule = RULES.get(str(row.get("gene") or ""))
    if rule is None:
        return None
    c = panel_rule_codes(row, rule)
    s = row.get(rule["score"]) if rule.get("score") else None
    parts = []
    if rule.get("score") is None:
        parts.append("no PP3/BP4 for missense variants")
    elif s is None:
        parts.append("PP3/BP4 not evaluated (no REVEL score)")
    elif c["PP3"]:
        parts.append(f"PP3_{NAMES[c['PP3']]} (REVEL {float(s):.3f})")
    elif c["BP4"]:
        parts.append(f"BP4_{NAMES[c['BP4']]} (REVEL {float(s):.3f})")
    else:
        parts.append(f"neither PP3 nor BP4 (REVEL {float(s):.3f})")
    if rule.get("pm2"):
        pm2_text = str(row.get("acmg_pm2") or "")
        af = row.get("gnomad_popmax_af")
        no_data = (af is None or (isinstance(af, float) and np.isnan(af))) and "absent from gnomAD" not in pm2_text
        low = "few reads there" in pm2_text
        if c["PM2"]:
            parts.append("PM2 met (absent from gnomAD, but few reads there: check)" if low else "PM2 met")
        elif no_data:
            parts.append("PM2 not assessed (no gnomAD data)")
        elif low and rule.get("pm2_coverage"):
            parts.append("PM2 not met (absent from gnomAD, but the panel requires 20x coverage)")
        else:
            parts.append("PM2 not met")
    if c["BA1"]:
        # the rule compares the grpmax allele frequency with a threshold the panel applies to the filtering allele
        # frequency, which is at most the grpmax: met here, not necessarily on the panel's own measure
        parts.append(
            "BA1 met on grpmax (the panel uses the filtering allele frequency: check)"
            if _ba1_on_faf(row, rule)
            else "BA1 met"
        )
    return ", ".join(parts)


def _ba1_on_faf(row: dict[str, Any], rule: dict[str, Any]) -> bool:
    """True when the panel's specification applies its BA1 threshold to the filtering allele frequency."""
    from .vcep_params import VCEP

    specs: list[dict[str, Any]] = cast(dict[str, Any], VCEP["genes"]).get(str(row.get("gene") or ""), [])
    mine = [x for x in specs if x["id"] == rule["spec"] or rule["panel"] in x["panel"]]
    return any((x.get("ba1") or {}).get("faf") for x in mine)
