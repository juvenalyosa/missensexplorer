"""ClinGen expert panel specifications (vcep.py, vcep_rules.py): what differs from the tool's evidence."""

import pytest

from missensexplorer import vcep, vcep_rules
from missensexplorer.gui import i18n


def test_a_gene_with_a_released_specification_is_named_with_its_disease_and_link():
    text = vcep.describe("MYH7")
    assert text.startswith("Cardiomyopathy VCEP v")
    assert "https://cspec.genome.network/" in text
    assert vcep.describe("NOTAGENE") is None


def test_several_specifications_of_a_gene_are_all_kept():
    ids = [s["id"] for s in vcep.specifications("RYR1")]
    assert ids == ["GN012", "GN150", "GN179"]  # malignant hyperthermia; dominant and recessive myopathy
    assert [s["id"] for s in vcep.specifications("VWF")] == ["GN081", "GN090"]
    notes = vcep.notes({"gene": "RYR1", "acmg_points": 4})
    assert len(notes) == 3
    assert any("autosomal recessive; GN179)" in n for n in notes)
    assert "x-linked" in vcep.describe("ABCD1")  # mode of inheritance from CSpec's structured field


def test_evidence_stronger_than_the_panel_allows_is_noted():
    notes = vcep.notes({"gene": "MYH7", "acmg_points": 4})  # the cardiomyopathy panel allows PP3 Supporting only
    assert len(notes) == 1
    assert "allows PP3 up to Supporting: the tool's PP3_Strong (4 points) is stronger" in notes[0]
    assert "sets its own predictor and thresholds for PP3/BP4" in notes[0]
    supporting = vcep.notes({"gene": "MYH7", "acmg_points": 1})[0]
    assert "allows PP3" not in supporting
    assert "own predictor" in supporting  # never silent on PP3/BP4
    assert vcep.notes({"gene": "NOTAGENE", "acmg_points": 4}) == []


def test_pp3_or_bp4_that_a_panel_does_not_use_for_missense_variants():
    note = vcep.notes({"gene": "ADA", "acmg_points": 1})[0]  # SCID panel: "Do not apply to missense variants"
    assert "does not use PP3 for missense variants" in note
    assert "allows PP3 up to" not in note
    assert "does not use BP4 for missense variants" in vcep.notes({"gene": "VHL", "acmg_points": -1})[0]
    excluded = {(g, c) for g, specs in vcep.VCEP["genes"].items() for s in specs for c in s["not_missense"]}
    assert {"ADA", "DCLRE1C", "IL7R", "JAK3", "RAG1", "RAG2", "IL2RG", "PALB2", "APC", "CDH1"} <= {
        g for g, c in excluded if c == "PP3"
    }
    assert {"PALB2", "VHL", "APC", "DDC", "CDH1", "PIK3CA", "MTOR"} <= {g for g, c in excluded if c == "BP4"}
    # splicing-only (CDH1) or non-coding-only (Brain Malformations) BP4: never "allows BP4 up to Supporting"
    for gene in ("CDH1", "PIK3CA"):
        note = vcep.notes({"gene": gene, "acmg_points": -2})[0]
        assert "does not use BP4 for missense variants" in note
        assert "allows BP4" not in note


def test_weaker_or_absent_strengths_and_the_three_point_level():
    ldlr = vcep.notes({"gene": "LDLR", "acmg_points": 0, "acmg_pm2": "PM2_Supporting (absent from gnomAD v4)"})[0]
    assert "uses PM2 only at Moderate: the tool's PM2_Supporting is weaker" in ldlr
    pah = vcep.notes({"gene": "PAH", "acmg_points": 3, "revel": 0.9})[0]  # the PKU panel allows 1, 2 and 4 points
    assert "has no 3-point PP3" in pah


def test_ba1_against_the_panel_threshold():
    common = {"gene": "VWF", "acmg_points": 0, "gnomad_popmax_af": 0.07, "population_frequency": "common: BA1 (7%)"}
    assert any("BA1 only > 10%: grpmax 7% does not meet it" in n for n in vcep.notes(common))
    rare = {"gene": "MYH7", "acmg_points": 0, "gnomad_popmax_af": 0.004, "population_frequency": "rare"}
    note = vcep.notes(rare)[0]  # the cardiomyopathy panel applies BA1 to the filtering allele frequency
    assert "is above its BA1 threshold (>= 0.1%), which it applies to the filtering allele frequency" in note
    assert "meets" not in note
    kcnq1 = {"gene": "KCNQ1", "acmg_points": 0, "gnomad_popmax_af": 0.02, "population_frequency": "frequent"}
    assert "meets its BA1 threshold (>= 0.4%)" in vcep.notes(kcnq1)[0]
    assert vcep.ba1_above_5pct("VWF") == pytest.approx(0.1)
    assert vcep.ba1_above_5pct("MYH7") is None


def test_encoded_rules_say_what_the_panel_would_give():
    note = vcep.notes({"gene": "LDLR", "acmg_points": 4, "revel": 0.81, "spliceai_max": 0.0,
                       "gnomad_popmax_af": None, "acmg_pm2": "PM2_Supporting (absent from gnomAD v4)"})[0]  # fmt: skip
    assert "its rules give PP3_Supporting (REVEL 0.810), PM2 met" in note
    ada = vcep_rules.describe({"gene": "ADA", "gnomad_popmax_af": 0.01})
    assert ada == "no PP3/BP4 for missense variants, PM2 not met, BA1 met"


def test_weaker_residue_evidence_is_not_an_invitation_to_upgrade():
    note = vcep.notes({"gene": "LDLR", "acmg_points": 0, "acmg_ps1_pm5": "PM5_Supporting (p.X: Likely pathogenic)"})[0]
    assert "does not use PM5 below Moderate: with a likely pathogenic or single-star reference" in note
    assert "is weaker" not in note.split("PM2")[0]
    assert "classified pathogenic (by these guidelines)" in note  # the panel's condition on the reference


def test_conditions_the_tool_cannot_check_are_quoted():
    note = vcep.notes({"gene": "DICER1", "acmg_points": 0, "acmg_ps1_pm5": "PS1_Moderate (same change: LP)"})[0]
    assert "on PS1: “Likely pathogenic changes do not apply" in note


@pytest.mark.parametrize("gene", ["MYH7", "RYR1"])
def test_vcep_notes_are_translated_part_by_part(gene):
    i18n.set_language("es")
    try:
        out = i18n.tr_text(vcep.SEP.join(vcep.notes({"gene": gene, "acmg_points": 4})))
    finally:
        i18n.set_language("en")
    assert "admite PP3 hasta Supporting" in out or "admite PP3 hasta Moderate" in out
    assert "allows" not in out  # every note of the cell, every part
    assert "own predictor" not in out
