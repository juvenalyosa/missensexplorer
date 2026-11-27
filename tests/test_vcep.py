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
    assert "X-linked" in vcep.describe("ABCD1")  # mode of inheritance from CSpec's structured field


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
    # the panels' further BA1 conditions are quoted whether they use the filtering allele frequency or not
    ryr1 = next(n for n in vcep.notes({"gene": "RYR1", "acmg_points": 0, "gnomad_popmax_af": 0.0004,
                                       "population_frequency": "rare"}) if "GN150" in n)  # fmt: skip
    assert "BA1/BS1 exclusion variants (NM_000540.3: c.14524G>A (p.Val4842Met), c.325C>T (p.Arg109Trp))" in ryr1
    mlh1 = vcep.notes({"gene": "MLH1", "acmg_points": 0, "gnomad_popmax_af": 0.002, "population_frequency": "rare"})[0]
    assert "excluded as founder pathogenic variant" in mlh1
    # a variant on the panel's own exclusion list, or on ClinGen's BA1 exception list, is never sent to BA1
    v4842m = {"gene": "RYR1", "acmg_points": 0, "gnomad_popmax_af": 0.004, "population_frequency": "rare",
              "mane_select": "NM_000540.3", "hgvsc": "ENST00000359596.8:c.14524G>A",
              "hgvsp": "ENSP00000352608.2:p.Val4842Met"}  # fmt: skip
    assert any("on the panel's BA1 exclusion list (well-known pathogenic variants): BA1 and BS1 do not apply" in n
               for n in vcep.notes(v4842m))  # fmt: skip
    # the lists are numbered on MANE Select: the same text on another transcript is another variant
    assert not any("exclusion list" in n for n in vcep.notes({**v4842m, "mane_select": None}))
    # NEB's list is numbered on NM_001271208.2: its p.Thr7417Pro is p.Thr7382Pro (c.22144A>C) on MANE Select
    neb = {
        "gene": "NEB",
        "acmg_points": 0,
        "gnomad_popmax_af": 0.01,
        "population_frequency": "frequent",
        "mane_select": "NM_001164508.2",
        "hgvsc": "ENST00000397345.8:c.22144A>C",
        "hgvsp": "ENSP1:p.Thr7382Pro",
    }
    assert any("on the panel's BA1 exclusion list" in n for n in vcep.notes(neb))
    assert not any(
        "exclusion list" in n for n in vcep.notes({**neb, "hgvsc": "x:c.22249A>C", "hgvsp": "x:p.Thr7417Pro"})
    )
    # 'p.Arg41*' and VEP's 'p.Arg41Ter'; 'p.Asp181fs' and any frameshift there (ACTA1)
    acta1 = {
        "gene": "ACTA1",
        "acmg_points": 0,
        "gnomad_popmax_af": 0.03,
        "population_frequency": "frequent",
        "mane_select": "NM_001100.4",
    }
    assert any("exclusion list" in n for n in vcep.notes({**acta1, "hgvsp": "x:p.Arg41Ter"}))
    assert any("exclusion list" in n for n in vcep.notes({**acta1, "hgvsp": "x:p.Asp181ThrfsTer3"}))
    gjb2 = {"gene": "GJB2", "acmg_points": 0, "gnomad_popmax_af": 0.071, "population_frequency": "frequent",
            "rsid": "rs72474224", "alt": "T"}  # fmt: skip
    assert "on ClinGen's BA1 exception list" in vcep.notes(gjb2)[0]
    assert "its own BA1 rule" not in vcep.notes(gjb2)[0]
    assert "meets" not in note
    kcnq1 = {"gene": "KCNQ1", "acmg_points": 0, "gnomad_popmax_af": 0.02, "population_frequency": "frequent"}
    assert (
        "grpmax 2% is at or above its BA1 threshold (>= 0.4%): check its other BA1 conditions" in vcep.notes(kcnq1)[0]
    )
    # never a flat "meets": the panels' BA1 texts carry exceptions (antithrombin Cambridge II is pathogenic)
    serpinc1 = {"gene": "SERPINC1", "acmg_points": 0, "gnomad_popmax_af": 0.0024, "population_frequency": "frequent"}
    note = vcep.notes(serpinc1)[0]
    assert "is at or above its BA1 threshold (>= 0.2%): check its other BA1 conditions" in note
    assert "meets" not in note
    # conditions beyond the threshold (allele counts, males only, other datasets) are quoted
    rs1 = {"gene": "RS1", "acmg_points": 0, "gnomad_popmax_af": 0.001, "population_frequency": "rare"}
    assert "in males in population databases" in vcep.notes(rs1)[0]
    hbb = {"gene": "HBB", "acmg_points": 0, "gnomad_popmax_af": 0.045, "population_frequency": "frequent"}
    assert "variant present in ≥5 alleles" in vcep.notes(hbb)[0]
    assert vcep.ba1_above_5pct("VWF") == pytest.approx(0.1)
    assert vcep.ba1_above_5pct("MYH7") is None


def test_encoded_rules_say_what_the_panel_would_give():
    note = vcep.notes({"gene": "LDLR", "acmg_points": 4, "revel": 0.81, "spliceai_max": 0.0,
                       "gnomad_popmax_af": None, "acmg_pm2": "PM2_Supporting (absent from gnomAD v4)"})[0]  # fmt: skip
    assert "its rules give PP3_Supporting (REVEL 0.810), PM2 met" in note
    # the panels' PM2 asks for absence, whatever the coverage: the rules count it, with a word of caution
    low = vcep_rules.describe({"gene": "SOS1", "revel": 0.734, "spliceai_max": 0.0, "gnomad_popmax_af": None,
                               "acmg_pm2": "not applicable (absent from gnomAD v4, but few reads there: >= 20 reads in "
                                           "48% of exomes)"})  # fmt: skip
    assert low == "PP3_Supporting (REVEL 0.734), PM2 met (absent from gnomAD, but few reads there: check)"
    # the Myeloid Malignancy panel asks for 20x coverage; with no gnomAD answer at all, PM2 is not assessed
    few_reads = "not applicable (absent from gnomAD v4, but few reads there)"
    runx1 = vcep_rules.describe({"gene": "RUNX1", "revel": 0.3, "spliceai_max": 0.0, "gnomad_popmax_af": None,
                                 "acmg_pm2": few_reads})  # fmt: skip
    assert "PM2 not met (absent from gnomAD, but the panel requires 20x coverage)" in runx1
    nodata = vcep_rules.describe({"gene": "LDLR", "revel": 0.81, "spliceai_max": 0.0, "gnomad_popmax_af": None,
                                  "acmg_pm2": None})  # fmt: skip
    assert "PM2 not assessed (no gnomAD data)" in nodata
    ada = vcep_rules.describe({"gene": "ADA", "gnomad_popmax_af": 0.01})
    assert ada == (
        "no PP3/BP4 for missense variants, PM2 not met, BA1 met on grpmax (the panel uses the filtering allele "
        "frequency: check)"
    )


def test_weaker_residue_evidence_is_not_an_invitation_to_upgrade():
    note = vcep.notes({"gene": "LDLR", "acmg_points": 0, "acmg_ps1_pm5": "PM5_Supporting (p.X: Likely pathogenic)"})[0]
    assert "has no PM5 below Moderate: the tool lowered PM5 because the ClinVar reference is likely pathogenic" in note
    assert "would not apply" not in note  # the panel judges the reference itself (it applied PM5 in 25 of 134)
    assert "is weaker" not in note.split("PM2")[0]
    assert "classified pathogenic (by these guidelines)" in note  # the panel's condition on the reference


def test_conditions_the_tool_cannot_check_are_quoted():
    note = vcep.notes({"gene": "DICER1", "acmg_points": 0, "acmg_ps1_pm5": "PS1_Moderate (same change: LP)"})[0]
    assert "on PS1: “Likely pathogenic changes do not apply" in note
    # a reference classified under the panel's own specification (TP53, mismatch repair, hearing loss, F8)
    for gene, crit, words in (("TP53", "PM5", "according to the TP53 VCEP"),
                              ("MLH1", "PS1", "established by this VCEP"),
                              ("GJB2", "PS1", "pathogenicity by the HL specification"),
                              ("F8", "PS1", "based on F8 gene rule specifications")):  # fmt: skip
        code = "PM5 (p.X: Pathogenic, 3★)" if crit == "PM5" else "PS1 (same change: Pathogenic, 3★)"
        assert words in vcep.notes({"gene": gene, "acmg_points": 0, "acmg_ps1_pm5": code})[0], gene


def test_a_panel_that_uses_neither_pp3_nor_bp4_is_not_said_to_set_its_own_thresholds():
    note = vcep.notes({"gene": "GFAP", "acmg_points": 2})[0]
    assert "does not use PP3" in note
    assert "own predictor" not in note


def test_the_disease_types_of_a_specification_are_listed_together():
    assert "von Willebrand disease types 2A, 2B, 2M" in vcep.describe("VWF")


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


def test_exclusion_lists_are_shown_and_matched_precisely():
    neb = next(
        n
        for n in vcep.notes({"gene": "NEB", "acmg_points": 0, "gnomad_popmax_af": 0.01, "population_frequency": "x"})
        if "exclusion variants" in n
    )
    # the panel's numbering, with the MANE Select numbering where it differs
    assert "c.22249A>C (p.Thr7417Pro) [MANE NM_001164508.2: c.22144A>C (p.Thr7382Pro)]" in neb
    assert "BA1 exclusion variants (NM_001271208.2:" in neb  # NEB excludes them from BA1 only
    acta1 = {"gene": "ACTA1", "acmg_points": 0, "gnomad_popmax_af": 0.03, "population_frequency": "frequent",
             "mane_select": "NM_001100.4"}  # fmt: skip
    assert not any("exclusion list" in n for n in vcep.notes({**acta1, "hgvsp": "x:p.Asp1810MetfsTer4"}))
    assert not any("exclusion list" in n for n in vcep.notes({**acta1, "hgvsp": "x:p.Asp18ThrfsTer3"}))
    assert any("exclusion list" in n for n in vcep.notes({**acta1, "hgvsp": "x:p.Asp181ThrfsTer3"}))
    i18n.set_language("es")
    try:
        assert "variantes excluidas de BA1/BS1 (NM_000540.3:" in i18n.tr_text(
            "BA1/BS1 exclusion variants (NM_000540.3: c.1A>G)"
        )
    finally:
        i18n.set_language("en")


def test_quoted_ba1_conditions_keep_their_parts():
    def ba1(gene):
        return next(x["ba1"] for x in vcep.specifications(gene) if x.get("ba1"))

    # a cited PMID does not make a clause noise (GAMT); a calculation step does, unless it sets an allele count (DDC)
    assert any(">2000 alleles" in c for c in ba1("GAMT")["conditions"])
    assert any(">2000 alleles" in c for c in ba1("DDC")["conditions"])
    assert not ba1("ABCD1")["conditions"]  # the Whiffen-Ware calculation alone
    # a clause is not cut inside a parenthesis (KCNQ1), and the stop codon of p.Arg41* survives the cleaning (ACTA1)
    assert any("(in one of the 5 continental populations; African" in c for c in ba1("KCNQ1")["conditions"])
    acta1 = next(x["ba1"] for x in vcep.specifications("ACTA1") if (x.get("ba1") or {}).get("excluded"))
    assert "(p.Arg41*)" in acta1["quote"]


def test_encoded_ba1_on_the_filtering_allele_frequency_is_hedged():
    hnf4a = {"gene": "HNF4A", "gnomad_popmax_af": 0.000208, "revel": 0.919, "acmg_pm2": ""}
    assert vcep_rules.describe(hnf4a).endswith(
        "BA1 met on grpmax (the panel uses the filtering allele frequency: check)"
    )
    pah = {"gene": "PAH", "gnomad_popmax_af": 0.02, "revel": 0.919, "acmg_pm2": ""}  # BA1 on the allele frequency
    assert vcep_rules.describe(pah).endswith(", BA1 met")
    i18n.set_language("es")
    try:
        assert i18n.tr_text(vcep_rules.describe(hnf4a)).endswith(
            "BA1 cumplido con el grpmax (el panel usa la frecuencia alélica filtrada: compruébelo)"
        )
    finally:
        i18n.set_language("en")


def test_panel_ba1_above_5pct_sets_the_frequency_label_and_the_common_threshold(monkeypatch):
    import importlib

    P = importlib.import_module("missensexplorer.predict")
    assert P._common_af({"gene": "VWF"}, 0.05) == pytest.approx(0.1)
    assert P._common_af({"gene": "BRAF"}, 0.05) == 0.05
    assert P._common_af({"gene": "VWF"}, None) is None  # no frequency rule asked for
    monkeypatch.setattr(P.acmg, "frequency", lambda r: "common: BA1 (grpmax 7%)")
    monkeypatch.setattr(P.acmg, "ba1_grpmax", lambda r: 0.07)
    assert P._population_frequency({"gene": "VWF"}) == (
        "frequent: BS1 candidate (grpmax 7%); BA1 not met: the gene's expert panel sets BA1 above 10%"
    )
    assert P._population_frequency({"gene": "BRAF"}) == "common: BA1 (grpmax 7%)"
