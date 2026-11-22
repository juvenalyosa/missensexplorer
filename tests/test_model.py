"""Missense model, ACMG PP3/BP4 evidence and the population-frequency rule (offline)."""

import json
from pathlib import Path

import pytest

from missensexplorer import acmg, model
from missensexplorer.predict import apply_frequency, classify, colocated_info

DATA = Path(__file__).parent / "data"


def test_scorer_matches_scikit_learn():
    vectors = json.loads((DATA / "model_vectors.json").read_text())
    assert len(vectors) >= 20
    for v in vectors:
        ms = model.score(v)
        if v["alphamissense"] is None and v["revel"] is None:
            assert ms is None
        else:
            assert ms.p == pytest.approx(v["expected"], abs=1e-4)


def test_model_needs_alphamissense_or_revel():
    assert model.score({"cadd_phred": 30, "sift_score": 0.0}) is None
    assert model.score({"revel": 0.9}) is not None
    assert model.score({"alphamissense": float("nan"), "revel": None}) is None


def test_model_direction_and_contributions():
    hi = model.score({"alphamissense": 0.99, "revel": 0.95, "cadd_phred": 32, "sift_score": 0.0, "polyphen_score": 1.0})
    lo = model.score({"alphamissense": 0.05, "revel": 0.05, "cadd_phred": 5, "sift_score": 0.8, "polyphen_score": 0.0})
    assert hi.p > 0.95 and lo.p < 0.05
    assert dict(hi.contributions)["REVEL"] > 0 > dict(lo.contributions)["REVEL"]
    assert abs(hi.contributions[0][1]) >= abs(hi.contributions[-1][1])   # largest first
    assert model.describe(hi).startswith(f"model score {hi.p:.2f} (")


def test_classify_uses_model_and_falls_back_to_votes():
    r = {"alphamissense": 0.99, "alphamissense_class": "likely_pathogenic", "revel": 0.95, "cadd_phred": 32}
    label, score, basis = classify(r, {"missense_variant"})
    assert label == "Likely damaging" and score > 0.9 and basis.startswith("model score")
    assert "3 damaging / 0 benign predictors" in basis
    # neither AlphaMissense nor REVEL: the v2.1 vote decides
    assert classify({"cadd_phred": 30.0}, {"missense_variant"})[0] == "Likely damaging (single predictor)"
    weak = {"sift": "deleterious", "polyphen": "probably_damaging"}
    assert classify(weak, {"missense_variant"})[0] == "Likely damaging (limited evidence: SIFT/PolyPhen only)"


@pytest.mark.parametrize("row,code,pts", [
    ({"revel": 0.932}, "PP3_Strong", 4), ({"revel": 0.931}, "PP3_Moderate", 3), ({"revel": 0.773}, "PP3_Moderate", 2),
    ({"revel": 0.644}, "PP3_Supporting", 1), ({"revel": 0.643}, "none", 0), ({"revel": 0.291}, "none", 0),
    ({"revel": 0.290}, "BP4_Supporting", -1), ({"revel": 0.183}, "BP4_Moderate", -2),
    ({"revel": 0.052}, "BP4_Moderate", -3), ({"revel": 0.016}, "BP4_Strong", -4),
    ({"alphamissense": 0.990}, "PP3_Strong", 4), ({"alphamissense": 0.989}, "PP3_Moderate", 3),
    ({"alphamissense": 0.906}, "PP3_Moderate", 2), ({"alphamissense": 0.792}, "PP3_Supporting", 1),
    ({"alphamissense": 0.791}, "none", 0), ({"alphamissense": 0.169}, "BP4_Supporting", -1),
    ({"alphamissense": 0.099}, "BP4_Moderate", -2), ({"alphamissense": 0.05}, "BP4_Moderate", -3),
    ({"revel": 0.95, "alphamissense": 0.1}, "PP3_Strong", 4),          # REVEL first, as ClinGen recommends one tool
])
def test_pp3_bp4_intervals(row, code, pts):
    text, p = acmg.computational(row)
    assert text.startswith(code) and p == pts


def test_pp3_bp4_without_scores():
    assert acmg.computational({"cadd_phred": 30}) == (None, None)
    assert acmg.computational({"revel": 0.5})[0] == "none (REVEL 0.500, indeterminate)"


def test_frequency_evidence():
    assert acmg.frequency({"gnomad_af": 0.04, "gnomad_popmax_af": 0.062}) == "common: BA1 (grpmax 6.2%)"
    assert acmg.frequency({"gnomad_af": 0.004, "gnomad_popmax_af": 0.013}) == "frequent: BS1 candidate (grpmax 1.3%)"
    assert acmg.frequency({"gnomad_af": 1e-5, "gnomad_popmax_af": 2e-5}) == "rare (grpmax 0.002%)"
    assert acmg.frequency({}) == "not found in gnomAD"
    assert acmg.frequency({"chrom": "MT"}) is None


def test_frequency_rule_keeps_the_protein_level_call_visible():
    label, basis = apply_frequency("Likely damaging", "model score 0.95", 0.08, 0.05)
    assert label == "Likely benign (common in population)"
    assert "protein-level call: Likely damaging" in basis and basis.endswith("| model score 0.95")
    assert apply_frequency("Likely damaging", "x", 0.02, 0.05) == ("Likely damaging", "x")
    assert apply_frequency("Likely damaging", "x", 0.08, None) == ("Likely damaging", "x")
    assert apply_frequency("Uncertain (intermediate scores)", "x", 0.2, 0.05)[0].startswith("Likely benign")
    assert apply_frequency("Non-coding / no protein change", "", 0.3, 0.05)[0].startswith("Non-coding")


def test_grpmax_excludes_bottlenecked_groups():
    coloc = [{"id": "rs1", "frequencies": {"T": {"gnomade": 0.001, "gnomade_asj": 0.2, "gnomade_fin": 0.1,
                                                 "gnomadg_ami": 0.3, "gnomade_remaining": 0.05, "gnomadg_mid": 0.04,
                                                 "gnomade_nfe": 0.002, "gnomadg_afr": 0.004}}}]
    info = colocated_info(coloc, "T", "T")
    assert info["gnomad_popmax_af"] == 0.004 and info["gnomad_af"] == 0.001


def test_ba1_exceptions_are_not_called_benign():
    c282y = {"rsid": "rs1800562", "alt": "A", "gnomad_af": 0.03, "gnomad_popmax_af": 0.071}
    assert acmg.ba1_exception(c282y).startswith("HFE")
    assert acmg.frequency(c282y) == "frequent, but on the ClinGen BA1 exception list (grpmax 7.1%)"
    assert acmg.ba1_exception({**c282y, "alt": "C"}) is None              # another allele of the same rsID
    assert acmg.ba1_exception({"rsid": "rs13078881", "alt": "C"}).startswith("BTD")   # D444H = D424H on MANE
    assert acmg.ba1_exception({"rsid": "rs1, rs387906242", "alt": "TAAG"}).startswith("ACAD9")
