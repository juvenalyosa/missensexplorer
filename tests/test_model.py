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


def test_score_with_frequency_matches_scikit_learn():
    assert model.FREQUENCY_FITTED, "the model was retrained: run benchmarks/train_frequency_score.py"
    vectors = json.loads((DATA / "frequency_vectors.json").read_text())
    assert len(vectors) >= 20
    for v in vectors:
        assert model.with_frequency(v["p"], v["grpmax"]) == pytest.approx(v["expected"], abs=1e-4)


def test_score_with_frequency_direction():
    f = model.with_frequency
    assert f(0.9, 1e-5) > f(0.5, 1e-5) > f(0.1, 1e-5)  # higher model score, higher
    assert f(0.9, None) > f(0.9, 1e-4) > f(0.9, 1e-2)  # more frequent, lower; absent counts as 0
    assert f(None, 1e-5) is None


def test_model_needs_alphamissense_or_revel():
    assert model.score({"cadd_phred": 30, "sift_score": 0.0}) is None
    assert model.score({"revel": 0.9}) is not None
    assert model.score({"alphamissense": float("nan"), "revel": None}) is None


def test_model_direction_and_contributions():
    hi = model.score({"alphamissense": 0.99, "revel": 0.95, "cadd_phred": 32, "sift_score": 0.0, "polyphen_score": 1.0})
    lo = model.score({"alphamissense": 0.05, "revel": 0.05, "cadd_phred": 5, "sift_score": 0.8, "polyphen_score": 0.0})
    assert hi.p > 0.95
    assert lo.p < 0.05
    assert dict(hi.contributions)["REVEL"] > 0 > dict(lo.contributions)["REVEL"]
    assert abs(hi.contributions[0][1]) >= abs(hi.contributions[-1][1])  # largest first
    assert model.describe(hi).startswith(f"model score {hi.p:.2f} (")


def test_classify_uses_model_and_falls_back_to_votes():
    r = {"alphamissense": 0.99, "alphamissense_class": "likely_pathogenic", "revel": 0.95, "cadd_phred": 32}
    label, score, basis = classify(r, {"missense_variant"})
    assert label == "Likely damaging"
    assert score > 0.9
    assert basis.startswith("model score")
    assert "3 damaging / 0 benign predictors" in basis
    # neither AlphaMissense nor REVEL: the v2.1 vote decides
    assert classify({"cadd_phred": 30.0}, {"missense_variant"})[0] == "Likely damaging (single predictor)"
    weak = {"sift": "deleterious", "polyphen": "probably_damaging"}
    assert classify(weak, {"missense_variant"})[0] == "Likely damaging (limited evidence: SIFT/PolyPhen only)"


@pytest.mark.parametrize(
    ("row", "code", "pts"),
    [
        ({"revel": 0.932}, "PP3_Strong", 4),
        ({"revel": 0.931}, "PP3_Moderate", 3),
        ({"revel": 0.773}, "PP3_Moderate", 2),
        ({"revel": 0.644}, "PP3_Supporting", 1),
        ({"revel": 0.643}, "none", 0),
        ({"revel": 0.291}, "none", 0),
        ({"revel": 0.290}, "BP4_Supporting", -1),
        ({"revel": 0.183}, "BP4_Moderate", -2),
        ({"revel": 0.052}, "BP4_Moderate", -3),
        ({"revel": 0.016}, "BP4_Strong", -4),
        ({"alphamissense": 0.990}, "PP3_Strong", 4),
        ({"alphamissense": 0.989}, "PP3_Moderate", 3),
        ({"alphamissense": 0.906}, "PP3_Moderate", 2),
        ({"alphamissense": 0.792}, "PP3_Supporting", 1),
        ({"alphamissense": 0.791}, "none", 0),
        ({"alphamissense": 0.169}, "BP4_Supporting", -1),
        ({"alphamissense": 0.099}, "BP4_Moderate", -2),
        ({"alphamissense": 0.05}, "BP4_Moderate", -3),
        ({"revel": 0.95, "alphamissense": 0.1}, "PP3_Strong", 4),  # REVEL first, as ClinGen recommends one tool
    ],
)
def test_pp3_bp4_intervals(row, code, pts):
    text, p = acmg.computational({**row, "spliceai_max": 0.01}, "revel")
    assert text.startswith(code)
    assert p == pts


def test_bp4_needs_spliceai_to_rule_out_splicing():
    text, p = acmg.computational({"revel": 0.016}, "revel")  # no SpliceAI score
    assert text == "none (REVEL 0.016; BP4_Strong withheld: no SpliceAI score to rule out a splicing effect)"
    assert p == 0
    assert acmg.computational({"revel": 0.016, "spliceai_max": 0.3}, "revel")[1] == 0
    assert acmg.computational({"revel": 0.016, "spliceai_max": 0.1}, "revel")[1] == -4


def test_pp3_bp4_without_scores():
    for source in acmg.SOURCES:
        assert acmg.computational({"cadd_phred": 30}, source) == (None, None)
    assert acmg.computational({"revel": 0.5}, "revel")[0] == "none (REVEL 0.500, indeterminate)"
    with pytest.raises(ValueError, match="unknown PP3/BP4 source"):
        acmg.computational({"revel": 0.5}, "cadd")


@pytest.mark.parametrize(
    ("evidence", "used"),
    [
        ("PP3_Strong (model 0.994)", "model"),
        ("BP4_Moderate (model 0.052)", "model"),
        ("PP3_Moderate (REVEL 0.930)", "clingen"),
        ("none (AlphaMissense 0.500, indeterminate)", "clingen"),
        (None, None),  # not missense: no PP3/BP4
        ("", None),
        (float("nan"), None),  # a missing cell of the results table
    ],
)
def test_thresholds_used_names_the_pp3_bp4_calibration(evidence, used):
    assert acmg.thresholds_used(evidence) == used


def test_model_pp3_bp4_follows_its_calibration():
    assert acmg.MODEL_CALIBRATED, "the model was retrained: run benchmarks/calibrate_acmg.py"
    cal = acmg.CALIBRATION
    assert [p for _, p in cal["pp3"]] == sorted((p for _, p in cal["pp3"]), reverse=True)
    assert cal["pp3"][-1][0] > cal["bp4"][-1][0]  # weakest PP3 cut-off above the weakest BP4 one
    seen = set()
    for v in json.loads((DATA / "model_vectors.json").read_text()):
        text, pts = acmg.computational({**v, "spliceai_max": 0.0})  # the model is the default source
        if v["alphamissense"] is None and v["revel"] is None:
            assert (text, pts) == (None, None)
            continue
        p = round(model.score(v).p, 3)
        assert f"model {p:.3f}" in text
        assert pts == acmg._points(p, cal["pp3"], cal["bp4"])
        seen.add((pts > 0) - (pts < 0))
    assert seen == {-1, 0, 1}  # the vectors reach PP3, BP4 and no evidence


def test_frequency_evidence():
    assert (
        acmg.frequency({"gnomad_af": 0.04, "gnomad_popmax_af": 0.062})
        == "common: BA1 (grpmax 6.2%; allele numbers not checked in gnomAD)"
    )
    checked = {"gnomad_af": 0.04, "gnomad_popmax_af": 0.062, "_grpmax_2000": (0.061, 30000)}
    assert acmg.frequency(checked) == "common: BA1 (grpmax 6.2%)"
    assert acmg.frequency({"gnomad_af": 0.004, "gnomad_popmax_af": 0.013}) == "frequent: BS1 candidate (grpmax 1.3%)"
    assert acmg.frequency({"gnomad_af": 1e-5, "gnomad_popmax_af": 2e-5}) == "rare (grpmax 0.002%)"
    assert acmg.frequency({}) == "no gnomAD frequency in Ensembl's data (not checked in gnomAD)"
    assert acmg.frequency({"_gnomad": {"found": False, "dataset": "gnomAD v4"}}) == "absent from gnomAD v4"
    assert acmg.frequency({"chrom": "MT"}) is None


def test_frequency_rule_keeps_the_protein_level_call_visible():
    label, basis = apply_frequency("Likely damaging", "model score 0.95", 0.08, 0.05)
    assert label == "Likely benign (common in population)"
    assert "protein-level call: Likely damaging" in basis
    assert basis.endswith("| model score 0.95")
    assert apply_frequency("Likely damaging", "x", 0.02, 0.05) == ("Likely damaging", "x")
    assert apply_frequency("Likely damaging", "x", 0.08, None) == ("Likely damaging", "x")
    assert apply_frequency("Uncertain (intermediate scores)", "x", 0.2, 0.05)[0].startswith("Likely benign")
    assert apply_frequency("Non-coding / no protein change", "", 0.3, 0.05)[0].startswith("Non-coding")


def test_grpmax_excludes_bottlenecked_groups():
    coloc = [
        {
            "id": "rs1",
            "frequencies": {
                "T": {
                    "gnomade": 0.001,
                    "gnomade_asj": 0.2,
                    "gnomade_fin": 0.1,
                    "gnomadg_ami": 0.3,
                    "gnomade_remaining": 0.05,
                    "gnomadg_mid": 0.04,
                    "gnomade_nfe": 0.002,
                    "gnomadg_afr": 0.004,
                }
            },
        }
    ]
    info = colocated_info(coloc, "T", "T")
    assert info["gnomad_popmax_af"] == 0.004
    assert info["gnomad_af"] == 0.001


def test_ba1_exceptions_are_not_called_benign():
    c282y = {"rsid": "rs1800562", "alt": "A", "gnomad_af": 0.03, "gnomad_popmax_af": 0.071}
    assert acmg.ba1_exception(c282y).startswith("HFE")
    assert acmg.frequency(c282y) == "frequent, but on the ClinGen BA1 exception list (grpmax 7.1%)"
    assert acmg.ba1_exception({**c282y, "alt": "C"}) is None  # another allele of the same rsID
    assert acmg.ba1_exception({"rsid": "rs13078881", "alt": "C"}).startswith("BTD")  # D444H = D424H on MANE
    assert acmg.ba1_exception({"rsid": "rs1, rs387906242", "alt": "TAAG"}).startswith("ACAD9")
