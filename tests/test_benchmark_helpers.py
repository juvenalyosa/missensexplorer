"""Helpers of the benchmark scripts that the paper's numbers depend on (offline)."""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

import erepo_clinvar_history as H
import make_vcep_table as V


@pytest.mark.parametrize(
    ("sig", "level"),
    [
        ("Pathogenic", "P"),
        ("Likely pathogenic", "LP"),
        ("Pathogenic/Likely pathogenic", "P"),
        ("Likely pathogenic;Pathogenic", "P"),  # releases before 2017 join assertions with ';'
        ("Pathogenic;not provided", "P"),
        ("Pathogenic;risk factor", None),
        ("Conflicting interpretations of pathogenicity", None),
        ("Benign", None),
    ],
)
def test_clinvar_history_reads_old_and_new_significance_strings(sig, level):
    assert H.level(sig) == level


def _refs(*rows):
    cols = ["nm", "pos", "ref", "alt", "vid", "level", "stars", "cdna"]
    return H.index(pd.DataFrame([dict(zip(cols, r, strict=True)) for r in rows]))


V600 = {"nm": "NM_004333", "pos": 600, "ref": "Val", "alt": "Glu", "vid": "13961", "cdna": "c.1799T>A",
        "spliceai": 0.01}  # fmt: skip


def test_clinvar_history_applies_the_tools_rule():
    assert H.call(V600, _refs(("NM_004333", 600, "Val", "Gly", "1", "P", 3, "c.1799T>G"))) == "PM5"
    assert H.call(V600, _refs(("NM_004333", 600, "Val", "Gly", "1", "P", 1, "c.1799T>G"))) == "PM5_Supporting"
    assert H.call(V600, _refs(("NM_004333", 600, "Val", "Gly", "1", "LP", 1, "c.1799T>G"))) is None
    assert H.call(V600, _refs(("NM_004333", 600, "Val", "Glu", "2", "P", 2, "c.1799_1800delinsAA"))) == "PS1"
    # the variant itself never counts, by ClinVar id or by c. change
    assert H.call(V600, _refs(("NM_004333", 600, "Val", "Glu", "13961", "P", 3, "c.1799T>A"))) is None
    assert H.call(V600, _refs(("NM_004333", 600, "Val", "Glu", "", "P", 3, "c.1799T>A"))) is None
    assert H.call({**V600, "spliceai": 0.5}, _refs(("NM_004333", 600, "Val", "Gly", "1", "P", 3, "c.1799T>G"))) is None


@pytest.mark.parametrize(
    ("sid", "text", "value", "op", "faf"),
    [
        ("GN013", "Variant has a PopMax FAF ≥0.005 (0.5%) in gnomAD.", 0.005, ">=", True),
        ("GN081", "Appropriate to use for variants with a Popmax MAF of >0.1 in gnomAD.", 0.1, ">", False),
        ("GN017", "gnomAD Grpmax FAF ≥ 1:10,000 (≥ 0.01% or 0.0001)", 0.0001, ">=", True),
        (
            "GN122",
            "Grpmax Filtering Allele Frequency greater than or equal to 1.11 x 10<sup>-5</sup>",
            1.11e-5,
            ">=",
            True,
        ),
        ("GN020", "Grpmax Filtering AF **\\>.5%** in gnomAD", 0.005, ">", True),  # corrected by hand
    ],
)
def test_ba1_thresholds_are_read_from_the_specification_text(sid, text, value, op, faf):
    b = V.ba1_threshold(sid, text)
    assert b["value"] == pytest.approx(value)
    assert (b["op"], b["faf"]) == (op, faf)


@pytest.mark.parametrize(
    ("text", "excluded"),
    [
        ("Only applicable to synonymous or intronic variants. **Do not apply to missense variants.**", True),
        ("Missense: Do not use.", True),
        ("Due to some pathogenic variants having low REVEL scores, BP4 is not applied for missense variants.", True),
        ("For synonymous variants, CADD ≤22.7 is sufficient since REVEL and AlphaMissense are not applicable.", False),
        ("For missense variants, award PP3 if REVEL score is ≥0.7", False),
    ],
)
def test_pp3_bp4_not_for_missense_is_recognised(text, excluded):
    assert bool(V.NOT_MISSENSE.search(V._clean(text))) is excluded
