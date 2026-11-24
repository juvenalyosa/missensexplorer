"""Edge cases on real VEP responses recorded by tests/data/record_fixtures.py (offline)."""

import json
from pathlib import Path

import pytest

from missensexplorer.inputs import parse_text, read_vcf
from missensexplorer.predict import _rows_for_variant, classify

DATA = Path(__file__).parent / "data"


def _load(name):
    return {r["input"]: r for r in json.loads((DATA / f"{name}.json").read_text())}


HGVS = _load("vep_edge_hgvs")
REGION = _load("vep_edge_region")
GRCH37 = _load("vep_edge_grch37")


def _row(text, results, gene, assembly="GRCh38"):
    v = parse_text(text).variants[0]
    rows = [r for r in _rows_for_variant(v, results, assembly, None, False) if r["gene"] == gene]
    assert len(rows) == 1, rows
    r = rows[0]
    label, score, _ = classify(r, r["_terms"])
    return r, label, score


@pytest.mark.parametrize(
    ("text", "key", "gene", "change", "label"),
    [
        ("BRCA1:c.68_69del", "BRCA1:c.68_69del", "BRCA1", "E23Vfs*17", "Loss of function (truncating)"),
        ("TP53:c.1A>G", "TP53:c.1A>G", "TP53", "M1?", "Loss of function (start codon lost)"),
        ("TP53:c.375+1G>A", "TP53:c.375+1G>A", "TP53", None, "Loss of function (canonical splice site)"),
        ("TP53:c.375G>A", "TP53:c.375G>A", "TP53", "T125=", "Likely splice-altering"),
    ],
)
def test_hgvs_edge_cases(text, key, gene, change, label):
    r, got, score = _row(text, [HGVS[key]], gene)
    assert r["protein_change"] == change
    assert got == label
    if label.startswith("Loss"):
        assert score == 1.0


def test_frameshift_nmd_estimate():
    r, _, _ = _row("BRCA1:c.68_69del", [HGVS["BRCA1:c.68_69del"]], "BRCA1")
    assert r["nmd"].startswith("NMD likely")


def test_synonymous_splice_score_is_used():
    r, _, _ = _row("TP53:c.375G>A", [HGVS["TP53:c.375G>A"]], "TP53")
    assert r["spliceai_max"] == pytest.approx(0.58)


def test_mitochondrial_missense():
    r, label, _ = _row("MT 3460 G A", [REGION["MT 3460 . G A . . ."]], "MT-ND1")
    assert r["protein_change"] == "A52T"
    assert r["uniprot"] == "P03886"
    assert label.startswith("Likely damaging")  # SIFT + PolyPhen only
    assert r["grantham"] == 58


def test_chromosome_x_hemizygous_from_vcf():
    vcf = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBOY",
        "X\t31478302\t.\tG\tA\t.\tPASS\t.\tGT\t1",
    ]
    v = read_vcf(vcf).variants[0]
    assert v.zygosity == "hemizygous"
    rows = [
        r for r in _rows_for_variant(v, [REGION["X 31478302 . G A . . ."]], "GRCh38", None, False) if r["gene"] == "DMD"
    ]
    assert rows[0]["protein_change"] == "E2914V"
    assert rows[0]["zygosity"] == "hemizygous"


def test_multiallelic_record_gives_one_row_per_allele():
    vcf = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        "17\t7673802\t.\tC\tT,A\t.\tPASS\t.\tGT\t1/2",
    ]
    parsed = read_vcf(vcf)
    assert [v.alt for v in parsed.variants] == ["T", "A"]
    assert {v.zygosity for v in parsed.variants} == {"heterozygous"}
    both = [REGION["17 7673802 . C T . . ."], REGION["17 7673802 . C A . . ."]]
    changes = []
    for v in parsed.variants:  # even if VEP returned both alleles, each row keeps its own
        rows = [r for r in _rows_for_variant(v, both, "GRCh38", None, False) if r["gene"] == "TP53"]
        changes.append([r["protein_change"] for r in rows])
    assert changes == [["R273H"], ["R273L"]]


def test_grch37_braf():
    r, label, _ = _row("7 140453136 A T", [GRCH37["7 140453136 . A T . . ."]], "BRAF", "GRCh37")
    assert r["protein_change"] == "V600E"
    assert label == "Likely damaging"
    assert r["alphamissense"] > 0.9
    assert r["assembly"] == "GRCh37"
