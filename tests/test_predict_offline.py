import json
from pathlib import Path

from missensexplorer import clinvar
from missensexplorer.inputs import parse_text
from missensexplorer.predict import _rows_for_variant, classify, colocated_info, short_protein_change
from missensexplorer.protein import context_window, features_at, mutate_protein

DATA = Path(__file__).parent / "data"


def test_rows_from_vep_region_fixture():
    vep = json.loads((DATA / "vep_braf_v600e_region.json").read_text())
    v = parse_text("7 140753336 A T").variants[0]
    rows = _rows_for_variant(v, vep, "GRCh38", None, False)
    assert len(rows) == 1
    r = rows[0]
    assert (r["gene"], r["protein_change"], r["hgvsp"], r["hgvsc"]) == ("BRAF", "V600E", "p.Val600Glu", "c.1799T>A")
    assert r["uniprot"] == "P15056" and r["protein_position"] == 600
    assert r["alphamissense_class"] == "likely_pathogenic" and r["cadd_phred"] == 29.8
    assert r["rsid"] == "rs113488022" and r["gnomad_af"] > 0
    assert r["grantham"] == 121 and r["_spdi"] == ("7", 140753335, "A", "T")
    label, score, basis = classify(r, r["_terms"])
    assert label == "Likely damaging" and score >= 0.8 and basis.startswith("model score") and "AlphaMissense" in basis


def test_gene_filter_drops_other_genes():
    vep = json.loads((DATA / "vep_braf_v600e_region.json").read_text())
    v = parse_text("7 140753336 A T").variants[0]
    assert _rows_for_variant(v, vep, "GRCh38", {"TP53"}, False) == []


def test_minus_strand_hgvs_alleles():
    vep = json.loads((DATA / "vep_hgvs_mixed.json").read_text())
    r213 = [x for x in vep if x["input"] == "TP53:p.R213*"]
    v = parse_text("TP53 R213*").variants[0]
    rows = _rows_for_variant(v, r213, "GRCh38", None, False)
    assert rows[0]["protein_change"] == "R213*"
    assert rows[0]["ref"] == "G" and rows[0]["alt"] == "A"  # reverse-complemented to the forward strand


def test_classify_categories():
    lof, s, _ = classify({"exon": "6/11", "nmd": "NMD likely"}, {"stop_gained"})
    assert lof.startswith("Loss of function") and s == 1.0
    benign = {"alphamissense": 0.07, "alphamissense_class": "likely_benign", "revel": 0.2, "cadd_phred": 5,
              "sift": "tolerated", "polyphen": "benign"}
    assert classify(benign, {"missense_variant"})[0] == "Likely benign"
    mixed = {"alphamissense": 0.99, "alphamissense_class": "likely_pathogenic", "revel": 0.1}
    assert classify(mixed, {"missense_variant"})[0].startswith("Uncertain")
    assert classify({"spliceai_max": 0.8}, {"synonymous_variant"})[0] == "Likely splice-altering"
    assert classify({}, {"intron_variant"})[0].startswith("Non-coding")


def test_short_protein_change():
    assert short_protein_change("ENSP1.1:p.Val600Glu") == "V600E"
    assert short_protein_change("p.Arg213Ter") == "R213*"
    assert short_protein_change("p.Thr125ProfsTer23") == "T125Pfs*23"
    assert short_protein_change("p.Leu12%3D") == "L12="


def test_colocated_strand_keys():
    coloc = [{"id": "rs1", "clin_sig_allele": "T:pathogenic;C:benign", "frequencies": {"A": {"gnomade": 1e-5}},
              "var_synonyms": {"ClinVar": ["VCV000013961", "RCV1"]}, "pubmed": [1, 2]}]
    info = colocated_info(coloc, "T", "A")
    assert info["_clinvar_vep"] == "pathogenic" and info["gnomad_af"] == 1e-5
    assert info["_vcvs"] == ["VCV000013961"] and info["pubmed_count"] == 2


def test_protein_mutation_helpers():
    seq = "MKVLAGRW"
    assert mutate_protein(seq, 3, 3, "V/E", "missense_variant") == ("MKELAGRW", "ok")
    mut, note = mutate_protein(seq, 4, 4, "L/*", "stop_gained")
    assert mut == "MKV" and note.startswith("truncated at residue 4")
    assert mutate_protein(seq, 3, 4, "VL/-", "inframe_deletion")[0] == "MKAGRW"
    assert mutate_protein(seq, 3, 3, "V/X", "frameshift_variant")[0] is None
    assert mutate_protein(seq, 3, 3, "Q/E", "missense_variant")[0] is None  # reference mismatch
    assert context_window(seq, 3, "V", "E", flank=2) == "MK[V>E]LA"


def test_features_at():
    entry = {"features": [
        {"type": "Domain", "start": 1, "end": 50, "description": "Kinase"},
        {"type": "Binding site", "start": 10, "end": 10, "description": "ATP"},
        {"type": "Disulfide bond", "start": 5, "end": 40, "description": ""},
        {"type": "Natural variant", "start": 10, "end": 10, "description": "in disease X", "original": "K",
         "alternatives": ["E"]},
    ]}
    f = features_at(entry, 10)
    assert f["key_site"] == "Binding site: ATP (10)" and "Domain: Kinase (1-50)" in f["protein_features"]
    assert f["uniprot_variants_here"] == "K10E: in disease X"
    assert "Disulfide" in features_at(entry, 40)["key_site"] and features_at(entry, 20)["key_site"] is None


def test_clinvar_match_by_spdi_and_cdna():
    recs = {"VCV1": {"spdi": "NC_000007.14:140753335:A:T", "title": "NM_004333.6(BRAF):c.1799T>A (p.Val600Glu)"},
            "VCV2": {"spdi": "NC_000007.14:140753335:A:C", "title": "NM_004333.6(BRAF):c.1799T>G (p.Val600Gly)"}}
    assert clinvar.match(recs, ["VCV2", "VCV1"], ("7", 140753335, "A", "T"), None, None) is recs["VCV1"]
    assert clinvar.match(recs, ["VCV1", "VCV2"], None, "NM_004333.6", "c.1799T>G") is recs["VCV2"]
    assert clinvar.stars("reviewed by expert panel") == 3


def test_sift_polyphen_only_as_fallback():
    # modern predictor present: SIFT/PolyPhen do not vote (benchmark: they add wrong calls)
    row = {"alphamissense": 0.1, "alphamissense_class": "likely_benign", "revel": 0.2,
           "sift": "deleterious", "polyphen": "probably_damaging"}
    assert classify(row, {"missense_variant"})[0] == "Likely benign"
    # nothing modern: SIFT/PolyPhen decide, and the label says the evidence is limited
    weak = {"sift": "deleterious", "polyphen": "probably_damaging"}
    assert classify(weak, {"missense_variant"})[0] == "Likely damaging (limited evidence: SIFT/PolyPhen only)"
    # modern scores in their intermediate range: uncertain, not "no scores"
    grey = {"revel": 0.5, "cadd_phred": 24.0, "sift": "deleterious"}
    label, score, _ = classify(grey, {"missense_variant"})
    assert label == "Uncertain (intermediate scores)" and 0.3 < score < 0.7
