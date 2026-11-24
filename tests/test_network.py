"""Live end-to-end tests against the public APIs:  pytest -m network"""

import pytest

from missensexplorer import PredictOptions, predict

pytestmark = pytest.mark.network


def test_end_to_end_mixed_inputs(tmp_path):
    res = predict(
        text="BRAF V600E\nchr17:7676154 G>C\nFAKEGENE V600E",
        options=PredictOptions(cache_path=str(tmp_path / "c.sqlite")),
    )
    t = res.table.set_index("protein_change")
    assert t.loc["V600E", "prediction"] == "Likely damaging"
    assert t.loc["P72R", "prediction"] == "Likely benign"
    assert t.loc["P72R", "clinvar"] == "Benign"
    assert t.loc["V600E", "plddt"] > 0
    assert t.loc["V600E", "protein_name"]
    assert len(res.errors) == 1
    assert 0 in res.sequences


def test_grch37_autodetect(tmp_path):
    # BRAF V600E in GRCh37 coordinates; no header -> build detected from REF alleles
    res = predict(text="7:140453136 A>T\n17:7577120 C>T", options=PredictOptions(cache_path=str(tmp_path / "c.sqlite")))
    assert res.assembly == "GRCh37"
    assert set(res.table["gene"]) == {"BRAF", "TP53"}


def test_ps1_pm5_and_pm2_live(tmp_path):
    # BRAF V600L by c.1798G>T: c.1798G>C makes the same change and is pathogenic (PS1); absent from gnomAD
    # (PM2). PKD1L1 P885R is a gnomAD v4 singleton that Ensembl gives no frequency: found by the direct check.
    res = predict(
        text="7:140753337 C>A\n7:47890563 G>C",
        options=PredictOptions(cache_path=str(tmp_path / "c.sqlite"), protein_context=False),
    )
    t = res.table.set_index("protein_change")
    assert t.loc["V600L", "acmg_ps1_pm5"].startswith("PS1 (same change by c.1798G>C")
    assert t.loc["V600L", "acmg_pm2"].startswith("PM2_Supporting (absent from gnomAD v4")
    assert t.loc["P885R", "population_frequency"].startswith("rare (grpmax")
    assert t.loc["P885R", "acmg_pm2"].startswith("depends on the disease")


# --------------------------------------------------------------------------- API contracts
# These check that the external services still return the fields MissensExplorer reads.
# They run every night on GitHub (.github/workflows/nightly.yml) to catch API changes early.


def test_contract_vep_plugins():
    from missensexplorer.api import ENSEMBL, Http
    from missensexplorer.vep import BASE_PARAMS, PREDICTOR_PARAMS

    http = Http()
    data = http.post_json(
        ENSEMBL["GRCh38"] + "/vep/human/region",
        {"variants": ["17 7673802 . C T . . ."]},
        params={**BASE_PARAMS, **PREDICTOR_PARAMS},
        headers={"Content-Type": "application/json"},
    )
    r = data[0]
    tc = next(t for t in r["transcript_consequences"] if t.get("gene_symbol") == "TP53")
    # PolyPhen-2 is not computed for every transcript, so it is optional (as in the pipeline)
    for key in (
        "hgvsp",
        "hgvsc",
        "amino_acids",
        "protein_start",
        "protein_id",
        "mane_select",
        "swissprot",
        "sift_prediction",
        "cadd_phred",
        "revel_score",
        "eve_class",
        "popeve_score",
    ):
        assert key in tc, key
    assert {"am_class", "am_pathogenicity"} <= set(tc["alphamissense"])
    assert {"DS_AG", "DS_AL", "DS_DG", "DS_DL"} <= set(tc["spliceai"])
    coloc = r["colocated_variants"]
    assert any("frequencies" in c for c in coloc)
    assert any("ClinVar" in (c.get("var_synonyms") or {}) for c in coloc)


def test_contract_clinvar_uniprot_alphafold():
    from missensexplorer.api import ALPHAFOLD, UNIPROT, Http
    from missensexplorer.clinvar import EUTILS
    from missensexplorer.protein import UNIPROT_FIELDS

    http = Http()
    cv = http.request("POST", EUTILS, data={"db": "clinvar", "id": "12347", "retmode": "json"}).json()["result"]
    rec = cv["12347"]
    assert "germline_classification" in rec
    assert "canonical_spdi" in rec["variation_set"][0]
    up = http.get_json(
        f"{UNIPROT}/uniprotkb/accessions", params={"accessions": "P04637", "fields": UNIPROT_FIELDS, "format": "json"}
    )["results"][0]
    assert up["sequence"]["length"] == 393
    assert any(f["type"] == "DNA binding" for f in up["features"])
    af = http.get_json(f"{ALPHAFOLD}/api/prediction/P04637")[0]
    assert af["pdbUrl"].endswith(".pdb")
    assert "plddtDocUrl" in af
    assert len(http.get_json(af["plddtDocUrl"])["confidenceScore"]) == 393
