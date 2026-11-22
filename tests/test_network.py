"""Live end-to-end tests against the public APIs:  pytest -m network"""
import pytest

from missensexplorer import PredictOptions, predict

pytestmark = pytest.mark.network


def test_end_to_end_mixed_inputs(tmp_path):
    res = predict(text="BRAF V600E\nchr17:7676154 G>C\nFAKEGENE V600E",
                  options=PredictOptions(cache_path=str(tmp_path / "c.sqlite")))
    t = res.table.set_index("protein_change")
    assert t.loc["V600E", "prediction"] == "Likely damaging"
    assert t.loc["P72R", "prediction"] == "Likely benign"
    assert t.loc["P72R", "clinvar"] == "Benign"
    assert t.loc["V600E", "plddt"] > 0 and t.loc["V600E", "protein_name"]
    assert len(res.errors) == 1 and 0 in res.sequences


def test_grch37_autodetect(tmp_path):
    # BRAF V600E in GRCh37 coordinates; no header -> build detected from REF alleles
    res = predict(text="7:140453136 A>T\n17:7577120 C>T",
                  options=PredictOptions(cache_path=str(tmp_path / "c.sqlite")))
    assert res.assembly == "GRCh37"
    assert set(res.table["gene"]) == {"BRAF", "TP53"}
