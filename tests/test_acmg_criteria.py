"""ACMG PS1/PM5 (ClinVar, same residue) and PM2 (direct gnomAD check), offline.

tests/data/clinvar_braf_codon600.json holds the real ClinVar records (E-utilities summaries) of the
pathogenic/likely pathogenic variants around BRAF codon 600, plus V600E (VCV000013961).
"""

import importlib
import json
from pathlib import Path

import httpx
import pytest

from missensexplorer import acmg, gnomad, residue
from missensexplorer import clinvar as cv
from missensexplorer.api import Http, SourceStatus
from missensexplorer.cache import Cache

P = importlib.import_module("missensexplorer.predict")
DATA = Path(__file__).parent / "data"
RECORDS = json.loads((DATA / "clinvar_braf_codon600.json").read_text())


def braf(pos, ref, alt, alt_aa, cdna, **kw):
    return {"consequence": "missense_variant", "protein_position": 600, "ref_aa": "V", "alt_aa": alt_aa,
            "chrom": "7", "pos": pos, "ref": ref, "alt": alt, "mane_select": "NM_004333.6",
            "hgvsc": f"ENST00000646891.2:{cdna}", "spliceai_max": 0.01, **kw}


V600E = braf(140753336, "A", "T", "E", "c.1799T>A")
V600L_GT = braf(140753337, "C", "A", "L", "c.1798G>T")   # likely pathogenic in ClinVar itself
V600L_GC = braf(140753337, "C", "G", "L", "c.1798G>C")   # pathogenic/likely pathogenic itself
V600M = braf(140753337, "C", "T", "M", "c.1798G>A")


def evidence(row, records=RECORDS):
    return residue.evidence(row, residue.references(row, records, "GRCh38"))


# --------------------------------------------------------------------------- PS1 / PM5

def test_pm5_from_other_pathogenic_changes_at_the_residue():
    refs = residue.references(V600E, RECORDS, "GRCh38")
    assert [x.protein for x in refs] == ["p.Val600Gly", "p.Val600Leu", "p.Val600Met", "p.Val600Leu"]
    # V600E itself (conflicting), Lys601 changes, in-frame indels and a record without criteria are left out
    assert not {"VCV000013961", "VCV000013966", "VCV000044817", "VCV000162793", "VCV002664272"} & {
        x.accession for x in refs}
    text, pts = residue.evidence(V600E, refs)
    assert text.startswith("PM5 (p.Val600Gly: Pathogenic, 3★, VCV000040389") and "and 2 more" in text
    assert pts == 2


def test_ps1_from_the_same_change_by_another_nucleotide():
    text, pts = evidence(V600L_GT)
    assert text == ("PS1 (same change by c.1798G>C: Pathogenic, 2★, VCV000376289)") and pts == 4


def test_ps1_moderate_when_the_same_change_is_only_likely_pathogenic():
    text, pts = evidence(V600L_GC)          # its own record does not count; c.1798G>T is likely pathogenic
    assert text.startswith("PS1_Moderate (same change by c.1798G>T: Likely pathogenic") and pts == 2


def test_pm5_supporting_when_only_likely_pathogenic_changes_are_known():
    lp_only = {k: v for k, v in RECORDS.items() if v["classification"] == "Likely pathogenic"}
    text, pts = evidence(V600E, lp_only)
    assert text.startswith("PM5_Supporting (p.Val600Met: Likely pathogenic, 1★") and pts == 1


def test_no_reference_other_transcript_and_splicing():
    assert evidence(V600E, {}) == ("none (no pathogenic missense variant at this residue in ClinVar)", 0)
    assert evidence({**V600E, "mane_select": "NM_001354609.2"})[1] == 0   # ClinVar names use NM_004333
    text, pts = evidence({**V600E, "spliceai_max": 0.35})
    assert text.startswith("not applied (SpliceAI 0.35") and pts == 0
    assert residue.evidence(V600E, None) == (None, None)


def test_window_only_for_missense_snvs():
    assert residue.window(V600E) == ("7", 140753334, 140753338)
    assert residue.window({**V600E, "consequence": "synonymous_variant"}) is None
    assert residue.window({**V600E, "ref": "AC", "alt": "TT"}) is None


def test_clinvar_window_search(tmp_path, monkeypatch):
    queries = []

    def handler(request):
        queries.append(request.content.decode())
        return httpx.Response(200, json={"esearchresult": {"idlist": [k[3:].lstrip("0") for k in RECORDS]}})

    http = Http()
    http.client = httpx.Client(transport=httpx.MockTransport(handler))
    monkeypatch.setattr(cv, "fetch", lambda http, cache, ids, status=None: {k: RECORDS[k] for k in ids})
    cache = Cache(tmp_path / "c.sqlite")
    w, far = ("7", 140753334, 140753338), ("7", 1, 5)
    out = cv.search_windows(http, cache, [w, far], "GRCh38")
    assert "chrpos38" in queries[0] and "clinsig_pathogenic" in queries[0]
    assert "VCV000040389" in out[w] and "VCV002664272" in out[w] and out[far] == []
    assert cv.search_windows(http, cache, [w], "GRCh38")[w] == out[w] and len(queries) == 1   # cached


# --------------------------------------------------------------------------- PM2 / gnomAD

def _gnomad_http(calls):
    def handler(request):
        q = json.loads(request.content)["query"]
        calls.append(q)
        if "variant(" in q:
            present = {"ac": 1, "an": 1461808, "filters": [], "populations": [{"id": "nfe", "ac": 1, "an": 1112010}]}
            return httpx.Response(200, json={"errors": [{"message": "Variant not found"}],
                                             "data": {"v0": {"exome": present, "genome": None}, "v1": None}})
        return httpx.Response(200, json={"data": {"c0": {"coverage": {
            "exome": [{"pos": 47890562, "over_20": 0.2}, {"pos": 47890563, "over_20": 0.93}], "genome": []}}}})

    http = Http()
    http.client = httpx.Client(transport=httpx.MockTransport(handler))
    return http


@pytest.fixture(autouse=True)
def no_pacing(monkeypatch):
    monkeypatch.setattr(gnomad, "_pace", lambda: None)


def test_gnomad_lookup_present_absent_and_coverage(tmp_path):
    calls = []
    st = SourceStatus()
    single, novel = ("7", 47890563, "G", "C"), ("7", 47890563, "G", "T")
    out = gnomad.lookup(_gnomad_http(calls), Cache(tmp_path / "c.sqlite"), [single, novel], "GRCh38", st)
    assert out[single]["found"] and out[single]["ac"] == 1
    assert out[single]["grpmax"] == pytest.approx(1 / 1112010)
    assert not out[novel]["found"] and out[novel]["over_20"] == {"exome": 0.93}
    assert gnomad.covered(out[novel]) and len(calls) == 2 and st.summary()["gnomAD"] == "ok"


def test_gnomad_failure_is_recorded(tmp_path):
    http = Http(max_retries=1)
    http.client = httpx.Client(transport=httpx.MockTransport(lambda r: httpx.Response(400, json={"errors": [{}]})))
    st = SourceStatus()
    assert gnomad.lookup(http, Cache(tmp_path / "c.sqlite"), [("1", 10, "A", "G")], "GRCh38", st) == {}
    assert st.failed("gnomAD", "1-10-A-G")


def test_pm2():
    absent = {"found": False, "dataset": "gnomAD v4", "over_20": {"exome": 0.93}}
    assert acmg.pm2({}, absent) == ("PM2_Supporting (absent from gnomAD v4; >= 20 reads in 93% of exomes)", 1)
    text, pts = acmg.pm2({}, {**absent, "over_20": {"exome": 0.1}})
    assert text.startswith("not applicable (absent from gnomAD v4, but few reads") and pts == 0
    assert acmg.pm2({"gnomad_af": 1e-6, "gnomad_popmax_af": 2e-6}, None)[0].startswith("depends on the disease")
    assert acmg.pm2({"gnomad_af": 0.02, "gnomad_popmax_af": 0.03}, None) == ("not met (frequent in gnomAD)", 0)
    assert acmg.pm2({}, None) == (None, None)


def test_codes_summary():
    row = {"acmg_pp3_bp4": "PP3_Strong (model 0.990)", "acmg_points": 4,
           "acmg_ps1_pm5": "PM5 (p.Val600Gly: Pathogenic, 3★, VCV000040389)",
           "acmg_pm2": "PM2_Supporting (absent from gnomAD v4; >= 20 reads in 93% of exomes)"}
    assert acmg.codes(row) == "PP3_Strong +4, PM5 +2, PM2_Supporting +1"
    assert acmg.codes({"acmg_pp3_bp4": "none (model 0.5, indeterminate)", "acmg_ps1_pm5": "none (...)"}) is None
    assert acmg.codes({"population_frequency": "common: BA1 (grpmax 6.2%)"}) == "BA1 (stand-alone benign)"


# --------------------------------------------------------------------------- in the pipeline

def test_pipeline_adds_ps1_pm5_and_pm2(monkeypatch):
    from test_robustness import BRAF

    class Vep:
        def __init__(self, *a, **k):
            pass

        def annotate(self, queries, params, progress=None, stage="", source="Ensembl VEP"):
            return {(kind, q): {"results": BRAF} for kind, qs in queries.items() for q in set(qs)}

    monkeypatch.setattr(P, "ensembl_release", lambda http, cache, asm: 116)
    monkeypatch.setattr(P, "_match_rate", lambda *a, **k: (1.0, {}))
    monkeypatch.setattr(P, "Vep", Vep)
    monkeypatch.setattr(P.cv, "fetch", lambda http, cache, ids, status=None: {k: RECORDS[k] for k in ids if k in RECORDS})
    monkeypatch.setattr(P.cv, "search_windows", lambda http, cache, ws, asm, status=None: {w: list(RECORDS) for w in ws})
    seen = []

    def lookup(http, cache, variants, assembly, status=None):
        seen.extend(variants)
        return {v: {"found": False, "dataset": "gnomAD v4", "over_20": {"exome": 0.9}} for v in variants}

    monkeypatch.setattr(P.gnomad, "lookup", lookup)
    res = P.predict(text="7 140753336 A T", options=P.PredictOptions(cache=False, protein_context=False, sequences=False))
    r = res.table.iloc[0]
    assert r["acmg_ps1_pm5"].startswith("PM5 (p.Val600Gly")
    assert seen == []                       # Ensembl already gave V600E a gnomAD frequency: no direct check
    assert r["acmg_pm2"].startswith("depends on the disease")
    assert "PM5 +2" in r["acmg_codes"]
