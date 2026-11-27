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
    return {
        "consequence": "missense_variant",
        "protein_position": 600,
        "ref_aa": "V",
        "alt_aa": alt_aa,
        "chrom": "7",
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "mane_select": "NM_004333.6",
        "hgvsc": f"ENST00000646891.2:{cdna}",
        "spliceai_max": 0.01,
        **kw,
    }


V600E = braf(140753336, "A", "T", "E", "c.1799T>A")
V600L_GT = braf(140753337, "C", "A", "L", "c.1798G>T")  # likely pathogenic in ClinVar itself
V600L_GC = braf(140753337, "C", "G", "L", "c.1798G>C")  # pathogenic/likely pathogenic itself
V600M = braf(140753337, "C", "T", "M", "c.1798G>A")


def evidence(row, records=RECORDS):
    return residue.evidence(row, residue.references(row, records, "GRCh38"))


# --------------------------------------------------------------------------- PS1 / PM5


def test_pm5_from_other_pathogenic_changes_at_the_residue():
    refs = residue.references(V600E, RECORDS, "GRCh38")
    assert [x.protein for x in refs] == ["p.Val600Gly", "p.Val600Leu", "p.Val600Met", "p.Val600Leu"]
    # V600E itself (conflicting), Lys601 changes, in-frame indels and a record without criteria are left out
    assert not {"VCV000013961", "VCV000013966", "VCV000044817", "VCV000162793", "VCV002664272"} & {
        x.accession for x in refs
    }
    text, pts = residue.evidence(V600E, refs)
    # the single-star likely pathogenic references are not used
    assert text == "PM5 (p.Val600Gly: Pathogenic, 3★, VCV000040389; p.Val600Leu: Pathogenic, 2★, VCV000376289)"
    assert pts == 2


def test_ps1_from_the_same_change_by_another_nucleotide():
    text, pts = evidence(V600L_GT)
    assert text == ("PS1 (same change by c.1798G>C: Pathogenic, 2★, VCV000376289)")
    assert pts == 4


def test_a_single_star_likely_pathogenic_same_change_is_not_used():
    text, pts = evidence(V600L_GC)  # its own record does not count; c.1798G>T is likely pathogenic, 1 star
    assert text.startswith("PM5 (p.Val600Gly: Pathogenic, 3★")
    assert pts == 2


def ref(level, stars, same, protein="p.Val600Leu"):
    return residue.Reference("VCV000000001", "c.1798G>T", protein, level, stars, same)


@pytest.mark.parametrize(
    ("level", "stars", "same", "code"),
    [
        ("P", 2, True, "PS1"),
        ("P", 1, True, "PS1_Moderate"),  # one review star: one step lower
        ("LP", 2, True, "PS1_Moderate"),
        ("LP", 1, True, "none"),  # likely pathogenic with one star: not used
        ("P", 3, False, "PM5"),
        ("P", 1, False, "PM5_Supporting"),
        ("LP", 2, False, "PM5_Supporting"),
        ("LP", 1, False, "none"),
    ],
)
def test_reference_strength_depends_on_level_and_review_stars(level, stars, same, code):
    text, _ = residue.evidence(V600E, [ref(level, stars, same, "p.Val600Glu" if same else "p.Val600Leu")])
    assert text.split(" ", 1)[0] == code


def test_none_when_only_single_star_likely_pathogenic_changes_are_known():
    lp_only = {k: v for k, v in RECORDS.items() if v["classification"] == "Likely pathogenic"}
    text, pts = evidence(V600E, lp_only)
    assert text.startswith("none (only single-star likely pathogenic variants at this residue: p.Val600Met")
    assert pts == 0


def test_no_reference_other_transcript_and_splicing():
    assert evidence(V600E, {}) == ("none (no pathogenic missense variant at this residue in ClinVar)", 0)
    assert evidence({**V600E, "mane_select": "NM_001354609.2"})[1] == 0  # ClinVar names use NM_004333
    text, pts = evidence({**V600E, "spliceai_max": 0.35})
    assert text.startswith("not applied (SpliceAI 0.35")
    assert pts == 0
    assert residue.evidence(V600E, None) == (None, None)


def test_window_only_for_missense_snvs():
    assert residue.window(V600E) == ("7", 140753334, 140753338)
    assert residue.window({**V600E, "consequence": "synonymous_variant"}) is None
    assert residue.window({**V600E, "ref": "AC", "alt": "TT"}) is None


def test_clinvar_window_search(tmp_path, monkeypatch, request):
    queries = []

    def handler(request):
        queries.append(request.content.decode())
        return httpx.Response(200, json={"esearchresult": {"idlist": [k[3:].lstrip("0") for k in RECORDS]}})

    http = Http()
    http.client = httpx.Client(transport=httpx.MockTransport(handler))
    monkeypatch.setattr(cv, "fetch", lambda http, cache, ids, status=None: {k: RECORDS[k] for k in ids})
    cache = Cache(tmp_path / "c.sqlite")
    request.addfinalizer(cache.close)
    w, far = ("7", 140753334, 140753338), ("7", 1, 5)
    out = cv.search_windows(http, cache, [w, far], "GRCh38")
    assert "chrpos38" in queries[0]
    assert "clinsig_pathogenic" in queries[0]
    assert "VCV000040389" in out[w]
    assert "VCV002664272" in out[w]
    assert out[far] == []
    assert cv.search_windows(http, cache, [w], "GRCh38")[w] == out[w]
    assert len(queries) == 1


# --------------------------------------------------------------------------- PM2 / gnomAD


def _gnomad_http(calls):
    def handler(request):
        q = json.loads(request.content)["query"]
        calls.append(q)
        if "variant(" in q:
            present = {"ac": 1, "an": 1461808, "filters": [], "populations": [{"id": "nfe", "ac": 1, "an": 1112010}]}
            return httpx.Response(
                200,
                json={
                    "errors": [{"message": "Variant not found"}],
                    "data": {"v0": {"exome": present, "genome": None}, "v1": None},
                },
            )
        return httpx.Response(
            200,
            json={
                "data": {
                    "c0": {
                        "coverage": {
                            "exome": [{"pos": 47890562, "over_20": 0.2}, {"pos": 47890563, "over_20": 0.93}],
                            "genome": [],
                        }
                    }
                }
            },
        )

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
    with Cache(tmp_path / "c.sqlite") as cache:
        out = gnomad.lookup(_gnomad_http(calls), cache, [single, novel], "GRCh38", st)
    assert out[single]["found"]
    assert out[single]["ac"] == 1
    assert out[single]["grpmax"] == pytest.approx(1 / 1112010)
    assert not out[novel]["found"]
    assert out[novel]["over_20"] == {"exome": 0.93}
    assert gnomad.covered(out[novel])
    assert len(calls) == 2
    assert st.summary()["gnomAD"] == "ok"


def test_gnomad_failure_is_recorded(tmp_path):
    http = Http(max_retries=1)
    http.client = httpx.Client(transport=httpx.MockTransport(lambda r: httpx.Response(400, json={"errors": [{}]})))
    st = SourceStatus()
    with Cache(tmp_path / "c.sqlite") as cache:
        assert gnomad.lookup(http, cache, [("1", 10, "A", "G")], "GRCh38", st) == {}
    assert st.failed("gnomAD", "1-10-A-G")


def test_pm2():
    absent = {"found": False, "dataset": "gnomAD v4", "over_20": {"exome": 0.93}}
    assert acmg.pm2({}, absent) == ("PM2_Supporting (absent from gnomAD v4; >= 20 reads in 93% of exomes)", 1)
    text, pts = acmg.pm2({}, {**absent, "over_20": {"exome": 0.1}})
    assert text.startswith("not applicable (absent from gnomAD v4, but few reads")
    assert pts == 0
    assert acmg.pm2({"gnomad_af": 1e-6, "gnomad_popmax_af": 2e-6}, None)[0].startswith("depends on the disease")
    assert acmg.pm2({"gnomad_af": 0.02, "gnomad_popmax_af": 0.03}, None) == ("not met (frequent in gnomAD)", 0)
    assert acmg.pm2({}, None) == (None, None)


def test_codes_summary():
    row = {
        "acmg_pp3_bp4": "PP3_Strong (model 0.990)",
        "acmg_points": 4,
        "acmg_ps1_pm5": "PM5 (p.Val600Gly: Pathogenic, 3★, VCV000040389)",
        "acmg_pm2": "PM2_Supporting (absent from gnomAD v4; >= 20 reads in 93% of exomes)",
    }
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

    monkeypatch.setattr(P, "ensembl_release", lambda http, cache, asm, **k: 116)
    monkeypatch.setattr(P, "_match_rate", lambda *a, **k: (1.0, {}))
    monkeypatch.setattr(P, "Vep", Vep)
    monkeypatch.setattr(
        P.cv, "fetch", lambda http, cache, ids, status=None: {k: RECORDS[k] for k in ids if k in RECORDS}
    )
    monkeypatch.setattr(
        P.cv, "search_windows", lambda http, cache, ws, asm, status=None: {w: list(RECORDS) for w in ws}
    )
    seen = []

    def lookup(http, cache, variants, assembly, status=None):
        seen.extend(variants)
        return {v: {"found": False, "dataset": "gnomAD v4", "over_20": {"exome": 0.9}} for v in variants}

    monkeypatch.setattr(P.gnomad, "lookup", lookup)
    res = P.predict(
        text="7 140753336 A T", options=P.PredictOptions(cache=False, protein_context=False, sequences=False)
    )
    r = res.table.iloc[0]
    assert r["acmg_ps1_pm5"].startswith("PM5 (p.Val600Gly")
    assert seen == []  # Ensembl already gave V600E a gnomAD frequency: no direct check
    assert r["acmg_pm2"].startswith("depends on the disease")
    assert "PM5 +2" in r["acmg_codes"]
    assert 0 < r["score_with_frequency"] < 1


def test_overlap_note_gives_whole_points_and_is_translated():
    import numpy as np

    from missensexplorer.gui import i18n

    note = acmg.overlap_note({"acmg_points": np.float64(4.0), "acmg_ps1_pm5": "PM5 (p.R1C: Pathogenic, 3★)"})
    assert note.startswith("PP3 (4 points) and PM5 may partly reflect")
    i18n.set_language("es")
    try:
        assert i18n.tr_text(note).startswith("PP3 (4 puntos) y PM5")
    finally:
        i18n.set_language("en")


def test_a_common_variant_that_clinvar_calls_pathogenic_carries_a_caution():
    from missensexplorer.gui import i18n

    hbs = {"clinvar": "Pathogenic", "clinvar_stars": 2, "population_frequency": "frequent: BS1 candidate (grpmax 5.7%)"}
    note = acmg.common_pathogenic_note(hbs)
    assert note.startswith("ClinVar: Pathogenic (2★): a common pathogenic variant; review BA1 and BS1")
    assert acmg.common_pathogenic_note({**hbs, "clinvar_stars": 1}) is None  # not reviewed by anyone else
    assert acmg.common_pathogenic_note({**hbs, "clinvar": "Conflicting classifications of pathogenicity"}) is None
    assert acmg.common_pathogenic_note({**hbs, "clinvar": "Uncertain significance"}) is None
    assert acmg.common_pathogenic_note({**hbs, "clinvar": "Pathogenic (submissions)"}) is None  # VEP's, not ClinVar's
    for cls in (  # ClinVar's qualifiers (SERPINA1 PiZ and PiS, GBA1 p.Asn409Ser, RYR1 malignant hyperthermia)
        "Pathogenic/Likely pathogenic; risk factor",
        "Pathogenic/Pathogenic, low penetrance; other",
        "Pathogenic; drug response",
        "Likely pathogenic, low penetrance",
    ):
        assert acmg.common_pathogenic_note({**hbs, "clinvar": cls}).startswith(f"ClinVar: {cls} (2★): a common"), cls
    rare = {**hbs, "population_frequency": "rare (grpmax 0.0002%)"}
    assert acmg.common_pathogenic_note(rare) is None
    at_ba1 = acmg.common_pathogenic_note(rare, ba1_reached=True)  # a panel's BA1 reached on grpmax (SERPINC1)
    assert "a pathogenic variant at or above an expert panel's BA1 threshold; review" in at_ba1
    assert "common" not in at_ba1
    i18n.set_language("es")
    try:
        assert i18n.tr_text(note).startswith("ClinVar: patogénica (2★): una variante patogénica común")
        assert i18n.tr_text(at_ba1).startswith("ClinVar: patogénica (2★): una variante patogénica en el umbral")
        piz = acmg.common_pathogenic_note({**hbs, "clinvar": "Pathogenic/Likely pathogenic; risk factor"})
        assert i18n.tr_text(piz).startswith(
            "ClinVar: patogénica/probablemente patogénica; factor de riesgo (2★): una variante patogénica común"
        )
    finally:
        i18n.set_language("en")


def _pipeline(monkeypatch, frequency, panel_notes, classification="Pathogenic; drug response"):
    """BRAF V600E through the pipeline with its own ClinVar record, frequency label and panel notes replaced."""
    from test_robustness import BRAF

    class Vep:
        def __init__(self, *a, **k):
            pass

        def annotate(self, queries, params, progress=None, stage="", source="Ensembl VEP"):
            return {(kind, q): {"results": BRAF} for kind, qs in queries.items() for q in set(qs)}

    own = {**RECORDS["VCV000013961"], "classification": classification,
           "review_status": "criteria provided, multiple submitters, no conflicts"}  # fmt: skip
    records = {**RECORDS, "VCV000013961": own}
    monkeypatch.setattr(P, "ensembl_release", lambda http, cache, asm, **k: 116)
    monkeypatch.setattr(P, "_match_rate", lambda *a, **k: (1.0, {}))
    monkeypatch.setattr(P, "Vep", Vep)
    monkeypatch.setattr(
        P.cv, "fetch", lambda http, cache, ids, status=None: {k: records[k] for k in ids if k in records}
    )
    monkeypatch.setattr(P.cv, "search_windows", lambda http, cache, ws, asm, status=None: {w: [] for w in ws})
    monkeypatch.setattr(P.acmg, "frequency", lambda r: frequency)
    monkeypatch.setattr(P.vcep, "notes", lambda r: list(panel_notes))
    res = P.predict(
        text="7 140753336 A T", options=P.PredictOptions(cache=False, protein_context=False, sequences=False)
    )
    return str(res.table.iloc[0]["acmg_notes"] or "")


def test_pipeline_adds_the_common_pathogenic_caution(monkeypatch):
    af = "X VCEP (BRAF): grpmax 0.2% is at or above its BA1 threshold (> 0.1%): check its other BA1 conditions"
    faf = (
        "X VCEP (BRAF): grpmax allele frequency 0.0002% is above its BA1 threshold (> 0.0001%), which it applies "
        "to the filtering allele frequency (not computed here): check the filtering allele frequency before BA1"
    )
    rule = "X VCEP (BRAF): its own BA1 rule: “MAF ≥ 0.5% for autosomal recessive”"
    notes = _pipeline(monkeypatch, "frequent: BS1 candidate (grpmax 1.2%)", [])
    assert notes.startswith("ClinVar: Pathogenic; drug response (2★): a common pathogenic variant")
    notes = _pipeline(monkeypatch, "rare (grpmax 0.2%)", [af])
    assert notes.startswith("ClinVar: Pathogenic; drug response (2★): a pathogenic variant at or above an expert")
    for other in (faf, rule):  # neither says that the panel's BA1 is reached
        assert _pipeline(monkeypatch, "rare (grpmax 0.2%)", [other]) == other
    assert _pipeline(monkeypatch, "frequent: BS1 candidate (grpmax 1.2%)", [], "Benign") == ""
