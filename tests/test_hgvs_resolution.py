"""How typed HGVS is resolved: ambiguous protein changes, indel alignment, readable errors (offline).

tests/data/vep_ambiguous_kras_q61h.json is VEP's real answer to KRAS:p.Gln61His with ambiguous_hgvs,
trimmed to the fields the code reads (see tests/data/record_fixtures.py).
"""

import importlib
import json
from pathlib import Path
from types import SimpleNamespace

import httpx
import pandas as pd
import pytest

from missensexplorer import protein as prot
from missensexplorer import vep as V
from missensexplorer.cache import Cache, NullCache
from missensexplorer.inputs import region_query

P = importlib.import_module("missensexplorer.predict")  # the package re-exports a function with the same name
DATA = Path(__file__).parent / "data"
Q61H = json.loads((DATA / "vep_ambiguous_kras_q61h.json").read_text())
VEP_Q61H_ERROR = (
    "Unable to parse HGVS notation 'KRAS:p.Gln61His':  : Sequence translated from reference (GGG -> G) does not "
    "match input sequence (Q): Sequence translated from reference (ACC -> T) does not match input sequence (Q)"
)


# --------------------------------------------------------------------------- readable errors


@pytest.mark.parametrize(
    ("raw", "short"),
    [
        (VEP_Q61H_ERROR, V.AA_MISMATCH),
        (
            "Unable to parse HGVS notation 'TP53:p.Arg999His': Can't locate object method \"strand\" via package "
            '"Bio::EnsEMBL::Mapper::Gap" at /usr/share/perl5/Bio/EnsEMBL/TranscriptMapper.pm line 380.',
            V.OUTSIDE_PROTEIN,
        ),
        (
            "Unable to parse HGVS notation 'BRCA1:p.Glu23fs':  : Sequence translated from reference (TGG -> W) does "
            "not match input sequence (E)",
            V.PROTEIN_FRAMESHIFT,
        ),
        ("Could not get a Transcript object for 'BRAFF'", "unknown gene or transcript symbol"),
        (V.AA_MISMATCH, V.AA_MISMATCH),  # already short: unchanged
    ],
)
def test_vep_errors_are_readable(raw, short):
    assert V._short_error(raw) == short


# --------------------------------------------------------------------------- ambiguous protein changes


def test_protein_matches_keeps_each_allele_that_gives_the_change_once():
    hits = V.protein_matches("KRAS:p.Gln61His", Q61H)
    assert len(Q61H) == 30  # one result per transcript and allele, with alleles that change another residue
    alleles = sorted((r["seq_region_name"], r["start"], r["allele_string"]) for r in hits)
    assert alleles == [("12", 25227341, "A/C"), ("12", 25227341, "A/T")]


def test_protein_matches_needs_the_requested_position_on_the_reported_transcript():
    def result(start, aa, pos):
        tc = {"gene_symbol": "BRAF", "amino_acids": aa, "protein_start": pos}
        return {"seq_region_name": "7", "start": start, "end": start, "allele_string": "C/A", "strand": 1,
                "transcript_consequences": [tc]}  # fmt: skip

    # BRAF:p.Ala600Glu (residue 600 is a valine): an isoform has an alanine there, which is A557E on MANE
    assert V.protein_matches("BRAF:p.Ala600Glu", [result(140776936, "A/E", 557)]) == []
    both = [result(100, "G/D", 368), result(200, "G/D", 396), result(300, "S/T", 396)]
    assert V.protein_matches("BRAF:p.Gly396Asp", both) == [both[1]]
    assert V.protein_matches("BRAF:p.Phe508del", both) == []  # only substitutions are resolved this way


class FakeVepHttp:
    """Answers like the VEP REST API: the ambiguous KRAS query fails unless ambiguous_hgvs is set."""

    def __init__(self):
        self.calls = []

    def request(self, method, url, params=None, json=None, headers=None):
        batch = json["hgvs_notations"]
        self.calls.append((dict(params), list(batch)))
        req = httpx.Request(method, url)
        if batch == ["BRAFF:p.Val600Glu"]:
            return httpx.Response(400, json={"error": "Could not get a Transcript object for 'BRAFF'"}, request=req)
        if params.get("ambiguous_hgvs"):
            return httpx.Response(200, json=Q61H, request=req)
        return httpx.Response(400, json={"error": VEP_Q61H_ERROR}, request=req)


def test_annotate_asks_again_for_an_ambiguous_protein_change():
    http = FakeVepHttp()
    vep = V.Vep(http, NullCache(), "GRCh38", workers=1)
    out = vep.annotate({"hgvs": ["KRAS:p.Gln61His", "BRAFF:p.Val600Glu"]}, V.BASE_PARAMS)
    q61h = out[("hgvs", "KRAS:p.Gln61His")]
    assert q61h["ambiguous"] == 2
    assert len(q61h["results"]) == 2
    assert out[("hgvs", "BRAFF:p.Val600Glu")]["error"] == "unknown gene or transcript symbol"  # not asked again
    retried = [batch for params, batch in http.calls if params.get("ambiguous_hgvs")]
    assert retried == [["KRAS:p.Gln61His"]]


def test_an_error_cached_in_vep_wording_is_reworded_and_retried(tmp_path):
    http = FakeVepHttp()
    with Cache(tmp_path / "c.sqlite") as cache:
        vep = V.Vep(http, cache, "GRCh38", workers=1)
        cache.set(vep._key("hgvs", V.BASE_PARAMS, "KRAS:p.Gln61His"), {"error": VEP_Q61H_ERROR, "cache_error": True})
        out = vep.annotate({"hgvs": ["KRAS:p.Gln61His"]}, V.BASE_PARAMS)
    assert out[("hgvs", "KRAS:p.Gln61His")]["ambiguous"] == 2
    assert [p.get("ambiguous_hgvs") for p, _ in http.calls] == [1]  # the cached error itself is not asked again


def test_a_protein_change_reported_at_another_position_gets_a_note():
    from missensexplorer.inputs import parse_text

    v = parse_text("MUTYH G396D").variants[0]
    row = {"gene": "MUTYH", "protein_position": 368, "protein_change": "G368D", "ref_aa": "G", "alt_aa": "D"}
    assert P._renumbered(v, row) == (
        "protein change given as p.Gly396Asp (another isoform's numbering); on the reported transcript it is G368D"
    )
    assert P._renumbered(v, {**row, "protein_position": 396}) is None
    assert P._renumbered(v, {**row, "gene": "TOE1"}) is None  # an overlapping gene
    assert P._renumbered(parse_text("chr1:45331556 C>T").variants[0], row) is None  # not a protein change
    # TP53 p.Ala175His: residue 175 is an arginine; Ensembl finds an isoform with an alanine there, whose
    # SNV is C182S on the reported transcript, so the amino acids differ too
    wrong = {"gene": "TP53", "protein_position": 182, "protein_change": "C182S", "ref_aa": "C", "alt_aa": "S"}
    assert P._renumbered(parse_text("TP53 A175H").variants[0], wrong) == (
        "protein change given as p.Ala175His exists only on another isoform; on the reported transcript this DNA "
        "change is C182S: check the input"
    )


# --------------------------------------------------------------------------- indel alignment

#            1234567890
GENOME = (1, "CCATCTTTGG")  # like CFTR around F508: ...ATC TTT GG...


@pytest.mark.parametrize(
    ("start", "deleted", "inserted", "vcf"),
    [
        (5, "CTT", "-", ("7", 3, "ATCT", "A")),  # c.1521_1523del: shifts left to the classic F508del
        (4, "TCT", "-", ("7", 3, "ATCT", "A")),  # already left-aligned
        (6, "TTT", "-", ("7", 5, "CTTT", "C")),  # another deletion, which cannot shift past the C
        (9, "-", "T", ("7", 5, "C", "CT")),  # T inserted after the run of Ts: moves to its start
        (4, "-", "G", ("7", 3, "A", "AG")),  # nothing to shift
    ],
)
def test_left_align(start, deleted, inserted, vcf):
    assert P.left_align("7", start, deleted, inserted, GENOME) == vcf


def test_left_align_refuses_a_deletion_the_genome_does_not_hold():
    assert P.left_align("7", 5, "GGG", "-", GENOME) is None


def test_indel_location_reads_vep_alleles_on_the_forward_strand():
    minus = {"seq_region_name": "17", "start": 43124028, "allele_string": "AG/-", "strand": -1}
    assert P._indel_location(minus) == ("17", 43124028, "CT", "-")
    assert P._indel_location({"seq_region_name": "7", "start": 5, "allele_string": "A/T"}) is None  # an SNV
    assert P._indel_location({"seq_region_name": "7", "start": 5, "allele_string": "AT/GC"}) is None  # delins
    assert P._indel_location({"seq_region_name": "GL000220.1", "start": 5, "allele_string": "A/-"}) is None


def test_hgvs_indels_are_annotated_in_their_left_aligned_form(monkeypatch):
    hgvs = ("hgvs", "CFTR:c.1521_1523del")
    answers = {hgvs: {"results": [{"input": hgvs[1], "seq_region_name": "7", "start": 5, "allele_string": "CTT/-"}]}}
    norm = region_query("7", 3, "ATCT", "A")
    asked = {}

    def fake_annotate(ctx, queries):
        asked.update(queries)
        return {("region", norm): {"results": [{"input": norm, "most_severe_consequence": "inframe_deletion"}]}}, {
            ("region", norm)
        }

    monkeypatch.setattr(P, "_annotate", fake_annotate)
    monkeypatch.setattr(P, "_genome_windows", lambda ctx, locs: dict.fromkeys(locs, GENOME))
    no_predictors = set()
    P._normalise_indels(SimpleNamespace(net=True), answers, no_predictors)
    assert asked == {"region": [norm]}
    assert answers[hgvs]["results"] == [{"input": hgvs[1], "most_severe_consequence": "inframe_deletion"}]
    assert no_predictors == {hgvs}  # the region answer had no predictor scores: neither has the input


def test_indel_alignment_needs_the_network():
    answers = {("hgvs", "X"): {"results": [{"seq_region_name": "7", "start": 5, "allele_string": "CTT/-"}]}}
    before = repr(answers)
    P._normalise_indels(SimpleNamespace(net=False), answers, set())
    assert repr(answers) == before


# --------------------------------------------------------------------------- GRCh37: UniProt entry names, PS1/PM5


def test_uniprot_entry_names_are_turned_into_accessions(tmp_path):
    seen = []

    def handler(request):
        seen.append(request.url.params["query"])
        return httpx.Response(
            200, json={"results": [{"primaryAccession": "P15056", "uniProtkbId": "BRAF_HUMAN"}]}, request=request
        )

    from missensexplorer.api import Http

    http = Http()
    http.client = httpx.Client(transport=httpx.MockTransport(handler))
    with Cache(tmp_path / "c.sqlite") as cache:
        ids = {"BRAF_HUMAN", "P04637", "NOPE_HUMAN", None}
        assert prot.uniprot_accessions(http, cache, ids) == {"BRAF_HUMAN": "P15056"}
        assert prot.uniprot_accessions(http, cache, {"BRAF_HUMAN"}) == {"BRAF_HUMAN": "P15056"}  # cached
    assert seen == ["id:BRAF_HUMAN OR id:NOPE_HUMAN"]  # accessions are not looked up


def test_same_residue_is_not_searched_without_a_mane_transcript(monkeypatch):
    searched = []
    monkeypatch.setattr(
        P.cv, "search_windows", lambda http, cache, windows, assembly, status: searched.extend(windows) or {}
    )
    monkeypatch.setattr(P.cv, "fetch", lambda *a, **k: {})
    snv = {"consequence": "missense_variant", "protein_position": 600, "chrom": "7", "ref": "A", "alt": "T"}
    rows = [{**snv, "pos": 140753336, "mane_select": "NM_004333.6"}, {**snv, "pos": 140453136, "mane_select": None}]
    ctx = SimpleNamespace(opt=SimpleNamespace(clinvar=True), net=True, http=None, cache=None, assembly="GRCh38",
                          status=None, say=lambda *a: None)  # fmt: skip
    windows, refs = P._search_same_residue(ctx, rows)
    assert set(windows) == {0, 1}
    assert searched == [windows[0]]
    found = P._Lookups()
    found.windows, found.residue_refs = windows, refs
    assert P._ps1_pm5(ctx, 1, rows[1], found).startswith("not checked (no MANE Select transcript")


# --------------------------------------------------------------------------- records dropped by FILTER


def test_terminal_summary_counts_records_that_failed_filter():
    from rich.console import Console

    res = P.PredictionResult(
        table=pd.DataFrame(columns=P.COLUMNS),
        errors=pd.DataFrame(columns=["input", "source", "error"]),
        assembly="GRCh38",
        stats={"inputs": 4, "rows": 0, "file_skipped_filter": 1, "file_skipped_not_carried": 2},
    )
    con = Console(record=True, width=200)
    from missensexplorer.report import print_result

    print_result(res, con)
    assert "Skipped file records: 2 not carried by the sample(s) · 1 not passing FILTER." in con.export_text()
