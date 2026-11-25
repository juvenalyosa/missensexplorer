"""Somatic oncogenicity points (ClinGen/CGC/VICC), offline: criteria, exclusions, data loading, pipeline."""

import importlib
import io
import json
import zipfile
from pathlib import Path

import httpx
import pytest

from missensexplorer import clinvar as cv
from missensexplorer import oncogenicity as onco
from missensexplorer.api import Http, SourceStatus
from missensexplorer.cache import Cache, NullCache

P = importlib.import_module("missensexplorer.predict")
DATA = Path(__file__).parent / "data"

HOTSPOTS = [
    {"hugoSymbol": "TP53", "residue": "R175", "type": "single residue", "tumorCount": 416,
     "variantAminoAcid": {"H": 386, "G": 9}, "transcriptId": "ENST00000269305"},
    {"hugoSymbol": "EGFR", "residue": "L858", "type": "single residue", "tumorCount": 144,
     "variantAminoAcid": {"R": 144}, "transcriptId": "ENST00000275493"},
    {"hugoSymbol": "EGFR", "residue": "G719", "type": "single residue", "tumorCount": 33,
     "variantAminoAcid": {"A": 18, "C": 7}, "transcriptId": "ENST00000275493"},
    {"hugoSymbol": "EGFR", "residue": "745-759", "type": "in-frame indel", "tumorCount": 155,
     "variantAminoAcid": {"E746_A750del": 123}, "transcriptId": "ENST00000275493"},
]  # fmt: skip
ROLES = {"TP53": onco.Role("LoF", 161, 212), "EGFR": onco.Role("Act", 22, 22)}


def data(**kw):
    residues, indels = onco.parse_hotspots(HOTSPOTS)
    return onco.Data(residues=residues, indels=indels, roles=ROLES, hotspots_ok=True, roles_ok=True, **kw)


def row(gene="TP53", change="R175H", pos=175, ref="R", alt="H", tx="ENST00000269305", **kw):
    return {
        "gene": gene,
        "protein_change": change,
        "protein_position": pos,
        "ref_aa": ref,
        "alt_aa": alt,
        "transcript": f"{tx}.9",
        "canonical": True,
        "mane_select": "NM_000546.6",
        "gnomad_af": None,
        "gnomad_popmax_af": None,
        **kw,
    }


def codes(a):
    return [c for c, _ in a.ordered()]


def ref(protein, alt_aa, level="O", same=False, cdna="c.1A>G"):
    return onco.Reference("VCV000000001", cdna, protein, alt_aa, level, same)


# --------------------------------------------------------------------------- classes and criteria


@pytest.mark.parametrize(
    ("points", "tier"),
    [(10, "Oncogenic"), (6, "Likely oncogenic"), (5, "VUS"), (0, "VUS"), (-1, "Likely benign"), (-7, "Benign")],
)
def test_point_totals_map_to_the_five_classes(points, tier):
    assert onco.tier(points) == tier


def test_tumour_suppressor_missense_hotspot_with_damaging_model_call():
    a = onco.assess(row(), {"missense_variant"}, "damaging", data(), [])
    assert codes(a) == ["OS3", "OP1", "OP4"]  # 416 tumours at R175, 386 with H; absent from gnomAD
    assert a.points == 6
    assert onco.tier(a.points) == "Likely oncogenic"
    assert "loss-of-function in 161 of 212 IntOGen cohorts" in a.evidence()


def test_in_an_activating_gene_the_model_call_does_not_count():
    r = row("EGFR", "L858R", 858, "L", "R", "ENST00000275493")
    a = onco.assess(r, {"missense_variant"}, "damaging", data(), [])
    assert codes(a) == ["OS3", "OP4"]
    assert "OP1/SBP1 not applied (activating gene" in a.evidence()
    b = onco.assess(r, {"missense_variant"}, "benign", data(), [])
    assert "SBP1" not in codes(b)


def test_hotspot_counts_set_its_strength():
    weak = row("TP53", "R175G", 175, "R", "G")  # 9 tumours with G
    assert "OP3" in codes(onco.assess(weak, {"missense_variant"}, None, data(), []))
    moderate = row("EGFR", "G719A", 719, "G", "A", "ENST00000275493")  # < 50 at residue, 18 with A
    assert "OM3" in codes(onco.assess(moderate, {"missense_variant"}, None, data(), []))


def test_a_hotspot_on_another_transcript_is_not_used():
    r = row(tx="ENST00000999999", canonical=False, mane_select=None)
    assert not {"OS3", "OM3", "OP3"} & set(codes(onco.assess(r, {"missense_variant"}, None, data(), [])))


def test_same_change_in_clinvar_gives_os1_and_removes_the_hotspot():
    refs = [ref("p.Arg175His", "H", same=True, cdna="c.524_525delinsAT")]
    a = onco.assess(row(), {"missense_variant"}, "damaging", data(), refs)
    assert codes(a) == ["OS1", "OP1", "OP4"]
    assert "hotspot not counted" in a.evidence()
    refs = [ref("p.Arg175His", "H", level="LO", same=True)]
    assert codes(onco.assess(row(), {"missense_variant"}, None, data(), refs))[0] == "OS1_Moderate"


def test_om4_needs_a_change_at_least_as_different_and_replaces_om3():
    g719c = row("EGFR", "G719C", 719, "G", "C", "ENST00000275493")
    a = onco.assess(g719c, {"missense_variant"}, None, data(), [ref("p.Gly719Ala", "A")])
    assert codes(a) == ["OM4", "OP3", "OP4"]  # Grantham G>C 159 >= 0.8 x G>A 60
    g719a = row("EGFR", "G719A", 719, "G", "A", "ENST00000275493")
    b = onco.assess(g719a, {"missense_variant"}, None, data(), [ref("p.Gly719Cys", "C")])
    assert codes(b) == ["OM3", "OP4"]  # G>A (60) is milder than G>C (159): no OM4, the hotspot stays
    assert "OM4 not applied" in b.evidence()
    c = onco.assess(g719c, {"missense_variant"}, None, data(), [ref("p.Gly719Ala", "A", level="LO")])
    assert "OM4_Supporting" in codes(c)
    d = onco.assess(g719a, {"missense_variant"}, None, data(), [ref("p.Gly719Ser", "S", level="LO")])
    assert codes(d) == ["OM3", "OP4"]  # OM4_Supporting (1) would be weaker than the hotspot's OM3 (2)
    assert "OM4_Supporting not used" in d.evidence()


def test_clinvar_not_searched_is_said():
    a = onco.assess(row(), {"missense_variant"}, None, data(), None)
    assert "OS1/OM4 not checked" in a.evidence()


def test_null_variant_in_a_loss_of_function_gene_unless_decay_is_escaped():
    r = row(change="R213*", pos=213, alt="*", nmd="NMD likely (premature stop upstream of last exon)")
    assert codes(onco.assess(r, {"stop_gained"}, None, data(), None)) == ["OVS1", "OP4"]
    last = {**r, "nmd": "NMD escape likely (last exon)"}
    a = onco.assess(last, {"stop_gained"}, None, data(), None)
    assert "OVS1" not in codes(a)
    assert "OVS1 not applied (NMD escape likely (last exon)" in a.evidence()
    egfr = row("EGFR", "W731*", 731, "W", "*")
    assert "OVS1" not in codes(onco.assess(egfr, {"stop_gained"}, None, data(), None))


def test_inframe_deletion_in_a_cancer_gene_with_an_indel_hotspot():
    r = row("EGFR", "E746_A750del", 746, "ELREA", None, "ENST00000275493")
    a = onco.assess(r, {"inframe_deletion"}, None, data(), None)
    assert codes(a) == ["OS3", "OM2", "OP4"]  # 155 tumours in 745-759, 123 with this deletion
    assert a.points == 7
    other = onco.assess({**r, "gene": "ABC1"}, {"inframe_deletion"}, None, data(), None)
    assert codes(other) == ["OP4"]
    assert "OM2 not applied" in other.evidence()


@pytest.mark.parametrize(
    ("af", "grpmax", "info", "code"),
    [
        (None, None, None, "OP4"),
        (9e-7, 2e-5, None, "OP4"),  # overall frequency below 1e-5
        (3e-5, 4e-5, None, None),
        (0.004, 0.02, None, "SBS1"),
        (0.3, 0.75, None, "SBVS1"),
        (None, None, {"found": True, "af": 2e-4}, None),  # gnomAD's own answer counts
        (None, None, {"found": False}, "OP4"),
    ],
)
def test_frequency_criteria(af, grpmax, info, code):
    a = onco.assess(row(gnomad_af=af, gnomad_popmax_af=grpmax), {"missense_variant"}, None, data(), [], info)
    got = set(codes(a)) & {"OP4", "SBS1", "SBVS1"}
    assert got == ({code} if code else set())


def test_codes_summary_is_strongest_first():
    a = onco.Assessment(codes=[("OP4", ""), ("SBP1", ""), ("OS3", ""), ("OM2", "")])
    assert a.summary() == "OS3 +4, OM2 +2, OP4 +1, SBP1 -1"
    assert a.points == 6


# --------------------------------------------------------------------------- CIViC (OS2/SBS2, clinical evidence)

CIVIC_COLUMNS = [
    "molecular_profile", "evidence_type", "evidence_direction", "evidence_level", "significance", "rating",
    "therapies", "disease", "citation_id", "source_type", "evidence_id", "variant_origin", "evidence_status",
]  # fmt: skip
CIVIC_A = {
    "type": "Predictive", "direction": "Supports", "level": "A", "significance": "Sensitivity/Response",
    "therapies": "Vemurafenib", "disease": "Skin Melanoma", "rating": 4, "pmid": "1", "eid": "1409",
    "origin": "Somatic",
}  # fmt: skip


def civic_tsv(*rows):
    lines = ["\t".join(CIVIC_COLUMNS)]
    for mp, typ, direction, level, sig, rating, *rest in rows:
        status = rest[0] if rest else "accepted"
        lines.append(
            "\t".join(
                [mp, typ, direction, level, sig, str(rating), "", "Cancer", "123", "PubMed", "9", "Somatic", status]
            )
        )
    return "\n".join(lines) + "\n"


def item(typ="Functional", sig="Gain of Function", direction="Supports", rating=3, level="D"):
    return {"type": typ, "direction": direction, "level": level, "significance": sig, "therapies": "", "disease": "X",
            "rating": rating, "pmid": "123", "eid": "9", "origin": "Somatic"}  # fmt: skip


def with_civic(gene, change, *items):
    return data(civic={gene: {change: list(items)}}, civic_ok=True)


def test_civic_index_keeps_single_protein_changes_only():
    idx = onco.parse_civic(
        civic_tsv(
            ["EGFR L858R", "Predictive", "Supports", "A", "Sensitivity/Response", 4],
            ["VHL L178Q (c.533T>A)", "Functional", "Supports", "D", "Loss of Function", 3],
            ["EGFR L858R AND EGFR T790M", "Predictive", "Supports", "B", "Resistance", 3],
            ["EML4::ALK Fusion", "Predictive", "Supports", "A", "Sensitivity/Response", 5],
            ["TP53 Mutation", "Prognostic", "Supports", "B", "Poor Outcome", 3],
            ["KRAS G12C", "Predictive", "Supports", "A", "Sensitivity/Response", 5, "submitted"],
        )
    )
    assert set(idx) == {"EGFR", "VHL"}
    assert set(idx["VHL"]) == {"L178Q"}
    assert idx["EGFR"]["L858R"][0]["level"] == "A"


def test_os2_from_civic_follows_the_gene_role():
    tp53 = row()
    a = onco.assess(tp53, {"missense_variant"}, None, with_civic("TP53", "R175H", item(sig="Dominant Negative")), [])
    assert set(codes(a)[:2]) == {"OS2", "OS3"}
    assert "CIViC: Dominant Negative, 3 stars, EID9 (PMID 123)" in a.evidence()
    egfr = row("EGFR", "L858R", 858, "L", "R", "ENST00000275493")
    lof = with_civic("EGFR", "L858R", item(sig="Loss of Function"))
    assert "OS2" not in codes(onco.assess(egfr, {"missense_variant"}, None, lof, []))  # not oncogenic in an oncogene
    gof = with_civic("EGFR", "L858R", item(sig="Gain of Function"))
    assert "OS2" in codes(onco.assess(egfr, {"missense_variant"}, None, gof, []))
    onc = with_civic("EGFR", "L858R", item(typ="Oncogenic", sig="Oncogenicity"))
    assert "OS2" in codes(onco.assess(egfr, {"missense_variant"}, None, onc, []))


def test_sbs2_conflicts_low_ratings_and_os1():
    r = row()
    unaltered = with_civic("TP53", "R175H", item(sig="Unaltered Function"))
    assert "SBS2" in codes(onco.assess(r, {"missense_variant"}, None, unaltered, []))
    against = with_civic("TP53", "R175H", item(typ="Oncogenic", sig="Oncogenicity", direction="Does Not Support"))
    assert "SBS2" in codes(onco.assess(r, {"missense_variant"}, None, against, []))
    both = with_civic("TP53", "R175H", item(sig="Loss of Function"), item(sig="Unaltered Function"))
    a = onco.assess(r, {"missense_variant"}, None, both, [])
    assert not {"OS2", "SBS2"} & set(codes(a))
    assert "CIViC functional evidence conflicts" in a.evidence()
    weak = with_civic("TP53", "R175H", item(sig="Loss of Function", rating=2))
    assert "rated below 3 stars not used" in onco.assess(r, {"missense_variant"}, None, weak, []).evidence()
    refs = [ref("p.Arg175His", "H", same=True)]
    b = onco.assess(r, {"missense_variant"}, None, with_civic("TP53", "R175H", item(sig="Loss of Function")), refs)
    assert "OS2" not in codes(b)
    assert "OS2 not applied (OS1 already uses" in b.evidence()
    assert (
        "OS2/SBS2 not checked (CIViC unavailable)" in onco.assess(r, {"missense_variant"}, None, data(), []).evidence()
    )


def test_clinical_evidence_summary_takes_the_best_level():
    d = with_civic(
        "BRAF",
        "V600E",
        item(typ="Prognostic", sig="Poor Outcome", level="B"),
        CIVIC_A,
        item(typ="Predictive", sig="Resistance", level="C"),
        {**CIVIC_A, "origin": "Rare Germline"},
        {**CIVIC_A, "direction": "Does Not Support"},
        item(),  # functional evidence is not clinical evidence
    )
    level, text = onco.clinical({"gene": "BRAF", "protein_change": "V600E"}, d)
    assert level == "A"
    assert text.startswith("3 CIViC items (level A 1, B 1, C 1); best: level A predictive, Sensitivity/Response")
    assert onco.clinical({"gene": "BRAF", "protein_change": "V600K"}, d) == (None, None)


# --------------------------------------------------------------------------- ClinVar references


def test_references_keep_reviewed_oncogenic_missense_at_the_residue_but_not_the_variant_itself():
    def rec(title, onc="Oncogenic", review="criteria provided, single submitter"):
        return {"title": title, "oncogenicity": onc, "oncogenicity_review": review}

    records = {
        "VCV1": rec("NM_000546.6(TP53):c.524G>A (p.Arg175His)"),  # the variant itself
        "VCV2": rec("NM_000546.6(TP53):c.523C>G (p.Arg175Gly)", "Likely oncogenic"),
        "VCV3": rec("NM_000546.6(TP53):c.523C>T (p.Arg175Cys)", review="no assertion criteria provided"),
        "VCV4": rec("NM_000546.6(TP53):c.523C>A (p.Arg175Ser)", "Uncertain significance"),
        "VCV5": rec("NM_000546.6(TP53):c.527G>A (p.Cys176Tyr)"),
    }
    r = {**row(), "chrom": "17", "pos": 7675088, "ref": "C", "alt": "T", "hgvsc": "ENST00000269305.9:c.524G>A"}
    refs = onco.references(r, records, "GRCh38")
    assert [(x.accession, x.level, x.same_change) for x in refs] == [("VCV2", "LO", False)]


def test_clinvar_records_keep_the_oncogenicity_review_status():
    rec = cv._compact(
        {
            "oncogenicity_classification": {
                "description": "Oncogenic",
                "review_status": "criteria provided, single submitter",
            }
        }
    )
    assert rec["oncogenicity"] == "Oncogenic"
    assert rec["oncogenicity_review"] == "criteria provided, single submitter"


def test_oncogenicity_window_search_uses_its_own_filter_and_cache(tmp_path, monkeypatch, request):
    queries = []

    def handler(request):
        queries.append(request.content.decode())
        return httpx.Response(200, json={"esearchresult": {"idlist": ["12374"]}})

    http = Http()
    http.client = httpx.Client(transport=httpx.MockTransport(handler))
    rec = {"loc": {"GRCh38": ["17", 7675088, 7675088]}}
    monkeypatch.setattr(cv, "fetch", lambda http, cache, ids, status=None: dict.fromkeys(ids, rec))
    cache = Cache(tmp_path / "c.sqlite")
    request.addfinalizer(cache.close)
    w = ("17", 7675086, 7675090)
    assert cv.search_oncogenicity_windows(http, cache, [w], "GRCh38") == {w: ["VCV000012374"]}
    assert "oncogenicity_classification" in queries[0]
    assert "clinsig_pathogenic" not in queries[0]
    assert cv.search_oncogenicity_windows(http, cache, [w], "GRCh38") == {w: ["VCV000012374"]}
    assert len(queries) == 1  # cached
    assert cache.get(f"clinvar-window|GRCh38|17:{w[1]}-{w[2]}") is None  # apart from the PS1/PM5 searches


# --------------------------------------------------------------------------- data


def intogen_zip(rows):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        text = "SYMBOL\tCOHORT\tROLE\n" + "".join(f"{g}\tC{i}\t{r}\n" for i, (g, r) in enumerate(rows))
        z.writestr("2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", text)
    return buf.getvalue()


def test_intogen_roles_are_the_majority_across_cohorts():
    roles = onco.parse_intogen(intogen_zip([("TP53", "LoF"), ("TP53", "LoF"), ("TP53", "Act"), ("EGFR", "Act")]))
    assert roles["TP53"] == onco.Role("LoF", 2, 3)
    assert roles["EGFR"].text() == "activating in 1 of 1 IntOGen cohorts"


def test_load_downloads_once_then_uses_the_cache(tmp_path, request):
    calls = []

    def handler(request):
        calls.append(str(request.url))
        if "cancerhotspots" in str(request.url):
            return httpx.Response(200, json=HOTSPOTS)
        if "civicdb" in str(request.url):
            return httpx.Response(
                200, text=civic_tsv(["TP53 R175H", "Functional", "Supports", "D", "Dominant Negative", 3])
            )
        return httpx.Response(200, content=intogen_zip([("TP53", "LoF")]))

    http = Http()
    http.client = httpx.Client(transport=httpx.MockTransport(handler))
    cache = Cache(tmp_path / "c.sqlite")
    request.addfinalizer(cache.close)
    st = SourceStatus()
    d = onco.load(http, cache, st)
    assert d.hotspots_ok
    assert d.roles_ok
    assert d.residues["TP53"]["R175"].tumours == 416
    assert d.roles["TP53"].role == "LoF"
    assert d.civic["TP53"]["R175H"][0]["significance"] == "Dominant Negative"
    assert st.summary() == {"CIViC": "ok", "IntOGen": "ok", "cancerhotspots.org": "ok"}
    onco.load(http, cache, SourceStatus())
    assert len(calls) == 3


def test_load_records_a_failed_download():
    http = Http()
    http.client = httpx.Client(transport=httpx.MockTransport(lambda request: httpx.Response(404)))
    st = SourceStatus()
    d = onco.load(http, NullCache(), st)
    assert not d.hotspots_ok
    assert not d.roles_ok
    assert st.failed("cancerhotspots.org")
    assert st.failed("IntOGen")
    assert st.failed("CIViC")
    assert not d.civic_ok


# --------------------------------------------------------------------------- pipeline


@pytest.fixture
def braf_pipeline(monkeypatch):
    braf = json.loads((DATA / "vep_braf_v600e_region.json").read_text())

    class Vep:
        def __init__(self, *a, **k):
            pass

        def annotate(self, queries, params, progress=None, stage="", source="Ensembl VEP"):
            return {(kind, q): {"results": braf} for kind, qs in queries.items() for q in set(qs)}

    monkeypatch.setattr(P, "ensembl_release", lambda http, cache, asm: 116)
    monkeypatch.setattr(P, "_match_rate", lambda *a, **k: (1.0, {}))
    monkeypatch.setattr(P, "Vep", Vep)
    monkeypatch.setattr(P.cv, "fetch", lambda http, cache, ids, status=None: {})
    monkeypatch.setattr(P.cv, "search_windows", lambda http, cache, ws, asm, status=None: {w: [] for w in ws})
    monkeypatch.setattr(P.gnomad, "lookup", lambda *a, **k: {})
    searched = []

    def onc_search(http, cache, ws, asm, status=None):
        searched.extend(ws)
        return {w: [] for w in ws}

    monkeypatch.setattr(P.cv, "search_oncogenicity_windows", onc_search)
    hot = [
        {"hugoSymbol": "BRAF", "residue": "V600", "type": "single residue", "tumorCount": 897,
         "variantAminoAcid": {"E": 833}, "transcriptId": "ENST00000288602"},
    ]  # fmt: skip
    residues, indels = onco.parse_hotspots(hot)
    civic = {"BRAF": {"V600E": [CIVIC_A]}}
    loaded = onco.Data(
        residues=residues,
        indels=indels,
        roles={"BRAF": onco.Role("Act", 54, 57)},
        civic=civic,
        hotspots_ok=True,
        roles_ok=True,
        civic_ok=True,
    )
    monkeypatch.setattr(P.onco, "load", lambda http, cache, status=None: loaded)
    return searched


def test_pipeline_adds_the_oncogenicity_columns_only_when_asked(braf_pipeline):
    opts = {"cache": False, "protein_context": False, "sequences": False}
    res = P.predict(text="7 140753336 A T", options=P.PredictOptions(oncogenicity=True, **opts))
    r = res.table.iloc[0]
    assert r["oncogenicity_codes"].startswith("OS3 +4")
    assert "hotspot V600: 897 tumours, 833 with E" in r["oncogenicity_evidence"]
    assert "OP1/SBP1 not applied" in r["oncogenicity_evidence"]
    assert r["oncogenicity"] in ("VUS", "Likely oncogenic")
    assert r["clinical_evidence_level"] == "A"
    assert "best: level A predictive, Sensitivity/Response: Vemurafenib in Skin Melanoma" in r["clinical_evidence"]
    assert braf_pipeline  # the residue was searched in ClinVar
    plain = P.predict(text="7 140753336 A T", options=P.PredictOptions(**opts)).table.iloc[0]
    assert P.pd.isna(plain["oncogenicity"])


def test_pipeline_says_when_the_oncogenicity_data_are_missing(braf_pipeline, monkeypatch):
    monkeypatch.setattr(P.onco, "load", lambda http, cache, status=None: onco.Data())
    res = P.predict(
        text="7 140753336 A T",
        options=P.PredictOptions(oncogenicity=True, cache=False, protein_context=False, sequences=False),
    )
    assert res.table.iloc[0]["oncogenicity"] == "not assessed (cancerhotspots.org or IntOGen unavailable)"
