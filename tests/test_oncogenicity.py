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
ROLES = {"TP53": onco.Role("LoF", "UniProt", 0.96, 212), "EGFR": onco.Role("Act", "UniProt", 0.0, 22)}


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
        "gnomad_af": 2e-6,  # rare in gnomAD: OP4 by frequency
        "gnomad_popmax_af": 4e-6,
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
    assert "loss-of-function (UniProt: tumour suppressor)" in a.evidence()


def test_the_computational_evidence_counts_in_every_gene_with_a_note_in_activating_genes():
    r = row("EGFR", "L858R", 858, "L", "R", "ENST00000275493")
    a = onco.assess(r, {"missense_variant"}, "damaging", data(), [])
    assert codes(a) == ["OS3", "OP1", "OP4"]  # 6 points: likely oncogenic, as the standard has it
    assert "does not tell activation from inactivation" in a.evidence()
    b = onco.assess(r, {"missense_variant"}, "benign", data(), [])
    assert "SBP1" in codes(b)


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
    assert codes(a) == ["OS1", "OP1", "OP4"]  # OS1 (4) against OS3 (4): the standard keeps OS1
    assert "OS3 not counted (OS1 uses this amino-acid change)" in a.evidence()
    # a likely oncogenic reference (2 points) never replaces a stronger hotspot (4 points)
    lo = [ref("p.Arg175His", "H", level="LO", same=True)]
    b = onco.assess(row(), {"missense_variant"}, "damaging", data(), lo)
    assert codes(b) == ["OS3", "OP1", "OP4"]
    assert "OS1_Moderate not counted" in b.evidence()
    # without a hotspot, the likely oncogenic reference gives OS1_Moderate
    c = onco.assess(row(pos=176, change="C176F", ref="C", alt="F"), {"missense_variant"}, None, data(),
                    [ref("p.Cys176Phe", "F", level="LO", same=True)])  # fmt: skip
    assert codes(c)[0] == "OS1_Moderate"


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


def test_ovs1_is_graded_like_pvs1_by_decay_and_the_share_of_protein_removed():
    r = row(change="R213*", pos=213, alt="*", nmd="NMD likely (premature stop upstream of last exon)")
    assert codes(onco.assess(r, {"stop_gained"}, None, data(), None)) == ["OVS1", "OP4"]
    last = {**r, "nmd": "NMD escape likely (last exon)", "protein_position": 300, "protein_length": 393}
    a = onco.assess(last, {"stop_gained"}, None, data(), None)
    assert codes(a) == ["OVS1_Strong", "OP4"]  # 24% of the protein removed
    assert "removes 24% of the protein" in a.evidence()
    tail = {**last, "protein_position": 380}
    assert codes(onco.assess(tail, {"stop_gained"}, None, data(), None))[0] == "OVS1_Moderate"  # 4%: extreme 3' end
    unknown = {**last, "protein_length": None}
    assert codes(onco.assess(unknown, {"stop_gained"}, None, data(), None))[0] == "OVS1_Moderate"
    structure = {"exons": [], "cds": [1, 1], "length": 393}
    assert codes(
        onco.assess({**last, "protein_length": None}, {"stop_gained"}, None, data(), None, structure=structure)
    )[0] == ("OVS1_Strong")
    start = row(change="M1?", pos=1, alt=None)
    assert codes(onco.assess(start, {"start_lost"}, None, data(), None))[0] == "OVS1_Moderate"
    egfr = row("EGFR", "W731*", 731, "W", "*")
    b = onco.assess(egfr, {"stop_gained"}, None, data(), None)
    assert not {"OVS1", "OVS1_Strong", "OVS1_Moderate"} & set(codes(b))
    assert "OVS1 not applied (not a loss-of-function gene" in b.evidence()


def test_splice_site_ovs1_depends_on_whether_skipping_the_exon_keeps_the_frame():
    # exons in transcript order (plus strand); the CDS runs from 105 to 400
    structure = {"exons": [[1, 110], [201, 260], [301, 400]], "cds": [105, 400], "length": 100}
    donor = row(change=None, pos=None, alt=None, intron="2/2")  # donor of intron 2: exon 2 (60 coding bases)
    a = onco.assess(donor, {"splice_donor_variant"}, None, data(), None, structure=structure)
    assert codes(a)[0] == "OVS1_Strong"  # in frame (60 % 3 == 0), 60 of 166 coding bases
    acceptor = row(change=None, pos=None, alt=None, intron="1/2")  # acceptor of intron 1: exon 2 too
    assert (
        codes(onco.assess(acceptor, {"splice_acceptor_variant"}, None, data(), None, structure=structure))[0]
        == "OVS1_Strong"
    )
    shifted = {"exons": [[1, 110], [201, 261], [301, 400]], "cds": [105, 400], "length": 100}
    assert codes(onco.assess(donor, {"splice_donor_variant"}, None, data(), None, structure=shifted))[0] == "OVS1"
    assert codes(onco.assess(donor, {"splice_donor_variant"}, None, data(), None, structure=None))[0] == "OVS1_Strong"


def test_inframe_deletion_in_a_cancer_gene_with_an_indel_hotspot():
    r = row("EGFR", "E746_A750del", 746, "ELREA", None, "ENST00000275493")
    a = onco.assess(r, {"inframe_deletion"}, None, data(), None)
    assert codes(a) == ["OS3", "OM2", "OP4"]  # 155 tumours in 745-759, 123 with this deletion
    assert a.points == 7
    other = onco.assess({**r, "gene": "ABC1"}, {"inframe_deletion"}, None, data(), None)
    assert codes(other) == ["OP4"]
    assert "OM2 not applied (not a known cancer gene" in other.evidence()


@pytest.mark.parametrize(
    ("af", "grpmax", "info", "code"),
    [
        (None, None, None, None),  # absence not confirmed by gnomAD itself: no OP4
        (9e-7, 2e-5, None, "OP4"),  # overall frequency below 1e-5
        (3e-5, 4e-5, None, None),
        (0.004, 0.02, None, "SBS1"),
        (0.3, 0.75, None, "SBVS1"),
        (None, None, {"found": True, "af": 2e-4}, None),  # gnomAD's own answer counts
        (None, None, {"found": False, "dataset": "gnomAD v4", "over_20": {"exome": 0.95}}, "OP4"),
        (None, None, {"found": False, "dataset": "gnomAD v4", "over_20": {"exome": 0.1}}, None),  # few reads
    ],
)
def test_frequency_criteria(af, grpmax, info, code):
    a = onco.assess(row(gnomad_af=af, gnomad_popmax_af=grpmax), {"missense_variant"}, None, data(), [], info)
    got = set(codes(a)) & {"OP4", "SBS1", "SBVS1"}
    assert got == ({code} if code else set())
    if af is None and info is None:
        assert "OP4 not checked (absence from gnomAD not confirmed)" in a.evidence()


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
    # a same-change ClinVar reference never lowers the points: the hotspot + CIViC path (OS3 + OS2) is stronger
    refs = [ref("p.Arg175His", "H", same=True)]
    b = onco.assess(r, {"missense_variant"}, None, with_civic("TP53", "R175H", item(sig="Loss of Function")), refs)
    assert codes(b) == ["OS3", "OS2", "OP4"]
    assert "OS1 not counted" in b.evidence()
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
        rows = [(g, r, *rest) for g, r, *rest in rows]
        text = "SYMBOL\tCOHORT\tROLE\tSAMPLES\n" + "".join(
            f"{g}\tC{i}\t{r}\t{(rest or [10])[0]}.0\n" for i, (g, r, *rest) in enumerate(rows)
        )
        z.writestr("2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", text)
    return buf.getvalue()


def test_gene_roles_come_from_intogen_weighted_by_samples_then_uniprot_keywords():
    counts = onco.parse_intogen(
        intogen_zip([("SMARCA4", "Act", 10), ("SMARCA4", "Act", 10), ("SMARCA4", "LoF", 80), ("EGFR", "Act", 50),
                     ("MIXED", "LoF", 50), ("MIXED", "Act", 50), ("ETV6", "LoF", 80), ("ETV6", "LoF", 10),
                     ("ETV6", "Act", 10), ("THIN", "Act", 5), ("FUSION", "Act", 3), ("PMS2", "ambiguous", 26),
                     ("SMB1", "LoF", 5), ("SMB1", "Act", 4), ("SMB1", "Act", 4), ("AMBX", "LoF", 4),
                     ("AMBX", "ambiguous", 10), ("MAML", "LoF", 14), ("MAML", "ambiguous", 5),
                     ("MAML", "ambiguous", 6)])
    )  # fmt: skip
    assert counts["SMARCA4"] == [80.0, 20.0, 0.0, 3, 0]
    assert counts["MAML"] == [14.0, 0.0, 11.0, 1, 2]  # cohorts with a role, then cohorts of ambiguous role
    roles = onco.combine_roles(counts, suppressors={"MLH1", "PMS2", "SMB1"}, oncogenes={"ETV6", "FUSION", "MIXED"})
    assert roles["SMARCA4"].role == "LoF"  # most cohorts say Act, but 80% of the mutated samples say LoF
    assert (
        roles["SMARCA4"].text()
        == "loss-of-function (IntOGen: 80% loss-of-function in 100 samples of cohorts with a role, 3 cohorts)"
    )
    # samples in cohorts of ambiguous role do not count against loss of function (PMS2 was "activating")
    assert (roles["PMS2"].role, roles["PMS2"].source) == ("LoF", "UniProt")
    assert roles["AMBX"].role == "ambiguous"
    assert "in cohorts of ambiguous role" in roles["AMBX"].text()
    # a few IntOGen samples do not override UniProt's keyword (SMARCB1 was "activating" on 13 samples)
    assert roles["SMB1"].role == "ambiguous"
    assert "they disagree" in roles["SMB1"].text()
    assert roles["EGFR"].role == "Act"
    assert roles["MLH1"].role == "LoF"  # UniProt only
    # strong IntOGen evidence (>= 3 cohorts, >= 50 samples) wins over UniProt's Proto-oncogene keyword
    assert (roles["ETV6"].role, roles["ETV6"].source) == ("LoF", "IntOGen")
    # an ambiguous IntOGen role, or one from a single small cohort, gives way to UniProt's keyword
    assert (roles["MIXED"].role, roles["MIXED"].source) == ("Act", "UniProt")
    assert (roles["FUSION"].role, roles["FUSION"].source) == ("Act", "UniProt")
    # a single small cohort and no UniProt keyword: the role is left undecided
    assert roles["THIN"].role == "ambiguous"
    assert not roles["THIN"].decisive
    assert roles["THIN"].text() == "ambiguous (IntOGen: 5 mutated samples in 1 cohort with a role: too few to decide)"
    # cohorts of ambiguous role do not count towards the two cohorts IntOGen needs (MAML2: one cohort of 14)
    assert roles["MAML"].role == "ambiguous"
    assert not roles["MAML"].decisive
    assert onco.parse_uniprot("Entry\tGene Names (primary)\nP1\tTP53\nP2\t\n") == {"TP53"}


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
        if "uniprot" in str(request.url):
            gene = "TP53" if "KW-0043" in str(request.url) else "KRAS"
            return httpx.Response(200, text=f"Entry\tGene Names (primary)\nP1\t{gene}\n")
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
    assert d.roles["KRAS"].role == "Act"
    assert st.summary() == {"CIViC": "ok", "IntOGen": "ok", "UniProt": "ok", "cancerhotspots.org": "ok"}
    onco.load(http, cache, SourceStatus())
    assert len(calls) == 5


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
        roles={"BRAF": onco.Role("Act", "UniProt", 0.01, 57)},
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
    assert (
        "does not tell activation from inactivation" in r["oncogenicity_evidence"]
        or "OP1" not in r["oncogenicity_codes"]
    )
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


def test_hotspot_gene_symbols_are_updated_to_current_hgnc_symbols():
    residues, _ = onco.parse_hotspots(
        [{"hugoSymbol": "H3F3A", "residue": "K28", "type": "single residue", "tumorCount": 33,
          "variantAminoAcid": {"M": 33}, "transcriptId": "ENST00000366696"}]
    )  # fmt: skip
    assert "H3-3A" in residues
    assert "H3F3A" not in residues


def test_os2_needs_preclinical_civic_evidence():
    r = row()
    clinical = with_civic("TP53", "R175H", item(typ="Oncogenic", sig="Oncogenicity", level="B"))
    a = onco.assess(r, {"missense_variant"}, None, clinical, [])
    assert "OS2" not in codes(a)
    assert "from clinical observations (not functional studies) not used" in a.evidence()


def egfr(change, pos=719, ref_aa="G", alt="C"):
    return row("EGFR", change, pos, ref_aa, alt, "ENST00000275493")


def test_adding_a_clinvar_reference_never_lowers_the_points():
    """For every combination of references (same change or other change, oncogenic or likely oncogenic),
    with hotspots of each strength (OS3, OM3, OP3, none), OM4, CIViC evidence and model calls, adding a
    reference never lowers the points."""
    cases = [  # (row, CIViC significance, milder other changes at the residue)
        (row(), "Loss of Function", [("p.Arg175Lys", "K"), ("p.Arg175Lys", "K")]),  # OS3 (OM4 excluded)
        (row(change="R175G", alt="G"), "Loss of Function", [("p.Arg175His", "H"), ("p.Arg175Lys", "K")]),  # OP3
        (row(pos=176, change="R176H"), "Loss of Function", [("p.Arg176Lys", "K"), ("p.Arg176Lys", "K")]),  # none
        (egfr("G719A", alt="A"), "Gain of Function", [("p.Gly719Ser", "S"), ("p.Gly719Ser", "S")]),  # OM3, OM4
        (egfr("G719C"), "Gain of Function", [("p.Gly719Ala", "A"), ("p.Gly719Ser", "S")]),  # OP3 with OM4
    ]
    seen: set[str] = set()
    for r, sig, milder in cases:
        aa3 = {"H": "His", "G": "Gly", "A": "Ala", "C": "Cys"}[r["alt_aa"]]
        same = f"p.{'Arg' if r['ref_aa'] == 'R' else 'Gly'}{r['protein_position']}{aa3}"
        refs = {
            "same O": ref(same, r["alt_aa"], level="O", same=True, cdna="c.1_2delinsAT"),
            "same LO": ref(same, r["alt_aa"], level="LO", same=True, cdna="c.1_2delinsAC"),
            "other O": ref(milder[0][0], milder[0][1], level="O"),
            "other LO": ref(milder[1][0], milder[1][1], level="LO"),
        }
        civic = with_civic(r["gene"], r["protein_change"], item(sig=sig))
        for d in (data(), civic, onco.Data(roles=ROLES, roles_ok=True, hotspots_ok=True)):
            for call in ("damaging", None):
                names = list(refs)
                for k in range(1 << len(names)):
                    chosen = [refs[n] for i, n in enumerate(names) if k >> i & 1]
                    a = onco.assess(r, {"missense_variant"}, call, d, chosen)
                    seen |= set(codes(a))
                    for extra in names:
                        if refs[extra] in chosen:
                            continue
                        more = onco.assess(r, {"missense_variant"}, call, d, [*chosen, refs[extra]]).points
                        assert more >= a.points, (r["protein_change"], [x.level for x in chosen], extra)
    # the cases exercise every criterion that a same-change reference is weighed against
    assert {"OS1", "OS1_Moderate", "OS3", "OM3", "OP3", "OS2", "OM4", "OM4_Supporting"} <= seen


def test_os1_combines_with_om3_and_op3_as_the_standard_does():
    same = ref("p.Gly719Ala", "A", same=True)
    # an oncogenic same-change reference and an OM3 hotspot: OS1 + OM3 (the standard does not exclude OM3)
    a = onco.assess(egfr("G719A", alt="A"), {"missense_variant"}, "damaging", data(), [same])
    assert codes(a) == ["OS1", "OM3", "OP1", "OP4"]
    assert a.points == 8
    assert "standard's rule" not in a.evidence()
    # with CIViC evidence too: OS1 + OM3 (6) ties OS2 + OM3 (6); either way 8 points, as the standard gives
    b = onco.assess(egfr("G719A", alt="A"), {"missense_variant"}, "damaging",
                    with_civic("EGFR", "G719A", item(sig="Gain of Function")), [same])  # fmt: skip
    assert b.points == 8
    assert "standard's rule" not in b.evidence()
    # an OP3 hotspot stays with OS1
    c = onco.assess(egfr("G719C"), {"missense_variant"}, None, data(), [ref("p.Gly719Cys", "C", same=True)])
    assert codes(c) == ["OS1", "OP3", "OP4"]
    # OM4 had replaced OM3; when OS1 replaces OM4, OM3 comes back
    d = onco.assess(egfr("G719A", alt="A"), {"missense_variant"}, None, data(), [same, ref("p.Gly719Ser", "S")])
    assert codes(d) == ["OS1", "OM3", "OP4"]
    assert "OM4 not counted (OS1 uses this amino-acid change)" in d.evidence()


def test_the_standards_rule_for_os1_is_reported_when_it_differs():
    refs = [ref("p.Arg175His", "H", same=True)]
    a = onco.assess(row(), {"missense_variant"}, "damaging", with_civic("TP53", "R175H", item(sig="Loss of Function")),
                    refs)  # fmt: skip
    assert a.points == 10  # OS3 + OS2 (stronger than OS1) + OP1 + OP4
    assert "the standard's rule for OS1 would give 6 points" in a.evidence()
    # a likely oncogenic reference gives no OS1 under the standard: the hotspot evidence stands
    lo = onco.assess(row(), {"missense_variant"}, None, data(), [ref("p.Arg175His", "H", level="LO", same=True)])
    assert "standard's rule" not in lo.evidence()


def test_unavailable_uniprot_keywords_are_said():
    a = onco.assess(row(), {"missense_variant"}, None, data(roles_partial=True), [])
    assert "UniProt keywords unavailable: gene roles from IntOGen only" in a.evidence()


def test_civic_release_names_are_english_whatever_the_locale():
    import datetime as dt
    import locale

    old = locale.setlocale(locale.LC_TIME)
    try:
        for name in ("es_ES.UTF-8", "de_DE.UTF-8", "fr_FR.UTF-8", "es_ES", "de_DE", "fr_FR"):
            try:
                locale.setlocale(locale.LC_TIME, name)
            except locale.Error:
                continue
            if dt.date(2026, 3, 1).strftime("%b") != "Mar":  # the locale really names months otherwise
                break
        else:
            pytest.skip("no non-English locale is installed")
        assert onco.civic_releases(dt.date(2026, 1, 15)) == ["01-Jan-2026", "01-Dec-2025"]
        assert onco.civic_releases(dt.date(2026, 8, 3)) == ["01-Aug-2026", "01-Jul-2026"]
    finally:
        locale.setlocale(locale.LC_TIME, old)
