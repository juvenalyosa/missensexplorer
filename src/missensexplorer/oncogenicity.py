"""Somatic oncogenicity points (ClinGen/CGC/VICC, Horak et al., Genet Med 2022) from open data.

The standard scores how likely a somatic variant is to drive a tumour. Each criterion adds points
(very strong 8, strong 4, moderate 2, supporting 1; benign ones subtract), and the total gives the
class: Oncogenic >= 10, Likely oncogenic 6 to 9, VUS 0 to 5, Likely benign -1 to -6, Benign <= -7.

Criteria applied automatically, and where their data come from:

OVS1  null variant (stop gained, frameshift, canonical splice site, start lost) in a gene IntOGen
      calls loss-of-function, unless nonsense-mediated decay is likely escaped (last exon).
OS1   same amino-acid change as a ClinVar Oncogenic variant made by another nucleotide change
      (OS1_Moderate, 2 points, when it is only Likely oncogenic).
OS3   cancerhotspots.org hotspot with >= 50 tumours at the residue and >= 10 with this change;
OM3   fewer than 50 at the residue and >= 10 with this change; OP3 fewer than 10 with this change.
      Single-residue hotspots for missense variants, in-frame indel hotspots for in-frame indels.
OM2   in-frame deletion or insertion in a known cancer gene (IntOGen), or stop lost in a
      loss-of-function gene.
OM4   another missense change at the residue is ClinVar Oncogenic (OM4_Supporting, 1 point, when
      only Likely oncogenic), and this change is at least about as different (Grantham distance at
      least 80% of the reference change's).
OP1   the MissensExplorer model calls the missense variant damaging; SBP1 (-1) benign. Not applied in
      genes IntOGen calls activating: there, a damaging call does not mean activation (in EGFR it
      pointed the wrong way, see docs/model.md).
OP4   absent from gnomAD or below an allele frequency of 1e-5 (the standard asks for "absent or
      extremely low" without a number);
SBS1  above 1% in a gnomAD continental population (grpmax); SBVS1 above 5%.
OS2   functional evidence for this protein change curated in CIViC (rated >= 3 of 5): CIViC
      "Oncogenic" evidence supporting oncogenicity, or "Functional" evidence of an effect that drives
      tumours in that kind of gene (gain of function or neomorphic in any driver; also loss of function
      or dominant negative in a loss-of-function gene). SBS2 (-4): CIViC evidence of unaltered
      function, or Oncogenic evidence that does not support oncogenicity. Neither when the evidence
      conflicts, and not with OS1 (the studies concern the amino-acid change, not one nucleotide change).

Not assessed (they need a curator): functional studies beyond CIViC, OM1 critical functional domain,
OP2 single-aetiology tumours, SBP2 synonymous variants. The variant's own ClinVar classification
never counts as evidence for itself: it is shown separately (``clinvar_oncogenicity``).

Clinical evidence (not points): :func:`clinical` summarises CIViC's predictive, diagnostic and
prognostic evidence for the protein change (the AMP/ASCO/CAP axis: what the variant means for
treatment), with the best evidence level (A validated or guideline ... E inferential).

Exclusions follow the standard: OS1 removes the hotspot criteria (its count mixes nucleotide
changes), OS1 and OS3 remove OM4, and OM4 and OM3 are not combined (the stronger one is kept; OM4
wins a tie, as the standard has it).

Data: cancerhotspots.org (Chang et al., Cancer Discov 2018), the IntOGen driver compendium
(Martinez-Jimenez et al., Nat Rev Cancer 2020; CC0) and CIViC's nightly evidence (Griffith et al.,
Nat Genet 2017; CC0), downloaded on first use and cached; ClinVar oncogenicity classifications through
NCBI E-utilities.
"""

from __future__ import annotations

import csv
import io
import re
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any

from ._types import KeyValueCache, Row
from .aa import grantham, to_one_letter
from .api import ApiError, Http, OfflineError, SourceStatus

HOTSPOTS_URL = "https://www.cancerhotspots.org/api/hotspots/single"
INTOGEN_RELEASE = "IntOGen-Drivers-20240920"
INTOGEN_URL = f"https://www.intogen.org/download?file={INTOGEN_RELEASE}.zip"
DATA_MAX_AGE = 90 * 86400  # both resources are released rarely
CIVIC_URL = "https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv"
CIVIC_MAX_AGE = 7 * 86400  # CIViC is updated nightly
CIVIC_MIN_RATING = 3  # of 5 stars: "well-established" functional studies
LEVELS = "ABCDE"
OP4_MAX_AF = 1e-5
SBS1_AF, SBVS1_AF = 0.01, 0.05
OM4_GRANTHAM_SHARE = 0.8

POINTS = {
    "OVS1": 8,
    "OS1": 4,
    "OS1_Moderate": 2,
    "OS2": 4,
    "OS3": 4,
    "OM2": 2,
    "OM3": 2,
    "OM4": 2,
    "OM4_Supporting": 1,
    "OP1": 1,
    "OP3": 1,
    "OP4": 1,
    "SBVS1": -8,
    "SBS1": -4,
    "SBS2": -4,
    "SBP1": -1,
}
NULL_TERMS = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost"}
INFRAME_TERMS = {"inframe_deletion", "inframe_insertion", "protein_altering_variant"}
ONCOGENIC = {"Oncogenic": "O", "Oncogenic/Likely oncogenic": "O", "Likely oncogenic": "LO"}
REVIEWED = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
}
ROLE_NAMES = {"Act": "activating", "LoF": "loss-of-function", "ambiguous": "ambiguous"}
TITLE = re.compile(r"^(N[MR]_\d+)(?:\.\d+)?\([^)]*\):(c\.\S+) \(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)$")
NOT_ASSESSED = "not assessed automatically: functional studies beyond CIViC, OM1 functional domain, OP2"
ONCOGENIC_FUNCTION = {  # CIViC functional significance that supports an oncogenic effect, by gene role
    "Act": {"Gain of Function", "Neomorphic"},
    "LoF": {"Gain of Function", "Neomorphic", "Loss of Function", "Dominant Negative"},
    "ambiguous": {"Gain of Function", "Neomorphic"},
}
CHANGE = re.compile(r"^[A-Z*]\d+(?:_[A-Z*]\d+)?(?:[A-Z*]|del[A-Z]*|ins[A-Z]+|delins[A-Z]+|dup[A-Z]*|fs\S*)?$")


def tier(points: int) -> str:
    """Return the class of a point total (Horak et al. 2022)."""
    if points >= 10:
        return "Oncogenic"
    if points >= 6:
        return "Likely oncogenic"
    if points >= 0:
        return "VUS"
    return "Likely benign" if points >= -6 else "Benign"


# --------------------------------------------------------------------------- data


@dataclass
class Role:
    """A gene's role in IntOGen: the role most of its cohorts give, and how many do."""

    role: str  # Act | LoF | ambiguous
    cohorts: int
    total: int

    def text(self) -> str:
        """'activating in 22 of 22 IntOGen cohorts'."""
        return f"{ROLE_NAMES.get(self.role, self.role)} in {self.cohorts} of {self.total} IntOGen cohorts"


@dataclass
class Hotspot:
    """A cancerhotspots.org hotspot: tumours with a variant at the residue (or in the indel range), by change."""

    tumours: int
    changes: dict[str, int]
    transcript: str


@dataclass
class Data:
    """Hotspots by gene and the IntOGen role by gene; empty when they could not be obtained."""

    residues: dict[str, dict[str, Hotspot]] = field(default_factory=dict)  # gene -> 'L858' -> hotspot
    indels: dict[str, list[tuple[int, int, Hotspot]]] = field(default_factory=dict)  # gene -> [(start, end, h)]
    roles: dict[str, Role] = field(default_factory=dict)
    civic: dict[str, dict[str, list[dict[str, Any]]]] = field(default_factory=dict)  # gene -> change -> items
    hotspots_ok: bool = False
    roles_ok: bool = False
    civic_ok: bool = False


def parse_hotspots(items: list[dict[str, Any]]) -> tuple[dict[str, dict[str, Hotspot]], dict[str, list]]:
    """Index cancerhotspots.org's single-residue and in-frame indel hotspots by gene."""
    residues: dict[str, dict[str, Hotspot]] = defaultdict(dict)
    indels: dict[str, list[tuple[int, int, Hotspot]]] = defaultdict(list)
    for h in items:
        gene, kind = h.get("hugoSymbol"), h.get("type")
        hs = Hotspot(
            int(h.get("tumorCount") or 0), dict(h.get("variantAminoAcid") or {}), str(h.get("transcriptId") or "")
        )
        if not gene:
            continue
        if kind == "single residue":
            residues[gene][str(h.get("residue"))] = hs
        elif kind == "in-frame indel" and re.fullmatch(r"\d+-\d+", str(h.get("residue"))):
            a, b = (int(x) for x in str(h["residue"]).split("-"))
            indels[gene].append((a, b, hs))
    return dict(residues), dict(indels)


def parse_intogen(zip_bytes: bytes) -> dict[str, Role]:
    """Return each gene's majority role across the cohorts of IntOGen's Compendium_Cancer_Genes.tsv."""
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as z:
        name = next(n for n in z.namelist() if n.endswith("Compendium_Cancer_Genes.tsv"))
        text = z.read(name).decode("utf-8")
    by_gene: dict[str, Counter[str]] = defaultdict(Counter)
    for rec in csv.DictReader(io.StringIO(text), delimiter="\t"):
        if rec.get("SYMBOL") and rec.get("ROLE"):
            by_gene[rec["SYMBOL"]][rec["ROLE"]] += 1
    out = {}
    for gene, c in by_gene.items():
        (role, n), total = c.most_common(1)[0], sum(c.values())
        out[gene] = Role(role, n, total)
    return out


def parse_civic(text: str) -> dict[str, dict[str, list[dict[str, Any]]]]:
    """Index CIViC's accepted evidence items by gene and protein change ('EGFR' -> 'L858R' -> [items]).

    Only molecular profiles of one variant written as a protein change are kept ('EGFR L858R',
    'VHL L178Q (c.533T>A)'); combinations, fusions, expression and 'Mutation' profiles are left out.
    """
    out: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for rec in csv.DictReader(io.StringIO(text), delimiter="\t"):
        if rec.get("evidence_status") != "accepted":
            continue
        parts = str(rec.get("molecular_profile") or "").split()
        if len(parts) < 2 or not CHANGE.match(parts[1]) or (len(parts) > 2 and not parts[2].startswith(("(c.", "c."))):
            continue
        rating = _num(rec.get("rating"))
        out[parts[0]][parts[1]].append(
            {
                "type": rec.get("evidence_type") or "",
                "direction": rec.get("evidence_direction") or "",
                "level": rec.get("evidence_level") or "",
                "significance": rec.get("significance") or "",
                "therapies": rec.get("therapies") or "",
                "disease": rec.get("disease") or "",
                "rating": int(rating) if rating is not None else 0,
                "pmid": rec.get("citation_id") if rec.get("source_type") == "PubMed" else None,
                "eid": rec.get("evidence_id") or "",
                "origin": rec.get("variant_origin") or "",
            }
        )
    return {g: dict(v) for g, v in out.items()}


def load(http: Http, cache: KeyValueCache, status: SourceStatus | None = None) -> Data:
    """Return the hotspots, gene roles and CIViC evidence, from the cache or downloaded.

    The hotspots and roles are kept 90 days, CIViC a week. A failed download is recorded in `status`
    (sources 'cancerhotspots.org', 'IntOGen' and 'CIViC') and leaves that part of the data empty.
    """
    data = Data()
    hot = cache.get("oncogenicity|cancerhotspots-single", max_age=DATA_MAX_AGE)
    if hot is None:
        try:
            r = http.request("GET", HOTSPOTS_URL)
            if r.status_code != 200:
                raise ApiError(f"HTTP {r.status_code}", r.status_code)
            hot = r.json()
            cache.set("oncogenicity|cancerhotspots-single", hot)
        except (ApiError, OfflineError, ValueError) as e:
            if status is not None:
                status.fail("cancerhotspots.org", str(e))
            hot = None
    if hot is not None:
        if status is not None:
            status.use("cancerhotspots.org")
        data.residues, data.indels = parse_hotspots(hot)
        data.hotspots_ok = True
    roles = cache.get(f"oncogenicity|{INTOGEN_RELEASE}", max_age=DATA_MAX_AGE)
    if roles is None:
        try:
            r = http.request("GET", INTOGEN_URL, follow_redirects=True)
            if r.status_code != 200:
                raise ApiError(f"HTTP {r.status_code}", r.status_code)
            roles = {g: [x.role, x.cohorts, x.total] for g, x in parse_intogen(r.content).items()}
            cache.set(f"oncogenicity|{INTOGEN_RELEASE}", roles)
        except (ApiError, OfflineError, ValueError, KeyError, StopIteration, zipfile.BadZipFile) as e:
            if status is not None:
                status.fail("IntOGen", str(e))
            roles = None
    if roles is not None:
        if status is not None:
            status.use("IntOGen")
        data.roles = {g: Role(*v) for g, v in roles.items()}
        data.roles_ok = True
    civic = cache.get("oncogenicity|civic-evidence", max_age=CIVIC_MAX_AGE)
    if civic is None:
        try:
            r = http.request("GET", CIVIC_URL, follow_redirects=True)
            if r.status_code != 200:
                raise ApiError(f"HTTP {r.status_code}", r.status_code)
            civic = parse_civic(r.text)
            cache.set("oncogenicity|civic-evidence", civic)
        except (ApiError, OfflineError, ValueError, KeyError) as e:
            if status is not None:
                status.fail("CIViC", str(e))
            civic = None
    if civic is not None:
        if status is not None:
            status.use("CIViC")
        data.civic = civic
        data.civic_ok = True
    return data


# --------------------------------------------------------------------------- ClinVar references


@dataclass
class Reference:
    """A ClinVar Oncogenic or Likely oncogenic missense variant at the row's residue (not the row itself)."""

    accession: str
    cdna: str
    protein: str
    alt_aa: str
    level: str  # O | LO
    same_change: bool


def references(row: Row, records: dict[str, dict[str, Any]], assembly: str) -> list[Reference]:
    """Return the ClinVar oncogenic missense variants at the row's residue, with review criteria (>= 1 star).

    Matched like ACMG PS1/PM5 (see residue.py): on the MANE Select transcript named in ClinVar's title,
    the protein position and the reference amino acid; the variant itself never counts.
    """
    from .residue import _same_allele

    nm = str(row.get("mane_select") or "").split(".")[0]
    pos, ref_aa, alt_aa = row.get("protein_position"), row.get("ref_aa"), row.get("alt_aa")
    own = (str(row.get("chrom")), int(row.get("pos") or 0), str(row.get("ref")), str(row.get("alt")))
    own_cdna = str(row.get("hgvsc") or "").split(":")[-1]
    out = []
    for acc, rec in records.items():
        level = ONCOGENIC.get(rec.get("oncogenicity") or "")
        m = TITLE.match(rec.get("title") or "")
        if not level or not REVIEWED.get(rec.get("oncogenicity_review") or "") or not m:
            continue
        if m.group(1) != nm or int(m.group(4)) != pos:
            continue
        r_aa, a_aa = to_one_letter(m.group(3)), to_one_letter(m.group(5))
        if r_aa != ref_aa or not a_aa or a_aa in (r_aa, "*"):
            continue
        if m.group(2) == own_cdna or _same_allele(rec, own, assembly):
            continue
        out.append(Reference(acc, m.group(2), f"p.{m.group(3)}{m.group(4)}{m.group(5)}", a_aa, level, a_aa == alt_aa))
    return sorted(out, key=lambda x: (not x.same_change, x.level != "O", x.accession))


# --------------------------------------------------------------------------- evidence


@dataclass
class Assessment:
    """The criteria met by one row, with their points, and the notes on criteria left out."""

    codes: list[tuple[str, str]] = field(default_factory=list)  # (code, why)
    notes: list[str] = field(default_factory=list)

    def ordered(self) -> list[tuple[str, str]]:
        """The criteria from the strongest oncogenic to the benign ones."""
        return sorted(self.codes, key=lambda x: (POINTS[x[0]] < 0, -abs(POINTS[x[0]])))

    @property
    def points(self) -> int:
        """Sum of the points of the criteria met."""
        return sum(POINTS[c] for c, _ in self.codes)

    def summary(self) -> str | None:
        """'OS3 +4, OP4 +1' (None when no criterion is met)."""
        return ", ".join(f"{c} {POINTS[c]:+d}" for c, _ in self.ordered()) or None

    def evidence(self) -> str:
        """Every criterion with its reason, then the notes."""
        return "; ".join([f"{c} ({why})" for c, why in self.ordered()] + self.notes)


def assess(
    row: Row,
    terms: set[str],
    model_call: str | None,
    data: Data,
    refs: list[Reference] | None,
    gnomad_info: dict[str, Any] | None = None,
) -> Assessment:
    """Return the oncogenicity criteria a row meets.

    Args:
        row: One results row; reads ``gene``, ``transcript``, ``canonical``, ``protein_position``,
            ``ref_aa``, ``alt_aa``, ``protein_change``, ``nmd``, ``gnomad_af`` and ``gnomad_popmax_af``.
        terms: The row's VEP consequence terms.
        model_call: ``"damaging"``, ``"benign"`` or None: the model's call before the frequency rule.
        data: Hotspots and gene roles, from :func:`load`.
        refs: ClinVar oncogenic variants at the residue (:func:`references`); None when not searched.
        gnomad_info: gnomAD's own answer for a variant Ensembl gave no frequency (``found``, ``af``).
    """
    a = Assessment()
    gene = str(row.get("gene") or "")
    role = data.roles.get(gene)
    if data.roles_ok:
        a.notes.append(f"{gene}: {role.text()}" if role else f"{gene} is not an IntOGen driver gene")
    missense = "missense_variant" in terms
    if terms & NULL_TERMS:
        _null(a, row, role)
    elif missense:
        _missense(a, row, model_call, data, refs, role)
    elif "stop_lost" in terms or (terms & INFRAME_TERMS and re.search(r"del|ins|dup", str(row.get("protein_change")))):
        _inframe(a, row, terms, data, role)
    _functional(a, row, data, role)
    _frequency(a, row, gnomad_info)
    a.notes.append(NOT_ASSESSED)
    return a


def _null(a: Assessment, row: Row, role: Role | None) -> None:
    """OVS1: a null variant in a loss-of-function gene, unless decay is likely escaped."""
    if not role or role.role != "LoF":
        a.notes.append("OVS1 not applied (the gene is not a loss-of-function driver in IntOGen)")
    elif "escape" in str(row.get("nmd") or ""):
        a.notes.append(f"OVS1 not applied ({row['nmd']}: the truncated protein may keep its function)")
    else:
        a.codes.append(("OVS1", f"null variant in a {ROLE_NAMES['LoF']} gene"))


def _missense(
    a: Assessment, row: Row, model_call: str | None, data: Data, refs: list[Reference] | None, role: Role | None
) -> None:
    """OS1/OM4 (ClinVar), OS3/OM3/OP3 (hotspots) and OP1/SBP1 (the model) for a missense variant."""
    ref_aa, alt_aa, pos = row.get("ref_aa"), row.get("alt_aa"), row.get("protein_position")
    same = [x for x in refs or [] if x.same_change]
    other = [x for x in refs or [] if not x.same_change]
    if refs is None:
        a.notes.append("OS1/OM4 not checked (ClinVar oncogenicity search unavailable)")
    if same:
        x = same[0]
        code = "OS1" if x.level == "O" else "OS1_Moderate"
        a.codes.append((code, f"same change by {x.cdna}: {_level(x)}, {x.accession}"))
    hot = _residue_hotspot(row, data)
    if hot and not same:
        n_aa = hot.changes.get(str(alt_aa), 0)
        code = "OS3" if hot.tumours >= 50 and n_aa >= 10 else "OM3" if n_aa >= 10 else "OP3"
        a.codes.append(
            (code, f"hotspot {ref_aa}{pos}: {hot.tumours} tumours, {n_aa} with {alt_aa}; cancerhotspots.org")
        )
    elif hot and same:
        a.notes.append("hotspot not counted (OS1 already uses this amino-acid change)")
    has = {c for c, _ in a.codes}
    if other and not has & {"OS1", "OS1_Moderate", "OS3"}:
        g_new = grantham(str(ref_aa), str(alt_aa)) or 0
        ok = [x for x in other if g_new >= OM4_GRANTHAM_SHARE * (grantham(str(ref_aa), x.alt_aa) or 0)]
        if ok:
            x = sorted(ok, key=lambda x: x.level != "O")[0]
            code = "OM4" if x.level == "O" else "OM4_Supporting"
            if "OM3" in has and POINTS[code] < POINTS["OM3"]:  # the two cannot be combined: keep the stronger
                a.notes.append(f"{code} not used ({x.protein}: {_level(x)}; the hotspot criterion OM3 is stronger)")
            else:
                a.codes = [(c, w) for c, w in a.codes if c != "OM3"]
                a.codes.append((code, f"{x.protein}: {_level(x)}, {x.accession}"))
        else:
            a.notes.append("OM4 not applied (this change is milder than the oncogenic ones at the residue)")
    if role and role.role == "Act":
        a.notes.append(
            "OP1/SBP1 not applied (activating gene: the model measures damage to the protein, not activation)"
        )
    elif model_call == "damaging":
        a.codes.append(("OP1", "the MissensExplorer model calls it damaging"))
    elif model_call == "benign":
        a.codes.append(("SBP1", "the MissensExplorer model calls it benign"))


def _inframe(a: Assessment, row: Row, terms: set[str], data: Data, role: Role | None) -> None:
    """OM2 and the in-frame indel hotspots."""
    gene = str(row.get("gene") or "")
    change = str(row.get("protein_change") or "")
    if "stop_lost" in terms:
        if role and role.role == "LoF":
            a.codes.append(("OM2", "stop lost in a loss-of-function gene"))
        return
    if role:
        a.codes.append(("OM2", f"in-frame change in a cancer gene ({role.text()})"))
    else:
        a.notes.append("OM2 not applied (not a known cancer gene in IntOGen)")
    m = re.match(r"^[A-Z*](\d+)(?:_[A-Z*](\d+))?", change)
    if not m:
        return
    start, end = int(m.group(1)), int(m.group(2) or m.group(1))
    for lo, hi, h in data.indels.get(gene, []):
        if lo <= start and end <= hi + 1:
            n = h.changes.get(change, 0)
            code = "OS3" if h.tumours >= 50 and n >= 10 else "OM3" if n >= 10 else "OP3"
            a.codes.append((code, f"in-frame indel hotspot {lo}-{hi}: {h.tumours} tumours, {n} with {change}"))
            break


def _civic_items(row: Row, data: Data) -> list[dict[str, Any]]:
    """CIViC evidence items for the row's gene and protein change."""
    return data.civic.get(str(row.get("gene") or ""), {}).get(str(row.get("protein_change") or ""), [])


def _functional(a: Assessment, row: Row, data: Data, role: Role | None) -> None:
    """OS2 / SBS2 from CIViC's Oncogenic and Functional evidence for this protein change."""
    if not data.civic_ok:
        a.notes.append("OS2/SBS2 not checked (CIViC unavailable)")
        return
    items = [x for x in _civic_items(row, data) if x["type"] in ("Oncogenic", "Functional")]
    if not items:
        return
    good = [x for x in items if x["rating"] >= CIVIC_MIN_RATING]
    if not good:
        a.notes.append(f"CIViC functional evidence rated below {CIVIC_MIN_RATING} stars not used")
        return
    onc_sig = ONCOGENIC_FUNCTION.get(role.role if role else "ambiguous", ONCOGENIC_FUNCTION["ambiguous"])
    pro: list[dict[str, Any]] = []
    con: list[dict[str, Any]] = []
    for x in good:
        supports = x["direction"] == "Supports"
        if x["type"] == "Oncogenic":
            (pro if supports else con).append(x)
        elif supports and x["significance"] in onc_sig:
            pro.append(x)
        elif supports and x["significance"] == "Unaltered Function":
            con.append(x)
    if pro and con:
        a.notes.append("OS2/SBS2 not applied (CIViC functional evidence conflicts)")
    elif pro and any(c.startswith("OS1") for c, _ in a.codes):
        a.notes.append("OS2 not applied (OS1 already uses this amino-acid change)")
    elif pro:
        a.codes.append(("OS2", "CIViC: " + _civic_refs(pro)))
    elif con:
        a.codes.append(("SBS2", "CIViC: " + _civic_refs(con)))


def _civic_refs(items: list[dict[str, Any]]) -> str:
    """'Gain of Function, 3 stars, EID1234 (PMID 123)' for the best-rated items (at most two)."""
    best = sorted(items, key=lambda x: -x["rating"])
    out = []
    for x in best[:2]:
        what = x["significance"] if x["type"] == "Functional" else f"oncogenicity {x['direction'].lower()}"
        pmid = f" (PMID {x['pmid']})" if x.get("pmid") else ""
        out.append(f"{what}, {x['rating']} stars, EID{x['eid']}{pmid}")
    more = f"; and {len(items) - 2} more" if len(items) > 2 else ""
    return "; ".join(out) + more


def clinical(row: Row, data: Data) -> tuple[str | None, str | None]:
    """CIViC's clinical evidence (predictive, diagnostic, prognostic) for the row's protein change.

    Returns:
        The best evidence level (A-E) and a summary: the number of items by level, then the best item
        (type, significance, therapies, disease, EID); (None, None) when CIViC has none.
    """
    items = [
        x
        for x in _civic_items(row, data)
        if x["type"] in ("Predictive", "Diagnostic", "Prognostic")
        and x["direction"] == "Supports"
        and x["origin"] not in ("Rare Germline", "Common Germline")
    ]
    if not items:
        return None, None
    items.sort(key=lambda x: (LEVELS.find(x["level"]) % 9, -x["rating"]))
    counts = Counter(x["level"] for x in items)
    by_level = ", ".join(f"{lv} {counts[lv]}" for lv in LEVELS if counts.get(lv))
    b = items[0]
    what = b["type"].lower() + (f", {b['significance']}" if b["significance"] else "")
    drugs = f": {b['therapies'].replace(',', ', ')}" if b["therapies"] else ""
    return b["level"], (
        f"{len(items)} CIViC items (level {by_level}); best: level {b['level']} {what}{drugs} in {b['disease']} "
        f"(EID{b['eid']})"
    )


def _frequency(a: Assessment, row: Row, gnomad_info: dict[str, Any] | None) -> None:
    """SBVS1/SBS1 from the gnomAD grpmax (continental populations); otherwise OP4 from the overall frequency."""
    grpmax = _num(row.get("gnomad_popmax_af"))
    af = _num(row.get("gnomad_af"))
    if af is None and gnomad_info and gnomad_info.get("found"):
        af = _num(gnomad_info.get("af"))
    if af is None:
        af = grpmax
    if grpmax is not None and grpmax > SBVS1_AF:
        a.codes.append(("SBVS1", f"gnomAD grpmax {grpmax:.3g} > 5%"))
    elif grpmax is not None and grpmax > SBS1_AF:
        a.codes.append(("SBS1", f"gnomAD grpmax {grpmax:.3g} > 1%"))
    elif not af:
        a.codes.append(("OP4", "absent from gnomAD"))
    elif af < OP4_MAX_AF:
        a.codes.append(("OP4", f"gnomAD frequency {af:.2g} < 1e-5"))


def _residue_hotspot(row: Row, data: Data) -> Hotspot | None:
    """The hotspot at the row's residue, when its transcript is the hotspot's (or the row's is canonical)."""
    hs = data.residues.get(str(row.get("gene") or ""), {})
    h = hs.get(f"{row.get('ref_aa')}{row.get('protein_position')}")
    if not h:
        return None
    tx = str(row.get("transcript") or "").split(".")[0]
    canonical = str(row.get("canonical")).lower() in ("true", "1", "yes")
    return h if tx == h.transcript or canonical or row.get("mane_select") else None


def _level(x: Reference) -> str:
    return "Oncogenic" if x.level == "O" else "Likely oncogenic"


def _num(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v
