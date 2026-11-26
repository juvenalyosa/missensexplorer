"""Somatic oncogenicity points (ClinGen/CGC/VICC, Horak et al., Genet Med 2022) from open data.

The standard scores how likely a somatic variant is to drive a tumour. Each criterion adds points
(very strong 8, strong 4, moderate 2, supporting 1; benign ones subtract), and the total gives the
class: Oncogenic >= 10, Likely oncogenic 6 to 9, VUS 0 to 5, Likely benign -1 to -6, Benign <= -7.
The class here comes from the criteria that can be computed automatically, with no tumour type: it is
a starting point for a curator, not an oncogenicity classification.

Gene roles: IntOGen's role weighted by mutated samples (loss-of-function when >= 60% of the samples in
loss-of-function or activating cohorts are in loss-of-function ones, activating when <= 40%), when it
rests on enough samples (Role); otherwise, or when IntOGen is ambiguous, UniProt's curated keyword
(Tumor suppressor / Proto-oncogene) when it names one role. When the two disagree, IntOGen decides only
with strong evidence; otherwise the role is ambiguous.

Criteria applied automatically, and where their data come from:

OVS1  null variant (stop gained, frameshift, canonical splice site, start lost) in a loss-of-function
      gene, graded as ClinGen's loss-of-function guidance grades PVS1 (Abou Tayoun et al. 2018), which
      the standard recommends, using its option to downgrade a criterion: full when nonsense-mediated
      decay is expected or skipping the exon shifts the frame; OVS1_Strong (4) when decay is escaped
      (last exon, single-exon gene) and more than 10% of the protein is removed, or an in-frame skip
      removes more than 10% of the coding sequence; OVS1_Moderate (2) otherwise, and for a lost start.
OS1   same amino-acid change as a ClinVar Oncogenic variant made by another nucleotide change
      (OS1_Moderate, 2 points, when it is only Likely oncogenic).
OS3   cancerhotspots.org hotspot with >= 50 tumours at the residue and >= 10 with this change;
OM3   fewer than 50 at the residue and >= 10 with this change; OP3 fewer than 10 with this change.
      Single-residue hotspots for missense variants, in-frame indel hotspots for in-frame indels.
OM2   in-frame deletion or insertion in a known cancer gene, or stop lost in a loss-of-function gene.
OM4   another missense change at the residue is ClinVar Oncogenic (OM4_Supporting, 1 point, when
      only Likely oncogenic), and this change is at least about as different (Grantham distance at
      least 80% of the reference change's), as the standard asks.
OP1   the calibrated computational evidence reaches PP3 (the model score, or REVEL, at the PP3/BP4
      thresholds of acmg.py, BP4 withheld when SpliceAI does not rule out a splicing effect);
      SBP1 (-1) when it reaches BP4. Applied in every gene, as the standard does; in activating genes a
      damage prediction does not tell activation from inactivation, which is noted.
OP4   absent from gnomAD or below an allele frequency of 1e-5 (the standard asks for "absent or
      extremely low" without a number);
SBS1  above 1% in a gnomAD continental population (grpmax); SBVS1 above 5%.
OS2   functional studies for this protein change curated in CIViC (preclinical evidence level D, rated
      >= 3 of 5): "Oncogenic" evidence supporting oncogenicity, or "Functional" evidence of an effect
      that drives tumours in that kind of gene (gain of function or neomorphic in any driver; also loss
      of function or dominant negative in a loss-of-function gene). SBS2 (-4): evidence of unaltered
      function, or Oncogenic evidence that does not support oncogenicity. Neither when the evidence
      conflicts.

Not assessed (they need a curator): functional studies beyond CIViC, OM1 critical functional domain,
OP2 single-aetiology tumours, SBP2 synonymous variants; nor the tumour type, which the standard asks
classifications to consider. The variant's own ClinVar classification never counts as evidence for
itself: it is shown separately (``clinvar_oncogenicity``).

Clinical evidence (not points): :func:`clinical` summarises CIViC's predictive, diagnostic and
prognostic evidence for the protein change (the AMP/ASCO/CAP axis: what the variant means for
treatment), with the best evidence level (A validated or guideline ... E inferential).

Exclusions: the standard does not use OS3 with OS1 (its count mixes nucleotide changes), counts OS2 with
OS1 only for studies of the particular nucleotide change, does not use OM4 with OS1 or OS3, and does not
use OM3 with OM4; OM3 and OP3 combine with OS1. Here OS1 (with OM3 or OP3) and the evidence it excludes
(OS3, OS2, OM4, with OM3 or OP3) are alternatives and the stronger is kept, so a reference never lowers
the points; when the standard's rule would give other points, the evidence says so. OS1_Moderate, the
tool's own level for a likely oncogenic reference, is not combined with OM3 or OP3 (that reference may
itself rest on the same hotspot) and gives way to the standard's criteria on a tie. OM4 and OM3 are not
combined (the stronger one is kept; OM4 wins a tie), and OM3 comes back when OS1 replaces OM4.

Data: cancerhotspots.org (Chang et al., Cancer Discov 2018; gene symbols updated to HGNC's current
ones), the IntOGen driver compendium (Martinez-Jimenez et al., Nat Rev Cancer 2020; CC0), UniProt's
keywords (CC BY 4.0) and CIViC's monthly release (Griffith et al., Nat Genet 2017; CC0), downloaded
on first use and cached; ClinVar oncogenicity classifications through NCBI E-utilities; transcript
structures from Ensembl.
"""

from __future__ import annotations

import csv
import datetime as dt
import io
import re
import time
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any

from ._types import KeyValueCache, Row
from .aa import grantham, to_one_letter
from .api import ENSEMBL, ApiError, Http, OfflineError, SourceStatus

HOTSPOTS_URL = "https://www.cancerhotspots.org/api/hotspots/single"
INTOGEN_RELEASE = "IntOGen-Drivers-20240920"
INTOGEN_URL = f"https://www.intogen.org/download?file={INTOGEN_RELEASE}.zip"
INTOGEN_KEY = f"oncogenicity|{INTOGEN_RELEASE}|samples-and-cohorts-by-role"  # parse_intogen's counts
DATA_MAX_AGE = 90 * 86400  # both resources are released rarely
CIVIC_URL = "https://civicdb.org/downloads/{d}/{d}-ClinicalEvidenceSummaries.tsv"  # monthly release, citable
CIVIC_MAX_AGE = 35 * 86400  # a new monthly release appears on the first of each month
CIVIC_MIN_RATING = 3  # of 5 stars (CIViC: 3 = "convincing evidence")
CIVIC_FUNCTIONAL_LEVELS = {"D"}  # OS2 needs functional studies: CIViC's preclinical level (cell or animal models)
LEVELS = "ABCDE"
# UniProt's curated keywords (Swiss-Prot, human; CC BY 4.0): Tumor suppressor (KW-0043), Proto-oncogene (KW-0656)
UNIPROT_KEYWORD_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?query=%28keyword%3A{kw}%29%20AND%20%28organism_id%3A9606%29"
    "%20AND%20%28reviewed%3Atrue%29&fields=accession%2Cgene_primary&format=tsv"
)
UNIPROT_KEYWORDS = {"LoF": "KW-0043", "Act": "KW-0656"}
LOF_SHARE, ACT_SHARE = 0.6, 0.4  # IntOGen: share of mutated samples in loss-of-function cohorts
INTOGEN_MIN_COHORTS, INTOGEN_MIN_SAMPLES = 2, 20  # IntOGen decides with 2 cohorts, or with 20 samples in one
INTOGEN_FEWEST_SAMPLES = 10  # and never on fewer samples (POLE was "activating" on 4, CDKN2C on 8)
INTOGEN_OVERRIDE_COHORTS, INTOGEN_OVERRIDE_SAMPLES = 3, 50  # IntOGen against UniProt's keyword needs both
MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")  # CIViC's names
TRUNCATION_SHARE = 0.10  # ClinGen LoF guidance (Abou Tayoun 2018): removing > 10% of the protein
# cancerhotspots.org names genes by the symbols of its release (2017); the current HGNC symbols (HGNC REST,
# fetch/prev_symbol, 2026-11-25), so that hotspots match the gene names Ensembl VEP reports.
HOTSPOT_SYMBOLS = {
    "FAM58A": "CCNQ",
    "H3F3A": "H3-3A",
    "HIST1H1C": "H1-2",
    "HIST1H3B": "H3C2",
    "HIST1H3C": "H3C3",
    "HIST1H3H": "H3C10",
    "PAK7": "PAK5",
    "STK19": "WHR1",
    "TCEB1": "ELOC",
    "WHSC1": "NSD2",
}
OP4_MAX_AF = 1e-5
SBS1_AF, SBVS1_AF = 0.01, 0.05
OM4_GRANTHAM_SHARE = 0.8

POINTS = {
    "OVS1": 8,
    "OVS1_Strong": 4,
    "OVS1_Moderate": 2,
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
NOT_ASSESSED = "not assessed automatically: functional studies beyond CIViC, OM1 functional domain, OP2, tumour type"
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
    """A gene's role in cancer and where it comes from.

    IntOGen decides when its evidence is enough: at least 10 mutated samples in the cohorts that give the
    gene a role (loss-of-function or activating), from at least 2 such cohorts unless there are 20, and at
    least as many samples there as in its cohorts of ambiguous role. The role is loss-of-function when at
    least 60% of those samples are in loss-of-function cohorts, activating when at most 40%. Counting samples
    rather than cohorts keeps many small cohorts from outvoting large ones (SMARCA4 is "Act" in 27 of 45
    cohorts but loss-of-function in 49% of the mutated samples of cohorts with a role; KEAP1 7 of 14 but
    69%). UniProt's curated keyword (Tumor suppressor or Proto-oncogene, not both) decides when IntOGen is
    absent, ambiguous or too thin; its Proto-oncogene keyword also marks genes known from fusions, such as
    RUNX1, NSD1 and CBFB, whose point mutations inactivate them. When IntOGen and UniProt disagree, IntOGen
    decides only with at least 3 cohorts and 50 samples (so RUNX1 is loss-of-function); otherwise the role
    is ambiguous (SMARCB1, on 13 samples).
    """

    role: str  # Act | LoF | ambiguous
    source: str  # UniProt | IntOGen
    lof_share: float | None = None  # among the samples of cohorts with a role
    cohorts: int = 0  # cohorts that give the gene a role (loss-of-function or activating)
    samples: float = 0.0  # in cohorts with a role
    why: str = "decided"  # decided | thin | ambiguous cohorts | conflict
    other: str | None = None  # UniProt's keyword when it disagrees with IntOGen
    ambiguous_samples: float = 0.0

    @property
    def decisive(self) -> bool:
        """False when IntOGen's evidence was too thin, mostly ambiguous or contradicted by UniProt."""
        return self.why == "decided"

    def text(self) -> str:
        """'loss-of-function (UniProt: tumour suppressor)' or '... (IntOGen: 69% loss-of-function ...)'."""
        name = ROLE_NAMES.get(self.role, self.role)
        kw = {"LoF": "tumour suppressor", "Act": "proto-oncogene"}
        if self.source == "UniProt":
            return f"{name} (UniProt: {kw.get(self.role, 'tumour suppressor and proto-oncogene')})"
        k = f"{self.cohorts} cohort{'s' if self.cohorts != 1 else ''}"  # with a role
        share = f"{self.lof_share:.0%}" if self.lof_share is not None else "?"
        if self.why == "thin":
            return f"{name} (IntOGen: {self.samples:.0f} mutated samples in {k} with a role: too few to decide)"
        if self.why == "ambiguous cohorts":
            total = self.samples + self.ambiguous_samples
            return (
                f"{name} (IntOGen: {self.ambiguous_samples:.0f} of {total:.0f} mutated samples in cohorts of "
                "ambiguous role)"
            )
        if self.why == "conflict":
            return (
                f"{name} (IntOGen: {share} loss-of-function in {self.samples:.0f} samples, {k}; "
                f"UniProt: {kw.get(str(self.other), str(self.other))}: they disagree)"
            )
        return f"{name} (IntOGen: {share} loss-of-function in {self.samples:.0f} samples of cohorts with a role, {k})"


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
    roles_partial: bool = False  # UniProt's keywords could not be obtained: roles from IntOGen only
    versions: dict[str, str] = field(default_factory=dict)  # source -> release or retrieval date


def parse_hotspots(items: list[dict[str, Any]]) -> tuple[dict[str, dict[str, Hotspot]], dict[str, list]]:
    """Index cancerhotspots.org's single-residue and in-frame indel hotspots by gene."""
    residues: dict[str, dict[str, Hotspot]] = defaultdict(dict)
    indels: dict[str, list[tuple[int, int, Hotspot]]] = defaultdict(list)
    for h in items:
        gene, kind = h.get("hugoSymbol"), h.get("type")
        gene = HOTSPOT_SYMBOLS.get(str(gene), gene) if gene else gene
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


def parse_intogen(zip_bytes: bytes) -> dict[str, list[float]]:
    """Return, for each gene of IntOGen's Compendium_Cancer_Genes.tsv, the mutated samples by cohort role.

    gene -> [samples in LoF cohorts, in Act cohorts, in ambiguous cohorts, cohorts with a role (LoF or
    Act), cohorts of ambiguous role].
    """
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as z:
        name = next(n for n in z.namelist() if n.endswith("Compendium_Cancer_Genes.tsv"))
        text = z.read(name).decode("utf-8")
    out: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0, 0, 0])
    for rec in csv.DictReader(io.StringIO(text), delimiter="\t"):
        gene, role = rec.get("SYMBOL"), rec.get("ROLE")
        if not gene or role not in ("LoF", "Act", "ambiguous"):
            continue
        c = out[gene]
        c[("LoF", "Act", "ambiguous").index(role)] += _num(rec.get("SAMPLES")) or 0.0
        c[4 if role == "ambiguous" else 3] += 1
    return dict(out)


def parse_uniprot(tsv: str) -> set[str]:
    """Gene symbols (primary names) in a UniProt TSV download with a 'Gene Names (primary)' column."""
    rows = list(csv.DictReader(io.StringIO(tsv), delimiter="\t"))
    col = next((k for k in (rows[0] if rows else {}) if "primary" in k.lower()), None)
    return {g.strip() for r in rows for g in str(r.get(col) or "").split(";") if g.strip()} if col else set()


def combine_roles(intogen: dict[str, list[float]], suppressors: set[str], oncogenes: set[str]) -> dict[str, Role]:
    """Each gene's role: IntOGen's sample-weighted role when it has enough evidence, else UniProt's keyword."""
    out: dict[str, Role] = {}
    for gene in set(intogen) | suppressors | oncogenes:
        ts, og = gene in suppressors, gene in oncogenes
        uni = ("LoF" if ts else "Act") if ts != og else None
        c = intogen.get(gene)
        if not c:
            out[gene] = Role(uni or "ambiguous", "UniProt")
            continue
        lof, act, amb, cohorts = float(c[0]), float(c[1]), float(c[2]), int(c[3])  # cohorts with a role
        decided = lof + act
        share = lof / decided if decided else None

        def make(role: str, source: str, why: str = "decided", other: str | None = None) -> Role:
            return Role(role, source, share, cohorts, decided, why, other, amb)  # noqa: B023 (used at once)

        if share is None or decided < amb:  # mostly in cohorts where IntOGen could not tell the role
            out[gene] = make(uni, "UniProt") if uni else make("ambiguous", "IntOGen", "ambiguous cohorts")
        elif decided < INTOGEN_FEWEST_SAMPLES or (cohorts < INTOGEN_MIN_COHORTS and decided < INTOGEN_MIN_SAMPLES):
            out[gene] = make(uni, "UniProt") if uni else make("ambiguous", "IntOGen", "thin")
        elif ACT_SHARE < share < LOF_SHARE:
            out[gene] = make(uni, "UniProt") if uni else make("ambiguous", "IntOGen")
        else:
            role = "LoF" if share >= LOF_SHARE else "Act"
            strong = cohorts >= INTOGEN_OVERRIDE_COHORTS and decided >= INTOGEN_OVERRIDE_SAMPLES
            out[gene] = (
                make("ambiguous", "IntOGen", "conflict", uni)
                if uni and uni != role and not strong
                else make(role, "IntOGen")
            )
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

    The hotspots and roles are kept 90 days, CIViC 35 days (it is released monthly). A failed download is
    recorded in `status` (sources 'cancerhotspots.org', 'IntOGen', 'UniProt' and 'CIViC') and leaves that
    part of the data empty.
    """
    data = Data()
    today = time.strftime("%Y-%m-%d", time.gmtime())
    stored: dict[str, str] = dict(cache.get("oncogenicity|versions") or {})  # what each cached download was
    versions: dict[str, str] = {}  # the sources used in this run

    def from_cache(source: str, key: str) -> None:
        when = getattr(cache, "stored", lambda k: None)(key)
        day = time.strftime("%Y-%m-%d", time.gmtime(when)) if when else "an earlier run"  # UTC, as the run record
        versions[source] = stored.get(source) or f"cached copy retrieved {day}"

    hot = cache.get("oncogenicity|cancerhotspots-single", max_age=DATA_MAX_AGE)
    if hot is None:
        try:
            r = http.request("GET", HOTSPOTS_URL)
            if r.status_code != 200:
                raise ApiError(f"HTTP {r.status_code}", r.status_code)
            hot = r.json()
            cache.set("oncogenicity|cancerhotspots-single", hot)
            versions["cancerhotspots.org"] = stored["cancerhotspots.org"] = f"API, retrieved {today}"
        except (ApiError, OfflineError, ValueError) as e:
            if status is not None:
                status.fail("cancerhotspots.org", str(e))
            hot = None
    if hot is not None:
        if "cancerhotspots.org" not in versions:
            from_cache("cancerhotspots.org", "oncogenicity|cancerhotspots-single")
        if status is not None:
            status.use("cancerhotspots.org")
        data.residues, data.indels = parse_hotspots(hot)
        data.hotspots_ok = True
    intogen = cache.get(INTOGEN_KEY, max_age=DATA_MAX_AGE)
    if intogen is None:
        try:
            r = http.request("GET", INTOGEN_URL, follow_redirects=True)
            if r.status_code != 200:
                raise ApiError(f"HTTP {r.status_code}", r.status_code)
            intogen = parse_intogen(r.content)
            cache.set(INTOGEN_KEY, intogen)
        except (ApiError, OfflineError, ValueError, KeyError, StopIteration, zipfile.BadZipFile) as e:
            if status is not None:
                status.fail("IntOGen", str(e))
            intogen = None
    keywords = cache.get("oncogenicity|uniprot-keywords", max_age=DATA_MAX_AGE)
    if keywords is None:
        try:
            keywords = {}
            for role, kw in UNIPROT_KEYWORDS.items():
                r = http.request("GET", UNIPROT_KEYWORD_URL.format(kw=kw))
                if r.status_code != 200:
                    raise ApiError(f"HTTP {r.status_code}", r.status_code)
                keywords[role] = sorted(parse_uniprot(r.text))
                release = r.headers.get("x-uniprot-release") if hasattr(r, "headers") else None
            cache.set("oncogenicity|uniprot-keywords", keywords)
            versions["UniProt"] = stored["UniProt"] = f"release {release}" if release else f"retrieved {today}"
        except (ApiError, OfflineError, ValueError) as e:
            if status is not None:
                status.fail("UniProt", str(e))
            keywords = None
    if intogen is not None:
        versions["IntOGen"] = INTOGEN_RELEASE
        if keywords is not None and "UniProt" not in versions:
            from_cache("UniProt", "oncogenicity|uniprot-keywords")
        if status is not None:
            status.use("IntOGen")
            if keywords is not None:
                status.use("UniProt")
        found = keywords or {}
        data.roles = combine_roles(intogen, set(found.get("LoF", [])), set(found.get("Act", [])))
        data.roles_ok = True
        data.roles_partial = keywords is None
    civic = cache.get("oncogenicity|civic-evidence-monthly", max_age=CIVIC_MAX_AGE)
    if civic is None:
        try:
            text, code = None, 0
            for d in civic_releases():
                r = http.request("GET", CIVIC_URL.format(d=d), follow_redirects=True)
                code = r.status_code
                if code == 200:
                    versions["CIViC"] = stored["CIViC"] = f"monthly release {d}"
                    text = r.text
                    break
            if text is None:
                raise ApiError(f"HTTP {code}", code)
            civic = parse_civic(text)
            cache.set("oncogenicity|civic-evidence-monthly", civic)
        except (ApiError, OfflineError, ValueError, KeyError) as e:
            if status is not None:
                status.fail("CIViC", str(e))
            civic = None
    if civic is not None:
        if "CIViC" not in versions:
            from_cache("CIViC", "oncogenicity|civic-evidence-monthly")
        if status is not None:
            status.use("CIViC")
        data.civic = civic
        data.civic_ok = True
    cache.set("oncogenicity|versions", stored)
    data.versions = versions
    return data


def civic_releases(today: dt.date | None = None) -> list[str]:
    """CIViC's monthly release names to try, newest first: ['01-Nov-2026', '01-Oct-2026']."""
    d = today or dt.date.today()
    first = d.replace(day=1)
    prev = (first - dt.timedelta(days=1)).replace(day=1)
    # English month names whatever the locale (the desktop app sets the user's locale)
    return [f"01-{MONTHS[x.month - 1]}-{x.year}" for x in (first, prev)]


def transcript_structures(
    http: Http, cache: KeyValueCache, transcripts: set[str], assembly: str, status: SourceStatus | None = None
) -> dict[str, dict[str, Any]]:
    """Exons (transcript order), translated region and protein length of Ensembl transcripts (cached).

    For the strength of OVS1: the share of the protein a truncation removes and whether skipping an exon
    keeps the frame. transcript id (no version) -> {"exons": [[start, end], ...], "cds": [start, end],
    "length": amino acids}; transcripts whose lookup failed are missing.
    """
    ids = sorted({t.split(".")[0] for t in transcripts if t})
    keys = {t: f"txstruct|{assembly}|{t}" for t in ids}
    hit = cache.get_many(list(keys.values()))
    out = {t: hit[keys[t]] for t in ids if keys[t] in hit}
    todo = [t for t in ids if t not in out]
    if status is not None and todo:
        status.use("Ensembl transcripts", len(todo))
    for i in range(0, len(todo), 50):
        chunk = todo[i : i + 50]
        try:
            data = http.post_json(
                f"{ENSEMBL[assembly]}/lookup/id",
                {"ids": chunk, "expand": 1},
                headers={"Content-Type": "application/json"},
            )
        except ApiError as e:
            if status is not None:
                status.fail("Ensembl transcripts", str(e), chunk)
            continue
        res = {}
        for t in chunk:
            d = (data or {}).get(t) or {}
            tr, exons = d.get("Translation") or {}, d.get("Exon") or []
            if not tr or not exons:
                continue
            ex = sorted([int(e["start"]), int(e["end"])] for e in exons)
            if d.get("strand") == -1:
                ex.reverse()
            res[t] = {"exons": ex, "cds": [int(tr["start"]), int(tr["end"])], "length": int(tr.get("length") or 0)}
        out.update(res)
        cache.set_many({keys[t]: v for t, v in res.items()})
    return out


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
    same_change: tuple[str, str] | None = None  # OS1 candidate, kept only if stronger than the alternatives
    hotspot: tuple[str, str] | None = None  # the residue's hotspot criterion, even when OM4 replaced OM3
    sop_delta: int = 0  # points the standard's rule for OS1 would add (negative: remove) to this assessment

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
    structure: dict[str, Any] | None = None,
) -> Assessment:
    """Return the oncogenicity criteria a row meets.

    Args:
        row: One results row; reads ``gene``, ``transcript``, ``canonical``, ``protein_position``,
            ``ref_aa``, ``alt_aa``, ``protein_change``, ``nmd``, ``gnomad_af`` and ``gnomad_popmax_af``.
        terms: The row's VEP consequence terms.
        model_call: ``"damaging"`` when the calibrated computational evidence reaches PP3 (supporting or
            stronger), ``"benign"`` when it reaches BP4, otherwise None.
        data: Hotspots and gene roles, from :func:`load`.
        refs: ClinVar oncogenic variants at the residue (:func:`references`); None when not searched.
        gnomad_info: gnomAD's own answer for a variant Ensembl gave no frequency (``found``, ``af``).
        structure: The transcript's exons and translation (:func:`transcript_structures`), for the
            strength of OVS1; None when unknown.
    """
    a = Assessment()
    gene = str(row.get("gene") or "")
    role = data.roles.get(gene)
    if data.roles_ok:
        a.notes.append(f"{gene}: {role.text()}" if role else f"{gene} is not a known cancer gene (UniProt, IntOGen)")
        if data.roles_partial:
            a.notes.append("UniProt keywords unavailable: gene roles from IntOGen only")
    missense = "missense_variant" in terms
    if terms & NULL_TERMS:
        _null(a, row, role, terms, structure)
    elif missense:
        _missense(a, row, model_call, data, refs, role)
    elif "stop_lost" in terms or (terms & INFRAME_TERMS and re.search(r"del|ins|dup", str(row.get("protein_change")))):
        _inframe(a, row, terms, data, role)
    _functional(a, row, data, role)
    _same_change(a)
    _frequency(a, row, gnomad_info)
    if a.sop_delta:
        a.notes.append(f"the standard's rule for OS1 would give {a.points + a.sop_delta} points")
    a.notes.append(NOT_ASSESSED)
    return a


def _null(a: Assessment, row: Row, role: Role | None, terms: set[str], structure: dict[str, Any] | None) -> None:
    """OVS1 for a null variant in a loss-of-function gene, graded as ClinGen's LoF guidance grades PVS1.

    The standard asks for caution with truncations at the extreme 3' end, splice variants that skip an
    exon in frame and alternative isoforms, recommends ClinGen's guidance on loss-of-function variants
    (Abou Tayoun et al. 2018), and lets a criterion be downgraded when the evidence falls short. So:
    a truncation expected to undergo nonsense-mediated decay gets OVS1; one that escapes it (last exon,
    single-exon gene) OVS1_Strong when it removes more than 10% of the protein and OVS1_Moderate
    otherwise; a canonical splice variant OVS1 when skipping its exon shifts the frame, and OVS1_Strong or
    OVS1_Moderate (more or less than 10% of the coding sequence) when the skip is in frame or the
    structure is unknown; a lost start codon OVS1_Moderate (another start codon may be used).
    """
    if not role or role.role != "LoF":
        why = role.text() if role else "not a known cancer gene"
        a.notes.append(f"OVS1 not applied (not a loss-of-function gene: {why})")
        return
    a.codes.append(_ovs1(row, terms, structure))


def _ovs1(row: Row, terms: set[str], structure: dict[str, Any] | None) -> tuple[str, str]:
    """The OVS1 code (full, Strong or Moderate) of a null variant in a loss-of-function gene, with its reason."""
    if "start_lost" in terms and not terms & {"stop_gained", "frameshift_variant"}:
        return "OVS1_Moderate", "start codon lost: another in-frame start codon may be used (ClinGen LoF guidance)"
    if terms & {"splice_donor_variant", "splice_acceptor_variant"} and not terms & {
        "stop_gained",
        "frameshift_variant",
    }:
        skip = skipped_exon(row, terms, structure)
        if skip is None:
            return "OVS1_Strong", "canonical splice site; whether skipping the exon keeps the frame is unknown"
        exon, coding, total = skip
        if coding % 3:
            return "OVS1", f"canonical splice site; skipping exon {exon} would shift the reading frame"
        share = coding / total if total else 0.0
        code = "OVS1_Strong" if share > TRUNCATION_SHARE else "OVS1_Moderate"
        return (
            code,
            f"canonical splice site; skipping exon {exon} would remove {share:.0%} of the coding sequence in frame",
        )
    nmd = str(row.get("nmd") or "")
    # NMD expected, or uncertain at the penultimate exon; a change across an exon-intron boundary is graded
    # like an escape, by the share of the protein removed (its splice effect is unknown)
    if "escape" not in nmd and "boundary" not in nmd:
        return "OVS1", "truncating variant expected to undergo nonsense-mediated decay" + (
            f" ({nmd})" if "uncertain" in nmd else ""
        )
    length = _num(row.get("protein_length")) or _num((structure or {}).get("length"))
    pos = _num(row.get("protein_position"))
    if not length or not pos:
        return "OVS1_Moderate", f"{nmd}; the share of the protein removed is unknown"
    share = max(0.0, (length - pos + 1) / length)
    code = "OVS1_Strong" if share > TRUNCATION_SHARE else "OVS1_Moderate"
    return code, f"{nmd}; the truncation removes {share:.0%} of the protein (ClinGen LoF guidance: > 10% strong)"


def skipped_exon(row: Row, terms: set[str], structure: dict[str, Any] | None) -> tuple[int, int, int] | None:
    """The exon a canonical splice variant would skip: (exon number, its coding bases, all coding bases).

    VEP numbers introns in transcript order: a donor variant in intron i acts on exon i, an acceptor
    variant on exon i + 1. Coding bases are the exon's overlap with the translated region.
    """
    if not structure or not structure.get("exons") or not structure.get("cds"):
        return None
    try:
        intron = int(str(row.get("intron") or "").split("/")[0])
    except ValueError:
        return None
    exon = intron + (1 if "splice_acceptor_variant" in terms and "splice_donor_variant" not in terms else 0)
    exons = structure["exons"]
    if not 1 <= exon <= len(exons):
        return None
    lo, hi = structure["cds"]

    def coding(e: list[int]) -> int:
        return max(0, min(e[1], hi) - max(e[0], lo) + 1)

    total = sum(coding(e) for e in exons)
    return exon, coding(exons[exon - 1]), total


def _missense(
    a: Assessment, row: Row, model_call: str | None, data: Data, refs: list[Reference] | None, role: Role | None
) -> None:
    """OS1/OM4 (ClinVar), OS3/OM3/OP3 (hotspots) and OP1/SBP1 (the model) for a missense variant."""
    ref_aa, alt_aa, pos = row.get("ref_aa"), row.get("alt_aa"), row.get("protein_position")
    same = [x for x in refs or [] if x.same_change]
    other = [x for x in refs or [] if not x.same_change]
    if refs is None:
        a.notes.append("OS1/OM4 not checked (ClinVar oncogenicity search unavailable)")
    if same:  # resolved in _same_change, once the hotspot and CIViC evidence are known
        x = same[0]
        a.same_change = (
            "OS1" if x.level == "O" else "OS1_Moderate",
            f"same change by {x.cdna}: {_level(x)}, {x.accession}",
        )
    hot = _residue_hotspot(row, data)
    if hot:
        n_aa = hot.changes.get(str(alt_aa), 0)
        code = "OS3" if hot.tumours >= 50 and n_aa >= 10 else "OM3" if n_aa >= 10 else "OP3"
        a.hotspot = (code, f"hotspot {ref_aa}{pos}: {hot.tumours} tumours, {n_aa} with {alt_aa}; cancerhotspots.org")
        a.codes.append(a.hotspot)
    has = {c for c, _ in a.codes}
    if other and "OS3" not in has:  # with a same-change reference, OM4 is weighed against OS1 in _same_change
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
    if model_call == "damaging":
        a.codes.append(("OP1", "calibrated computational evidence reaches PP3"))
    elif model_call == "benign":
        a.codes.append(("SBP1", "calibrated computational evidence reaches BP4"))
    if model_call and role and role.role == "Act":
        a.notes.append("in an activating gene, a damage prediction does not tell activation from inactivation")


SAME_CHANGE = ("OS3", "OS2", "OM4", "OM4_Supporting")  # the evidence that the standard does not use with OS1
WITH_OS1 = ("OM3", "OP3")  # hotspot criteria that the standard combines with OS1


def _same_change(a: Assessment) -> None:
    """OS1 from a ClinVar classification of the same amino-acid change, or the evidence it excludes.

    The standard does not use OS3 with OS1 (unless the hotspot count is specific to the nucleotide change),
    counts OS2 with OS1 only for studies of the particular nucleotide change, and does not use OM4 with OS1;
    OM3 and OP3 combine with OS1, and OM3, which cannot be used with OM4, comes back when OS1 replaces OM4
    (Horak et al. 2022, Table 1). Here the two ways of counting the change, OS1 with OM3 or OP3 and the
    evidence it excludes with OM3 or OP3, are alternatives and the stronger is kept, so that adding a
    reference never lowers the points: this departs from the standard when the excluded evidence is
    stronger. A likely oncogenic reference (OS1_Moderate, two points, the tool's own level) is weighed
    alone against all the evidence for the change, OM3 and OP3 included, because its classification may
    itself rest on the same hotspot; it is kept only when it is stronger, never on a tie.
    """
    if a.same_change is None:
        return
    code, why = a.same_change
    alt = [(c, w) for c, w in a.codes if c in SAME_CHANGE]
    alt_points = sum(POINTS[c] for c, _ in alt) + sum(POINTS[c] for c, _ in a.codes if c in WITH_OS1)
    # OM3/OP3 kept (or restored) with OS1, as the standard does; not with the tool's own OS1_Moderate
    minor = a.hotspot if code == "OS1" and a.hotspot and a.hotspot[0] in WITH_OS1 else None
    os1_points = POINTS[code] + (POINTS[minor[0]] if minor else 0)
    if os1_points > alt_points or (os1_points == alt_points and code == "OS1"):
        # OM3/OP3 are left out only with OS1_Moderate; with OS1 they stay (or come back as `minor`)
        dropped = [c for c, _ in a.codes if c in SAME_CHANGE or (c in WITH_OS1 and minor is None)]
        a.codes = [(c, w) for c, w in a.codes if c not in SAME_CHANGE and c not in WITH_OS1]
        a.codes.append((code, why))
        if minor:
            a.codes.append(minor)
        if dropped:
            a.notes.append(f"{', '.join(dropped)} not counted ({code} uses this amino-acid change)")
        chosen = os1_points
    else:
        parts = [c for c, _ in alt] + [c for c, _ in a.codes if c in WITH_OS1]
        than = "is stronger" if alt_points > os1_points else "counts as much and is the standard's"
        a.notes.append(f"{code} not counted ({why}; {' + '.join(parts)} for this change {than})")
        chosen = alt_points
    # the standard's rule: OS1 (with OM3 or OP3) from an oncogenic reference replaces the evidence it
    # excludes; a likely oncogenic reference gives no OS1, so the other evidence stands
    standard = POINTS["OS1"] + (POINTS[minor[0]] if minor else 0) if code == "OS1" else alt_points
    a.sop_delta = standard - chosen


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
        a.notes.append("OM2 not applied (not a known cancer gene: UniProt, IntOGen)")
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
    clinical_only = [x for x in items if x["level"] not in CIVIC_FUNCTIONAL_LEVELS]
    items = [x for x in items if x["level"] in CIVIC_FUNCTIONAL_LEVELS]
    if clinical_only and not items:
        a.notes.append("CIViC oncogenic evidence from clinical observations (not functional studies) not used")
        return
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
    elif af is None:
        from .gnomad import covered

        if gnomad_info and not gnomad_info.get("found") and covered(gnomad_info):
            a.codes.append(("OP4", f"absent from {gnomad_info.get('dataset') or 'gnomAD'}"))
        elif gnomad_info and not gnomad_info.get("found"):
            a.notes.append("OP4 not applied (absent from gnomAD, but few reads there)")
        else:
            a.notes.append("OP4 not checked (absence from gnomAD not confirmed)")
    elif af == 0 or af < OP4_MAX_AF:
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
