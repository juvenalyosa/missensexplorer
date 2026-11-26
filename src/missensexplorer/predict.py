"""End-to-end pipeline: inputs -> VEP -> predictors -> protein context -> one table."""

from __future__ import annotations

import gzip
import hashlib
import io
import logging
import os
import platform
import re
import time
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any
from urllib.parse import unquote

import pandas as pd

from . import __version__, acmg, gnomad, model, residue, vcep
from . import clinvar as cv
from . import oncogenicity as onco
from . import protein as prot
from ._types import JSON, KeyValueCache, Progress, Row
from .aa import AA3_TO_1, substitution_profile
from .api import ENSEMBL, ApiError, Http, OfflineError, SourceStatus, ensembl_release
from .cache import Cache, NullCache
from .inputs import CHROMS, ParsedInput, VariantInput, _assembly_from_header, parse_text, read_file, region_query
from .vep import BASE_PARAMS, DETAIL_CONSEQUENCES, PREDICTOR_PARAMS, Vep, on_canonical
from .vep_local import LocalVep

log = logging.getLogger(__name__)

LOF_TERMS = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost"}
INFRAME_TERMS = {"inframe_insertion", "inframe_deletion", "protein_altering_variant", "stop_lost"}
SPLICE_NEAR_TERMS = {
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
}
PROTEIN_ALTERING = LOF_TERMS | INFRAME_TERMS | {"missense_variant"}

CLINVAR_ORDER = [
    "pathogenic",
    "pathogenic/likely_pathogenic",
    "likely_pathogenic",
    "risk_factor",
    "drug_response",
    "association",
    "affects",
    "protective",
    "uncertain_significance",
    "conflicting_interpretations_of_pathogenicity",
    "conflicting_classifications_of_pathogenicity",
    "likely_benign",
    "benign/likely_benign",
    "benign",
    "other",
    "not_provided",
]

# Calibrated thresholds (Pejaver et al., AJHG 2022) for PP3/BP4 "supporting" evidence.
REVEL_DAMAGING, REVEL_BENIGN = 0.644, 0.290
CADD_DAMAGING, CADD_BENIGN = 25.3, 22.7

COLUMNS = [
    "input",
    "source",
    "sample",
    "zygosity",
    "genotype",
    "carriers",
    "gene",
    "protein_change",
    "hgvsp",
    "hgvsc",
    "consequence",
    "impact",
    "prediction",
    "effect_score",
    "score_with_frequency",
    "prediction_basis",
    "acmg_codes",
    "acmg_pp3_bp4",
    "acmg_points",
    "acmg_ps1_pm5",
    "acmg_pm2",
    "population_frequency",
    "vcep_specification",
    "acmg_notes",
    "oncogenicity",
    "oncogenicity_points",
    "oncogenicity_codes",
    "oncogenicity_evidence",
    "clinical_evidence_level",
    "clinical_evidence",
    "alphamissense",
    "alphamissense_class",
    "revel",
    "cadd_phred",
    "eve",
    "eve_class",
    "popeve",
    "sift",
    "sift_score",
    "polyphen",
    "polyphen_score",
    "spliceai_max",
    "spliceai_effect",
    "grantham",
    "grantham_class",
    "blosum62",
    "hydropathy_change",
    "volume_change",
    "charge_change",
    "aa_class_change",
    "physchem_notes",
    "protein_name",
    "uniprot",
    "protein_position",
    "protein_length",
    "ref_aa",
    "alt_aa",
    "protein_window",
    "mutant_protein",
    "protein_features",
    "key_site",
    "uniprot_variants_here",
    "uniprot_mutagenesis_here",
    "plddt",
    "plddt_class",
    "alphafold_model",
    "protein_function",
    "clinvar",
    "clinvar_review",
    "clinvar_stars",
    "clinvar_conditions",
    "clinvar_oncogenicity",
    "clinvar_ids",
    "gnomad_af",
    "gnomad_popmax_af",
    "rsid",
    "cosmic",
    "pubmed_count",
    "pubmed_ids",
    "chrom",
    "pos",
    "ref",
    "alt",
    "assembly",
    "variant_class",
    "transcript",
    "mane_select",
    "canonical",
    "protein_id",
    "exon",
    "intron",
    "codons",
    "cds_position",
    "nmd",
    "vcf_id",
    "vcf_filter",
    "notes",
]


@dataclass
class PredictOptions:
    """Options of one :func:`predict` run.

    Attributes:
        assembly: auto | GRCh38 | GRCh37
        default_gene: gene used for bare changes like "V600E"
        genes: keep only variants in these genes
        sample: VCF sample to report (default: all carriers)
        pass_only: VCF: keep FILTER=PASS only
        include_noncoding: also report intronic/UTR/intergenic consequences
        predictors: AlphaMissense, REVEL, CADD, EVE, SpliceAI
        protein_context: UniProt features + AlphaFold pLDDT
        clinvar: exact-allele ClinVar classification (NCBI)
        sequences: mutant protein sequences
        vep: rest (public Ensembl API) | local (VEP on this computer)
        vep_command: local VEP command (default: $MISSENSEXPLORER_VEP or 'vep --cache --offline')
        vep_args: extra local VEP arguments, e.g. plugin data files
        offline: never use the internet (local VEP only; no ClinVar/UniProt/AlphaFold)
        common_af: gnomAD grpmax at/above which damaging/uncertain calls become benign (None = off)
        acmg_source: PP3/BP4 from: model (own ClinGen-method calibration) | revel (ClinGen's)
        gnomad_check: ask gnomAD directly about variants Ensembl gives no frequency (PM2)
        oncogenicity: somatic oncogenicity points (ClinGen/CGC/VICC) from cancerhotspots.org, the IntOGen
            driver genes, CIViC and ClinVar oncogenicity classifications, and CIViC's clinical evidence
            (see oncogenicity.py)
        workers: parallel requests to Ensembl VEP, UniProt and AlphaFold
        cache: keep the web answers in an on-disk cache (False: ask the services every time)
        cache_path: cache file (None: the default cache directory)
    """

    assembly: str = "auto"
    default_gene: str | None = None
    genes: list[str] | None = None
    sample: str | None = None
    pass_only: bool = False
    include_noncoding: bool = False
    predictors: bool = True
    protein_context: bool = True
    clinvar: bool = True
    sequences: bool = True
    vep: str = "rest"
    vep_command: str | None = None
    vep_args: str | None = None
    offline: bool = False
    common_af: float | None = 0.05
    acmg_source: str = "model"
    gnomad_check: bool = True
    oncogenicity: bool = False
    workers: int = 4
    cache: bool = True
    cache_path: str | None = None


@dataclass
class PredictionResult:
    """Everything one :func:`predict` run produced.

    Attributes:
        table: the results, one row per variant and transcript, with the columns of :data:`COLUMNS`
        errors: inputs that could not be annotated (columns input, source, error)
        assembly: genome build used (GRCh38 or GRCh37)
        warnings: messages for the user (genome build, failed services, multi-sample VCF, ...)
        stats: counts for the summary (inputs, rows, damaging calls, errors, seconds, VCF counts)
        sequences: row index -> {wt, mut, header, note}
        alphafold: uniprot -> AlphaFold info
        uniprot: UniProt accession -> UniProt entry used for the protein context
        run: provenance: versions, options, inputs, sources
        sources: source -> ok | failed: ... | skipped: ...
    """

    table: pd.DataFrame
    errors: pd.DataFrame
    assembly: str
    warnings: list[str] = field(default_factory=list)
    stats: dict[str, Any] = field(default_factory=dict)
    sequences: dict[int, dict[str, str]] = field(default_factory=dict)
    alphafold: dict[str, dict[str, Any]] = field(default_factory=dict)
    uniprot: dict[str, dict[str, Any]] = field(default_factory=dict)
    run: dict[str, Any] = field(default_factory=dict)
    sources: dict[str, str] = field(default_factory=dict)


# --------------------------------------------------------------------------- small helpers

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    """Return the reverse complement of a DNA sequence; '-' and '' (no bases) are returned unchanged."""
    return s if s in ("-", "") else s.translate(_COMP)[::-1]


def _f(x: Any) -> float | None:
    """Return ``x`` as a float, or None when it is None, '' or not a number."""
    try:
        return None if x is None or x == "" else float(x)
    except (TypeError, ValueError):
        return None


def _max_listed(s: Any) -> float | None:
    """Return the largest number of a per-transcript score list, or None.

    DbNSFP returns per-transcript lists like '.,.,0.931,.'.
    """
    if s is None:
        return None
    if isinstance(s, (int, float)):
        return float(s)
    vals = [v for v in (_f(x) for x in str(s).split(",")) if v is not None]
    return max(vals) if vals else None


def short_protein_change(hgvsp: str | None) -> str | None:
    """Return the one-letter form of an HGVS protein change.

    'p.Val600Glu' -> 'V600E', 'p.Arg213Ter' -> 'R213*', 'p.Thr125ProfsTer23' -> 'T125Pfs*23'.

    Args:
        hgvsp: HGVS protein notation, with or without the protein ID prefix ('ENSP...:p.Val600Glu');
            URL-encoded text (as VEP returns it) is decoded.

    Returns:
        The one-letter notation, or None when ``hgvsp`` is empty.
    """
    if not hgvsp:
        return None
    s = unquote(hgvsp).split(":")[-1]
    s = s.removeprefix("p.").strip("()")

    def rep(m: re.Match[str]) -> str:
        tok = m.group(0)
        return AA3_TO_1.get(tok.upper(), tok)

    return re.sub(r"[A-Z][a-z]{2}", rep, s)


def _strip_ref(h: str | None) -> str | None:
    """Return an HGVS name without its sequence prefix ('ENST...:c.1799T>A' -> 'c.1799T>A'), URL-decoded."""
    return unquote(h.split(":", 1)[1]) if h and ":" in h else h


def _clinvar_label(sigs: list[str]) -> str | None:
    """Join the distinct ClinVar significances, most severe first; 'not_provided' only when alone."""
    uniq = list(dict.fromkeys(s.strip() for s in sigs if s and s.strip()))
    if len(uniq) > 1:
        uniq = [s for s in uniq if s != "not_provided"]
    uniq.sort(key=lambda s: CLINVAR_ORDER.index(s) if s in CLINVAR_ORDER else len(CLINVAR_ORDER))
    return ", ".join(uniq) or None


def colocated_info(colocated: list[dict[str, Any]], alt_plus: str | None, alt_raw: str | None = None) -> Row:
    """Collect ClinVar/gnomAD/dbSNP/COSMIC/PubMed data for one allele.

    clin_sig_allele is keyed by the forward-strand allele, while VEP keys
    `frequencies` by the allele as given in the input (reverse strand for HGVS
    on minus-strand genes), so both spellings are needed.

    Args:
        colocated: VEP's ``colocated_variants`` for one result.
        alt_plus: the ALT allele on the forward strand; None takes the ClinVar significances of every allele
            (and no frequency).
        alt_raw: the ALT allele as given in the input (VEP's ``variant_allele``).

    Returns:
        Row fields ``rsid``, ``gnomad_af``, ``gnomad_popmax_af`` (grpmax), ``cosmic``, ``pubmed_count`` and
        ``pubmed_ids`` (at most 25), plus the private ``_clinvar_vep`` (ClinVar significance from VEP's data)
        and ``_vcvs`` (ClinVar VCV accessions), which the pipeline removes later.
    """
    out: Row = {
        "rsid": None,
        "_clinvar_vep": None,
        "_vcvs": [],
        "gnomad_af": None,
        "gnomad_popmax_af": None,
        "cosmic": None,
        "pubmed_count": 0,
        "pubmed_ids": None,
    }
    sigs: list[str] = []
    vcv: list[str] = []
    cosmic: list[str] = []
    pubmed: list[JSON] = []
    rs: list[str] = []
    af: float | None = None
    popmax: float | None = None
    for c in colocated or []:
        cid = c.get("id") or ""
        if cid.startswith("rs"):
            rs.append(cid)
        if cid.startswith("COSV"):
            cosmic.append(cid)
        for part in (c.get("clin_sig_allele") or "").split(";"):
            if ":" in part:
                a, sig = part.split(":", 1)
                if alt_plus is None or a == alt_plus:
                    sigs.append(sig)
        syn = c.get("var_synonyms") or {}
        vcv += [x for x in syn.get("ClinVar", []) if str(x).startswith("VCV")]
        pubmed += c.get("pubmed") or []
        freqs = c.get("frequencies") or {}
        fr = freqs.get(alt_raw) if alt_raw in freqs else freqs.get(alt_plus) if alt_plus else None
        if fr:
            overall = [fr[k] for k in ("gnomade", "gnomadg", "af") if fr.get(k) is not None]
            if overall:
                af = max(overall + ([af] if af is not None else []))
            # grpmax: main genetic ancestry groups only (bottlenecked groups excluded, as gnomAD does)
            pops = [
                v
                for k, v in fr.items()
                if k.startswith(("gnomade_", "gnomadg_")) and v is not None and k.split("_", 1)[1] in acmg.GRPMAX_GROUPS
            ]
            if pops:
                popmax = max(pops + ([popmax] if popmax is not None else []))
    out["rsid"] = ", ".join(dict.fromkeys(rs)) or None
    out["_clinvar_vep"] = _clinvar_label(sigs)
    out["_vcvs"] = list(dict.fromkeys(vcv))
    out["cosmic"] = ", ".join(dict.fromkeys(cosmic)) or None
    pm = list(dict.fromkeys(str(p) for p in pubmed))
    out["pubmed_count"] = len(pm)
    out["pubmed_ids"] = ", ".join(pm[:25]) or None
    out["gnomad_af"] = af
    out["gnomad_popmax_af"] = popmax
    return out


# --------------------------------------------------------------------------- verdict

# Vote weights: modern, well-calibrated predictors count fully; SIFT/PolyPhen
# (older, conservation-only) count half.
WEIGHTS = {"AlphaMissense": 1.0, "REVEL": 1.0, "EVE": 1.0, "CADD": 0.75, "PolyPhen": 0.5, "SIFT": 0.5}


def _num(x: Any) -> float | None:
    """Return ``x`` as a float, or None for missing values (None, NaN, '', non-numeric)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _str(x: Any) -> str | None:
    """Return ``x`` if it is a non-empty string, else None (so NaN in table rows counts as missing)."""
    return x if isinstance(x, str) and x else None


MODERN = ("AlphaMissense", "REVEL", "EVE", "CADD")


def _modern_present(row: Row) -> bool:
    """Return True when any of the MODERN predictors has a value for the row."""
    return bool(
        _str(row.get("alphamissense_class"))
        or _num(row.get("revel")) is not None
        or _num(row.get("cadd_phred")) is not None
        or _str(row.get("eve_class"))
    )


def _votes(row: Row) -> tuple[list[str], list[str], list[str]]:
    """Return the damaging/benign votes per predictor.

    Accepts pipeline dicts and table rows (NaN = missing).

    Returns:
        The predictors that vote damaging, those that vote benign, and one description per available
        score (e.g. 'REVEL 0.93').
    """
    dmg: list[str] = []
    ben: list[str] = []
    detail: list[str] = []
    am, amc = _num(row.get("alphamissense")), _str(row.get("alphamissense_class"))
    if amc:
        detail.append(f"AlphaMissense {am:.2f} ({amc.replace('_', ' ')})" if am is not None else f"AlphaMissense {amc}")
        (dmg if amc == "likely_pathogenic" else ben if amc == "likely_benign" else []).append("AlphaMissense")
    rv = _num(row.get("revel"))
    if rv is not None:
        detail.append(f"REVEL {rv:.2f}")
        (dmg if rv >= REVEL_DAMAGING else ben if rv <= REVEL_BENIGN else []).append("REVEL")
    cd = _num(row.get("cadd_phred"))
    if cd is not None:
        detail.append(f"CADD {cd:.1f}")
        (dmg if cd >= CADD_DAMAGING else ben if cd <= CADD_BENIGN else []).append("CADD")
    ev = _str(row.get("eve_class"))
    if ev:
        detail.append(f"EVE {ev}")
        (dmg if ev.lower() == "pathogenic" else ben if ev.lower() == "benign" else []).append("EVE")
    sf = _str(row.get("sift"))
    if sf:
        detail.append(f"SIFT {sf.replace('_', ' ')}")
        (dmg if sf.startswith("deleterious") else ben if sf.startswith("tolerated") else []).append("SIFT")
    pp = _str(row.get("polyphen"))
    if pp:
        detail.append(f"PolyPhen {pp.replace('_', ' ')}")
        (dmg if pp == "probably_damaging" else ben if pp == "benign" else []).append("PolyPhen")
    return dmg, ben, detail


def classify(row: Row, terms: set[str]) -> tuple[str, float | None, str]:
    """Return (prediction label, effect score 0-1, human-readable basis).

    Args:
        row: a pipeline row or a results-table row (NaN counts as missing).
        terms: the row's VEP consequence terms.

    Returns:
        The label (see :func:`prediction_status` for how labels are grouped), the effect score (None when
        no score applies) and the evidence behind the call.
    """
    sp = _num(row.get("spliceai_max"))
    if terms & LOF_TERMS:
        if terms & {"stop_gained", "frameshift_variant"}:
            kind = "truncating"
        elif terms & {"splice_acceptor_variant", "splice_donor_variant"}:
            kind = "canonical splice site"
        else:
            kind = "start codon lost"
        basis = [x for x in (row.get("nmd"), row.get("_truncation")) if x]
        return f"Loss of function ({kind})", 1.0, "; ".join(basis) or "HIGH impact consequence"
    if "missense_variant" in terms:
        return _classify_missense(row, sp)
    if terms & INFRAME_TERMS:
        cd = _num(row.get("cadd_phred"))
        if cd is not None and cd >= CADD_DAMAGING:
            return "Protein-altering in-frame (likely damaging)", 0.75, f"CADD {cd:.1f}"
        return (
            "Protein-altering in-frame",
            0.5 if cd is None else 0.25,
            f"CADD {cd:.1f}" if cd is not None else "no scores available",
        )
    if sp is not None and sp >= 0.5:
        return "Likely splice-altering", round(sp, 2), f"SpliceAI {sp:.2f}"
    if sp is not None and sp >= 0.2:
        return "Possible splice effect", round(sp, 2), f"SpliceAI {sp:.2f}"
    if terms & {"synonymous_variant", "stop_retained_variant", "start_retained_variant"}:
        return "Silent (no amino-acid change)", 0.0, f"SpliceAI {sp:.2f}" if sp is not None else ""
    if terms & SPLICE_NEAR_TERMS:
        return (
            "Splice region (no strong splice signal)",
            0.1 if sp is not None else None,
            f"SpliceAI {sp:.2f}" if sp is not None else "",
        )
    return "Non-coding / no protein change", None, ""


def _classify_missense(row: Row, sp: float | None) -> tuple[str, float | None, str]:
    """Classify a missense variant with the model, or with a weighted predictor vote when it has no score."""
    dmg, ben, detail = _votes(row)
    basis = "; ".join(detail)
    if sp is not None and sp >= 0.2:
        basis += f"; SpliceAI {sp:.2f} (may also affect splicing)"
    # Main path: logistic regression over the scores (model.py), trained on ClinVar; it beat
    # the vote on held-out genes (benchmarks/results/clinvar_benchmark.md).
    ms = model.score(row)
    if ms is not None:
        dm, bm = [x for x in dmg if x in MODERN], [x for x in ben if x in MODERN]
        if ms.p >= model.DAMAGING:
            label = "Likely damaging"
        elif ms.p <= model.BENIGN:
            label = "Likely benign"
        elif dm and bm:
            label = "Uncertain (conflicting predictors)"
        else:
            label = "Uncertain (intermediate scores)"
        summary = f"{model.describe(ms)}; {len(dm)} damaging / {len(bm)} benign predictors"
        return label, round(ms.p, 2), f"{summary} | {basis}"
    return _classify_by_vote(row, dmg, ben, basis)


def _classify_by_vote(row: Row, dmg: list[str], ben: list[str], basis: str) -> tuple[str, float | None, str]:
    """Classify a missense variant by the weighted vote of the predictors in ``dmg`` and ``ben``."""
    # Fallback without AlphaMissense and REVEL: weighted vote. SIFT and PolyPhen-2 only vote when
    # no modern predictor has a score: on expert-reviewed ClinVar variants they added wrong calls.
    modern = _modern_present(row)
    limited = not modern and bool(dmg or ben)
    if modern:
        dmg = [x for x in dmg if x in MODERN]
        ben = [x for x in ben if x in MODERN]
    n = len(dmg) + len(ben)
    if n == 0:
        if modern:
            return "Uncertain (intermediate scores)", 0.5, basis
        return "Missense (no predictor scores)", None, basis or "no scores available"
    wd, wb = sum(WEIGHTS[x] for x in dmg), sum(WEIGHTS[x] for x in ben)
    frac = wd / (wd + wb)
    summary = f"{len(dmg)} damaging / {len(ben)} benign predictors (weighted {frac:.0%} damaging)"
    suffix = " (limited evidence: SIFT/PolyPhen only)" if limited else ""
    if n == 1 and not limited:
        label = "Likely damaging (single predictor)" if dmg else "Likely benign (single predictor)"
    elif frac >= 0.7:
        label = "Likely damaging" + suffix
    elif frac <= 0.3:
        label = "Likely benign" + suffix
    else:
        label = "Uncertain (conflicting predictors)"
    return label, round(frac, 2), f"{summary} | {basis}"


def is_damaging(label: str | None) -> bool:
    """Return True for the labels that count as damaging (likely damaging, loss of function, splice-altering)."""
    label = label or ""
    return (
        label.startswith(("Likely damaging", "Loss of function", "Likely splice-altering"))
        or label == "Protein-altering in-frame (likely damaging)"
    )


STATUSES = ["Damaging", "Uncertain", "Benign / silent", "Not assessed"]


def prediction_status(label: Any) -> str:
    """Collapse a prediction label into one of STATUSES (for colour, filters and charts).

    Args:
        label: a prediction label; anything that is not a string (None, NaN) counts as "Not assessed".
    """
    label = label if isinstance(label, str) else ""
    if is_damaging(label):
        return "Damaging"
    if "benign" in label.lower() or label.startswith("Silent"):
        return "Benign / silent"
    if label.startswith(("Uncertain", "Possible", "Protein-altering", "Splice region")):
        return "Uncertain"
    return "Not assessed"


def apply_frequency(label: str, basis: str, grpmax: float | None, common_af: float | None) -> tuple[str, str]:
    """BA1-style rule: a damaging or uncertain call becomes benign when the variant is common in gnomAD.

    Kept separate from the score so it stays visible: the basis still shows the protein-level call.

    Args:
        label: the protein-level prediction label.
        basis: the evidence behind it.
        grpmax: gnomAD grpmax allele frequency (0-1); None or NaN when unknown.
        common_af: frequency at/above which the rule applies; None or 0 turns it off.

    Returns:
        The (label, basis) pair, changed only when the rule applies.
    """
    af = _num(grpmax)
    if not common_af or af is None or af < common_af or prediction_status(label) not in ("Damaging", "Uncertain"):
        return label, basis
    note = (
        f"common in gnomAD: grpmax {af * 100:.2g}% >= {common_af * 100:g}% (protein-level call: {label}); "
        "frequent variants can still cause recessive or founder disease"
    )
    return "Likely benign (common in population)", f"{note} | {basis}"


def _nmd(exon: str | None, terms: set[str]) -> str | None:
    """Return the expected nonsense-mediated decay outcome of a stop-gain or frameshift, from VEP's exon ('3/10').

    A change that spans exons ('9-10/10') also crosses an exon-intron boundary and removes splice sites, so the
    fate of the transcript cannot be told from the position of the new stop codon: NMD is uncertain.
    """
    if not exon or not terms & {"stop_gained", "frameshift_variant"}:
        return None
    try:
        first, total = exon.split("/")[0:2]
        a, b = int(first.split("-")[-1]), int(total)
    except ValueError:
        return None
    if "-" in first:
        return "NMD uncertain (the change spans an exon-intron boundary)"
    if b == 1:
        return "NMD escape likely (single-exon gene)"
    if a == b:
        return "NMD escape likely (last exon)"
    if a == b - 1:
        return "NMD uncertain (penultimate exon; depends on distance to last junction)"
    return "NMD likely (premature stop upstream of last exon)"


# --------------------------------------------------------------------------- rows


def _rows_for_variant(
    v: VariantInput, results: list[dict[str, Any]], assembly: str, genes: set[str] | None, include_noncoding: bool
) -> list[Row]:
    """Turn VEP's results for one input into table rows: one per transcript consequence kept.

    Args:
        v: the input variant.
        results: VEP's results for it (one per genomic allele).
        assembly: genome build, recorded in each row.
        genes: keep only consequences in these genes (upper case); None keeps all.
        include_noncoding: also keep consequences outside DETAIL_CONSEQUENCES, and report results
            without any transcript.

    Returns:
        Rows with the private keys ``_terms`` (consequence terms), ``_protein_end``, ``_spdi``, ``_vcvs``
        and ``_clinvar_vep``, which the later stages use and remove.
    """
    rows: list[Row] = []
    for r in results:
        strand = r.get("strand", 1)
        alleles = (r.get("allele_string") or "").split("/")
        g_ref = alleles[0] if alleles else None
        if strand == -1 and g_ref:
            g_ref = revcomp(g_ref)
        tcs = r.get("transcript_consequences") or []
        coloc = r.get("colocated_variants") or []
        base: Row = {
            "input": v.raw,
            "source": v.source,
            "sample": v.sample,
            "zygosity": v.zygosity,
            "genotype": v.genotype,
            "carriers": v.carriers,
            "vcf_id": v.vcf_id,
            "vcf_filter": v.vcf_filter,
            "assembly": assembly,
            "variant_class": r.get("variant_class"),
            "chrom": v.chrom or r.get("seq_region_name"),
            "pos": v.pos or r.get("start"),
            "ref": v.ref or g_ref,
        }
        kept = []
        for tc in tcs:
            terms = set(tc.get("consequence_terms") or [])
            sym = tc.get("gene_symbol") or tc.get("gene_id")
            if genes and (sym or "").upper() not in genes:
                continue
            if not include_noncoding and not terms & DETAIL_CONSEQUENCES:
                continue
            kept.append((tc, terms))
        if not kept and include_noncoding and not tcs and not genes:
            row = dict(base)
            row["alt"] = v.alt or (alleles[1] if len(alleles) > 1 else None)
            row["consequence"] = r.get("most_severe_consequence")
            row.update(colocated_info(coloc, row["alt"]))
            row["_terms"] = {r.get("most_severe_consequence")}
            rows.append(row)
        snv = bool(v.ref and v.alt and len(v.ref) == 1 and len(v.alt) == 1)
        for tc, terms in kept:
            row = dict(base)
            va = tc.get("variant_allele")
            alt_plus = revcomp(va) if (va and strand == -1) else va
            if snv and alt_plus and alt_plus != v.alt:  # a multi-allelic result: keep only this allele
                continue
            row["alt"] = v.alt or alt_plus
            row.update(colocated_info(coloc, alt_plus, va))
            if r.get("start") and g_ref is not None and alt_plus is not None:
                row["_spdi"] = (
                    str(r.get("seq_region_name")),
                    int(r["start"]) - 1,
                    g_ref.replace("-", ""),
                    alt_plus.replace("-", ""),
                )
            _add_transcript(row, tc, terms)
            rows.append(row)
    return rows


def _add_transcript(row: Row, tc: dict[str, Any], terms: set[str]) -> None:
    """Add the fields of one VEP transcript consequence (HGVS, predictor scores, ...) to ``row``."""
    aa = tc.get("amino_acids") or ""
    ref_aa, _, alt_aa = aa.partition("/")
    hgvsp = tc.get("hgvsp")
    # REST nests plugin output; a local VEP may report flat keys (am_class, spliceai_pred_ds_ag, ...)
    am = tc.get("alphamissense") or {
        "am_class": tc.get("am_class"),
        "am_pathogenicity": tc.get("am_pathogenicity"),
    }
    amc: JSON = am.get("am_class")
    amc = {"pathogenic": "likely_pathogenic", "benign": "likely_benign"}.get(amc, amc)
    sp = tc.get("spliceai") or {k: tc.get(f"spliceai_pred_{k.lower()}") for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")}
    sp_vals = [x for x in (_f(sp.get(k)) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")) if x is not None]
    sp_max = max(sp_vals) if sp_vals else None
    acc = (tc.get("swissprot") or tc.get("trembl") or [None])[0]
    row.update(
        {
            "gene": tc.get("gene_symbol") or tc.get("gene_id"),
            "consequence": ",".join(tc.get("consequence_terms") or []),
            "impact": tc.get("impact"),
            "hgvsc": _strip_ref(tc.get("hgvsc")),
            "hgvsp": _strip_ref(hgvsp),
            "protein_change": short_protein_change(hgvsp)
            or (f"{ref_aa}{tc.get('protein_start')}{alt_aa}" if ref_aa and alt_aa else None),
            "transcript": tc.get("transcript_id"),
            "mane_select": tc.get("mane_select"),
            "canonical": bool(tc.get("canonical")),
            "protein_id": tc.get("protein_id"),
            "uniprot": acc.split(".")[0] if acc else None,
            "exon": tc.get("exon"),
            "intron": tc.get("intron"),
            "codons": tc.get("codons"),
            "cds_position": tc.get("cds_start"),
            "protein_position": tc.get("protein_start"),
            "_protein_end": tc.get("protein_end"),
            "ref_aa": ref_aa or None,
            "alt_aa": alt_aa or None,
            "sift": tc.get("sift_prediction"),
            "sift_score": _f(tc.get("sift_score")),
            "polyphen": tc.get("polyphen_prediction"),
            "polyphen_score": _f(tc.get("polyphen_score")),
            "alphamissense": _f(am.get("am_pathogenicity")),
            "alphamissense_class": amc,
            "revel": _f(tc.get("revel")) if tc.get("revel") is not None else _max_listed(tc.get("revel_score")),
            "cadd_phred": _f(tc.get("cadd_phred")),
            "eve": _f(tc.get("eve_score")),
            "eve_class": tc.get("eve_class"),
            "popeve": _f(tc.get("popeve_score")),
            "spliceai_max": sp_max,
            "_terms": terms,
        }
    )
    row["nmd"] = _nmd(row["exon"], terms)
    if "missense_variant" in terms and len(ref_aa) == 1 and len(alt_aa) == 1:
        row.update(substitution_profile(ref_aa, alt_aa))


# --------------------------------------------------------------------------- genome build


def _ref_bases(
    http: Http, assembly: str, variants: list[VariantInput], cache: KeyValueCache | None = None
) -> dict[str, str]:
    """Return query -> genome sequence at the REF position(s), from the Ensembl REST API.

    Only nuclear variants with a REF of 1-50 bases are checked. Answers are cached; a failed request
    is logged and its variants are left out.
    """
    cache = cache or NullCache()
    vs = [
        (v, v.pos, v.ref)
        for v in variants
        if v.chrom in CHROMS and v.pos is not None and v.ref and v.ref != "-" and len(v.ref) <= 50
    ]
    key = {v.query: f"refseq|{assembly}|{v.chrom}:{pos}:{len(ref)}" for v, pos, ref in vs}
    hit = cache.get_many(list(set(key.values())))
    out = {q: hit[k] for q, k in key.items() if k in hit}
    vs = [x for x in vs if x[0].query not in out]
    for i in range(0, len(vs), 50):
        chunk = vs[i : i + 50]
        regs = {f"{v.chrom}:{pos}..{pos + len(ref) - 1}:1": v for v, pos, ref in chunk}
        try:
            data = http.post_json(
                f"{ENSEMBL[assembly]}/sequence/region/human",
                {"regions": list(regs)},
                headers={"Content-Type": "application/json"},
            )
        except ApiError as e:
            log.warning("reference check failed: %s", e)
            continue
        for item in data:
            v = regs.get(item.get("query"))
            if v is not None:
                out[v.query] = (item.get("seq") or "").upper()
                cache.set(key[v.query], out[v.query])
    return out


def _match_rate(
    http: Http, assembly: str, variants: list[VariantInput], cache: KeyValueCache | None = None
) -> tuple[float | None, dict[str, str]]:
    """Return the share of checked REF alleles that match ``assembly`` (None: none checked), and the genome bases."""
    seqs = _ref_bases(http, assembly, variants, cache)
    checked = [v for v in variants if v.query in seqs]
    if not checked:
        return None, seqs
    ok = sum(seqs[v.query] == v.ref for v in checked)
    return ok / len(checked), seqs


def gene_regions(
    http: Http, cache: KeyValueCache, genes: list[str], assembly: str
) -> tuple[list[tuple[str, int, int, str]], list[str]]:
    """Look up the genomic span of each gene in Ensembl (cached).

    Args:
        http: HTTP client for the Ensembl REST API.
        cache: cache for the answers.
        genes: gene symbols.
        assembly: GRCh38 or GRCh37.

    Returns:
        ``[(chrom, start, end, gene), ...]`` for the genes found, and the gene symbols not found.

    Raises:
        ApiError: the Ensembl lookup failed.
    """
    regions: list[tuple[str, int, int, str]] = []
    missing: list[str] = []
    todo = []
    for g in genes:
        hit = cache.get(f"generegion|{assembly}|{g}")
        if hit:
            regions.append(tuple(hit))
        else:
            todo.append(g)
    if todo:
        data = http.post_json(
            f"{ENSEMBL[assembly]}/lookup/symbol/homo_sapiens",
            {"symbols": todo},
            headers={"Content-Type": "application/json"},
        )
        for g in todo:
            d = data.get(g) or data.get(g.upper())
            if d and d.get("seq_region_name"):
                reg = (d["seq_region_name"], d["start"], d["end"], g)
                regions.append(reg)
                cache.set(f"generegion|{assembly}|{g}", list(reg))
            else:
                missing.append(g)
    return regions, missing


def _peek_assembly(name: str, data: bytes | None) -> str | None:
    """Return the genome build named in a VCF header (first 2 MB of the file), or None."""
    try:
        if data is not None:
            raw = data
        else:
            with open(name, "rb") as fh:
                raw = fh.read(2_000_000)
        if raw[:2] == b"\x1f\x8b":
            with gzip.open(io.BytesIO(raw)) as gz:
                raw = gz.read(2_000_000)
        for line in raw[:2_000_000].decode("utf-8", "replace").splitlines():
            if not line.startswith("##"):
                break
            hint = _assembly_from_header(line)
            if hint:
                return hint
    except (OSError, EOFError):
        return None
    return None


# --------------------------------------------------------------------------- main entry


def _file_info(name: str, data: bytes | None) -> dict[str, str | int]:
    """Name, size and SHA-256 of an input file (for the run record)."""
    h = hashlib.sha256()
    if data is not None:
        h.update(data)
        size = len(data)
    else:
        size = os.path.getsize(name)
        with open(name, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                h.update(block)
    return {"name": os.path.basename(name), "bytes": size, "sha256": h.hexdigest()}


def _vcf_alleles(r: Row) -> bool:
    """Return True for SNVs and VCF-style indels (no '-'), which gnomAD can look up; not mitochondrial."""
    ref, alt = str(r.get("ref") or ""), str(r.get("alt") or "")
    return (
        bool(ref)
        and bool(alt)
        and "-" not in ref + alt
        and ref.isalpha()
        and alt.isalpha()
        and str(r.get("chrom") or "").upper() not in ("MT", "M", "CHRM")
    )


def _source_warnings(status: SourceStatus) -> list[str]:
    """Return one warning per web service that failed, telling the user how to retry."""
    out = []
    for src, msgs in status.errors.items():
        n = len(status.failed_items.get(src, ()))
        out.append(
            f"{src} did not respond{f' for {n} item(s)' if n else ''} ({msgs[0][:160]}). Affected rows are "
            "marked in 'notes'. Run the analysis again to retry: everything already obtained is cached."
        )
    return out


@dataclass
class _Context:
    """State shared by the stages of one :func:`predict` run.

    Attributes:
        opt: the run's options.
        text: free-text input.
        files: input files as (name, bytes or None to read the file from disk).
        genes: gene filter, upper case (empty: no filter).
        cache: cache for the web answers.
        http: HTTP client (closed by :func:`predict`).
        status: which sources were used and which requests failed.
        say: progress callback.
        local: annotate with the local VEP.
        net: use the internet (turned off when Ensembl cannot be reached in local-VEP mode).
        t0: start time (``time.time()``), for the duration.
        started: start time, UTC ISO 8601, for the run record.
        assembly: genome build in use (may change after the REF allele check).
        hint: genome build named in a VCF header, if any.
        release: Ensembl release of the REST server, if known.
        cache_only: Ensembl could not be reached and the cache held earlier answers: the run uses them only.
        local_version: local VEP version, if known.
        n_unique: unique VEP queries, for the statistics.
        n_no_protein: valid inputs that gave no row, for the statistics.
        warnings: messages for the user, in the order they arose.
        data_versions: releases (or retrieval dates) of the downloaded data sets, for the run record.
    """

    opt: PredictOptions
    text: str
    files: list[tuple[str, bytes | None]]
    genes: list[str]
    cache: KeyValueCache
    http: Http
    status: SourceStatus
    say: Progress
    local: bool
    net: bool
    t0: float
    started: str
    assembly: str = "GRCh38"
    hint: str | None = None
    release: int | None = None
    cache_only: bool = False
    local_version: str | None = None
    n_unique: int = 0
    n_no_protein: int = 0
    warnings: list[str] = field(default_factory=list)
    data_versions: dict[str, str] = field(default_factory=dict)


@dataclass
class _Lookups:
    """Answers of the lookup stages that the last, per-row stage combines.

    Attributes:
        windows: row index -> genomic window (chrom, start, end) searched in ClinVar for the same residue.
        residue_refs: row index -> ClinVar variants at the same residue, for the rows whose search succeeded.
        gnomad_need: row index -> (chrom, pos, ref, alt) of the rows without a frequency in Ensembl's data.
        gnomad_info: row index -> gnomAD's own answer, for the rows of ``gnomad_need`` it answered.
        uniprot: UniProt accession -> UniProt entry.
        alphafold: UniProt accession -> AlphaFold info.
        protein_seqs: Ensembl protein ID (without version) -> reference protein sequence.
        onco_data: hotspots and gene roles for the oncogenicity points (None: not asked for or offline).
        onco_refs: row index -> ClinVar oncogenic variants at the same residue, for the rows searched.
        onco_structures: transcript id -> exons and protein length, for the null variants (OVS1 strength).
    """

    windows: dict[int, tuple[str, int, int]] = field(default_factory=dict)
    residue_refs: dict[int, list[residue.Reference]] = field(default_factory=dict)
    gnomad_need: dict[int, tuple[str, int, str, str]] = field(default_factory=dict)
    gnomad_info: dict[int, dict[str, Any]] = field(default_factory=dict)
    uniprot: dict[str, dict[str, Any]] = field(default_factory=dict)
    alphafold: dict[str, dict[str, Any]] = field(default_factory=dict)
    protein_seqs: dict[str, str] = field(default_factory=dict)
    onco_data: onco.Data | None = None
    onco_refs: dict[int, list[onco.Reference]] = field(default_factory=dict)
    onco_structures: dict[str, dict[str, Any]] = field(default_factory=dict)


def predict(
    text: str = "",
    files: Iterable[str | tuple[str, bytes | None]] = (),
    options: PredictOptions | None = None,
    progress: Progress | None = None,
) -> PredictionResult:
    """Annotate variants given as free text and/or files [(name, bytes_or_None), ...].

    Stages, in order: read the inputs, check the genome build, annotate with Ensembl VEP (with the
    predictor plugins), build the rows, add ClinVar (exact allele, then the same residue for PS1/PM5),
    ask gnomAD about variants without a frequency, fetch the protein context, then classify each row
    and add its ACMG evidence. A web service that fails does not stop the run: the affected rows are
    flagged in ``notes`` and the failure is reported in ``warnings`` and ``sources``.

    Args:
        text: variants as free text, one or more per line (see :func:`missensexplorer.inputs.parse_text`).
        files: input files (VCF, VCF.gz, CSV, TSV, XLSX, TXT), each a path or a ``(name, bytes)`` pair;
            bytes None means "read the file from disk".
        options: run options (default: ``PredictOptions()``). Network calls, none of them with
            ``offline=True``: the Ensembl REST API (connectivity, genome build check, release, the gene
            coordinates of ``genes`` for files, and the VEP annotation itself with ``vep="rest"``, or
            only for HGVS/rsID inputs with ``vep="local"``); ``predictors`` adds the VEP plugins to
            those requests; ``clinvar`` queries NCBI E-utilities; ``gnomad_check`` queries the gnomAD
            API; ``protein_context`` queries UniProt and the AlphaFold DB; ``sequences`` fetches protein
            sequences from Ensembl. Answers are stored in the cache unless ``cache=False``.
        progress: called as ``progress(stage, done, total)`` as the run advances.

    Returns:
        The results table, the inputs that failed, warnings, statistics, mutant protein sequences, the
        UniProt and AlphaFold data used, and a run record (provenance).

    Raises:
        ValueError: ``offline`` is set without the local VEP.
        ConnectionError: the Ensembl REST API cannot be reached and the REST VEP is in use.
    """
    t0 = time.time()
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    opt = options or PredictOptions()
    local = opt.vep == "local"
    if opt.offline and not local:
        raise ValueError("offline mode needs the local VEP (vep='local')")
    cache: Cache | NullCache = Cache(opt.cache_path) if opt.cache else NullCache()
    http = Http()
    say: Progress = progress or (lambda *a, **k: None)
    try:
        ctx = _Context(
            opt=opt,
            text=text,
            genes=[g.strip().upper() for g in (opt.genes or []) if g and g.strip()],
            files=[(f, None) if isinstance(f, str) else f for f in files],
            cache=cache,
            http=http,
            status=SourceStatus(),
            say=say,
            local=local,
            net=not opt.offline,
            t0=t0,
            started=started,
        )
        # ---- genome build (first guess: explicit > VCF header > GRCh38)
        say("Reading input", 0, 1)
        ctx.hint = _files_assembly(ctx.files)
        ctx.assembly = opt.assembly if opt.assembly != "auto" else (ctx.hint or "GRCh38")
        _check_connection(ctx)
        parsed = _load_inputs(ctx)
        say("Reading input", 1, 1)
        parsed, ref_note = _check_genome_build(ctx, parsed)
        ctx.release = _release(ctx)

        # ---- VEP (local and/or REST), then one row per variant and transcript
        queries = _vep_queries(parsed)
        ctx.n_unique = sum(len(set(q)) for q in queries.values())
        answers, no_predictors = _annotate(ctx, queries)
        _normalise_indels(ctx, answers, no_predictors)
        rows, errors = _build_rows(ctx, parsed, answers, no_predictors, ref_note)

        # ---- ClinVar, gnomAD and protein context, then the verdict per row
        _add_clinvar(ctx, rows)
        found = _Lookups()
        found.windows, found.residue_refs = _search_same_residue(ctx, rows)
        found.onco_data, found.onco_refs, found.onco_structures = _oncogenicity_lookups(ctx, rows, found.windows)
        if found.onco_data is not None:
            ctx.data_versions.update(found.onco_data.versions)
        found.gnomad_need, found.gnomad_info = _check_gnomad(ctx, rows)
        found.uniprot, found.alphafold, found.protein_seqs = _fetch_protein_context(ctx, rows)
        sequences = _finish_rows(ctx, rows, found)
        return _result(ctx, parsed, rows, errors, sequences, found)
    finally:
        http.close()
        cache.close()


# --------------------------------------------------------------------------- pipeline stages


def _files_assembly(files: list[tuple[str, bytes | None]]) -> str | None:
    """Return the genome build named in the first VCF header that names one, or None."""
    hint = None
    for name, data in files:
        hint = hint or _peek_assembly(name, data)
    return hint


def _check_connection(ctx: _Context) -> None:
    """Check that Ensembl is reachable; in local-VEP mode, carry on offline when it is not.

    Raises:
        ConnectionError: Ensembl cannot be reached and the REST VEP is in use.
    """
    if not ctx.net:
        return
    # one quick request tells us whether Ensembl is reachable at all
    try:
        ensembl_release(ctx.http, ctx.cache, "GRCh38")
    except OfflineError as e:
        if not ctx.local and any(ctx.cache.get(f"ensembl_release|{a}", max_age=0) for a in ("GRCh38", "GRCh37")):
            # a cache filled by earlier runs (a rerun of a published analysis): use what it holds, whatever its age
            ctx.cache.freeze()
            ctx.cache_only = True
            ctx.warnings.append("Ensembl cannot be reached: only answers already in the cache are used")
            return
        if not ctx.local:
            raise ConnectionError(
                "Cannot reach the Ensembl REST API (rest.ensembl.org). Check the internet connection, "
                "or use the local VEP mode for VCF/genomic input."
            ) from e
        ctx.net = False
        ctx.warnings.append("No internet connection: using local VEP only; ClinVar, UniProt and AlphaFold are skipped.")
    except ApiError as e:
        log.warning("Ensembl release lookup failed: %s", e)


def _load_inputs(ctx: _Context) -> ParsedInput:
    """Parse the text and read the files; with a gene filter, VCF lines outside the genes are skipped."""
    regions = None
    if ctx.genes and ctx.files:
        if ctx.net:
            try:
                regions, missing = gene_regions(ctx.http, ctx.cache, ctx.genes, ctx.assembly)
                if missing:
                    ctx.warnings.append(f"Gene(s) not found in Ensembl: {', '.join(missing)}")
            except ApiError as e:
                ctx.status.fail("Ensembl gene lookup", str(e), ctx.genes)
        if regions is None:
            ctx.warnings.append(
                "Gene coordinates unavailable: the gene filter is applied after annotation (slower for large files)."
            )
    opt = ctx.opt
    parsed = ParsedInput()
    if ctx.text and ctx.text.strip():
        parsed.extend(parse_text(ctx.text, gene=opt.default_gene))
    for name, data in ctx.files:
        parsed.extend(
            read_file(name, data, gene=opt.default_gene, sample=opt.sample, regions=regions, pass_only=opt.pass_only)
        )
    return parsed


def _region_variants(parsed: ParsedInput) -> list[VariantInput]:
    """Return the valid inputs given as genomic coordinates."""
    return [v for v in parsed.variants if v.ok and v.kind == "region"]


def _check_genome_build(ctx: _Context, parsed: ParsedInput) -> tuple[ParsedInput, dict[str, str]]:
    """Verify REF alleles against the genome; auto-switch build if clearly wrong.

    All genomic-coordinate inputs are checked, or a sample of 50 when there are more than 200. The
    build is switched only when it was not chosen (option or VCF header) and the other build matches
    better and at least 90% of the checked REF alleles.

    Returns:
        The inputs (read again when a switch changes the gene-filter coordinates), and query -> note for
        each REF allele that does not match the genome.
    """
    region_vars = _region_variants(parsed)
    ref_seqs: dict[str, str] = {}
    if region_vars and ctx.net:
        ctx.say("Checking genome build", 0, 1)
        small = len(region_vars) <= 200
        sample = region_vars if small else region_vars[:: max(1, len(region_vars) // 50)][:50]
        try:
            rate, ref_seqs = _match_rate(ctx.http, ctx.assembly, sample, ctx.cache)
            if rate is not None and rate < 0.9:
                other = "GRCh37" if ctx.assembly == "GRCh38" else "GRCh38"
                rate2, seqs2 = _match_rate(ctx.http, other, sample, ctx.cache)
                if ctx.opt.assembly == "auto" and not ctx.hint and rate2 is not None and rate2 > rate and rate2 >= 0.9:
                    ctx.warnings.append(
                        f"REF alleles match {other} ({rate2:.0%}) better than {ctx.assembly} ({rate:.0%}); "
                        f"using {other}."
                    )
                    ctx.assembly, ref_seqs = other, seqs2
                    if ctx.genes and ctx.files:
                        parsed = _load_inputs(ctx)
                        region_vars = _region_variants(parsed)
                else:
                    ctx.warnings.append(
                        f"Only {rate:.0%} of checked REF alleles match {ctx.assembly}"
                        + (f" ({other}: {rate2:.0%})" if rate2 is not None else "")
                        + ". Check the genome build and coordinates."
                    )
        except ApiError as e:
            ctx.status.fail("Genome build check", str(e))
        ctx.say("Checking genome build", 1, 1)
    elif region_vars and ctx.opt.assembly == "auto" and not ctx.hint:
        ctx.warnings.append(f"Genome build not verified (offline); assuming {ctx.assembly}.")
    ref_note = {}
    for v in region_vars:
        g = ref_seqs.get(v.query)
        if g and g != v.ref:
            ref_note[v.query] = f"REF {v.ref} does not match {ctx.assembly} ({g})"
    return parsed, ref_note


def _release(ctx: _Context) -> int | None:
    """Return the Ensembl release of the REST server for the build in use; None offline or on failure."""
    if not ctx.net:
        return None
    if ctx.cache_only:  # the release of the cached answers (their keys name it)
        rel = ctx.cache.get(f"ensembl_release|{ctx.assembly}", max_age=0)
        return int(rel) if rel else None
    try:
        return ensembl_release(ctx.http, ctx.cache, ctx.assembly)
    except ApiError as e:
        log.warning("Ensembl release lookup failed: %s", e)
        return None


def _vep_queries(parsed: ParsedInput) -> dict[str, list[str]]:
    """Return the VEP queries of the valid inputs, grouped by kind (region, hgvs, id)."""
    queries: dict[str, list[str]] = {}
    for v in parsed.variants:
        if v.ok:
            queries.setdefault(v.kind, []).append(v.query)
    return queries


def _annotate(
    ctx: _Context, queries: dict[str, list[str]]
) -> tuple[dict[tuple[str, str], dict[str, Any]], set[tuple[str, str]]]:
    """Annotate the queries with the local VEP and/or the Ensembl VEP REST API.

    The local VEP takes the genomic-coordinate queries; the REST API the others (all of them in REST mode).
    Offline, the local VEP takes everything and reports the queries it cannot annotate as errors.

    Returns:
        (kind, query) -> VEP answer (``{"results": [...]}`` or ``{"error": ...}``), and the (kind, query)
        pairs annotated without predictor scores because the predictor plugins failed.
    """
    answers: dict[tuple[str, str], dict[str, Any]] = {}
    no_predictors: set[tuple[str, str]] = set()
    if ctx.local:
        lv = LocalVep(ctx.cache, ctx.assembly, ctx.opt.vep_command, ctx.opt.vep_args, status=ctx.status)
        ctx.local_version = lv.version()
        answers.update(
            lv.annotate({"region": queries.get("region", [])}, None, ctx.say, "Annotating variants (local VEP)")
        )
        rest_queries = {k: v for k, v in queries.items() if k != "region"}
        if rest_queries and not ctx.net:
            answers.update(lv.annotate(rest_queries))
            rest_queries = {}
    else:
        rest_queries = queries
    if rest_queries:
        vep = Vep(ctx.http, ctx.cache, ctx.assembly, workers=ctx.opt.workers, release=ctx.release, status=ctx.status)
        n_rest = sum(len(set(q)) for q in rest_queries.values())
        if not ctx.opt.predictors:
            answers.update(vep.annotate(rest_queries, BASE_PARAMS, ctx.say, "Annotating variants (VEP)"))
        elif n_rest <= 200:
            part, no_predictors = _annotate_one_pass(ctx, vep, rest_queries)
            answers.update(part)
        else:
            part, no_predictors = _annotate_two_passes(ctx, vep, rest_queries)
            answers.update(part)
    return answers, no_predictors


def _annotate_one_pass(
    ctx: _Context, vep: Vep, queries: dict[str, list[str]]
) -> tuple[dict[tuple[str, str], dict[str, Any]], set[tuple[str, str]]]:
    """Annotate with the predictor plugins; retry without them the queries whose request failed."""
    full = {**BASE_PARAMS, **PREDICTOR_PARAMS}
    no_predictors: set[tuple[str, str]] = set()
    part = vep.annotate(queries, full, ctx.say, "Annotating + predictors (VEP)", source="Ensembl VEP + predictors")
    retry: dict[str, list[str]] = {}
    for (kind, q), val in part.items():
        if val.get("service_error"):
            retry.setdefault(kind, []).append(q)
    if retry:  # plugins are the usual culprit: try again without them
        base = vep.annotate(retry, BASE_PARAMS, ctx.say, "Retrying without predictors")
        for k, val in base.items():
            if "results" in val:
                part[k] = val
                no_predictors.add(k)
    return part, no_predictors


def _annotate_two_passes(
    ctx: _Context, vep: Vep, queries: dict[str, list[str]]
) -> tuple[dict[tuple[str, str], dict[str, Any]], set[tuple[str, str]]]:
    """Annotate without the predictor plugins, then score only the queries with a coding consequence.

    Used for more than 200 queries. A query whose scoring request fails keeps its base annotation.
    """
    full = {**BASE_PARAMS, **PREDICTOR_PARAMS}
    no_predictors: set[tuple[str, str]] = set()
    part = vep.annotate(queries, BASE_PARAMS, ctx.say, "Annotating variants (VEP, fast pass)")
    detail: dict[str, list[str]] = {}
    for (kind, q), val in part.items():
        if any(_needs_predictors(r, ctx.genes) for r in val.get("results", [])):
            detail.setdefault(kind, []).append(q)
    if detail:
        scored = vep.annotate(
            detail, full, ctx.say, "Scoring coding variants (predictors)", source="VEP predictor plugins"
        )
        for k, val in scored.items():
            if "results" in val:
                part[k] = val
            else:  # keep the base annotation instead of losing the variant
                no_predictors.add(k)
    return part, no_predictors


OUTSIDE_CODING = (
    "Ensembl places this protein change outside the coding sequence of the gene's main transcript "
    "(check the position and the isoform)"
)


def _protein_change_of_kept_gene(v: VariantInput, gene_set: set[str] | None) -> bool:
    """Whether `v` is a protein change written with a gene symbol that the gene filter keeps.

    Such an input must give a protein-coding row; when it gives none, Ensembl placed it outside the
    protein, and the input is reported as an error instead of being counted among the non-coding ones.
    """
    gene, _, change = v.query.partition(":")
    if v.kind != "hgvs" or not change.startswith("p."):
        return False
    return gene_set is None or gene.upper() in {g.upper() for g in gene_set}


_PROTEIN_POS = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})?")


def _renumbered(v: VariantInput, row: Row) -> str | None:
    """Say so when a protein change was given at another position than the reported transcript's.

    Ensembl resolves a protein change on any isoform of the gene, so a change written in an older or
    another isoform's numbering (MUTYH p.Gly396Asp) is reported on the MANE transcript at its position
    there (G368D). The note keeps that change of numbering from going unnoticed.
    """
    gene, _, protein = v.query.partition(":")
    m = _PROTEIN_POS.match(protein)
    pos = row.get("protein_position")
    if v.kind != "hgvs" or not m or not pos or str(row.get("gene") or "").upper() != gene.upper():
        return None
    if int(m.group(2)) == int(pos) or _bare(row.get("hgvsp")) == _bare(protein):
        return None  # same position; or a deletion across two codons, numbered by its first codon
    shown = row.get("protein_change") or pos
    given = (AA3_TO_1.get(m.group(1).upper()), AA3_TO_1.get((m.group(3) or "").upper()))
    if m.group(3) and given != (row.get("ref_aa"), row.get("alt_aa")):  # not even the same amino acids
        return (
            f"protein change given as {protein} exists only on another isoform; on the reported transcript this "
            f"DNA change is {shown}: check the input"
        )
    return f"protein change given as {protein} (another isoform's numbering); on the reported transcript it is {shown}"


def _bare(hgvsp: Any) -> str:
    """Return a protein HGVS without its sequence id and parentheses: ``ENSP…:p.(Leu747del)`` -> ``p.Leu747del``."""
    return str(hgvsp or "").rpartition(":")[2].replace("(", "").replace(")", "")


_PROTEIN_INDEL = re.compile(r"^p\.\(?[A-Z][a-z]{2}\d+(?:_[A-Z][a-z]{2}\d+)?(?:del|dup|ins)(?!.*fs)", re.I)
INDEL_WINDOW = 100  # bases of genome fetched before an insertion or deletion, to left-align it


def _normalise_indels(
    ctx: _Context, answers: dict[tuple[str, str], dict[str, Any]], no_predictors: set[tuple[str, str]]
) -> None:
    """Re-annotate the insertions and deletions resolved from HGVS in their left-aligned VCF form.

    A deletion inside a repeat can be placed at several equivalent positions, and Ensembl resolves an
    HGVS notation to the position the notation implies. CADD's indel scores, gnomAD's frequencies and
    ClinVar's records are keyed by the left-aligned VCF form that variant callers produce, so another
    placement misses them: ``CFTR:c.1521_1523del`` got no CADD score, while ``chr7:117559590 ATCT>A``
    (the same deletion) did. Each such result is replaced, in ``answers``, by the annotation of its
    left-aligned form; the input keeps what the user wrote. Needs the Ensembl REST API (not offline).
    """
    if not ctx.net:
        return
    todo: dict[tuple[str, str], list[tuple[int, tuple[str, int, str, str]]]] = {}
    for key, val in answers.items():
        if key[0] != "hgvs":
            continue
        for j, r in enumerate(val.get("results") or []):
            loc = _indel_location(r)
            if loc:
                todo.setdefault(key, []).append((j, loc))
    if not todo:
        return
    locs = {loc for items in todo.values() for _, loc in items}
    windows = _genome_windows(ctx, locs)
    region_of: dict[tuple[str, int, str, str], str] = {}
    for loc in locs:
        norm = left_align(*loc, windows[loc]) if loc in windows else None
        if norm:
            region_of[loc] = region_query(*norm)
    if not region_of:
        return
    extra, extra_np = _annotate(ctx, {"region": sorted(set(region_of.values()))})
    for key, items in todo.items():
        results = list(answers[key]["results"])
        for j, loc in items:
            rq = region_of.get(loc)
            got = extra.get(("region", rq), {}).get("results") if rq else None
            if got:
                if on_canonical(key[1], results[j]):  # aligning moves the start, not the residues
                    answers[key] = {**answers[key], "given_on_canonical": True}
                results[j] = {**got[0], "input": key[1]}
                if ("region", rq) in extra_np:
                    no_predictors.add(key)
        answers[key] = {**answers[key], "results": results}


def _indel_location(result: dict[str, Any]) -> tuple[str, int, str, str] | None:
    """Return a VEP result's insertion or deletion on the forward strand: (chrom, start, deleted, inserted).

    ``start`` is VEP's: the first deleted base, or the base after an insertion. The empty allele is "-".
    None for anything else (substitutions, delins, non-nuclear sequences).
    """
    deleted, _, inserted = str(result.get("allele_string") or "").partition("/")
    chrom, start = str(result.get("seq_region_name") or ""), result.get("start")
    if "/" in inserted or chrom not in CHROMS or not start or (deleted == "-") == (inserted == "-"):
        return None
    if result.get("strand", 1) == -1:
        deleted, inserted = (revcomp(a) if a != "-" else a for a in (deleted, inserted))
    if not re.fullmatch(r"[ACGT]+|-", deleted) or not re.fullmatch(r"[ACGT]+|-", inserted):
        return None
    return chrom, int(start), deleted, inserted


def _genome_windows(
    ctx: _Context, locs: Iterable[tuple[str, int, str, str]]
) -> dict[tuple[str, int, str, str], tuple[int, str]]:
    """Return, per indel location, (first position, genome sequence) from INDEL_WINDOW bases before it to just after it.

    Answers are cached; a failed request is logged and its locations are left out.
    """
    want = {}
    for loc in locs:
        chrom, start, deleted, _ = loc
        length = 0 if deleted == "-" else len(deleted)
        want[loc] = (chrom, max(1, start - INDEL_WINDOW), start + length + 1)
    key = {loc: f"refwin|{ctx.assembly}|{c}:{a}-{b}" for loc, (c, a, b) in want.items()}
    hit = ctx.cache.get_many(list(set(key.values())))
    out = {loc: (want[loc][1], hit[k]) for loc, k in key.items() if k in hit}
    missing = [loc for loc in want if loc not in out]
    for i in range(0, len(missing), 50):
        chunk = {f"{c}:{a}..{b}:1": loc for loc in missing[i : i + 50] for c, a, b in [want[loc]]}
        try:
            data = ctx.http.post_json(
                f"{ENSEMBL[ctx.assembly]}/sequence/region/human",
                {"regions": list(chunk)},
                headers={"Content-Type": "application/json"},
            )
        except ApiError as e:
            log.warning("genome sequence for indel alignment failed: %s", e)
            continue
        for item in data:
            found = chunk.get(item.get("query"))
            if found is not None and item.get("seq"):
                seq = str(item["seq"]).upper()
                out[found] = (want[found][1], seq)
                ctx.cache.set(key[found], seq)
    return out


def left_align(
    chrom: str, start: int, deleted: str, inserted: str, window: tuple[int, str]
) -> tuple[str, int, str, str] | None:
    """Shift an insertion or deletion as far left as it goes, and write it in VCF form.

    Args:
        chrom: Chromosome.
        start: VEP's start: the first deleted base, or the base after the insertion point.
        deleted: Deleted bases, or "-" for an insertion.
        inserted: Inserted bases, or "-" for a deletion.
        window: (position of the first base, genome sequence) around the variant.

    Returns:
        ``(chrom, pos, ref, alt)`` with the preceding base included, as in a VCF; None when the window
        does not cover the variant or the genome does not hold the deleted bases.
    """
    first, seq = window

    def base(x: int) -> str | None:
        return seq[x - first] if 0 <= x - first < len(seq) else None

    if inserted == "-":  # deletion of `deleted` at start..start+len-1
        if seq[start - first : start - first + len(deleted)] != deleted:
            return None
        while (b := base(start - 1)) and base(start - 2) and b == deleted[-1]:
            deleted, start = b + deleted[:-1], start - 1
        anchor = base(start - 1)
        return (chrom, start - 1, anchor + deleted, anchor) if anchor else None
    at = start - 1  # insertion after this base
    while (b := base(at)) and base(at - 1) and b == inserted[-1]:
        inserted, at = b + inserted[:-1], at - 1
    anchor = base(at)
    return (chrom, at, anchor, anchor + inserted) if anchor else None


def _needs_predictors(result: dict[str, Any], genes: list[str]) -> bool:
    """Return True when a VEP result has a consequence in DETAIL_CONSEQUENCES (in ``genes``, if given)."""
    return any(
        (set(tc.get("consequence_terms") or []) & DETAIL_CONSEQUENCES)
        and (not genes or (tc.get("gene_symbol") or "").upper() in genes)
        for tc in result.get("transcript_consequences") or []
    )


def _build_rows(
    ctx: _Context,
    parsed: ParsedInput,
    answers: dict[tuple[str, str], dict[str, Any]],
    no_predictors: set[tuple[str, str]],
    ref_note: dict[str, str],
) -> tuple[list[Row], list[dict[str, str]]]:
    """Turn the VEP answers into rows, one per variant and transcript; counts the inputs without rows.

    Returns:
        The rows (with ``notes`` on REF mismatches and missing predictor scores), and one error
        ``{input, source, error}`` per input that could not be parsed or annotated.
    """
    rows: list[Row] = []
    errors: list[dict[str, str]] = []
    gene_set = set(ctx.genes) if ctx.genes else None
    for v in parsed.variants:
        if not v.ok:
            errors.append({"input": v.raw, "source": v.source, "error": v.error or "could not parse"})
            continue
        val = answers.get((v.kind, v.query), {})
        if "error" in val:
            errors.append({"input": v.raw, "source": v.source, "error": val["error"]})
            continue
        vrows = _rows_for_variant(v, val.get("results", []), ctx.assembly, gene_set, ctx.opt.include_noncoding)
        if not vrows and _protein_change_of_kept_gene(v, gene_set):
            errors.append({"input": v.raw, "source": v.source, "error": OUTSIDE_CODING})
            continue
        if not vrows:
            ctx.n_no_protein += 1
        for r in vrows:
            notes = []
            if val.get("ambiguous", 0) > 1:
                notes.append(f"{val['ambiguous']} nucleotide changes give this protein change; each is reported")
            renumbered = None if val.get("given_on_canonical") else _renumbered(v, r)
            if renumbered:
                notes.append(renumbered)
            if v.kind == "hgvs" and _PROTEIN_INDEL.match(v.query.split(":", 1)[-1]):
                notes.append(
                    f"protein notation: other DNA changes can give the same {r.get('protein_change') or 'change'}; "
                    f"this is the one Ensembl maps ({r.get('hgvsc') or 'see HGVS c.'}). Give the c. or genomic "
                    "change to analyse a specific one"
                )
            if v.query in ref_note:
                notes.append(ref_note[v.query])
            if (v.kind, v.query) in no_predictors:
                notes.append("predictors unavailable (service error; run again to retry)")
            r["notes"] = "; ".join(notes) or None
        rows.extend(vrows)
    return rows, errors


def _add_clinvar(ctx: _Context, rows: list[Row]) -> None:
    """Add the ClinVar classification of the exact allele to each row, and drop the private ClinVar keys.

    Rows whose record cannot be matched keep the significances in VEP's data, marked "(submissions)".
    """
    records: dict[str, dict[str, Any]] = {}
    if ctx.opt.clinvar and rows and ctx.net:
        ctx.say("ClinVar classification", 0, 1)
        records = cv.fetch(ctx.http, ctx.cache, {x for r in rows for x in r.get("_vcvs", [])}, status=ctx.status)
        ctx.say("ClinVar classification", 1, 1)
    elif ctx.opt.clinvar and rows:
        ctx.status.skip("ClinVar", "offline: classification taken from the VEP data (all submissions)")
    for r in rows:
        rec = cv.match(records, r.get("_vcvs", []), r.get("_spdi"), r.get("mane_select"), r.get("hgvsc"))
        if rec:
            r["clinvar"] = rec["classification"]
            r["clinvar_review"] = rec["review_status"]
            r["clinvar_stars"] = cv.stars(rec["review_status"])
            r["clinvar_conditions"] = "; ".join(rec["conditions"][:6]) or None
            r["clinvar_oncogenicity"] = rec.get("oncogenicity") or rec.get("somatic_impact")
            r["clinvar_ids"] = rec["accession"]
        elif r.get("_clinvar_vep"):
            r["clinvar"] = f"{r['_clinvar_vep']} (submissions)"
        if any(ctx.status.failed("ClinVar", x) for x in r.get("_vcvs", [])):
            r["notes"] = "; ".join(x for x in (r.get("notes"), "ClinVar unavailable (service error)") if x)
        for k in ("_vcvs", "_spdi", "_clinvar_vep"):
            r.pop(k, None)


def _search_same_residue(
    ctx: _Context, rows: list[Row]
) -> tuple[dict[int, tuple[str, int, int]], dict[int, list[residue.Reference]]]:
    """ACMG PS1/PM5: find pathogenic missense variants at the same residue (ClinVar).

    Returns:
        Row index -> genomic window, for every row PS1/PM5 apply to; and row index -> the ClinVar variants
        found in it, for the rows whose search and records were obtained.
    """
    windows = {i: w for i, r in enumerate(rows) if (w := residue.window(r))}
    # ClinVar names its variants on the MANE Select transcript: without one (GRCh37) nothing can be compared
    searched = {i: w for i, w in windows.items() if rows[i].get("mane_select")}
    residue_refs: dict[int, list[residue.Reference]] = {}
    if searched and ctx.opt.clinvar and ctx.net:
        ctx.say("ClinVar: same residue", 0, 1)
        found = cv.search_windows(ctx.http, ctx.cache, list(searched.values()), ctx.assembly, status=ctx.status)
        recs = cv.fetch(ctx.http, ctx.cache, {v for w in found.values() for v in w}, status=ctx.status)
        for i, w in searched.items():
            if w in found and all(v in recs for v in found[w]):
                residue_refs[i] = residue.references(rows[i], {v: recs[v] for v in found[w]}, ctx.assembly)
        ctx.say("ClinVar: same residue", 1, 1)
    return windows, residue_refs


def _oncogenicity_lookups(
    ctx: _Context, rows: list[Row], windows: dict[int, tuple[str, int, int]]
) -> tuple[onco.Data | None, dict[int, list[onco.Reference]], dict[str, dict[str, Any]]]:
    """Somatic oncogenicity: the hotspots and gene roles, and the ClinVar oncogenic variants at each residue.

    Returns:
        The data (None when the points were not asked for, or offline), row index -> ClinVar oncogenic
        missense variants at its residue, for the rows whose search and records were obtained, and the
        exons and protein length of the transcripts of null variants.
    """
    if not ctx.opt.oncogenicity or not rows:
        return None, {}, {}
    if not ctx.net:
        for src in ("cancerhotspots.org", "IntOGen", "CIViC"):
            ctx.status.skip(src, "offline: oncogenicity points not computed")
        return None, {}, {}
    ctx.say("Oncogenicity (hotspots, IntOGen, CIViC)", 0, 1)
    data = onco.load(ctx.http, ctx.cache, ctx.status)
    refs: dict[int, list[onco.Reference]] = {}
    searched = {i: w for i, w in windows.items() if rows[i].get("mane_select")}
    if searched and ctx.opt.clinvar:
        found = cv.search_oncogenicity_windows(ctx.http, ctx.cache, list(searched.values()), ctx.assembly, ctx.status)
        recs = cv.fetch(ctx.http, ctx.cache, {v for w in found.values() for v in w}, status=ctx.status)
        for i, w in searched.items():
            if w in found and all(v in recs for v in found[w]):
                refs[i] = onco.references(rows[i], {v: recs[v] for v in found[w]}, ctx.assembly)
        _own_oncogenicity(ctx, rows)
    null_tx = {str(r.get("transcript") or "") for r in rows if set(r.get("_terms") or ()) & onco.NULL_TERMS}
    structures = onco.transcript_structures(ctx.http, ctx.cache, null_tx, ctx.assembly, ctx.status) if null_tx else {}
    ctx.say("Oncogenicity (hotspots, IntOGen, CIViC)", 1, 1)
    return data, refs, structures


def _check_gnomad(
    ctx: _Context, rows: list[Row]
) -> tuple[dict[int, tuple[str, int, str, str]], dict[int, dict[str, Any]]]:
    """Ask gnomAD itself about variants without a frequency in Ensembl's data (singletons, PM2).

    Rows found in gnomAD get its frequencies.

    Returns:
        Row index -> (chrom, pos, ref, alt) of the rows that need the check, and row index -> gnomAD's answer
        for those it answered.
    """
    need = {
        i: (str(r["chrom"]), int(r["pos"]), str(r["ref"]), str(r["alt"]))
        for i, r in enumerate(rows)
        if r.get("gnomad_af") is None and r.get("pos") and _vcf_alleles(r)
    }
    # BA1 candidates: Ensembl gives the frequency but not the allele number, and BA1 needs a group of at
    # least 2,000 alleles (Ghosh et al. 2018), so gnomAD itself is asked for these too.
    common = {
        i: (str(r["chrom"]), int(r["pos"]), str(r["ref"]), str(r["alt"]))
        for i, r in enumerate(rows)
        if i not in need and (_num(r.get("gnomad_popmax_af")) or 0) >= acmg.BA1_AF and r.get("pos") and _vcf_alleles(r)
    }
    gnomad_info: dict[int, dict[str, Any]] = {}
    if (need or common) and ctx.opt.gnomad_check and ctx.net:
        ctx.say("gnomAD", 0, 1)
        got = gnomad.lookup(
            ctx.http,
            ctx.cache,
            [*need.values(), *common.values()],
            ctx.assembly,
            status=ctx.status,
            progress=lambda k, n: ctx.say("gnomAD", k, n),
        )
        for i, k in need.items():
            if k in got:
                gnomad_info[i] = got[k]
                if got[k]["found"]:
                    rows[i]["gnomad_af"], rows[i]["gnomad_popmax_af"] = got[k]["af"], got[k]["grpmax"]
        for i, k in common.items():
            if k in got and got[k].get("grpmax_2000") is not None:
                rows[i]["_grpmax_2000"] = (got[k]["grpmax_2000"], got[k].get("grpmax_2000_an"))
        ctx.say("gnomAD", 1, 1)
    elif need and ctx.opt.gnomad_check:
        ctx.status.skip("gnomAD", "offline")
    return need, gnomad_info


def _fetch_protein_context(
    ctx: _Context, rows: list[Row]
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]], dict[str, str]]:
    """Fetch UniProt entries, AlphaFold models and Ensembl protein sequences for the protein-level rows.

    Returns:
        UniProt accession -> entry, UniProt accession -> AlphaFold info, and Ensembl protein ID (without
        version) -> sequence; each empty when its option is off or the run is offline.
    """
    uni: dict[str, dict[str, Any]] = {}
    afd: dict[str, dict[str, Any]] = {}
    seqs: dict[str, str] = {}
    prot_rows = [r for r in rows if r.get("protein_position")]
    if prot_rows and ctx.opt.protein_context and ctx.net:
        ctx.say("Protein context (UniProt, AlphaFold)", 0, 1)
        # Ensembl's GRCh37 server gives UniProt entry names (BRAF_HUMAN) where GRCh38 gives accessions
        names = prot.uniprot_accessions(ctx.http, ctx.cache, {r.get("uniprot") for r in prot_rows}, status=ctx.status)
        for r in prot_rows:
            if r.get("uniprot") in names:
                r["uniprot"] = names[r["uniprot"]]
        accs = {r["uniprot"] for r in prot_rows if r.get("uniprot")}
        uni = prot.fetch_uniprot(ctx.http, ctx.cache, accs, workers=ctx.opt.workers, status=ctx.status)
        afd = prot.fetch_alphafold(ctx.http, ctx.cache, accs, workers=ctx.opt.workers, status=ctx.status)
        ctx.say("Protein context (UniProt, AlphaFold)", 1, 1)
    elif prot_rows and ctx.opt.protein_context:
        ctx.status.skip("UniProt", "offline")
        ctx.status.skip("AlphaFold DB", "offline")
    if prot_rows and ctx.opt.sequences and ctx.net:
        seqs = prot.fetch_protein_sequences(
            ctx.http, ctx.cache, {r.get("protein_id") for r in prot_rows}, ctx.assembly, status=ctx.status
        )
    return uni, afd, seqs


def _finish_rows(ctx: _Context, rows: list[Row], found: _Lookups) -> dict[int, dict[str, str]]:
    """Add the protein context, the prediction and the ACMG evidence to each row; drop the private keys.

    Returns:
        Row index -> {wt, mut, header, note}, for the rows whose mutant protein could be built.
    """
    sequences: dict[int, dict[str, str]] = {}
    for i, r in enumerate(rows):
        terms: set[str] = r.pop("_terms", set()) or set()
        pend: int | None = r.pop("_protein_end", None)
        notes = [r["notes"]] if r.get("notes") else []
        pos: int | None = r.get("protein_position")
        pid = (r.get("protein_id") or "").split(".")[0]
        wt = found.protein_seqs.get(pid)
        if pid and not wt and ctx.status.failed("Ensembl protein sequences", pid):
            notes.append("protein sequence unavailable (service error)")
        if wt:
            r["protein_length"] = len(wt)
        if pos and wt and r.get("ref_aa"):
            seq = _add_mutant_protein(r, wt, pos, pend, terms)
            if seq is not None:
                sequences[i] = seq
        _add_protein_features(ctx, r, pos, pend, wt, notes, found)
        _add_verdict(ctx, i, r, terms, notes, found)
        r["notes"] = "; ".join(dict.fromkeys(notes)) or None
    return sequences


def _add_mutant_protein(r: Row, wt: str, pos: int, pend: int | None, terms: set[str]) -> dict[str, str] | None:
    """Add the sequence window and the mutant-protein note to ``r``; return its sequences entry, if any."""
    r["protein_window"] = prot.context_window(wt, pos, r["ref_aa"], r.get("alt_aa"))
    mut, note = prot.mutate_protein(wt, pos, pend or pos, f"{r['ref_aa']}/{r.get('alt_aa') or '-'}", ",".join(terms))
    r["mutant_protein"] = note
    if note.startswith("truncated"):
        r["_truncation"] = note
    if mut is None:
        return None
    return {
        "wt": wt,
        "mut": mut,
        "note": note,
        "header": f"{r.get('gene')}|{r.get('protein_change')}|{r.get('protein_id')}|{r.get('input')}",
    }


def _add_protein_features(
    ctx: _Context, r: Row, pos: int | None, pend: int | None, wt: str | None, notes: list[str], found: _Lookups
) -> None:
    """Add the UniProt features and AlphaFold pLDDT at the variant, when UniProt's numbering fits the transcript."""
    acc = r.get("uniprot") or ""
    if acc and pos and ctx.status.failed("UniProt", acc):
        notes.append("UniProt unavailable (service error)")
    if acc and pos and ctx.status.failed("AlphaFold DB", acc):
        notes.append("AlphaFold unavailable (service error)")
    entry = found.uniprot.get(acc)
    if not (entry and pos):
        return
    r["protein_name"] = entry.get("name")
    r["protein_function"] = entry.get("function")
    r["protein_length"] = r.get("protein_length") or entry.get("length")
    useq = entry.get("sequence") or ""
    same = bool(wt) and useq == wt
    local_ok = (
        r.get("ref_aa")
        and useq[pos - 1 : pos - 1 + len(r["ref_aa"])] == r["ref_aa"]
        and (not wt or useq[max(0, pos - 6) : pos + 5] == wt[max(0, pos - 6) : pos + 5])
    )
    if same or local_ok:
        r.update(prot.features_at(entry, pos, pend or pos))
        af = found.alphafold.get(r["uniprot"])
        r["plddt"] = prot.plddt_at(af, pos)
        r["plddt_class"] = prot.plddt_class(r["plddt"])
        if af and af.get("available"):
            r["alphafold_model"] = af.get("entry_id")
        if not same:
            notes.append("UniProt canonical isoform differs from transcript; features mapped by local sequence match")
    else:
        notes.append("UniProt canonical numbering differs from this transcript; UniProt/AlphaFold features skipped")


def _add_verdict(ctx: _Context, i: int, r: Row, terms: set[str], notes: list[str], found: _Lookups) -> None:
    """Add the prediction, the ACMG evidence and the population frequency to row ``i``."""
    label, score, basis = classify(r, terms)
    r.pop("_truncation", None)
    if "missense_variant" in terms:
        r["acmg_pp3_bp4"], r["acmg_points"] = acmg.computational(r, ctx.opt.acmg_source)
    if i in found.windows:
        r["acmg_ps1_pm5"] = _ps1_pm5(ctx, i, r, found)
        chrom, start, end = found.windows[i]
        if ctx.status.failed("ClinVar", f"{ctx.assembly}:{chrom}:{start}-{end}"):
            notes.append("ClinVar same-residue search unavailable (service error)")
    r["_gnomad"] = found.gnomad_info.get(i)
    r["acmg_pm2"] = acmg.pm2(r, r["_gnomad"])[0]
    if i in found.gnomad_need and i not in found.gnomad_info and ctx.opt.gnomad_check and ctx.net:
        notes.append("gnomAD check unavailable (service error)")
    r["population_frequency"] = acmg.frequency(r)
    if str(r["population_frequency"] or "").startswith("common: BA1") and (t := vcep.ba1_above_5pct(r.get("gene"))):
        # the gene's expert panels set BA1 above 5% (Ghosh et al.: BA1 at 5% only without a gene-specific rule)
        g = acmg.ba1_grpmax(r) or 0.0
        if not g > t:
            r["population_frequency"] = (
                f"frequent: BS1 candidate (grpmax {acmg._pct(g)}); BA1 not met: the gene's expert panel sets BA1 "
                f"above {100 * t:.3g}%"
            )
    r.pop("_gnomad")
    if "missense_variant" in terms and (ms := model.score(r)) is not None:
        r["score_with_frequency"] = model.with_frequency(ms.p, r.get("gnomad_popmax_af"))
    r["acmg_codes"] = acmg.codes(r)
    acmg_notes = [n for n in (acmg.overlap_note(r),) if n]
    if r.get("acmg_codes") or r.get("acmg_pp3_bp4"):
        r["vcep_specification"] = vcep.describe(r.get("gene"))
        acmg_notes.extend(vcep.notes(r))
    r["acmg_notes"] = vcep.SEP.join(acmg_notes) or None
    if found.onco_data is not None:
        _add_oncogenicity(r, i, terms, label, found)
    if not acmg.ba1_exception(r):  # common but pathogenic/uncertain (ClinGen list): no frequency rule
        common = ctx.opt.common_af
        if common and (t := vcep.ba1_above_5pct(r.get("gene"))) and t > common:
            common = t  # the gene's expert panels set BA1 higher (VWF, ABCA4): a variant below it is not "common"
        label, basis = apply_frequency(label, basis, acmg.ba1_grpmax(r), common)
    r.pop("_grpmax_2000", None)
    r["prediction"], r["effect_score"], r["prediction_basis"] = label, score, basis
    r["spliceai_effect"] = (
        "likely splice-altering"
        if (r.get("spliceai_max") or 0) >= 0.5
        else "possible"
        if (r.get("spliceai_max") or 0) >= 0.2
        else None
    )


def _own_oncogenicity(ctx: _Context, rows: list[Row]) -> None:
    """Fill each row's own ClinVar oncogenicity classification by searching its position.

    VEP links a variant to ClinVar only through the variation database's identifiers, which miss many
    recent somatic records; the oncogenicity search at the variant's own position finds them. Shown
    only: the variant's own classification is never evidence for itself.
    """
    need: dict[int, tuple[str, int, int]] = {}
    for i, r in enumerate(rows):
        if r.get("clinvar_oncogenicity") or not r.get("pos") or not _vcf_alleles(r):
            continue
        pos = int(r["pos"])
        need[i] = (str(r["chrom"]), pos, pos + max(len(str(r.get("ref") or "")), 1) - 1)
    if not need:
        return
    found = cv.search_oncogenicity_windows(ctx.http, ctx.cache, sorted(set(need.values())), ctx.assembly, ctx.status)
    recs = cv.fetch(ctx.http, ctx.cache, {v for w in found.values() for v in w}, status=ctx.status)
    for i, w in need.items():
        r = rows[i]
        own = (str(r["chrom"]), int(r["pos"]), str(r["ref"]), str(r["alt"]))
        cands = [v for v in found.get(w, []) if v in recs]
        hit = next((recs[v] for v in cands if residue._same_allele(recs[v], own, ctx.assembly)), None)
        hit = hit or cv.match(recs, cands, None, r.get("mane_select"), str(r.get("hgvsc") or "").split(":")[-1])
        if hit and (hit.get("oncogenicity") or hit.get("somatic_impact")):
            r["clinvar_oncogenicity"] = hit.get("oncogenicity") or hit.get("somatic_impact")


def _add_oncogenicity(r: Row, i: int, terms: set[str], label: str, found: _Lookups) -> None:
    """Add the somatic oncogenicity class, points, criteria and evidence to row ``i``."""
    data = found.onco_data
    if data is None:
        return
    if not (data.hotspots_ok and data.roles_ok):
        r["oncogenicity"] = "not assessed (cancerhotspots.org or IntOGen unavailable)"
        return
    # OP1/SBP1: the calibrated computational evidence (PP3 or BP4 at any strength), as the standard asks
    # that "all utilized lines of computational evidence" agree, splicing impact included
    pts = r.get("acmg_points")
    call = "damaging" if isinstance(pts, int) and pts > 0 else "benign" if isinstance(pts, int) and pts < 0 else None
    refs = found.onco_refs.get(i) if i in found.windows else []
    structure = found.onco_structures.get(str(r.get("transcript") or "").split(".")[0])
    a = onco.assess(r, terms, call, data, refs, found.gnomad_info.get(i), structure)
    r["oncogenicity"] = onco.tier(a.points)
    r["oncogenicity_points"] = a.points
    r["oncogenicity_codes"] = a.summary()
    r["oncogenicity_evidence"] = a.evidence()
    if data.civic_ok:
        r["clinical_evidence_level"], r["clinical_evidence"] = onco.clinical(r, data)


def _ps1_pm5(ctx: _Context, i: int, r: Row, found: _Lookups) -> str | None:
    """Return the PS1/PM5 evidence of row ``i``, or why the same residue was not checked."""
    if i in found.residue_refs:
        return residue.evidence(r, found.residue_refs[i])[0]
    if not r.get("mane_select"):
        return "not checked (no MANE Select transcript, on which ClinVar names its variants)"
    reason = "ClinVar turned off" if not ctx.opt.clinvar else "offline" if not ctx.net else "ClinVar did not respond"
    return f"not checked ({reason})"


def _result(
    ctx: _Context,
    parsed: ParsedInput,
    rows: list[Row],
    errors: list[dict[str, str]],
    sequences: dict[int, dict[str, str]],
    found: _Lookups,
) -> PredictionResult:
    """Assemble the tables, statistics, warnings and run record into the result."""
    table = pd.DataFrame(rows, columns=COLUMNS) if rows else pd.DataFrame(columns=COLUMNS)
    errs = pd.DataFrame(errors, columns=["input", "source", "error"])
    stats = _stats(ctx, parsed, table, errs)
    if parsed.samples and len(parsed.samples) > 1 and not ctx.opt.sample:
        ctx.warnings.append(
            f"Multi-sample VCF ({len(parsed.samples)} samples): carriers listed per variant. "
            "Choose one sample to get its zygosity."
        )
    ctx.warnings.extend(parsed.warnings)
    ctx.warnings.extend(_source_warnings(ctx.status))
    run = _run_record(ctx, stats["seconds"])
    ctx.say("Done", 1, 1)
    return PredictionResult(
        table=table,
        errors=errs,
        assembly=ctx.assembly,
        warnings=list(dict.fromkeys(ctx.warnings)),  # a build switch can load the inputs twice
        stats=stats,
        sequences=sequences,
        alphafold=found.alphafold,
        uniprot=found.uniprot,
        run=run,
        sources=ctx.status.summary(),
    )


def _stats(ctx: _Context, parsed: ParsedInput, table: pd.DataFrame, errs: pd.DataFrame) -> dict[str, Any]:
    """Return the counts for the summary, the run time and the file counts (VCF records, skipped lines)."""
    return {
        "inputs": len(parsed.variants),
        "unique_queries": ctx.n_unique,
        "rows": len(table),
        "protein_altering": int(table["consequence"].fillna("").str.contains("|".join(PROTEIN_ALTERING)).sum())
        if len(table)
        else 0,
        "likely_damaging": int(table["prediction"].map(is_damaging).sum()) if len(table) else 0,
        "no_protein_effect": ctx.n_no_protein,
        "errors": len(errs),
        "seconds": round(time.time() - ctx.t0, 1),
        **{f"file_{k}": v for k, v in parsed.stats.items() if k.startswith(("vcf_", "skipped_"))},
    }


def _run_record(ctx: _Context, seconds: float) -> dict[str, Any]:
    """Return the provenance of the run: versions, options, inputs (with SHA-256) and sources."""
    return {
        "tool": f"MissensExplorer {__version__}",
        "started_utc": ctx.started,
        "seconds": seconds,
        "assembly": ctx.assembly,
        "annotation": (
            f"local VEP ({ctx.local_version or 'version unknown'})"
            if ctx.local
            else f"Ensembl VEP REST, release {ctx.release or 'unknown'}"
        ),
        "ensembl_release": ctx.release,
        "cache_only": ctx.cache_only,  # Ensembl could not be reached: every answer came from the cache
        "options": {k: v for k, v in asdict(ctx.opt).items() if k not in ("cache_path",)},
        "inputs": {
            "text_lines": len([x for x in (ctx.text or "").splitlines() if x.strip()]),
            "files": [_file_info(n, d) for n, d in ctx.files],
        },
        "sources": ctx.status.summary(),
        "data_versions": {
            **ctx.data_versions,
            "ClinGen CSpec (expert panel specifications)": f"accessed {vcep.VCEP['accessed']}",
        },
        "missense_model": {
            "version": model.MODEL["version"],
            "trained": model.MODEL["trained"],
            "training_variants": model.MODEL["training"]["n"],
            "cut_offs": [model.BENIGN, model.DAMAGING],
        },
        "python": platform.python_version(),
    }


def predict_text(text: str, **kwargs: Any) -> pd.DataFrame:
    """Return the results table for variants given as text, e.g. ``predict_text("BRAF V600E")``.

    Args:
        text: variants as free text (see :func:`predict`).
        **kwargs: :class:`PredictOptions` fields.
    """
    return predict(text=text, options=PredictOptions(**kwargs)).table
