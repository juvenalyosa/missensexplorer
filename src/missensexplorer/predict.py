"""End-to-end pipeline: inputs -> VEP -> predictors -> protein context -> one table."""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass, field
from urllib.parse import unquote

import pandas as pd

from . import acmg
from . import clinvar as cv
from . import gnomad, residue
from . import model
from . import protein as prot
from .aa import AA3_TO_1, substitution_profile
from .api import ENSEMBL, ApiError, Http, OfflineError, SourceStatus, ensembl_release
from .cache import Cache, NullCache
from .inputs import CHROMS, ParsedInput, VariantInput, parse_text, read_file
from .vep import BASE_PARAMS, DETAIL_CONSEQUENCES, PREDICTOR_PARAMS, Vep
from .vep_local import LocalVep

log = logging.getLogger(__name__)

LOF_TERMS = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost"}
INFRAME_TERMS = {"inframe_insertion", "inframe_deletion", "protein_altering_variant", "stop_lost"}
SPLICE_NEAR_TERMS = {"splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant",
                     "splice_polypyrimidine_tract_variant"}
PROTEIN_ALTERING = LOF_TERMS | INFRAME_TERMS | {"missense_variant"}

CLINVAR_ORDER = [
    "pathogenic", "pathogenic/likely_pathogenic", "likely_pathogenic", "risk_factor", "drug_response",
    "association", "affects", "protective", "uncertain_significance",
    "conflicting_interpretations_of_pathogenicity", "conflicting_classifications_of_pathogenicity",
    "likely_benign", "benign/likely_benign", "benign", "other", "not_provided",
]

# Calibrated thresholds (Pejaver et al., AJHG 2022) for PP3/BP4 "supporting" evidence.
REVEL_DAMAGING, REVEL_BENIGN = 0.644, 0.290
CADD_DAMAGING, CADD_BENIGN = 25.3, 22.7

COLUMNS = [
    "input", "source", "sample", "zygosity", "genotype", "carriers",
    "gene", "protein_change", "hgvsp", "hgvsc", "consequence", "impact",
    "prediction", "effect_score", "prediction_basis", "acmg_codes", "acmg_pp3_bp4", "acmg_points", "acmg_ps1_pm5",
    "acmg_pm2", "population_frequency",
    "alphamissense", "alphamissense_class", "revel", "cadd_phred", "eve", "eve_class", "popeve",
    "sift", "sift_score", "polyphen", "polyphen_score", "spliceai_max", "spliceai_effect",
    "grantham", "grantham_class", "blosum62", "hydropathy_change", "volume_change", "charge_change",
    "aa_class_change", "physchem_notes",
    "protein_name", "uniprot", "protein_position", "protein_length", "ref_aa", "alt_aa", "protein_window",
    "mutant_protein", "protein_features", "key_site", "uniprot_variants_here", "uniprot_mutagenesis_here",
    "plddt", "plddt_class", "alphafold_model", "protein_function",
    "clinvar", "clinvar_review", "clinvar_stars", "clinvar_conditions", "clinvar_oncogenicity", "clinvar_ids",
    "gnomad_af", "gnomad_popmax_af", "rsid", "cosmic", "pubmed_count", "pubmed_ids",
    "chrom", "pos", "ref", "alt", "assembly", "variant_class", "transcript", "mane_select", "canonical",
    "protein_id", "exon", "intron", "codons", "cds_position", "nmd", "vcf_id", "vcf_filter", "notes",
]


@dataclass
class PredictOptions:
    assembly: str = "auto"                  # auto | GRCh38 | GRCh37
    default_gene: str | None = None         # gene used for bare changes like "V600E"
    genes: list[str] | None = None          # keep only variants in these genes
    sample: str | None = None               # VCF sample to report (default: all carriers)
    pass_only: bool = False                 # VCF: keep FILTER=PASS only
    include_noncoding: bool = False         # also report intronic/UTR/intergenic consequences
    predictors: bool = True                 # AlphaMissense, REVEL, CADD, EVE, SpliceAI
    protein_context: bool = True            # UniProt features + AlphaFold pLDDT
    clinvar: bool = True                    # exact-allele ClinVar classification (NCBI)
    sequences: bool = True                  # mutant protein sequences
    vep: str = "rest"                       # rest (public Ensembl API) | local (VEP on this computer)
    vep_command: str | None = None          # local VEP command (default: $MISSENSEXPLORER_VEP or 'vep --cache --offline')
    vep_args: str | None = None             # extra local VEP arguments, e.g. plugin data files
    offline: bool = False                   # never use the internet (local VEP only; no ClinVar/UniProt/AlphaFold)
    common_af: float | None = 0.05          # gnomAD grpmax at/above which damaging/uncertain calls become benign (None = off)
    acmg_source: str = "model"              # PP3/BP4 from: model (own ClinGen-method calibration) | revel (ClinGen's)
    gnomad_check: bool = True               # ask gnomAD directly about variants Ensembl gives no frequency (PM2)
    workers: int = 4
    cache: bool = True
    cache_path: str | None = None


@dataclass
class PredictionResult:
    table: pd.DataFrame
    errors: pd.DataFrame
    assembly: str
    warnings: list[str] = field(default_factory=list)
    stats: dict = field(default_factory=dict)
    sequences: dict[int, dict] = field(default_factory=dict)   # row index -> {wt, mut, header, note}
    alphafold: dict[str, dict] = field(default_factory=dict)   # uniprot -> AlphaFold info
    uniprot: dict[str, dict] = field(default_factory=dict)
    run: dict = field(default_factory=dict)                     # provenance: versions, options, inputs, sources
    sources: dict[str, str] = field(default_factory=dict)       # source -> ok | failed: ... | skipped: ...


# --------------------------------------------------------------------------- small helpers

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s if s in ("-", "") else s.translate(_COMP)[::-1]


def _f(x):
    try:
        return None if x is None or x == "" else float(x)
    except (TypeError, ValueError):
        return None


def _max_listed(s) -> float | None:
    """dbNSFP returns per-transcript lists like '.,.,0.931,.'"""
    if s is None:
        return None
    if isinstance(s, (int, float)):
        return float(s)
    vals = [_f(x) for x in str(s).split(",")]
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None


def short_protein_change(hgvsp: str | None) -> str | None:
    """'p.Val600Glu' -> 'V600E', 'p.Arg213Ter' -> 'R213*', 'p.Thr125ProfsTer23' -> 'T125Pfs*23'."""
    if not hgvsp:
        return None
    s = unquote(hgvsp).split(":")[-1]
    s = s.removeprefix("p.").strip("()")

    def rep(m):
        tok = m.group(0)
        return AA3_TO_1.get(tok.upper(), tok)

    return re.sub(r"[A-Z][a-z]{2}", rep, s)


def _strip_ref(h: str | None) -> str | None:
    return unquote(h.split(":", 1)[1]) if h and ":" in h else h


def _clinvar_label(sigs: list[str]) -> str | None:
    uniq = list(dict.fromkeys(s.strip() for s in sigs if s and s.strip()))
    if len(uniq) > 1:
        uniq = [s for s in uniq if s != "not_provided"]
    uniq.sort(key=lambda s: CLINVAR_ORDER.index(s) if s in CLINVAR_ORDER else len(CLINVAR_ORDER))
    return ", ".join(uniq) or None


def colocated_info(colocated: list[dict], alt_plus: str | None, alt_raw: str | None = None) -> dict:
    """ClinVar/gnomAD/dbSNP/COSMIC/PubMed for one allele.

    clin_sig_allele is keyed by the forward-strand allele, while VEP keys
    `frequencies` by the allele as given in the input (reverse strand for HGVS
    on minus-strand genes), so both spellings are needed.
    """
    out = {"rsid": None, "_clinvar_vep": None, "_vcvs": [], "gnomad_af": None, "gnomad_popmax_af": None,
           "cosmic": None, "pubmed_count": 0, "pubmed_ids": None}
    sigs, vcv, cosmic, pubmed, rs = [], [], [], [], []
    af, popmax = None, None
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
            pops = [v for k, v in fr.items() if k.startswith(("gnomade_", "gnomadg_")) and v is not None
                    and k.split("_", 1)[1] in acmg.GRPMAX_GROUPS]
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


def _num(x) -> float | None:
    """Float, or None for missing values (None, NaN, '', non-numeric)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _str(x) -> str | None:
    return x if isinstance(x, str) and x else None


MODERN = ("AlphaMissense", "REVEL", "EVE", "CADD")


def _modern_present(row: dict) -> bool:
    return bool(_str(row.get("alphamissense_class")) or _num(row.get("revel")) is not None
                or _num(row.get("cadd_phred")) is not None or _str(row.get("eve_class")))


def _votes(row: dict) -> tuple[list[str], list[str], list[str]]:
    """Damaging/benign votes per predictor. Accepts pipeline dicts and table rows (NaN = missing)."""
    dmg, ben, detail = [], [], []
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


def classify(row: dict, terms: set[str]) -> tuple[str, float | None, str]:
    """Return (prediction label, effect score 0-1, human-readable basis)."""
    sp = row.get("spliceai_max")
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
    if terms & INFRAME_TERMS:
        cd = row.get("cadd_phred")
        if cd is not None and cd >= CADD_DAMAGING:
            return "Protein-altering in-frame (likely damaging)", 0.75, f"CADD {cd:.1f}"
        return "Protein-altering in-frame", 0.5 if cd is None else 0.25, f"CADD {cd:.1f}" if cd is not None else "no scores available"
    if sp is not None and sp >= 0.5:
        return "Likely splice-altering", round(sp, 2), f"SpliceAI {sp:.2f}"
    if sp is not None and sp >= 0.2:
        return "Possible splice effect", round(sp, 2), f"SpliceAI {sp:.2f}"
    if terms & {"synonymous_variant", "stop_retained_variant", "start_retained_variant"}:
        return "Silent (no amino-acid change)", 0.0, f"SpliceAI {sp:.2f}" if sp is not None else ""
    if terms & SPLICE_NEAR_TERMS:
        return "Splice region (no strong splice signal)", 0.1 if sp is not None else None, f"SpliceAI {sp:.2f}" if sp is not None else ""
    return "Non-coding / no protein change", None, ""


def is_damaging(label) -> bool:
    label = label or ""
    return (label.startswith(("Likely damaging", "Loss of function", "Likely splice-altering"))
            or label == "Protein-altering in-frame (likely damaging)")


STATUSES = ["Damaging", "Uncertain", "Benign / silent", "Not assessed"]


def prediction_status(label) -> str:
    """Collapse a prediction label into one of STATUSES (for colour, filters and charts)."""
    label = label if isinstance(label, str) else ""
    if is_damaging(label):
        return "Damaging"
    if "benign" in label.lower() or label.startswith("Silent"):
        return "Benign / silent"
    if label.startswith(("Uncertain", "Possible", "Protein-altering", "Splice region")):
        return "Uncertain"
    return "Not assessed"


def apply_frequency(label: str, basis: str, grpmax, common_af: float | None) -> tuple[str, str]:
    """BA1-style rule: a damaging or uncertain call becomes benign when the variant is common in gnomAD.

    Kept separate from the score so it stays visible: the basis still shows the protein-level call."""
    af = _num(grpmax)
    if not common_af or af is None or af < common_af or prediction_status(label) not in ("Damaging", "Uncertain"):
        return label, basis
    note = (f"common in gnomAD: grpmax {af * 100:.2g}% >= {common_af * 100:g}% (protein-level call: {label}); "
            "frequent variants can still cause recessive or founder disease")
    return "Likely benign (common in population)", f"{note} | {basis}"


def _nmd(exon: str | None, terms: set[str]) -> str | None:
    if not exon or not terms & {"stop_gained", "frameshift_variant"}:
        return None
    try:
        a, b = (int(x) for x in exon.split("/")[0:2])
    except ValueError:
        return None
    if b == 1:
        return "NMD escape likely (single-exon gene)"
    if a == b:
        return "NMD escape likely (last exon)"
    if a == b - 1:
        return "NMD uncertain (penultimate exon; depends on distance to last junction)"
    return "NMD likely (premature stop upstream of last exon)"


# --------------------------------------------------------------------------- rows

def _rows_for_variant(v: VariantInput, results: list[dict], assembly: str, genes: set[str] | None,
                      include_noncoding: bool) -> list[dict]:
    rows = []
    for r in results:
        strand = r.get("strand", 1)
        alleles = (r.get("allele_string") or "").split("/")
        g_ref = alleles[0] if alleles else None
        if strand == -1 and g_ref:
            g_ref = revcomp(g_ref)
        tcs = r.get("transcript_consequences") or []
        coloc = r.get("colocated_variants") or []
        base = {
            "input": v.raw, "source": v.source, "sample": v.sample, "zygosity": v.zygosity,
            "genotype": v.genotype, "carriers": v.carriers, "vcf_id": v.vcf_id, "vcf_filter": v.vcf_filter,
            "assembly": assembly, "variant_class": r.get("variant_class"),
            "chrom": v.chrom or r.get("seq_region_name"), "pos": v.pos or r.get("start"),
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
                row["_spdi"] = (str(r.get("seq_region_name")), int(r["start"]) - 1,
                                g_ref.replace("-", ""), alt_plus.replace("-", ""))
            aa = tc.get("amino_acids") or ""
            ref_aa, _, alt_aa = aa.partition("/")
            hgvsp = tc.get("hgvsp")
            # REST nests plugin output; a local VEP may report flat keys (am_class, spliceai_pred_ds_ag, ...)
            am = tc.get("alphamissense") or {"am_class": tc.get("am_class"), "am_pathogenicity": tc.get("am_pathogenicity")}
            amc = am.get("am_class")
            amc = {"pathogenic": "likely_pathogenic", "benign": "likely_benign"}.get(amc, amc)
            sp = tc.get("spliceai") or {k: tc.get(f"spliceai_pred_{k.lower()}") for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")}
            sp_vals = [_f(sp.get(k)) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")]
            sp_vals = [x for x in sp_vals if x is not None]
            sp_max = max(sp_vals) if sp_vals else None
            acc = (tc.get("swissprot") or tc.get("trembl") or [None])[0]
            row.update({
                "gene": tc.get("gene_symbol") or tc.get("gene_id"),
                "consequence": ",".join(tc.get("consequence_terms") or []),
                "impact": tc.get("impact"),
                "hgvsc": _strip_ref(tc.get("hgvsc")),
                "hgvsp": _strip_ref(hgvsp),
                "protein_change": short_protein_change(hgvsp) or (
                    f"{ref_aa}{tc.get('protein_start')}{alt_aa}" if ref_aa and alt_aa else None),
                "transcript": tc.get("transcript_id"),
                "mane_select": tc.get("mane_select"),
                "canonical": bool(tc.get("canonical")),
                "protein_id": tc.get("protein_id"),
                "uniprot": acc.split(".")[0] if acc else None,
                "exon": tc.get("exon"), "intron": tc.get("intron"),
                "codons": tc.get("codons"),
                "cds_position": tc.get("cds_start"),
                "protein_position": tc.get("protein_start"),
                "_protein_end": tc.get("protein_end"),
                "ref_aa": ref_aa or None, "alt_aa": alt_aa or None,
                "sift": tc.get("sift_prediction"), "sift_score": _f(tc.get("sift_score")),
                "polyphen": tc.get("polyphen_prediction"), "polyphen_score": _f(tc.get("polyphen_score")),
                "alphamissense": _f(am.get("am_pathogenicity")), "alphamissense_class": amc,
                "revel": _f(tc.get("revel")) if tc.get("revel") is not None else _max_listed(tc.get("revel_score")),
                "cadd_phred": _f(tc.get("cadd_phred")),
                "eve": _f(tc.get("eve_score")), "eve_class": tc.get("eve_class"),
                "popeve": _f(tc.get("popeve_score")),
                "spliceai_max": sp_max,
                "_terms": terms,
            })
            row["nmd"] = _nmd(row["exon"], terms)
            if "missense_variant" in terms and len(ref_aa) == 1 and len(alt_aa) == 1:
                row.update(substitution_profile(ref_aa, alt_aa))
            rows.append(row)
    return rows


# --------------------------------------------------------------------------- genome build

def _ref_bases(http: Http, assembly: str, variants: list[VariantInput], cache=None) -> dict[str, str]:
    """query -> genome sequence at the REF position(s)."""
    cache = cache or NullCache()
    vs = [v for v in variants if v.chrom in CHROMS and v.ref and v.ref != "-" and len(v.ref) <= 50]
    key = {v.query: f"refseq|{assembly}|{v.chrom}:{v.pos}:{len(v.ref)}" for v in vs}
    hit = cache.get_many(list(set(key.values())))
    out = {q: hit[k] for q, k in key.items() if k in hit}
    vs = [v for v in vs if v.query not in out]
    for i in range(0, len(vs), 50):
        chunk = vs[i : i + 50]
        regs = {f"{v.chrom}:{v.pos}..{v.pos + len(v.ref) - 1}:1": v for v in chunk}
        try:
            data = http.post_json(f"{ENSEMBL[assembly]}/sequence/region/human", {"regions": list(regs)},
                                  headers={"Content-Type": "application/json"})
        except ApiError as e:
            log.warning("reference check failed: %s", e)
            continue
        for item in data:
            v = regs.get(item.get("query"))
            if v is not None:
                out[v.query] = (item.get("seq") or "").upper()
                cache.set(key[v.query], out[v.query])
    return out


def _match_rate(http, assembly, variants, cache=None) -> tuple[float | None, dict]:
    seqs = _ref_bases(http, assembly, variants, cache)
    checked = [v for v in variants if v.query in seqs]
    if not checked:
        return None, seqs
    ok = sum(seqs[v.query] == v.ref for v in checked)
    return ok / len(checked), seqs


def gene_regions(http: Http, cache, genes: list[str], assembly: str) -> tuple[list, list[str]]:
    regions, missing = [], []
    todo = []
    for g in genes:
        hit = cache.get(f"generegion|{assembly}|{g}")
        if hit:
            regions.append(tuple(hit))
        else:
            todo.append(g)
    if todo:
        data = http.post_json(f"{ENSEMBL[assembly]}/lookup/symbol/homo_sapiens", {"symbols": todo},
                              headers={"Content-Type": "application/json"})
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
    from .inputs import _assembly_from_header
    import gzip
    import io

    try:
        raw = data if data is not None else open(name, "rb").read(2_000_000)
        if raw[:2] == b"\x1f\x8b":
            with gzip.open(io.BytesIO(raw)) as fh:
                raw = fh.read(2_000_000)
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

def _file_info(name: str, data: bytes | None) -> dict:
    """Name, size and SHA-256 of an input file (for the run record)."""
    import hashlib
    import os

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


def _vcf_alleles(r: dict) -> bool:
    """SNVs and VCF-style indels (no '-'), which gnomAD can look up; not mitochondrial."""
    ref, alt = str(r.get("ref") or ""), str(r.get("alt") or "")
    return (bool(ref) and bool(alt) and "-" not in ref + alt and ref.isalpha() and alt.isalpha()
            and str(r.get("chrom") or "").upper() not in ("MT", "M", "CHRM"))


def _source_warnings(status: SourceStatus) -> list[str]:
    out = []
    for src, msgs in status.errors.items():
        n = len(status.failed_items.get(src, ()))
        out.append(f"{src} did not respond{f' for {n} item(s)' if n else ''} ({msgs[0][:160]}). Affected rows are "
                   "marked in 'notes'. Run the analysis again to retry: everything already obtained is cached.")
    return out


def predict(text: str = "", files=(), options: PredictOptions | None = None, progress=None) -> PredictionResult:
    """Annotate variants given as free text and/or files [(name, bytes_or_None), ...]."""
    import platform
    from dataclasses import asdict
    from datetime import datetime, timezone

    from . import __version__

    t0 = time.time()
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    opt = options or PredictOptions()
    cache = Cache(opt.cache_path) if opt.cache else NullCache()
    http = Http()
    status = SourceStatus()
    warnings: list[str] = []
    say = progress or (lambda *a, **k: None)
    local = opt.vep == "local"
    net = not opt.offline
    if opt.offline and not local:
        raise ValueError("offline mode needs the local VEP (vep='local')")
    try:
        genes = [g.strip().upper() for g in (opt.genes or []) if g and g.strip()]
        files = [(f, None) if isinstance(f, str) else f for f in files]

        # ---- genome build (first guess: explicit > VCF header > GRCh38)
        say("Reading input", 0, 1)
        hint = None
        for name, data in files:
            hint = hint or _peek_assembly(name, data)
        assembly = opt.assembly if opt.assembly != "auto" else (hint or "GRCh38")

        # ---- connectivity: one quick request tells us whether Ensembl is reachable at all
        if net:
            try:
                ensembl_release(http, cache, "GRCh38")
            except OfflineError as e:
                if not local:
                    raise ConnectionError(
                        "Cannot reach the Ensembl REST API (rest.ensembl.org). Check the internet connection, "
                        "or use the local VEP mode for VCF/genomic input.") from e
                net = False
                warnings.append("No internet connection: using local VEP only; ClinVar, UniProt and AlphaFold "
                                "are skipped.")
            except ApiError as e:
                log.warning("Ensembl release lookup failed: %s", e)

        def load(asm):
            regions = None
            if genes and files:
                if net:
                    try:
                        regions, missing = gene_regions(http, cache, genes, asm)
                        if missing:
                            warnings.append(f"Gene(s) not found in Ensembl: {', '.join(missing)}")
                    except ApiError as e:
                        status.fail("Ensembl gene lookup", str(e), genes)
                        regions = None
                if regions is None:
                    warnings.append("Gene coordinates unavailable: the gene filter is applied after annotation "
                                    "(slower for large files).")
            parsed = ParsedInput()
            if text and text.strip():
                parsed.extend(parse_text(text, gene=opt.default_gene))
            for name, data in files:
                parsed.extend(read_file(name, data, gene=opt.default_gene, sample=opt.sample,
                                        regions=regions, pass_only=opt.pass_only))
            return parsed

        parsed = load(assembly)
        say("Reading input", 1, 1)

        # ---- verify REF alleles against the genome; auto-switch build if clearly wrong
        region_vars = [v for v in parsed.variants if v.ok and v.kind == "region"]
        ref_seqs = {}
        if region_vars and net:
            say("Checking genome build", 0, 1)
            small = len(region_vars) <= 200
            sample = region_vars if small else region_vars[:: max(1, len(region_vars) // 50)][:50]
            try:
                rate, ref_seqs = _match_rate(http, assembly, sample, cache)
                if rate is not None and rate < 0.9:
                    other = "GRCh37" if assembly == "GRCh38" else "GRCh38"
                    rate2, seqs2 = _match_rate(http, other, sample, cache)
                    if opt.assembly == "auto" and not hint and rate2 is not None and rate2 > rate and rate2 >= 0.9:
                        warnings.append(f"REF alleles match {other} ({rate2:.0%}) better than {assembly} ({rate:.0%}); "
                                        f"using {other}.")
                        assembly, ref_seqs = other, seqs2
                        if genes and files:
                            parsed = load(assembly)
                            region_vars = [v for v in parsed.variants if v.ok and v.kind == "region"]
                    else:
                        warnings.append(f"Only {rate:.0%} of checked REF alleles match {assembly}"
                                        + (f" ({other}: {rate2:.0%})" if rate2 is not None else "")
                                        + ". Check the genome build and coordinates.")
            except ApiError as e:
                status.fail("Genome build check", str(e))
            say("Checking genome build", 1, 1)
        elif region_vars and opt.assembly == "auto" and not hint:
            warnings.append(f"Genome build not verified (offline); assuming {assembly}.")
        ref_note = {}
        for v in region_vars:
            g = ref_seqs.get(v.query)
            if g and g != v.ref:
                ref_note[v.query] = f"REF {v.ref} does not match {assembly} ({g})"

        release = None
        if net:
            try:
                release = ensembl_release(http, cache, assembly)
            except ApiError as e:
                log.warning("Ensembl release lookup failed: %s", e)

        # ---- VEP (local and/or REST)
        ok_vars = [v for v in parsed.variants if v.ok]
        queries: dict[str, list[str]] = {}
        for v in ok_vars:
            queries.setdefault(v.kind, []).append(v.query)
        n_unique = sum(len(set(q)) for q in queries.values())
        full = {**BASE_PARAMS, **PREDICTOR_PARAMS}
        no_predictors: set = set()   # (kind, query) annotated, but the predictor plugins failed
        res: dict = {}
        local_version = None
        if local:
            lv = LocalVep(cache, assembly, opt.vep_command, opt.vep_args, status=status)
            local_version = lv.version()
            res.update(lv.annotate({"region": queries.get("region", [])}, None, say, "Annotating variants (local VEP)"))
            rest_queries = {k: v for k, v in queries.items() if k != "region"}
            if rest_queries and not net:
                res.update(lv.annotate(rest_queries))
                rest_queries = {}
        else:
            rest_queries = queries
        if rest_queries:
            vep = Vep(http, cache, assembly, workers=opt.workers, release=release, status=status)
            n_rest = sum(len(set(q)) for q in rest_queries.values())
            if not opt.predictors:
                res.update(vep.annotate(rest_queries, BASE_PARAMS, say, "Annotating variants (VEP)"))
            elif n_rest <= 200:
                part = vep.annotate(rest_queries, full, say, "Annotating + predictors (VEP)",
                                    source="Ensembl VEP + predictors")
                retry: dict[str, list[str]] = {}
                for (kind, q), val in part.items():
                    if val.get("service_error"):
                        retry.setdefault(kind, []).append(q)
                if retry:  # plugins are the usual culprit: try again without them
                    base = vep.annotate(retry, BASE_PARAMS, say, "Retrying without predictors")
                    for k, val in base.items():
                        if "results" in val:
                            part[k] = val
                            no_predictors.add(k)
                res.update(part)
            else:
                part = vep.annotate(rest_queries, BASE_PARAMS, say, "Annotating variants (VEP, fast pass)")
                detail: dict[str, list[str]] = {}
                for (kind, q), val in part.items():
                    for r in val.get("results", []):
                        if any((set(tc.get("consequence_terms") or []) & DETAIL_CONSEQUENCES)
                               and (not genes or (tc.get("gene_symbol") or "").upper() in genes)
                               for tc in r.get("transcript_consequences") or []):
                            detail.setdefault(kind, []).append(q)
                            break
                if detail:
                    scored = vep.annotate(detail, full, say, "Scoring coding variants (predictors)",
                                          source="VEP predictor plugins")
                    for k, val in scored.items():
                        if "results" in val:
                            part[k] = val
                        else:  # keep the base annotation instead of losing the variant
                            no_predictors.add(k)
                res.update(part)

        # ---- rows
        rows, errors = [], []
        n_no_protein = 0
        gene_set = set(genes) if genes else None
        for v in parsed.variants:
            if not v.ok:
                errors.append({"input": v.raw, "source": v.source, "error": v.error or "could not parse"})
                continue
            val = res.get((v.kind, v.query), {})
            if "error" in val:
                errors.append({"input": v.raw, "source": v.source, "error": val["error"]})
                continue
            vrows = _rows_for_variant(v, val.get("results", []), assembly, gene_set, opt.include_noncoding)
            if not vrows:
                n_no_protein += 1
            for r in vrows:
                notes = []
                if v.query in ref_note:
                    notes.append(ref_note[v.query])
                if (v.kind, v.query) in no_predictors:
                    notes.append("predictors unavailable (service error; run again to retry)")
                r["notes"] = "; ".join(notes) or None
            rows.extend(vrows)

        # ---- ClinVar (exact allele)
        records = {}
        if opt.clinvar and rows and net:
            say("ClinVar classification", 0, 1)
            records = cv.fetch(http, cache, {x for r in rows for x in r.get("_vcvs", [])}, status=status)
            say("ClinVar classification", 1, 1)
        elif opt.clinvar and rows:
            status.skip("ClinVar", "offline: classification taken from the VEP data (all submissions)")
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
            if any(status.failed("ClinVar", x) for x in r.get("_vcvs", [])):
                r["notes"] = "; ".join(x for x in (r.get("notes"), "ClinVar unavailable (service error)") if x)
            for k in ("_vcvs", "_spdi", "_clinvar_vep"):
                r.pop(k, None)

        # ---- ACMG PS1/PM5: pathogenic missense variants at the same residue (ClinVar)
        windows = {i: w for i, r in enumerate(rows) if (w := residue.window(r))}
        residue_refs: dict[int, list] = {}
        if windows and opt.clinvar and net:
            say("ClinVar: same residue", 0, 1)
            found = cv.search_windows(http, cache, list(windows.values()), assembly, status=status)
            recs = cv.fetch(http, cache, {v for w in found.values() for v in w}, status=status)
            for i, w in windows.items():
                if w in found and all(v in recs for v in found[w]):
                    residue_refs[i] = residue.references(rows[i], {v: recs[v] for v in found[w]}, assembly)
            say("ClinVar: same residue", 1, 1)

        # ---- gnomAD itself, for variants without a frequency in Ensembl's data (singletons, PM2)
        need = {i: (str(r["chrom"]), int(r["pos"]), str(r["ref"]), str(r["alt"])) for i, r in enumerate(rows)
                if r.get("gnomad_af") is None and r.get("pos") and _vcf_alleles(r)}
        gnomad_info: dict[int, dict] = {}
        if need and opt.gnomad_check and net:
            say("gnomAD", 0, 1)
            got = gnomad.lookup(http, cache, list(need.values()), assembly, status=status,
                                progress=lambda k, n: say("gnomAD", k, n))
            for i, k in need.items():
                if k in got:
                    gnomad_info[i] = got[k]
                    if got[k]["found"]:
                        rows[i]["gnomad_af"], rows[i]["gnomad_popmax_af"] = got[k]["af"], got[k]["grpmax"]
            say("gnomAD", 1, 1)
        elif need and opt.gnomad_check:
            status.skip("gnomAD", "offline")

        # ---- protein context
        uni, afd, seqs = {}, {}, {}
        prot_rows = [r for r in rows if r.get("protein_position")]
        if prot_rows and opt.protein_context and net:
            say("Protein context (UniProt, AlphaFold)", 0, 1)
            accs = {r["uniprot"] for r in prot_rows if r.get("uniprot")}
            uni = prot.fetch_uniprot(http, cache, accs, workers=opt.workers, status=status)
            afd = prot.fetch_alphafold(http, cache, accs, workers=opt.workers, status=status)
            say("Protein context (UniProt, AlphaFold)", 1, 1)
        elif prot_rows and opt.protein_context:
            status.skip("UniProt", "offline")
            status.skip("AlphaFold DB", "offline")
        if prot_rows and opt.sequences and net:
            seqs = prot.fetch_protein_sequences(http, cache, {r.get("protein_id") for r in prot_rows}, assembly,
                                                status=status)

        sequences = {}
        for i, r in enumerate(rows):
            terms = r.pop("_terms", set()) or set()
            pend = r.pop("_protein_end", None)
            notes = [r["notes"]] if r.get("notes") else []
            pos = r.get("protein_position")
            pid = (r.get("protein_id") or "").split(".")[0]
            wt = seqs.get(pid)
            if pid and not wt and status.failed("Ensembl protein sequences", pid):
                notes.append("protein sequence unavailable (service error)")
            if wt:
                r["protein_length"] = len(wt)
            if pos and wt and r.get("ref_aa"):
                r["protein_window"] = prot.context_window(wt, pos, r["ref_aa"], r.get("alt_aa"))
                mut, note = prot.mutate_protein(wt, pos, pend or pos, f"{r['ref_aa']}/{r.get('alt_aa') or '-'}",
                                                ",".join(terms))
                r["mutant_protein"] = note
                if note.startswith("truncated"):
                    r["_truncation"] = note
                if mut is not None:
                    sequences[i] = {"wt": wt, "mut": mut, "note": note,
                                    "header": f"{r.get('gene')}|{r.get('protein_change')}|{r.get('protein_id')}|{r.get('input')}"}
            acc = r.get("uniprot") or ""
            if acc and pos and status.failed("UniProt", acc):
                notes.append("UniProt unavailable (service error)")
            if acc and pos and status.failed("AlphaFold DB", acc):
                notes.append("AlphaFold unavailable (service error)")
            entry = uni.get(acc)
            if entry and pos:
                r["protein_name"] = entry.get("name")
                r["protein_function"] = entry.get("function")
                r["protein_length"] = r.get("protein_length") or entry.get("length")
                useq = entry.get("sequence") or ""
                same = bool(wt) and useq == wt
                local_ok = r.get("ref_aa") and useq[pos - 1 : pos - 1 + len(r["ref_aa"])] == r["ref_aa"] and (
                    not wt or useq[max(0, pos - 6) : pos + 5] == wt[max(0, pos - 6) : pos + 5])
                if same or local_ok:
                    r.update(prot.features_at(entry, pos, pend or pos))
                    af = afd.get(r["uniprot"])
                    r["plddt"] = prot.plddt_at(af, pos)
                    r["plddt_class"] = prot.plddt_class(r["plddt"])
                    if af and af.get("available"):
                        r["alphafold_model"] = af.get("entry_id")
                    if not same:
                        notes.append("UniProt canonical isoform differs from transcript; features mapped by local sequence match")
                else:
                    notes.append("UniProt canonical numbering differs from this transcript; UniProt/AlphaFold features skipped")
            label, score, basis = classify(r, terms)
            r.pop("_truncation", None)
            if "missense_variant" in terms:
                r["acmg_pp3_bp4"], r["acmg_points"] = acmg.computational(r, opt.acmg_source)
            if i in windows:
                r["acmg_ps1_pm5"] = residue.evidence(r, residue_refs[i])[0] if i in residue_refs else (
                    "not checked (" + ("ClinVar turned off" if not opt.clinvar else "offline" if not net else
                                       "ClinVar did not respond") + ")")
                if status.failed("ClinVar", f"{assembly}:{windows[i][0]}:{windows[i][1]}-{windows[i][2]}"):
                    notes.append("ClinVar same-residue search unavailable (service error)")
            r["_gnomad"] = gnomad_info.get(i)
            r["acmg_pm2"] = acmg.pm2(r, r["_gnomad"])[0]
            if i in need and i not in gnomad_info and opt.gnomad_check and net:
                notes.append("gnomAD check unavailable (service error)")
            r["population_frequency"] = acmg.frequency(r)
            r.pop("_gnomad")
            r["acmg_codes"] = acmg.codes(r)
            if not acmg.ba1_exception(r):     # common but pathogenic/uncertain (ClinGen list): no frequency rule
                label, basis = apply_frequency(label, basis, r.get("gnomad_popmax_af"), opt.common_af)
            r["prediction"], r["effect_score"], r["prediction_basis"] = label, score, basis
            r["spliceai_effect"] = ("likely splice-altering" if (r.get("spliceai_max") or 0) >= 0.5 else
                                    "possible" if (r.get("spliceai_max") or 0) >= 0.2 else None)
            r["notes"] = "; ".join(dict.fromkeys(notes)) or None

        table = pd.DataFrame(rows, columns=COLUMNS) if rows else pd.DataFrame(columns=COLUMNS)
        errs = pd.DataFrame(errors, columns=["input", "source", "error"])
        stats = {
            "inputs": len(parsed.variants),
            "unique_queries": n_unique,
            "rows": len(table),
            "protein_altering": int(table["consequence"].fillna("").str.contains(
                "|".join(PROTEIN_ALTERING)).sum()) if len(table) else 0,
            "likely_damaging": int(table["prediction"].map(is_damaging).sum()) if len(table) else 0,
            "no_protein_effect": n_no_protein,
            "errors": len(errs),
            "seconds": round(time.time() - t0, 1),
            **{f"file_{k}": v for k, v in parsed.stats.items() if k.startswith(("vcf_", "skipped_"))},
        }
        if parsed.samples and len(parsed.samples) > 1 and not opt.sample:
            warnings.append(f"Multi-sample VCF ({len(parsed.samples)} samples): carriers listed per variant. "
                            "Choose one sample to get its zygosity.")
        warnings.extend(parsed.warnings)
        warnings.extend(_source_warnings(status))
        run = {
            "tool": f"MissensExplorer {__version__}",
            "started_utc": started,
            "seconds": stats["seconds"],
            "assembly": assembly,
            "annotation": (f"local VEP ({local_version or 'version unknown'})" if local else
                           f"Ensembl VEP REST, release {release or 'unknown'}"),
            "ensembl_release": release,
            "options": {k: v for k, v in asdict(opt).items() if k not in ("cache_path",)},
            "inputs": {"text_lines": len([x for x in (text or "").splitlines() if x.strip()]),
                       "files": [_file_info(n, d) for n, d in files]},
            "sources": status.summary(),
            "missense_model": {"version": model.MODEL["version"], "trained": model.MODEL["trained"],
                               "training_variants": model.MODEL["training"]["n"],
                               "cut_offs": [model.BENIGN, model.DAMAGING]},
            "python": platform.python_version(),
        }
        say("Done", 1, 1)
        return PredictionResult(table=table, errors=errs, assembly=assembly, warnings=warnings, stats=stats,
                                sequences=sequences, alphafold=afd, uniprot=uni, run=run, sources=status.summary())
    finally:
        http.close()


def predict_text(text: str, **kwargs) -> pd.DataFrame:
    """Convenience: predict_text('BRAF V600E') -> DataFrame."""
    return predict(text=text, options=PredictOptions(**kwargs)).table
