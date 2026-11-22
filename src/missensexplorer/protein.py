"""Protein-level context: UniProt features, AlphaFold confidence, Ensembl protein sequences."""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor

from .api import ALPHAFOLD, ENSEMBL, UNIPROT, ApiError, Http

log = logging.getLogger(__name__)

UNIPROT_FIELDS = ",".join([
    "accession", "protein_name", "gene_primary", "cc_function", "length", "sequence",
    "ft_domain", "ft_region", "ft_motif", "ft_act_site", "ft_binding", "ft_site", "ft_mod_res",
    "ft_lipid", "ft_carbohyd", "ft_disulfid", "ft_crosslnk", "ft_transmem", "ft_intramem",
    "ft_topo_dom", "ft_signal", "ft_transit", "ft_propep", "ft_zn_fing", "ft_dna_bind",
    "ft_repeat", "ft_coiled", "ft_compbias", "ft_variant", "ft_mutagen",
])

# Feature types reported in "protein_features" (order = display priority).
CONTEXT_TYPES = [
    "Active site", "Binding site", "Site", "Disulfide bond", "Modified residue", "Lipidation",
    "Glycosylation", "Cross-link", "Zinc finger", "DNA binding", "Motif", "Domain",
    "Transmembrane", "Intramembrane", "Signal", "Transit peptide", "Propeptide",
    "Topological domain", "Repeat", "Coiled coil", "Region", "Compositional bias",
]
# Feature types that mark residues with a direct functional role.
KEY_SITE_TYPES = {"Active site", "Binding site", "Site", "Disulfide bond", "Modified residue",
                  "Lipidation", "Glycosylation", "Cross-link", "Zinc finger", "DNA binding"}


# --------------------------------------------------------------------------- UniProt

def _compact_uniprot(entry: dict) -> dict:
    desc = entry.get("proteinDescription", {})
    name = (desc.get("recommendedName") or (desc.get("submissionNames") or [{}])[0]).get("fullName", {}).get("value")
    genes = entry.get("genes") or []
    gene = genes[0].get("geneName", {}).get("value") if genes else None
    func = None
    for c in entry.get("comments", []):
        if c.get("commentType") == "FUNCTION" and c.get("texts"):
            func = c["texts"][0].get("value")
            break
    feats = []
    for f in entry.get("features", []):
        loc = f.get("location", {})
        s, e = loc.get("start", {}).get("value"), loc.get("end", {}).get("value")
        if s is None or e is None:
            continue
        alt = f.get("alternativeSequence") or {}
        feats.append({
            "type": f.get("type"), "start": s, "end": e,
            "description": f.get("description") or "",
            "original": alt.get("originalSequence"),
            "alternatives": alt.get("alternativeSequences"),
        })
    seq = entry.get("sequence", {})
    return {"accession": entry.get("primaryAccession"), "name": name, "gene": gene,
            "function": func, "length": seq.get("length"), "sequence": seq.get("value"),
            "features": feats}


def fetch_uniprot(http: Http, cache, accessions, workers: int = 4) -> dict[str, dict]:
    accs = sorted({a for a in accessions if a})
    out = {}
    keys = {a: f"uniprot|{a}" for a in accs}
    hit = cache.get_many(list(keys.values()))
    missing = [a for a in accs if keys[a] not in hit]
    out.update({a: hit[keys[a]] for a in accs if keys[a] in hit})

    def run(chunk):
        try:
            data = http.get_json(f"{UNIPROT}/uniprotkb/accessions",
                                 params={"accessions": ",".join(chunk), "fields": UNIPROT_FIELDS,
                                         "format": "json", "size": 500})
        except ApiError as e:
            log.warning("UniProt request failed: %s", e)
            return {}
        return {r["primaryAccession"]: _compact_uniprot(r) for r in data.get("results", [])}

    chunks = [missing[i : i + 25] for i in range(0, len(missing), 25)]
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(run, chunks):
            out.update(res)
            cache.set_many({keys[a]: v for a, v in res.items() if a in keys})
    return out


def _feature_label(f: dict) -> str:
    span = f"{f['start']}" if f["start"] == f["end"] else f"{f['start']}-{f['end']}"
    desc = f" {f['description']}" if f["description"] else ""
    return f"{f['type']}:{desc} ({span})".replace(": (", " (")


def features_at(entry: dict, start: int, end: int | None = None) -> dict:
    """UniProt annotations overlapping residues start..end."""
    end = end or start
    ctx, sites, variants, mutagen = [], [], [], []
    for f in entry.get("features", []):
        t = f["type"]
        if t == "Disulfide bond":
            hit = start <= f["start"] <= end or start <= f["end"] <= end
        else:
            hit = f["start"] <= end and f["end"] >= start
        if not hit:
            continue
        if t == "Natural variant":
            alts = "/".join(f.get("alternatives") or ["?"])
            variants.append(f"{f.get('original') or ''}{f['start']}{alts}: {f['description']}")
        elif t == "Mutagenesis":
            alts = "/".join(f.get("alternatives") or ["?"])
            mutagen.append(f"{f.get('original') or ''}{f['start']}{alts}: {f['description']}")
        elif t in CONTEXT_TYPES:
            ctx.append((CONTEXT_TYPES.index(t), _feature_label(f)))
            # only residue-level features count as a "functional site" (not long regions)
            if t in KEY_SITE_TYPES and (t == "Disulfide bond" or f["end"] - f["start"] <= 10):
                sites.append(_feature_label(f))
    ctx.sort()
    return {
        "protein_features": "; ".join(lbl for _, lbl in ctx) or None,
        "key_site": "; ".join(sites) or None,
        "uniprot_variants_here": " | ".join(variants) or None,
        "uniprot_mutagenesis_here": " | ".join(mutagen) or None,
    }


# --------------------------------------------------------------------------- AlphaFold

def fetch_alphafold(http: Http, cache, accessions, workers: int = 4) -> dict[str, dict]:
    accs = sorted({a for a in accessions if a})
    out = {}

    def run(acc):
        key = f"alphafold|{acc}"
        val = cache.get(key)
        if val is not None:
            return acc, val
        val = {"available": False}
        try:
            r = http.request("GET", f"{ALPHAFOLD}/api/prediction/{acc}", ok_statuses=(200, 404, 400, 422))
            if r.status_code == 200 and r.json():
                entries = r.json()
                e = next((x for x in entries if x.get("uniprotAccession") == acc), entries[0])
                conf = http.get_json(e["plddtDocUrl"]) if e.get("plddtDocUrl") else {}
                val = {
                    "available": True,
                    "entry_id": e.get("modelEntityId") or e.get("entryId"),
                    "version": e.get("latestVersion"),
                    "pdb_url": e.get("pdbUrl"), "cif_url": e.get("cifUrl"),
                    "start": e.get("uniprotStart", 1), "end": e.get("uniprotEnd"),
                    "global_plddt": e.get("globalMetricValue"),
                    "plddt": conf.get("confidenceScore"),
                    "sequence": e.get("uniprotSequence") or e.get("sequence"),
                }
        except (ApiError, ValueError, KeyError) as ex:
            log.warning("AlphaFold lookup failed for %s: %s", acc, ex)
            return acc, None  # do not cache transient failures
        cache.set(key, val)
        return acc, val

    with ThreadPoolExecutor(max_workers=workers) as ex:
        for acc, val in ex.map(run, accs):
            if val is not None:
                out[acc] = val
    return out


def plddt_at(af: dict | None, pos: int | None) -> float | None:
    if not af or not af.get("available") or not pos or not af.get("plddt"):
        return None
    i = pos - (af.get("start") or 1)
    scores = af["plddt"]
    return round(float(scores[i]), 1) if 0 <= i < len(scores) else None


def plddt_class(v: float | None) -> str | None:
    if v is None:
        return None
    if v >= 90:
        return "very high"
    if v >= 70:
        return "confident"
    if v >= 50:
        return "low"
    return "very low (likely disordered)"


# --------------------------------------------------------------------------- sequences

def fetch_protein_sequences(http: Http, cache, protein_ids, assembly: str = "GRCh38") -> dict[str, str]:
    ids = sorted({p.split(".")[0] for p in protein_ids if p})
    keys = {p: f"protseq|{assembly}|{p}" for p in ids}
    hit = cache.get_many(list(keys.values()))
    out = {p: hit[keys[p]] for p in ids if keys[p] in hit}
    missing = [p for p in ids if p not in out]
    for i in range(0, len(missing), 50):
        chunk = missing[i : i + 50]
        try:
            data = http.post_json(f"{ENSEMBL[assembly]}/sequence/id", {"ids": chunk, "type": "protein"},
                                  headers={"Content-Type": "application/json"})
        except ApiError as e:
            log.warning("Ensembl protein sequence request failed: %s", e)
            continue
        res = {}
        for item in data:
            pid = (item.get("query") or item.get("id") or "").split(".")[0]
            if pid and item.get("seq"):
                res[pid] = item["seq"]
        out.update(res)
        cache.set_many({keys[p]: s for p, s in res.items()})
    return out


def mutate_protein(seq: str, start: int, end: int, amino_acids: str, consequence: str) -> tuple[str | None, str]:
    """Apply a VEP 'amino_acids' change (e.g. 'V/E', 'R/*', 'KL/-', 'S/SG') to a protein.

    Returns (mutant_sequence or None, note).
    """
    if "/" not in (amino_acids or ""):
        return (seq, "no amino-acid change") if "synonymous" in consequence else (None, "no amino-acid change reported")
    ref, alt = amino_acids.split("/", 1)
    ref = "" if ref == "-" else ref
    alt = "" if alt == "-" else alt
    if "frameshift" in consequence:
        return None, "frameshift: downstream sequence changes, not reconstructed"
    if "X" in alt:
        return None, "unknown residue in alternate allele"
    s0 = start - 1
    if ref and seq[s0 : s0 + len(ref)] != ref:
        return None, f"reference residue mismatch ({seq[s0:s0 + len(ref)] or '?'} in sequence vs {ref})"
    if not ref:  # pure insertion between start-1 and start (VEP: end = start - 1)
        s0 = max(end, start - 1) if end < start else start
        mut = seq[:s0] + alt + seq[s0:]
    else:
        mut = seq[:s0] + alt + seq[s0 + len(ref):]
    if "*" in mut:
        stop = mut.index("*")
        lost = len(seq) - stop
        return mut[:stop], f"truncated at residue {stop + 1}: {lost} of {len(seq)} residues lost ({100 * lost / len(seq):.0f}%)"
    return mut, "ok"


def context_window(seq: str, start: int, ref: str, alt: str, flank: int = 7) -> str | None:
    """e.g. 'GLATV[V>E]KSRWS' around the changed residue(s)."""
    if not seq or not start:
        return None
    s0 = start - 1
    ref = "" if ref in ("-", None) else ref
    left = seq[max(0, s0 - flank) : s0]
    right = seq[s0 + len(ref) : s0 + len(ref) + flank]
    return f"{left}[{ref or '-'}>{alt or '-'}]{right}"
