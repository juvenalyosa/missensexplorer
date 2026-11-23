"""ClinVar through NCBI E-utilities (batched): the classification of an exact allele, and the
pathogenic variants near a position (for ACMG PS1/PM5, see residue.py)."""

from __future__ import annotations

import logging
import re
import time

from .api import ApiError, Http

log = logging.getLogger(__name__)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PATHOGENIC_FILTER = "(clinsig_pathogenic[prop] OR clinsig_likely_pathogenic[prop])"
WINDOWS_PER_SEARCH = 25

STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
}

_NC_CHROM = {f"NC_{i:06d}": str(i) for i in range(1, 23)}
_NC_CHROM.update({"NC_000023": "X", "NC_000024": "Y", "NC_012920": "MT"})


def _compact(r: dict) -> dict:
    g = r.get("germline_classification") or {}
    onc = (r.get("oncogenicity_classification") or {}).get("description") or None
    som = (r.get("clinical_impact_classification") or {}).get("description") or None
    vsets = r.get("variation_set") or []
    spdi, cdna, loc = None, None, {}
    if len(vsets) == 1:
        spdi = vsets[0].get("canonical_spdi")
        cdna = vsets[0].get("cdna_change")
        for x in vsets[0].get("variation_loc") or []:
            if x.get("assembly_name") in ("GRCh37", "GRCh38") and str(x.get("start") or "").isdigit():
                loc[x["assembly_name"]] = [x.get("chr"), int(x["start"]), int(x.get("stop") or x["start"])]
    return {
        "accession": r.get("accession"),
        "title": r.get("title"),
        "classification": g.get("description") or None,
        "review_status": g.get("review_status") or None,
        "last_evaluated": (g.get("last_evaluated") or "").split(" ")[0] or None,
        "conditions": [t.get("trait_name") for t in g.get("trait_set", []) if t.get("trait_name")],
        "oncogenicity": onc,
        "somatic_impact": som,
        "spdi": spdi,
        "cdna": cdna,
        "loc": loc,                                    # assembly -> [chrom, start, stop] (1-based)
    }


CLINVAR_MAX_AGE = 7 * 86400  # ClinVar is updated weekly


def fetch(http: Http, cache, vcv_ids, status=None) -> dict[str, dict]:
    """{'VCV000013961': {...}, ...}. Failed requests are recorded in `status` (source 'ClinVar')."""
    ids = sorted({v for v in vcv_ids if v and v.startswith("VCV")})
    if status is not None and ids:
        status.use("ClinVar", len(ids))
    keys = {v: f"clinvar2|{v}" for v in ids}    # clinvar2: records with positions
    hit = cache.get_many(list(keys.values()), max_age=CLINVAR_MAX_AGE)
    out = {v: hit[keys[v]] for v in ids if keys[v] in hit}
    missing = [v for v in ids if v not in out]
    for i in range(0, len(missing), 200):
        chunk = missing[i : i + 200]
        uid = {str(int(v[3:])): v for v in chunk}
        try:
            r = http.request("POST", EUTILS, data={"db": "clinvar", "id": ",".join(uid), "retmode": "json",
                                                   "tool": "missensexplorer"})
            if r.status_code != 200:
                raise ApiError(f"HTTP {r.status_code}", r.status_code)
            data = r.json().get("result", {})
        except (ApiError, ValueError) as e:
            log.warning("ClinVar request failed: %s", e)
            if status is not None:
                status.fail("ClinVar", str(e), chunk)
            continue
        res = {}
        for u in data.get("uids", []):
            if u in uid and isinstance(data.get(u), dict):
                res[uid[u]] = _compact(data[u])
        out.update(res)
        cache.set_many({keys[v]: val for v, val in res.items()})
        if i + 200 < len(missing):
            time.sleep(0.4)  # NCBI etiquette without an API key: <= 3 requests/s
    return out


def _parse_spdi(spdi: str | None):
    if not spdi:
        return None
    try:
        seq, pos, dele, ins = spdi.split(":")
        return _NC_CHROM.get(seq.split(".")[0]), int(pos), dele, ins
    except ValueError:
        return None


def match(records: dict[str, dict], candidates: list[str], spdi_key: tuple | None,
          mane: str | None, hgvsc: str | None) -> dict | None:
    """Pick the ClinVar record describing exactly this allele."""
    for v in candidates:
        rec = records.get(v)
        if not rec:
            continue
        p = _parse_spdi(rec.get("spdi"))
        if p and spdi_key and p == spdi_key:
            return rec
    if mane and hgvsc:
        nm = mane.split(".")[0]
        for v in candidates:
            rec = records.get(v)
            if not rec or not rec.get("title"):
                continue
            m = re.match(r"^(N[MR]_\d+)(?:\.\d+)?(?:\([^)]*\))?:(c\.[^\s]+)", rec["title"])
            if m and m.group(1) == nm and m.group(2) == hgvsc:
                return rec
    return None


def stars(review_status: str | None) -> int | None:
    return None if not review_status else STARS.get(review_status, 0)


def vcv(uid: str) -> str:
    return f"VCV{int(uid):09d}"


def search_windows(http: Http, cache, windows, assembly: str, status=None) -> dict[tuple, list[str]]:
    """VCV accessions of the pathogenic or likely pathogenic ClinVar variants that overlap each genomic
    window (chrom, start, end; 1-based, on `assembly`), {window: [VCV...]}. Searches are batched
    (WINDOWS_PER_SEARCH windows each) and the result of every window is cached for a week. Windows of
    a failed search are left out of the result and recorded in `status` (source 'ClinVar')."""
    field = "chrpos38" if assembly == "GRCh38" else "chrpos37"
    windows = sorted(set(windows))
    keys = {w: f"clinvar-window|{assembly}|{w[0]}:{w[1]}-{w[2]}" for w in windows}
    hit = cache.get_many(list(keys.values()), max_age=CLINVAR_MAX_AGE)
    out = {w: hit[keys[w]] for w in windows if keys[w] in hit}
    todo = [w for w in windows if w not in out]
    for i in range(0, len(todo), WINDOWS_PER_SEARCH):
        chunk = todo[i : i + WINDOWS_PER_SEARCH]
        term = " OR ".join(f"({c}[chr] AND {a}:{b}[{field}])" for c, a, b in chunk)
        try:
            r = http.request("POST", ESEARCH, data={"db": "clinvar", "term": f"({term}) AND {PATHOGENIC_FILTER}",
                                                    "retmode": "json", "retmax": 10000, "tool": "missensexplorer"})
            if r.status_code != 200:
                raise ApiError(f"HTTP {r.status_code}", r.status_code)
            ids = [vcv(u) for u in r.json()["esearchresult"].get("idlist", [])]
        except (ApiError, ValueError, KeyError) as e:
            log.warning("ClinVar search failed: %s", e)
            if status is not None:
                status.fail("ClinVar", f"search: {e}", [f"{assembly}:{c}:{a}-{b}" for c, a, b in chunk])
            continue
        records = fetch(http, cache, ids, status=status)
        if len(records) < len(ids):
            continue                                   # incomplete: do not cache these windows
        for w in chunk:
            out[w] = [v for v, rec in records.items() if _overlaps(rec, w, assembly)]
        cache.set_many({keys[w]: out[w] for w in chunk})
    return out


def _overlaps(rec: dict, window: tuple, assembly: str) -> bool:
    c, a, b = window
    loc = (rec.get("loc") or {}).get(assembly)
    if loc:
        return str(loc[0]) == str(c) and loc[1] <= b and loc[2] >= a
    p = _parse_spdi(rec.get("spdi")) if assembly == "GRCh38" else None
    return bool(p) and p[0] == str(c) and a <= p[1] + 1 <= b

