"""Authoritative ClinVar classification for an exact allele (NCBI E-utilities, batched)."""

from __future__ import annotations

import logging
import re
import time

from .api import ApiError, Http

log = logging.getLogger(__name__)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

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
    spdi, cdna = None, None
    if len(vsets) == 1:
        spdi = vsets[0].get("canonical_spdi")
        cdna = vsets[0].get("cdna_change")
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
    }


CLINVAR_MAX_AGE = 7 * 86400  # ClinVar is updated weekly


def fetch(http: Http, cache, vcv_ids, status=None) -> dict[str, dict]:
    """{'VCV000013961': {...}, ...}. Failed requests are recorded in `status` (source 'ClinVar')."""
    ids = sorted({v for v in vcv_ids if v and v.startswith("VCV")})
    if status is not None and ids:
        status.use("ClinVar", len(ids))
    keys = {v: f"clinvar|{v}" for v in ids}
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
