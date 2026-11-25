"""ClinVar through NCBI E-utilities (batched).

Provides the classification of an exact allele, the pathogenic variants near a position (for ACMG
PS1/PM5, see residue.py) and the variants with an oncogenicity classification near a position (for
the somatic OS1/OM4, see oncogenicity.py).
"""

from __future__ import annotations

import logging
import re
import time
from collections.abc import Iterable
from typing import Any

from ._types import JSON, KeyValueCache
from .api import ApiError, Http, SourceStatus

log = logging.getLogger(__name__)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PATHOGENIC_FILTER = "(clinsig_pathogenic[prop] OR clinsig_likely_pathogenic[prop])"
ONCOGENICITY_FILTER = "oncogenicity_classification[prop]"  # any oncogenicity classification
WINDOWS_PER_SEARCH = 25
CLINVAR_MAX_AGE = 7 * 86400  # ClinVar is updated weekly

# Review status -> the number of gold stars ClinVar shows for it.
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


def _compact(r: JSON) -> dict[str, Any]:
    """Return the fields used here from one E-utilities esummary record (the ClinVar record kept in the cache)."""
    g = r.get("germline_classification") or {}
    oc = r.get("oncogenicity_classification") or {}
    onc = oc.get("description") or None
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
        "oncogenicity_review": oc.get("review_status") or None,
        "somatic_impact": som,
        "spdi": spdi,
        "cdna": cdna,
        "loc": loc,  # assembly -> [chrom, start, stop] (1-based)
    }


def fetch(
    http: Http, cache: KeyValueCache, vcv_ids: Iterable[str], status: SourceStatus | None = None
) -> dict[str, dict[str, Any]]:
    """Fetch ClinVar records by VCV accession, from the cache (one week) or E-utilities.

    {'VCV000013961': {...}, ...}. Failed requests are recorded in `status` (source 'ClinVar').

    Args:
        http: HTTP client (network calls, 200 accessions per request).
        cache: Response cache; records fetched are written to it.
        vcv_ids: VCV accessions; anything else is ignored.
        status: Records which sources were used and which requests failed; None to not record.

    Returns:
        The records found, by accession: ``accession``, ``title``, ``classification``, ``review_status``,
        ``last_evaluated``, ``conditions``, ``oncogenicity``, ``oncogenicity_review``, ``somatic_impact``,
        ``spdi``, ``cdna`` and ``loc`` (assembly -> [chrom, start, stop], 1-based). Accessions of failed
        requests are missing.
    """
    ids = sorted({v for v in vcv_ids if v and v.startswith("VCV")})
    if status is not None and ids:
        status.use("ClinVar", len(ids))
    keys = {v: f"clinvar3|{v}" for v in ids}  # clinvar3: records with positions and oncogenicity review
    hit = cache.get_many(list(keys.values()), max_age=CLINVAR_MAX_AGE)
    out = {v: hit[keys[v]] for v in ids if keys[v] in hit}
    missing = [v for v in ids if v not in out]
    for i in range(0, len(missing), 200):
        chunk = missing[i : i + 200]
        uid = {str(int(v[3:])): v for v in chunk}
        try:
            r = http.request(
                "POST",
                EUTILS,
                data={"db": "clinvar", "id": ",".join(uid), "retmode": "json", "tool": "missensexplorer"},
            )
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


def _parse_spdi(spdi: str | None) -> tuple[str | None, int, str, str] | None:
    """Split a canonical SPDI into (chrom, 0-based position, deleted, inserted); None if empty or malformed.

    The chromosome is None when the sequence is not a human chromosome (NC_) accession.
    """
    if not spdi:
        return None
    try:
        seq, pos, dele, ins = spdi.split(":")
        return _NC_CHROM.get(seq.split(".")[0]), int(pos), dele, ins
    except ValueError:
        return None


def match(
    records: dict[str, dict[str, Any]],
    candidates: list[str],
    spdi_key: tuple[str, int, str, str] | None,
    mane: str | None,
    hgvsc: str | None,
) -> dict[str, Any] | None:
    """Pick the ClinVar record describing exactly this allele.

    Args:
        records: ClinVar records by accession, from :func:`fetch`.
        candidates: The row's VCV accessions, in order of preference.
        spdi_key: The allele as (chrom, 0-based position, deleted, inserted), compared with each record's
            canonical SPDI first; None if unknown.
        mane: MANE Select RefSeq transcript (NM_...), compared with the transcript in the record's title.
        hgvsc: The allele's c. change on that transcript.

    Returns:
        The first record matching by SPDI, else the first matching by transcript and c. change, else None.
    """
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
    """Return ClinVar's gold stars for a review status: 0 if unknown, None if there is no status."""
    return None if not review_status else STARS.get(review_status, 0)


def vcv(uid: str) -> str:
    """Return the VCV accession of an E-utilities ClinVar uid: '13961' -> 'VCV000013961'."""
    return f"VCV{int(uid):09d}"


def search_windows(
    http: Http,
    cache: KeyValueCache,
    windows: Iterable[tuple[str, int, int]],
    assembly: str,
    status: SourceStatus | None = None,
) -> dict[tuple[str, int, int], list[str]]:
    """Find the pathogenic ClinVar variants in genomic windows.

    VCV accessions of the pathogenic or likely pathogenic ClinVar variants that overlap each genomic
    window (chrom, start, end; 1-based, on `assembly`), {window: [VCV...]}. Searches are batched
    (WINDOWS_PER_SEARCH windows each) and the result of every window is cached for a week. Windows of
    a failed search are left out of the result and recorded in `status` (source 'ClinVar').

    Args:
        http: HTTP client (network calls to esearch, and to esummary through :func:`fetch`).
        cache: Response cache; windows and records found are written to it.
        windows: Genomic windows, as from :func:`missensexplorer.residue.window`.
        assembly: ``"GRCh38"`` or ``"GRCh37"``.
        status: Records which sources were used and which requests failed; None to not record.
    """
    return _search(http, cache, windows, assembly, PATHOGENIC_FILTER, "clinvar-window", status)


def search_oncogenicity_windows(
    http: Http,
    cache: KeyValueCache,
    windows: Iterable[tuple[str, int, int]],
    assembly: str,
    status: SourceStatus | None = None,
) -> dict[tuple[str, int, int], list[str]]:
    """Find the ClinVar variants with an oncogenicity classification in genomic windows.

    Like :func:`search_windows`, for any oncogenicity classification (the caller keeps the oncogenic
    ones); results are cached for a week under their own keys.
    """
    return _search(http, cache, windows, assembly, ONCOGENICITY_FILTER, "clinvar-onc-window", status)


def _search(
    http: Http,
    cache: KeyValueCache,
    windows: Iterable[tuple[str, int, int]],
    assembly: str,
    term_filter: str,
    key_prefix: str,
    status: SourceStatus | None,
) -> dict[tuple[str, int, int], list[str]]:
    """Search ClinVar for the variants matching `term_filter` in each window (batched, cached a week)."""
    field = "chrpos38" if assembly == "GRCh38" else "chrpos37"
    windows = sorted(set(windows))
    keys = {w: f"{key_prefix}|{assembly}|{w[0]}:{w[1]}-{w[2]}" for w in windows}
    hit = cache.get_many(list(keys.values()), max_age=CLINVAR_MAX_AGE)
    out = {w: hit[keys[w]] for w in windows if keys[w] in hit}
    todo = [w for w in windows if w not in out]
    for i in range(0, len(todo), WINDOWS_PER_SEARCH):
        chunk = todo[i : i + WINDOWS_PER_SEARCH]
        term = " OR ".join(f"({c}[chr] AND {a}:{b}[{field}])" for c, a, b in chunk)
        try:
            r = http.request(
                "POST",
                ESEARCH,
                data={
                    "db": "clinvar",
                    "term": f"({term}) AND {term_filter}",
                    "retmode": "json",
                    "retmax": 10000,
                    "tool": "missensexplorer",
                },
            )
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
            continue  # incomplete: do not cache these windows
        for w in chunk:
            out[w] = [v for v, rec in records.items() if _overlaps(rec, w, assembly)]
        cache.set_many({keys[w]: out[w] for w in chunk})
    return out


def _overlaps(rec: dict[str, Any], window: tuple[str, int, int], assembly: str) -> bool:
    """Whether a ClinVar record lies in a window (its location on `assembly`, else its GRCh38 SPDI)."""
    c, a, b = window
    loc = (rec.get("loc") or {}).get(assembly)
    if loc:
        return str(loc[0]) == str(c) and loc[1] <= b and loc[2] >= a
    p = _parse_spdi(rec.get("spdi")) if assembly == "GRCh38" else None
    return p is not None and p[0] == str(c) and a <= p[1] + 1 <= b
