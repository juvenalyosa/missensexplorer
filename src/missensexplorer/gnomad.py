"""Direct gnomAD check (GraphQL API) for variants without a frequency in Ensembl's annotation.

Ensembl attaches gnomAD frequencies only to variants in its variation database (dbSNP), so many rare
gnomAD variants, singletons above all, come back without one. For those, this module asks gnomAD
itself (v4 on GRCh38, v2.1 on GRCh37): allele count and number in exomes and genomes, and the
highest frequency among the main genetic ancestry groups (grpmax, as acmg.GRPMAX_GROUPS). A variant
gnomAD does not know is "absent" only where gnomAD has reads: the share of exome samples with at least
20 reads at the position is fetched too (the exomes are 730,947 of gnomAD v4's 807,162 individuals).

The public API limits the cost of a query (25: one variant costs 1, one position's exome coverage 6)
and the request rate (about 10 per minute), so lookups are batched to that cost and paced.
"""

from __future__ import annotations

import logging
import threading
import time
from collections import deque
from collections.abc import Callable, Iterable
from typing import Any

from ._types import JSON, KeyValueCache
from .api import ApiError, Http, SourceStatus

log = logging.getLogger(__name__)

API = "https://gnomad.broadinstitute.org/api"
DATASETS = {"GRCh38": ("gnomad_r4", "gnomAD v4"), "GRCh37": ("gnomad_r2_1", "gnomAD v2.1")}
GROUPS = ("afr", "amr", "eas", "nfe", "sas")
PER_REQUEST = 25  # variants per query (cost 1 each; the API allows 25)
COVERAGE_PER_REQUEST = 4  # positions per coverage query (exome coverage costs 6 each)
RATE = (10, 60.0)  # at most 10 requests per 60 s
MAX_AGE = 90 * 86400
_sent: deque[float] = deque()  # monotonic times of the requests sent in the last RATE window
_lock = threading.Lock()
_FIELDS = "ac an filters populations { id ac an }"


def _vid(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Return gnomAD's variant id, such as '7-140753336-A-T'."""
    return f"{str(chrom).removeprefix('chr')}-{pos}-{ref}-{alt}"


def _pace() -> None:
    """Wait as needed to stay within the API's request rate."""
    n, window = RATE
    with _lock:
        now = time.monotonic()
        while _sent and now - _sent[0] > window:
            _sent.popleft()
        if len(_sent) >= n:
            time.sleep(window - (now - _sent[0]) + 0.5)
            _sent.popleft()
        _sent.append(time.monotonic())


def _post(http: Http, query: str) -> dict[str, JSON]:
    """Send one GraphQL query (paced) and return its ``data``; raise ApiError if the API returned none."""
    _pace()
    r = http.request("POST", API, json={"query": query}, ok_statuses=(200, 400))
    if r.status_code != 200:
        raise ApiError(f"HTTP {r.status_code}", r.status_code)
    d = r.json()
    if d.get("data") is None:
        raise ApiError(f"gnomAD API error: {(d.get('errors') or [{}])[0].get('message')}")
    return d["data"]  # "Variant not found" errors leave the absent variants at null


def _summary(v: JSON | None) -> dict[str, Any]:
    """Summarise one ``variant`` answer (exome and genome together); ``found`` False when gnomAD lacks it."""
    if not v:
        return {"found": False}
    ac = an = 0
    group = {g: [0, 0] for g in GROUPS}
    passed = False
    for part in (v.get("exome"), v.get("genome")):
        if not part:
            continue
        ac, an = ac + (part.get("ac") or 0), an + (part.get("an") or 0)
        passed = passed or not part.get("filters")
        for p in part.get("populations") or []:
            if p.get("id") in group:
                group[p["id"]][0] += p.get("ac") or 0
                group[p["id"]][1] += p.get("an") or 0
    afs = [a / n for a, n in group.values() if n]
    return {
        "found": ac > 0,
        "ac": ac,
        "an": an,
        "af": ac / an if an else None,
        "grpmax": max(afs) if afs and ac else None,
        "passes_filters": passed,
    }


def _coverage(region: JSON | None, pos: int) -> dict[str, float | None]:
    """Return the share of exome samples with >= 20 reads at ``pos``, from the nearest bin within 10 bases."""
    bins = ((region or {}).get("coverage") or {}).get("exome") or []
    near = min(bins, key=lambda b: abs(b["pos"] - pos), default=None)
    return {"exome": None if near is None or abs(near["pos"] - pos) > 10 else near.get("over_20")}


def _variants_query(chunk: list[tuple[str, int, str, str]], dataset: str) -> str:
    """Return one GraphQL query for the exome and genome counts of several variants (aliases v0, v1, ...)."""
    return (
        "{ "
        + " ".join(
            f'v{j}: variant(variantId: "{_vid(*v)}", dataset: {dataset}) {{ exome {{ {_FIELDS} }} '
            f"genome {{ {_FIELDS} }} }}"
            for j, v in enumerate(chunk)
        )
        + " }"
    )


def _coverage_query(part: list[tuple[str, int, str, str]], dataset: str, ref_genome: str) -> str:
    """Return one GraphQL query for the exome coverage at several variants' positions (aliases c0, c1, ...)."""
    return (
        "{ "
        + " ".join(
            f'c{j}: region(chrom: "{str(v[0]).removeprefix("chr")}", start: {v[1]}, stop: {v[1]}, '
            f"reference_genome: {ref_genome}) "
            f"{{ coverage(dataset: {dataset}) {{ exome {{ pos over_20 }} }} }}"
            for j, v in enumerate(part)
        )
        + " }"
    )


def lookup(
    http: Http,
    cache: KeyValueCache,
    variants: Iterable[tuple[str, int, str, str]],
    assembly: str,
    status: SourceStatus | None = None,
    progress: Callable[[int, int], None] | None = None,
) -> dict[tuple[str, int, str, str], dict[str, Any]]:
    """Look up variants in gnomAD itself, from the cache (90 days) or the GraphQL API (paced).

    {(chrom, pos, ref, alt): {found, ac, an, af, grpmax, passes_filters, dataset, and for absent
    variants over_20: {exome}}} for SNVs and VCF-style indels. Failures are recorded in `status`
    (source 'gnomAD'); `progress(done, total)` is called after each batch.

    Args:
        http: HTTP client (network calls).
        cache: Response cache; the answers are written to it.
        variants: VCF-style alleles ``(chrom, pos, ref, alt)`` on `assembly`.
        assembly: ``"GRCh38"`` (gnomAD v4) or ``"GRCh37"`` (gnomAD v2.1); any other gives ``{}``.
        status: Records which sources were used and which requests failed; None to not record.
        progress: Called as ``progress(done, total)`` after each batch; None for no calls.

    Returns:
        The answers by variant; variants of a failed batch are missing.
    """
    if assembly not in DATASETS:
        return {}
    dataset, label = DATASETS[assembly]
    keys = {v: f"gnomad|{dataset}|{_vid(*v)}" for v in set(variants)}
    hit = cache.get_many(list(keys.values()), max_age=MAX_AGE)
    out = {v: hit[k] for v, k in keys.items() if k in hit}
    todo = sorted(v for v in keys if v not in out)
    if status is not None and todo:
        status.use("gnomAD", len(todo))
    for i in range(0, len(todo), PER_REQUEST):
        chunk = todo[i : i + PER_REQUEST]
        try:
            data = _post(http, _variants_query(chunk, dataset))
            res = {v: {**_summary(data.get(f"v{j}")), "dataset": label} for j, v in enumerate(chunk)}
            absent = [v for v in chunk if not res[v]["found"]]
            ref_genome = "GRCh38" if assembly == "GRCh38" else "GRCh37"
            for k in range(0, len(absent), COVERAGE_PER_REQUEST):
                part = absent[k : k + COVERAGE_PER_REQUEST]
                cov = _post(http, _coverage_query(part, dataset, ref_genome))
                for j, v in enumerate(part):
                    res[v]["over_20"] = _coverage(cov.get(f"c{j}"), v[1])
        except (ApiError, ValueError, KeyError, TypeError) as e:
            log.warning("gnomAD request failed: %s", e)
            if status is not None:
                status.fail("gnomAD", str(e), [_vid(*v) for v in chunk])
            continue
        out.update(res)
        cache.set_many({keys[v]: res[v] for v in chunk})
        if progress:
            progress(min(i + PER_REQUEST, len(todo)), len(todo))
    return out


def covered(info: dict[str, Any], share: float = 0.5) -> bool:
    """At least `share` of the exome samples have >= 20 reads at the position.

    Args:
        info: A variant's entry from :func:`lookup`.
        share: Minimum share of exome samples (0-1).

    Returns:
        False also when the coverage is unknown.
    """
    return ((info.get("over_20") or {}).get("exome") or 0) >= share
