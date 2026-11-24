"""Batched, cached, parallel calls to the Ensembl VEP REST API."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import re
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

import httpx

from ._types import JSON, KeyValueCache, Progress
from .api import ENSEMBL, ApiError, Http, SourceStatus

log = logging.getLogger(__name__)

# Always requested: transcript choice, HGVS, protein ids, UniProt mapping, exon numbers.
# pick_allele_gene => one consequence block per allele per overlapping gene, chosen with
# VEP's ranking (MANE Select > canonical > APPRIS > TSL > biotype > ...).
BASE_PARAMS: dict[str, int] = {
    "pick_allele_gene": 1,
    "hgvs": 1,
    "mane": 1,
    "canonical": 1,
    "uniprot": 1,
    "numbers": 1,
    "variant_class": 1,
    "protein": 1,
}
# Pathogenicity predictors served by the Ensembl REST plugins (only requested for
# protein-altering/splice variants on large inputs because they are slower).
PREDICTOR_PARAMS: dict[str, int | str] = {
    "AlphaMissense": 1,
    "CADD": 1,
    "EVE": 1,
    "SpliceAI": 1,
    "dbNSFP": "REVEL_score",
}

ENDPOINTS = {  # kind -> (path, payload key, batch size)
    "region": ("/vep/human/region", "variants", 200),
    "hgvs": ("/vep/human/hgvs", "hgvs_notations", 100),
    "id": ("/vep/human/id", "ids", 200),
}

# Consequences for which the predictor pass is worth running.
DETAIL_CONSEQUENCES = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "synonymous_variant",
    "stop_retained_variant",
    "start_retained_variant",
    "incomplete_terminal_codon_variant",
    "coding_sequence_variant",
}


def _sig(params: Mapping[str, int | str]) -> str:
    """Return a short hash of the request parameters, so cache keys differ per parameter set."""
    return hashlib.sha1(json.dumps(params, sort_keys=True).encode()).hexdigest()[:10]


class Vep:
    """Ensembl VEP REST client.

    Cache keys include the Ensembl release, so answers from an
    older release are not reused after the server is updated.

    Args:
        http: Client used for the requests.
        cache: Cache of the answers per query (failed requests are not cached).
        assembly: "GRCh38" or "GRCh37" (a key of `ENSEMBL`).
        workers: Number of batches sent in parallel.
        release: Ensembl release of the server (see `ensembl_release`); None when unknown.
        status: Where failed requests are recorded; None to not record them.

    Attributes:
        name: Display name of the annotation source.
        base: Base URL of the Ensembl REST server for `assembly`.

    Raises:
        ValueError: `assembly` is not supported.
    """

    name = "Ensembl VEP (REST)"

    def __init__(
        self,
        http: Http,
        cache: KeyValueCache,
        assembly: str = "GRCh38",
        workers: int = 4,
        release: int | None = None,
        status: SourceStatus | None = None,
    ) -> None:
        if assembly not in ENSEMBL:
            raise ValueError(f"assembly must be one of {list(ENSEMBL)}")
        self.http, self.cache, self.assembly, self.workers = http, cache, assembly, workers
        self.base = ENSEMBL[assembly]
        self.release = release
        self.status = status

    def _key(self, kind: str, params: Mapping[str, int | str], q: str) -> str:
        """Return the cache key of one query."""
        return f"vep|e{self.release or '?'}|{self.assembly}|{kind}|{_sig(params)}|{q}"

    def annotate(
        self,
        queries: dict[str, list[str]],
        params: Mapping[str, int | str],
        progress: Progress | None = None,
        stage: str = "VEP",
        source: str = "Ensembl VEP",
    ) -> dict[tuple[str, str], dict[str, Any]]:
        """Annotate variant queries with VEP, in parallel batches, reusing cached answers.

        queries: {kind: [query, ...]} -> {(kind, query): {"results": [...]} | {"error": str}}.

        A batch that still fails after the retries gives {"error": ..., "service_error": True}
        for its queries (not cached, so the next run asks again) and is recorded in `status`.

        New answers with results, and errors that will not change on a retry
        ({"error": ..., "cache_error": True}, e.g. a notation VEP rejects), are written to
        the cache. Duplicate queries are sent once.

        Args:
            queries: Kind ("region", "hgvs" or "id", see `ENDPOINTS`) -> query strings.
            params: VEP request parameters, e.g. `BASE_PARAMS` or `BASE_PARAMS` plus `PREDICTOR_PARAMS`.
            progress: Called with (stage, queries done, queries to send).
            stage: Stage name passed to `progress`.
            source: Source name used in `status`.

        Returns:
            (kind, query) -> {"results": [VEP JSON object per consequence block]}, or
            {"error": message} with "cache_error" or "service_error" set to True.
        """
        out: dict[tuple[str, str], dict[str, Any]] = {}
        todo: list[tuple[str, list[str]]] = []
        for kind, qs in queries.items():
            qs = list(dict.fromkeys(qs))
            keys = [self._key(kind, params, q) for q in qs]
            hit = self.cache.get_many(keys)
            missing = []
            for q, k in zip(qs, keys, strict=True):
                if k in hit:
                    out[(kind, q)] = hit[k]
                else:
                    missing.append(q)
            # small jobs: spread over the workers so requests run in parallel
            size = min(ENDPOINTS[kind][2], max(1, math.ceil(len(missing) / self.workers)))
            todo += [(kind, missing[i : i + size]) for i in range(0, len(missing), size)]

        total = sum(len(b) for _, b in todo)
        if self.status is not None:
            self.status.use(source, sum(len(set(q)) for q in queries.values()))
        done = 0
        if progress:
            progress(stage, 0, total)
        if not todo:
            return out
        with ThreadPoolExecutor(max_workers=self.workers) as ex:
            futs = {ex.submit(self._run_batch, kind, batch, params): (kind, batch) for kind, batch in todo}
            for fut in as_completed(futs):
                kind, batch = futs[fut]
                try:
                    res = fut.result()
                except Exception as e:  # network failure after retries
                    log.warning("VEP batch failed: %s", e)
                    res = {q: {"error": f"VEP request failed: {e}", "service_error": True} for q in batch}
                    if self.status is not None:
                        self.status.fail(source, f"{type(e).__name__}: {e}", batch)
                cacheable = {}
                for q, val in res.items():
                    out[(kind, q)] = val
                    if "results" in val or val.get("cache_error"):
                        cacheable[self._key(kind, params, q)] = val
                self.cache.set_many(cacheable)
                done += len(batch)
                if progress:
                    progress(stage, done, total)
        return out

    # ------------------------------------------------------------------ internals
    def _post(self, kind: str, batch: list[str], params: Mapping[str, int | str]) -> httpx.Response:
        """POST one batch of queries to the VEP endpoint of `kind`."""
        path, key, _ = ENDPOINTS[kind]
        return self.http.request(
            "POST", self.base + path, params=params, json={key: batch}, headers={"Content-Type": "application/json"}
        )

    def _run_batch(
        self, kind: str, batch: list[str], params: Mapping[str, int | str], _attempt: int = 0
    ) -> dict[str, dict[str, Any]]:
        """Annotate one batch, retrying unreadable answers and splitting the batch to isolate bad queries.

        Returns:
            Query -> answer, as in `annotate`.

        Raises:
            ApiError: The request failed (other than with HTTP 400).
        """
        r = self._post(kind, batch, params)
        if r.status_code == 200:
            items = _parse_json(r.text)
            if items is None:  # truncated/garbled response: retry, then split the batch
                log.warning("VEP returned an unreadable response (%d variants, attempt %d)", len(batch), _attempt + 1)
                if _attempt < 2:
                    return self._run_batch(kind, batch, params, _attempt + 1)
                if len(batch) > 1:
                    mid = len(batch) // 2
                    res = self._run_batch(kind, batch[:mid], params)
                    res.update(self._run_batch(kind, batch[mid:], params))
                    return res
                # a single variant whose answer cannot be read: lose only this one
                if self.status is not None:
                    self.status.fail("Ensembl VEP", "unreadable response", batch)
                return {
                    batch[0]: {"error": "VEP returned an unreadable response for this variant", "service_error": True}
                }
        if r.status_code == 400:
            if len(batch) == 1:
                try:
                    body = r.json()
                except ValueError:
                    body = None
                msg = body.get("error", r.text) if isinstance(body, dict) else r.text
                return {batch[0]: {"error": _short_error(msg), "cache_error": True}}
            mid = len(batch) // 2  # isolate the offending entry
            res = self._run_batch(kind, batch[:mid], params)
            res.update(self._run_batch(kind, batch[mid:], params))
            return res
        if r.status_code != 200:
            raise ApiError(f"HTTP {r.status_code}: {r.text[:200]}", r.status_code)

        grouped: dict[str, list[JSON]] = {}
        for item in items:
            grouped.setdefault(item.get("input"), []).append(item)
        res = {}
        lost = []
        for q in batch:
            if q in grouped:
                res[q] = {"results": grouped[q]}
            else:
                lost.append(q)
        # VEP silently drops inputs it cannot parse/map; ask individually for the reason.
        for q in lost[:25] if len(batch) > 1 else []:
            res.update(self._run_batch(kind, [q], params))
        for q in lost:
            res.setdefault(
                q, {"error": "VEP could not map this variant (check notation/assembly)", "cache_error": True}
            )
        return res


_BAD_NUMBER = re.compile(r"(?<=[:\[,])(\s*)-?(?:nan|NaN|inf|Inf|infinity|Infinity)(?=\s*[,}\]])")


def _parse_json(text: str) -> JSON:
    """json.loads, tolerating bare nan/inf tokens that some VEP plugins print; None if unreadable."""
    try:
        return json.loads(text)
    except ValueError:
        pass
    try:
        return json.loads(_BAD_NUMBER.sub(r"\1null", text))
    except ValueError:
        return None


def _short_error(msg: str) -> str:
    """Return a short, readable version of a VEP error message for the results table."""
    msg = str(msg).strip()
    if "Could not get a Transcript object" in msg:
        return "unknown gene or transcript symbol"
    if "does not match reference" in msg.lower() or "reference allele" in msg.lower():
        return "reference allele does not match the genome/transcript"
    return msg.split(" at /")[0][:240]
