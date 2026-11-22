"""Batched, cached, parallel calls to the Ensembl VEP REST API."""

from __future__ import annotations

import hashlib
import json
import logging
import math
from concurrent.futures import ThreadPoolExecutor, as_completed

from .api import ENSEMBL, ApiError, Http

log = logging.getLogger(__name__)

# Always requested: transcript choice, HGVS, protein ids, UniProt mapping, exon numbers.
# pick_allele_gene => one consequence block per allele per overlapping gene, chosen with
# VEP's ranking (MANE Select > canonical > APPRIS > TSL > biotype > ...).
BASE_PARAMS = {
    "pick_allele_gene": 1, "hgvs": 1, "mane": 1, "canonical": 1, "uniprot": 1,
    "numbers": 1, "variant_class": 1, "protein": 1,
}
# Pathogenicity predictors served by the Ensembl REST plugins (only requested for
# protein-altering/splice variants on large inputs because they are slower).
PREDICTOR_PARAMS = {
    "AlphaMissense": 1, "CADD": 1, "EVE": 1, "SpliceAI": 1, "dbNSFP": "REVEL_score",
}

ENDPOINTS = {  # kind -> (path, payload key, batch size)
    "region": ("/vep/human/region", "variants", 200),
    "hgvs": ("/vep/human/hgvs", "hgvs_notations", 100),
    "id": ("/vep/human/id", "ids", 200),
}

# Consequences for which the predictor pass is worth running.
DETAIL_CONSEQUENCES = {
    "missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
    "inframe_insertion", "inframe_deletion", "protein_altering_variant",
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "splice_donor_5th_base_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant", "synonymous_variant", "stop_retained_variant",
    "start_retained_variant", "incomplete_terminal_codon_variant", "coding_sequence_variant",
}


def _sig(params: dict) -> str:
    return hashlib.sha1(json.dumps(params, sort_keys=True).encode()).hexdigest()[:10]


class Vep:
    def __init__(self, http: Http, cache, assembly: str = "GRCh38", workers: int = 4):
        if assembly not in ENSEMBL:
            raise ValueError(f"assembly must be one of {list(ENSEMBL)}")
        self.http, self.cache, self.assembly, self.workers = http, cache, assembly, workers
        self.base = ENSEMBL[assembly]

    def _key(self, kind, params, q):
        return f"vep|{self.assembly}|{kind}|{_sig(params)}|{q}"

    def annotate(self, queries: dict[str, list[str]], params: dict, progress=None,
                 stage: str = "VEP") -> dict[tuple[str, str], dict]:
        """queries: {kind: [query, ...]} -> {(kind, query): {"results": [...]} | {"error": str}}."""
        out: dict[tuple[str, str], dict] = {}
        todo: list[tuple[str, list[str]]] = []
        for kind, qs in queries.items():
            qs = list(dict.fromkeys(qs))
            keys = [self._key(kind, params, q) for q in qs]
            hit = self.cache.get_many(keys)
            missing = []
            for q, k in zip(qs, keys):
                if k in hit:
                    out[(kind, q)] = hit[k]
                else:
                    missing.append(q)
            # small jobs: spread over the workers so requests run in parallel
            size = min(ENDPOINTS[kind][2], max(1, math.ceil(len(missing) / self.workers)))
            todo += [(kind, missing[i : i + size]) for i in range(0, len(missing), size)]

        total = sum(len(b) for _, b in todo)
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
                    res = {q: {"error": f"VEP request failed: {e}"} for q in batch}
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
    def _post(self, kind, batch, params):
        path, key, _ = ENDPOINTS[kind]
        r = self.http.request("POST", self.base + path, params=params, json={key: batch},
                              headers={"Content-Type": "application/json"})
        return r

    def _run_batch(self, kind, batch, params) -> dict[str, dict]:
        r = self._post(kind, batch, params)
        if r.status_code == 400:
            if len(batch) == 1:
                try:
                    msg = r.json().get("error", r.text)
                except ValueError:
                    msg = r.text
                return {batch[0]: {"error": _short_error(msg), "cache_error": True}}
            mid = len(batch) // 2  # isolate the offending entry
            res = self._run_batch(kind, batch[:mid], params)
            res.update(self._run_batch(kind, batch[mid:], params))
            return res
        if r.status_code != 200:
            raise ApiError(f"HTTP {r.status_code}: {r.text[:200]}", r.status_code)

        grouped: dict[str, list] = {}
        for item in r.json():
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
            res.setdefault(q, {"error": "VEP could not map this variant (check notation/assembly)", "cache_error": True})
        return res


def _short_error(msg: str) -> str:
    msg = str(msg).strip()
    if "Could not get a Transcript object" in msg:
        return "unknown gene or transcript symbol"
    if "does not match reference" in msg.lower() or "reference allele" in msg.lower():
        return "reference allele does not match the genome/transcript"
    return msg.split(" at /")[0][:240]
