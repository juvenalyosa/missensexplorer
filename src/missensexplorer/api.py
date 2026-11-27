"""HTTP access to Ensembl, UniProt and AlphaFold DB with retries and rate-limit handling."""

from __future__ import annotations

import logging
import random
import threading
import time
from collections.abc import Container, Iterable
from typing import Any
from urllib.parse import urlsplit

import httpx

from . import __version__
from ._types import JSON, KeyValueCache

log = logging.getLogger(__name__)

ENSEMBL = {
    "GRCh38": "https://rest.ensembl.org",
    "GRCh37": "https://grch37.rest.ensembl.org",
}
UNIPROT = "https://rest.uniprot.org"
ALPHAFOLD = "https://alphafold.ebi.ac.uk"


class ApiError(RuntimeError):
    """A web API request failed after its retries, or the API answered with an error.

    Args:
        message: Human-readable description of the failure.
        status: HTTP status code of the failed response; None when no response was received.

    Attributes:
        status: HTTP status code of the failed response; None when no response was received.
    """

    def __init__(self, message: str, status: int | None = None) -> None:
        super().__init__(message)
        self.status = status


class OfflineError(ApiError):
    """A host could not be reached at all (no network, DNS failure, firewall)."""


class Http:
    """Shared, thread-safe HTTP client (httpx.Client pools connections across threads).

    Server errors and rate limits are retried with backoff. Network errors are retried
    only a few times; after that the host is marked as down and later requests to it
    fail immediately, so an offline computer gets a clear error instead of long waits.
    A host whose requests keep failing after all retries (rate limits or server errors on
    `failures_before_down` requests in a row) is marked as down the same way: the analysis
    then finishes without that source, reported, instead of waiting on it for every item.

    Args:
        timeout: Read/write timeout of one request, in seconds (connecting times out after 20 s).
        max_retries: Attempts per request on rate limits (HTTP 429) and server errors (5xx).
        transport_retries: Attempts per request on network errors before giving up.
        failures_before_down: The host is marked as down once this many requests in a row
            have failed after all retries, the last one with an HTTP status (429 or 5xx).
        max_wait: Total back-off per request, in seconds; the request fails once the next
            wait would exceed it.

    Attributes:
        down: Host -> reason, for hosts marked as down (requests to them raise `OfflineError`).
        failures: Host -> number of requests in a row that failed after all retries.
        client: The underlying connection pool.
    """

    def __init__(
        self,
        timeout: float = 300,
        max_retries: int = 6,
        transport_retries: int = 3,
        failures_before_down: int = 3,
        max_wait: float = 90,
    ) -> None:
        self.max_retries = max_retries
        self.transport_retries = transport_retries
        self.failures_before_down = failures_before_down
        self.max_wait = max_wait  # total back-off per request, in seconds
        self.down: dict[str, str] = {}
        self.failures: dict[str, int] = {}
        self._lock = threading.Lock()
        self.client = httpx.Client(
            timeout=httpx.Timeout(timeout, connect=20),
            headers={
                "User-Agent": f"missensexplorer/{__version__}",
                "Accept": "application/json",
            },
            follow_redirects=True,
            limits=httpx.Limits(max_connections=16, max_keepalive_connections=8),
        )

    def request(self, method: str, url: str, *, ok_statuses: Container[int] = (200,), **kwargs: Any) -> httpx.Response:
        """Return the response for a 2xx/expected status; retry 429/5xx/network errors.

        4xx client errors (other than 429) are returned to the caller without
        retrying, because they usually carry a useful message (e.g. a bad HGVS).

        Args:
            method: HTTP method, e.g. "GET" or "POST".
            url: Full URL; its host is used for the down/failure bookkeeping.
            ok_statuses: Status codes returned without retrying, in addition to the 4xx codes
                other than 429 (e.g. a 404 that means "no entry").
            **kwargs: Passed to `httpx.Client.request` (params, json, data, headers, ...).

        Returns:
            The response, whose status is in `ok_statuses` or is a 4xx other than 429.

        Raises:
            OfflineError: The host is marked as down, or could not be connected to after
                `transport_retries` attempts (the host is then marked as down).
            ApiError: The request still failed after all retries or `max_wait` seconds of back-off.
        """
        host = urlsplit(url).netloc
        if host in self.down:
            raise OfflineError(f"{host} is unreachable ({self.down[host]})")
        last_err = ""
        transport_failures = 0
        waited = 0.0
        for attempt in range(self.max_retries):
            try:
                r = self.client.request(method, url, **kwargs)
            except httpx.TransportError as e:
                last_err = f"{type(e).__name__}: {e}"
                transport_failures += 1
                if transport_failures >= self.transport_retries:
                    if isinstance(e, (httpx.ConnectError, httpx.ConnectTimeout)):
                        with self._lock:
                            self.down[host] = last_err
                        raise OfflineError(f"cannot connect to {host} ({last_err})") from e
                    break
                self._sleep(transport_failures - 1)
                continue
            if r.status_code in ok_statuses or (r.status_code < 500 and r.status_code != 429):
                with self._lock:
                    self.failures[host] = 0
                return r  # expected status, or another 4xx for the caller to read
            if r.status_code == 429:
                last_err = "HTTP 429, rate limited"
                try:
                    wait = float(r.headers.get("Retry-After") or 0) or 2**attempt
                except ValueError:
                    wait = 2**attempt
            else:
                last_err = f"HTTP {r.status_code}"
                wait = min(60, 2**attempt + random.random())
            if attempt == self.max_retries - 1 or waited + wait > self.max_wait:
                break  # no point waiting when no attempt follows
            log.info("%s from %s, waiting %.1fs", last_err, url, wait)
            self._pause(min(wait, 120))
            waited += wait
        with self._lock:
            self.failures[host] = self.failures.get(host, 0) + 1
            if self.failures[host] >= self.failures_before_down and last_err.startswith("HTTP"):
                self.down[host] = f"keeps failing: {last_err}"
        raise ApiError(f"{method} {url} failed after retries ({last_err or 'no response'})")

    def get_json(self, url: str, **kwargs: Any) -> JSON:
        """GET `url` and return its parsed JSON body.

        Args:
            url: Full URL.
            **kwargs: Passed to `request` (e.g. params, headers, ok_statuses).

        Raises:
            ApiError: The status is not 200 (the message is the API's error text), or the
                request failed (see `request`).
        """
        r = self.request("GET", url, **kwargs)
        if r.status_code != 200:
            raise ApiError(_error_text(r), r.status_code)
        return r.json()

    def post_json(self, url: str, payload: dict[str, Any], **kwargs: Any) -> JSON:
        """POST `payload` as JSON to `url` and return the parsed JSON answer.

        Args:
            url: Full URL.
            payload: JSON-serialisable request body.
            **kwargs: Passed to `request` (e.g. params, headers, ok_statuses).

        Raises:
            ApiError: The status is not 200 (the message is the API's error text), or the
                request failed (see `request`).
        """
        r = self.request("POST", url, json=payload, **kwargs)
        if r.status_code != 200:
            raise ApiError(_error_text(r), r.status_code)
        return r.json()

    @staticmethod
    def _pause(seconds: float) -> None:
        """Sleep before retrying a rate-limited or failed request (tests replace it)."""
        time.sleep(seconds)

    @staticmethod
    def _sleep(attempt: int) -> None:
        """Sleep with exponential back-off and jitter after a network error (tests replace it)."""
        time.sleep(min(60, (2**attempt) + random.random()))

    def close(self) -> None:
        """Close the pooled connections."""
        self.client.close()


class SourceStatus:
    """Which external sources were used, and which requests to them failed.

    A failed source is not the same as "no data": rows affected by a failure are
    flagged so the user knows the gap comes from the service and can retry.

    `use`, `skip` and `fail` may be called from several threads at once.

    Attributes:
        used: Source -> number of items looked up.
        errors: Source -> first distinct error messages (at most 5).
        failed_items: Source -> items (variant queries, accessions, ...) whose request failed.
        skipped: Source -> reason the source was not queried at all (e.g. offline).
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self.used: dict[str, int] = {}
        self.errors: dict[str, list[str]] = {}
        self.failed_items: dict[str, set[str]] = {}
        self.skipped: dict[str, str] = {}

    def use(self, source: str, n: int = 1) -> None:
        """Record that `n` items were looked up in `source`."""
        with self._lock:
            self.used[source] = self.used.get(source, 0) + n

    def skip(self, source: str, reason: str) -> None:
        """Record that `source` was not queried, and why."""
        with self._lock:
            self.skipped[source] = reason

    def fail(self, source: str, message: str, items: Iterable[str | None] = ()) -> None:
        """Record a failed request to `source`.

        Args:
            source: Source name, as shown in the run summary.
            message: Error message; only the first 5 distinct messages are kept.
            items: Items the failed request was for; None entries are ignored. With no items
                the whole source counts as failed (see `failed`).
        """
        with self._lock:
            msgs = self.errors.setdefault(source, [])
            if len(msgs) < 5 and str(message) not in msgs:
                msgs.append(str(message))
            self.failed_items.setdefault(source, set()).update(i for i in items if i is not None)

    def failed(self, source: str, item: str | None = None) -> bool:
        """Return whether `source` had a failure (for `item`, when given)."""
        if source not in self.errors:
            return False
        return item is None or item in self.failed_items.get(source, ())

    def summary(self) -> dict[str, str]:
        """Return source -> "ok", "failed (n items): first error" or "skipped: reason", sorted by source."""
        out = {}
        for src in sorted(set(self.used) | set(self.errors) | set(self.skipped)):
            if src in self.errors:
                n = len(self.failed_items.get(src, ()))
                out[src] = f"failed{f' ({n} items)' if n else ''}: {self.errors[src][0]}"
            elif src in self.skipped:
                out[src] = f"skipped: {self.skipped[src]}"
            else:
                out[src] = "ok"
        return out


def ensembl_release(http: Http, cache: KeyValueCache, assembly: str, timeout: float | None = None) -> int | None:
    """Return the current Ensembl release of the REST server (cached for a day).

    Args:
        http: Client used for the request (GET /info/software).
        cache: Cache read and written under "ensembl_release|<assembly>".
        assembly: "GRCh38" or "GRCh37" (a key of `ENSEMBL`).
        timeout: Read timeout of the request, in seconds; None for the client's.

    Returns:
        The release number, or None when the server does not report one.

    Raises:
        ApiError: The request failed (OfflineError when the server cannot be reached).
    """
    key = f"ensembl_release|{assembly}"
    hit = cache.get(key, max_age=86400)
    if hit:
        return hit
    extra = {"timeout": httpx.Timeout(timeout, connect=20)} if timeout else {}
    data = http.get_json(f"{ENSEMBL[assembly]}/info/software", params={"content-type": "application/json"}, **extra)
    rel = data.get("release")
    if rel:
        cache.set(key, rel)
    return rel


def _error_text(r: httpx.Response) -> str:
    """Return the API's own error message from an error response, else the status and body start."""
    try:
        data = r.json()
        if isinstance(data, dict) and "error" in data:
            return str(data["error"])
        if isinstance(data, dict) and "messages" in data:
            return "; ".join(map(str, data["messages"]))
    except ValueError:
        pass
    return f"HTTP {r.status_code}: {r.text[:300]}"
