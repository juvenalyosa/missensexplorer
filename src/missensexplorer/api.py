"""HTTP access to Ensembl, UniProt and AlphaFold DB with retries and rate-limit handling."""

from __future__ import annotations

import logging
import random
import threading
import time
from urllib.parse import urlsplit

import httpx

from . import __version__

log = logging.getLogger(__name__)

ENSEMBL = {
    "GRCh38": "https://rest.ensembl.org",
    "GRCh37": "https://grch37.rest.ensembl.org",
}
UNIPROT = "https://rest.uniprot.org"
ALPHAFOLD = "https://alphafold.ebi.ac.uk"


class ApiError(RuntimeError):
    def __init__(self, message: str, status: int | None = None):
        super().__init__(message)
        self.status = status


class OfflineError(ApiError):
    """A host could not be reached at all (no network, DNS failure, firewall)."""


class Http:
    """Shared, thread-safe HTTP client (httpx.Client pools connections across threads).

    Server errors and rate limits are retried with backoff. Network errors are retried
    only a few times; after that the host is marked as down and later requests to it
    fail immediately, so an offline computer gets a clear error instead of long waits.
    """

    def __init__(self, timeout: float = 300, max_retries: int = 6, transport_retries: int = 3):
        self.max_retries = max_retries
        self.transport_retries = transport_retries
        self.down: dict[str, str] = {}
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

    def request(self, method: str, url: str, *, ok_statuses=(200,), **kwargs) -> httpx.Response:
        """Return the response for a 2xx/expected status; retry 429/5xx/network errors.

        4xx client errors (other than 429) are returned to the caller without
        retrying, because they usually carry a useful message (e.g. a bad HGVS).
        """
        host = urlsplit(url).netloc
        if host in self.down:
            raise OfflineError(f"{host} is unreachable ({self.down[host]})")
        last_err = ""
        transport_failures = 0
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
            if r.status_code in ok_statuses:
                return r
            if r.status_code == 429:
                wait = float(r.headers.get("Retry-After", 0) or 0) or (2 ** attempt)
                log.info("Rate limited by %s, waiting %.1fs", url, wait)
                time.sleep(min(wait, 120))
                continue
            if r.status_code >= 500:
                last_err = f"HTTP {r.status_code}"
                self._sleep(attempt)
                continue
            return r  # other 4xx: let caller decide
        raise ApiError(f"{method} {url} failed after {self.max_retries} attempts ({last_err})")

    def get_json(self, url: str, **kwargs):
        r = self.request("GET", url, **kwargs)
        if r.status_code != 200:
            raise ApiError(_error_text(r), r.status_code)
        return r.json()

    def post_json(self, url: str, payload: dict, **kwargs):
        r = self.request("POST", url, json=payload, **kwargs)
        if r.status_code != 200:
            raise ApiError(_error_text(r), r.status_code)
        return r.json()

    @staticmethod
    def _sleep(attempt: int) -> None:
        time.sleep(min(60, (2 ** attempt) + random.random()))

    def close(self) -> None:
        self.client.close()


class SourceStatus:
    """Which external sources were used, and which requests to them failed.

    A failed source is not the same as "no data": rows affected by a failure are
    flagged so the user knows the gap comes from the service and can retry.
    """

    def __init__(self):
        self._lock = threading.Lock()
        self.used: dict[str, int] = {}
        self.errors: dict[str, list[str]] = {}
        self.failed_items: dict[str, set] = {}
        self.skipped: dict[str, str] = {}

    def use(self, source: str, n: int = 1) -> None:
        with self._lock:
            self.used[source] = self.used.get(source, 0) + n

    def skip(self, source: str, reason: str) -> None:
        with self._lock:
            self.skipped[source] = reason

    def fail(self, source: str, message: str, items=()) -> None:
        with self._lock:
            msgs = self.errors.setdefault(source, [])
            if len(msgs) < 5 and str(message) not in msgs:
                msgs.append(str(message))
            self.failed_items.setdefault(source, set()).update(i for i in items if i is not None)

    def failed(self, source: str, item=None) -> bool:
        if source not in self.errors:
            return False
        return item is None or item in self.failed_items.get(source, ())

    def summary(self) -> dict[str, str]:
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


def ensembl_release(http: Http, cache, assembly: str) -> int | None:
    """Current Ensembl release of the REST server (cached for a day)."""
    key = f"ensembl_release|{assembly}"
    hit = cache.get(key, max_age=86400)
    if hit:
        return hit
    data = http.get_json(f"{ENSEMBL[assembly]}/info/software", params={"content-type": "application/json"})
    rel = data.get("release")
    if rel:
        cache.set(key, rel)
    return rel


def _error_text(r: httpx.Response) -> str:
    try:
        data = r.json()
        if isinstance(data, dict) and "error" in data:
            return str(data["error"])
        if isinstance(data, dict) and "messages" in data:
            return "; ".join(map(str, data["messages"]))
    except ValueError:
        pass
    return f"HTTP {r.status_code}: {r.text[:300]}"
