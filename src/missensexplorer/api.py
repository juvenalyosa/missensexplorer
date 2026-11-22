"""HTTP access to Ensembl, UniProt and AlphaFold DB with retries and rate-limit handling."""

from __future__ import annotations

import logging
import random
import time

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


class Http:
    """Shared, thread-safe HTTP client (httpx.Client pools connections across threads)."""

    def __init__(self, timeout: float = 300, max_retries: int = 6):
        self.max_retries = max_retries
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
        last_err = ""
        for attempt in range(self.max_retries):
            try:
                r = self.client.request(method, url, **kwargs)
            except httpx.TransportError as e:
                last_err = f"{type(e).__name__}: {e}"
                self._sleep(attempt)
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
