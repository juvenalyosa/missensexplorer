"""Tiny thread-safe SQLite key/value cache for API responses.

Re-running the same variants (or variants in the same proteins) is then
instant and does not hit the public APIs again.
"""

from __future__ import annotations

import json
import os
import sqlite3
import threading
import time
from pathlib import Path
from typing import Any


def default_cache_dir() -> Path:
    """Return the cache directory (not created here).

    The first that applies: ``$MISSENSEXPLORER_CACHE``, ``$XDG_CACHE_HOME/missensexplorer``,
    ``~/.cache/missensexplorer``.
    """
    base = os.environ.get("MISSENSEXPLORER_CACHE")
    if base:
        return Path(base)
    xdg = os.environ.get("XDG_CACHE_HOME")
    return Path(xdg) / "missensexplorer" if xdg else Path.home() / ".cache" / "missensexplorer"


class Cache:
    """SQLite-backed key/value store of JSON values, with an age limit per read.

    Opening the cache creates the database file (and its directory) if needed. One
    connection is shared by all threads, serialised by a lock. Use it as a context manager,
    or call `close`, to release the connection.

    Args:
        path: Database file; None means ``default_cache_dir() / "cache.sqlite"``.
        ttl_days: Default maximum age of a value returned by `get`/`get_many`, in days
            (0 means no limit).

    Attributes:
        path: The database file.
        ttl: Default maximum age, in seconds.
    """

    def __init__(self, path: str | Path | None = None, ttl_days: float = 30) -> None:
        path = Path(path) if path else default_cache_dir() / "cache.sqlite"
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.ttl = ttl_days * 86400
        self._lock = threading.Lock()
        self._db = sqlite3.connect(str(path), check_same_thread=False, timeout=30)
        self._db.execute("PRAGMA journal_mode=WAL")
        self._db.execute("CREATE TABLE IF NOT EXISTS kv (k TEXT PRIMARY KEY, v TEXT NOT NULL, t REAL NOT NULL)")
        self._db.commit()

    def get(self, key: str, max_age: float | None = None) -> Any:
        """Return the cached value, or None when missing or older than `max_age` seconds (default: the TTL)."""
        limit = max_age if max_age is not None else self.ttl
        with self._lock:
            row = self._db.execute("SELECT v, t FROM kv WHERE k = ?", (key,)).fetchone()
        if row is None or (limit and time.time() - row[1] > limit):
            return None
        return json.loads(row[0])

    def stored(self, key: str) -> float | None:
        """Return when `key` was stored (epoch seconds), or None when it is not in the cache."""
        with self._lock:
            row = self._db.execute("SELECT t FROM kv WHERE k = ?", (key,)).fetchone()
        return None if row is None else float(row[0])

    def get_many(self, keys: list[str], max_age: float | None = None) -> dict[str, Any]:
        """Return key -> value for the `keys` that are stored and not older than `max_age` seconds.

        Args:
            keys: Keys to look up (queried 500 at a time).
            max_age: Maximum age in seconds; None means the TTL, 0 means no limit.
        """
        out: dict[str, Any] = {}
        now = time.time()
        limit = max_age if max_age is not None else self.ttl
        with self._lock:
            for i in range(0, len(keys), 500):
                chunk = keys[i : i + 500]
                q = f"SELECT k, v, t FROM kv WHERE k IN ({','.join('?' * len(chunk))})"
                for k, v, t in self._db.execute(q, chunk):
                    if not limit or now - t <= limit:
                        out[k] = json.loads(v)
        return out

    def set(self, key: str, value: Any) -> None:
        """Store a JSON-serialisable `value` under `key` (replacing any older value)."""
        self.set_many({key: value})

    def set_many(self, items: dict[str, Any]) -> None:
        """Store several JSON-serialisable values in one transaction, all with the current time."""
        if not items:
            return
        now = time.time()
        rows = [(k, json.dumps(v), now) for k, v in items.items()]
        with self._lock:
            self._db.executemany("INSERT OR REPLACE INTO kv (k, v, t) VALUES (?, ?, ?)", rows)
            self._db.commit()

    def clear(self) -> None:
        """Delete every stored value."""
        with self._lock:
            self._db.execute("DELETE FROM kv")
            self._db.commit()

    def close(self) -> None:
        """Close the database connection. The cache cannot be used afterwards."""
        with self._lock:
            self._db.close()

    def __enter__(self) -> Cache:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


class NullCache:
    """Drop-in replacement when caching is disabled."""

    def get(self, key: str, max_age: float | None = None) -> Any:
        """Return None: nothing is ever cached."""
        return None

    def stored(self, key: str) -> float | None:
        """Return None: nothing is ever cached."""
        return None

    def get_many(self, keys: list[str], max_age: float | None = None) -> dict[str, Any]:
        """Return an empty dict: nothing is ever cached."""
        return {}

    def set(self, key: str, value: Any) -> None:
        """Discard the value."""

    def set_many(self, items: dict[str, Any]) -> None:
        """Discard the values."""

    def clear(self) -> None:
        """Do nothing: there is nothing to delete."""

    def close(self) -> None:
        """Do nothing: there is no connection."""
