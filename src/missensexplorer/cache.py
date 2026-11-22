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


def default_cache_dir() -> Path:
    base = os.environ.get("MISSENSEXPLORER_CACHE")
    if base:
        return Path(base)
    xdg = os.environ.get("XDG_CACHE_HOME")
    return Path(xdg) / "missensexplorer" if xdg else Path.home() / ".cache" / "missensexplorer"


class Cache:
    def __init__(self, path: str | Path | None = None, ttl_days: float = 30):
        path = Path(path) if path else default_cache_dir() / "cache.sqlite"
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.ttl = ttl_days * 86400
        self._lock = threading.Lock()
        self._db = sqlite3.connect(str(path), check_same_thread=False, timeout=30)
        self._db.execute("PRAGMA journal_mode=WAL")
        self._db.execute(
            "CREATE TABLE IF NOT EXISTS kv (k TEXT PRIMARY KEY, v TEXT NOT NULL, t REAL NOT NULL)"
        )
        self._db.commit()

    def get(self, key: str):
        with self._lock:
            row = self._db.execute("SELECT v, t FROM kv WHERE k = ?", (key,)).fetchone()
        if row is None or (self.ttl and time.time() - row[1] > self.ttl):
            return None
        return json.loads(row[0])

    def get_many(self, keys: list[str]) -> dict:
        out = {}
        now = time.time()
        with self._lock:
            for i in range(0, len(keys), 500):
                chunk = keys[i : i + 500]
                q = f"SELECT k, v, t FROM kv WHERE k IN ({','.join('?' * len(chunk))})"
                for k, v, t in self._db.execute(q, chunk):
                    if not self.ttl or now - t <= self.ttl:
                        out[k] = json.loads(v)
        return out

    def set(self, key: str, value) -> None:
        self.set_many({key: value})

    def set_many(self, items: dict) -> None:
        if not items:
            return
        now = time.time()
        rows = [(k, json.dumps(v), now) for k, v in items.items()]
        with self._lock:
            self._db.executemany("INSERT OR REPLACE INTO kv (k, v, t) VALUES (?, ?, ?)", rows)
            self._db.commit()

    def clear(self) -> None:
        with self._lock:
            self._db.execute("DELETE FROM kv")
            self._db.commit()


class NullCache:
    """Drop-in replacement when caching is disabled."""

    def get(self, key):
        return None

    def get_many(self, keys):
        return {}

    def set(self, key, value):
        pass

    def set_many(self, items):
        pass

    def clear(self):
        pass
