"""Type aliases and protocols shared across the package."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, Protocol, TypeAlias

Row: TypeAlias = dict[str, Any]
"""One row of the results table: column name -> value (str, int, float or None)."""

JSON: TypeAlias = Any
"""Parsed JSON from a web API. Its shape is documented where it is read."""

Progress: TypeAlias = Callable[[str, int, int], None]
"""Progress callback: ``(stage name, items done, items total)``."""


class KeyValueCache(Protocol):
    """What the pipeline needs from a cache: :class:`~missensexplorer.cache.Cache` or ``NullCache``."""

    def get(self, key: str, max_age: float | None = None) -> Any:
        """Return the value stored under ``key``, or ``None`` if absent or older than ``max_age`` seconds."""
        ...

    def get_many(self, keys: list[str], max_age: float | None = None) -> dict[str, Any]:
        """Return the stored values of ``keys`` that exist and are fresh enough."""
        ...

    def set(self, key: str, value: Any) -> None:
        """Store a JSON-serialisable ``value`` under ``key``."""
        ...

    def set_many(self, items: dict[str, Any]) -> None:
        """Store several values at once."""
        ...

    def freeze(self) -> None:
        """Return stored values whatever their age (a rerun from the cache without the internet)."""
        ...
