"""Shared test set-up."""

import pytest


@pytest.fixture(autouse=True)
def _isolated_cache(tmp_path_factory, monkeypatch):
    """Every test reads and writes its own cache, never the user's (whatever MISSENSEXPLORER_CACHE says)."""
    monkeypatch.setenv("MISSENSEXPLORER_CACHE", str(tmp_path_factory.mktemp("cache")))
