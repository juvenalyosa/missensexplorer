"""Error log for the desktop app (the packaged app has no console to print to)."""

from __future__ import annotations

import logging
import logging.handlers
import sys
import threading
from pathlib import Path


def log_dir() -> Path:
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Logs" / "MissensExplorer"
    from ..cache import default_cache_dir

    return default_cache_dir() / "logs"


def log_file() -> Path:
    return log_dir() / "missensexplorer.log"


def setup(on_unexpected=None) -> Path:
    """Log warnings and errors to a rotating file; route uncaught exceptions there too."""
    path = log_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.handlers.RotatingFileHandler(path, maxBytes=1_000_000, backupCount=3, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    handler.setLevel(logging.INFO)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    if not any(isinstance(h, logging.handlers.RotatingFileHandler) for h in root.handlers):
        root.addHandler(handler)
    for noisy in ("httpx", "httpcore"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    log = logging.getLogger("missensexplorer.gui")

    def hook(exc_type, exc, tb):
        log.error("uncaught exception", exc_info=(exc_type, exc, tb))
        if on_unexpected is not None:
            try:
                on_unexpected(exc)
            except Exception:  # never let the error dialog itself crash the app
                pass

    sys.excepthook = hook
    threading.excepthook = lambda a: log.error("uncaught exception in thread %s", a.thread,
                                               exc_info=(a.exc_type, a.exc_value, a.exc_traceback))
    log.info("MissensExplorer started")
    return path
