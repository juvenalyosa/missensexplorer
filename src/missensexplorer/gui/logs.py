"""Error log for the desktop app (the packaged app has no console to print to)."""

from __future__ import annotations

import contextlib
import logging
import logging.handlers
import sys
import threading
from collections.abc import Callable
from pathlib import Path
from types import TracebackType


def log_dir() -> Path:
    """Return the folder of the log file: ~/Library/Logs/MissensExplorer on macOS, else ``<cache>/logs``."""
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Logs" / "MissensExplorer"
    from ..cache import default_cache_dir

    return default_cache_dir() / "logs"


def log_file() -> Path:
    """Return the path of the app's log file (it may not exist yet)."""
    return log_dir() / "missensexplorer.log"


def setup(on_unexpected: Callable[[BaseException], object] | None = None) -> Path:
    """Log warnings and errors to a rotating file; route uncaught exceptions there too.

    Args:
        on_unexpected: Called with the exception after an uncaught exception in the main thread is logged
            (the app shows an error dialog); errors raised by it are ignored.

    Returns:
        The path of the log file.
    """
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

    def hook(exc_type: type[BaseException], exc: BaseException, tb: TracebackType | None) -> None:
        log.error("uncaught exception", exc_info=(exc_type, exc, tb))
        if on_unexpected is not None:
            with contextlib.suppress(Exception):  # never let the error dialog itself crash the app
                on_unexpected(exc)

    def thread_hook(a: threading.ExceptHookArgs) -> None:
        log.error(
            "uncaught exception in thread %s",
            a.thread,
            # typeshed allows exc_value=None here, logging's stub does not; logging handles it at runtime.
            exc_info=(a.exc_type, a.exc_value, a.exc_traceback),  # type: ignore[arg-type]
        )

    sys.excepthook = hook
    threading.excepthook = thread_hook
    log.info("MissensExplorer started")
    return path
