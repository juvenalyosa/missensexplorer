"""Where the desktop app keeps its preferences."""

from __future__ import annotations

import os

from PySide6.QtCore import QSettings

ENV = "MISSENSEXPLORER_SETTINGS"


def app_settings() -> QSettings:
    """The user's preferences, or an INI file named by $MISSENSEXPLORER_SETTINGS (tests and the packaged
    self-test use one, so they neither read nor overwrite the user's choices)."""
    path = os.environ.get(ENV)
    if path:
        return QSettings(path, QSettings.Format.IniFormat)
    return QSettings("MissensExplorer", "MissensExplorer")
