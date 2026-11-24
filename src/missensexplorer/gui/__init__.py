"""MissensExplorer desktop app (Qt / PySide6).

Run with:  missensexplorer gui     (or)     missensexplorer-gui
"""

from __future__ import annotations

import contextlib
import sys
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from .window import MainWindow


def prepare_web_engine() -> None:
    """QtWebEngine must be loaded before the QApplication exists (3D structure viewer)."""
    from PySide6.QtCore import QCoreApplication, Qt

    QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
    import importlib

    with contextlib.suppress(ImportError):
        importlib.import_module("PySide6.QtWebEngineWidgets")


def main(argv: list[str] | None = None) -> int:
    """Run the desktop app until its window is closed.

    Args:
        argv: Command-line arguments for Qt, without the program name; None uses ``sys.argv``.

    Returns:
        The exit code: Qt's, or 1 when PySide6 is not installed.
    """
    try:
        from PySide6.QtCore import QLocale, Qt
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("The desktop app needs Qt for Python:  pip install 'missensexplorer[gui]'", file=sys.stderr)
        return 1

    from . import i18n
    from . import theme as th

    argv = list(sys.argv if argv is None else [sys.argv[0], *argv])
    prepare_web_engine()
    QApplication.setHighDpiScaleFactorRoundingPolicy(Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    app = cast(QApplication, QApplication.instance() or QApplication(argv))
    app.setApplicationName("MissensExplorer")
    app.setOrganizationName("MissensExplorer")
    app.setApplicationDisplayName("MissensExplorer")
    app.setStyle("Fusion")
    app.setPalette(th.palette())
    app.setFont(th.font(13))
    app.setStyleSheet(th.stylesheet())
    app.setWindowIcon(th.app_icon())

    from . import logs

    def unexpected(exc: BaseException) -> None:
        """Tell the user about an uncaught exception (already written to the log)."""
        from PySide6.QtWidgets import QMessageBox

        QMessageBox.warning(
            holder.get("win"), i18n.t("unexpected"), i18n.t("unexpected_body", path=str(logs.log_file())) + f"\n\n{exc}"
        )

    holder: dict[str, MainWindow] = {}
    logs.setup(unexpected)
    from .settings import app_settings

    settings = app_settings()
    i18n.set_language(str(settings.value("language", QLocale.system().name())))

    from .window import MainWindow

    def build(old: MainWindow | None = None) -> MainWindow:
        """Create and show a main window, carrying over the inputs and results of ``old``."""
        win = MainWindow()
        win.language_changed.connect(switch)
        if old is not None:
            win.transfer_state(old)
        win.show()
        win.raise_()
        holder["win"] = win
        return win

    def switch(lang: str) -> None:
        """Save ``lang`` as the preferred language and rebuild the window in it."""
        settings.setValue("language", lang)
        if lang == i18n.LANG:
            return
        old = holder["win"]
        i18n.set_language(lang)
        build(old)
        old._save_settings()
        old.close()
        old.deleteLater()

    build()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
