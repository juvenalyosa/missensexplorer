"""MissensExplorer desktop app (Qt / PySide6).

Run with:  missensexplorer gui     (or)     missensexplorer-gui
"""

from __future__ import annotations

import sys


def prepare_web_engine() -> None:
    """QtWebEngine must be loaded before the QApplication exists (3D structure viewer)."""
    from PySide6.QtCore import QCoreApplication, Qt

    QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
    import importlib

    try:
        importlib.import_module("PySide6.QtWebEngineWidgets")
    except ImportError:
        pass


def main(argv: list[str] | None = None) -> int:
    try:
        from PySide6.QtCore import QLocale, QSettings, Qt
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("The desktop app needs Qt for Python:  pip install 'missensexplorer[gui]'", file=sys.stderr)
        return 1

    from . import i18n
    from . import theme as th

    argv = list(sys.argv if argv is None else [sys.argv[0], *argv])
    prepare_web_engine()
    QApplication.setHighDpiScaleFactorRoundingPolicy(Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    app = QApplication.instance() or QApplication(argv)
    app.setApplicationName("MissensExplorer")
    app.setOrganizationName("MissensExplorer")
    app.setApplicationDisplayName("MissensExplorer")
    app.setStyle("Fusion")
    app.setPalette(th.palette())
    app.setFont(th.font(13))
    app.setStyleSheet(th.stylesheet())
    app.setWindowIcon(th.app_icon())

    from . import logs

    def unexpected(exc):
        from PySide6.QtWidgets import QMessageBox

        QMessageBox.warning(holder.get("win"), i18n.t("unexpected"),
                            i18n.t("unexpected_body", path=str(logs.log_file())) + f"\n\n{exc}")

    holder: dict = {}
    logs.setup(unexpected)
    settings = QSettings("MissensExplorer", "MissensExplorer")
    i18n.set_language(str(settings.value("language", QLocale.system().name())))

    from .window import MainWindow

    def build(old=None):
        win = MainWindow()
        win.language_changed.connect(switch)
        if old is not None:
            win.transfer_state(old)
        win.show()
        win.raise_()
        holder["win"] = win
        return win

    def switch(lang: str):
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
