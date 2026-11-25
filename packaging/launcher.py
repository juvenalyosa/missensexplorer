"""Entry point of the packaged desktop app (PyInstaller: macOS, Windows and Linux).

Set MISSENSEXPLORER_SMOKE=/path/shot.png to run a self-test of a packaged build: it predicts
BRAF V600E, renders the 3D structure, saves a screenshot and a log next to it (shot.png.log) and
exits. Exit codes: 0 all fine, 1 no result, 2 no 3D viewer, 3 the 3D viewer did not draw.
"""

import os
import sys

from missensexplorer.gui import main


def smoke(out: str) -> int:
    from missensexplorer.gui import prepare_web_engine

    def say(*parts: object) -> None:  # a windowed app may have no console: log to a file too
        line = " ".join(str(p) for p in parts)
        print(line, flush=True)
        with open(out + ".log", "a", encoding="utf-8") as fh:
            fh.write(line + "\n")

    prepare_web_engine()
    import tempfile

    from PySide6.QtCore import QTimer

    # self-test with default options, independent of (and without touching) the user's preferences
    os.environ["MISSENSEXPLORER_SETTINGS"] = os.path.join(tempfile.mkdtemp(prefix="mx-smoke-"), "settings.ini")
    from PySide6.QtWidgets import QApplication

    from missensexplorer.gui import i18n
    from missensexplorer.gui import theme as th

    app = QApplication(sys.argv[:1])
    app.setStyle("Fusion")
    app.setPalette(th.palette())
    app.setStyleSheet(th.stylesheet())
    i18n.set_language("es")
    from missensexplorer.gui.window import MainWindow

    w = MainWindow()
    w.resize(1500, 940)
    w.show()
    w.text.setPlainText("BRAF V600E")
    w.run_prediction()
    state = {"code": 1, "ticks": 0}

    def finish():
        w.grab().save(out)
        app.quit()

    def js(result):
        say("3Dmol:", result)
        state["code"] = 0 if str(result).startswith("object|1") else 3
        QTimer.singleShot(500, finish)

    def poll():
        state["ticks"] += 1
        if w.thread is None and w.res is not None and state["ticks"] > 40:
            say("rows:", w.proxy.rowCount(), "prediction:", w.res.table["prediction"].tolist())
            if w.web is None:
                state["code"] = 2
                return finish()
            w.web.page().runJavaScript("[typeof $3Dmol, document.querySelectorAll('canvas').length].join('|')", 0, js)
        elif state["ticks"] > 600:
            say("timeout")
            finish()
        else:
            QTimer.singleShot(250, poll)

    QTimer.singleShot(250, poll)
    app.exec()
    return state["code"]


if __name__ == "__main__":
    if sys.platform.startswith("linux"):
        # Ubuntu 24.04+ blocks the user namespaces Chromium's sandbox needs, and Qt's web engine
        # then fails to start. The web view only shows the app's own page, 3Dmol and AlphaFold models.
        os.environ.setdefault("QTWEBENGINE_DISABLE_SANDBOX", "1")
    if os.environ.get("MISSENSEXPLORER_SMOKE"):
        sys.exit(smoke(os.environ["MISSENSEXPLORER_SMOKE"]))
    sys.exit(main([]))
