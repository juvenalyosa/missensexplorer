"""Entry point for the bundled macOS app (PyInstaller).

Set MISSENSEXPLORER_SMOKE=/path/shot.png to run a self-test of a packaged build: it predicts
BRAF V600E, renders the 3D structure, saves a screenshot and exits (code 0 on success).
"""

import os
import sys

from missensexplorer.gui import main


def smoke(out: str) -> int:
    from missensexplorer.gui import prepare_web_engine

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
        print("3Dmol:", result, flush=True)
        state["code"] = 0 if str(result).startswith("object|1") else 3
        QTimer.singleShot(500, finish)

    def poll():
        state["ticks"] += 1
        if w.thread is None and w.res is not None and state["ticks"] > 40:
            print("rows:", w.proxy.rowCount(), "prediction:", w.res.table["prediction"].tolist(), flush=True)
            if w.web is None:
                state["code"] = 2
                return finish()
            w.web.page().runJavaScript("[typeof $3Dmol, document.querySelectorAll('canvas').length].join('|')", 0, js)
        elif state["ticks"] > 600:
            print("timeout", flush=True)
            finish()
        else:
            QTimer.singleShot(250, poll)

    QTimer.singleShot(250, poll)
    app.exec()
    return state["code"]


if __name__ == "__main__":
    if os.environ.get("MISSENSEXPLORER_SMOKE"):
        sys.exit(smoke(os.environ["MISSENSEXPLORER_SMOKE"]))
    sys.exit(main([]))
