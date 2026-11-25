"""Write the app icon as PNG files (Linux, and the source of the Windows .ico), from the app's own icon.

    python packaging/make_icons.py   ->  packaging/icons/missensexplorer-<size>.png
"""

import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtCore import QSize  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from missensexplorer.gui import theme  # noqa: E402

OUT = Path(__file__).resolve().parent / "icons"


def main() -> int:
    app = QApplication(sys.argv[:1])  # noqa: F841 - needed to paint
    OUT.mkdir(exist_ok=True)
    icon = theme.app_icon()
    for size in (16, 32, 48, 64, 128, 256, 512):
        path = OUT / f"missensexplorer-{size}.png"
        icon.pixmap(QSize(size, size)).save(str(path))
        print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
