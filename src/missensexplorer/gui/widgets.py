"""Reusable widgets: cards, switches, drop zone, predictor bars, lollipop chart, 3D viewer."""

from __future__ import annotations

import html
import json
from pathlib import Path

from PySide6.QtCore import (QEasingCurve, QPoint, QPointF, QPropertyAnimation, QRect, QRectF, QSize, Qt, QUrl,
                            Property, Signal)
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen
from PySide6.QtWidgets import (QAbstractButton, QFileDialog, QFrame, QHBoxLayout, QLabel, QLayout, QSizePolicy,
                               QToolButton, QToolTip, QVBoxLayout, QWidget)

from . import theme as th
from .i18n import t, tr_call, tr_prediction
from .models import isnum, paint_status_mark

INPUT_SUFFIXES = (".vcf", ".vcf.gz", ".gz", ".bgz", ".csv", ".tsv", ".tab", ".txt", ".xlsx", ".xls")


def label(text: str = "", obj: str | None = None, wrap: bool = False, selectable: bool = False) -> QLabel:
    w = QLabel(text)
    if obj:
        w.setObjectName(obj)
    w.setWordWrap(wrap)
    if selectable:
        w.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
    return w


class Card(QFrame):
    def __init__(self, parent=None, margins=(18, 16, 18, 16), spacing=10):
        super().__init__(parent)
        self.setObjectName("Card")
        self.lay = QVBoxLayout(self)
        self.lay.setContentsMargins(*margins)
        self.lay.setSpacing(spacing)


class Kpi(QFrame):
    def __init__(self, caption: str, status: str | None = None):
        super().__init__()
        self.setObjectName("Kpi")
        self.status = status
        lay = QVBoxLayout(self)
        lay.setContentsMargins(18, 14, 18, 14)
        lay.setSpacing(2)
        self.value = label("—", "KpiValue")
        self.caption = label(caption, "KpiLabel")
        lay.addWidget(self.value)
        lay.addWidget(self.caption)
        self.setMinimumWidth(150)

    def set(self, v) -> None:
        self.value.setText(f"{v:,}".replace(",", " ") if isinstance(v, int) else str(v))

    def paintEvent(self, e):
        super().paintEvent(e)
        if self.status:
            p = QPainter(self)
            p.setRenderHint(QPainter.RenderHint.Antialiasing)
            paint_status_mark(p, self.status, QPointF(self.width() - 22, 22), 5.5)
            p.end()


class Toggle(QAbstractButton):
    """iOS-like switch with a label, painted to match the sidebar."""

    def __init__(self, text: str, checked: bool = True, tooltip: str | None = None, dark: bool = True):
        super().__init__()
        self.setCheckable(True)
        self.setChecked(checked)
        self.setText(text)
        self.dark = dark
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        if tooltip:
            self.setToolTip(tooltip)
        self._pos = 1.0 if checked else 0.0
        self._anim = QPropertyAnimation(self, b"knob", self)
        self._anim.setDuration(140)
        self._anim.setEasingCurve(QEasingCurve.Type.OutCubic)
        self.toggled.connect(self._animate)
        self.setFont(th.font(12.5))

    def _animate(self, on):
        self._anim.stop()
        self._anim.setStartValue(self._pos)
        self._anim.setEndValue(1.0 if on else 0.0)
        self._anim.start()

    def get_knob(self):
        return self._pos

    def set_knob(self, v):
        self._pos = v
        self.update()

    knob = Property(float, get_knob, set_knob)

    def sizeHint(self):
        fm = QFontMetrics(self.font())
        return QSize(44 + fm.horizontalAdvance(self.text()) + 4, 26)

    def paintEvent(self, e):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        h = 18
        track = QRectF(1, (self.height() - h) / 2, 32, h)
        on = QColor(th.MINERAL_300 if self.dark else th.PETROL_700)
        off = QColor(th.PETROL_700 if self.dark else th.BORDER_STRONG)
        c = QColor(off.red() + (on.red() - off.red()) * self._pos, off.green() + (on.green() - off.green()) * self._pos,
                   off.blue() + (on.blue() - off.blue()) * self._pos)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(c)
        p.drawRoundedRect(track, h / 2, h / 2)
        r = h / 2 - 2.5
        cx = track.left() + h / 2 + (track.width() - h) * self._pos
        p.setBrush(QColor(th.PETROL_900 if (self.dark and self._pos > 0.5) else th.WHITE))
        p.drawEllipse(QPointF(cx, track.center().y()), r, r)
        p.setPen(QColor(th.MINERAL_100 if self.dark else th.INK))
        p.setFont(self.font())
        p.drawText(QRectF(44, 0, self.width() - 44, self.height()),
                   Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, self.text())
        p.end()


class FileChip(QFrame):
    removed = Signal(str)

    def __init__(self, path: str):
        super().__init__()
        self.setObjectName("FileChip")
        self.path = path
        lay = QHBoxLayout(self)
        lay.setContentsMargins(10, 6, 6, 6)
        lay.setSpacing(8)
        ic = QLabel()
        ic.setPixmap(th.icon("file", th.MINERAL_300, 16).pixmap(16, 16))
        lay.addWidget(ic)
        p = Path(path)
        size = p.stat().st_size if p.exists() else 0
        name = label(p.name)
        name.setToolTip(path)
        fm = QFontMetrics(name.font())
        name.setText(fm.elidedText(p.name, Qt.TextElideMode.ElideMiddle, 200))
        lay.addWidget(name, 1)
        sz = label(_human(size), "Hint")
        lay.addWidget(sz)
        x = QToolButton()
        x.setIcon(th.icon("x", th.MINERAL_300, 14))
        x.setCursor(Qt.CursorShape.PointingHandCursor)
        x.clicked.connect(lambda: self.removed.emit(self.path))
        lay.addWidget(x)


def _human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return ""


class DropZone(QFrame):
    files_dropped = Signal(list)

    def __init__(self):
        super().__init__()
        self.setObjectName("DropZone")
        self.setAcceptDrops(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(14, 12, 14, 12)
        lay.setSpacing(12)
        ic = QLabel()
        ic.setPixmap(th.icon("upload", th.MINERAL_300, 22).pixmap(22, 22))
        lay.addWidget(ic)
        self.text = label(t("drop"), "Hint", wrap=True)
        lay.addWidget(self.text, 1)
        self.start_dir = str(Path.home())

    def _set_hover(self, on: bool):
        self.setProperty("hover", "true" if on else "false")
        self.style().unpolish(self)
        self.style().polish(self)

    def mouseReleaseEvent(self, e):
        self.browse()

    def browse(self):
        files, _ = QFileDialog.getOpenFileNames(
            self, t("m_open"), self.start_dir,
            "Variant files (*.vcf *.vcf.gz *.gz *.bgz *.csv *.tsv *.tab *.txt *.xlsx);;All files (*)")
        if files:
            self.start_dir = str(Path(files[0]).parent)
            self.files_dropped.emit(files)

    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls():
            e.acceptProposedAction()
            self._set_hover(True)

    def dragLeaveEvent(self, e):
        self._set_hover(False)

    def dropEvent(self, e):
        self._set_hover(False)
        files = [u.toLocalFile() for u in e.mimeData().urls() if u.isLocalFile()]
        if files:
            self.files_dropped.emit(files)
            e.acceptProposedAction()


class FlowLayout(QLayout):
    """Wraps child widgets onto new lines (for chips)."""

    def __init__(self, parent=None, spacing=8):
        super().__init__(parent)
        self.items = []
        self._sp = spacing
        self.setContentsMargins(0, 0, 0, 0)

    def addItem(self, item):
        self.items.append(item)

    def count(self):
        return len(self.items)

    def itemAt(self, i):
        return self.items[i] if 0 <= i < len(self.items) else None

    def takeAt(self, i):
        return self.items.pop(i) if 0 <= i < len(self.items) else None

    def hasHeightForWidth(self):
        return True

    def heightForWidth(self, w):
        return self._layout(QRect(0, 0, w, 0), True)

    def setGeometry(self, rect):
        super().setGeometry(rect)
        self._layout(rect, False)

    def sizeHint(self):
        return self.minimumSize()

    def minimumSize(self):
        s = QSize()
        for it in self.items:
            s = s.expandedTo(it.minimumSize())
        return s

    def _layout(self, rect, test):
        x, y, line = rect.x(), rect.y(), 0
        for it in self.items:
            hint = it.sizeHint()
            if x + hint.width() > rect.right() + 1 and line > 0:
                x, y = rect.x(), y + line + self._sp
                line = 0
            if not test:
                it.setGeometry(QRect(QPoint(x, y), hint))
            x += hint.width() + self._sp
            line = max(line, hint.height())
        return y + line - rect.y()


# --------------------------------------------------------------------------- predictor bars

# name, column, call column, (min, max), thresholds [(value, label)], higher_is_damaging
PREDICTORS = [
    ("AlphaMissense", "alphamissense", "alphamissense_class", (0, 1), [(0.34, ""), (0.564, "")], True),
    ("REVEL", "revel", None, (0, 1), [(0.290, ""), (0.644, "")], True),
    ("CADD", "cadd_phred", None, (0, 50), [(22.7, ""), (25.3, "")], True),
    ("EVE", "eve", "eve_class", (0, 1), [(0.5, "")], True),
    ("SIFT", "sift_score", "sift", (0, 1), [(0.05, "")], False),
    ("PolyPhen-2", "polyphen_score", "polyphen", (0, 1), [(0.446, ""), (0.908, "")], True),
    ("SpliceAI", "spliceai_max", "spliceai_effect", (0, 1), [(0.2, ""), (0.5, "")], True),
]


def predictor_call(name: str, v, call) -> str | None:
    """'Damaging' | 'Benign / silent' | 'Uncertain' for the bar colour (same rules as predict._votes)."""
    c = (call or "").lower() if isinstance(call, str) else ""
    if name == "AlphaMissense":
        return {"likely_pathogenic": "Damaging", "likely_benign": "Benign / silent"}.get(c, "Uncertain" if c else None)
    if name == "REVEL" and isnum(v):
        return "Damaging" if v >= 0.644 else "Benign / silent" if v <= 0.290 else "Uncertain"
    if name == "CADD" and isnum(v):
        return "Damaging" if v >= 25.3 else "Benign / silent" if v <= 22.7 else "Uncertain"
    if name == "EVE" and c:
        return "Damaging" if c == "pathogenic" else "Benign / silent" if c == "benign" else "Uncertain"
    if name == "SIFT" and c:
        return "Damaging" if c.startswith("deleterious") else "Benign / silent"
    if name == "PolyPhen-2" and c:
        return "Damaging" if c == "probably_damaging" else "Benign / silent" if c == "benign" else "Uncertain"
    if name == "SpliceAI" and isnum(v):
        return "Damaging" if v >= 0.5 else "Uncertain" if v >= 0.2 else "Benign / silent"
    return None


class PredictorRow(QWidget):
    def __init__(self, name, v, call, rng, thresholds, higher_bad):
        super().__init__()
        self.name, self.v, self.call, self.rng = name, v, call, rng
        self.thresholds = thresholds
        self.higher_bad = higher_bad
        self.status = predictor_call(name, v, call)
        self.setMinimumHeight(30)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        tip = f"{name}: {v:.3f}" if isnum(v) else name
        if isinstance(call, str) and call:
            tip += f" · {tr_call(call)}"
        self.setToolTip(tip)

    def sizeHint(self):
        return QSize(380, 30)

    def paintEvent(self, e):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        h = self.height()
        name_w, val_w = 104, 46
        p.setFont(th.font(12.5, QFont.Weight.Medium))
        p.setPen(QColor(th.INK))
        p.drawText(QRectF(0, 0, name_w, h), Qt.AlignmentFlag.AlignVCenter, self.name)
        call_txt = ""
        if isinstance(self.call, str) and self.call:
            call_txt = tr_call(self.call)
        elif self.status:
            call_txt = t({"Damaging": "call_damaging", "Benign / silent": "call_benign"}.get(self.status, "call_neutral"))
        p.setFont(th.font(11.5))
        call_w = 168
        bar = QRectF(name_w, h / 2 - 3, max(40, self.width() - name_w - val_w - call_w - 16), 6)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor(th.MINERAL_100))
        p.drawRoundedRect(bar, 3, 3)
        lo, hi = self.rng
        if isnum(self.v):
            frac = max(0.0, min(1.0, (self.v - lo) / (hi - lo)))
            color = QColor(th.STATUS_COLOR.get(self.status or "Not assessed"))
            p.setBrush(color)
            p.drawRoundedRect(QRectF(bar.left(), bar.top(), max(6, bar.width() * frac), bar.height()), 3, 3)
        for tv, _ in self.thresholds:
            x = bar.left() + bar.width() * (tv - lo) / (hi - lo)
            p.setPen(QPen(QColor(th.INK_2), 1.2))
            p.drawLine(QPointF(x, bar.top() - 3), QPointF(x, bar.bottom() + 3))
        p.setPen(QColor(th.INK))
        p.setFont(th.font(12.5, mono=False))
        vt = (f"{self.v:.1f}" if self.name == "CADD" else f"{self.v:.2f}") if isnum(self.v) else "—"
        p.drawText(QRectF(bar.right() + 8, 0, val_w, h), Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignRight, vt)
        if call_txt:
            x0 = bar.right() + 8 + val_w + 8
            if self.status:
                paint_status_mark(p, self.status, QPointF(x0 + 6, h / 2 + 0.5), 4.2)
            p.setPen(QColor(th.INK_2))
            p.setFont(th.font(11.5))
            tr = QRectF(x0 + 16, 0, self.width() - x0 - 16, h)
            p.drawText(tr, Qt.AlignmentFlag.AlignVCenter,
                       QFontMetrics(p.font()).elidedText(call_txt, Qt.TextElideMode.ElideRight, int(tr.width())))
        p.end()


# --------------------------------------------------------------------------- lollipop chart

class LollipopChart(QWidget):
    """Variants along the protein: x = residue, y = effect score; colour + shape = prediction status."""

    picked = Signal(object)  # dataframe index label

    def __init__(self):
        super().__init__()
        self.setMouseTracking(True)
        self.points = []       # dicts: idx, pos, score, status, label, prediction
        self.domains = []      # dicts: name, start, end
        self.length = 0
        self.selected = None
        self.hover = None
        self.title = ""
        self.setMinimumHeight(260)

    def set_data(self, points, domains, length, selected, title=""):
        self.points, self.domains, self.length, self.selected, self.title = points, domains, length, selected, title
        self.hover = None
        self.update()

    def _geom(self):
        left, right, top, bottom = 54, 18, 44, 64
        return QRectF(left, top, max(50, self.width() - left - right), max(50, self.height() - top - bottom))

    def _xy(self, plot, pos, score):
        L = max(self.length or 1, 1)
        x = plot.left() + plot.width() * min(pos, L) / L
        y = plot.bottom() - plot.height() * (0.5 if score is None else score)
        return QPointF(x, y)

    def paintEvent(self, e):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.fillRect(self.rect(), QColor(th.WHITE))
        plot = self._geom()
        # legend (top-left), shape + label
        x = plot.left()
        p.setFont(th.font(11.5))
        present = [s for s in th.STATUS_ORDER if any(pt["status"] == s for pt in self.points)]
        for s in present:
            paint_status_mark(p, s, QPointF(x + 6, 16), 5)
            p.setPen(QColor(th.INK_2))
            txt = t(s)
            p.drawText(QPointF(x + 16, 20), txt)
            x += 16 + QFontMetrics(p.font()).horizontalAdvance(txt) + 18
        # grid + y axis (recessive)
        p.setFont(th.font(10.5))
        for yv in (0, 0.5, 1):
            y = plot.bottom() - plot.height() * yv
            p.setPen(QPen(QColor("#EDF2F0"), 1))
            p.drawLine(QPointF(plot.left(), y), QPointF(plot.right(), y))
            p.setPen(QColor(th.MUTED))
            p.drawText(QRectF(0, y - 8, plot.left() - 8, 16), Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                       f"{yv:g}")
        p.save()
        p.translate(14, plot.center().y())
        p.rotate(-90)
        p.drawText(QRectF(-80, -8, 160, 16), Qt.AlignmentFlag.AlignCenter, t("score_axis"))
        p.restore()
        # protein backbone + domains
        band = QRectF(plot.left(), plot.bottom() + 12, plot.width(), 16)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor(th.MINERAL_100))
        p.drawRoundedRect(QRectF(band.left(), band.center().y() - 3, band.width(), 6), 3, 3)
        L = max(self.length or 1, 1)
        p.setFont(th.font(10))
        for d in self.domains:
            x0 = band.left() + band.width() * d["start"] / L
            x1 = band.left() + band.width() * d["end"] / L
            r = QRectF(x0, band.top(), max(3, x1 - x0 - 2), band.height())
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(th.MINERAL_300))
            p.drawRoundedRect(r, 4, 4)
            name = d["name"]
            fm = QFontMetrics(p.font())
            if fm.horizontalAdvance(name) + 8 < r.width():
                p.setPen(QColor(th.PETROL_900))
                p.drawText(r, Qt.AlignmentFlag.AlignCenter, name)
        # x ticks
        p.setPen(QColor(th.MUTED))
        p.setFont(th.font(10.5))
        step = _nice_step(L)
        tick = 0
        while tick <= L:
            xx = band.left() + band.width() * tick / L
            p.drawText(QRectF(xx - 30, band.bottom() + 4, 60, 14), Qt.AlignmentFlag.AlignCenter, f"{tick}")
            tick += step
        p.drawText(QRectF(plot.left(), band.bottom() + 18, plot.width(), 14), Qt.AlignmentFlag.AlignCenter,
                   t("residue_axis") + (f" · {self.title}" if self.title else ""))
        # stems then heads (selected last)
        pts = sorted(self.points, key=lambda pt: pt["idx"] == self.selected)
        for pt in pts:
            c = self._xy(plot, pt["pos"], pt["score"])
            col = QColor(th.STATUS_COLOR[pt["status"]])
            col.setAlphaF(0.55)
            p.setPen(QPen(col, 2))
            p.drawLine(QPointF(c.x(), plot.bottom()), c)
        for pt in pts:
            c = self._xy(plot, pt["pos"], pt["score"])
            sel = pt["idx"] == self.selected
            r = 7.0 if sel else 5.5
            if sel:
                p.setPen(QPen(QColor(th.INK), 1.5))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawEllipse(c, r + 5, r + 5)
            paint_status_mark(p, pt["status"], c, r, ring=th.WHITE)
            if sel or pt["idx"] == self.hover:
                p.setPen(QColor(th.INK))
                p.setFont(th.font(11.5, QFont.Weight.DemiBold))
                fm = QFontMetrics(p.font())
                w = fm.horizontalAdvance(pt["label"]) + 4
                lx = min(max(c.x() - w / 2, plot.left()), plot.right() - w)
                p.drawText(QPointF(lx, c.y() - r - 9), pt["label"])
        p.end()

    def _nearest(self, pos):
        plot = self._geom()
        best, bd = None, 14.0
        for pt in self.points:
            c = self._xy(plot, pt["pos"], pt["score"])
            d = ((c.x() - pos.x()) ** 2 + (c.y() - pos.y()) ** 2) ** 0.5
            if d < bd:
                best, bd = pt, d
        return best

    def mouseMoveEvent(self, e):
        pt = self._nearest(e.position())
        idx = pt["idx"] if pt else None
        if idx != self.hover:
            self.hover = idx
            self.update()
        if pt:
            self.setCursor(Qt.CursorShape.PointingHandCursor)
            score = f"{pt['score']:.2f}" if pt["score"] is not None else "—"
            QToolTip.showText(e.globalPosition().toPoint(),
                              f"<b>{html.escape(pt['label'])}</b><br>{html.escape(tr_prediction(pt['prediction']))}"
                              f"<br>{t('col_effect')}: {score}<br>{html.escape(pt['input'])}", self)
        else:
            self.unsetCursor()
            QToolTip.hideText()

    def leaveEvent(self, e):
        self.hover = None
        self.update()

    def mouseReleaseEvent(self, e):
        pt = self._nearest(e.position())
        if pt:
            self.picked.emit(pt["idx"])


def _nice_step(L: int) -> int:
    for s in (10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 5000):
        if L / s <= 8:
            return s
    return 10000


# --------------------------------------------------------------------------- sequence

def sequence_html(wt: str, mut: str, pos: int, n_ref: int, n_alt: int, per_line: int = 50) -> str:
    """Wild type and mutant protein around the change (up to 4 lines of 50), change highlighted."""
    mono = th.mono_family()
    s0 = pos - 1
    lo = max(0, (s0 // per_line - 1) * per_line)

    def block(seq: str, hl_len: int, color: str) -> str:
        part = seq[lo: lo + per_line * 4]
        lines = []
        for i in range(0, len(part), per_line):
            cells = []
            for j, aa in enumerate(part[i: i + per_line]):
                k = lo + i + j
                ch = html.escape(aa)
                if s0 <= k < s0 + hl_len:
                    ch = f"<span style='background:{color};color:white;'>{ch}</span>"
                if j and j % 10 == 0:
                    cells.append(" ")
                cells.append(ch)
            lines.append(f"<span style='color:{th.MUTED}'>{lo + i + 1:>5}</span>  " + "".join(cells))
        return "<br>".join(lines)

    def head(txt: str) -> str:
        return (f"<div style='font-size:11px;letter-spacing:1px;color:{th.MUTED};margin:8px 0 4px'>{txt}</div>")

    css = f"font-family:{mono};font-size:12px;line-height:170%;color:{th.INK};white-space:pre"
    return (head(f"{t('wild_type').upper()} · {len(wt)} aa") + f"<div style='{css}'>{block(wt, n_ref, th.PETROL_600)}</div>"
            + head(f"{t('mutant').upper()} · {len(mut)} aa") + f"<div style='{css}'>{block(mut, max(n_alt, 1), th.ACCENT)}</div>")


# --------------------------------------------------------------------------- 3D structure

THREEDMOL_URL = "https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"
THREEDMOL_FILE = "3Dmol-2.4.2-min.js"


def ensure_3dmol(tmp_dir: Path) -> bool:
    """Download 3Dmol.js once (Python side) so the embedded browser never needs the network."""
    import httpx

    f = tmp_dir / THREEDMOL_FILE
    if f.exists() and f.stat().st_size > 100_000:
        return True
    try:
        r = httpx.get(THREEDMOL_URL, timeout=60, follow_redirects=True)
    except httpx.HTTPError:
        return False
    if r.status_code != 200 or len(r.content) < 100_000:
        return False
    tmp_dir.mkdir(parents=True, exist_ok=True)
    f.write_bytes(r.content)
    return True


def _js(value) -> str:
    """JSON for inline <script>: '</' must not appear, or the text could close the script tag."""
    return json.dumps(value).replace("</", "<\\/")


def structure_page(pdb: str, pos: int, label_text: str) -> str:
    return f"""<!doctype html><html><head><meta charset="utf-8">
<style>html,body{{margin:0;height:100%;background:#fff;overflow:hidden;font-family:Helvetica Neue,Helvetica,Arial}}
#mol{{position:absolute;inset:0}}
#legend{{position:absolute;left:12px;bottom:10px;font-size:11px;color:{th.INK_2};background:rgba(255,255,255,.9);
padding:6px 9px;border-radius:8px;border:1px solid {th.BORDER}}}
#legend i{{display:inline-block;width:10px;height:10px;border-radius:2px;margin:0 4px -1px 8px}}
#legend i:first-child{{margin-left:0}}</style>
<script src="{THREEDMOL_FILE}"></script></head>
<body><div id="mol"></div>
<div id="legend">pLDDT {''.join(f'<i style="background:{c}"></i>{lbl}' for _, c, lbl in th.PLDDT_COLORS)}
<i style="background:{th.ACCENT};border-radius:50%"></i>{html.escape(label_text)}</div>
<script>
  const v = $3Dmol.createViewer(document.getElementById("mol"), {{backgroundColor: "white"}});
  v.addModel({_js(pdb)}, "pdb");
  const plddt = a => a.b >= 90 ? "#0053D6" : a.b >= 70 ? "#65CBF3" : a.b >= 50 ? "#FFDB13" : "#FF7D45";
  v.setStyle({{}}, {{cartoon: {{colorfunc: plddt}}}});
  v.addStyle({{resi: {pos}}}, {{stick: {{color: "{th.ACCENT}", radius: 0.28}}}});
  v.addStyle({{resi: {pos}, atom: "CA"}}, {{sphere: {{color: "{th.ACCENT}", radius: 1.15}}}});
  v.addLabel({_js(label_text)}, {{backgroundColor: "{th.PETROL_900}", backgroundOpacity: 0.92,
      fontColor: "white", fontSize: 13, inFront: true, showBackground: true, borderRadius: 6}},
      {{resi: {pos}, atom: "CA"}});
  v.zoomTo({{within: {{distance: 18, sel: {{resi: {pos}}}}}}});
  v.render();
  v.spin(false);
</script></body></html>"""


def make_web_view():
    """QWebEngineView for local pages (3Dmol.js and the PDB are fetched by Python), or None."""
    import os

    from PySide6.QtGui import QGuiApplication

    # Chromium needs a real display: skip it for headless runs (tests, screenshots) or on request.
    if os.environ.get("MISSENSEXPLORER_NO_WEBENGINE") or QGuiApplication.platformName() in ("offscreen", "minimal"):
        return None
    try:
        from PySide6.QtWebEngineCore import QWebEngineSettings
        from PySide6.QtWebEngineWidgets import QWebEngineView
    except ImportError:
        return None
    view = QWebEngineView()
    s = view.settings()
    s.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls, False)
    s.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessFileUrls, True)
    view.page().setBackgroundColor(QColor(th.WHITE))
    return view


def load_page(view, page: str, tmp_dir: Path) -> None:
    tmp_dir.mkdir(parents=True, exist_ok=True)
    f = tmp_dir / "structure.html"
    f.write_text(page, encoding="utf-8")
    view.setUrl(QUrl.fromLocalFile(str(f)))

