"""Table model, filtering and cell painting for the results grid."""

from __future__ import annotations

import math
import numbers

import pandas as pd
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QPointF, QRectF, QSortFilterProxyModel, Qt
from PySide6.QtGui import QBrush, QColor, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QStyle, QStyledItemDelegate

from ..predict import prediction_status
from . import theme as th
from .i18n import t, tr_prediction, tr_text

STATUS_ROLE = Qt.ItemDataRole.UserRole + 1
SORT_ROLE = Qt.ItemDataRole.UserRole + 2


def isnum(v) -> bool:
    return isinstance(v, numbers.Real) and not isinstance(v, bool) and not (isinstance(v, float) and math.isnan(v))


def val(v, default=None):
    if v is None or (isinstance(v, float) and math.isnan(v)) or (isinstance(v, str) and not v.strip()):
        return default
    try:
        if pd.isna(v):
            return default
    except (TypeError, ValueError):
        pass
    return v


def consequence_text(c) -> str:
    c = val(c, "")
    return ", ".join(x.replace("_variant", "").replace("_", " ") for x in str(c).split(",") if x)


def clinvar_text(r) -> str:
    c = val(r.get("clinvar"), "")
    stars = r.get("clinvar_stars")
    if c and isnum(stars) and stars > 0:
        c = f"{c}  {'★' * int(stars)}"
    return c


# (key, header string key, kind)
COLUMNS = [
    ("input", "col_input", "text"),
    ("gene", "col_gene", "gene"),
    ("protein_change", "col_change", "change"),
    ("consequence", "col_consequence", "consequence"),
    ("prediction", "col_prediction", "status"),
    ("effect_score", "col_effect", "bar"),
    ("alphamissense", "AlphaMissense", "bar"),
    ("revel", "REVEL", "bar"),
    ("cadd_phred", "CADD", "num1"),
    ("plddt", "pLDDT", "num0"),
    ("clinvar", "col_clinvar", "clinvar"),
    ("gnomad_af", "col_gnomad", "sci"),
    ("zygosity", "col_zyg", "zyg"),
]


class ResultsModel(QAbstractTableModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.df = pd.DataFrame()
        self.rows: list[dict] = []
        self.index_labels: list = []
        self.zyg_from_carriers = False

    def set_frame(self, df: pd.DataFrame) -> None:
        self.beginResetModel()
        self.df = df
        self.rows = df.to_dict("records") if len(df) else []
        self.index_labels = list(df.index)
        self.zyg_from_carriers = bool(len(df)) and df["zygosity"].isna().all() and df["carriers"].notna().any()
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if orientation == Qt.Orientation.Horizontal and role == Qt.ItemDataRole.DisplayRole:
            key = COLUMNS[section][1]
            return t(key).upper() if key.startswith("col_") else key.upper()
        if orientation == Qt.Orientation.Horizontal and role == Qt.ItemDataRole.TextAlignmentRole:
            kind = COLUMNS[section][2]
            return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter) if kind in ("num1", "num0", "sci") \
                else int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        return None

    def raw(self, row: int, key: str):
        r = self.rows[row]
        if key == "zygosity" and self.zyg_from_carriers:
            return val(r.get("carriers"))
        return val(r.get(key))

    def column_empty(self, col: int) -> bool:
        key = COLUMNS[col][0]
        return all(self.raw(i, key) is None for i in range(len(self.rows)))

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        key, _, kind = COLUMNS[index.column()]
        r = self.rows[index.row()]
        v = self.raw(index.row(), key)
        if role == Qt.ItemDataRole.DisplayRole:
            if kind == "status":
                return tr_prediction(v)
            if kind == "consequence":
                return consequence_text(v)
            if kind == "clinvar":
                return clinvar_text(r)
            if kind == "zyg":
                return tr_text(v) if v else ""
            if v is None:
                return ""
            if kind == "bar":
                return f"{v:.2f}" if isnum(v) else ""
            if kind == "num1":
                return f"{v:.1f}" if isnum(v) else str(v)
            if kind == "num0":
                return f"{v:.0f}" if isnum(v) else str(v)
            if kind == "sci":
                return f"{v:.1e}" if isnum(v) else str(v)
            return str(v)
        if role == SORT_ROLE:
            if kind == "status":
                return th.STATUS_ORDER.index(prediction_status(v)) * 10 - (r.get("effect_score") or 0)
            if kind in ("bar", "num1", "num0", "sci"):
                return float(v) if isnum(v) else -1.0
            return str(v or "").lower()
        if role == STATUS_ROLE:
            return prediction_status(r.get("prediction"))
        if role == Qt.ItemDataRole.TextAlignmentRole:
            if kind in ("num1", "num0", "sci"):
                return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            return int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        if role == Qt.ItemDataRole.FontRole:
            if kind == "gene":
                return th.font(13, th.QFont.Weight.DemiBold)
            if kind in ("change", "text"):
                return th.font(12, mono=True)
        if role == Qt.ItemDataRole.ForegroundRole:
            if kind == "text":
                return QBrush(QColor(th.INK_2))
            if kind in ("consequence", "zyg", "sci"):
                return QBrush(QColor(th.INK_2))
        if role == Qt.ItemDataRole.ToolTipRole:
            if kind == "status":
                return tr_text(val(r.get("prediction_basis"), "")) or None
            if kind == "clinvar" and val(r.get("clinvar_review")):
                return r["clinvar_review"]
            if kind == "text":
                return str(v or "")
        return None


class FilterProxy(QSortFilterProxyModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.statuses: set[str] | None = None
        self.query = ""
        self.setSortRole(SORT_ROLE)

    def _refilter(self, change):
        if hasattr(self, "beginFilterChange"):  # Qt >= 6.10
            self.beginFilterChange()
            change()
            self.endFilterChange(QSortFilterProxyModel.Direction.Rows)
        else:
            change()
            self.invalidateFilter()

    def set_statuses(self, s):
        self._refilter(lambda: setattr(self, "statuses", set(s) if s is not None else None))

    def set_query(self, q: str):
        self._refilter(lambda: setattr(self, "query", q.strip().lower()))

    def filterAcceptsRow(self, row, parent):
        m: ResultsModel = self.sourceModel()
        r = m.rows[row]
        if self.statuses is not None and prediction_status(r.get("prediction")) not in self.statuses:
            return False
        if self.query:
            hay = " ".join(str(val(r.get(k), "")) for k in ("gene", "protein_change", "input", "hgvsp", "hgvsc",
                                                             "rsid", "clinvar", "sample"))
            return self.query in hay.lower()
        return True


# --------------------------------------------------------------------------- painting helpers

def shape_path(shape: str, c: QPointF, r: float) -> QPainterPath:
    p = QPainterPath()
    if shape == "triangle":
        p.moveTo(c.x(), c.y() - r * 1.1)
        p.lineTo(c.x() + r * 1.05, c.y() + r * 0.75)
        p.lineTo(c.x() - r * 1.05, c.y() + r * 0.75)
        p.closeSubpath()
    elif shape == "diamond":
        p.moveTo(c.x(), c.y() - r * 1.15)
        p.lineTo(c.x() + r * 1.0, c.y())
        p.lineTo(c.x(), c.y() + r * 1.15)
        p.lineTo(c.x() - r * 1.0, c.y())
        p.closeSubpath()
    elif shape == "square":
        p.addRoundedRect(QRectF(c.x() - r * 0.85, c.y() - r * 0.85, r * 1.7, r * 1.7), 1.5, 1.5)
    else:
        p.addEllipse(c, r, r)
    return p


def paint_status_mark(p: QPainter, status: str, c: QPointF, r: float = 5.0, ring: str | None = None) -> None:
    path = shape_path(th.STATUS_SHAPE[status], c, r)
    if ring:
        p.setPen(QPen(QColor(ring), 2))
    else:
        p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(QColor(th.STATUS_COLOR[status]))
    p.drawPath(path)


class StatusDelegate(QStyledItemDelegate):
    """Shape + colour + label: the status is never encoded by colour alone."""

    def paint(self, p: QPainter, opt, index):
        self.initStyleOption(opt, index)
        text = opt.text
        opt.text = ""
        opt.widget.style().drawControl(QStyle.ControlElement.CE_ItemViewItem, opt, p, opt.widget)
        status = index.data(STATUS_ROLE) or "Not assessed"
        p.save()
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        rect = opt.rect
        cx = rect.left() + 16
        paint_status_mark(p, status, QPointF(cx, rect.center().y() + 0.5), 5.2)
        p.setPen(QColor(th.INK))
        p.setFont(th.font(12.5))
        tr = rect.adjusted(28, 0, -6, 0)
        p.drawText(tr, Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                   p.fontMetrics().elidedText(text, Qt.TextElideMode.ElideRight, tr.width()))
        p.restore()


class BarDelegate(QStyledItemDelegate):
    """Thin 0–1 bar with the value as text (one hue: magnitude)."""

    def paint(self, p: QPainter, opt, index):
        self.initStyleOption(opt, index)
        text = opt.text
        opt.text = ""
        opt.widget.style().drawControl(QStyle.ControlElement.CE_ItemViewItem, opt, p, opt.widget)
        if not text:
            return
        try:
            v = max(0.0, min(1.0, float(text)))
        except ValueError:
            return
        p.save()
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        rect = opt.rect.adjusted(8, 0, -10, 0)
        p.setFont(th.font(12))
        p.setPen(QColor(th.INK))
        tw = 34
        p.drawText(QRectF(rect.left(), rect.top(), tw, rect.height()),
                   Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, text)
        track = QRectF(rect.left() + tw + 4, rect.center().y() - 2, max(10, rect.width() - tw - 4), 4)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor(th.MINERAL_100))
        p.drawRoundedRect(track, 2, 2)
        if v > 0:
            fill = QRectF(track.left(), track.top(), max(4, track.width() * v), track.height())
            p.setBrush(QColor(th.PETROL_500))
            p.drawRoundedRect(fill, 2, 2)
        p.restore()
