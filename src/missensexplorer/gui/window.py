"""Main window of the MissensExplorer desktop app."""

from __future__ import annotations

import functools
import html
import logging
import traceback
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any, cast
from urllib.parse import quote

import httpx
import pandas as pd
from PySide6.QtCore import (
    QItemSelectionModel,
    QModelIndex,
    QObject,
    QPointF,
    QRunnable,
    Qt,
    QThread,
    QThreadPool,
    QTimer,
    QUrl,
    Signal,
)
from PySide6.QtGui import (
    QAction,
    QActionGroup,
    QCloseEvent,
    QDesktopServices,
    QDragEnterEvent,
    QDropEvent,
    QGuiApplication,
    QIcon,
    QKeySequence,
    QPainter,
    QPaintEvent,
    QPixmap,
    QResizeEvent,
    QStandardItem,
    QStandardItemModel,
)
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLayout,
    QLineEdit,
    QMainWindow,
    QMenu,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QStackedWidget,
    QTableView,
    QTabWidget,
    QTextBrowser,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from .. import __version__, acmg, model
from .._types import Row
from ..cache import Cache, default_cache_dir
from ..inputs import vcf_samples
from ..predict import PredictionResult, PredictOptions, predict, prediction_status
from ..report import to_fasta, write
from . import i18n
from . import theme as th
from .i18n import aa_class, aa_name, t, tr_prediction, tr_stage, tr_text
from .models import COLUMNS, BarDelegate, FilterProxy, ResultsModel, StatusDelegate, isnum, paint_status_mark, val
from .settings import app_settings
from .widgets import (
    PREDICTORS,
    Card,
    ContributionChart,
    Domain,
    DropZone,
    FileChip,
    FlowLayout,
    Kpi,
    LollipopChart,
    LollipopPoint,
    PredictorRow,
    Toggle,
    ensure_3dmol,
    label,
    load_page,
    make_web_view,
    sequence_html,
    structure_page,
)

log = logging.getLogger(__name__)

EXAMPLE = """BRAF V600E
TP53: R175H, R248Q, R273H, P72R
KRAS:p.G12D
BRCA1:c.5096G>A
chr17:7675088 C>T
rs113488022
CFTR p.Phe508del"""

QUICK_EXAMPLES = [
    "BRAF V600E",
    "TP53: R175H, R248Q, R273H",
    "KRAS G12D",
    "CFTR p.Phe508del",
    "BRCA1:c.5096G>A",
    "rs28934578",
]


# --------------------------------------------------------------------------- background work


class Cancelled(Exception):
    """Raised from the progress callback to stop a prediction the user cancelled."""


class PredictWorker(QObject):
    """Run :func:`~missensexplorer.predict.predict` in a worker thread and report back through signals.

    Signals: ``progress(stage, done, total)``, ``finished(PredictionResult)``, ``failed(message)`` and
    ``cancelled()``.

    Attributes:
        stop: Set to True (from the GUI thread) to cancel at the next progress report.
    """

    progress = Signal(str, int, int)
    finished = Signal(object)
    failed = Signal(str)
    cancelled = Signal()

    def __init__(self, text: str, files: list[str], options: PredictOptions) -> None:
        super().__init__()
        self.text, self.files, self.options = text, files, options
        self.stop = False

    def _cb(self, stage: str, done: int, total: int) -> None:
        """Forward the progress of :func:`predict`, or raise :class:`Cancelled` once ``stop`` is set."""
        if self.stop:
            raise Cancelled
        self.progress.emit(stage, done, total)

    def run(self) -> None:
        """Run the prediction (network calls); emit ``finished``, ``failed`` or ``cancelled``."""
        try:
            res = predict(text=self.text, files=self.files, options=self.options, progress=self._cb)
        except Cancelled:
            self.cancelled.emit()
        except Exception as e:  # report, never crash the UI
            log.error("prediction failed\n%s", traceback.format_exc())
            self.failed.emit(str(e) or e.__class__.__name__)
        else:
            self.finished.emit(res)


class _PdbSignals(QObject):
    """Signals of :class:`PdbLoader` (a QRunnable is not a QObject and cannot emit them itself)."""

    done = Signal(str, object)  # url, pdb text or None


class PdbLoader(QRunnable):
    """Fetch an AlphaFold model in the thread pool and emit ``signals.done(url, pdb text or None)``.

    The PDB text is cached; 3Dmol.js is downloaded first if needed (see :func:`ensure_3dmol`).
    """

    def __init__(self, url: str, signals: _PdbSignals, tmp_dir: Path) -> None:
        super().__init__()
        self.url, self.signals, self.tmp_dir = url, signals, tmp_dir

    def run(self) -> None:
        """Load the PDB file (cache, else network) and emit it; None when 3Dmol.js or the file is unavailable."""
        key = f"pdb|{self.url}"
        if not ensure_3dmol(self.tmp_dir):
            self.signals.done.emit(self.url, None)
            return
        try:
            with Cache() as cache:
                txt = cache.get(key)
                if txt is None:
                    r = httpx.get(self.url, timeout=60, follow_redirects=True)
                    txt = r.text if r.status_code == 200 else None
                    if txt:
                        cache.set(key, txt)
        except Exception:  # network problems: show the fallback message
            txt = None
        self.signals.done.emit(self.url, txt)


# --------------------------------------------------------------------------- small pieces


class StatusDot(QWidget):
    """The shape and colour of a prediction status, as a small square widget of ``size`` pixels."""

    def __init__(self, status: str, size: int = 16) -> None:
        super().__init__()
        self.status = status
        self.setFixedSize(size, size)

    def paintEvent(self, e: QPaintEvent) -> None:
        """Paint the status mark in the centre."""
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        paint_status_mark(p, self.status, QPointF(self.width() / 2, self.height() / 2 + 0.5), self.width() / 3)
        p.end()


class Toast(QLabel):
    """Short message shown at the bottom centre of its parent, hidden after a few seconds."""

    def __init__(self, parent: QWidget) -> None:
        super().__init__(parent)
        self.setStyleSheet(
            f"background:{th.PETROL_900};color:white;border-radius:10px;padding:10px 16px;font-size:12.5px;"
        )
        self.hide()
        self.timer = QTimer(self, singleShot=True)
        self.timer.timeout.connect(self.hide)

    def show_message(self, text: str, ms: int = 2600) -> None:
        """Show ``text`` for ``ms`` milliseconds."""
        self.setText(text)
        self.adjustSize()
        self.place()
        self.raise_()
        self.show()
        self.timer.start(ms)

    def place(self) -> None:
        """Move to the bottom centre of the parent (also called when the window is resized)."""
        par = self.parentWidget()
        assert par is not None  # always created with a parent
        self.move((par.width() - self.width()) // 2, par.height() - self.height() - 28)


class Collapsible(QWidget):
    """A title button that shows or hides a rich-text body (closed at first)."""

    def __init__(self, title: str, body: str) -> None:
        super().__init__()
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        self.btn = QToolButton()
        self.btn.setText(title)
        self.btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.btn.setIcon(th.icon("chevron-right", th.INK_2, 14))
        self.btn.setCheckable(True)
        self.btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.btn.setStyleSheet(
            f"QToolButton{{border:0;background:transparent;color:{th.PETROL_700};font-weight:500;"
            f"padding:2px 0}} QToolButton:hover{{color:{th.PETROL_950}}}"
        )
        self.body = label(body, "Body", wrap=True, selectable=True)
        self.body.setTextFormat(Qt.TextFormat.RichText)
        self.body.setContentsMargins(18, 0, 0, 4)
        self.body.hide()
        self.btn.toggled.connect(self._toggle)
        lay.addWidget(self.btn)
        lay.addWidget(self.body)

    def _toggle(self, on: bool) -> None:
        """Show or hide the body and turn the chevron."""
        self.body.setVisible(on)
        self.btn.setIcon(th.icon("chevron-down" if on else "chevron-right", th.INK_2, 14))


def section(title: str) -> QLabel:
    """Return a section heading of the detail pane (style ``H3``, with space above)."""
    w = label(title, "H3")
    w.setContentsMargins(0, 10, 0, 0)
    return w


def divider() -> QFrame:
    """Return a 1-pixel horizontal rule."""
    f = QFrame()
    f.setObjectName("Divider")
    f.setFixedHeight(1)
    return f


def clear_layout(lay: QLayout) -> None:
    """Remove and delete every widget of ``lay``, emptying nested layouts too."""
    while lay.count():
        it = lay.takeAt(0)
        assert it is not None  # count() > 0
        w = it.widget()
        if w is not None:
            w.setParent(None)
            w.deleteLater()
        elif (sub := it.layout()) is not None:
            clear_layout(sub)


def _rich(text: str) -> QLabel:
    """Return a wrapped, selectable rich-text label in the ``Body`` style."""
    lw = label(text, "Body", wrap=True, selectable=True)
    lw.setTextFormat(Qt.TextFormat.RichText)
    return lw


def _getter(r: Row) -> Callable[..., Any]:
    """Return ``g(key, default=None)``: the value of ``r[key]``, or ``default`` when it is missing (see ``val``)."""

    def g(key: str, default: Any = None) -> Any:
        return val(r.get(key), default)

    return g


def _code_points(text: object) -> int | None:
    """Return the points of the ACMG code that ``text`` starts with (``"PM5 (p.Val600Gly ...)"`` -> 2), or None."""
    return acmg.CODE_POINTS.get(str(text or "").split(" ")[0])


# UI string of the note under the ACMG evidence, by the thresholds behind PP3/BP4 (acmg.thresholds_used).
ACMG_NOTE_KEY = {"model": "acmg_note_model", "clingen": "acmg_note", None: "acmg_note_general"}


def _fact_box(value: str, caption: str) -> QFrame:
    """Return a tile of the physico-chemical grid: a large value over a small caption."""
    box = QFrame()
    box.setStyleSheet(f"QFrame{{background:{th.MINERAL_50};border-radius:8px;}}")
    bl = QVBoxLayout(box)
    bl.setContentsMargins(10, 7, 10, 7)
    bl.setSpacing(0)
    kv = label(value)
    kv.setStyleSheet(f"font-size:16px;color:{th.PETROL_900};background:transparent;")
    bl.addWidget(kv)
    kl = label(caption, "Muted")
    kl.setStyleSheet("background:transparent;font-size:11px;")
    bl.addWidget(kl)
    return box


# --------------------------------------------------------------------------- main window


class MainWindow(QMainWindow):
    """The app window: inputs in the sidebar, results table, detail pane and structure / lollipop / sequence tabs.

    Emits ``language_changed(code)`` when the user picks another language; the app then builds a new window
    and moves the state over with :meth:`transfer_state`.

    Attributes:
        settings: The user's preferences (see :func:`~missensexplorer.gui.settings.app_settings`).
        res: The result shown, or None before the first run.
        df: ``res.table`` plus a ``status`` column (see :func:`~missensexplorer.predict.prediction_status`).
        files: Input files, in the order they were added.
        thread: The thread of the running prediction, or None when idle.
        worker: The running prediction, or None when idle.
        current_idx: Index label (in ``df``) of the variant shown in the detail pane, or None.
        pdb_cache: AlphaFold PDB text by URL (None: the download failed).
        pending_structure: (PDB URL, residue, label) of the structure to show once its PDB file is loaded.
        tmp_dir: Folder of the 3D viewer page and of 3Dmol.js.
    """

    language_changed = Signal(str)

    def __init__(self) -> None:
        super().__init__()
        self.settings = app_settings()
        self.res: PredictionResult | None = None
        self.df = pd.DataFrame()
        self.files: list[str] = []
        # Shadows QObject.thread(); kept: the packaged app's self-test (packaging/macos/launcher.py) reads it.
        self.thread: QThread | None = None  # type: ignore[assignment]
        self.worker: PredictWorker | None = None
        self.current_idx: int | None = None
        self.pdb_cache: dict[str, str | None] = {}
        self.pending_structure: tuple[str, int, str] | None = None
        self.pdb_signals = _PdbSignals()
        self.pdb_signals.done.connect(self._pdb_ready)
        self.tmp_dir = default_cache_dir() / "gui"
        self._sized = False  # splitter sizes set (on the first result)
        self._af_acc: str | None = None  # UniProt accession of the "Open in AlphaFold DB" button
        self._fasta: str | None = None  # FASTA of the variant shown in the sequence tab

        self.setWindowTitle("MissensExplorer")
        self.setWindowIcon(th.app_icon())
        self.setAcceptDrops(True)
        self.resize(1500, 940)
        self.setMinimumSize(1180, 740)

        central = QWidget()
        central.setObjectName("Main")
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(self._build_sidebar())
        main = QWidget()
        mv = QVBoxLayout(main)
        mv.setContentsMargins(0, 0, 0, 0)
        mv.setSpacing(0)
        mv.addWidget(self._build_header())
        self.stack = QStackedWidget()
        self.stack.addWidget(self._build_welcome())
        self.stack.addWidget(self._build_results())
        mv.addWidget(self.stack, 1)
        root.addWidget(main, 1)
        self.setCentralWidget(central)
        self.toast = Toast(central)
        self._build_menus()
        self._restore_settings()
        self._show_header_tools(False)

    # ------------------------------------------------------------------ sidebar
    def _build_sidebar(self) -> QWidget:
        """Build the left sidebar: brand, variant and file inputs, options, and the predict button."""
        side = QFrame()
        side.setObjectName("Sidebar")
        side.setFixedWidth(372)
        v = QVBoxLayout(side)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)
        v.addWidget(self._build_brand())
        v.addWidget(divider())

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        body = QWidget()
        body.setObjectName("SidebarBody")
        b = QVBoxLayout(body)
        b.setContentsMargins(24, 18, 24, 18)
        b.setSpacing(10)
        self._add_variant_inputs(b)
        b.addSpacing(8)
        self._add_file_inputs(b)
        b.addSpacing(8)
        self._add_options(b)
        b.addStretch(1)
        scroll.setWidget(body)
        v.addWidget(scroll, 1)

        foot = self._build_sidebar_footer()
        v.addWidget(divider())
        v.addWidget(foot)
        return side

    def _build_brand(self) -> QWidget:
        """Build the logo, app name and tagline at the top of the sidebar."""
        brand = QWidget()
        bl = QHBoxLayout(brand)
        bl.setContentsMargins(24, 22, 24, 18)
        bl.setSpacing(12)
        logo = QLabel()
        logo.setPixmap(th.app_icon().pixmap(38, 38))
        bl.addWidget(logo)
        bt = QVBoxLayout()
        bt.setSpacing(1)
        bt.addWidget(label("MissensExplorer", "Brand"))
        bt.addWidget(label(t("brand_sub"), "BrandSub"))
        bl.addLayout(bt, 1)
        return brand

    def _add_variant_inputs(self, b: QVBoxLayout) -> None:
        """Add the variants text box, its hint, the default gene, and the example and formats buttons to ``b``."""
        b.addWidget(label(t("variants"), "Section"))
        self.text = QPlainTextEdit()
        self.text.setPlaceholderText(EXAMPLE)
        self.text.setMinimumHeight(118)
        self.text.setMaximumHeight(170)
        self.text.setTabChangesFocus(True)
        b.addWidget(self.text)
        b.addWidget(label(t("variants_hint"), "Hint", wrap=True))
        row = QHBoxLayout()
        row.setSpacing(8)
        self.gene = QLineEdit(placeholderText=t("default_gene"))
        self.gene.setToolTip(t("default_gene_tip"))
        row.addWidget(self.gene, 1)
        ex = QPushButton(t("example"))
        ex.setObjectName("Ghost")
        ex.setCursor(Qt.CursorShape.PointingHandCursor)
        ex.clicked.connect(lambda: self.text.setPlainText(EXAMPLE))
        row.addWidget(ex)
        fm = QPushButton()
        fm.setObjectName("Ghost")
        fm.setIcon(th.icon("help", th.MINERAL_200, 16))
        fm.setToolTip(t("formats"))
        fm.setCursor(Qt.CursorShape.PointingHandCursor)
        fm.clicked.connect(self.show_formats)
        row.addWidget(fm)
        b.addLayout(row)

    def _add_file_inputs(self, b: QVBoxLayout) -> None:
        """Add the file drop zone, the file list, the gene filter and the VCF sample / PASS options to ``b``."""
        b.addWidget(label(t("files"), "Section"))
        self.drop = DropZone()
        self.drop.files_dropped.connect(self.add_files)
        b.addWidget(self.drop)
        self.file_list = QVBoxLayout()
        self.file_list.setSpacing(6)
        b.addLayout(self.file_list)
        self.genes = QLineEdit(placeholderText=t("only_genes"))
        self.genes.setToolTip(t("only_genes_tip"))
        b.addWidget(self.genes)
        self.sample = QComboBox()
        self.sample.addItem(t("all_carriers"))
        self.sample.setEnabled(False)
        self.sample.setToolTip(t("sample"))
        self.sample.setMinimumWidth(10)
        b.addWidget(self.sample)
        self.pass_only = Toggle(t("pass_only"), False)
        b.addWidget(self.pass_only)
        self.sample.hide()
        self.pass_only.hide()

    def _add_options(self, b: QVBoxLayout) -> None:
        """Add the collapsible options (button, one-line summary and panel) to ``b``."""
        self.opt_btn = QToolButton()
        self.opt_btn.setCheckable(True)
        self.opt_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.opt_btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.opt_btn.setText(t("options"))
        self.opt_btn.setStyleSheet(
            f"QToolButton{{border:0;background:transparent;color:{th.MINERAL_400};font-size:10.5px;"
            f"font-weight:600;letter-spacing:1.4px;padding:0}} QToolButton:hover{{color:{th.WHITE}}}"
        )
        b.addWidget(self.opt_btn, 0, Qt.AlignmentFlag.AlignLeft)
        self.opt_summary = label("", "Hint")
        b.addWidget(self.opt_summary)
        self.opt_panel = self._build_options_panel()
        b.addWidget(self.opt_panel)
        self.opt_btn.toggled.connect(self._toggle_options)
        self._engine_changed()
        self._toggle_options(False)

    def _build_options_panel(self) -> QWidget:
        """Build the options: genome build, VEP engine, data sources, frequency rule, ACMG source, workers."""
        panel = QWidget()
        op = QVBoxLayout(panel)
        op.setContentsMargins(0, 4, 0, 0)
        op.setSpacing(10)
        grow = QHBoxLayout()
        grow.addWidget(label(t("build")))
        self.build = QComboBox()
        self.build.addItem(t("build_auto"), "auto")
        self.build.addItem("GRCh38 / hg38", "GRCh38")
        self.build.addItem("GRCh37 / hg19", "GRCh37")
        self.build.setMinimumWidth(10)
        grow.addWidget(self.build, 1)
        op.addLayout(grow)
        erow = QHBoxLayout()
        erow.addWidget(label(t("engine")))
        self.engine = QComboBox()
        self.engine.addItem(t("engine_rest"), "rest")
        self.engine.addItem(t("engine_local"), "local")
        self.engine.setToolTip(t("engine_tip"))
        self.engine.setMinimumWidth(10)
        erow.addWidget(self.engine, 1)
        op.addLayout(erow)
        self.vep_cmd = QLineEdit(placeholderText=t("vep_cmd"))
        self.vep_cmd.setToolTip(t("engine_tip"))
        self.vep_args = QLineEdit(placeholderText=t("vep_args"))
        self.offline = Toggle(t("offline"), False, t("offline_tip"))
        for w in (self.vep_cmd, self.vep_args, self.offline):
            op.addWidget(w)
        self.engine.currentIndexChanged.connect(self._engine_changed)
        self.engine.currentIndexChanged.connect(self._options_summary)
        self.opt_pred = Toggle(t("opt_predictors"), True, t("opt_predictors_tip"))
        self.opt_ctx = Toggle(t("opt_context"), True)
        self.opt_cv = Toggle(t("opt_clinvar"), True)
        self.opt_nc = Toggle(t("opt_noncoding"), False)
        self.opt_onc = Toggle(t("opt_oncogenicity"), False, t("opt_oncogenicity_tip"))
        for tg in (self.opt_pred, self.opt_ctx, self.opt_cv, self.opt_nc, self.opt_onc):
            op.addWidget(tg)
            tg.toggled.connect(self._options_summary)
        self.build.currentIndexChanged.connect(self._options_summary)
        frow = QHBoxLayout()
        frow.addWidget(label(t("freq_rule")))
        self.freq_rule = QComboBox()
        self.freq_rule.addItem(t("freq_ba1"), 0.05)
        self.freq_rule.addItem(t("freq_bs1"), 0.01)
        self.freq_rule.addItem(t("freq_off"), 0.0)
        self.freq_rule.setToolTip(t("freq_tip"))
        self.freq_rule.setMinimumWidth(10)
        frow.addWidget(self.freq_rule, 1)
        op.addLayout(frow)
        arow = QHBoxLayout()
        arow.addWidget(label(t("acmg_src")))
        self.acmg_src = QComboBox()
        self.acmg_src.addItem(t("acmg_src_model"), "model")
        self.acmg_src.addItem(t("acmg_src_revel"), "revel")
        self.acmg_src.setToolTip(t("acmg_src_tip"))
        self.acmg_src.setMinimumWidth(10)
        arow.addWidget(self.acmg_src, 1)
        op.addLayout(arow)
        wrow = QHBoxLayout()
        wrow.addWidget(label(t("workers")))
        wrow.addStretch(1)
        self.workers = QSpinBox(minimum=1, maximum=8, value=4)
        self.workers.setFixedWidth(70)
        wrow.addWidget(self.workers)
        op.addLayout(wrow)
        return panel

    def _build_sidebar_footer(self) -> QWidget:
        """Build the predict button, the progress bar with its stage text and cancel link, and the disclaimer."""
        foot = QWidget()
        f = QVBoxLayout(foot)
        f.setContentsMargins(24, 14, 24, 20)
        f.setSpacing(8)
        self.run_btn = QPushButton(f"  {t('predict')}")
        self.run_btn.setObjectName("Primary")
        self.run_btn.setIcon(th.icon("play", th.PETROL_900, 16))
        self.run_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.run_btn.setToolTip("⌘↩" if QGuiApplication.platformName() == "cocoa" else "Ctrl+Enter")
        self.run_btn.clicked.connect(self.run_prediction)
        f.addWidget(self.run_btn)
        self.progress = QProgressBar(textVisible=False, maximum=1000)
        self.progress.hide()
        f.addWidget(self.progress)
        prow = QHBoxLayout()
        self.stage = label("", "Stage")
        self.stage.hide()
        prow.addWidget(self.stage, 1)
        self.cancel_btn = QPushButton(t("cancel"))
        self.cancel_btn.setObjectName("Link")
        self.cancel_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.cancel_btn.clicked.connect(self.cancel_prediction)
        self.cancel_btn.hide()
        prow.addWidget(self.cancel_btn)
        f.addLayout(prow)
        f.addWidget(label(t("research_only"), "Hint", wrap=True))
        return foot

    def _toggle_options(self, on: bool) -> None:
        """Show the options panel (``on``) or its one-line summary."""
        self.opt_panel.setVisible(on)
        self.opt_summary.setVisible(not on)
        self.opt_btn.setIcon(th.icon("chevron-down" if on else "chevron-right", th.MINERAL_400, 14))
        self._options_summary()

    def _options_summary(self, *_: object) -> None:
        """Update the one-line options summary (genome build, number of data sources, local VEP)."""
        n = sum(w.isChecked() for w in (self.opt_pred, self.opt_ctx, self.opt_cv, self.opt_nc, self.opt_onc))
        txt = t("options_summary", build=self.build.currentText().split(" /")[0], n=n)
        if self.engine.currentData() == "local":
            txt += f" · {t('engine_local')}"
        self.opt_summary.setText(txt)

    def _engine_changed(self, *_: object) -> None:
        """Show the local-VEP fields only for the local engine (offline mode needs it)."""
        local = self.engine.currentData() == "local"
        for w in (self.vep_cmd, self.vep_args, self.offline):
            w.setVisible(local)
        if not local:
            self.offline.setChecked(False)

    # ------------------------------------------------------------------ header
    def _build_header(self) -> QWidget:
        """Build the header: title and subtitle, errors button, search box and export menu."""
        h = QFrame()
        h.setObjectName("Header")
        h.setFixedHeight(72)
        lay = QHBoxLayout(h)
        lay.setContentsMargins(28, 10, 24, 10)
        lay.setSpacing(10)
        tv = QVBoxLayout()
        tv.setSpacing(0)
        self.title = label(t("results"), "Title")
        self.subtitle = label("", "Subtitle")
        tv.addWidget(self.title)
        tv.addWidget(self.subtitle)
        lay.addLayout(tv)
        lay.addStretch(1)
        self.err_btn = QPushButton()
        self.err_btn.setIcon(th.icon("alert", th.STATUS_COLOR["Uncertain"], 16))
        self.err_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.err_btn.clicked.connect(self.show_errors)
        lay.addWidget(self.err_btn)
        self.search = QLineEdit(placeholderText=t("search"))
        self.search.setObjectName("Search")
        self.search.addAction(th.icon("search", th.MUTED, 16), QLineEdit.ActionPosition.LeadingPosition)
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self._filter_changed)
        lay.addWidget(self.search)
        self.export_btn = QPushButton(f" {t('export')}")
        self.export_btn.setObjectName("Accent")
        self.export_btn.setIcon(th.icon("download", th.WHITE, 16))
        self.export_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        menu = QMenu(self.export_btn)
        for fmt, key in (
            ("xlsx", "m_export_xlsx"),
            ("csv", "m_export_csv"),
            ("json", "m_export_json"),
            ("fasta", "m_export_fasta"),
        ):
            menu.addAction(t(key), functools.partial(self.export, fmt))
        self.export_btn.setMenu(menu)
        lay.addWidget(self.export_btn)
        return h

    def _show_header_tools(self, on: bool) -> None:
        """Show the search box, export and errors buttons (results page), or the app name (welcome page)."""
        for w in (self.search, self.export_btn):
            w.setVisible(on)
        self.err_btn.setVisible(on and self.res is not None and len(self.res.errors) > 0)
        if not on:
            self.title.setText("MissensExplorer")
            self.subtitle.setText(t("brand_sub"))

    # ------------------------------------------------------------------ welcome
    def _build_welcome(self) -> QWidget:
        """Build the welcome page: introduction, three steps and one-click examples."""
        page = QScrollArea()
        page.setWidgetResizable(True)
        inner = QWidget()
        outer = QHBoxLayout(inner)
        outer.addStretch(1)
        col = QWidget()
        col.setMaximumWidth(940)
        c = QVBoxLayout(col)
        c.setContentsMargins(32, 64, 32, 40)
        c.setSpacing(18)
        c.addWidget(label(t("welcome_title"), "H1", wrap=True))
        c.addWidget(label(t("welcome_lead"), "Lead", wrap=True))
        c.addSpacing(10)
        steps = QHBoxLayout()
        steps.setSpacing(14)
        for tk, bk, ic in (
            ("step1_t", "step1_b", "upload"),
            ("step2_t", "step2_b", "sliders"),
            ("step3_t", "step3_b", "dna"),
        ):
            card = Card(margins=(18, 18, 18, 18), spacing=8)
            icl = QLabel()
            icl.setPixmap(th.icon(ic, th.PETROL_700, 24).pixmap(24, 24))
            card.lay.addWidget(icl)
            card.lay.addWidget(label(t(tk), "H2"))
            card.lay.addWidget(label(t(bk), "Body", wrap=True))
            card.lay.addStretch(1)
            steps.addWidget(card, 1)
        c.addLayout(steps)
        c.addSpacing(14)
        c.addWidget(label(t("try"), "H3"))
        chips = QWidget()
        flow = FlowLayout(chips, spacing=8)
        for ex in QUICK_EXAMPLES:
            b = QPushButton(ex)
            b.setObjectName("ExampleChip")
            b.setCursor(Qt.CursorShape.PointingHandCursor)
            b.clicked.connect(lambda _=False, e=ex: self._quick_run(e))
            flow.addWidget(b)
        vcf = Path(__file__).resolve().parents[3] / "examples" / "example.vcf"
        if vcf.exists():

            def run_example_vcf() -> None:
                self.add_files([str(vcf)])
                self.run_prediction()

            b = QPushButton("example.vcf")
            b.setObjectName("ExampleChip")
            b.setIcon(th.icon("file", th.PETROL_700, 14))
            b.setCursor(Qt.CursorShape.PointingHandCursor)
            b.clicked.connect(run_example_vcf)
            flow.addWidget(b)
        c.addWidget(chips)
        c.addStretch(1)
        outer.addWidget(col, 10)
        outer.addStretch(1)
        page.setWidget(inner)
        return page

    def _quick_run(self, text: str) -> None:
        """Predict ``text`` (an example chip of the welcome page)."""
        self.text.setPlainText(text)
        self.run_prediction()

    # ------------------------------------------------------------------ results
    def _build_results(self) -> QWidget:
        """Build the results page: KPIs, warnings, status filters, the table, the detail pane and the tabs."""
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(24, 20, 24, 20)
        v.setSpacing(14)
        v.addLayout(self._build_kpis())

        self.banners = QVBoxLayout()
        self.banners.setSpacing(8)
        v.addLayout(self.banners)
        v.addLayout(self._build_filter_row())

        self.vsplit = QSplitter(Qt.Orientation.Vertical)
        self.vsplit.setChildrenCollapsible(False)
        self.vsplit.addWidget(self._build_table_card())

        self.hsplit = QSplitter(Qt.Orientation.Horizontal)
        self.hsplit.setChildrenCollapsible(False)
        self.hsplit.addWidget(self._build_detail_card())
        self.hsplit.addWidget(self._build_tabs_card())
        self.hsplit.setStretchFactor(0, 11)
        self.hsplit.setStretchFactor(1, 10)
        self.vsplit.addWidget(self.hsplit)
        self.vsplit.setStretchFactor(0, 4)
        self.vsplit.setStretchFactor(1, 6)
        v.addWidget(self.vsplit, 1)
        return page

    def _build_kpis(self) -> QHBoxLayout:
        """Build the row of key figures: inputs, results, likely damaging, pathogenic in ClinVar."""
        kpis = QHBoxLayout()
        kpis.setSpacing(14)
        self.k_inputs = Kpi(t("kpi_inputs"))
        self.k_rows = Kpi(t("kpi_rows"))
        self.k_dmg = Kpi(t("kpi_damaging"), "Damaging")
        self.k_cv = Kpi(t("kpi_clinvar"))
        for k in (self.k_inputs, self.k_rows, self.k_dmg, self.k_cv):
            kpis.addWidget(k, 1)
        return kpis

    def _build_filter_row(self) -> QHBoxLayout:
        """Build the row of status filter buttons (filled per result) and the skipped-variants note."""
        frow = QHBoxLayout()
        frow.setSpacing(8)
        self.segments: dict[str, QPushButton] = {}
        self.seg_box = QHBoxLayout()
        self.seg_box.setSpacing(8)
        frow.addLayout(self.seg_box)
        frow.addStretch(1)
        self.skipped = label("", "Muted")
        frow.addWidget(self.skipped)
        return frow

    def _build_table_card(self) -> Card:
        """Build the results table (model, filter proxy, delegates, column widths) and its empty-state text."""
        tcard = Card(margins=(2, 4, 2, 2), spacing=0)
        self.model = ResultsModel(self)
        self.proxy = FilterProxy(self)
        self.proxy.setSourceModel(self.model)
        self.table = QTableView()
        self.table.setModel(self.proxy)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.verticalHeader().hide()
        self.table.verticalHeader().setDefaultSectionSize(36)
        self.table.setShowGrid(False)
        self.table.setWordWrap(False)
        self.table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        hh = self.table.horizontalHeader()
        hh.setHighlightSections(False)
        hh.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        hh.setStretchLastSection(True)
        hh.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.table.setItemDelegateForColumn(4, StatusDelegate(self.table))
        for c in (5, 6, 7):
            self.table.setItemDelegateForColumn(c, BarDelegate(self.table))
        for c, w in enumerate([170, 68, 90, 130, 240, 112, 112, 112, 60, 60, 220, 84, 110]):
            self.table.setColumnWidth(c, w)
        self.table.selectionModel().currentRowChanged.connect(self._row_changed)
        tcard.lay.addWidget(self.table)
        self.empty_rows = label(t("no_rows"), "Muted", wrap=True)
        self.empty_rows.setContentsMargins(18, 18, 18, 18)
        self.empty_rows.hide()
        tcard.lay.addWidget(self.empty_rows)
        return tcard

    def _build_detail_card(self) -> Card:
        """Build the scrollable detail pane; :meth:`show_detail` fills its layout ``dl``."""
        dcard = Card(margins=(0, 0, 0, 0), spacing=0)
        dscroll = QScrollArea()
        dscroll.setWidgetResizable(True)
        dscroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.detail = QWidget()
        self.detail.setStyleSheet("background: white;")
        self.dl = QVBoxLayout(self.detail)
        self.dl.setContentsMargins(22, 18, 22, 22)
        self.dl.setSpacing(8)
        dscroll.setWidget(self.detail)
        dcard.lay.addWidget(dscroll)
        return dcard

    def _build_tabs_card(self) -> Card:
        """Build the right-hand tabs: 3D structure, variants along the protein, mutant protein sequence."""
        rcard = Card(margins=(0, 4, 0, 0), spacing=0)
        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.addTab(self._build_structure_tab(), t("tab_structure"))
        self.tabs.addTab(self._build_lollipop_tab(), t("tab_lollipop"))
        self.tabs.addTab(self._build_sequence_tab(), t("tab_sequence"))
        rcard.lay.addWidget(self.tabs)
        return rcard

    def _build_structure_tab(self) -> QWidget:
        """Build the 3D viewer (or a message when QtWebEngine is unavailable) with its caption and AlphaFold link."""
        sw = QWidget()
        sl = QVBoxLayout(sw)
        sl.setContentsMargins(0, 0, 0, 0)
        sl.setSpacing(0)
        self.web = make_web_view()
        self.struct_msg = label("", "Muted", wrap=True)
        self.struct_msg.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.struct_stack = QStackedWidget()
        self.struct_stack.addWidget(self.struct_msg)
        if self.web is not None:
            self.struct_stack.addWidget(self.web)
        sl.addWidget(self.struct_stack, 1)
        cap = QHBoxLayout()
        cap.setContentsMargins(14, 6, 10, 8)
        self.struct_caption = label("", "Muted", wrap=True)
        cap.addWidget(self.struct_caption, 1)
        self.af_btn = QPushButton(t("open_alphafold"))
        self.af_btn.setObjectName("LinkButton")
        self.af_btn.setIcon(th.icon("external", th.PETROL_700, 14))
        self.af_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.af_btn.clicked.connect(self._open_alphafold)
        cap.addWidget(self.af_btn)
        sl.addLayout(cap)
        return sw

    def _build_lollipop_tab(self) -> QWidget:
        """Build the lollipop chart; clicking a variant selects its row."""
        self.lolli = LollipopChart()
        self.lolli.picked.connect(self.select_index)
        lw = QWidget()
        ll = QVBoxLayout(lw)
        ll.setContentsMargins(10, 10, 10, 10)
        ll.addWidget(self.lolli)
        return lw

    def _build_sequence_tab(self) -> QWidget:
        """Build the wild-type / mutant sequence view with its note and a copy-FASTA button."""
        qw = QWidget()
        ql = QVBoxLayout(qw)
        ql.setContentsMargins(16, 10, 16, 12)
        self.seq_view = QTextBrowser()
        self.seq_view.setOpenExternalLinks(False)
        ql.addWidget(self.seq_view, 1)
        qrow = QHBoxLayout()
        self.seq_note = label("", "Muted", wrap=True)
        qrow.addWidget(self.seq_note, 1)
        cp = QPushButton(f" {t('copy_fasta')}")
        cp.setIcon(th.icon("copy", th.INK_2, 14))
        cp.clicked.connect(self._copy_fasta)
        qrow.addWidget(cp)
        ql.addLayout(qrow)
        return qw

    # ------------------------------------------------------------------ menus & settings
    def _build_menus(self) -> None:
        """Build the File, View (language) and Help menus and the run / find shortcuts."""
        mb = self.menuBar()
        m = mb.addMenu(t("m_file"))
        a = m.addAction(t("m_open"), self.drop.browse)
        a.setShortcut(QKeySequence.StandardKey.Open)
        m.addSeparator()
        for fmt, key, sc in (
            ("xlsx", "m_export_xlsx", "Ctrl+E"),
            ("csv", "m_export_csv", None),
            ("json", "m_export_json", None),
            ("fasta", "m_export_fasta", None),
        ):
            act = m.addAction(t(key), functools.partial(self.export, fmt))
            if sc:
                act.setShortcut(QKeySequence(sc))
        m.addSeparator()
        q = m.addAction(t("m_quit"), self.close)
        q.setShortcut(QKeySequence.StandardKey.Quit)
        v = mb.addMenu(t("m_view"))
        lang = v.addMenu(t("m_language"))
        grp = QActionGroup(self)
        for code, name in (("es", "Español"), ("en", "English")):
            act = QAction(name, self, checkable=True, checked=code == i18n.LANG)
            act.triggered.connect(lambda _=False, c=code: self.language_changed.emit(c))
            grp.addAction(act)
            lang.addAction(act)
        h = mb.addMenu(t("m_help"))
        h.addAction(t("m_formats"), self.show_formats)
        h.addAction(t("m_run_info"), self.show_run_info)
        h.addAction(t("m_open_log"), self.open_log)
        h.addAction(t("m_clear_cache"), self.clear_cache)
        h.addSeparator()
        h.addAction(t("m_about"), self.show_about)
        run = QAction(self)
        run.setShortcuts([QKeySequence("Ctrl+Return"), QKeySequence("Ctrl+Enter")])
        run.triggered.connect(self.run_prediction)
        self.addAction(run)
        find = QAction(self)
        find.setShortcut(QKeySequence.StandardKey.Find)
        find.triggered.connect(lambda: self.search.setFocus())
        self.addAction(find)

    def _restore_settings(self) -> None:
        """Restore the window geometry and the options from the user's preferences."""
        s = self.settings

        def value(key: str, default: object = None) -> Any:
            # Any: QSettings returns what its backend stored (strings, in INI files), hence the conversions below.
            return s.value(key, default)

        geo = value("geometry")
        if geo is not None:
            self.restoreGeometry(geo)
        b = value("build", "auto")
        self.build.setCurrentIndex(max(0, self.build.findData(b)))
        for key, w, default in (
            ("predictors", self.opt_pred, True),
            ("context", self.opt_ctx, True),
            ("clinvar", self.opt_cv, True),
            ("noncoding", self.opt_nc, False),
            ("oncogenicity", self.opt_onc, False),
        ):
            w.setChecked(str(value(key, default)).lower() == "true")
        self.workers.setValue(int(value("workers", 4)))
        self.freq_rule.setCurrentIndex(max(0, self.freq_rule.findData(float(value("common_af", 0.05) or 0))))
        self.acmg_src.setCurrentIndex(max(0, self.acmg_src.findData(value("acmg_source", "model"))))
        self.opt_btn.setChecked(str(value("options_open", False)).lower() == "true")
        self.engine.setCurrentIndex(max(0, self.engine.findData(value("engine", "rest"))))
        self.vep_cmd.setText(value("vep_cmd", "") or "")
        self.vep_args.setText(value("vep_args", "") or "")
        self.offline.setChecked(str(value("offline", False)).lower() == "true")
        self._engine_changed()
        self.drop.start_dir = value("last_dir", str(Path.home()))

    def _save_settings(self) -> None:
        """Save the window geometry and the options to the user's preferences."""
        s = self.settings
        s.setValue("geometry", self.saveGeometry())
        s.setValue("build", self.build.currentData())
        s.setValue("predictors", self.opt_pred.isChecked())
        s.setValue("context", self.opt_ctx.isChecked())
        s.setValue("clinvar", self.opt_cv.isChecked())
        s.setValue("noncoding", self.opt_nc.isChecked())
        s.setValue("oncogenicity", self.opt_onc.isChecked())
        s.setValue("workers", self.workers.value())
        s.setValue("common_af", self.freq_rule.currentData())
        s.setValue("acmg_source", self.acmg_src.currentData())
        s.setValue("options_open", self.opt_btn.isChecked())
        s.setValue("engine", self.engine.currentData())
        s.setValue("vep_cmd", self.vep_cmd.text())
        s.setValue("vep_args", self.vep_args.text())
        s.setValue("offline", self.offline.isChecked())
        s.setValue("last_dir", self.drop.start_dir)

    def closeEvent(self, e: QCloseEvent) -> None:
        """Save the preferences and stop a running prediction (waiting up to 3 s) before closing."""
        self._save_settings()
        if self.worker:
            self.worker.stop = True
        if self.thread and self.thread.isRunning():
            self.thread.quit()
            self.thread.wait(3000)
        super().closeEvent(e)

    def transfer_state(self, other: MainWindow) -> None:
        """Carry inputs and results into a freshly built window (language switch)."""
        self.text.setPlainText(other.text.toPlainText())
        self.gene.setText(other.gene.text())
        self.genes.setText(other.genes.text())
        self.add_files(other.files)
        if other.res is not None:
            self.show_result(other.res)
            if other.current_idx is not None:
                self.select_index(other.current_idx)
        self.setGeometry(other.geometry())

    # ------------------------------------------------------------------ drag & drop on the window
    def dragEnterEvent(self, e: QDragEnterEvent) -> None:
        """Accept drags that carry files or links."""
        if e.mimeData().hasUrls():
            e.acceptProposedAction()

    def dropEvent(self, e: QDropEvent) -> None:
        """Add the dropped local files to the inputs."""
        files = [u.toLocalFile() for u in e.mimeData().urls() if u.isLocalFile()]
        if files:
            self.add_files(files)

    # ------------------------------------------------------------------ files
    def add_files(self, files: Iterable[str]) -> None:
        """Add the existing files of ``files`` that are not yet in the list (the next file dialog opens there)."""
        for f in files:
            if f and Path(f).is_file() and f not in self.files:
                self.files.append(f)
                self.drop.start_dir = str(Path(f).parent)
        self._refresh_files()

    def remove_file(self, f: str) -> None:
        """Remove ``f`` from the input files."""
        if f in self.files:
            self.files.remove(f)
        self._refresh_files()

    def _refresh_files(self) -> None:
        """Rebuild the file chips and the sample menu (the samples of the VCF files, read from their header)."""
        clear_layout(self.file_list)
        samples: list[str] = []
        for f in self.files:
            chip = FileChip(f)
            chip.removed.connect(self.remove_file)
            self.file_list.addWidget(chip)
            try:
                with open(f, "rb") as fh:
                    head = fh.read(4_000_000)
                for s in vcf_samples(f, head):
                    if s not in samples:
                        samples.append(s)
            except (OSError, ValueError):
                pass
        cur = self.sample.currentText()
        self.sample.clear()
        self.sample.addItem(t("all_carriers"))
        self.sample.addItems(samples)
        self.sample.setEnabled(bool(samples))
        self.sample.setVisible(bool(samples))
        self.pass_only.setVisible(bool(self.files))
        if cur in samples:
            self.sample.setCurrentText(cur)
        elif len(samples) == 1:
            self.sample.setCurrentIndex(1)

    # ------------------------------------------------------------------ run
    def _options(self) -> PredictOptions:
        """Return the prediction options chosen in the sidebar."""
        genes = [g.strip() for g in self.genes.text().replace(";", ",").replace(" ", ",").split(",") if g.strip()]
        return PredictOptions(
            assembly=self.build.currentData(),
            default_gene=self.gene.text().strip() or None,
            genes=genes or None,
            sample=self.sample.currentText() if self.sample.currentIndex() > 0 else None,
            pass_only=self.pass_only.isChecked(),
            include_noncoding=self.opt_nc.isChecked(),
            predictors=self.opt_pred.isChecked(),
            protein_context=self.opt_ctx.isChecked(),
            clinvar=self.opt_cv.isChecked(),
            workers=self.workers.value(),
            vep=self.engine.currentData(),
            vep_command=self.vep_cmd.text().strip() or None,
            vep_args=self.vep_args.text().strip() or None,
            offline=self.offline.isChecked() and self.engine.currentData() == "local",
            common_af=self.freq_rule.currentData() or None,
            acmg_source=self.acmg_src.currentData(),
            oncogenicity=self.opt_onc.isChecked(),
        )

    def run_prediction(self) -> None:
        """Start predicting the typed variants and the input files in a worker thread (unless one is running)."""
        if self.thread is not None:
            return
        text = self.text.toPlainText()
        if not text.strip() and not self.files:
            self.toast.show_message(t("need_input"))
            self.text.setFocus()
            return
        self.worker = PredictWorker(text, list(self.files), self._options())
        self.thread = QThread(self)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.failed.connect(self._on_failed)
        self.worker.cancelled.connect(self._on_cancelled)
        for sig in (self.worker.finished, self.worker.failed, self.worker.cancelled):
            sig.connect(self.thread.quit)
        self.thread.finished.connect(self._thread_done)
        self._busy(True)
        self.thread.start()

    def cancel_prediction(self) -> None:
        """Ask the running prediction to stop (at its next progress report)."""
        if self.worker:
            self.worker.stop = True
            self.stage.setText(t("cancelling"))

    def _busy(self, on: bool) -> None:
        """Switch the sidebar footer between the predict button and the progress display."""
        self.run_btn.setEnabled(not on)
        self.progress.setVisible(on)
        self.stage.setVisible(on)
        self.cancel_btn.setVisible(on)
        self.progress.setValue(0)
        if on:
            self.stage.setText(tr_stage("Reading input"))
            self.progress.setMaximum(0)  # indeterminate until the first count arrives

    def _on_progress(self, stage: str, done: int, total: int) -> None:
        """Show the stage of the prediction and its progress."""
        self.progress.setMaximum(1000)
        self.progress.setValue(int(1000 * done / max(total, 1)))
        suffix = f" · {done}/{total}" if total > 1 else ""
        self.stage.setText(tr_stage(stage) + suffix)

    def _on_finished(self, res: PredictionResult) -> None:
        """Show the finished prediction."""
        self.show_result(res)

    def _on_failed(self, msg: str) -> None:
        """Tell the user the prediction failed, and where the log is."""
        from .logs import log_file

        QMessageBox.warning(self, t("failed"), f"{tr_text(msg)}\n\n{log_file()}")

    def _on_cancelled(self) -> None:
        """Confirm the cancellation."""
        self.toast.show_message(t("cancel"))

    def _thread_done(self) -> None:
        """Clean up after the worker thread has finished."""
        assert self.worker is not None  # both set by run_prediction
        assert self.thread is not None
        self._busy(False)
        self.worker.deleteLater()
        self.thread.deleteLater()
        self.worker = None
        self.thread = None

    # ------------------------------------------------------------------ show results
    def show_result(self, res: PredictionResult) -> None:
        """Show a prediction: header, KPIs, warnings, status filters, the table and the first row's details."""
        self.res = res
        df = res.table.copy()
        df["status"] = df["prediction"].map(prediction_status) if len(df) else []
        self.df = df
        self.stack.setCurrentIndex(1)
        self._show_summary(res, df)
        self._show_banners(res)
        self._show_segments(df)
        self._show_table(df)

    def _show_summary(self, res: PredictionResult, df: pd.DataFrame) -> None:
        """Fill the header (counts, genome build, annotation source, run time), the KPIs and the errors button."""
        s = res.stats
        self._show_header_tools(True)
        self.title.setText(t("results"))
        run = res.run or {}
        engine = (
            t("local_label")
            if str(run.get("annotation", "")).startswith("local")
            else f"Ensembl {run.get('ensembl_release')}"
            if run.get("ensembl_release")
            else "Ensembl"
        )
        self.subtitle.setText(
            t(
                "subtitle",
                inputs=s.get("inputs", 0),
                rows=s.get("rows", 0),
                build=res.assembly,
                engine=engine,
                secs=s.get("seconds", 0),
            )
        )
        self.k_inputs.set(s.get("inputs", 0))
        self.k_rows.set(s.get("rows", 0))
        self.k_dmg.set(s.get("likely_damaging", 0))
        cvp = 0
        if len(df):
            c = df["clinvar"].fillna("").str.lower()
            cvp = int(
                (c.str.contains("pathogenic") & ~c.str.contains("conflicting") & ~c.str.contains("submissions")).sum()
            )
        self.k_cv.set(cvp)
        n_err = len(res.errors)
        self.err_btn.setText(f" {t('not_analysed', n=n_err)}")
        self.err_btn.setVisible(n_err > 0)

    def _show_banners(self, res: PredictionResult) -> None:
        """Show the run's warnings as banners (with a retry button for services that did not respond)."""
        s = res.stats
        clear_layout(self.banners)
        for w in res.warnings:
            self.banners.addWidget(self._banner(tr_text(w), retry="did not respond" in w))
        skipped = [
            t(k, n=s[key])
            for key, k in (
                ("no_protein_effect", "skip_noprot"),
                ("file_skipped_not_carried", "skip_carried"),
                ("file_skipped_outside_genes", "skip_genes"),
                ("file_skipped_symbolic", "skip_symbolic"),
                ("file_skipped_filter", "skip_filter"),
            )
            if s.get(key)
        ]
        self.skipped.setText(t("skipped", items=" · ".join(skipped)) if skipped else "")

    def _show_segments(self, df: pd.DataFrame) -> None:
        """Add one checkable filter button per prediction status present, with its count."""
        clear_layout(self.seg_box)
        self.segments = {}
        counts = df["status"].value_counts().to_dict() if len(df) else {}
        for st in th.STATUS_ORDER:
            if not counts.get(st):
                continue
            b = QPushButton(f"  {t(st)}  {counts[st]}")
            b.setObjectName("Segment")
            b.setCheckable(True)
            b.setChecked(True)
            b.setIcon(self._status_icon(st))
            b.setCursor(Qt.CursorShape.PointingHandCursor)
            b.toggled.connect(self._filter_changed)
            self.seg_box.addWidget(b)
            self.segments[st] = b

    def _show_table(self, df: pd.DataFrame) -> None:
        """Load ``df`` into the table (filters cleared, empty columns hidden, original order) and select row 1."""
        self.model.set_frame(df)
        self.proxy.set_statuses(None)
        self.search.blockSignals(True)
        self.search.clear()
        self.search.blockSignals(False)
        self.proxy.set_query("")
        for c in range(len(COLUMNS)):
            self.table.setColumnHidden(c, self.model.column_empty(c))
        self.table.sortByColumn(-1, Qt.SortOrder.AscendingOrder)
        self.proxy.sort(-1)
        has_rows = len(df) > 0
        if not self._sized:
            total_h = max(500, self.vsplit.height())
            self.vsplit.setSizes([int(total_h * 0.36), int(total_h * 0.64)])
            total_w = max(800, self.hsplit.width())
            self.hsplit.setSizes([int(total_w * 0.5), int(total_w * 0.5)])
            self._sized = True
        self.table.setVisible(has_rows)
        self.empty_rows.setVisible(not has_rows)
        self.hsplit.setVisible(has_rows)
        if has_rows:
            self.table.selectRow(0)
            self._row_changed(self.proxy.index(0, 0), QModelIndex())
        self.table.setFocus()

    def _status_icon(self, st: str) -> QIcon:
        """Return the mark of status ``st`` as a button icon."""
        pm = QPixmap(28, 28)
        pm.fill(Qt.GlobalColor.transparent)
        pm.setDevicePixelRatio(2)
        p = QPainter(pm)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        paint_status_mark(p, st, QPointF(7, 7.5), 4.6)
        p.end()
        return QIcon(pm)

    def _banner(self, text: str, retry: bool = False) -> QFrame:
        """Return a warning banner showing ``text``, with a button that runs the prediction again if ``retry``."""
        f = QFrame()
        f.setObjectName("Banner")
        lay = QHBoxLayout(f)
        lay.setContentsMargins(14, 9, 14, 9)
        ic = QLabel()
        ic.setPixmap(th.icon("alert", "#9A7414", 16).pixmap(16, 16))
        lay.addWidget(ic)
        lay.addWidget(label(text, wrap=True), 1)
        if retry:
            b = QPushButton(f" {t('retry')}")
            b.setIcon(th.icon("reset", th.INK_2, 14))
            b.setCursor(Qt.CursorShape.PointingHandCursor)
            b.clicked.connect(self.run_prediction)
            lay.addWidget(b)
        return f

    def _filter_changed(self, *_: object) -> None:
        """Apply the status buttons and the search text to the table (all statuses if none is checked)."""
        if self.segments:
            on = {k for k, b in self.segments.items() if b.isChecked()}
            self.proxy.set_statuses(on if on else set(self.segments))
        self.proxy.set_query(self.search.text())
        if self.proxy.rowCount() and not self.table.selectionModel().hasSelection():
            self.table.selectRow(0)

    def _row_changed(self, cur: QModelIndex, _prev: QModelIndex | None = None) -> None:
        """Show the details of the table row that became current."""
        if not cur.isValid():
            return
        src = self.proxy.mapToSource(cur)
        idx = self.model.index_labels[src.row()]
        self.show_detail(idx)

    def select_index(self, idx: int) -> None:
        """Select the table row of dataframe index `idx` (e.g. from the lollipop chart)."""
        if idx not in self.model.index_labels:
            return
        src_row = self.model.index_labels.index(idx)
        prow = self.proxy.mapFromSource(self.model.index(src_row, 0))
        if not prow.isValid():  # filtered out: clear filters first
            for b in self.segments.values():
                b.blockSignals(True)
                b.setChecked(True)
                b.blockSignals(False)
            self.search.clear()
            self.proxy.set_statuses(None)
            prow = self.proxy.mapFromSource(self.model.index(src_row, 0))
        self.table.selectionModel().setCurrentIndex(
            prow, QItemSelectionModel.SelectionFlag.ClearAndSelect | QItemSelectionModel.SelectionFlag.Rows
        )
        self.table.scrollTo(prow)

    # ------------------------------------------------------------------ detail
    def show_detail(self, idx: int) -> None:
        """Fill the detail pane and the right-hand tabs for the variant with dataframe index ``idx``."""
        self.current_idx = idx
        r = cast(Row, self.df.loc[idx].to_dict())
        clear_layout(self.dl)
        self._detail_header(r)
        self._detail_why(r)
        self._detail_acmg(r)
        self._detail_somatic(r)
        self._detail_predictors(r)
        self._detail_physchem(r)
        self._detail_sequence_window(r)
        self._detail_context(r)
        self._detail_clinical(r)
        self._detail_links(r)
        self.dl.addStretch(1)

        self._show_structure(r)
        self._show_lollipop(r, idx)
        self._show_sequence(r, idx)

    def _detail_header(self, r: Row) -> None:
        """Add the title (gene and HGVS p.), the prediction with its scores, the consequence and the input."""
        g = _getter(r)
        status = r.get("status") or "Not assessed"

        head = QHBoxLayout()
        head.setSpacing(10)
        title = f"{g('gene', '')} {g('hgvsp', '') or ''}".strip() or str(g("input", ""))
        h = label(title, "H2", selectable=True)
        h.setStyleSheet("font-size: 19px;")
        head.addWidget(h)
        if g("protein_change"):
            head.addWidget(label(g("protein_change"), "Chip"))
        head.addStretch(1)
        self.dl.addLayout(head)

        srow = QHBoxLayout()
        srow.setSpacing(8)
        srow.addWidget(StatusDot(status, 15))
        pl = label(tr_prediction(g("prediction", "")))
        pl.setStyleSheet(f"font-weight: 600; color: {th.INK}; font-size: 13.5px;")
        srow.addWidget(pl)
        if isnum(r.get("effect_score")):
            srow.addWidget(label(f"· {t('col_effect').lower()} {r['effect_score']:.2f}", "Muted"))
        if isnum(r.get("score_with_frequency")):
            fw = label(f"· {t('col_with_freq').lower()} {r['score_with_frequency']:.2f}", "Muted")
            fw.setToolTip(t("with_freq_tip"))
            srow.addWidget(fw)
        srow.addStretch(1)
        self.dl.addLayout(srow)
        meta = [
            x
            for x in (
                str(g("consequence", "")).replace("_variant", "").replace("_", " ").replace(",", ", "),
                g("transcript"),
                g("hgvsc"),
            )
            if x
        ]
        mrow = label(" · ".join(meta), "Muted", wrap=True, selectable=True)
        self.dl.addWidget(mrow)
        self.dl.addWidget(label(str(g("input", "")), "Mono", selectable=True))

    def _detail_why(self, r: Row) -> None:
        """Add the reasons for the prediction and, for missense variants, the model's predictor contributions."""
        g = _getter(r)
        basis = g("prediction_basis")
        if basis:
            summary = basis.split(" | ")[0]
            self.dl.addWidget(section(t("why")))
            self.dl.addWidget(label(tr_text(summary), "Body", wrap=True, selectable=True))
        ms = model.score(r) if "missense" in str(g("consequence", "")) else None
        if ms is not None and ms.contributions:
            self.dl.addWidget(ContributionChart(ms.contributions, t("contrib_damaging"), t("contrib_benign")))
            self.dl.addWidget(label(t("model_note"), "Muted", wrap=True))

    def _detail_acmg(self, r: Row) -> None:
        """Add the ACMG evidence codes with their points, and the notes on how they were set."""
        g = _getter(r)
        ev = [
            (t("acmg_codes"), g("acmg_codes"), None),
            (t("acmg_pp3"), g("acmg_pp3_bp4"), g("acmg_points")),
            (t("acmg_ps1pm5"), g("acmg_ps1_pm5"), _code_points(g("acmg_ps1_pm5"))),
            (t("acmg_pm2"), g("acmg_pm2"), _code_points(g("acmg_pm2"))),
            (t("acmg_freq"), g("population_frequency"), None),
        ]
        ev = [(k, v, pts) for k, v, pts in ev if v]
        if ev:
            self.dl.addWidget(section(t("acmg")))
            for k, v, pts in ev:
                extra = f" · {t('points', n=f'{int(pts):+d}')}" if isnum(pts) and pts else ""
                self.dl.addWidget(_rich(f"<b>{k}:</b> {html.escape(tr_text(v))}{extra}"))
            note = ACMG_NOTE_KEY[acmg.thresholds_used(g("acmg_pp3_bp4"))]
            self.dl.addWidget(label(t(note), "Muted", wrap=True))
            if g("acmg_ps1_pm5") or g("acmg_pm2"):
                self.dl.addWidget(label(t("acmg_note_residue"), "Muted", wrap=True))

    def _detail_somatic(self, r: Row) -> None:
        """Add the somatic oncogenicity class, its criteria and evidence, and ClinVar's call on this variant."""
        g = _getter(r)
        if not g("oncogenicity"):
            return
        self.dl.addWidget(section(t("somatic")))
        pts = r.get("oncogenicity_points")
        extra = f" · {t('points', n=f'{int(pts):+d}')}" if isnum(pts) else ""
        lines = [(t("somatic_class"), tr_text(g("oncogenicity")) + extra)]
        if g("oncogenicity_codes"):
            lines.append((t("somatic_codes"), g("oncogenicity_codes")))
        if g("clinvar_oncogenicity"):
            lines.append((t("somatic_clinvar"), g("clinvar_oncogenicity")))
        if g("clinical_evidence"):
            lines.append((t("somatic_clinical"), tr_text(g("clinical_evidence"))))
        for k, v in lines:
            self.dl.addWidget(_rich(f"<b>{k}:</b> {html.escape(v)}"))
        if g("oncogenicity_evidence"):
            self.dl.addWidget(label(tr_text(g("oncogenicity_evidence")), "Body", wrap=True, selectable=True))
        self.dl.addWidget(label(t("somatic_note"), "Muted", wrap=True))

    def _detail_predictors(self, r: Row) -> None:
        """Add one bar per predictor with a score or a call (and the popEVE score, which has no fixed scale)."""
        rows = []
        for name, col, call_col, rng, thr, hb in PREDICTORS:
            v, call = r.get(col), (r.get(call_col) if call_col else None)
            v = v if isnum(v) else None
            call = call if isinstance(call, str) and call else None
            if v is not None or call:
                rows.append(PredictorRow(name, v, call, rng, thr, hb))
        if rows:
            self.dl.addWidget(section(t("predictors")))
            for w in rows:
                self.dl.addWidget(w)
            if isnum(r.get("popeve")):
                self.dl.addWidget(
                    label(
                        f"popEVE {r['popeve']:.2f} · " + t("popeve_note"),
                        "Muted",
                    )
                )

    def _detail_physchem(self, r: Row) -> None:
        """Add the amino-acid change and its physico-chemical facts (Grantham, BLOSUM62, hydropathy, ...)."""
        g = _getter(r)
        ra, aa = g("ref_aa"), g("alt_aa")
        if isnum(r.get("grantham")) and ra and aa:
            self.dl.addWidget(section(t("physchem")))
            classes = str(g("aa_class_change", "")).split(" -> ")
            cls = " → ".join(aa_class(c) for c in classes) if len(classes) == 2 else ""
            self.dl.addWidget(
                label(
                    f"<b>{aa_name(ra)} → {aa_name(aa)}</b>"
                    + (f"  <span style='color:{th.MUTED}'>({cls})</span>" if cls else ""),
                    "Body",
                )
            )
            grid = QGridLayout()
            grid.setHorizontalSpacing(10)
            grid.setVerticalSpacing(8)
            facts = [
                (t("grantham"), f"{int(r['grantham'])}", tr_text(g("grantham_class", ""))),
                (t("blosum"), f"{int(r['blosum62']):+d}" if isnum(r.get("blosum62")) else "—", ""),
                (
                    t("hydropathy"),
                    f"{r['hydropathy_change']:+.1f}" if isnum(r.get("hydropathy_change")) else "—",
                    "Kyte–Doolittle",
                ),
                (t("volume"), f"{r['volume_change']:+.0f} Å³" if isnum(r.get("volume_change")) else "—", ""),
                (t("charge"), f"{int(r['charge_change']):+d}" if isnum(r.get("charge_change")) else "—", ""),
            ]
            for i, (k, v, sub) in enumerate(facts):
                grid.addWidget(_fact_box(v, k + (f" · {sub}" if sub else "")), i // 3, i % 3)
            self.dl.addLayout(grid)
            if g("physchem_notes"):
                self.dl.addWidget(label(tr_text(g("physchem_notes")), "Muted", wrap=True))

    def _detail_sequence_window(self, r: Row) -> None:
        """Add the protein sequence around the change (change highlighted) and the notes on the mutant protein."""
        g = _getter(r)
        if g("protein_window"):
            self.dl.addWidget(section(t("sequence")))
            w = html.escape(g("protein_window"))
            w = w.replace("[", f"<span style='background:{th.ACCENT};color:white;'>&nbsp;").replace(
                "]", "&nbsp;</span>"
            )
            seq = label(w)
            seq.setTextFormat(Qt.TextFormat.RichText)
            seq.setStyleSheet(f"font-family:'{th.mono_family()}';font-size:15px;letter-spacing:1.5px;color:{th.INK};")
            seq.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            self.dl.addWidget(seq)
            note = g("mutant_protein")
            if note and note != "ok":
                self.dl.addWidget(label(tr_text(note), "Muted", wrap=True))
            if g("nmd"):
                self.dl.addWidget(label(tr_text(g("nmd")), "Muted", wrap=True))

    def _detail_context(self, r: Row) -> None:
        """Add the protein context: protein and residue, key site, region, pLDDT, UniProt notes and function."""
        g = _getter(r)
        self.dl.addWidget(section(t("context")))
        facts = []
        if g("protein_name"):
            pos = int(r["protein_position"]) if isnum(r.get("protein_position")) else None
            where = (
                t("residue_of", pos=pos, len=int(r["protein_length"]))
                if pos and isnum(r.get("protein_length"))
                else t("residue", pos=pos)
                if pos
                else ""
            )
            facts.append(
                f"<b>{html.escape(g('protein_name'))}</b>"
                + (f" · {where}" if where else "")
                + (f" · <span style='color:{th.MUTED}'>{g('uniprot')}</span>" if g("uniprot") else "")
            )
        if g("key_site"):
            facts.append(f"<b>{t('key_site')}:</b> {html.escape(g('key_site'))}")
        if g("protein_features"):
            facts.append(f"<b>{t('region')}:</b> {html.escape(g('protein_features'))}")
        if isnum(r.get("plddt")):
            pc = next(c for lim, c, _ in th.PLDDT_COLORS if r["plddt"] >= lim)
            facts.append(
                f"<b>{t('plddt')}:</b> <span style='color:{pc}'>●</span> {r['plddt']:.0f} "
                f"({tr_text(g('plddt_class', ''))})"
            )
        if facts:
            for f in facts:
                self.dl.addWidget(_rich(f))
        else:
            self.dl.addWidget(label(t("no_context"), "Muted"))
        if g("notes"):
            self.dl.addWidget(label(tr_text(g("notes")), "Muted", wrap=True))
        for key, tk in (("uniprot_variants_here", "uniprot_variants"), ("uniprot_mutagenesis_here", "mutagenesis")):
            if g(key):
                items = "<br>".join("• " + html.escape(x) for x in str(g(key)).split(" | "))
                self.dl.addWidget(Collapsible(t(tk), items))
        if g("protein_function"):
            self.dl.addWidget(Collapsible(t("function"), html.escape(g("protein_function"))))

    def _detail_clinical(self, r: Row) -> None:
        """Add ClinVar (classification, stars, review, conditions), gnomAD frequency and the genotype."""
        g = _getter(r)
        self.dl.addWidget(section(t("clinical")))
        cl = []
        if g("clinvar"):
            stars = (
                "★" * int(r["clinvar_stars"]) + "☆" * (4 - int(r["clinvar_stars"]))
                if isnum(r.get("clinvar_stars"))
                else ""
            )
            txt = (
                f"<b>{t('clinvar')}:</b> {html.escape(g('clinvar'))} "
                f"<span style='color:{th.STATUS_COLOR['Uncertain']}'>{stars}</span>"
            )
            if g("clinvar_review"):
                txt += f"<br><span style='color:{th.MUTED}'>{html.escape(g('clinvar_review'))}</span>"
            cl.append(txt)
        if g("clinvar_oncogenicity"):
            cl.append(f"<b>{t('oncogenicity')}:</b> {html.escape(g('clinvar_oncogenicity'))}")
        if g("clinvar_conditions"):
            cl.append(f"<b>{t('conditions')}:</b> {html.escape(g('clinvar_conditions'))}")
        if isnum(r.get("gnomad_af")):
            popmax = r.get("gnomad_popmax_af")
            pm = f" · {t('popmax', v=format(popmax, '.2e'))}" if isnum(popmax) else ""
            cl.append(f"<b>{t('gnomad')}:</b> {r['gnomad_af']:.2e}{pm}")
        else:
            cl.append(t("gnomad_none"))
        if g("zygosity"):
            cl.append(
                f"<b>{t('zygosity')}:</b> {tr_text(g('zygosity'))} "
                f"({g('genotype', '')}, {html.escape(str(g('sample', '')))})"
            )
        elif g("carriers"):
            cl.append(f"<b>{t('carriers')}:</b> {html.escape(str(g('carriers')))}")
        for c in cl:
            self.dl.addWidget(_rich(c))

    def _detail_links(self, r: Row) -> None:
        """Add buttons that open the variant in external databases (see :meth:`_links`)."""
        links = self._links(r)
        if links:
            self.dl.addWidget(section(t("links")))
            box = QWidget()
            flow = FlowLayout(box, spacing=4)
            for name, url in links:
                b = QPushButton(name)
                b.setObjectName("LinkButton")
                b.setIcon(th.icon("external", th.PETROL_700, 13))
                b.setCursor(Qt.CursorShape.PointingHandCursor)
                b.setToolTip(url)
                b.clicked.connect(lambda _=False, u=url: QDesktopServices.openUrl(QUrl(u)))
                flow.addWidget(b)
            self.dl.addWidget(box)

    def _links(self, r: Row) -> list[tuple[str, str]]:
        """Return (name, URL) of the variant's pages: UniProt, AlphaFold, ClinVar, dbSNP, gnomAD, Ensembl, PubMed."""
        g = _getter(r)
        out = []
        if g("uniprot"):
            out += [
                ("UniProt", f"https://www.uniprot.org/uniprotkb/{g('uniprot')}/entry"),
                ("AlphaFold", f"https://alphafold.ebi.ac.uk/entry/{g('uniprot')}"),
            ]
        if g("clinvar_ids") and str(g("clinvar_ids")).startswith("VCV"):
            out.append(("ClinVar", f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{int(str(g('clinvar_ids'))[3:])}/"))
        if g("rsid"):
            out.append(("dbSNP", f"https://www.ncbi.nlm.nih.gov/snp/{str(g('rsid')).split(',')[0].strip()}"))
        ref, alt = g("ref"), g("alt")
        if (
            ref
            and alt
            and "-" not in (ref + alt)
            and self.res
            and self.res.assembly == "GRCh38"
            and isnum(r.get("pos"))
        ):
            out.append(
                (
                    "gnomAD",
                    f"https://gnomad.broadinstitute.org/variant/{g('chrom')}-{int(r['pos'])}-{ref}-{alt}"
                    "?dataset=gnomad_r4",
                )
            )
        if g("transcript"):
            host = "grch37.ensembl.org" if self.res and self.res.assembly == "GRCh37" else "www.ensembl.org"
            out.append(("Ensembl", f"https://{host}/Homo_sapiens/Transcript/Summary?t={g('transcript')}"))
        if g("pubmed_ids"):
            ids = str(g("pubmed_ids")).replace(" ", "").split(",")
            n = int(r["pubmed_count"]) if isnum(r.get("pubmed_count")) else len(ids)
            out.append((f"PubMed ({n})", "https://pubmed.ncbi.nlm.nih.gov/?term=" + quote(" OR ".join(ids[:25]))))
        return out

    # ------------------------------------------------------------------ right-hand tabs
    def _show_structure(self, r: Row) -> None:
        """Show the AlphaFold model with the residue marked (loading it in the background), or why it cannot."""
        acc = val(r.get("uniprot"))
        af = self.res.alphafold.get(acc) if (self.res and acc) else None
        pos = int(r["protein_position"]) if isnum(r.get("protein_position")) else None
        self.af_btn.setVisible(bool(acc))
        self._af_acc = acc
        if not (af and af.get("available") and pos and isnum(r.get("plddt"))):
            self.struct_msg.setText(t("no_structure"))
            self.struct_stack.setCurrentIndex(0)
            self.struct_caption.setText("")
            self.pending_structure = None
            return
        lbl = str(val(r.get("protein_change"), pos))
        self.struct_caption.setText(t("structure_caption", id=af.get("entry_id")))
        if self.web is None:
            self.struct_msg.setText(t("no_webengine"))
            self.struct_stack.setCurrentIndex(0)
            return
        url = af["pdb_url"]
        self.pending_structure = (url, pos, lbl)
        if url in self.pdb_cache:
            self._render_structure()
        else:
            self.struct_msg.setText(t("loading_structure"))
            self.struct_stack.setCurrentIndex(0)
            QThreadPool.globalInstance().start(PdbLoader(url, self.pdb_signals, self.tmp_dir))

    def _pdb_ready(self, url: str, txt: str | None) -> None:
        """Store a downloaded PDB file (None: failed) and show it if it is still the one wanted."""
        self.pdb_cache[url] = txt
        if self.pending_structure and self.pending_structure[0] == url:
            self._render_structure()

    def _render_structure(self) -> None:
        """Show the pending structure in the 3D viewer, or a message if its PDB file could not be loaded."""
        assert self.pending_structure is not None  # only set when there is a web view (see _show_structure)
        assert self.web is not None
        url, pos, lbl = self.pending_structure
        pdb = self.pdb_cache.get(url)
        if not pdb:
            self.struct_msg.setText(t("no_structure"))
            self.struct_stack.setCurrentIndex(0)
            return
        load_page(self.web, structure_page(pdb, pos, lbl), self.tmp_dir)
        self.struct_stack.setCurrentIndex(1)

    def _open_alphafold(self) -> None:
        """Open the AlphaFold DB entry of the current protein in the browser."""
        if self._af_acc:
            QDesktopServices.openUrl(QUrl(f"https://alphafold.ebi.ac.uk/entry/{self._af_acc}"))

    def _show_lollipop(self, r: Row, idx: int) -> None:
        """Show all results of the variant's gene along its protein, with its UniProt domains; select ``idx``."""
        gene = val(r.get("gene"))
        d = self.df[(self.df["gene"] == gene) & self.df["protein_position"].notna()] if gene else self.df.iloc[0:0]
        entry = (self.res.uniprot.get(val(r.get("uniprot"), "")) or {}) if self.res else {}
        length = int(r["protein_length"]) if isnum(r.get("protein_length")) else entry.get("length")
        if not length and len(d):
            length = int(d["protein_position"].max() * 1.05) + 1
        domains: list[Domain] = []
        if "skipped" not in str(val(r.get("notes"), "")):
            for f in entry.get("features", []):
                if f["type"] in ("Domain", "Zinc finger", "DNA binding", "Transmembrane"):
                    domains.append({"name": f["description"] or f["type"], "start": f["start"], "end": f["end"]})
        pts: list[LollipopPoint] = []
        for i, x in d.iterrows():
            pts.append(
                {
                    "idx": cast(int, i),  # the results table has a RangeIndex
                    "pos": int(x["protein_position"]),
                    "score": float(x["effect_score"]) if isnum(x["effect_score"]) else None,
                    "status": x["status"],
                    "label": str(val(x["protein_change"], "")),
                    "prediction": val(x["prediction"], ""),
                    "input": str(val(x["input"], "")),
                }
            )
        title = f"{gene} · {val(entry.get('name'), '')}" if gene else ""
        self.lolli.set_data(pts, domains, length or 1, idx, title.strip(" ·"))

    def _show_sequence(self, r: Row, idx: int) -> None:
        """Show the wild-type and mutant protein around the change, and keep their FASTA for the copy button."""
        seq = self.res.sequences.get(idx) if self.res else None
        pos = int(r["protein_position"]) if isnum(r.get("protein_position")) else None
        if not seq or not pos:
            self.seq_view.setHtml(f"<p style='color:{th.MUTED}'>{t('no_sequence')}</p>")
            self.seq_note.setText("")
            self._fasta = None
            return
        ra = str(val(r.get("ref_aa"), "") or "")
        aa = str(val(r.get("alt_aa"), "") or "")
        n_ref = len(ra.replace("-", "")) or 1
        n_alt = len(aa.replace("-", "").replace("*", "")) or (0 if "del" in str(r.get("protein_change")) else 1)
        self.seq_view.setHtml(sequence_html(seq["wt"], seq["mut"], pos, n_ref, n_alt))
        note = seq.get("note")
        self.seq_note.setText(tr_text(note) if note and note != "ok" else "")
        hdr = seq["header"]
        self._fasta = f">WT|{hdr}\n{_wrap(seq['wt'])}\n>MUT|{hdr}\n{_wrap(seq['mut'])}\n"

    def _copy_fasta(self) -> None:
        """Copy the wild-type and mutant FASTA of the current variant to the clipboard."""
        if self._fasta:
            QGuiApplication.clipboard().setText(self._fasta)
            self.toast.show_message(t("copied"))

    # ------------------------------------------------------------------ dialogs & actions
    def show_errors(self) -> None:
        """Show the inputs that could not be analysed, with the reason, in a dialog."""
        if not self.res or not len(self.res.errors):
            return
        dlg = QDialog(self)
        dlg.setWindowTitle(t("errors_title"))
        dlg.resize(720, 420)
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(20, 18, 20, 18)
        lay.addWidget(label(t("errors_title"), "H2"))
        tb = QTableView()
        m = QStandardItemModel(0, 3, tb)
        m.setHorizontalHeaderLabels([t("col_input").upper(), t("col_source").upper(), t("col_error").upper()])
        for _, e in self.res.errors.iterrows():
            m.appendRow(
                [
                    QStandardItem(str(e["input"])),
                    QStandardItem(str(e["source"])),
                    QStandardItem(tr_text(str(e["error"]))),
                ]
            )
        tb.setModel(m)
        tb.verticalHeader().hide()
        tb.setShowGrid(False)
        tb.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        tb.horizontalHeader().setStretchLastSection(True)
        tb.setColumnWidth(0, 200)
        tb.setColumnWidth(1, 120)
        lay.addWidget(tb)
        dlg.exec()

    def export(self, fmt: str) -> None:
        """Ask for a file name and save the results (``xlsx``, ``csv``, ``json``) or the proteins (``fasta``)."""
        if not self.res:
            return
        if fmt == "fasta" and not self.res.sequences:
            self.toast.show_message(t("no_sequence"))
            return
        ext = {"xlsx": "xlsx", "csv": "csv", "json": "json", "fasta": "fasta"}[fmt]
        name = "missensexplorer_proteins.fasta" if fmt == "fasta" else f"missensexplorer.{ext}"
        start = str(Path(self.drop.start_dir) / name)
        filters = {
            "xlsx": "Excel (*.xlsx)",
            "csv": "CSV (*.csv)",
            "json": "JSON (*.json)",
            "fasta": "FASTA (*.fasta *.fa)",
        }
        path, _ = QFileDialog.getSaveFileName(self, t("export"), start, filters[fmt])
        if not path:
            return
        if not Path(path).suffix:
            path += f".{ext}"
        try:
            if fmt == "fasta":
                Path(path).write_text(to_fasta(self.res))
            else:
                write(self.res, path)
        except OSError as e:
            QMessageBox.warning(self, t("export"), str(e))
            return
        self.drop.start_dir = str(Path(path).parent)
        self.toast.show_message(t("saved", path=Path(path).name))

    def show_run_info(self) -> None:
        """Show the provenance of the run (versions, options, sources) in a dialog."""
        if not self.res:
            return
        from ..report import run_info

        info = run_info(self.res)
        dlg = QDialog(self)
        dlg.setWindowTitle(t("run_info_title"))
        dlg.resize(820, 520)
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(20, 18, 20, 18)
        lay.addWidget(label(t("run_info_title"), "H2"))
        view = QTextBrowser()
        rows = "".join(
            f"<tr><td style='color:{th.MUTED};padding:4px 14px 4px 0;vertical-align:top'>{html.escape(str(k))}</td>"
            f"<td style='padding:4px 0'>{html.escape(tr_text(str(v)))}</td></tr>"
            for k, v in info.itertuples(index=False)
        )
        view.setHtml(f"<table cellspacing='0'>{rows}</table>")
        lay.addWidget(view)
        dlg.exec()

    def open_log(self) -> None:
        """Open the log file (or its folder, if there is no log yet) in the system's viewer."""
        from .logs import log_file

        f = log_file()
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(f if f.exists() else f.parent)))

    def show_formats(self) -> None:
        """Show the accepted input formats in a dialog."""
        box = QMessageBox(self)
        box.setWindowTitle(t("m_formats"))
        box.setTextFormat(Qt.TextFormat.RichText)
        box.setText(f"<h3 style='font-weight:500'>{t('m_formats')}</h3>{t('formats_html')}")
        box.exec()

    def show_about(self) -> None:
        """Show the About dialog (version, data sources)."""
        box = QMessageBox(self)
        box.setWindowTitle(t("m_about"))
        box.setIconPixmap(th.app_icon().pixmap(72, 72))
        box.setTextFormat(Qt.TextFormat.RichText)
        box.setText(t("about", v=__version__))
        box.exec()

    def clear_cache(self) -> None:
        """Delete the on-disk cache of web-service answers and the loaded structures."""
        with Cache() as cache:
            cache.clear()
        self.pdb_cache.clear()
        self.toast.show_message(t("cache_cleared"))

    def resizeEvent(self, e: QResizeEvent) -> None:
        """Keep a visible toast at the bottom centre."""
        super().resizeEvent(e)
        if self.toast.isVisible():
            self.toast.place()


def _wrap(s: str, n: int = 60) -> str:
    """Split ``s`` into lines of ``n`` characters (FASTA)."""
    return "\n".join(s[i : i + n] for i in range(0, len(s), n))
