"""Visual identity of the desktop app: desaturated petrol / mineral green on white, Helvetica."""

from __future__ import annotations

from PySide6.QtCore import QByteArray, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QFontDatabase, QIcon, QPainter, QPalette, QPixmap
from PySide6.QtSvg import QSvgRenderer

# --------------------------------------------------------------------------- colour tokens

PETROL_950 = "#10262A"
PETROL_900 = "#17333A"   # sidebar
PETROL_800 = "#21434A"   # sidebar inputs
PETROL_700 = "#2B5659"   # primary
PETROL_600 = "#376669"   # primary hover
PETROL_500 = "#4D7C7B"
MINERAL_400 = "#7C9F99"
MINERAL_300 = "#A7C1BB"
MINERAL_200 = "#CCDDD8"
MINERAL_100 = "#E5EEEB"
MINERAL_50 = "#F3F7F6"   # app background
WHITE = "#FFFFFF"
INK = "#1A2829"
INK_2 = "#4A5B5C"
MUTED = "#7E8D8C"
BORDER = "#DDE6E3"
BORDER_STRONG = "#C5D3CF"
ACCENT = "#AD3F78"       # the mutated residue (sequence window, 3D structure)

# Prediction status: reserved colours, always shown with a shape and a label.
# Validated (dataviz validate_palette.js, white surface): CVD ΔE ≥ 12.9, normal-vision ΔE ≥ 17.9.
STATUS_ORDER = ["Damaging", "Uncertain", "Benign / silent", "Not assessed"]
STATUS_COLOR = {
    "Damaging": "#B8433A",
    "Uncertain": "#D29A1E",
    "Benign / silent": "#2E8752",
    "Not assessed": "#9AA6A4",
}
STATUS_SHAPE = {"Damaging": "triangle", "Uncertain": "diamond", "Benign / silent": "circle",
                "Not assessed": "square"}

# AlphaFold pLDDT convention.
PLDDT_COLORS = [(90, "#0053D6", "≥ 90"), (70, "#65CBF3", "70–90"), (50, "#FFDB13", "50–70"), (0, "#FF7D45", "< 50")]

FONT_FAMILIES = ["Helvetica Neue", "Helvetica", "Arial"]
MONO_FAMILIES = ["Menlo", "SF Mono", "Monaco", "Consolas", "DejaVu Sans Mono"]


def font(size: float = 13, weight: QFont.Weight = QFont.Weight.Normal, mono: bool = False) -> QFont:
    f = QFont()
    f.setFamilies(MONO_FAMILIES if mono else FONT_FAMILIES)
    f.setPixelSize(max(8, round(size)))  # pixels, so painted text matches the px sizes of the stylesheet
    f.setWeight(weight)
    f.setHintingPreference(QFont.HintingPreference.PreferNoHinting)
    return f


def mono_family() -> str:
    fams = set(QFontDatabase.families())
    return next((f for f in MONO_FAMILIES if f in fams), "monospace")


def ui_family() -> str:
    fams = set(QFontDatabase.families())
    return next((f for f in FONT_FAMILIES if f in fams), "sans-serif")


# --------------------------------------------------------------------------- icons (inline SVG)

_ICONS = {
    "search": '<circle cx="11" cy="11" r="7"/><path d="M20 20l-3.5-3.5"/>',
    "upload": '<path d="M12 15V4"/><path d="M7 9l5-5 5 5"/><path d="M4 15v3a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2v-3"/>',
    "download": '<path d="M12 4v11"/><path d="M7 10l5 5 5-5"/><path d="M4 15v3a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2v-3"/>',
    "file": '<path d="M14 3H7a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h10a2 2 0 0 0 2-2V8z"/><path d="M14 3v5h5"/>',
    "x": '<path d="M6 6l12 12"/><path d="M18 6L6 18"/>',
    "play": '<path d="M8 5.5v13l10.5-6.5z" fill="currentColor"/>',
    "alert": '<path d="M12 3.5L2.5 20h19z"/><path d="M12 10v4.5"/><path d="M12 17.5v.01"/>',
    "info": '<circle cx="12" cy="12" r="9"/><path d="M12 11v5.5"/><path d="M12 7.5v.01"/>',
    "external": '<path d="M14 4h6v6"/><path d="M20 4l-9 9"/><path d="M18 14v5a1 1 0 0 1-1 1H5a1 1 0 0 1-1-1V7a1 1 0 0 1 1-1h5"/>',
    "copy": '<rect x="8" y="8" width="12" height="12" rx="2"/><path d="M16 8V5a1 1 0 0 0-1-1H5a1 1 0 0 0-1 1v10a1 1 0 0 0 1 1h3"/>',
    "sliders": '<path d="M4 7h10"/><path d="M18 7h2"/><circle cx="16" cy="7" r="2"/><path d="M4 17h4"/><path d="M12 17h8"/><circle cx="10" cy="17" r="2"/>',
    "chevron-down": '<path d="M6 9l6 6 6-6"/>',
    "chevron-right": '<path d="M9 6l6 6-6 6"/>',
    "help": '<circle cx="12" cy="12" r="9"/><path d="M9.5 9.5a2.5 2.5 0 1 1 3.5 2.3c-.6.3-1 .9-1 1.6v.6"/><path d="M12 17.3v.01"/>',
    "sparkle": '<path d="M12 3v4M12 17v4M3 12h4M17 12h4M6 6l2.5 2.5M15.5 15.5L18 18M6 18l2.5-2.5M15.5 8.5L18 6"/>',
    "reset": '<path d="M4 12a8 8 0 1 0 2.4-5.7"/><path d="M4 4v4h4"/>',
    "dna": '<path d="M7 3c0 6 10 6 10 12s-10 3-10 6"/><path d="M17 3c0 6-10 6-10 12s10 3 10 6"/><path d="M8.5 7h7M8.5 17h7M10 11.5h4"/>',
}


def icon(name: str, color: str = INK_2, size: int = 18, stroke: float = 1.8) -> QIcon:
    body = _ICONS[name].replace("currentColor", color)
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="{color}" '
           f'stroke-width="{stroke}" stroke-linecap="round" stroke-linejoin="round">{body}</svg>')
    renderer = QSvgRenderer(QByteArray(svg.encode()))
    ic = QIcon()
    for scale in (1, 2):
        pm = QPixmap(size * scale, size * scale)
        pm.fill(Qt.GlobalColor.transparent)
        p = QPainter(pm)
        renderer.render(p, QRectF(0, 0, size * scale, size * scale))
        p.end()
        pm.setDevicePixelRatio(scale)
        ic.addPixmap(pm)
    return ic


def app_icon() -> QIcon:
    """Rounded petrol tile with a white double helix."""
    ic = QIcon()
    for s in (16, 32, 64, 128, 256, 512):
        svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">
          <defs><linearGradient id="g" x1="0" y1="0" x2="1" y2="1">
            <stop offset="0" stop-color="{PETROL_600}"/><stop offset="1" stop-color="{PETROL_900}"/></linearGradient></defs>
          <rect x="2" y="2" width="60" height="60" rx="14" fill="url(#g)"/>
          <g fill="none" stroke="#FFFFFF" stroke-width="3.6" stroke-linecap="round">
            <path d="M22 12c0 14 20 14 20 20s-20 6-20 20"/>
            <path d="M42 12c0 14-20 14-20 20s20 6 20 20" stroke="{MINERAL_200}"/>
            <path d="M25 19h14M25 45h14" stroke-width="2.6"/>
          </g>
          <circle cx="32" cy="32" r="4.2" fill="{ACCENT}" stroke="#FFFFFF" stroke-width="2"/>
        </svg>'''
        r = QSvgRenderer(QByteArray(svg.encode()))
        pm = QPixmap(s, s)
        pm.fill(Qt.GlobalColor.transparent)
        p = QPainter(pm)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        r.render(p)
        p.end()
        ic.addPixmap(pm)
    return ic


# --------------------------------------------------------------------------- palette + stylesheet

def palette() -> QPalette:
    pal = QPalette()
    roles = {
        QPalette.ColorRole.Window: MINERAL_50, QPalette.ColorRole.WindowText: INK,
        QPalette.ColorRole.Base: WHITE, QPalette.ColorRole.AlternateBase: "#F7FAF9",
        QPalette.ColorRole.Text: INK, QPalette.ColorRole.Button: WHITE, QPalette.ColorRole.ButtonText: INK,
        QPalette.ColorRole.Highlight: MINERAL_200, QPalette.ColorRole.HighlightedText: INK,
        QPalette.ColorRole.ToolTipBase: PETROL_900, QPalette.ColorRole.ToolTipText: WHITE,
        QPalette.ColorRole.PlaceholderText: MUTED, QPalette.ColorRole.Link: PETROL_700,
        QPalette.ColorRole.BrightText: WHITE, QPalette.ColorRole.Mid: BORDER, QPalette.ColorRole.Dark: BORDER_STRONG,
    }
    for role, c in roles.items():
        pal.setColor(role, QColor(c))
    pal.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, QColor(MUTED))
    pal.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, QColor(MUTED))
    return pal


def stylesheet() -> str:
    fam = ", ".join(f'"{f}"' for f in FONT_FAMILIES)
    mono = ", ".join(f'"{f}"' for f in MONO_FAMILIES)
    return f"""
* {{ font-family: {fam}; }}
QWidget {{ color: {INK}; font-size: 13px; }}
QMainWindow, #Main {{ background: {MINERAL_50}; }}
QToolTip {{ background: {PETROL_900}; color: {WHITE}; border: 0; padding: 6px 8px; border-radius: 6px; }}

/* ---------- sidebar ---------- */
#Sidebar {{ background: {PETROL_900}; }}
#Sidebar QLabel {{ color: {MINERAL_200}; background: transparent; }}
#Sidebar QLabel#Brand {{ color: {WHITE}; font-size: 19px; font-weight: 500; letter-spacing: .2px; }}
#Sidebar QLabel#BrandSub {{ color: {MINERAL_400}; font-size: 11px; }}
#Sidebar QLabel#Section {{ color: {MINERAL_400}; font-size: 10.5px; font-weight: 600; letter-spacing: 1.4px; }}
#Sidebar QLabel#Hint {{ color: {MINERAL_400}; font-size: 11px; }}
#Sidebar QLabel#Stage {{ color: {MINERAL_200}; font-size: 11.5px; }}
#Sidebar QPlainTextEdit, #Sidebar QLineEdit, #Sidebar QComboBox, #Sidebar QSpinBox {{
    background: {PETROL_800}; color: {WHITE}; border: 1px solid {PETROL_600}; border-radius: 8px;
    selection-background-color: {PETROL_500}; selection-color: {WHITE};
}}
#Sidebar QPlainTextEdit {{ padding: 8px 6px; font-family: {mono}; font-size: 12.5px; }}
#Sidebar QLineEdit, #Sidebar QComboBox, #Sidebar QSpinBox {{ padding: 6px 10px; min-height: 20px; }}
#Sidebar QPlainTextEdit:focus, #Sidebar QLineEdit:focus, #Sidebar QComboBox:focus, #Sidebar QSpinBox:focus {{
    border: 1px solid {MINERAL_300}; }}
#Sidebar QComboBox::drop-down {{ border: 0; width: 22px; }}
#Sidebar QComboBox QAbstractItemView {{ background: {WHITE}; color: {INK}; selection-background-color: {MINERAL_100};
    selection-color: {INK}; border: 1px solid {BORDER}; outline: 0; padding: 4px; }}
#Sidebar QComboBox:disabled, #Sidebar QLineEdit:disabled {{ color: {MINERAL_400}; border-color: {PETROL_800}; }}
#Sidebar QScrollArea, #Sidebar #SidebarBody {{ background: transparent; border: 0; }}
#Sidebar QPushButton#Ghost {{ background: transparent; color: {MINERAL_200}; border: 1px solid {PETROL_600};
    border-radius: 8px; padding: 6px 12px; }}
#Sidebar QPushButton#Ghost:hover {{ background: {PETROL_800}; color: {WHITE}; }}
#Sidebar QPushButton#Primary {{ background: {WHITE}; color: {PETROL_900}; border: 0; border-radius: 10px;
    padding: 12px 16px; font-size: 14px; font-weight: 600; }}
#Sidebar QPushButton#Primary:hover {{ background: {MINERAL_100}; }}
#Sidebar QPushButton#Primary:pressed {{ background: {MINERAL_200}; }}
#Sidebar QPushButton#Primary:disabled {{ background: {PETROL_700}; color: {MINERAL_300}; }}
#Sidebar QPushButton#Link {{ background: transparent; border: 0; color: {MINERAL_300}; padding: 2px; text-align: left; }}
#Sidebar QPushButton#Link:hover {{ color: {WHITE}; text-decoration: underline; }}
#Sidebar QProgressBar {{ background: {PETROL_800}; border: 0; border-radius: 2px; max-height: 4px; }}
#Sidebar QProgressBar::chunk {{ background: {MINERAL_300}; border-radius: 2px; }}
#Sidebar #Divider {{ background: {PETROL_800}; max-height: 1px; min-height: 1px; }}
#DropZone {{ border: 1.5px dashed {PETROL_500}; border-radius: 10px; background: rgba(255,255,255,0.03); }}
#DropZone[hover="true"] {{ border-color: {MINERAL_300}; background: rgba(255,255,255,0.08); }}
#FileChip {{ background: {PETROL_800}; border-radius: 8px; }}
#FileChip QLabel {{ color: {WHITE}; }}
#FileChip QToolButton {{ border: 0; background: transparent; border-radius: 4px; }}
#FileChip QToolButton:hover {{ background: {PETROL_600}; }}

/* ---------- main area ---------- */
#Header {{ background: {WHITE}; border-bottom: 1px solid {BORDER}; }}
#Header QLabel#Title {{ font-size: 20px; font-weight: 500; color: {INK}; }}
#Header QLabel#Subtitle {{ font-size: 12px; color: {MUTED}; }}
QLabel#H1 {{ font-size: 30px; font-weight: 300; color: {PETROL_900}; }}
QLabel#H2 {{ font-size: 17px; font-weight: 500; color: {INK}; }}
QLabel#H3 {{ font-size: 11px; font-weight: 600; color: {MUTED}; letter-spacing: 1.2px; }}
QLabel#Lead {{ font-size: 14px; color: {INK_2}; }}
QLabel#Muted {{ color: {MUTED}; font-size: 12px; }}
QLabel#Body {{ color: {INK_2}; font-size: 12.5px; }}
QLabel#Mono {{ font-family: {mono}; font-size: 12px; color: {INK_2}; }}
QLabel#Chip {{ background: {MINERAL_100}; color: {PETROL_800}; border-radius: 10px; padding: 2px 9px;
    font-size: 11.5px; font-weight: 500; }}
#Card {{ background: {WHITE}; border: 1px solid {BORDER}; border-radius: 12px; }}
#Kpi {{ background: {WHITE}; border: 1px solid {BORDER}; border-radius: 12px; }}
#Kpi QLabel#KpiValue {{ font-size: 26px; font-weight: 300; color: {PETROL_900}; }}
#Kpi QLabel#KpiLabel {{ font-size: 11.5px; color: {MUTED}; }}
#Banner {{ background: #FBF4E4; border: 1px solid #EFDDB0; border-radius: 10px; }}
#Banner QLabel {{ color: #6B5318; }}
QPushButton {{ background: {WHITE}; border: 1px solid {BORDER_STRONG}; border-radius: 8px; padding: 6px 12px;
    color: {INK}; }}
QPushButton:hover {{ background: {MINERAL_50}; border-color: {MINERAL_300}; }}
QPushButton:pressed {{ background: {MINERAL_100}; }}
QPushButton#Accent {{ background: {PETROL_700}; color: {WHITE}; border: 0; font-weight: 500; padding: 7px 14px; }}
QPushButton#Accent:hover {{ background: {PETROL_600}; }}
QPushButton#Accent::menu-indicator {{ image: none; width: 0; }}
QPushButton#Segment {{ background: {WHITE}; border: 1px solid {BORDER}; border-radius: 14px; padding: 5px 12px;
    min-height: 18px;
    color: {INK_2}; }}
QPushButton#Segment:checked {{ background: {MINERAL_100}; border-color: {MINERAL_300}; color: {INK}; }}
QPushButton#Segment:hover {{ border-color: {MINERAL_300}; }}
QPushButton#ExampleChip {{ background: {WHITE}; border: 1px solid {BORDER}; border-radius: 15px; padding: 6px 14px;
    min-height: 18px;
    font-family: {mono}; font-size: 12px; color: {PETROL_800}; }}
QPushButton#ExampleChip:hover {{ background: {MINERAL_100}; border-color: {MINERAL_300}; }}
QPushButton#LinkButton {{ background: transparent; border: 0; color: {PETROL_700}; padding: 2px 4px; font-weight: 500; }}
QPushButton#LinkButton:hover {{ color: {PETROL_950}; text-decoration: underline; }}
QLineEdit#Search {{ background: {WHITE}; border: 1px solid {BORDER}; border-radius: 15px; padding: 5px 12px;
    min-width: 220px; }}
QLineEdit#Search:focus {{ border-color: {PETROL_500}; }}

QTableView {{ background: {WHITE}; border: 0; gridline-color: transparent; outline: 0;
    selection-background-color: {MINERAL_100}; selection-color: {INK}; alternate-background-color: #F8FBFA; }}
QTableView::item {{ padding: 0 8px; border-bottom: 1px solid #EEF3F1; }}
QTableView::item:selected {{ background: {MINERAL_100}; color: {INK}; }}
QHeaderView {{ background: {WHITE}; border: 0; }}
QHeaderView::section {{ background: {WHITE}; color: {MUTED}; border: 0; border-bottom: 1px solid {BORDER};
    padding: 8px 8px; font-size: 11px; font-weight: 600; letter-spacing: .6px; }}
QHeaderView::section:hover {{ color: {INK}; }}
QTableCornerButton::section {{ background: {WHITE}; border: 0; }}

QTabWidget::pane {{ border: 0; background: {WHITE}; border-top: 1px solid {BORDER}; }}
QTabBar {{ background: transparent; }}
QTabBar::tab {{ background: transparent; color: {MUTED}; padding: 10px 14px; margin-right: 4px; border: 0;
    border-bottom: 2px solid transparent; font-weight: 500; }}
QTabBar::tab:selected {{ color: {PETROL_900}; border-bottom: 2px solid {PETROL_700}; }}
QTabBar::tab:hover {{ color: {INK}; }}

QScrollArea {{ background: transparent; border: 0; }}
QScrollBar:vertical {{ background: transparent; width: 10px; margin: 2px; }}
QScrollBar::handle:vertical {{ background: {MINERAL_200}; border-radius: 3px; min-height: 30px; }}
QScrollBar::handle:vertical:hover {{ background: {MINERAL_300}; }}
QScrollBar:horizontal {{ background: transparent; height: 10px; margin: 2px; }}
QScrollBar::handle:horizontal {{ background: {MINERAL_200}; border-radius: 3px; min-width: 30px; }}
QScrollBar::add-line, QScrollBar::sub-line {{ width: 0; height: 0; }}
QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}
#Sidebar QScrollBar::handle:vertical {{ background: {PETROL_700}; }}
QSplitter::handle {{ background: transparent; }}
QSplitter::handle:vertical {{ height: 10px; }}
QSplitter::handle:horizontal {{ width: 10px; }}
QMenu {{ background: {WHITE}; border: 1px solid {BORDER}; border-radius: 8px; padding: 6px; }}
QMenu::item {{ padding: 7px 22px 7px 12px; border-radius: 6px; }}
QMenu::item:selected {{ background: {MINERAL_100}; color: {INK}; }}
QMenu::separator {{ height: 1px; background: {BORDER}; margin: 4px 6px; }}
QMessageBox, QDialog {{ background: {WHITE}; }}
QTextBrowser {{ background: {WHITE}; border: 0; }}
"""
