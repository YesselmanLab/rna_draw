"""Dark / light theme tokens + QSS builder for the desktop editor.

One palette-token dict per theme (mirrored from the loved editor mockup) feeds
a single `build_qss` template, so the app chrome for both themes is generated
from the same source and only the colors differ. `MainWindow.set_theme` swaps
the applied stylesheet live and pushes the matching canvas/grid colors into the
`RnaGraphicsView`; nothing here touches editor behavior or geometry.

The accent teal is kept consistent with the in-scene selection teal so a
selection reads the same in the panels, the toolbar, and on the canvas.
"""

from __future__ import annotations

# Palette tokens, straight from the mockup (editor_plan.html).
DARK = {
    "bg": "#0e1211",
    "panel": "#161b19",
    "sunk": "#101413",
    "ink": "#e8ecea",
    "muted": "#9aa6a1",
    "faint": "#6b756f",
    "line": "#232a27",
    "accent": "#2bd4c2",
    "accent_ink": "#7ff0e4",
    "accent_soft": "#123330",
    "good": "#4ad07f",
    "warn": "#e5674f",
    "grid": "#1b211f",
    "canvas": "#0a0d0c",  # near-black canvas from the mockup
    "hover": "#1f2623",
    "field": "#101413",
    "sel_fill": "rgba(43, 212, 194, 0.16)",
}

LIGHT = {
    "bg": "#eef1f0",
    "panel": "#ffffff",
    "sunk": "#f4f6f5",
    "ink": "#1a1f1d",
    "muted": "#5c6763",
    "faint": "#8a938f",
    "line": "#dfe4e2",
    "accent": "#0f9d8f",
    "accent_ink": "#0a6b61",
    "accent_soft": "#d5efec",
    "good": "#1f9d55",
    "warn": "#c23b2c",
    "grid": "#e4e8e6",
    "canvas": "#ffffff",
    "hover": "#e7ecea",
    "field": "#ffffff",
    "sel_fill": "rgba(15, 157, 143, 0.14)",
}

THEMES = {"dark": DARK, "light": LIGHT}

# A teal used for the on-canvas selection box + rotate handle, kept aligned
# with each theme's accent so selection reads consistently everywhere.
SELECT_TEAL = "#1fb6a6"

_MONO = '"SF Mono","JetBrains Mono",Menlo,Consolas,"DejaVu Sans Mono",monospace'


def build_qss(t: dict) -> str:
    """Render the full application stylesheet from a token dict `t`."""
    return f"""
/* ==== base ==== */
QMainWindow, QWidget {{
    background: {t['bg']};
    color: {t['ink']};
    font-size: 13px;
}}
QToolTip {{
    background: {t['panel']};
    color: {t['ink']};
    border: 1px solid {t['line']};
    padding: 4px 7px;
    border-radius: 5px;
}}

/* ==== menu bar ==== */
QMenuBar {{
    background: {t['sunk']};
    color: {t['ink']};
    border-bottom: 1px solid {t['line']};
    padding: 2px 4px;
}}
QMenuBar::item {{
    background: transparent;
    padding: 4px 10px;
    border-radius: 6px;
}}
QMenuBar::item:selected {{ background: {t['hover']}; }}
QMenu {{
    background: {t['panel']};
    color: {t['ink']};
    border: 1px solid {t['line']};
    border-radius: 8px;
    padding: 5px;
}}
QMenu::item {{ padding: 5px 22px 5px 14px; border-radius: 5px; }}
QMenu::item:selected {{ background: {t['accent_soft']}; color: {t['accent_ink']}; }}
QMenu::separator {{ height: 1px; background: {t['line']}; margin: 5px 8px; }}

/* ==== toolbars ==== */
QToolBar {{
    background: {t['sunk']};
    border-bottom: 1px solid {t['line']};
    padding: 6px 10px;
    spacing: 6px;
}}
QToolBar::separator {{
    width: 1px;
    background: {t['line']};
    margin: 2px 6px;
}}
QToolBar QToolButton {{
    padding: 5px 11px;
    border-radius: 7px;
    border: 1px solid {t['line']};
    background: {t['panel']};
    color: {t['ink']};
}}
QToolBar QToolButton:hover {{
    background: {t['hover']};
    border: 1px solid {t['accent']};
}}
QToolBar QToolButton:checked {{
    background: {t['accent_soft']};
    border: 1px solid {t['accent']};
    color: {t['accent_ink']};
}}
QToolBar QToolButton:disabled {{ color: {t['faint']}; }}
QToolBar QLabel {{ color: {t['muted']}; padding-left: 4px; }}
/* primary teal-filled action button (e.g. Checkpoint / Snapshot) */
QToolButton#primaryButton, QToolBar QToolButton#primaryButton {{
    background: {t['accent']};
    border: 1px solid {t['accent']};
    color: #08201d;
    font-weight: 600;
}}
QToolButton#primaryButton:hover, QToolBar QToolButton#primaryButton:hover {{
    background: {t['accent_ink']};
    border: 1px solid {t['accent_ink']};
}}
/* collapsible-section headers: muted, quiet, flat */
QToolButton#sectionHeader {{
    border: none;
    background: transparent;
    color: {t['muted']};
    padding: 5px 2px;
    text-align: left;
}}
QToolButton#sectionHeader:hover {{ color: {t['ink']}; }}

/* ==== docks / panels ==== */
QDockWidget {{
    color: {t['muted']};
    titlebar-close-icon: none;
    titlebar-normal-icon: none;
}}
QDockWidget::title {{
    background: {t['sunk']};
    padding: 7px 12px;
    border-bottom: 1px solid {t['line']};
}}
QScrollArea {{ background: {t['bg']}; border: none; }}
QScrollBar:vertical {{ background: transparent; width: 11px; margin: 2px; }}
QScrollBar::handle:vertical {{
    background: {t['line']};
    border-radius: 5px;
    min-height: 30px;
}}
QScrollBar::handle:vertical:hover {{ background: {t['faint']}; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; }}
QScrollBar:horizontal {{ background: transparent; height: 11px; margin: 2px; }}
QScrollBar::handle:horizontal {{
    background: {t['line']};
    border-radius: 5px;
    min-width: 30px;
}}

/* ==== group boxes ==== */
QGroupBox {{
    background: {t['panel']};
    border: 1px solid {t['line']};
    border-radius: 9px;
    margin-top: 14px;
    padding: 10px 9px 8px 9px;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 11px;
    padding: 0 5px;
    color: {t['faint']};
    font-weight: 600;
}}
QGroupBox QLabel {{ color: {t['muted']}; }}
QLabel {{ color: {t['ink']}; background: transparent; }}

/* ==== inputs ==== */
QLineEdit {{
    background: {t['field']};
    border: 1px solid {t['line']};
    border-radius: 7px;
    padding: 5px 9px;
    color: {t['ink']};
    selection-background-color: {t['accent']};
    selection-color: #08201d;
}}
QLineEdit:focus {{ border: 1px solid {t['accent']}; }}

QDoubleSpinBox, QSpinBox, QComboBox {{
    background: {t['field']};
    border: 1px solid {t['line']};
    border-radius: 7px;
    padding: 4px 8px;
    color: {t['ink']};
    min-height: 20px;
}}
QDoubleSpinBox:focus, QSpinBox:focus, QComboBox:focus {{ border: 1px solid {t['accent']}; }}
QComboBox::drop-down {{ border: none; width: 20px; }}
QComboBox::down-arrow {{
    width: 0; height: 0;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-top: 5px solid {t['muted']};
    margin-right: 7px;
}}
QComboBox QAbstractItemView {{
    background: {t['panel']};
    color: {t['ink']};
    border: 1px solid {t['line']};
    border-radius: 8px;
    selection-background-color: {t['accent_soft']};
    selection-color: {t['accent_ink']};
    outline: none;
    padding: 3px;
}}
QSpinBox::up-button, QSpinBox::down-button,
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {{
    background: transparent; width: 16px; border: none;
}}

/* ==== checkboxes ==== */
QCheckBox {{ spacing: 7px; color: {t['ink']}; background: transparent; }}
QCheckBox::indicator {{
    width: 16px; height: 16px;
    border: 1px solid {t['line']};
    border-radius: 5px;
    background: {t['field']};
}}
QCheckBox::indicator:hover {{ border: 1px solid {t['accent']}; }}
QCheckBox::indicator:checked {{
    background: {t['accent']};
    border: 1px solid {t['accent']};
}}

/* ==== sliders ==== */
QSlider::groove:horizontal {{
    height: 4px; border-radius: 2px; background: {t['line']};
}}
QSlider::sub-page:horizontal {{ background: {t['accent']}; border-radius: 2px; }}
QSlider::handle:horizontal {{
    background: {t['accent']};
    width: 14px; height: 14px;
    margin: -6px 0; border-radius: 7px;
}}

/* ==== buttons ==== */
QPushButton {{
    padding: 6px 13px;
    border-radius: 7px;
    background: {t['panel']};
    color: {t['ink']};
    border: 1px solid {t['line']};
}}
QPushButton:hover {{ background: {t['hover']}; border: 1px solid {t['accent']}; }}
QPushButton:pressed {{ background: {t['accent_soft']}; }}
QPushButton:disabled {{ color: {t['faint']}; background: {t['sunk']}; }}
QPushButton#primaryButton {{
    background: {t['accent']};
    border: 1px solid {t['accent']};
    color: #08201d;
    font-weight: 600;
}}
QPushButton#primaryButton:hover {{
    background: {t['accent_ink']};
    border: 1px solid {t['accent_ink']};
}}

/* ==== versions list: card-like rows, current outlined in teal ==== */
QListWidget {{
    background: {t['sunk']};
    border: 1px solid {t['line']};
    border-radius: 8px;
    padding: 4px;
    outline: none;
}}
QListWidget::item {{
    background: {t['panel']};
    border: 1px solid {t['line']};
    border-radius: 7px;
    padding: 7px 9px;
    margin: 3px 2px;
    color: {t['muted']};
}}
QListWidget::item:hover {{ border: 1px solid {t['faint']}; }}
QListWidget::item:selected {{
    background: {t['accent_soft']};
    border: 1px solid {t['accent']};
    color: {t['accent_ink']};
}}

/* ==== canvas + status bar ==== */
QGraphicsView {{ background: {t['canvas']}; border: none; }}
QStatusBar {{
    background: {t['sunk']};
    border-top: 1px solid {t['line']};
    color: {t['muted']};
    font-family: {_MONO};
    font-size: 11.5px;
}}
QStatusBar::item {{ border: none; }}
QStatusBar QLabel {{ padding: 2px 10px; color: {t['muted']}; }}
"""


def qss_for(name: str) -> str:
    """Full application stylesheet for theme `name` (falls back to dark)."""
    return build_qss(THEMES.get(name, DARK))


def canvas_colors(name: str) -> tuple[str, str]:
    """`(canvas_base, grid_dot)` hex colors for theme `name`."""
    t = THEMES.get(name, DARK)
    return t["canvas"], t["grid"]


__all__ = ["THEMES", "DARK", "LIGHT", "SELECT_TEAL", "build_qss", "qss_for", "canvas_colors"]
