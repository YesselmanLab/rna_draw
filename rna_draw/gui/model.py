"""`EditorModel`: the frontend-agnostic RNA-helix editing controller.

This holds the authoritative geometry (`_x`, `_y`, the structure tree) and
performs the only two rigid edits the editors expose -- rotate and translate
a selected helix about its junction pivot. It imports NO GUI toolkit, so the
editing logic is unit-testable headlessly (see `tests/test_desktop_mvp.py`);
the Qt view in `rna_draw.gui.desktop` and the Jupyter widget in
`rna_draw.gui.editor_widget` are thin shells over this same logic.

Never-silent contract (identical to the Jupyter widget): every committed move
is re-checked in-process by the REAL native overlap kernel
(`rna_draw.overlap_native.check_overlaps_native`, the exact twin of the frozen
arbiter) at the layout's own SCALED `node_r` -- the same radius the pipeline
gated at, NOT a fresh default. Overlaps are flagged (never presented as
clean); offending nucleotide indices come from the definitional checker's
witnesses so a view can tint them red. Coordinates are ENGINE space (y-up); a
view flips y for screen.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass

from rna_draw.document import Document
from rna_draw.document_render import SourceInputs, document_from_layout, resolve_preset
from rna_draw.gui.hinge import redistribute_loop
from rna_draw.layout import LayoutResult
from rna_draw.layout.base import empty_report, params_at_node_r
from rna_draw.layout.pipeline import layout_guaranteed, resolve_engine
from rna_draw.layout.postpass import (
    rotate_range,
    translate_range,
    witness_nucleotides,
)
from rna_draw.layout.structure_tree import Loop, build_structure_tree, loop_center
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.overlap_native import check_overlaps_native
from rna_draw.parameters import RenderType
from rna_draw.render_rna import get_pairmap_from_secstruct
from rna_draw.scene import _rgb_to_hex, build_scene
from rna_draw.style import StylePreset, default_preset
from rna_draw.validate import never_silent_gate

# The target render radius; the pipeline scales down from here when gating.
_TARGET_NODE_R = 10.0

# Rotation clamp: a selected helix may swing at most this far (radians) from
# the direction it was LAID OUT / SELECTED at, in either sense. 120 degrees is
# generous enough for real re-arrangement while making it impossible to flip a
# helix around into an unphysical pose (a full 180 degree flip is barred). The
# per-selection bound is the TIGHTER of this and the neighbouring-anchor limit
# (a helix cannot rotate past a sibling helix / the loop's closing pair).
MAX_HELIX_ROTATION = math.radians(120.0)

# Angular breathing room kept between a rotating helix's base and a neighbour
# anchor on the parent loop, so the clamp stops just short of a collision.
_NEIGHBOR_MARGIN = math.radians(8.0)


def _normalize_pi(angle: float) -> float:
    """Wrap an angle to the half-open interval ``(-pi, pi]``."""
    wrapped = math.fmod(angle, 2.0 * math.pi)
    if wrapped <= -math.pi:
        wrapped += 2.0 * math.pi
    elif wrapped > math.pi:
        wrapped -= 2.0 * math.pi
    return wrapped

# The demo structure shown on launch so the window is populated immediately: a
# three-way junction (three hairpins off a central multiloop).
DEMO_SS = "((((...((((....))))...((((....))))...))))"

# A valid A/C/G/U sequence exactly len(DEMO_SS) long, complementary across the
# three stems, so the demo shows real letters (and any res_type/paired scheme
# has something to color) the moment the window opens.
DEMO_SEQ = "GGGGAAAGGGGAAAACCCCAAAGGGGAAAACCCCAAACCCC"

# Render-type selector strings (panel + preset.extra["render_type"]) mapped to
# the colorer's `RenderType` enum. "none" means "no scheme -> default fill".
RENDER_TYPES = {"none": None, "res_type": RenderType.RES_TYPE, "paired": RenderType.PAIRED}

# Base (100%) reference values for the percentage sliders in the style panel.
# 100% reproduces the historical default look; `view_style` multiplies each base
# by its stored percent so the scene values below are byte-identical to before
# when every slider sits at 100%.
_BASE_EDGE_WIDTH = 1.0
_BASE_PAIR_WIDTH = 1.6
_BASE_CROSSING_WIDTH = 1.6
_BASE_ROUTED_WIDTH = 1.6
_BASE_BACKBONE_WIDTH = 2.0

# Non-geometry visual knobs that have no home in the canonical `StylePreset`
# fields. They are carried in `preset.extra["view"]` so they round-trip through
# `.rnastyle.json` save/load (via the preset's `extra` merge) without inventing
# a new file format. Sizes are stored as PERCENTAGES of the base defaults above
# (100% = the historical look); `view_style` converts them to concrete widths /
# scales. `sphere_scale` is a LIVE DISPLAY scale on the drawn disk radius only --
# it never touches the layout node_r the overlap checker gates on.
DEFAULT_VIEW = {
    "background": "#ffffff",
    "backbone_color": "#8a8f98",
    "nt_edge_color": "#2b2f36",
    "letter_color": "#101418",
    "routed_style": "dash",  # solid | dash | dot
    "render_type": "none",  # none | res_type | paired
    # Overall render STYLE (extensible string): "spheres" = the classic
    # disk/backbone/pair look; "letters" = the RFview-style colored-letters mode
    # (letters, gapped backbone, typed pair symbols). Display-only; add more
    # style strings here as they are implemented.
    "render_style": "spheres",  # spheres | letters
    "data_palette": "viridis",  # SHAPE/data colormap name (persisted; see report)
    # Percentage size sliders (100% = base default). Display-only, no relayout.
    "sphere_pct": 100,  # drawn disk radius as % of layout node_r
    "letter_pct": 100,  # letter font size as % (maps to letter_scale)
    "edge_pct": 100,  # nucleotide outline width
    "pair_pct": 100,  # base-pair line width
    "backbone_pct": 100,  # backbone line width
    "routed_pct": 100,  # PK routed line width
    "show_spheres": True,  # draw the nucleotide disks at all
}


def _load_overlap_note(report, dirty) -> str:
    """Human-readable never-silent note for an overlapping loaded layout."""
    parts = []
    if not report.passed:
        parts.append(f"{report.num_overlaps} disk/backbone/pair overlap(s)")
    if dirty:
        parts.append(f"{len(dirty)} routed line(s) crossing the layout")
    detail = "; ".join(parts) if parts else "overlaps"
    return f"loaded layout is NOT clean ({detail}) -- shown flagged, not silently clean"


def _bond_type(a: str, b: str) -> str:
    """Classify a base pair by its two bases into a drawn ``bond`` type.

    Normalizes ``T -> U`` and uppercases; returns one of ``"gc"`` ({G,C}),
    ``"au"`` ({A,U}), ``"gu"`` ({G,U}), or ``"other"`` (any non-canonical
    combination, or a blank/unknown base at either position). Purely a
    display classification -- it drives which pair SYMBOL the letters render
    style draws (=, -, wobble, LW placeholder); it never affects geometry.
    """
    a = (a or "").upper().replace("T", "U").strip()
    b = (b or "").upper().replace("T", "U").strip()
    if a not in ("A", "C", "G", "U") or b not in ("A", "C", "G", "U"):
        return "other"
    pair = {a, b}
    if pair == {"G", "C"}:
        return "gc"
    if pair == {"A", "U"}:
        return "au"
    if pair == {"G", "U"}:
        return "gu"
    return "other"


def _resolve_rgb(value, palette: dict) -> list[float]:
    """Resolve a color spec to an ``[r, g, b]`` 0..1 triple.

    Accepts a ``#rrggbb`` hex string, a single-letter palette key (looked up
    in ``palette``), or an existing ``[r, g, b]`` triple.
    """
    if isinstance(value, str):
        if value.startswith("#") and len(value) >= 7:
            v = value.lstrip("#")
            return [int(v[i : i + 2], 16) / 255.0 for i in (0, 2, 4)]
        if value in palette:
            return list(palette[value])
    elif isinstance(value, (list, tuple)) and len(value) >= 3:
        return [float(c) for c in value[:3]]
    return [0.588, 0.588, 0.588]


def _resolve_color(value, palette: dict) -> str:
    """Resolve a color spec (hex, palette key, or rgb triple) to ``#rrggbb``."""
    return _rgb_to_hex(_resolve_rgb(value, palette))


def view_style(preset: StylePreset) -> dict:
    """Resolve a `StylePreset` into the flat, Qt-free style dict the scene carries.

    Every color is a ``#rrggbb`` string and every size a float, so the value
    is pure data (no toolkit types) and lands directly in `Scene`'s ``style``
    block. Connector colors are resolved through the preset palette (a
    single-letter key) or taken verbatim if the panel wrote a hex string.
    """
    view = dict(DEFAULT_VIEW)
    view.update(preset.extra.get("view", {}))
    palette = preset.palette
    ld = preset.layout_defaults
    cc = preset.connector_colors

    def _pct(key: str) -> float:
        """Stored percentage for `key` as a 0..1+ multiplier (100 -> 1.0)."""
        return float(view.get(key, 100)) / 100.0

    return {
        "background": view["background"],
        "backbone_color": view["backbone_color"],
        "backbone_width": _BASE_BACKBONE_WIDTH * _pct("backbone_pct"),
        "default_fill": _resolve_color(preset.default_color, palette),
        "nt_edge_color": view["nt_edge_color"],
        "nt_edge_width": _BASE_EDGE_WIDTH * _pct("edge_pct"),
        "show_letters": bool(ld.render_in_letters),
        "show_spheres": bool(view.get("show_spheres", True)),
        # Live DISPLAY scale on the drawn disk radius only (not a relayout): the
        # scene draws node_r * sphere_scale, decoupled from the layout node_r the
        # overlap checker gates on. A big visual sphere may look overlapping while
        # the layout is clean -- that is cosmetic and correct.
        "sphere_scale": _pct("sphere_pct"),
        "letter_scale": _pct("letter_pct"),
        "letter_color": view["letter_color"],
        "pair_color": _resolve_color(cc.get("nested_pair", "e"), palette),
        "pair_width": _BASE_PAIR_WIDTH * _pct("pair_pct"),
        "crossing_color": _resolve_color(cc.get("pk_connector", "o"), palette),
        "crossing_width": _BASE_CROSSING_WIDTH,
        "routed_color": _resolve_color(cc.get("pk_line", "r"), palette),
        "routed_width": _BASE_ROUTED_WIDTH * _pct("routed_pct"),
        "routed_style": view["routed_style"],
        # Overall render style ("spheres" | "letters"); a plain string so more
        # styles can be added without touching the value model.
        "render_style": view.get("render_style", "spheres"),
    }


@dataclass
class Result:
    """Outcome of a committed rigid move.

    Args:
        flagged: True iff the native checker found >=1 overlap in the
            committed geometry (never presented as clean when True).
        overlaps: Offending nucleotide indices (for red tinting), empty
            when `flagged` is False.
    """

    flagged: bool
    overlaps: list[int]


@dataclass
class Selection:
    """A resolved click selection at a chosen granularity.

    Args:
        kind: One of ``"residue"`` (a single nucleotide), ``"helix"`` (the
            enclosing stem slice, the only rotatable selection), or
            ``"motif"`` (the enclosing loop / structural element).
        indices: The selected nucleotide indices (sorted).
        pivot: Rotation pivot in engine space; set only for a ``"helix"``
            selection (``None`` otherwise -- residues and motifs do not rotate).
        start: First index of the contiguous helix slice (``"helix"`` only).
        end: Last index of the contiguous helix slice (``"helix"`` only).
    """

    kind: str
    indices: list[int]
    pivot: tuple[float, float] | None = None
    start: int | None = None
    end: int | None = None


class EditorModel:
    """Authoritative click-to-select, rotate/translate RNA-helix editor.

    Pure logic, no Qt import. Holds engine-space coordinates plus the
    topological structure tree and mutates them via exactly-rigid moves,
    re-checking every commit with the native kernel.
    """

    def __init__(self) -> None:
        self._x: list[float] = []
        self._y: list[float] = []
        self._node_r: float = _TARGET_NODE_R
        self._pair_map: list[int] = []
        self._ss: str | None = None
        self._seq: str | None = None
        self._colors = None
        # The authoritative styling for scene rendering (palette, connector
        # colors, letters, geometry defaults). Mutated by `set_style`.
        self._preset: StylePreset = default_preset()
        self._crossing_pairs: list = []
        self._crossing_lines: list = []
        self._tree = None
        self._engine_name: str = ""
        # Current selection state.
        self._sel_start: int | None = None
        self._sel_end: int | None = None
        self._pivot: tuple[float, float] | None = None
        # Parent loop of the current selection, redistributed on rotation so
        # the loop stays a clean circle (VARNA hinge; see rna_draw.gui.hinge).
        self._sel_loop: Loop | None = None
        # Granular selection: the highlighted nucleotide set and its kind
        # ("residue" | "helix" | "motif"), independent of the rotatable helix
        # slice above. Drives the teal highlight in the scene.
        self._sel_indices: set[int] = set()
        self._sel_kind: str | None = None
        # Rotation clamp reference (set when a HELIX is selected): the helix
        # base direction at selection, and the allowed cumulative deviation
        # band [lo, hi] (radians) -- the tighter of +/-MAX_HELIX_ROTATION and
        # the neighbouring-anchor limits. See `_clamp_rotation`.
        self._sel_ref_angle: float | None = None
        self._sel_dev_lo: float = 0.0
        self._sel_dev_hi: float = 0.0
        # Cached checker verdict for the current geometry.
        self._flagged: bool = False
        self._overlaps: list[int] = []
        # Human-readable note from the last `from_document`/`_adopt_document`
        # load whose stored coords failed the live never-silent gate (empty
        # when the load was clean); a view surfaces it non-blockingly.
        self._load_note: str = ""
        # Residue-numbering (VARNA-style "10, 20, 30 ..." labels). Off by
        # default; the desktop toggle turns it on. Purely visual -- the
        # numbers are derived from geometry in `scene()`, never stored.
        self._show_numbers: bool = False
        self._number_interval: int = 10
        # Per-region VISUAL styling (no geometry, no overlap risk): explicit
        # per-nucleotide color overrides (highest fill precedence, layered on
        # top of the colorer output) and persistent highlight halos. Both
        # round-trip through the Document's `extra` band.
        self._color_overrides: dict[int, str] = {}
        self._highlights: list[dict] = []
        # Undo/redo history: each entry is a full editable-state snapshot from
        # `_capture_state`. `push_undo` is called by the frontend ONCE per user
        # gesture (a whole drag = one entry). Bounded so a long session cannot
        # grow the history without limit; the oldest entry is dropped past cap.
        self._undo: list[dict] = []
        self._redo: list[dict] = []
        self._undo_cap: int = 100

    # -- undo / redo (full-state snapshots) ---------------------------------

    def _capture_state(self) -> dict:
        """Snapshot the whole editable model state as an independent copy.

        Deep/independent copies of every mutable field so a later edit cannot
        corrupt a stored snapshot (selection is intentionally NOT captured --
        selection changes are not undoable, only STATE changes are).
        """
        return {
            "x": list(self._x),
            "y": list(self._y),
            "node_r": self._node_r,
            "pair_map": list(self._pair_map),
            "ss": self._ss,
            "seq": self._seq,
            "preset": copy.deepcopy(self._preset),
            "color_overrides": dict(self._color_overrides),
            "highlights": [
                {"indices": list(h["indices"]), "color": h["color"]} for h in self._highlights
            ],
            "crossing_pairs": copy.deepcopy(self._crossing_pairs),
            "crossing_lines": copy.deepcopy(self._crossing_lines),
            "engine_name": self._engine_name,
        }

    def _restore_state(self, state: dict) -> None:
        """Restore a `_capture_state` snapshot, rebuild derived data, re-check.

        Rebuilds the structure tree (pair_map/seq may have changed), recomputes
        the per-residue fills, clears the (non-undoable) selection, and RE-RUNS
        the never-silent native checker so the restored layout's flagged /
        overlap status is honest -- a restored pose is never silently clean.
        """
        self._x = list(state["x"])
        self._y = list(state["y"])
        self._node_r = float(state["node_r"])
        self._pair_map = list(state["pair_map"])
        self._ss = state["ss"]
        self._seq = state["seq"]
        self._preset = copy.deepcopy(state["preset"])
        self._color_overrides = dict(state["color_overrides"])
        self._highlights = [
            {"indices": list(h["indices"]), "color": h["color"]} for h in state["highlights"]
        ]
        self._crossing_pairs = copy.deepcopy(state["crossing_pairs"])
        self._crossing_lines = copy.deepcopy(state["crossing_lines"])
        self._engine_name = state["engine_name"]
        self._clear_selection()
        self._tree = build_structure_tree(self._pair_map)
        self._recompute_colors()
        # NEVER-SILENT: re-validate the restored coords with the real native
        # kernel at the layout's own gate radius; flag reflects reality.
        if len(self._x) >= 2:
            count = check_overlaps_native(
                self._x, self._y, self._pair_map, self._scaled_params()
            )
        else:
            count = 0
        self._flagged = count > 0
        self._overlaps = self._overlap_indices() if self._flagged else []

    def push_undo(self) -> None:
        """Record the current state as an undo point and clear the redo stack.

        The frontend calls this ONCE at the start of each user gesture, BEFORE
        the first mutation (a whole drag coalesces to a single entry because
        the per-frame mutations do not push). Any new edit invalidates the redo
        stack. History is capped; the oldest entry is dropped past the cap.
        """
        self._undo.append(self._capture_state())
        if len(self._undo) > self._undo_cap:
            self._undo.pop(0)
        self._redo.clear()

    def can_undo(self) -> bool:
        """Whether there is at least one state to undo to."""
        return bool(self._undo)

    def can_redo(self) -> bool:
        """Whether there is at least one undone state to redo to."""
        return bool(self._redo)

    def undo(self) -> bool:
        """Restore the previous state; push the current state onto redo.

        Safe no-op (returns False) when the undo stack is empty.
        """
        if not self._undo:
            return False
        self._redo.append(self._capture_state())
        self._restore_state(self._undo.pop())
        return True

    def redo(self) -> bool:
        """Re-apply the most recently undone state; push current onto undo.

        Safe no-op (returns False) when the redo stack is empty.
        """
        if not self._redo:
            return False
        self._undo.append(self._capture_state())
        self._restore_state(self._redo.pop())
        return True

    def reset_history(self) -> None:
        """Discard all undo/redo history (e.g. when a new document is loaded)."""
        self._undo.clear()
        self._redo.clear()

    # -- construction -------------------------------------------------------

    @classmethod
    def from_ss(cls, ss: str, seq: str | None = None, engine: str = "auto") -> "EditorModel":
        """Build a model by running the guaranteed layout on a dot-bracket.

        Args:
            ss: Dot-bracket secondary structure.
            seq: Optional sequence for per-nt labels; must match ``len(ss)``.
            engine: Layout-engine selection string for `resolve_engine`.

        Returns:
            A ready-to-edit `EditorModel`.
        """
        result = layout_guaranteed(
            ss, engine=resolve_engine(engine), params=OverlapParams(node_r=_TARGET_NODE_R)
        )
        pair_map = result.pair_map
        if pair_map is None:
            pair_map = get_pairmap_from_secstruct(ss)

        colors = None
        if seq is not None:
            try:
                from rna_draw.colorer import Colorer

                colors = Colorer().get_rgb_colors(seq, ss, default_color=None)
            except Exception:
                colors = None

        model = cls()
        model._x = [float(v) for v in result.x]
        model._y = [float(v) for v in result.y]
        model._node_r = float(result.node_r)
        model._pair_map = [int(v) for v in pair_map]
        model._ss = ss
        model._seq = seq
        model._colors = colors
        model._crossing_pairs = list(result.crossing_pairs or [])
        model._crossing_lines = list(result.crossing_lines or [])
        model._tree = build_structure_tree(model._pair_map)
        model._engine_name = result.engine_name
        # Show letters out of the box only when a sequence is present (matches
        # the pre-panel behaviour); the panel toggle overrides this.
        model._preset.layout_defaults.render_in_letters = seq is not None
        model._recompute_colors()
        # Honest initial verdict from the real checker.
        model._flagged = bool(result.flagged) or len(model._overlap_indices()) > 0
        model._overlaps = model._overlap_indices() if model._flagged else []
        return model

    # -- document (de)serialization -----------------------------------------

    def to_document(self, name: str | None = None) -> Document:
        """Capture the CURRENT editor state as a `Document`.

        The derived layout band is built from the model's LIVE, possibly
        hand-edited coordinates (`_x`/`_y`, engine space, UNSHIFTED) -- NOT by
        re-running layout -- so an arranged drawing is persisted exactly as the
        user left it. The source band records the structure the layout was
        built from (`_ss`/`_seq`) plus the live style preset (inlined, so it
        round-trips without a registry lookup). The advisory checker verdict is
        an honest snapshot of the current flag; the live never-silent gate is
        re-run on load, so it is never trusted as authoritative.

        Args:
            name: Optional label stored in ``Document.extra["name"]`` (used by
                the editor's version-history list); ``None`` omits it.

        Returns:
            A `Document` whose `derived.layout` mirrors the live coordinates.
        """
        n = len(self._x)
        ss = self._ss if self._ss is not None else "." * n
        seq = self._seq if self._seq is not None else ""
        view = {**DEFAULT_VIEW, **self._preset.extra.get("view", {})}
        source = SourceInputs(
            ss=ss,
            seq=seq,
            render_type=view.get("render_type", "none"),
            style_preset=self._preset.to_dict(),
            engine="auto",
        )
        report = (
            check_overlaps(self._x, self._y, self._pair_map, self._scaled_params())
            if n >= 2
            else empty_report()
        )
        result = LayoutResult(
            x=list(self._x),
            y=list(self._y),
            engine_name=self._engine_name,
            report=report,
            flagged=self._flagged,
            node_r=self._node_r,
            pair_map=list(self._pair_map),
            crossing_pairs=list(self._crossing_pairs),
            crossing_lines=list(self._crossing_lines),
        )
        doc = document_from_layout(source, result)
        if name is not None:
            doc.extra["name"] = name
        # Persist per-region visual styling in the Document's extra band so
        # save/load + version snapshots keep recolored regions and highlights.
        if self._color_overrides:
            doc.extra["color_overrides"] = {
                str(k): v for k, v in sorted(self._color_overrides.items())
            }
        if self._highlights:
            doc.extra["highlights"] = [
                {"indices": list(h["indices"]), "color": h["color"]} for h in self._highlights
            ]
        return doc

    @classmethod
    def from_document(cls, doc: Document) -> "EditorModel":
        """Build a model from a `Document`, re-running the never-silent gate.

        The stored coordinates/style are adopted verbatim, but the live
        `validate.never_silent_gate` ALWAYS re-validates them (the file's
        advisory checker block is not trusted): a clean layout is drawn as-is;
        an overlapping stored layout is adopted flagged (its `load_note` set),
        never presented as silently clean.

        Args:
            doc: The document to restore.

        Returns:
            A ready-to-edit `EditorModel`.
        """
        model = cls()
        model._adopt_document(doc)
        return model

    def _adopt_document(self, doc: Document) -> str:
        """Restore `doc`'s coords/seq/style into this model (never-silent).

        Returns the load note (empty when the stored layout is clean).
        """
        structure = doc.source.structure
        self._preset = resolve_preset(doc.source.style_preset)
        self._ss = structure.ss or None
        self._seq = structure.seq or None
        self._clear_selection()
        self._restore_region_styling(doc)
        layout = doc.derived.layout
        if layout is None:
            # Source-only document: no stored coords, lay out fresh (gated by
            # the guaranteed pipeline exactly as `from_ss` does).
            fresh = EditorModel.from_ss(structure.ss, seq=self._seq)
            self._x, self._y = fresh._x, fresh._y
            self._node_r = fresh._node_r
            self._pair_map = fresh._pair_map
            self._crossing_pairs = fresh._crossing_pairs
            self._crossing_lines = fresh._crossing_lines
            self._tree = fresh._tree
            self._engine_name = fresh._engine_name
            self._flagged = fresh._flagged
            self._overlaps = list(fresh._overlaps)
            self._recompute_colors()
            self._load_note = ""
            return ""

        xs = [float(px) for px, _ in layout.coords]
        ys = [float(py) for _, py in layout.coords]
        self._x, self._y = xs, ys
        self._node_r = float(layout.node_r)
        self._pair_map = [int(v) for v in layout.pair_map]
        self._crossing_pairs = list(layout.crossing_pairs)
        self._crossing_lines = list(layout.crossing_lines)
        self._engine_name = layout.engine_name
        self._tree = build_structure_tree(self._pair_map)
        self._recompute_colors()
        # NEVER-SILENT: re-validate the stored coords live at the model's own
        # gate radius; the file's advisory `checker` block is ignored.
        clean, report, dirty = never_silent_gate(
            xs, ys, self._pair_map, self._crossing_lines, self._node_r, _TARGET_NODE_R
        )
        self._flagged = not clean
        self._overlaps = self._overlap_indices() if not clean else []
        self._load_note = "" if clean else _load_overlap_note(report, dirty)
        return self._load_note

    # -- read-only accessors ------------------------------------------------

    @property
    def engine_name(self) -> str:
        """Name of the engine that produced the current geometry."""
        return self._engine_name

    @property
    def flagged(self) -> bool:
        """Whether the current committed geometry has overlaps."""
        return self._flagged

    @property
    def load_note(self) -> str:
        """Never-silent note from the last document load (``""`` if clean)."""
        return self._load_note

    @property
    def selection(self) -> tuple[int, int] | None:
        """The selected helix slice ``(start, end)`` or ``None``."""
        if self._sel_start is None or self._sel_end is None:
            return None
        return self._sel_start, self._sel_end

    @property
    def pivot(self) -> tuple[float, float] | None:
        """The current rotation pivot in engine space, or ``None``."""
        return self._pivot

    @property
    def sel_indices(self) -> list[int]:
        """The currently highlighted nucleotide indices (sorted)."""
        return sorted(self._sel_indices)

    @property
    def sel_kind(self) -> str | None:
        """The current selection kind (``"residue"``/``"helix"``/``"motif"``) or ``None``."""
        return self._sel_kind

    # -- checker plumbing ---------------------------------------------------

    def _scaled_params(self) -> OverlapParams:
        """Checker params at the layout's actual ``node_r`` (the gate radius)."""
        return params_at_node_r(OverlapParams(node_r=_TARGET_NODE_R), self._node_r)

    def _overlap_indices(self) -> list[int]:
        """Offending nucleotide indices from the definitional checker's witnesses."""
        if len(self._x) < 2:
            return []
        report = check_overlaps(self._x, self._y, self._pair_map, self._scaled_params())
        idx: set[int] = set()
        for w in report.witnesses:
            ka, kb = witness_nucleotides(w, self._pair_map)
            idx.add(ka)
            idx.add(kb)
        return sorted(idx)

    def _centroid(self) -> tuple[float, float]:
        """Fallback pivot: the whole-structure centroid."""
        n = len(self._x)
        if n == 0:
            return 0.0, 0.0
        return sum(self._x) / n, sum(self._y) / n

    # -- selection ----------------------------------------------------------

    def select(self, index: int) -> tuple[int, int, tuple[float, float]] | None:
        """Select the innermost helix containing ``index`` and set its pivot.

        Back-compat entry point (helix granularity); `select_at` is the
        general form.

        Args:
            index: Clicked nucleotide index.

        Returns:
            ``(start, end, (cx, cy))`` -- the contiguous helix slice and its
            parent-loop-center pivot -- or ``None`` if ``index`` is an
            exterior/unpaired nucleotide inside no helix (selection cleared).
        """
        if self._tree is None or not (0 <= index < len(self._pair_map)):
            self._clear_selection()
            return None
        sel = self._select_helix(index)
        if sel is None:
            return None
        return sel.start, sel.end, sel.pivot

    def select_at(self, index: int, granularity: str = "helix") -> Selection | None:
        """Resolve a click to a selection at the requested granularity.

        Args:
            index: Clicked nucleotide index.
            granularity: ``"residue"`` (just this nucleotide), ``"helix"``
                (the enclosing stem slice), or ``"motif"`` (the enclosing
                loop's members, or the helix when clicked on a stem).

        Returns:
            A `Selection`, or ``None`` when nothing resolves (e.g. a helix
            request on an exterior/unpaired nucleotide); selection cleared.
        """
        if self._tree is None or not (0 <= index < len(self._pair_map)):
            self._clear_selection()
            return None
        if granularity == "residue":
            return self._select_residue(index)
        if granularity == "motif":
            return self._select_motif(index)
        return self._select_helix(index)

    def select_indices(self, indices) -> Selection:
        """Select an arbitrary SET of nucleotide indices (box/marquee select).

        The generic set selection used by the view's rubber-band tool: the
        given indices become the highlighted set (kind ``"range"``), replacing
        any prior selection. Out-of-range indices are dropped. No pivot/handle
        is armed (a range is highlight-only, like a motif); the scene tints the
        set teal via the ``selected`` list. This is the clean selected-set API
        a per-region styling pass consumes (see `sel_indices`/`sel_kind`).

        Args:
            indices: Any iterable of nucleotide indices.

        Returns:
            A `Selection` of kind ``"range"`` over the valid indices.
        """
        n = len(self._pair_map)
        valid = {int(i) for i in indices if 0 <= int(i) < n}
        self._clear_selection()
        self._sel_indices = valid
        self._sel_kind = "range"
        return Selection("range", sorted(valid))

    def _select_helix(self, index: int) -> Selection | None:
        """Select the innermost enclosing helix slice and arm its rotation clamp."""
        stack = self._tree.enclosing_pairs.get(index, [])
        if not stack:
            self._clear_selection()
            return None
        start, end = stack[-1]
        parent_pair = stack[-2] if len(stack) >= 2 else None
        parent_loop = self._tree.loop_by_closing_pair.get(parent_pair)
        if parent_loop is not None and parent_loop.members:
            try:
                pivot = loop_center(parent_loop.members, self._x, self._y)
            except ValueError:
                pivot = self._centroid()
        else:
            pivot = self._centroid()
        self._clear_selection()
        self._sel_start, self._sel_end, self._pivot = start, end, pivot
        self._sel_loop = parent_loop
        self._sel_indices = set(range(start, end + 1))
        self._sel_kind = "helix"
        self._arm_rotation_clamp(start, end, pivot, parent_loop)
        return Selection("helix", sorted(self._sel_indices), pivot, start, end)

    def _select_residue(self, index: int) -> Selection:
        """Select a single nucleotide (free-form draggable, no rotation)."""
        self._clear_selection()
        self._sel_indices = {index}
        self._sel_kind = "residue"
        return Selection("residue", [index])

    def _helix_indices(self, index: int) -> set[int]:
        """Every index of the WHOLE physical helix that ``index`` pairs into.

        Walks outward from ``index``'s base pair while consecutive pairs stack
        (``pair_map[a-1] == b+1``) to the helix's outermost rung, then inward
        counting the stacked run -- so the result is BOTH strands of the entire
        contiguous stacked-pair helix (not just the one base pair, and not the
        branch-with-enclosed-loops the rotatable helix slice covers).
        """
        partner = self._pair_map[index]
        a, b = (index, partner) if index < partner else (partner, index)
        # Climb to the outermost stacked rung.
        while a - 1 >= 0 and b + 1 < len(self._pair_map) and self._pair_map[a - 1] == b + 1:
            a, b = a - 1, b + 1
        # Descend the stacked run, collecting both strands.
        indices: set[int] = set()
        while a < b and self._pair_map[a] == b:
            indices.add(a)
            indices.add(b)
            a, b = a + 1, b - 1
        return indices

    def _select_motif(self, index: int) -> Selection:
        """Select the structural element under the click.

        On a stem (paired nucleotide) this is the WHOLE HELIX -- both strands
        of the maximal contiguous stacked-pair run through ``index`` (see
        `_helix_indices`), highlighted but NOT rotatable. On a loop (unpaired
        nucleotide) it is the enclosing loop's members plus the base pairs of
        any child stems, so the whole junction ring highlights.
        """
        if self._pair_map[index] != -1:
            helix = self._helix_indices(index)
            self._clear_selection()
            self._sel_indices = set(helix)
            self._sel_kind = "motif"
            return Selection("motif", sorted(helix))
        stack = self._tree.enclosing_pairs.get(index, [])
        closing = stack[-1] if stack else None
        loop = self._tree.loop_by_closing_pair.get(closing) or self._tree.exterior
        members = set(loop.members)
        for child in loop.children:
            members.add(child.end)
        self._clear_selection()
        self._sel_loop = loop
        self._sel_indices = set(members)
        self._sel_kind = "motif"
        return Selection("motif", sorted(members))

    def _arm_rotation_clamp(
        self, start: int, end: int, pivot: tuple[float, float], loop: Loop | None
    ) -> None:
        """Record the helix's reference direction + allowed deviation band.

        The reference is the direction from ``pivot`` to the helix base
        midpoint at selection. The band is ``[-MAX, +MAX]`` tightened by each
        neighbouring anchor on the parent loop (sibling stems and the loop's
        closing pair), so the helix cannot rotate past a sibling.
        """
        cx, cy = pivot
        mx, my = (self._x[start] + self._x[end]) / 2.0, (self._y[start] + self._y[end]) / 2.0
        ref = math.atan2(my - cy, mx - cx)
        lo, hi = -MAX_HELIX_ROTATION, MAX_HELIX_ROTATION
        if loop is not None:
            anchors: set[int] = set()
            if loop.closing_pair is not None:
                anchors.update(loop.closing_pair)
            for child in loop.children:
                anchors.add(child.start)
                anchors.add(child.end)
            anchors.discard(start)
            anchors.discard(end)
            for a in anchors:
                if not (0 <= a < len(self._x)):
                    continue
                dev = _normalize_pi(math.atan2(self._y[a] - cy, self._x[a] - cx) - ref)
                if dev > 0.0:
                    hi = min(hi, max(0.0, dev - _NEIGHBOR_MARGIN))
                elif dev < 0.0:
                    lo = max(lo, min(0.0, dev + _NEIGHBOR_MARGIN))
        self._sel_ref_angle = ref
        self._sel_dev_lo, self._sel_dev_hi = lo, hi

    def deselect(self) -> None:
        """Clear the current selection and pivot."""
        self._clear_selection()

    def _clear_selection(self) -> None:
        self._sel_start = self._sel_end = None
        self._pivot = None
        self._sel_loop = None
        self._sel_indices = set()
        self._sel_kind = None
        self._sel_ref_angle = None
        self._sel_dev_lo = self._sel_dev_hi = 0.0

    # -- rigid moves (never-silent) -----------------------------------------

    def rotate_selection(self, angle: float, redistribute: bool = True) -> Result:
        """Rotate the selected helix about its pivot, redistribute, re-check, commit.

        The helix's rigid slice is rotated about the parent-loop center
        exactly as before; then, unless ``redistribute`` is False, the
        parent loop's unpaired members are re-spaced evenly on the loop
        circle so the loop stays a clean, compact ring instead of distorting
        (the VARNA hinge; see `rna_draw.gui.hinge.redistribute_loop`). Only
        the dragged helix moves -- sibling helices and the loop's closing
        pair are untouched. When the parent loop has no unpaired members
        (e.g. a stacked-pair continuation), redistribution is a no-op and
        this is a pure rigid rotation.

        Args:
            angle: Rotation angle in radians (counter-clockwise positive).
            redistribute: When True (default), re-space the parent loop's
                unpaired members after rotating. Pass False for a pure
                rigid-slice rotation (the pre-hinge behaviour).

        Returns:
            A `Result` reflecting the native checker's verdict on the
            committed geometry.
        """
        if self._sel_start is None or self._pivot is None:
            return Result(self._flagged, list(self._overlaps))
        angle = self._clamp_rotation(angle)
        cx, cy = self._pivot
        nx, ny = rotate_range(self._x, self._y, self._sel_start, self._sel_end, cx, cy, angle)
        if redistribute and self._sel_loop is not None:
            nx, ny = redistribute_loop(self._sel_loop, nx, ny, (cx, cy))
        return self._commit(nx, ny)

    def translate_selection(self, dx: float, dy: float) -> Result:
        """Rigidly translate the selected helix, re-check, commit.

        Args:
            dx: Translation along engine x.
            dy: Translation along engine y.

        Returns:
            A `Result` reflecting the native checker's verdict.
        """
        if self._sel_start is None:
            return Result(self._flagged, list(self._overlaps))
        nx, ny = translate_range(self._x, self._y, self._sel_start, self._sel_end, dx, dy)
        return self._commit(nx, ny)

    def _clamp_rotation(self, angle: float) -> float:
        """Clamp a requested rotation so the helix stays in its allowed band.

        Reads the helix base's CURRENT deviation from the reference direction
        (armed in `_arm_rotation_clamp`) and returns the largest rotation that
        keeps the cumulative deviation within ``[lo, hi]``. A request beyond
        the bound is clamped to stop AT the limit (the helix never flips or
        jumps); a request within the band is returned unchanged.
        """
        if self._sel_ref_angle is None or self._sel_start is None or self._pivot is None:
            return angle
        cx, cy = self._pivot
        mx = (self._x[self._sel_start] + self._x[self._sel_end]) / 2.0
        my = (self._y[self._sel_start] + self._y[self._sel_end]) / 2.0
        cur_dev = _normalize_pi(math.atan2(my - cy, mx - cx) - self._sel_ref_angle)
        target = min(self._sel_dev_hi, max(self._sel_dev_lo, cur_dev + angle))
        return target - cur_dev

    def move_residue(self, index: int, x: float, y: float) -> Result:
        """Free-form move a SINGLE nucleotide to ``(x, y)``, re-check, commit.

        The residue-granularity edit: only nucleotide ``index`` moves (to the
        given engine-space position); every other coordinate is untouched. The
        never-silent native re-check runs on the committed geometry exactly as
        for the rigid helix moves, so an overlap the drag creates is flagged.

        Args:
            index: Nucleotide index to move.
            x: New engine-space x-coordinate.
            y: New engine-space y-coordinate.

        Returns:
            A `Result` reflecting the native checker's verdict.
        """
        if not (0 <= index < len(self._x)):
            return Result(self._flagged, list(self._overlaps))
        nx, ny = list(self._x), list(self._y)
        nx[index], ny[index] = float(x), float(y)
        return self._commit(nx, ny)

    def _commit(self, nx: list[float], ny: list[float]) -> Result:
        """Commit new coordinates after the never-silent native re-check.

        The authoritative verdict always comes from
        `check_overlaps_native` at the scaled gate radius -- a committed
        layout is never reported clean without that check.
        """
        count = check_overlaps_native(nx, ny, self._pair_map, self._scaled_params())
        self._x, self._y = nx, ny
        self._flagged = count > 0
        self._overlaps = self._overlap_indices() if self._flagged else []
        return Result(self._flagged, list(self._overlaps))

    # -- styling (visual, no geometry change) -------------------------------

    @property
    def preset(self) -> StylePreset:
        """The live `StylePreset` driving this model's rendering."""
        return self._preset

    def set_style(self, preset: StylePreset) -> None:
        """Adopt a new style preset and recompute per-residue fills.

        This is a pure VISUAL update: it never touches coordinates, so it
        cannot introduce an overlap (nothing to re-check). GEOMETRY changes
        (node_r, spacing) must go through `relayout`, never here -- silently
        rescaling stored coordinates is the B3 violation this split avoids.
        Node radius shown on disks stays at the layout's gated `node_r`, so a
        bigger-disk request is a relayout, not a here-and-now resize.
        """
        self._preset = preset
        self._recompute_colors()

    def restyle(self) -> None:
        """Recompute fills from the current preset (after mutating it in place)."""
        self._recompute_colors()

    def _recompute_colors(self) -> None:
        """Recompute per-residue RGB fills via the shared colorer + preset palette.

        Reuses `colorer.Colorer.get_rgb_colors` verbatim -- the same
        `color_str > data > render_type > default` precedence the library
        uses everywhere -- so this never forks the coloring logic. The preset
        palette is applied by scoping the colorer's module `COLORS` table to
        the preset's palette for the duration of the call (restype/paired
        schemes read that table), then restoring it.
        """
        n = len(self._x)
        if n == 0:
            self._colors = None
            return
        from rna_draw import colorer as _colorer

        preset = self._preset
        view = {**DEFAULT_VIEW, **preset.extra.get("view", {})}
        render_type = RENDER_TYPES.get(view.get("render_type", "none"))
        # colorer needs seq/ss of matching length; synthesize placeholders when
        # absent (labels still come from the real seq, so this only feeds color).
        seq = self._seq if self._seq is not None else "N" * n
        ss = self._ss if self._ss is not None else "." * n
        default_rgb = _resolve_rgb(preset.default_color, preset.palette)

        saved = dict(_colorer.COLORS)
        try:
            for key, rgb in preset.palette.items():
                if len(key) == 1:
                    _colorer.COLORS[key] = list(rgb)
            self._colors = _colorer.Colorer().get_rgb_colors(
                seq, ss, render_type=render_type, default_color=default_rgb
            )
        finally:
            _colorer.COLORS.clear()
            _colorer.COLORS.update(saved)

    def set_sequence(self, seq: str) -> str:
        """Update the per-nucleotide sequence labels WITHOUT moving coordinates.

        Pure relabel/recolor: the geometry is untouched (no relayout, no
        re-check needed -- letters and fills can never introduce an overlap).
        The sequence is uppercased, whitespace-stripped, then padded with
        spaces or truncated so its length matches the structure's; a
        human-readable note is returned when that adjustment happened (empty
        string otherwise) so the caller can surface it non-blockingly.

        When a non-empty sequence is applied, letters are turned on so the
        labels are visible; the scene fills are recomputed so an active
        ``res_type``/``paired`` scheme repaints from the new sequence.

        Args:
            seq: The raw sequence text (any letters kept; case-normalized).

        Returns:
            A note string describing any length adjustment, or ``""``.
        """
        n = len(self._x)
        cleaned = "".join(ch for ch in (seq or "").upper() if not ch.isspace())
        note = ""
        if n and len(cleaned) != n:
            note = (
                f"sequence length {len(cleaned)} != structure length {n}; "
                f"{'truncated' if len(cleaned) > n else 'padded'} to fit"
            )
            cleaned = cleaned[:n] if len(cleaned) > n else cleaned + " " * (n - len(cleaned))
        self._seq = cleaned
        if cleaned.strip():
            self._preset.layout_defaults.render_in_letters = True
        self._recompute_colors()
        return note

    def relayout(self, preset: StylePreset | None = None) -> Result:
        """Re-run the guaranteed layout at the preset's geometry, never-silent.

        GEOMETRY knobs (node_r) change the checker's gate radius, so a fresh
        `layout_guaranteed` is run and its verdict adopted -- overlaps are
        never hidden by a visual-only path. Selection is cleared (coordinates
        are wholly replaced). Styling (palette/letters/connectors) carries
        over. Requires an `ss` (present for any `from_ss` model).

        Returns:
            A `Result` mirroring the fresh layout's checker verdict.
        """
        if preset is not None:
            self._preset = preset
        if self._ss is None:
            return Result(self._flagged, list(self._overlaps))
        node_r = float(self._preset.layout_defaults.node_r)
        result = layout_guaranteed(
            self._ss,
            engine=resolve_engine("auto"),
            params=OverlapParams(node_r=node_r),
        )
        pair_map = result.pair_map
        if pair_map is None:
            pair_map = get_pairmap_from_secstruct(self._ss)
        self._x = [float(v) for v in result.x]
        self._y = [float(v) for v in result.y]
        self._node_r = float(result.node_r)
        self._pair_map = [int(v) for v in pair_map]
        self._crossing_pairs = list(result.crossing_pairs or [])
        self._crossing_lines = list(result.crossing_lines or [])
        self._tree = build_structure_tree(self._pair_map)
        self._engine_name = result.engine_name
        self._clear_selection()
        self._recompute_colors()
        self._flagged = bool(result.flagged) or len(self._overlap_indices()) > 0
        self._overlaps = self._overlap_indices() if self._flagged else []
        return Result(self._flagged, list(self._overlaps))

    # -- residue numbering (visual, no geometry change) ---------------------

    def set_numbering(self, show: bool, interval: int = 10) -> None:
        """Toggle VARNA-style residue-position numbers and set their interval.

        Pure visual: the numbers are recomputed from the current geometry in
        `scene()` (never stored), so this can never move a coordinate or
        introduce an overlap.
        """
        self._show_numbers = bool(show)
        self._number_interval = max(1, int(interval))

    def _numbers(self) -> list[dict]:
        """Position labels for residue 1, every Nth residue, and the last.

        Placed OUTSIDE the drawing, in the direction of greatest clearance to
        neighbouring disks so a number never sits on top of a sphere: each label
        scans candidate directions and takes the one whose position is farthest
        from any nearby disk centre, biased outward (away from the structure
        centroid) to break ties. Indices are 0-based; the shown text is 1-based
        (VARNA convention). Coordinates are ENGINE space (y-up); the view flips y.
        """
        n = len(self._x)
        if n == 0:
            return []
        cx, cy = self._centroid()
        step = max(1, self._number_interval)
        indices = sorted({0, n - 1} | {i for i in range(n) if (i + 1) % step == 0})
        r = self._node_r
        offset = r * 2.2  # push the label clear of its own + neighbouring disks
        window = offset + r * 2.0  # only disks this close can collide with a label
        out: list[dict] = []
        for i in indices:
            xi, yi = self._x[i], self._y[i]
            near = [
                j
                for j in range(n)
                if j != i and abs(self._x[j] - xi) < window and abs(self._y[j] - yi) < window
            ]
            ox, oy = xi - cx, yi - cy
            od = math.hypot(ox, oy) or 1.0
            ox, oy = ox / od, oy / od  # preferred outward direction
            best_pos, best_score = (xi + ox * offset, yi + oy * offset), -1.0
            for k in range(16):
                a = 2.0 * math.pi * k / 16.0
                dx, dy = math.cos(a), math.sin(a)
                px, py = xi + dx * offset, yi + dy * offset
                clear = min(
                    (math.hypot(px - self._x[j], py - self._y[j]) for j in near),
                    default=window,
                )
                # clearance dominates; the small outward term only breaks ties.
                score = clear + r * 0.5 * (dx * ox + dy * oy)
                if score > best_score:
                    best_score, best_pos = score, (px, py)
            out.append({"index": i, "x": best_pos[0], "y": best_pos[1], "text": str(i + 1)})
        return out

    # -- per-region styling (visual overrides + highlights) -----------------

    def apply_color_to_selection(self, hex_color: str) -> None:
        """Recolor every currently-selected nucleotide with ``hex_color``.

        Sets a per-nt color OVERRIDE for each index in `sel_indices`. The
        override wins over the scheme/style fill in `scene()` (highest
        precedence); it is layered on top of the colorer output, never forking
        the coloring logic. Purely visual -- no coordinate is touched, so there
        is nothing to re-check (the never-silent contract is unaffected).
        """
        for i in self._sel_indices:
            self._color_overrides[int(i)] = hex_color

    def clear_color_on_selection(self) -> None:
        """Drop any color override on the currently-selected nucleotides.

        Those nts revert to their scheme/style fill (the colorer output).
        """
        for i in self._sel_indices:
            self._color_overrides.pop(int(i), None)

    def clear_all_overrides(self) -> None:
        """Remove every per-nt color override (whole structure)."""
        self._color_overrides = {}

    def highlight_selection(self, hex_color: str) -> None:
        """Add a persistent highlight halo over the currently-selected region.

        Appends a ``{indices, color}`` highlight for `sel_indices`; the scene
        emits it in a ``highlights`` list and the view draws a translucent halo
        behind those disks. Distinct from a color override (emphasis, not a
        repaint) and from the transient teal selection tint. No-op when nothing
        is selected.
        """
        if self._sel_indices:
            self._highlights.append(
                {"indices": sorted(self._sel_indices), "color": hex_color}
            )

    def clear_highlights(self) -> None:
        """Remove every persistent highlight."""
        self._highlights = []

    def remove_highlight(self, index: int) -> None:
        """Remove the highlight at ``index`` in insertion order (no-op if out of range)."""
        if 0 <= index < len(self._highlights):
            del self._highlights[index]

    @property
    def color_overrides(self) -> dict[int, str]:
        """A copy of the per-nt color overrides (index -> ``#rrggbb``)."""
        return dict(self._color_overrides)

    @property
    def highlights(self) -> list[dict]:
        """A copy of the persistent highlights (``{indices, color}``)."""
        return [{"indices": list(h["indices"]), "color": h["color"]} for h in self._highlights]

    def _restore_region_styling(self, doc: Document) -> None:
        """Restore color overrides + highlights from a Document's `extra` band."""
        self._color_overrides = {}
        for key, value in (doc.extra.get("color_overrides") or {}).items():
            try:
                self._color_overrides[int(key)] = str(value)
            except (TypeError, ValueError):
                continue
        self._highlights = []
        for h in doc.extra.get("highlights") or []:
            try:
                indices = [int(i) for i in h.get("indices", [])]
            except (TypeError, ValueError, AttributeError):
                continue
            self._highlights.append(
                {"indices": indices, "color": str(h.get("color", "#ffff00"))}
            )

    # -- scene --------------------------------------------------------------

    def scene(self) -> dict:
        """Serialize the current geometry + verdict to a portable scene dict."""
        scene = build_scene(
            self._x,
            self._y,
            self._pair_map,
            self._node_r,
            seq=self._seq,
            colors=self._colors,
            crossing_pairs=self._crossing_pairs,
            routed_lines=self._crossing_lines,
            flagged=self._flagged,
            overlaps=self._overlaps,
            pivot=self._pivot,
        )
        data = scene.to_dict()
        # Tag each drawn base pair with a `bond` type from the SEQUENCE bases at
        # (i, j): gc / au / gu / other. Crossing (PK) pairs keep kind="crossing"
        # but also carry a bond. Display-only -- consumed by the letters render
        # style to pick a pair symbol (=, -, wobble, LW placeholder).
        seq = self._seq
        for pair in data["pairs"]:
            i, j = pair["i"], pair["j"]
            if seq is not None:
                a = seq[i] if i < len(seq) else ""
                b = seq[j] if j < len(seq) else ""
                pair["bond"] = _bond_type(a, b)
            else:
                pair["bond"] = "other"
        # Per-nt color overrides win over the scheme/style fill (HIGHEST
        # precedence), layered on top of the colorer output in `build_scene`.
        if self._color_overrides:
            for nt in data["nucleotides"]:
                override = self._color_overrides.get(nt["id"])
                if override is not None:
                    nt["fill"] = override
        # Carry resolved visual styling so the view restyles without re-layout.
        data["style"] = view_style(self._preset)
        # The highlighted nucleotide set (any granularity) + its kind, so the
        # view can teal-tint non-contiguous selections (a motif ring) too.
        data["selected"] = sorted(self._sel_indices)
        data["selection_kind"] = self._sel_kind
        # Residue-position numbers (optional; absent for the Jupyter payload).
        if self._show_numbers:
            data["numbers"] = self._numbers()
        # Persistent highlight halos (optional; absent when none defined so the
        # Jupyter payload and unstyled scenes are unchanged).
        if self._highlights:
            data["highlights"] = [
                {"indices": list(h["indices"]), "color": h["color"]} for h in self._highlights
            ]
        return data


__all__ = [
    "EditorModel",
    "Result",
    "Selection",
    "MAX_HELIX_ROTATION",
    "DEMO_SS",
    "DEMO_SEQ",
    "RENDER_TYPES",
    "DEFAULT_VIEW",
    "view_style",
]
