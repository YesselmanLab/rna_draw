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

import math
from dataclasses import dataclass

from rna_draw.gui.hinge import redistribute_loop
from rna_draw.layout.base import params_at_node_r
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

# Non-geometry visual knobs that have no home in the canonical `StylePreset`
# fields. They are carried in `preset.extra["view"]` so they round-trip through
# `.rnastyle.json` save/load (via the preset's `extra` merge) without inventing
# a new file format. Defaults mirror the desktop `scene_view` module constants
# so an unstyled scene renders byte-identically to before.
DEFAULT_VIEW = {
    "background": "#ffffff",
    "backbone_color": "#8a8f98",
    "backbone_width": 2.0,
    "nt_edge_color": "#2b2f36",
    "nt_edge_width": 1.0,
    "letter_color": "#101418",
    "pair_width": 1.6,
    "crossing_width": 1.6,
    "routed_width": 1.6,
    "routed_style": "dash",  # solid | dash | dot
    "render_type": "none",  # none | res_type | paired
    "data_palette": "viridis",  # SHAPE/data colormap name (persisted; see report)
}


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
    return {
        "background": view["background"],
        "backbone_color": view["backbone_color"],
        "backbone_width": float(view["backbone_width"]),
        "default_fill": _resolve_color(preset.default_color, palette),
        "nt_edge_color": view["nt_edge_color"],
        "nt_edge_width": float(view["nt_edge_width"]),
        "show_letters": bool(ld.render_in_letters),
        "letter_scale": float(ld.text_size) / 50.0,
        "letter_color": view["letter_color"],
        "pair_color": _resolve_color(cc.get("nested_pair", "e"), palette),
        "pair_width": float(view["pair_width"]),
        "crossing_color": _resolve_color(cc.get("pk_connector", "o"), palette),
        "crossing_width": float(view["crossing_width"]),
        "routed_color": _resolve_color(cc.get("pk_line", "r"), palette),
        "routed_width": float(view["routed_width"]),
        "routed_style": view["routed_style"],
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
        # Residue-numbering (VARNA-style "10, 20, 30 ..." labels). Off by
        # default; the desktop toggle turns it on. Purely visual -- the
        # numbers are derived from geometry in `scene()`, never stored.
        self._show_numbers: bool = False
        self._number_interval: int = 10

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

    def _select_motif(self, index: int) -> Selection:
        """Select the structural element under the click.

        On a stem (paired nucleotide) this is the enclosing helix; on a loop
        (unpaired nucleotide) it is the enclosing loop's members plus the base
        pairs of any child stems, so the whole junction ring highlights.
        """
        if self._pair_map[index] != -1:
            helix = self._select_helix(index)
            if helix is not None:
                return helix
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

        Placed just OUTSIDE each disk, offset radially away from the whole-
        structure centroid so the number never sits on top of the drawing.
        Indices are 0-based; the shown text is 1-based (VARNA convention).
        Coordinates are ENGINE space (y-up); the view flips y for screen.
        """
        n = len(self._x)
        if n == 0:
            return []
        cx, cy = self._centroid()
        step = max(1, self._number_interval)
        indices = sorted({0, n - 1} | {i for i in range(n) if (i + 1) % step == 0})
        offset = self._node_r * 1.9
        out: list[dict] = []
        for i in indices:
            dx, dy = self._x[i] - cx, self._y[i] - cy
            dist = math.hypot(dx, dy) or 1.0
            out.append(
                {
                    "index": i,
                    "x": self._x[i] + (dx / dist) * offset,
                    "y": self._y[i] + (dy / dist) * offset,
                    "text": str(i + 1),
                }
            )
        return out

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
        # Carry resolved visual styling so the view restyles without re-layout.
        data["style"] = view_style(self._preset)
        # The highlighted nucleotide set (any granularity) + its kind, so the
        # view can teal-tint non-contiguous selections (a motif ring) too.
        data["selected"] = sorted(self._sel_indices)
        data["selection_kind"] = self._sel_kind
        # Residue-position numbers (optional; absent for the Jupyter payload).
        if self._show_numbers:
            data["numbers"] = self._numbers()
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
