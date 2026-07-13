"""Bridge between `Document`/`StylePreset` and the existing matplotlib renderer.

`draw_document` is the load/draw flow's never-silent gate: if a document
already carries a resolved layout (`derived.layout`), its stored coordinates
are re-validated live (`validate.never_silent_gate`, RC1/RC2) and drawn as-is
-- NEVER re-laid-out, and NEVER silently drawn if they fail. A source-only
document (no stored layout yet) falls through to `relayout`, which is also
the ONLY path an explicit geometry/preset change takes effect (B3): a preset
swap alone never rescales coordinates already stored in a `Document`.

Reuses, never forks: `colorer.Colorer.get_rgb_colors` (its own
`color_str > data > render_type > default_color` precedence, `colorer.py`),
`draw.py`'s `_drawn_pairs`/`_shift_routed_lines`/`_set_bounds_and_size`
helpers and `PK_LINE_COLOR`/`_get_render_type`, and `render_rna.RNARenderer`.

Scope note (R4): `StylePreset.rules` is a forward-compat placeholder for the
P2b CSS-like cascade; the matplotlib path can't consume arbitrary rules, so
`resolve_style` only applies `layout_defaults` (geometry) and the two
non-geometry knobs P1 understands (`default_color`, `render_in_letters`) --
full rule application is deferred to P2b, not silently dropped.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import matplotlib.cm
from matplotlib.figure import Figure

from rna_draw import colorer, parameters
from rna_draw.data import Data
from rna_draw.document import (
    ColoringIntent,
    DataIntent,
    DerivedBand,
    Document,
    LayoutBand,
    SourceBand,
    StructureIntent,
)
from rna_draw.draw import (
    PK_LINE_COLOR,
    _drawn_pairs,
    _get_render_type,
    _set_bounds_and_size,
    _shift_routed_lines,
)
from rna_draw.layout import LayoutResult, layout_guaranteed, resolve_engine
from rna_draw.layout.base import empty_report
from rna_draw.overlap import OverlapParams, OverlapReport
from rna_draw.render_rna import RNARenderer, get_pairmap_from_secstruct
from rna_draw.style import StylePreset, default_preset
from rna_draw.validate import LayoutOverlapError, never_silent_gate

_PRESET_REGISTRY: dict[str, Callable[[], StylePreset]] = {"default": default_preset}


@dataclass
class SourceInputs:
    """Ergonomic bundle of structure + coloring + style intent.

    Mirrors `SourceBand`'s fields without requiring its nested dataclasses
    to be pre-built -- the argument `document_from_layout`/`relayout` take.

    Args:
        ss: Dot-bracket secondary structure.
        seq: Sequence, same length as `ss`.
        render_type: `"res_type"`, `"paired"`, `"none"`, or `None`.
        color_str: Raw color-string spec, or `None`.
        default_color: Palette key used as the base color, or `None`.
        data: `DataIntent`, or `None`.
        style_preset: A registry name, a path, or an inline preset mapping.
        style_overrides: Doc-level overrides (`default_color`,
            `render_in_letters`).
        engine: Layout-engine selection string.
    """

    ss: str
    seq: str
    render_type: str | None = None
    color_str: str | None = None
    default_color: str | None = None
    data: DataIntent | None = None
    style_preset: str | dict[str, Any] = "default"
    style_overrides: dict[str, Any] = field(default_factory=dict)
    engine: str = "auto"


def resolve_preset(preset: StylePreset | str | dict[str, Any] | None) -> StylePreset:
    """Normalize a `draw_document`/`relayout` `preset` argument.

    Args:
        preset: A `StylePreset`, a registry name (currently just
            `"default"`, R6 -- no plugin system), a path to a
            `.rnastyle.json` file, an inline preset mapping, or `None`.

    Returns:
        A `StylePreset`. `None` and registry names resolve to a FRESH
        instance each call (the registry stores factories, not shared
        mutable objects, so one caller mutating its preset never leaks to
        another).
    """
    if preset is None:
        return default_preset()
    if isinstance(preset, StylePreset):
        return preset
    if isinstance(preset, dict):
        return StylePreset.from_dict(preset)
    if preset in _PRESET_REGISTRY:
        return _PRESET_REGISTRY[preset]()
    return StylePreset.load(preset)


def resolve_style(
    preset: StylePreset, doc: Document
) -> tuple[parameters.DrawParameters, str | None]:
    """Compose preset defaults + doc-level overrides into effective draw settings.

    Full precedence (low -> high, per the plan): preset defaults <
    preset.rules < doc.source.style_overrides < coloring intent <
    per-residue. P1 collapses this to the two non-geometry knobs the
    matplotlib path understands: `render_in_letters` (style_overrides can
    flip it) and `default_color` (style_overrides, then coloring intent,
    wins over the preset). Layout GEOMETRY is never touched here (B3) --
    only `relayout` changes it. `rules` is parsed/round-tripped but not
    applied (R4).

    Args:
        preset: The resolved `StylePreset`.
        doc: The `Document` being rendered.

    Returns:
        `(draw_params, default_color)`.
    """
    draw_params = preset.to_draw_parameters()
    overrides = doc.source.style_overrides
    if "render_in_letters" in overrides:
        draw_params.RENDER_IN_LETTERS = overrides["render_in_letters"]
    default_color = doc.source.coloring.default_color
    if default_color is None:
        default_color = overrides.get("default_color", preset.default_color)
    return draw_params, default_color


def _data_from_intent(intent: DataIntent) -> Data:
    """Build a `Data` object from a `DataIntent`.

    R1: resolves `intent.palette`'s NAME to a matplotlib colormap callable
    BEFORE constructing `Data` (`colorer.color_by_data` calls
    `data.palette(norm(d))`, so it needs a callable; the CLI's own
    `-data_palette` string path is broken upstream -- see `data.py`'s
    `Data.__init__` -- this fixes it for the document route without
    touching `data.py`'s frozen-for-this-plan semantics).

    Args:
        intent: The stored data-coloring intent.

    Returns:
        A `Data` object ready for `colorer.Colorer.get_rgb_colors`.
    """
    palette = matplotlib.cm.get_cmap(intent.palette)
    data_str = ";".join(str(value) for value in intent.values)
    return Data(
        data_str=data_str,
        palette=palette,
        vmin=intent.vmin,
        vmax=intent.vmax,
        ignore_restype=intent.ignore_restype,
    )


def resolve_colors(doc: Document, preset: StylePreset) -> list[list[float]]:
    """Build per-nucleotide RGB colors for a `Document`.

    Reuses `colorer.Colorer.get_rgb_colors` verbatim (never forked), which
    preserves its own `color_str > data > render_type > default_color`
    precedence.

    Args:
        doc: The document being colored.
        preset: The resolved `StylePreset` (only `default_color` feeds in
            here, via `resolve_style`; palette/connector-color fidelity is
            schema-only in P1, see the module docstring).

    Returns:
        One `[r, g, b]` per nucleotide, sequence order.
    """
    coloring = doc.source.coloring
    structure = doc.source.structure
    _, default_color_key = resolve_style(preset, doc)
    default_rgb = colorer.parse_color_code(default_color_key) if default_color_key else None
    render_type = _get_render_type(coloring.render_type)
    data_obj = _data_from_intent(coloring.data) if coloring.data is not None else None
    return colorer.Colorer().get_rgb_colors(
        structure.seq, structure.ss, coloring.color_str, data_obj, render_type, default_rgb
    )


def _source_band(source: SourceInputs) -> SourceBand:
    """Assemble a `SourceBand` from the ergonomic `SourceInputs` bundle."""
    return SourceBand(
        structure=StructureIntent(ss=source.ss, seq=source.seq),
        coloring=ColoringIntent(
            render_type=source.render_type,
            color_str=source.color_str,
            default_color=source.default_color,
            data=source.data,
        ),
        style_preset=source.style_preset,
        style_overrides=source.style_overrides,
        engine=source.engine,
    )


def document_from_layout(source: SourceInputs, result: LayoutResult) -> Document:
    """Build a `Document` from source intent + a PRE-`set_coords` layout snapshot.

    MUST be called with `result` before any `RNARenderer.set_coords` call --
    `set_coords` mutates `result.x`/`result.y` IN PLACE (see `document.py`'s
    "Coords must be stored UNSHIFTED").

    Args:
        source: Structure + coloring + style intent.
        result: `layout_guaranteed`'s result, unmutated.

    Returns:
        A `Document` with a populated `derived.layout` band.
    """
    pair_map = result.pair_map
    if pair_map is None:
        pair_map = get_pairmap_from_secstruct(source.ss)
    coords = [(float(px), float(py)) for px, py in zip(result.x, result.y)]
    checker = {
        "verdict": "passed" if result.report.passed else "overlap",
        "node_r": result.node_r,
    }
    layout = LayoutBand(
        node_r=result.node_r,
        coords=coords,
        pair_map=list(pair_map),
        crossing_pairs=list(result.crossing_pairs),
        crossing_lines=list(result.crossing_lines),
        engine_name=result.engine_name,
        flagged=result.flagged,
        checker=checker,
    )
    return Document(source=_source_band(source), derived=DerivedBand(layout=layout))


def relayout(doc: Document, preset: StylePreset | str | dict[str, Any] | None = None) -> Document:
    """Lay `doc.source` out fresh, dropping any existing `derived.layout` (B3).

    The ONLY path a geometry change takes effect: swapping presets on a
    document with existing stored coords never rescales them by itself.

    Args:
        doc: The document; its `source` band is reused, `derived.layout`
            is discarded.
        preset: Resolved the same way as `draw_document`'s `preset`.

    Returns:
        A NEW `Document` with a freshly computed `derived.layout` band.
    """
    resolved_preset = resolve_preset(preset)
    structure = doc.source.structure
    engine = resolve_engine(doc.source.engine)
    params = OverlapParams(node_r=resolved_preset.layout_defaults.node_r)
    result = layout_guaranteed(structure.ss, engine=engine, params=params)
    coloring = doc.source.coloring
    source = SourceInputs(
        ss=structure.ss,
        seq=structure.seq,
        render_type=coloring.render_type,
        color_str=coloring.color_str,
        default_color=coloring.default_color,
        data=coloring.data,
        style_preset=doc.source.style_preset,
        style_overrides=doc.source.style_overrides,
        engine=doc.source.engine,
    )
    return document_from_layout(source, result)


def _layout_result_from_band(layout: LayoutBand) -> LayoutResult:
    """Reconstruct a `LayoutResult` from a `LayoutBand` for reuse with
    `draw.py`'s existing pair/bounds helpers.

    `report` is a placeholder (`empty_report()`) -- the never-silent gate's
    own LIVE `OverlapReport` (from `validate.never_silent_gate`), not this
    one, is what decides whether the document may be drawn.

    Args:
        layout: The stored layout band.

    Returns:
        A `LayoutResult` equivalent to the one the pipeline originally
        produced (minus its transient `report`).
    """
    xs = [point[0] for point in layout.coords]
    ys = [point[1] for point in layout.coords]
    return LayoutResult(
        x=xs,
        y=ys,
        engine_name=layout.engine_name,
        report=empty_report(),
        flagged=layout.flagged,
        node_r=layout.node_r,
        pair_map=list(layout.pair_map),
        crossing_pairs=list(layout.crossing_pairs),
        crossing_lines=list(layout.crossing_lines),
    )


def _overlap_message(report: OverlapReport, dirty_lines: list) -> str:
    """Human-readable summary for a `LayoutOverlapError`."""
    parts = []
    if not report.passed:
        parts.append(f"{report.num_overlaps} disk/backbone/pair overlap(s)")
    if dirty_lines:
        parts.append(f"{len(dirty_lines)} routed line(s) crossing the layout")
    return "document's stored coordinates are not drawable: " + "; ".join(parts)


def _render_stored_layout(doc: Document, preset: StylePreset, filename: str) -> Figure:
    """Gate, then draw, a document's already-resolved `derived.layout`.

    Args:
        doc: A document whose `derived.layout` is not `None`.
        preset: The resolved `StylePreset`.
        filename: Output PNG path stem.

    Returns:
        The matplotlib `Figure`.

    Raises:
        LayoutOverlapError: If the live never-silent gate (RC1/RC2) fails.
    """
    layout = doc.derived.layout
    assert layout is not None
    xs = [point[0] for point in layout.coords]
    ys = [point[1] for point in layout.coords]
    clean, report, dirty = never_silent_gate(
        xs,
        ys,
        layout.pair_map,
        layout.crossing_lines,
        layout.node_r,
        preset.layout_defaults.node_r,
    )
    if not clean:
        raise LayoutOverlapError(_overlap_message(report, dirty))

    colors = resolve_colors(doc, preset)
    draw_params, _ = resolve_style(preset, doc)
    result = _layout_result_from_band(layout)
    pairs = _drawn_pairs(result, doc.source.structure.ss)

    renderer = RNARenderer()
    min_x = min(result.x) - result.node_r
    min_y = min(result.y) - result.node_r
    renderer.set_coords(result.x, result.y, result.node_r)
    lines = _shift_routed_lines(result.crossing_lines, min_x, min_y)

    renderer.ax.axis("off")
    _set_bounds_and_size(renderer, draw_params, lines)
    renderer.draw(
        draw_params.CELL_PADDING,
        draw_params.CELL_PADDING,
        colors,
        pairs,
        doc.source.structure.seq,
        draw_params.RENDER_IN_LETTERS,
    )
    renderer.draw_routed_lines(
        lines, draw_params.CELL_PADDING, draw_params.CELL_PADDING, PK_LINE_COLOR
    )
    renderer.fig.savefig(filename + ".png")
    return renderer.fig


def draw_document(
    doc: Document,
    preset: StylePreset | str | dict[str, Any] | None = None,
    filename: str = "secstruct",
) -> Figure:
    """Draw a `Document`, never silently drawing an overlap.

    If `doc.derived.layout` is present, its stored coordinates are
    re-validated live and rendered AS-IS (never re-laid-out); if it is
    absent (a source-only document), `relayout` runs first.

    Args:
        doc: The document to draw.
        preset: A `StylePreset`, registry name, path, inline mapping, or
            `None` (uses `default_preset()`).
        filename: Output PNG path stem (mirrors `RNADrawer.draw`).

    Returns:
        The matplotlib `Figure` -- P1 does not change the return type.

    Raises:
        LayoutOverlapError: If a stored layout fails the live never-silent
            gate (RC1/RC2).
    """
    resolved_preset = resolve_preset(preset)
    if doc.derived.layout is not None:
        return _render_stored_layout(doc, resolved_preset, filename)
    relaid = relayout(doc, resolved_preset)
    return _render_stored_layout(relaid, resolved_preset, filename)


__all__ = [
    "SourceInputs",
    "resolve_preset",
    "resolve_style",
    "resolve_colors",
    "document_from_layout",
    "relayout",
    "draw_document",
]
