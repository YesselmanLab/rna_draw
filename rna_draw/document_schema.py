"""Pure (de)serialization records for `.rnadoc.json` (split out of
`document.py` per the plan's 300-line guidance).

Every dataclass here is a plain schema record: `to_dict`/`from_dict` only,
no engine calls, no rendering. `document.py` composes these into `Document`
and adds the layout-provenance bridge (`document_from_layout`) plus
save/load. See `rna_draw/document.py`'s module docstring for the full
`.rnadoc.json` shape.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from rna_draw.layout.base import RoutedLine
from rna_draw.schema_utils import merge_extra, split_known

_STRUCTURE_FIELDS = {"ss", "seq"}
_DATA_INTENT_FIELDS = {"values", "palette", "vmin", "vmax", "ignore_restype"}
_COLORING_FIELDS = {"render_type", "color_str", "default_color", "data"}
_SOURCE_FIELDS = {"structure", "coloring", "style_preset", "style_overrides", "engine"}
_LAYOUT_FIELDS = {
    "node_r",
    "coords",
    "pair_map",
    "crossing_pairs",
    "crossing_lines",
    "engine_name",
    "flagged",
    "checker",
}


@dataclass
class StructureIntent:
    """The dot-bracket structure + sequence a `Document` was built from.

    Args:
        ss: Dot-bracket secondary structure.
        seq: Sequence, same length as `ss`.
        extra: Unrecognized fields carried through.
    """

    ss: str
    seq: str
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to the `source.structure` JSON object."""
        return merge_extra({"ss": self.ss, "seq": self.seq}, self.extra)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> StructureIntent:
        """Parse a `source.structure` JSON object."""
        return cls(ss=data["ss"], seq=data["seq"], extra=split_known(data, _STRUCTURE_FIELDS))


@dataclass
class DataIntent:
    """Experimental-data coloring intent (mirrors `data.Data`'s ctor args).

    Stored as INTENT, resolved lazily by `document_render.resolve_colors`
    (R1: a `palette` name is resolved to a matplotlib colormap there, since
    `data.Data`/`colorer.color_by_data` need a callable).

    Args:
        values: One data value per nucleotide.
        palette: Matplotlib colormap name (e.g. `"Reds"`).
        vmin: Values below this are clamped; `None` uses `min(values)`.
        vmax: Values above this are clamped; `None` uses `max(values)`.
        ignore_restype: Residue letters (e.g. `"GU"`) drawn with
            `default_color` instead of the data palette.
        extra: Unrecognized fields carried through.
    """

    values: list[float]
    palette: str = "Reds"
    vmin: float | None = None
    vmax: float | None = None
    ignore_restype: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to the `coloring.data` JSON object."""
        base = {
            "values": self.values,
            "palette": self.palette,
            "vmin": self.vmin,
            "vmax": self.vmax,
            "ignore_restype": self.ignore_restype,
        }
        return merge_extra(base, self.extra)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DataIntent:
        """Parse a `coloring.data` JSON object."""
        return cls(
            values=list(data["values"]),
            palette=data.get("palette", "Reds"),
            vmin=data.get("vmin"),
            vmax=data.get("vmax"),
            ignore_restype=data.get("ignore_restype"),
            extra=split_known(data, _DATA_INTENT_FIELDS),
        )


@dataclass
class ColoringIntent:
    """Coloring intent: mutually exclusive `render_type`/`data`, plus
    `color_str` (highest precedence) and `default_color` (base).

    Precedence mirrors `colorer.Colorer.get_rgb_colors` exactly (reused,
    never forked): `color_str > data > render_type > default_color`.

    Args:
        render_type: `"res_type"`, `"paired"`, `"none"`, or `None`.
        color_str: Raw color-string spec, or `None`.
        default_color: Palette key used as the base color, or `None`.
        data: `DataIntent`, or `None`.
        extra: Unrecognized fields carried through.
    """

    render_type: str | None = None
    color_str: str | None = None
    default_color: str | None = None
    data: DataIntent | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to the `source.coloring` JSON object."""
        base: dict[str, Any] = {
            "render_type": self.render_type,
            "color_str": self.color_str,
            "default_color": self.default_color,
            "data": self.data.to_dict() if self.data is not None else None,
        }
        return merge_extra(base, self.extra)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ColoringIntent:
        """Parse a `source.coloring` JSON object."""
        raw_data = data.get("data")
        return cls(
            render_type=data.get("render_type"),
            color_str=data.get("color_str"),
            default_color=data.get("default_color"),
            data=DataIntent.from_dict(raw_data) if raw_data is not None else None,
            extra=split_known(data, _COLORING_FIELDS),
        )


@dataclass
class SourceBand:
    """Everything a layout was (or will be) produced FROM.

    Args:
        structure: `StructureIntent`.
        coloring: `ColoringIntent`.
        style_preset: A registry name (e.g. `"default"`), a path, or an
            inline preset mapping.
        style_overrides: Doc-level overrides for the non-geometry knobs P1
            understands (`default_color`, `render_in_letters`); see
            `document_render.resolve_style`.
        engine: Layout-engine selection string (`"auto"`, `"legacy"`, ...).
        extra: Unrecognized fields carried through.
    """

    structure: StructureIntent
    coloring: ColoringIntent = field(default_factory=ColoringIntent)
    style_preset: str | dict[str, Any] = "default"
    style_overrides: dict[str, Any] = field(default_factory=dict)
    engine: str = "auto"
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to the `source` JSON object."""
        base = {
            "structure": self.structure.to_dict(),
            "coloring": self.coloring.to_dict(),
            "style_preset": self.style_preset,
            "style_overrides": self.style_overrides,
            "engine": self.engine,
        }
        return merge_extra(base, self.extra)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SourceBand:
        """Parse a `source` JSON object."""
        return cls(
            structure=StructureIntent.from_dict(data["structure"]),
            coloring=ColoringIntent.from_dict(data.get("coloring", {})),
            style_preset=data.get("style_preset", "default"),
            style_overrides=data.get("style_overrides", {}),
            engine=data.get("engine", "auto"),
            extra=split_known(data, _SOURCE_FIELDS),
        )


def _routed_line_to_dict(line: RoutedLine) -> dict[str, Any]:
    """Serialize one `RoutedLine` to its `crossing_lines` JSON object."""
    return {"i": line.i, "j": line.j, "points": [[px, py] for px, py in line.points]}


def _routed_line_from_dict(data: dict[str, Any]) -> RoutedLine:
    """Parse one `crossing_lines` JSON object into a `RoutedLine`."""
    points = [(float(px), float(py)) for px, py in data["points"]]
    return RoutedLine(i=data["i"], j=data["j"], points=points)


@dataclass
class LayoutBand:
    """The resolved, checker-gated layout provenance to persist.

    Coordinates are stored UNSHIFTED (see `document.py`'s module docstring
    for why); `node_r` is AUTHORITATIVE for these coords and is never
    rescaled by a preset swap (B3) -- only an explicit
    `document_render.relayout` changes geometry.

    Args:
        node_r: Disk radius these coords were checker-gated/drawn at
            (`LayoutResult.node_r`; may be `< target` if the pipeline's
            adaptive search shrank it -- see RC1).
        coords: One `(x, y)` per nucleotide, sequence order, UNSHIFTED.
        pair_map: The FULL drawn base-pair map (`LayoutResult.pair_map`),
            or an all-`-1` map if the source had none.
        crossing_pairs: PK-B in-plane crossing connectors.
        crossing_lines: PK-A routed polylines.
        engine_name: Which engine produced this layout.
        flagged: `LayoutResult.flagged`.
        checker: Advisory cache `{"verdict": ..., "node_r": ...}` -- NEVER
            trusted on load; `validate.never_silent_gate` always re-checks
            live (see RC1/RC2).
        extra: Unrecognized fields carried through (e.g. a future GUI's
            `edited_indices`; NOT populated by P1, see plan Risk R2).
    """

    node_r: float
    coords: list[tuple[float, float]]
    pair_map: list[int]
    crossing_pairs: list[tuple[int, int]] = field(default_factory=list)
    crossing_lines: list[RoutedLine] = field(default_factory=list)
    engine_name: str = ""
    flagged: bool = False
    checker: dict[str, Any] = field(default_factory=dict)
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to the `derived.layout` JSON object."""
        base = {
            "node_r": self.node_r,
            "coords": [[x, y] for x, y in self.coords],
            "pair_map": self.pair_map,
            "crossing_pairs": [[i, j] for i, j in self.crossing_pairs],
            "crossing_lines": [_routed_line_to_dict(line) for line in self.crossing_lines],
            "engine_name": self.engine_name,
            "flagged": self.flagged,
            "checker": self.checker,
        }
        return merge_extra(base, self.extra)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> LayoutBand:
        """Parse a `derived.layout` JSON object."""
        return cls(
            node_r=data["node_r"],
            coords=[(float(x), float(y)) for x, y in data["coords"]],
            pair_map=list(data["pair_map"]),
            crossing_pairs=[(i, j) for i, j in data.get("crossing_pairs", [])],
            crossing_lines=[
                _routed_line_from_dict(line) for line in data.get("crossing_lines", [])
            ],
            engine_name=data.get("engine_name", ""),
            flagged=data.get("flagged", False),
            checker=data.get("checker", {}),
            extra=split_known(data, _LAYOUT_FIELDS),
        )


@dataclass
class DerivedBand:
    """Everything derived FROM `SourceBand` by actually running layout.

    Args:
        layout: `LayoutBand`, or `None` for a source-only document (no
            layout has been run yet; `document_render.draw_document` then
            falls through to `relayout`).
    """

    layout: LayoutBand | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize to the `derived` JSON object."""
        return {"layout": self.layout.to_dict() if self.layout is not None else None}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DerivedBand:
        """Parse a `derived` JSON object."""
        raw_layout = data.get("layout")
        return cls(layout=LayoutBand.from_dict(raw_layout) if raw_layout is not None else None)


__all__ = [
    "StructureIntent",
    "DataIntent",
    "ColoringIntent",
    "SourceBand",
    "LayoutBand",
    "DerivedBand",
]
