"""P1 style preset: layout geometry + coloring, saved/loaded as `.rnastyle.json`.

`StylePreset` is the persisted counterpart of `parameters.DrawParameters` +
`colorer.COLORS` + the PK connector-color constants (`draw.py`'s
`PK_CONNECTOR_COLOR`/`PK_LINE_COLOR`). `default_preset()` reproduces the
library's built-in defaults exactly (`to_draw_parameters()` equals
`parameters.DrawParameters()` field-for-field).

Scope note (R4, mirrored from the plan): `rules` is a forward-compat
placeholder for the P2b CSS-like cascade. P1's matplotlib renderer can't
consume arbitrary rules, so they are parsed and round-tripped (`extra`-style
fidelity) but never applied; only `layout_defaults` and `default_color` (via
`document_render.resolve_style`) drive P1 rendering. This is a documented
limitation, not a silent gap -- see `document_render.resolve_style`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rna_draw import parameters
from rna_draw.colorer import COLORS
from rna_draw.schema_utils import merge_extra, split_known, validate_schema

SCHEMA = "rna_draw/style"
SCHEMA_VERSION = 1

# Mirrors `draw.py`'s `PK_CONNECTOR_COLOR = COLORS["o"]`, `PK_LINE_COLOR =
# COLORS["r"]`, and the plain nested-pair default `COLORS["e"]` -- persisted
# here for schema completeness/forward-compat; P1's actual render bridge
# still draws with `draw.py`'s own constants (see `document_render`'s
# docstring for why this isn't yet wired into rendering).
DEFAULT_CONNECTOR_COLORS = {"nested_pair": "e", "pk_connector": "o", "pk_line": "r"}

_LAYOUT_DEFAULTS_FIELDS = {
    "node_r",
    "primary_space",
    "pair_space",
    "cell_padding",
    "text_size",
    "render_in_letters",
    "output_width",
    "output_height",
}
_STYLE_PRESET_FIELDS = {
    "name",
    "layout_defaults",
    "palette",
    "connector_colors",
    "default_color",
    "rules",
}


@dataclass
class LayoutDefaults:
    """Layout INPUT geometry: used to lay out a new structure, or on an
    explicit `document_render.relayout` -- it never rescales coordinates
    already stored in a `Document` (B3; see `document_render`'s docstring).

    Args:
        node_r: Nucleotide disk radius, layout units.
        primary_space: Backbone step target, layout units.
        pair_space: Base-pair connector length target, layout units.
        cell_padding: Render offset added to every coordinate.
        text_size: Font size for `RENDER_IN_LETTERS` mode.
        render_in_letters: Draw nucleotide letters instead of disks.
        output_width: Target PNG width in pixels.
        output_height: Target PNG height in pixels.
        extra: Unrecognized fields carried through from a newer schema.
    """

    node_r: float = 10.0
    primary_space: float = 20.0
    pair_space: float = 23.0
    cell_padding: float = 40.0
    text_size: float = 50.0
    render_in_letters: bool = False
    output_width: int = 1280
    output_height: int = 960
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to the `layout_defaults` JSON object."""
        base = {
            "node_r": self.node_r,
            "primary_space": self.primary_space,
            "pair_space": self.pair_space,
            "cell_padding": self.cell_padding,
            "text_size": self.text_size,
            "render_in_letters": self.render_in_letters,
            "output_width": self.output_width,
            "output_height": self.output_height,
        }
        return merge_extra(base, self.extra)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> LayoutDefaults:
        """Parse a `layout_defaults` JSON object, defaulting missing fields."""
        defaults = cls()
        return cls(
            node_r=data.get("node_r", defaults.node_r),
            primary_space=data.get("primary_space", defaults.primary_space),
            pair_space=data.get("pair_space", defaults.pair_space),
            cell_padding=data.get("cell_padding", defaults.cell_padding),
            text_size=data.get("text_size", defaults.text_size),
            render_in_letters=data.get("render_in_letters", defaults.render_in_letters),
            output_width=data.get("output_width", defaults.output_width),
            output_height=data.get("output_height", defaults.output_height),
            extra=split_known(data, _LAYOUT_DEFAULTS_FIELDS),
        )


def _default_palette() -> dict[str, list[float]]:
    """Copy `colorer.COLORS` byte-for-byte into a plain JSON-able dict."""
    return {key: list(value) for key, value in COLORS.items()}


@dataclass
class StylePreset:
    """A named, persisted rendering style: geometry + palette + coloring.

    Args:
        name: Preset name (used as a registry key by `document_render`).
        layout_defaults: See `LayoutDefaults`.
        palette: Color-code -> `[r, g, b]` (0-1 floats), reproducing
            `colorer.COLORS` by default.
        connector_colors: Which `palette` key colors nested pairs, PK-B
            in-plane connectors, and PK-A routed lines (schema-complete;
            not yet wired into P1 rendering -- see the module docstring).
        default_color: `palette` key used when no other color is supplied.
        rules: Forward-compat CSS-like cascade placeholder (P2b); parsed
            and round-tripped, never applied in P1 (R4).
        extra: Unrecognized top-level fields carried through.
    """

    name: str = "default"
    layout_defaults: LayoutDefaults = field(default_factory=LayoutDefaults)
    palette: dict[str, list[float]] = field(default_factory=_default_palette)
    connector_colors: dict[str, str] = field(default_factory=lambda: dict(DEFAULT_CONNECTOR_COLORS))
    default_color: str = "e"
    rules: list[Any] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)

    def to_draw_parameters(self) -> parameters.DrawParameters:
        """Build a `DrawParameters` from `layout_defaults`.

        Returns:
            A `parameters.DrawParameters` whose fields mirror
            `layout_defaults`; `default_preset()`'s output is equal
            field-for-field to `parameters.DrawParameters()`.
        """
        draw_params = parameters.DrawParameters()
        ld = self.layout_defaults
        draw_params.NODE_R = ld.node_r
        draw_params.PRIMARY_SPACE = ld.primary_space
        draw_params.PAIR_SPACE = ld.pair_space
        draw_params.CELL_PADDING = ld.cell_padding
        draw_params.TEXT_SIZE = ld.text_size
        draw_params.RENDER_IN_LETTERS = ld.render_in_letters
        draw_params.output_width = ld.output_width
        draw_params.output_height = ld.output_height
        return draw_params

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a full `.rnastyle.json` document (with schema/version)."""
        base = {
            "schema": SCHEMA,
            "version": SCHEMA_VERSION,
            "name": self.name,
            "layout_defaults": self.layout_defaults.to_dict(),
            "palette": self.palette,
            "connector_colors": self.connector_colors,
            "default_color": self.default_color,
            "rules": self.rules,
        }
        return merge_extra(base, self.extra)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> StylePreset:
        """Parse a `.rnastyle.json`-shaped mapping.

        Args:
            data: A raw JSON-decoded mapping (schema/version already
                validated by the caller if it came straight off disk).

        Returns:
            A `StylePreset`; unrecognized top-level keys land in `extra`.
        """
        defaults = cls()
        return cls(
            name=data.get("name", defaults.name),
            layout_defaults=LayoutDefaults.from_dict(data.get("layout_defaults", {})),
            palette=data.get("palette", defaults.palette),
            connector_colors=data.get("connector_colors", defaults.connector_colors),
            default_color=data.get("default_color", defaults.default_color),
            rules=data.get("rules", []),
            extra=split_known(data, _STYLE_PRESET_FIELDS | {"schema", "version"}),
        )

    def save(self, path: str | Path) -> None:
        """Write this preset to `path` as pretty-printed `.rnastyle.json`."""
        Path(path).write_text(json.dumps(self.to_dict(), indent=2))

    @classmethod
    def load(cls, path: str | Path) -> StylePreset:
        """Load a `.rnastyle.json` file written by `save`.

        Raises:
            ValueError: If the file's `schema` isn't `"rna_draw/style"`.
        """
        data = json.loads(Path(path).read_text())
        validate_schema(data, SCHEMA)
        return cls.from_dict(data)


def default_preset() -> StylePreset:
    """Build the library's built-in style, byte-for-byte.

    Returns:
        A `StylePreset` whose `to_draw_parameters()` equals
        `parameters.DrawParameters()` field-for-field and whose `palette`
        equals `colorer.COLORS` field-for-field.
    """
    return StylePreset()


__all__ = [
    "SCHEMA",
    "SCHEMA_VERSION",
    "DEFAULT_CONNECTOR_COLORS",
    "LayoutDefaults",
    "StylePreset",
    "default_preset",
]
