"""`Document`: a saved/loaded `.rnadoc.json` drawing (structure + coloring
intent + resolved layout coords).

```json
{
  "schema": "rna_draw/document",
  "version": 1,
  "source": {
    "structure": {"ss": "((((....))))", "seq": "GGGGAAAACCCC"},
    "coloring": {"render_type": "res_type", "color_str": null,
                 "default_color": "e", "data": null},
    "style_preset": "default",
    "style_overrides": {},
    "engine": "auto"
  },
  "derived": {
    "layout": {
      "node_r": 10.0,
      "coords": [[x0, y0], [x1, y1], ...],
      "pair_map": [-1, 5, ...],
      "crossing_pairs": [[i, j], ...],
      "crossing_lines": [{"i": 3, "j": 20, "points": [[x, y], ...]}],
      "engine_name": "pseudoknot",
      "flagged": false,
      "checker": {"verdict": "passed", "node_r": 10.0}
    }
  }
}
```

**Coords are stored UNSHIFTED.** `RNADrawer.__render` (`draw.py`) computes a
shift and calls `RNARenderer.set_coords`, which mutates `result.x`/`result.y`
IN PLACE, subtracting that shift; `crossing_lines` are then separately
re-shifted for rendering. So a `Document`'s layout band must be built from a
snapshot taken BEFORE `set_coords` runs (see
`document_render.document_from_layout`). This is safe to re-validate at load
time because `overlap.check_overlaps` is translation-invariant -- every
predicate in `geometry.py` (`disks_overlap`, `point_segment_distance`,
`segment_segment_distance`) uses only coordinate *differences* -- so
validating the stored UNSHIFTED coords yields the same verdict the batch
pipeline got at shifted coords.

**The `checker` field is advisory only.** `derived.layout.checker` is a
cache written at save time; `validate.never_silent_gate` (called by
`document_render.draw_document`) ALWAYS re-validates live against the frozen
checker plus the routed-line test, overriding whatever this file claims (a
tampered `"passed"` is caught -- see the plan's "Advisory-cache override"
test).

**Unknown-field passthrough.** Every dataclass here (and in
`document_schema`) carries an `extra: dict` that `from_dict` populates with
any key it doesn't recognize and `to_dict` merges back in last, so
`save -> load -> save` is byte-idempotent even across schema versions.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rna_draw.document_schema import (
    ColoringIntent,
    DataIntent,
    DerivedBand,
    LayoutBand,
    SourceBand,
    StructureIntent,
)
from rna_draw.schema_utils import merge_extra, split_known, validate_schema

SCHEMA = "rna_draw/document"
SCHEMA_VERSION = 1

_DOCUMENT_FIELDS = {"source", "derived"}


@dataclass
class Document:
    """A complete, self-contained saveable/loadable drawing.

    Args:
        source: What this drawing was (or will be) produced FROM.
        derived: What running layout actually produced, if it has been run.
        version: Schema version; an unknown future version keeps every
            field this code recognizes and passes the rest through rather
            than hard-failing.
        extra: Unrecognized top-level fields carried through.
    """

    source: SourceBand
    derived: DerivedBand = field(default_factory=DerivedBand)
    version: int = SCHEMA_VERSION
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a full `.rnadoc.json` mapping."""
        base = {
            "schema": SCHEMA,
            "version": self.version,
            "source": self.source.to_dict(),
            "derived": self.derived.to_dict(),
        }
        return merge_extra(base, self.extra)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Document:
        """Parse a `.rnadoc.json`-shaped mapping (schema not re-checked here;
        `load` validates it before calling this).
        """
        return cls(
            source=SourceBand.from_dict(data["source"]),
            derived=DerivedBand.from_dict(data.get("derived", {})),
            version=data.get("version", SCHEMA_VERSION),
            extra=split_known(data, _DOCUMENT_FIELDS | {"schema", "version"}),
        )

    def save(self, path: str | Path) -> None:
        """Write this document to `path` as pretty-printed `.rnadoc.json`."""
        Path(path).write_text(json.dumps(self.to_dict(), indent=2))

    @classmethod
    def load(cls, path: str | Path) -> Document:
        """Load a `.rnadoc.json` file written by `save`.

        Raises:
            ValueError: If the file's `schema` isn't `"rna_draw/document"`.
        """
        data = json.loads(Path(path).read_text())
        validate_schema(data, SCHEMA)
        return cls.from_dict(data)


__all__ = [
    "SCHEMA",
    "SCHEMA_VERSION",
    "Document",
    "SourceBand",
    "ColoringIntent",
    "DataIntent",
    "DerivedBand",
    "LayoutBand",
    "StructureIntent",
]
