"""Shared JSON-schema (de)serialization helpers for `style.py`/`document.py`.

Every persisted dataclass in the P1 persistence layer (`style.LayoutDefaults`,
`style.StylePreset`, `document_schema.StructureIntent`/`ColoringIntent`/
`DataIntent`/`SourceBand`/`LayoutBand`/`DerivedBand`, `document.Document`)
carries an `extra: dict` passthrough field so `save -> load -> save` stays
byte-idempotent even across schema versions, and so a future field (e.g. a
GUI's `edited_indices`) survives a round trip through a version of this code
that doesn't understand it yet. Factored out once the pattern crossed three
dataclasses -- reuse over repetition.
"""

from __future__ import annotations

from typing import Any


def split_known(data: dict[str, Any], known: set[str]) -> dict[str, Any]:
    """Return the subset of `data` whose keys are NOT in `known`.

    Args:
        data: A raw JSON-decoded mapping.
        known: The field names a dataclass's `from_dict` already consumes.

    Returns:
        Every `data` entry `from_dict` doesn't recognize, to stash in the
        dataclass's own `extra` field.
    """
    return {key: value for key, value in data.items() if key not in known}


def merge_extra(base: dict[str, Any], extra: dict[str, Any]) -> dict[str, Any]:
    """Merge a dataclass's `extra` passthrough back into its `to_dict` output.

    Args:
        base: The fields `to_dict` builds from its own known attributes.
        extra: The dataclass's `extra` dict (unknown fields carried through
            from a prior `from_dict`).

    Returns:
        `base` with `extra` merged in last, so a round trip never loses an
        unrecognized field.
    """
    merged = dict(base)
    merged.update(extra)
    return merged


def validate_schema(data: dict[str, Any], expected_schema: str) -> None:
    """Guard a document/preset's `schema` tag before parsing its fields.

    Args:
        data: The raw JSON-decoded top-level mapping.
        expected_schema: The exact `schema` string this loader requires
            (e.g. `"rna_draw/document"`).

    Raises:
        ValueError: If `data["schema"]` is missing or does not match
            `expected_schema`. An unknown `version` is NOT an error (see
            each `from_dict`'s own version-tolerant field handling) -- only
            the wrong schema family is a hard failure.
    """
    schema = data.get("schema")
    if schema != expected_schema:
        raise ValueError(f"expected schema {expected_schema!r}, got {schema!r}")
