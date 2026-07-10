"""Regression baseline test locking the legacy layout engine's coordinates.

These tests protect ``rna_draw.render_rna.RNARenderer.setup_tree`` against
silent coordinate changes while the surrounding modernization work (M1) is
applied. Only per-nucleotide (x, y) coordinate arrays are asserted -- never
rendered pixels/PNG bytes.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import rna_draw as rd
from rna_draw.render_rna import RNARenderer

BASELINE_PATH = Path(__file__).parent / "resources" / "layout_baseline.json"
NODE_R = 10
PRIMARY_SPACE = 20
PAIR_SPACE = 23
ATOL = 1e-9


def _load_baseline() -> dict[str, Any]:
    """Load the committed coordinate baseline JSON."""
    with BASELINE_PATH.open() as handle:
        return json.load(handle)


def _compute_coords(secstruct: str) -> tuple[list[float], list[float]]:
    """Run ``setup_tree`` and return its resulting coordinate arrays."""
    renderer = RNARenderer()
    renderer.setup_tree(
        secstruct, NODE_R=NODE_R, PRIMARY_SPACE=PRIMARY_SPACE, PAIR_SPACE=PAIR_SPACE
    )
    return list(renderer.xarray_), list(renderer.yarray_)


def _assert_all_close(actual: list[float], expected: list[float]) -> None:
    """Assert two equal-length float lists match within the fixed tolerance."""
    assert len(actual) == len(expected)
    for got, want in zip(actual, expected):
        assert math.isclose(got, want, rel_tol=0.0, abs_tol=ATOL)


def _assert_matches_baseline(name: str) -> None:
    """Assert recomputed coordinates match the committed baseline for `name`."""
    fixture = _load_baseline()[name]
    xarray, yarray = _compute_coords(fixture["ss"])
    _assert_all_close(xarray, fixture["xarray"])
    _assert_all_close(yarray, fixture["yarray"])


def test_baseline_hairpin() -> None:
    """Coords for a simple hairpin match the committed baseline."""
    _assert_matches_baseline("hairpin")


def test_baseline_dangling() -> None:
    """Coords for a structure with dangling ends match (guards the dead-code fix)."""
    _assert_matches_baseline("dangling")


def test_baseline_multiloop() -> None:
    """Coords for a two-branch multiloop junction match the committed baseline."""
    _assert_matches_baseline("multiloop")


def test_baseline_trna_midsize() -> None:
    """Coords for a mid-size tRNA-like structure match the committed baseline."""
    _assert_matches_baseline("trna_midsize")


def test_coord_array_lengths() -> None:
    """Every fixture's coordinate arrays are the same length as its secstruct."""
    baseline = _load_baseline()
    for fixture in baseline.values():
        xarray, yarray = _compute_coords(fixture["ss"])
        assert len(xarray) == len(fixture["ss"])
        assert len(yarray) == len(fixture["ss"])


def test_smoke_render_writes_requested_path(tmp_path: Path) -> None:
    """The public ``rna_draw`` API writes a PNG at the caller-requested path."""
    out = tmp_path / "smoke"
    rd.rna_draw(ss="((((....))))", seq="GGGGAAAACCCC", out=str(out))
    assert (tmp_path / "smoke.png").exists()
