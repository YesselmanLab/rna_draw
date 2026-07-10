"""Tests for `LegacyEngine` (baseline parity) and `SafeFallbackEngine`
(always-clean-by-construction).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pytest
from conftest import random_structure

from rna_draw.layout.base import EngineUnavailableError
from rna_draw.layout.fallback import SafeFallbackEngine
from rna_draw.layout.legacy import LegacyEngine
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.render_rna import get_pairmap_from_secstruct

BASELINE_PATH = Path(__file__).parent / "resources" / "layout_baseline.json"
ATOL = 1e-9


def _load_baseline() -> dict[str, Any]:
    """Load the committed coordinate baseline JSON."""
    with BASELINE_PATH.open() as handle:
        return json.load(handle)


def _assert_all_close(actual: list[float], expected: list[float]) -> None:
    """Assert two equal-length float lists match within the fixed tolerance."""
    assert len(actual) == len(expected)
    for got, want in zip(actual, expected):
        assert math.isclose(got, want, rel_tol=0.0, abs_tol=ATOL)


class TestLegacyEngineMatchesBaseline:
    """`LegacyEngine` must reproduce the M1 coordinate baseline exactly."""

    @pytest.mark.parametrize("name", ["hairpin", "dangling", "multiloop", "trna_midsize"])
    def test_matches_committed_baseline(self, name: str) -> None:
        fixture = _load_baseline()[name]
        x, y = LegacyEngine().layout(fixture["ss"])
        _assert_all_close(x, fixture["xarray"])
        _assert_all_close(y, fixture["yarray"])

    def test_name_is_legacy(self) -> None:
        assert LegacyEngine().name == "legacy"


class TestLegacyEngineEdgeCases:
    def test_empty_structure_returns_empty_lists(self) -> None:
        assert LegacyEngine().layout("") == ([], [])

    def test_single_unpaired_nucleotide(self) -> None:
        x, y = LegacyEngine().layout(".")
        assert len(x) == 1
        assert len(y) == 1

    def test_empty_loop_raises_engine_unavailable(self) -> None:
        with pytest.raises(EngineUnavailableError):
            LegacyEngine().layout("()")


class TestSafeFallbackEngine:
    def test_name_is_fallback(self) -> None:
        assert SafeFallbackEngine().name == "fallback"

    def test_empty_structure_returns_empty_lists(self) -> None:
        assert SafeFallbackEngine().layout("") == ([], [])

    @pytest.mark.parametrize(
        "secstruct",
        ["(....)", "((..))((..))", "....", ".", "()", "((((....))))"],
    )
    def test_curated_cases_are_clean(self, secstruct: str) -> None:
        x, y = SafeFallbackEngine().layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map).passed

    @pytest.mark.parametrize("n", [2, 5, 20, 60, 120])
    @pytest.mark.parametrize("seed", range(10))
    def test_random_pseudoknot_free_structures_are_always_clean(self, seed: int, n: int) -> None:
        secstruct = random_structure(seed, n)
        x, y = SafeFallbackEngine().layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map).passed

    def test_clean_under_a_tighter_node_r(self) -> None:
        """A tighter node_r exercises the scaling loop doing >0 doublings."""
        secstruct = "((((....))))((((....))))"
        params = OverlapParams(node_r=20.0)
        x, y = SafeFallbackEngine(params=params).layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map, params).passed
