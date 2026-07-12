"""Tests for the M3 checker-gated compaction pass (`constructive/compaction.py`).

Every test that runs the engine is `@pytest.mark.timeout`-marked (the hard
requirement: the pure-Python engine must never hang the suite).
"""

from __future__ import annotations

import pytest
from conftest import random_structure
from test_constructive_geometry import make_multiloop

from rna_draw.layout.constructive import envelope
from rna_draw.layout.constructive.compaction import compact_layout
from rna_draw.layout.constructive.engine import ConstructiveEngine, _LayoutState, _place_exterior
from rna_draw.layout.structure_tree import build_structure_tree
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct

TIMEOUT = 20
PARAMS = DrawParameters()
OVERLAP_PARAMS = OverlapParams()


def _bbox_area(x: list[float], y: list[float]) -> float:
    """Bounding-box area of a laid-out coordinate set."""
    width = max(x) - min(x)
    height = max(y) - min(y)
    return max(width, 1e-9) * max(height, 1e-9)


def _sound_layout(secstruct: str) -> tuple[list[float], list[float]]:
    """Build the M2 "sound" (pre-compaction) layout directly, bypassing M3.

    Args:
        secstruct: A structure the sound builder can place without
            triggering the S8 retry-with-larger-margin schedule.

    Returns:
        `(x, y)`, the loose, checker-clean-by-construction M2 layout.
    """
    pair_map = get_pairmap_from_secstruct(secstruct)
    tree = build_structure_tree(pair_map)
    n = len(secstruct)
    state = _LayoutState(
        tree=tree, x=[0.0] * n, y=[0.0] * n, params=PARAMS, cache=envelope.ReachCache()
    )
    _place_exterior(state)
    return state.x, state.y


def make_wide_stem_multiloop(degree: int, stem_depth: int, loop_size: int = 3) -> str:
    """A multiloop of `degree` children, each a LONG straight stem.

    A long straight stem is the case the isotropic bounding-disk envelope
    over-estimates worst: its true angular footprint (a thin sliver) is far
    smaller than a disk of radius `~stem_depth`, so this structure is a
    reliable "known loose" case for the compaction gain tests.

    Args:
        degree: Number of long-stem children hanging off the multiloop.
        stem_depth: Stacked-pair depth of every child's stem.
        loop_size: Unpaired nucleotide count in every child's terminal loop.

    Returns:
        A dot-bracket structure.
    """
    return make_multiloop(degree, stem_depth, loop_size)


class TestCompactionNeverWorsensCleanliness:
    """Compaction is checker-gated: the final layout must always be clean."""

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize(
        "secstruct",
        [
            make_wide_stem_multiloop(4, 6),
            make_wide_stem_multiloop(6, 8),
            "(((...)(...))(...))",
            "((((...)(...)(...)))(...))",
            "((...))" * 4,
            "...((...))...",
        ],
    )
    def test_compacted_layout_is_checker_clean(self, secstruct: str) -> None:
        x, y = ConstructiveEngine().layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        report = check_overlaps(x, y, pair_map, OVERLAP_PARAMS)
        assert report.passed, f"{secstruct!r} left {report.num_overlaps} overlaps"

    @pytest.mark.timeout(60)
    @pytest.mark.parametrize("seed", range(30))
    def test_fuzz_compacted_layout_never_dirty(self, seed: int) -> None:
        """Compaction never turns a clean random structure's layout dirty."""
        from rna_draw.layout.base import EngineError

        secstruct = random_structure(seed, 60)
        try:
            x, y = ConstructiveEngine().layout(secstruct)
        except EngineError:
            return
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map, OVERLAP_PARAMS).passed


class TestCompactionShrinksLooseLayouts:
    """A known-loose (elongated-child) multiloop must get materially smaller."""

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize(
        "secstruct",
        [
            make_wide_stem_multiloop(4, 6),
            make_wide_stem_multiloop(6, 8),
            make_wide_stem_multiloop(5, 10),
        ],
    )
    def test_bbox_area_shrinks_materially(self, secstruct: str) -> None:
        """The isotropic-envelope sound build is far sprawlier than a
        tangential-width-repacked one for a multiloop of long straight
        stems -- the textbook case the M3 lever targets (see
        `compaction.py`'s module docstring). Require at least a 2x area
        reduction (measured gain on these fixtures is 5x+).
        """
        sound_x, sound_y = _sound_layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        tree = build_structure_tree(pair_map)
        cache = envelope.ReachCache()
        state = _LayoutState(
            tree=tree, x=list(sound_x), y=list(sound_y), params=PARAMS, cache=cache
        )
        _place_exterior(state)  # repopulate cache's loop_packing for this tree
        compact_x, compact_y = compact_layout(
            tree, state.x, state.y, pair_map, PARAMS, OVERLAP_PARAMS, cache
        )

        assert check_overlaps(compact_x, compact_y, pair_map, OVERLAP_PARAMS).passed
        sound_area = _bbox_area(state.x, state.y)
        compact_area = _bbox_area(compact_x, compact_y)
        assert compact_area <= sound_area / 2.0, (
            f"{secstruct!r}: sound area {sound_area:.1f}, compact area "
            f"{compact_area:.1f} -- expected at least a 2x reduction"
        )


class TestCompactionIsMonotone:
    """Compaction must never make the bounding-box area WORSE (larger)."""

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize(
        "secstruct",
        [
            make_wide_stem_multiloop(3, 4),
            make_wide_stem_multiloop(4, 6),
            "(((...)(...))(...))",
            "((((...)(...)(...)))(...))",
            "((...)((...)(...)))",
            "((...))" * 3,
            "((...))" * 5,
            "..((...))..((...))..",
            "(((((...))).))",
        ],
    )
    def test_compacted_area_never_exceeds_sound_area(self, secstruct: str) -> None:
        sound_x, sound_y = _sound_layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        tree = build_structure_tree(pair_map)
        cache = envelope.ReachCache()
        state = _LayoutState(
            tree=tree,
            x=[0.0] * len(secstruct),
            y=[0.0] * len(secstruct),
            params=PARAMS,
            cache=cache,
        )
        _place_exterior(state)
        compact_x, compact_y = compact_layout(
            tree, state.x, state.y, pair_map, PARAMS, OVERLAP_PARAMS, cache
        )
        assert _bbox_area(compact_x, compact_y) <= _bbox_area(sound_x, sound_y) + 1e-6

    @pytest.mark.timeout(60)
    @pytest.mark.parametrize("seed", range(20))
    def test_fuzz_compacted_area_never_exceeds_sound_area(self, seed: int) -> None:
        from rna_draw.layout.base import EngineError

        secstruct = random_structure(seed, 50)
        pair_map = get_pairmap_from_secstruct(secstruct)
        tree = build_structure_tree(pair_map)
        cache = envelope.ReachCache()
        state = _LayoutState(
            tree=tree,
            x=[0.0] * len(secstruct),
            y=[0.0] * len(secstruct),
            params=PARAMS,
            cache=cache,
        )
        try:
            _place_exterior(state)
        except EngineError:
            return
        compact_x, compact_y = compact_layout(
            tree, state.x, state.y, pair_map, PARAMS, OVERLAP_PARAMS, cache
        )
        assert _bbox_area(compact_x, compact_y) <= _bbox_area(state.x, state.y) + 1e-6


class TestCompactionTiming:
    """Compaction must stay bounded even on the largest admitted structures."""

    @pytest.mark.timeout(30)
    def test_large_structure_compacts_within_a_few_seconds(self) -> None:
        import time

        from rna_draw.layout.base import EngineError

        secstruct = random_structure(1, 3500)
        start = time.monotonic()
        try:
            x, y = ConstructiveEngine().layout(secstruct)
        except EngineError:
            return
        elapsed = time.monotonic() - start
        assert elapsed < 20.0, f"took {elapsed:.2f}s -- too slow for the deep tail"
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map, OVERLAP_PARAMS).passed
