"""Tests for the `constructive` layout engine's M1 scope (see the
constructive-engine runbook, milestone M1): the straight-helix ladder
placer, the round-loop exact angular-interval packer, and their
composition into a single multiloop of hairpins -- proven checker-clean
BY CONSTRUCTION on 50 synthetic structures.
"""

from __future__ import annotations

import math
import statistics

import pytest

from rna_draw.layout.base import EngineError
from rna_draw.layout.constructive import ConstructiveEngine
from rna_draw.layout.constructive.geometry_helpers import (
    pack_loop_angles,
    place_stem,
)
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct

PARAMS = DrawParameters()
COLLINEAR_ATOL = 1e-9
SPACING_ATOL = 1e-6


def make_isolated_stem(depth: int, loop_size: int = 3) -> str:
    """A `depth`-pair stem closing a plain `loop_size`-nt hairpin loop."""
    return "(" * depth + "." * loop_size + ")" * depth


def make_multiloop(degree: int, hairpin_depth: int, loop_size: int = 3) -> str:
    """A single multiloop of `degree` hairpins, each `hairpin_depth` pairs deep.

    Args:
        degree: Number of hairpin children hanging off the multiloop.
        hairpin_depth: Stacked-pair depth of every child's stem.
        loop_size: Unpaired nucleotide count in every child's terminal loop.

    Returns:
        A dot-bracket structure: one closing stem (a single pair) wrapping
        `degree` consecutive hairpins.
    """
    hairpin = "(" * hairpin_depth + "." * loop_size + ")" * hairpin_depth
    return "(" + hairpin * degree + ")"


def _cross(o: tuple[float, float], a: tuple[float, float], b: tuple[float, float]) -> float:
    """2D cross product of `oa` and `ob`, used to test three-point collinearity."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


class TestPlaceStem:
    """`place_stem` must build an exactly-straight, constantly-spaced ladder."""

    @pytest.mark.parametrize("depth", range(2, 11))
    def test_strands_are_collinear(self, depth: int) -> None:
        ladder = place_stem(depth, (1.0, -2.0), (0.6, 0.8), PARAMS.PRIMARY_SPACE, PARAMS.PAIR_SPACE)
        for strand in (ladder.strand_a, ladder.strand_b):
            for k in range(2, len(strand)):
                assert abs(_cross(strand[0], strand[1], strand[k])) < COLLINEAR_ATOL

    @pytest.mark.parametrize("depth", range(2, 11))
    def test_constant_rail_separation(self, depth: int) -> None:
        ladder = place_stem(depth, (0.0, 0.0), (1.0, 0.0), PARAMS.PRIMARY_SPACE, PARAMS.PAIR_SPACE)
        for a, b in zip(ladder.strand_a, ladder.strand_b):
            dist = math.hypot(a[0] - b[0], a[1] - b[1])
            assert math.isclose(dist, PARAMS.PAIR_SPACE, rel_tol=0.0, abs_tol=SPACING_ATOL)

    @pytest.mark.parametrize("depth", range(2, 11))
    def test_constant_rise_along_axis(self, depth: int) -> None:
        ladder = place_stem(depth, (0.0, 0.0), (1.0, 0.0), PARAMS.PRIMARY_SPACE, PARAMS.PAIR_SPACE)
        for strand in (ladder.strand_a, ladder.strand_b):
            for k in range(1, depth):
                step = math.hypot(strand[k][0] - strand[k - 1][0], strand[k][1] - strand[k - 1][1])
                assert math.isclose(step, PARAMS.PRIMARY_SPACE, rel_tol=0.0, abs_tol=SPACING_ATOL)

    def test_rejects_nonpositive_depth(self) -> None:
        with pytest.raises(ValueError):
            place_stem(0, (0.0, 0.0), (1.0, 0.0), PARAMS.PRIMARY_SPACE, PARAMS.PAIR_SPACE)

    @pytest.mark.parametrize("depth", range(2, 11))
    def test_isolated_stem_is_checker_clean(self, depth: int) -> None:
        secstruct = make_isolated_stem(depth)
        x, y = ConstructiveEngine().layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map, OverlapParams()).passed


class TestPackLoopAngles:
    """`pack_loop_angles` must exactly pack disjoint, ordered intervals."""

    @pytest.mark.parametrize("k", range(3, 13))
    def test_intervals_disjoint_ordered_and_fit(self, k: int) -> None:
        half_widths = [12.0] * k
        reserved = 15.0
        packing = pack_loop_angles(half_widths, reserved, radius_floor=1.0, radius_step=0.5)

        assert packing.angles == sorted(packing.angles)
        half_angles = [math.asin(hw / packing.radius) for hw in half_widths]
        reserved_half_angle = math.asin(reserved / packing.radius)
        lower_edges = [a - h for a, h in zip(packing.angles, half_angles)]
        upper_edges = [a + h for a, h in zip(packing.angles, half_angles)]
        assert lower_edges[0] >= reserved_half_angle - 1e-9
        for upper, next_lower in zip(upper_edges, lower_edges[1:]):
            assert next_lower >= upper - 1e-9
        assert upper_edges[-1] <= 2 * math.pi - reserved_half_angle + 1e-9

        total_required = 2 * reserved_half_angle + sum(2 * h for h in half_angles)
        assert total_required <= 2 * math.pi + 1e-9

    @pytest.mark.parametrize("k", range(3, 13))
    def test_radius_is_smallest_fitting_grid_point(self, k: int) -> None:
        half_widths = [10.0] * k
        reserved = 10.0
        radius_floor, radius_step = 1.0, 0.25
        packing = pack_loop_angles(half_widths, reserved, radius_floor, radius_step)

        n_steps = round((packing.radius - radius_floor) / radius_step)
        assert math.isclose(packing.radius, radius_floor + n_steps * radius_step, abs_tol=1e-9)
        if n_steps > 0:
            smaller_radius = radius_floor + (n_steps - 1) * radius_step
            total = 2 * math.asin(reserved / smaller_radius)
            total += sum(2 * math.asin(hw / smaller_radius) for hw in half_widths)
            assert total > 2 * math.pi

    def test_larger_children_need_non_decreasing_radius(self) -> None:
        radii = []
        for half_width in (5.0, 10.0, 20.0, 40.0):
            packing = pack_loop_angles(
                [half_width] * 6, reserved_half_width=5.0, radius_floor=1.0, radius_step=0.5
            )
            radii.append(packing.radius)
        assert radii == sorted(radii)

    def test_rejects_nonpositive_radius_step(self) -> None:
        with pytest.raises(ValueError):
            pack_loop_angles([10.0], 5.0, radius_floor=1.0, radius_step=0.0)


class TestConstructiveEngineSingleMultiloopOfHairpins:
    """M1's Checkpoint A deliverable: clean-by-construction on 50 synthetics
    spanning `degree in 3..12` x `hairpin_depth in 2..6` (the runbook's
    required >= 40), each with a varied terminal-loop size.
    """

    @pytest.mark.parametrize("degree", range(3, 13))
    @pytest.mark.parametrize("hairpin_depth", range(2, 7))
    def test_synthetic_multiloop_is_checker_clean(self, degree: int, hairpin_depth: int) -> None:
        loop_size = 3 + (degree + hairpin_depth) % 4
        secstruct = make_multiloop(degree, hairpin_depth, loop_size)
        pair_map = get_pairmap_from_secstruct(secstruct)

        x, y = ConstructiveEngine().layout(secstruct)

        report = check_overlaps(x, y, pair_map, OverlapParams())
        assert report.passed, f"{secstruct!r} left {report.num_overlaps} overlaps"

    @pytest.mark.parametrize("degree", range(3, 13))
    @pytest.mark.parametrize("hairpin_depth", range(2, 7))
    def test_median_backbone_step_within_5_percent(self, degree: int, hairpin_depth: int) -> None:
        secstruct = make_multiloop(degree, hairpin_depth)
        x, y = ConstructiveEngine().layout(secstruct)

        steps = [math.hypot(x[i + 1] - x[i], y[i + 1] - y[i]) for i in range(len(x) - 1)]
        median_step = statistics.median(steps)

        assert abs(median_step - PARAMS.PRIMARY_SPACE) / PARAMS.PRIMARY_SPACE <= 0.05


class TestConstructiveEngineScope:
    """M1's explicit scope guardrails: unsupported shapes raise `EngineError`,
    never a silently returned overlap.
    """

    def test_name_is_constructive(self) -> None:
        assert ConstructiveEngine().name == "constructive"

    def test_empty_structure_returns_empty_lists(self) -> None:
        assert ConstructiveEngine().layout("") == ([], [])

    def test_single_unpaired_nucleotide(self) -> None:
        assert ConstructiveEngine().layout(".") == ([0.0], [0.0])

    def test_pseudoknot_raises(self) -> None:
        with pytest.raises(EngineError):
            ConstructiveEngine().layout("([)]")

    def test_dangling_tail_raises(self) -> None:
        with pytest.raises(EngineError):
            ConstructiveEngine().layout(".((...))")

    def test_multiple_top_level_branches_raises(self) -> None:
        with pytest.raises(EngineError):
            ConstructiveEngine().layout("((...))((...))")

    def test_all_unpaired_raises(self) -> None:
        with pytest.raises(EngineError):
            ConstructiveEngine().layout("....")

    def test_nested_multiloop_child_raises(self) -> None:
        with pytest.raises(EngineError):
            ConstructiveEngine().layout("((((...))(...)).(...))")

    def test_plain_hairpin_is_clean(self) -> None:
        secstruct = "(((...)))"
        x, y = ConstructiveEngine().layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map, OverlapParams()).passed

    def test_shallow_depth_one_children_are_clean(self) -> None:
        secstruct = "((...)(...))"
        x, y = ConstructiveEngine().layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map, OverlapParams()).passed
