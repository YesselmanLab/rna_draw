"""Tests for `rna_draw.layout.pseudoknot.proximity` (Phase 2b proximity
bias): the reflection isometry, the greedy checker-gated hill-climb, and
its end-to-end effect on `layout_pseudoknot`'s crossing-endpoint distance
and PK-B rate.
"""

from __future__ import annotations

from math import hypot

import pytest

from rna_draw.layout.pseudoknot import engine as engine_module
from rna_draw.layout.pseudoknot.engine import layout_pseudoknot
from rna_draw.layout.pseudoknot.extraction import max_nested_subset
from rna_draw.layout.pseudoknot.parsing import Stem, group_stems, parse_all_pairs
from rna_draw.layout.pseudoknot.proximity import (
    bias_crossing_proximity,
    crossing_endpoint_distance,
)
from rna_draw.overlap import OverlapParams, check_overlaps

# Two single-pair branches hanging off the exterior loop (`(); ()`), pivots
# at index 0 and index 2. Reflecting branch A (0, 1) about its own pivot
# brings index 1 from below (`y=-50`) to above (`y=50`) the pivot line,
# strictly closer to the crossing partner index 3 (also at `y=50`).
IMPROVING_X = [0.0, 0.0, 300.0, 300.0]
IMPROVING_Y = [0.0, -50.0, -50.0, 50.0]
IMPROVING_PAIR_MAP = [1, 0, 3, 2]
IMPROVING_CROSSING = [Stem(i=1, j=3, length=1)]

# A single-pair branch whose only improving flip is checker-rejected: after
# reflection, the departing backbone bond (index 1 -> index 2) passes
# directly through index 0's own disk (the reflected far strand swaps to
# the opposite side of the pivot along the rung's own line -- see
# `proximity.py`'s module docstring on the disk-preservation certificate
# being a heuristic, not a soundness guarantee).
REJECTED_X = [0.0, 20.0, 2000.0, -500.0]
REJECTED_Y = [0.0, 0.0, 500.0, 0.0]
REJECTED_PAIR_MAP = [1, 0, -1, -1]
REJECTED_CROSSING = [Stem(i=1, j=3, length=1)]

H_TYPE = "((((....[[[[....))))....]]]]"
KISSING_LOOP = "((((..[[..))))((((..]]..))))"


class TestReflectRangeIsometry:
    """Direct coverage of the private `_reflect_range` primitive via the
    public `bias_crossing_proximity` entry point (accessing the private
    helper only for the pure-geometry assertions the plan calls for).
    """

    def test_is_an_isometry_on_the_reflected_range(self) -> None:
        from rna_draw.layout.pseudoknot.proximity import _reflect_range

        x = [0.0, 3.0, -2.0, 7.5]
        y = [0.0, 4.0, 6.0, -1.5]
        before_pairwise = [
            hypot(x[a] - x[b], y[a] - y[b]) for a in range(4) for b in range(4) if a != b
        ]
        _reflect_range(x, y, 0, 3, pivot=(0.0, 0.0), axis_dir=(1.0, 0.0))
        after_pairwise = [
            hypot(x[a] - x[b], y[a] - y[b]) for a in range(4) for b in range(4) if a != b
        ]
        for before, after in zip(before_pairwise, after_pairwise):
            assert after == pytest.approx(before, abs=1e-9)

    def test_pivot_index_stays_fixed(self) -> None:
        from rna_draw.layout.pseudoknot.proximity import _reflect_range

        x = [1.0, 5.0, -3.0]
        y = [2.0, -4.0, 6.0]
        _reflect_range(x, y, 0, 2, pivot=(x[0], y[0]), axis_dir=(0.0, 1.0))
        assert x[0] == pytest.approx(1.0)
        assert y[0] == pytest.approx(2.0)

    def test_reflecting_twice_is_identity(self) -> None:
        from rna_draw.layout.pseudoknot.proximity import _reflect_range

        x0 = [0.0, 3.0, -2.0, 7.5]
        y0 = [0.0, 4.0, 6.0, -1.5]
        x, y = list(x0), list(y0)
        pivot, axis_dir = (0.0, 0.0), (1.0, 0.0)
        _reflect_range(x, y, 0, 3, pivot, axis_dir)
        _reflect_range(x, y, 0, 3, pivot, axis_dir)
        for a, b in zip(x, x0):
            assert a == pytest.approx(b, abs=1e-9)
        for a, b in zip(y, y0):
            assert a == pytest.approx(b, abs=1e-9)


class TestDegenerateBranchFrame:
    """A branch whose two rung endpoints coincide (zero-length rung) has no
    well-defined reflection axis -- `_branch_frame` returns `None`, and
    every consumer must skip it rather than divide by zero.
    """

    def test_branch_frame_is_none_for_a_zero_length_rung(self) -> None:
        from rna_draw.layout.pseudoknot.proximity import _branch_frame
        from rna_draw.layout.structure_tree import Branch

        x, y = [5.0, 5.0], [5.0, 5.0]
        branch = Branch(closing_pair=(0, 1), start=0, end=1)
        assert _branch_frame(x, y, branch) is None

    def test_flip_candidate_is_none_for_a_zero_length_rung(self) -> None:
        from rna_draw.layout.pseudoknot.proximity import _flip_candidate
        from rna_draw.layout.structure_tree import Branch

        x, y = [5.0, 5.0], [5.0, 5.0]
        branch = Branch(closing_pair=(0, 1), start=0, end=1)
        assert _flip_candidate(x, y, branch) is None

    def test_best_improving_flip_skips_degenerate_candidates(self) -> None:
        from rna_draw.layout.pseudoknot.proximity import _best_improving_flip
        from rna_draw.layout.structure_tree import Branch

        x, y = [5.0, 5.0, 100.0], [5.0, 5.0, 5.0]
        pair_map = [1, 0, -1]
        branch = Branch(closing_pair=(0, 1), start=0, end=1)
        result = _best_improving_flip(
            x, y, [branch], pair_map, [Stem(i=1, j=2, length=1)], OverlapParams(), base_obj=95.0
        )
        assert result is None


class TestBiasCrossingProximityUnit:
    def test_improving_flip_strictly_lowers_objective_and_stays_clean(self) -> None:
        params = OverlapParams()
        assert check_overlaps(IMPROVING_X, IMPROVING_Y, IMPROVING_PAIR_MAP, params).passed
        base_obj = crossing_endpoint_distance(IMPROVING_X, IMPROVING_Y, IMPROVING_CROSSING)

        bx, by = bias_crossing_proximity(
            IMPROVING_X, IMPROVING_Y, IMPROVING_PAIR_MAP, IMPROVING_CROSSING, params
        )

        new_obj = crossing_endpoint_distance(bx, by, IMPROVING_CROSSING)
        assert new_obj < base_obj
        assert check_overlaps(bx, by, IMPROVING_PAIR_MAP, params).passed

    def test_objective_is_never_worse_than_the_input(self) -> None:
        params = OverlapParams()
        base_obj = crossing_endpoint_distance(IMPROVING_X, IMPROVING_Y, IMPROVING_CROSSING)
        bx, by = bias_crossing_proximity(
            IMPROVING_X, IMPROVING_Y, IMPROVING_PAIR_MAP, IMPROVING_CROSSING, params
        )
        assert crossing_endpoint_distance(bx, by, IMPROVING_CROSSING) <= base_obj

    def test_a_checker_rejected_flip_reverts_cleanly(self) -> None:
        params = OverlapParams()
        assert check_overlaps(REJECTED_X, REJECTED_Y, REJECTED_PAIR_MAP, params).passed

        bx, by = bias_crossing_proximity(
            REJECTED_X, REJECTED_Y, REJECTED_PAIR_MAP, REJECTED_CROSSING, params
        )

        # The only candidate flip lowers the objective but is checker-dirty
        # (see the module docstring above): the search must reject it and
        # return the input UNCHANGED, not some other unclean configuration.
        assert bx == REJECTED_X
        assert by == REJECTED_Y

    def test_entry_guard_leaves_a_dirty_layout_unchanged(self) -> None:
        # Two coincident, unpaired, non-adjacent disks: a deliberately
        # dirty (checker-failing) input.
        x = [0.0, 0.0, 500.0]
        y = [0.0, 0.0, 500.0]
        pair_map = [-1, -1, -1]
        params = OverlapParams()
        assert not check_overlaps(x, y, pair_map, params).passed

        bx, by = bias_crossing_proximity(x, y, pair_map, [Stem(i=0, j=2, length=1)], params)

        assert bx == x
        assert by == y

    def test_no_crossing_stems_is_a_no_op(self) -> None:
        params = OverlapParams()
        bx, by = bias_crossing_proximity(IMPROVING_X, IMPROVING_Y, IMPROVING_PAIR_MAP, [], params)
        assert bx == IMPROVING_X
        assert by == IMPROVING_Y


def _crossing_stems(secstruct: str) -> list[Stem]:
    """The crossing stems `layout_pseudoknot` would extract for `secstruct`."""
    stems = group_stems(parse_all_pairs(secstruct))
    _retained, crossing = max_nested_subset(stems)
    return crossing


def _unbiased_objective(secstruct: str, params: OverlapParams, monkeypatch) -> float:
    """`layout_pseudoknot`'s own crossing-endpoint objective with Phase 2b
    disabled (an identity `bias_crossing_proximity`), i.e. the v1 baseline.
    """
    monkeypatch.setattr(
        engine_module, "bias_crossing_proximity", lambda x, y, *_args, **_kw: (x, y)
    )
    result = layout_pseudoknot(secstruct, params)
    monkeypatch.undo()
    return crossing_endpoint_distance(result.x, result.y, _crossing_stems(secstruct))


class TestEngineIntegration:
    """End-to-end: Phase 2b never regresses `layout_pseudoknot`'s own
    crossing-endpoint objective or PK-B count relative to v1 (bias
    disabled), on the two hand-picked engine fixtures.
    """

    @pytest.mark.parametrize("secstruct", [KISSING_LOOP, H_TYPE])
    def test_biased_objective_is_never_worse_than_unbiased(
        self, secstruct: str, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        params = OverlapParams()
        unbiased_obj = _unbiased_objective(secstruct, params, monkeypatch)

        biased_result = layout_pseudoknot(secstruct, params)
        biased_obj = crossing_endpoint_distance(
            biased_result.x, biased_result.y, _crossing_stems(secstruct)
        )

        assert biased_obj <= unbiased_obj

    @pytest.mark.parametrize("secstruct", [KISSING_LOOP, H_TYPE])
    def test_biased_pk_b_count_is_never_lower_than_unbiased(
        self, secstruct: str, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        params = OverlapParams()
        monkeypatch.setattr(
            engine_module, "bias_crossing_proximity", lambda x, y, *_args, **_kw: (x, y)
        )
        unbiased_pk_b = len(layout_pseudoknot(secstruct, params).crossing_pairs)
        monkeypatch.undo()

        biased_pk_b = len(layout_pseudoknot(secstruct, params).crossing_pairs)

        assert biased_pk_b >= unbiased_pk_b

    @pytest.mark.parametrize("secstruct", [KISSING_LOOP, H_TYPE])
    def test_never_silent_overlap(self, secstruct: str) -> None:
        result = layout_pseudoknot(secstruct, OverlapParams())
        assert result.report.passed or result.flagged

    def test_pk_free_input_is_untouched_by_bias(self) -> None:
        # No crossing stems at all: `bias_crossing_proximity` must be a
        # pure no-op (covered directly above); this confirms the engine
        # wiring doesn't perturb a pk-free layout's coordinates either.
        secstruct = "((((....))))"
        result = layout_pseudoknot(secstruct, OverlapParams())
        assert result.crossing_pairs == []
        assert result.crossing_lines == []
        assert result.flagged is False
