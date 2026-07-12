"""Tests for the provable area-minimization pass (`constructive/area_min.py`).

Covers the one provably-clean DOF (per-branch rotation about its own
attachment pivot) plus the certificate helpers it relies on
(`envelope.branch_disk`/`disks_disjoint`, `area_min._capsule_clear`).

Every test that runs the engine is `@pytest.mark.timeout`-marked (the hard
requirement: the pure-Python engine must never hang the suite).
"""

from __future__ import annotations

import pytest
from conftest import random_structure
from test_constructive_engine import make_degree2_chain
from test_constructive_geometry import make_multiloop

from rna_draw.layout.base import EngineError
from rna_draw.layout.constructive import area_min, envelope
from rna_draw.layout.constructive.engine import ConstructiveEngine, _LayoutState, _place_exterior
from rna_draw.layout.structure_tree import Loop, StructureTree, build_structure_tree, collapse_stem
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


def _sound_layout(
    secstruct: str,
) -> tuple[StructureTree, list[int], envelope.ReachCache, list[float], list[float]]:
    """Build the M2 "sound" (pre-compaction) layout, bypassing M3 and area_min.

    Retries `_MARGIN_SCALES`, exactly mirroring `engine._build_verified`:
    at the default margin scale a structure is not always already
    checker-clean (S8's bounded re-inflation schedule exists precisely
    because of this), so a test helper that skips the retry can hand
    downstream code (compaction, area_min) an UNSOUND input and blame it
    for the resulting dirt.

    Returns:
        `(tree, pair_map, cache, x, y)` for the first margin scale that
        comes back checker-clean.

    Raises:
        EngineError: If no margin scale in `_MARGIN_SCALES` is clean.
    """
    from rna_draw.layout.constructive.engine import _MARGIN_SCALES

    pair_map = get_pairmap_from_secstruct(secstruct)
    tree = build_structure_tree(pair_map)
    n = len(secstruct)
    for margin_scale in _MARGIN_SCALES:
        cache = envelope.ReachCache()
        state = _LayoutState(
            tree=tree,
            x=[0.0] * n,
            y=[0.0] * n,
            params=PARAMS,
            cache=cache,
            margin_scale=margin_scale,
        )
        _place_exterior(state)
        if check_overlaps(state.x, state.y, pair_map, OVERLAP_PARAMS).passed:
            return tree, pair_map, cache, state.x, state.y
    raise EngineError(f"no margin scale produced a clean sound layout for {secstruct!r}")


class TestDisksDisjointCertificate:
    """`envelope.branch_disk`/`disks_disjoint`: rotation's own certificate."""

    def test_disks_disjoint_true_on_sound_build(self) -> None:
        secstruct = make_multiloop(5, 6, 3)
        tree, _pair_map, cache, x, y = _sound_layout(secstruct)
        outer = tree.exterior.children[0]
        _depth, loop = collapse_stem(tree, outer.closing_pair)
        disks = [envelope.branch_disk(child.closing_pair, x, y, cache) for child in loop.children]
        assert len(disks) == 5
        assert envelope.disks_disjoint(disks)

    def test_disks_disjoint_detects_forced_overlap(self) -> None:
        same_center = [((0.0, 0.0), 5.0), ((0.0, 0.0), 5.0)]
        assert not envelope.disks_disjoint(same_center)
        too_close = [((0.0, 0.0), 5.0), ((6.0, 0.0), 5.0)]
        assert not envelope.disks_disjoint(too_close)
        touching_ok = [((0.0, 0.0), 5.0), ((11.0, 0.0), 5.0)]
        assert envelope.disks_disjoint(touching_ok)


class TestCapsuleClear:
    """`area_min._capsule_clear`: the single-capsule-vs-all certificate test."""

    def test_capsule_clear_flags_known_clip(self) -> None:
        x = [0.0, 10.0, 5.0]
        y = [0.0, 0.0, 0.0]
        pair_map = [-1, -1, -1]
        assert not area_min._capsule_clear(x, y, pair_map, OVERLAP_PARAMS, 0)

    def test_capsule_clear_clean_passes(self) -> None:
        x = [0.0, 10.0, 500.0]
        y = [0.0, 0.0, 0.0]
        pair_map = [-1, -1, -1]
        assert area_min._capsule_clear(x, y, pair_map, OVERLAP_PARAMS, 0)

    def test_capsule_clear_on_a_real_sound_layout(self) -> None:
        secstruct = make_multiloop(4, 6, 3)
        _tree, pair_map, _cache, x, y = _sound_layout(secstruct)
        for k in range(len(x) - 1):
            assert area_min._capsule_clear(x, y, pair_map, OVERLAP_PARAMS, k)


class TestRotateBranches:
    """DOF 1: per-branch rotation about its own attachment pivot."""

    @pytest.mark.timeout(TIMEOUT)
    def test_orientation_reduces_bbox_on_elongated_branch(self) -> None:
        """A multiloop of long straight stems is the textbook case an
        isotropic-disk sound build over-orients: reorienting a child so its
        far tip lands closer to the rest of the structure should shrink the
        bbox (never grow it)."""
        secstruct = make_multiloop(5, 10, 3)
        tree, pair_map, cache, x, y = _sound_layout(secstruct)
        before = _bbox_area(x, y)
        rx, ry = area_min.rotate_branches(tree, x, y, pair_map, PARAMS, OVERLAP_PARAMS, cache)
        report = check_overlaps(rx, ry, pair_map, OVERLAP_PARAMS)
        assert report.passed
        assert _bbox_area(rx, ry) <= before + 1e-6

    def test_rejects_rotation_that_would_clip_departing_capsule(self) -> None:
        """A hand-built collision, at a realistic scale (a branch's own
        reach is always much larger than the checker's disk radius):
        rotating index range `[0, 1]` 90 degrees about the origin would
        swing index 1 from `(50, 0)` to `(0, -50)`, sending the departing
        capsule `(1, 2)` (to a fixed point at `(0, -100)`) straight through
        an unrelated fixed disk at `(0, -70)`.
        """
        pair_map = [-1, -1, -1, -1]
        tree = build_structure_tree(pair_map)
        state = area_min._AreaMinState(
            tree=tree,
            x=[0.0, 50.0, 0.0, 0.0],
            y=[0.0, 0.0, -100.0, -70.0],
            pair_map=pair_map,
            params=PARAMS,
            overlap_params=OVERLAP_PARAMS,
            margin_scale=1.0,
            deadline=None,
        )
        saved_x, saved_y = list(state.x), list(state.y)
        accepted = area_min._apply_and_verify_rotation(
            state, lo=0, hi=1, j_b=1, pivot=(0.0, 0.0), old_dir=(1.0, 0.0), new_dir=(0.0, -1.0)
        )
        assert not accepted
        assert state.x == saved_x
        assert state.y == saved_y

    def test_accepts_rotation_when_departing_capsule_stays_clear(self) -> None:
        """The same shape of rotation as above, at a realistic scale (a
        branch's own reach is always much larger than the checker's disk
        radius) with the obstacle moved far away."""
        pair_map = [-1, -1, -1, -1]
        tree = build_structure_tree(pair_map)
        state = area_min._AreaMinState(
            tree=tree,
            x=[0.0, 50.0, 0.0, 1000.0],
            y=[0.0, 0.0, -100.0, 1000.0],
            pair_map=pair_map,
            params=PARAMS,
            overlap_params=OVERLAP_PARAMS,
            margin_scale=1.0,
            deadline=None,
        )
        accepted = area_min._apply_and_verify_rotation(
            state, lo=0, hi=1, j_b=1, pivot=(0.0, 0.0), old_dir=(1.0, 0.0), new_dir=(0.0, -1.0)
        )
        assert accepted
        assert state.x[1] == pytest.approx(0.0, abs=1e-9)
        assert state.y[1] == pytest.approx(-50.0, abs=1e-9)

    @pytest.mark.timeout(TIMEOUT)
    def test_pinned_degree2_left_straight(self) -> None:
        """A degree-2 chain past `_DEGREE2_CHAIN_LENGTH_FLOOR` is pinned
        collinear by the sound build (M2b(i)); `rotate_branches` must not
        touch those loops' children."""
        secstruct = make_degree2_chain(10)
        tree, pair_map, cache, x, y = _sound_layout(secstruct)
        pinned_loops = [
            loop
            for loop in tree.loops
            if loop.closing_pair is not None and cache.degree2_pinned.get(loop.closing_pair, False)
        ]
        assert pinned_loops, "fixture should contain at least one pinned degree-2 loop"
        for loop in pinned_loops:
            assert not area_min._loop_allows_rotation(cache, loop)

        # `rotate_branches` may still rigidly move a pinned loop's whole
        # range as part of rotating some OUTER (unpinned/exterior) branch
        # that contains it -- a rigid transform preserves collinearity, so
        # the right invariant is the dominant child's axis staying PARALLEL
        # to the loop's own axis (the pinned-at-pi, straight-continuation
        # relationship `envelope._degree2_packing` establishes), not that
        # absolute coordinates never move.
        rx, ry = area_min.rotate_branches(tree, x, y, pair_map, PARAMS, OVERLAP_PARAMS, cache)
        for loop in pinned_loops:
            assert _dom_axis_parallel_to_loop_axis(tree, loop, cache, rx, ry)

    @pytest.mark.timeout(TIMEOUT)
    def test_bulge_child_left_straight(self) -> None:
        """A bulge/interior loop's sole child must keep continuing the
        parent stem's own axis (M2b(i))."""
        secstruct = "((.(...)))"
        tree, _pair_map, cache, _x, _y = _sound_layout(secstruct)
        outer = tree.exterior.children[0]
        _depth, loop = collapse_stem(tree, outer.closing_pair)
        assert len(loop.children) == 1
        assert not area_min._loop_allows_rotation(cache, loop)

    @pytest.mark.timeout(60)
    @pytest.mark.parametrize("seed", range(20))
    def test_fuzz_rotation_never_introduces_overlap(self, seed: int) -> None:
        secstruct = random_structure(seed, 150)
        try:
            tree, pair_map, cache, x, y = _sound_layout(secstruct)
        except EngineError:
            return
        rx, ry = area_min.rotate_branches(tree, x, y, pair_map, PARAMS, OVERLAP_PARAMS, cache)
        assert check_overlaps(rx, ry, pair_map, OVERLAP_PARAMS).passed

    @pytest.mark.timeout(60)
    @pytest.mark.parametrize("seed", range(20))
    def test_fuzz_rotation_never_grows_bbox(self, seed: int) -> None:
        secstruct = random_structure(seed, 150)
        try:
            tree, pair_map, cache, x, y = _sound_layout(secstruct)
        except EngineError:
            return
        before = _bbox_area(x, y)
        rx, ry = area_min.rotate_branches(tree, x, y, pair_map, PARAMS, OVERLAP_PARAMS, cache)
        assert _bbox_area(rx, ry) <= before + 1e-6


def _dom_axis_parallel_to_loop_axis(
    tree: StructureTree,
    loop: Loop,
    cache: envelope.ReachCache,
    x: list[float],
    y: list[float],
) -> bool:
    """Whether a pinned degree-2 loop's dominant child still continues
    straight through (its own axis parallel to the loop's own axis)."""
    from rna_draw.layout.constructive import compaction
    from rna_draw.layout.constructive import envelope as env

    dom, _side = env._dominant_and_side(tree, loop, PARAMS, cache, 1.0)
    assert loop.closing_pair is not None  # caller only passes pinned (closing-paired) loops
    i, j = loop.closing_pair
    loop_axis = compaction._axis_dir_from_rung((x[i], y[i]), (x[j], y[j]), PARAMS.PAIR_SPACE)
    di, dj = dom.closing_pair
    dom_axis = compaction._axis_dir_from_rung((x[di], y[di]), (x[dj], y[dj]), PARAMS.PAIR_SPACE)
    cos_angle = loop_axis[0] * dom_axis[0] + loop_axis[1] * dom_axis[1]
    return cos_angle == pytest.approx(1.0, abs=1e-6)


class TestEngineIntegration:
    """`ConstructiveEngine` end-to-end, with area_min wired in."""

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize(
        "secstruct",
        [
            make_multiloop(4, 6, 3),
            make_multiloop(6, 8, 3),
            make_degree2_chain(10),
            "(((...)(...))(...))",
            "((...))" * 4,
        ],
    )
    def test_engine_returns_clean_layout(self, secstruct: str) -> None:
        x, y = ConstructiveEngine().layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        report = check_overlaps(x, y, pair_map, OVERLAP_PARAMS)
        assert report.passed, f"{secstruct!r} left {report.num_overlaps} overlaps"

    @pytest.mark.timeout(TIMEOUT)
    def test_engine_never_worse_than_compaction_alone(self) -> None:
        """The area-min pass's own `_smaller_bbox` safety net: the engine's
        final rendered area must never exceed what plain M3 compaction
        (the pre-area_min baseline) alone would have produced."""
        from rna_draw.layout.constructive import compaction
        from rna_draw.layout.constructive.engine import _rescale_or_keep_clean

        secstruct = make_multiloop(5, 10, 3)
        tree, pair_map, cache, x, y = _sound_layout(secstruct)
        cx, cy = compaction.compact_layout(tree, x, y, pair_map, PARAMS, OVERLAP_PARAMS, cache)
        if check_overlaps(cx, cy, pair_map, OVERLAP_PARAMS).passed:
            baseline_x, baseline_y = cx, cy
        else:
            baseline_x, baseline_y = x, y
        baseline_x, baseline_y = _rescale_or_keep_clean(
            baseline_x, baseline_y, pair_map, PARAMS.PRIMARY_SPACE
        )
        baseline_area = _bbox_area(baseline_x, baseline_y)

        engine_x, engine_y = ConstructiveEngine().layout(secstruct)
        assert _bbox_area(engine_x, engine_y) <= baseline_area + 1e-6

    @pytest.mark.timeout(60)
    @pytest.mark.parametrize("seed", range(15))
    def test_fuzz_engine_never_silent_overlap_area_min(self, seed: int) -> None:
        secstruct = random_structure(seed, 200)
        try:
            x, y = ConstructiveEngine().layout(secstruct)
        except EngineError:
            return
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map, OVERLAP_PARAMS).passed
