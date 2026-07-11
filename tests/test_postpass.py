"""Tests for `rna_draw.layout.postpass`.

Three properties matter most here, in order of importance:

1. **Rigidity** (`TestRigidRangeTransforms`): `rotate_range`/`translate_range`
   preserve every pairwise distance within the moved slice and leave
   everything outside it untouched -- this is what keeps helices straight
   and loops round; bending would violate the whole point of a *rigid*
   post-pass.
2. **Monotone / never-worse** (`TestRemoveOverlapsMonotone`): the guarantee
   that makes the pass safe to run unconditionally on any layout --
   `report_after.num_overlaps <= report_before.num_overlaps` for every
   input, checked on synthetic layouts and on real dirty structures laid
   out by `EscalatingClearanceEngine`.
3. **Efficacy** (`TestRemoveOverlapsEfficacy`): a softer property -- at
   least the easy near-clean cases actually get fixed.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import pytest

from rna_draw.geometry import PrimitiveId
from rna_draw.layout.legacy import LegacyEngine
from rna_draw.layout.postpass import (
    POSTPASS_PARAMS,
    PostPassConfig,
    _alt_witness_nucleotides,
    _resolve_witness,
    inflate_loop,
    remove_overlaps,
    rotate_range,
    translate_range,
    witness_nucleotides,
)
from rna_draw.layout.structure_tree import build_structure_tree, loop_center
from rna_draw.overlap import OverlapKind, Witness, check_overlaps
from rna_draw.render_rna import get_pairmap_from_secstruct

WORST_SET_JSON = Path(__file__).parent.parent / "benchmarks" / "worst_set.json"
TWO_HAIRPIN_MULTILOOP = "((..)((..))..)"
CROWDED_HAIRPIN = "((........))"


def _witness(kind: OverlapKind, a: PrimitiveId, b: PrimitiveId) -> Witness:
    """Build a `Witness` with placeholder geometry (only the ids matter here)."""
    return Witness(kind=kind, id_a=a, id_b=b, separation=1.0, overlap_depth=1.0)


class TestWitnessNucleotides:
    """`witness_nucleotides` must recover the documented representative
    index per the reuse map: the disk's own index, a backbone capsule's
    lower endpoint, or a pair capsule's lower endpoint.
    """

    def test_nt_kind_uses_its_own_index(self) -> None:
        w = _witness(OverlapKind.DISK_DISK, PrimitiveId("nt", 2), PrimitiveId("nt", 7))
        assert witness_nucleotides(w, [-1] * 10) == (2, 7)

    def test_bb_kind_uses_the_lower_backbone_endpoint(self) -> None:
        w = _witness(OverlapKind.CAPSULE_CAPSULE, PrimitiveId("bb", 4), PrimitiveId("bb", 9))
        assert witness_nucleotides(w, [-1] * 11) == (4, 9)

    def test_pair_kind_uses_the_lower_of_the_two_partners(self) -> None:
        pair_map = [-1] * 12
        pair_map[5], pair_map[10] = 10, 5
        pair_map[6], pair_map[9] = 9, 6
        w = _witness(OverlapKind.CAPSULE_CAPSULE, PrimitiveId("pair", 5), PrimitiveId("pair", 6))
        assert witness_nucleotides(w, pair_map) == (5, 6)


class TestAltWitnessNucleotides:
    def test_nt_kind_alt_equals_the_only_endpoint(self) -> None:
        w = _witness(OverlapKind.DISK_DISK, PrimitiveId("nt", 2), PrimitiveId("nt", 7))
        assert _alt_witness_nucleotides(w, [-1] * 10) == (2, 7)

    def test_bb_kind_alt_uses_the_upper_backbone_endpoint(self) -> None:
        w = _witness(OverlapKind.CAPSULE_CAPSULE, PrimitiveId("bb", 4), PrimitiveId("bb", 9))
        assert _alt_witness_nucleotides(w, [-1] * 11) == (5, 10)

    def test_pair_kind_alt_uses_the_upper_partner(self) -> None:
        pair_map = [-1] * 12
        pair_map[5], pair_map[10] = 10, 5
        pair_map[6], pair_map[9] = 9, 6
        w = _witness(OverlapKind.CAPSULE_CAPSULE, PrimitiveId("pair", 5), PrimitiveId("pair", 6))
        assert _alt_witness_nucleotides(w, pair_map) == (10, 9)


class TestResolveWitnessR6EndpointRetry:
    """R6: if the lower-index representative resolves to no branch on
    EITHER side, retry with each side's other capsule endpoint.
    """

    def test_retry_finds_a_branch_when_lower_endpoints_do_not(self) -> None:
        pair_map = get_pairmap_from_secstruct(TWO_HAIRPIN_MULTILOOP)
        tree = build_structure_tree(pair_map)
        # bb(0) spans nt(0, 1): nt0 is the multiloop's own boundary member.
        # bb(11) spans nt(11, 12): both bare multiloop members. The lower
        # reps (0, 11) are both bare loop members (branch is None on both
        # sides); the upper reps (1, 12) resolve nt1 into branch (1, 4).
        w = _witness(OverlapKind.CAPSULE_CAPSULE, PrimitiveId("bb", 0), PrimitiveId("bb", 11))
        loop, branch_a, ka, branch_b, kb = _resolve_witness(tree, w, pair_map)
        assert loop.closing_pair == (0, 13)
        assert branch_a is not None
        assert branch_a.closing_pair == (1, 4)
        assert ka == 1
        assert branch_b is None
        assert kb == 12

    def test_no_branch_anywhere_stays_unresolved(self) -> None:
        # Two bare members of the same hairpin loop: no branch on either
        # side, and the retry endpoints do not change that (genuinely
        # nothing to move -- overlap is *within* one rigid branch).
        pair_map = get_pairmap_from_secstruct(TWO_HAIRPIN_MULTILOOP)
        tree = build_structure_tree(pair_map)
        w = _witness(OverlapKind.DISK_DISK, PrimitiveId("nt", 2), PrimitiveId("nt", 3))
        loop, branch_a, _ka, branch_b, _kb = _resolve_witness(tree, w, pair_map)
        assert loop.closing_pair == (1, 4)
        assert branch_a is None
        assert branch_b is None


class TestRigidRangeTransforms:
    """Exact rigidity is the load-bearing property: it is what keeps
    helices straight and loops round as a pure corollary (a transform
    that preserves every pairwise distance within a slice cannot bend it).
    """

    def _coords(self) -> tuple[list[float], list[float]]:
        x = [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0]
        y = [0.0, -2.0, 4.0, 1.0, -5.0, 6.0, 2.0, -1.0, 3.0, 0.0]
        return x, y

    def test_rotate_preserves_internal_pairwise_distances(self) -> None:
        x, y = self._coords()
        nx, ny = rotate_range(x, y, 3, 7, cx=10.0, cy=-1.5, angle_rad=0.83)
        for i in range(3, 8):
            for j in range(3, 8):
                before = math.hypot(x[i] - x[j], y[i] - y[j])
                after = math.hypot(nx[i] - nx[j], ny[i] - ny[j])
                assert after == pytest.approx(before, abs=1e-9)

    def test_rotate_leaves_outside_slice_untouched(self) -> None:
        x, y = self._coords()
        nx, ny = rotate_range(x, y, 3, 7, cx=10.0, cy=-1.5, angle_rad=0.83)
        assert nx[:3] == x[:3]
        assert ny[:3] == y[:3]
        assert nx[8:] == x[8:]
        assert ny[8:] == y[8:]

    def test_translate_preserves_internal_pairwise_distances(self) -> None:
        x, y = self._coords()
        nx, ny = translate_range(x, y, 2, 6, dx=5.0, dy=-3.0)
        for i in range(2, 7):
            for j in range(2, 7):
                before = math.hypot(x[i] - x[j], y[i] - y[j])
                after = math.hypot(nx[i] - nx[j], ny[i] - ny[j])
                assert after == pytest.approx(before, abs=1e-9)

    def test_translate_leaves_outside_slice_untouched(self) -> None:
        x, y = self._coords()
        nx, ny = translate_range(x, y, 2, 6, dx=5.0, dy=-3.0)
        assert nx[:2] == x[:2]
        assert ny[:2] == y[:2]
        assert nx[7:] == x[7:]
        assert ny[7:] == y[7:]

    def test_translate_shifts_slice_by_exactly_dx_dy(self) -> None:
        x, y = self._coords()
        nx, ny = translate_range(x, y, 2, 6, dx=5.0, dy=-3.0)
        for i in range(2, 7):
            assert nx[i] == pytest.approx(x[i] + 5.0)
            assert ny[i] == pytest.approx(y[i] - 3.0)


def _nudged_multiloop(fraction: float) -> tuple[list[float], list[float], list[int]]:
    """A near-clean dirty layout: `TWO_HAIRPIN_MULTILOOP` laid out cleanly,
    then branch (5, 10) nudged toward branch (1, 4) by `fraction` of the
    distance between their centroids.
    """
    pair_map = get_pairmap_from_secstruct(TWO_HAIRPIN_MULTILOOP)
    x, y = LegacyEngine().layout(TWO_HAIRPIN_MULTILOOP)
    ax = sum(x[i] for i in range(1, 5)) / 4
    ay = sum(y[i] for i in range(1, 5)) / 4
    bx = sum(x[i] for i in range(5, 11)) / 6
    by = sum(y[i] for i in range(5, 11)) / 6
    nx, ny = translate_range(x, y, 5, 10, (ax - bx) * fraction, (ay - by) * fraction)
    return nx, ny, pair_map


def _crowded_hairpin_loop(shrink: float) -> tuple[list[float], list[float], list[int]]:
    """`CROWDED_HAIRPIN` laid out cleanly, then its hairpin loop shrunk by
    `shrink` (a factor < 1 fed to `inflate_loop`, reusing it in reverse) so
    its 8 unpaired nts crowd together. The hairpin has no child branches,
    so `_resolve_witness` finds no branch on either side of any resulting
    witness -- this is fixable ONLY by `inflate_loop`, not by rotating or
    translating a branch, and so isolates the inflation move's efficacy.
    """
    pair_map = get_pairmap_from_secstruct(CROWDED_HAIRPIN)
    x, y = LegacyEngine().layout(CROWDED_HAIRPIN)
    tree = build_structure_tree(pair_map)
    hairpin_loop = tree.loop_by_closing_pair[(1, 10)]
    nx, ny = inflate_loop(hairpin_loop, x, y, pair_map, shrink)
    return nx, ny, pair_map


def _scattered_exterior_clash() -> tuple[list[float], list[float], list[int]]:
    """Two bare, unpaired exterior nucleotides forced to overlap: no branch
    exists on either side, so this is a case `remove_overlaps` cannot fix
    (a genuine stuck case), included to prove monotonicity holds trivially
    even when zero moves are available.
    """
    x = [0.0, 0.05, 40.0, 80.0]
    y = [0.0, 0.0, 0.0, 0.0]
    pair_map = [-1, -1, -1, -1]
    return x, y, pair_map


def _load_worst_structures(names: list[str]) -> list[str]:
    """Look up dot-bracket structures by `.dbn` filename in the worst set."""
    entries = json.loads(WORST_SET_JSON.read_text())
    by_name = {e["name"]: e["structure"] for e in entries}
    return [by_name[name] for name in names]


class TestRemoveOverlapsMonotone:
    """The hard guarantee: never return a layout with more overlaps than
    the input had, for any input.
    """

    def test_clean_input_returned_unchanged(self) -> None:
        x = [0.0, 20.0, 40.0]
        y = [0.0, 0.0, 0.0]
        pair_map = [-1, -1, -1]
        result = remove_overlaps(x, y, pair_map)
        assert result.x == x
        assert result.y == y
        assert result.moves_applied == 0
        assert result.report_after.num_overlaps == 0
        assert result.report_before is result.report_after

    def test_fewer_than_two_nucleotides_returned_unchanged(self) -> None:
        result = remove_overlaps([0.0], [0.0], [-1])
        assert result.x == [0.0]
        assert result.moves_applied == 0

    def test_empty_input_returned_unchanged(self) -> None:
        result = remove_overlaps([], [], [])
        assert result.x == []
        assert result.report_before.num_overlaps == 0

    @pytest.mark.parametrize("fraction", [0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.9])
    def test_synthetic_nudged_multiloop_never_worsens(self, fraction: float) -> None:
        x, y, pair_map = _nudged_multiloop(fraction)
        before = check_overlaps(x, y, pair_map, POSTPASS_PARAMS).num_overlaps
        result = remove_overlaps(x, y, pair_map)
        assert result.report_before.num_overlaps == before
        assert result.report_after.num_overlaps <= before

    @pytest.mark.parametrize("shrink", [0.3, 0.5, 0.7, 0.9])
    def test_crowded_hairpin_inflation_never_worsens(self, shrink: float) -> None:
        x, y, pair_map = _crowded_hairpin_loop(shrink)
        before = check_overlaps(x, y, pair_map, POSTPASS_PARAMS).num_overlaps
        result = remove_overlaps(x, y, pair_map)
        assert result.report_before.num_overlaps == before
        assert result.report_after.num_overlaps <= before

    def test_unfixable_exterior_clash_stays_monotone(self) -> None:
        x, y, pair_map = _scattered_exterior_clash()
        before = check_overlaps(x, y, pair_map, POSTPASS_PARAMS).num_overlaps
        assert before > 0
        result = remove_overlaps(x, y, pair_map)
        assert result.report_after.num_overlaps <= before

    @pytest.mark.parametrize("name", ["bpRNA_CRW_11829.dbn", "bpRNA_RFAM_15398.dbn"])
    def test_real_worst_set_structures_never_worsen(self, name: str) -> None:
        pytest.importorskip("rna_draw._vienna_layout")
        from benchmarks.engines import EscalatingClearanceEngine

        (structure,) = _load_worst_structures([name])
        x, y = EscalatingClearanceEngine().layout(structure)
        pair_map = get_pairmap_from_secstruct(structure)
        before = check_overlaps(x, y, pair_map, POSTPASS_PARAMS).num_overlaps
        assert before > 0, "fixture expected to be dirty before the post-pass"

        result = remove_overlaps(x, y, pair_map)
        assert result.report_after.num_overlaps <= before

    @pytest.mark.parametrize("time_budget_s", [0.0, 1e-6])
    def test_tiny_time_budget_never_worsens_and_returns_quickly(
        self, time_budget_s: float
    ) -> None:
        # A tiny budget should cut the move loop off almost immediately --
        # the monotone guarantee holds regardless (see docstring): every
        # accepted move already strictly reduced the count, so best-so-far
        # is at most the input's overlap count even if the loop exits after
        # zero or one moves.
        x, y, pair_map = _nudged_multiloop(0.5)
        before = check_overlaps(x, y, pair_map, POSTPASS_PARAMS).num_overlaps
        assert before > 0, "fixture expected to be dirty before the post-pass"

        start = time.monotonic()
        result = remove_overlaps(x, y, pair_map, PostPassConfig(time_budget_s=time_budget_s))
        elapsed = time.monotonic() - start

        assert result.report_after.num_overlaps <= before
        assert elapsed < 2.0

    def test_time_budget_none_reproduces_prior_unbounded_behavior(self) -> None:
        x, y, pair_map = _nudged_multiloop(0.15)
        unbounded = remove_overlaps(x, y, pair_map, PostPassConfig(time_budget_s=None))
        default = remove_overlaps(x, y, pair_map, PostPassConfig())
        assert unbounded.x == default.x
        assert unbounded.y == default.y
        assert unbounded.moves_applied == default.moves_applied


class TestRemoveOverlapsEfficacy:
    """Softer property: at least the easy near-clean cases get fixed."""

    def test_small_nudge_is_driven_to_zero(self) -> None:
        x, y, pair_map = _nudged_multiloop(0.15)
        before = check_overlaps(x, y, pair_map, POSTPASS_PARAMS).num_overlaps
        assert 0 < before <= 6, "fixture should be near-clean, not deeply dirty"

        result = remove_overlaps(x, y, pair_map)
        assert result.report_after.num_overlaps == 0
        assert result.moves_applied > 0

    def test_config_max_moves_bounds_the_search(self) -> None:
        x, y, pair_map = _nudged_multiloop(0.1)
        result = remove_overlaps(x, y, pair_map, PostPassConfig(max_moves=0))
        assert result.moves_applied == 0
        assert result.report_after.num_overlaps == result.report_before.num_overlaps

    def test_crowded_hairpin_is_fixed_by_inflation_alone(self) -> None:
        """No branch exists on either side of any witness here (see
        `_crowded_hairpin_loop`), so a fix proves the inflation move works
        on its own, independent of the rotate/translate candidates.
        """
        x, y, pair_map = _crowded_hairpin_loop(0.3)
        before = check_overlaps(x, y, pair_map, POSTPASS_PARAMS).num_overlaps
        assert before > 0, "fixture expected to be dirty before the post-pass"

        result = remove_overlaps(x, y, pair_map)
        assert result.report_after.num_overlaps == 0
        assert result.moves_applied > 0


class TestInflateLoop:
    """`inflate_loop` is the local-crowding move: it must scale a loop's own
    unpaired members and rigidly translate each child branch, without ever
    distorting a stem or touching anything outside the loop's scope.
    """

    def _multiloop_fixture(self) -> tuple[list[float], list[float], list[int]]:
        pair_map = get_pairmap_from_secstruct(TWO_HAIRPIN_MULTILOOP)
        x, y = LegacyEngine().layout(TWO_HAIRPIN_MULTILOOP)
        return x, y, pair_map

    def test_child_branches_stay_internally_rigid(self) -> None:
        x, y, pair_map = self._multiloop_fixture()
        tree = build_structure_tree(pair_map)
        multiloop = tree.loop_by_closing_pair[(0, 13)]
        nx, ny = inflate_loop(multiloop, x, y, pair_map, 1.5)
        for branch in multiloop.children:
            indices = range(branch.start, branch.end + 1)
            for i in indices:
                for j in indices:
                    before = math.hypot(x[i] - x[j], y[i] - y[j])
                    after = math.hypot(nx[i] - nx[j], ny[i] - ny[j])
                    assert after == pytest.approx(before, abs=1e-9)

    def test_leaves_own_closing_pair_untouched(self) -> None:
        # The multiloop's own closing pair (0, 13) belongs to its *parent*
        # stem, not this loop -- inflating this loop must not move it.
        x, y, pair_map = self._multiloop_fixture()
        tree = build_structure_tree(pair_map)
        multiloop = tree.loop_by_closing_pair[(0, 13)]
        nx, ny = inflate_loop(multiloop, x, y, pair_map, 1.5)
        for k in (0, 13):
            assert nx[k] == x[k]
            assert ny[k] == y[k]

    def test_unpaired_members_scale_radially_from_loop_center(self) -> None:
        x, y, pair_map = self._multiloop_fixture()
        tree = build_structure_tree(pair_map)
        multiloop = tree.loop_by_closing_pair[(0, 13)]
        cx, cy = loop_center(multiloop.members, x, y)
        nx, ny = inflate_loop(multiloop, x, y, pair_map, 1.5)
        for m in (11, 12):  # the multiloop's only bare, unpaired members
            assert nx[m] == pytest.approx(cx + 1.5 * (x[m] - cx))
            assert ny[m] == pytest.approx(cy + 1.5 * (y[m] - cy))

    def test_hairpin_loop_with_no_children_only_spreads_its_members(self) -> None:
        x, y, pair_map = self._multiloop_fixture()
        tree = build_structure_tree(pair_map)
        hairpin_a = tree.loop_by_closing_pair[(1, 4)]
        assert hairpin_a.children == []
        nx, ny = inflate_loop(hairpin_a, x, y, pair_map, 1.5)
        # everything outside this hairpin's own range is untouched
        for k in [0, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
            assert nx[k] == x[k]
            assert ny[k] == y[k]
        # its unpaired members (2, 3) did move
        assert (nx[2], ny[2]) != (x[2], y[2])
        assert (nx[3], ny[3]) != (x[3], y[3])
