"""Tests for the `constructive` layout engine's M2 scope (see the
constructive-engine runbook, milestone M2): general bottom-up composition
over the WHOLE structure tree -- bulges/interior loops (2a), 2-level
nesting (2b), deeper/mixed nesting (2c), and the exterior loop's open
boundary with tails and multiple top-level branches (2d) -- plus the S5
envelope containment invariant and the S8 never-silent-overlap contract.

Every test that runs the engine is `@pytest.mark.timeout`-marked (the hard
requirement: the pure-Python engine must never hang the suite).
"""

from __future__ import annotations

import math
import time

import pytest
from conftest import random_structure

from rna_draw.layout.base import EngineError
from rna_draw.layout.constructive import envelope
from rna_draw.layout.constructive.engine import _MAX_NUCLEOTIDES, ConstructiveEngine
from rna_draw.layout.structure_tree import build_structure_tree, collapse_stem
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.render_rna import get_pairmap_from_secstruct

TIMEOUT = 20
PARAMS_OVERLAP = OverlapParams()


def _assert_clean(secstruct: str) -> None:
    """Lay `secstruct` out and assert the frozen checker finds 0 witnesses."""
    x, y = ConstructiveEngine().layout(secstruct)
    pair_map = get_pairmap_from_secstruct(secstruct)
    report = check_overlaps(x, y, pair_map, PARAMS_OVERLAP)
    assert report.passed, f"{secstruct!r} left {report.num_overlaps} overlaps: {report.witnesses}"


class TestBulgesAndInteriorLoops:
    """2a: a stem with unpaired bases on one/both strands before a hairpin."""

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize(
        "secstruct",
        [
            "((.(...)))",  # bulge on the 5' strand only
            "((...).)",  # bulge on the 3' strand only
            "((.(...).))",  # symmetric interior loop
            "((..(...)))",  # 2-nt bulge, 5' side
            "((.(...)..))",  # asymmetric interior loop
            "(((((...))).))",  # deep collapsed run then a 3' bulge
            "(((((...)).)))",  # deep collapsed run then a 5' bulge
        ],
    )
    def test_bulge_or_interior_loop_is_clean(self, secstruct: str) -> None:
        _assert_clean(secstruct)

    @pytest.mark.timeout(TIMEOUT)
    def test_bulge_does_not_collapse_the_stem(self) -> None:
        """A bulge loop must NOT be treated as a stack continuation: the
        (0, 9)/(1, 8) run collapses (no unpaired members between them), but
        collapsing must STOP at (1, 8) -- the loop with the bulge -- rather
        than continuing through it to (3, 7).
        """
        secstruct = "((.(...)))"
        pair_map = get_pairmap_from_secstruct(secstruct)
        tree = build_structure_tree(pair_map)
        depth, loop = collapse_stem(tree, (0, 9))
        assert depth == 2
        assert loop.closing_pair == (1, 8)
        assert len(loop.children) == 1


def make_bulge_chain(levels: int, n_before: int = 1, n_after: int = 0, loop_size: int = 3) -> str:
    """A chain of `levels` nested single-child bulge/interior loops.

    Each level has `n_before` unpaired nts before its child stem opens and
    `n_after` after it closes, ending in a `loop_size`-nt hairpin.

    Args:
        levels: Number of nested bulge/interior loops in the chain.
        n_before: Unpaired count on the near (5') strand at every level.
        n_after: Unpaired count on the far (3') strand at every level.
        loop_size: Unpaired nucleotide count in the terminal hairpin loop.

    Returns:
        A dot-bracket structure.
    """
    open_part = ("." * n_before + "(") * levels
    close_part = (")" + "." * n_after) * levels
    return "(" + open_part + "." * loop_size + close_part + ")"


class TestBulgeChainStraightPlacement:
    """The bulge-chain envelope-compounding fix: a single-child bulge or
    interior loop is placed as a straight continuation of the parent
    stem's axis (`engine._place_bulge`) instead of a fresh circular
    envelope, so a CHAIN of such loops has LINEAR (not exponential) reach.
    """

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize(
        ("n_before", "n_after"),
        [(1, 0), (0, 1), (2, 2), (1, 4), (4, 1), (3, 0), (0, 3)],
    )
    def test_asymmetric_bulge_chain_is_clean(self, n_before: int, n_after: int) -> None:
        _assert_clean(make_bulge_chain(15, n_before, n_after))

    @pytest.mark.timeout(TIMEOUT)
    def test_bulge_chain_feeding_a_multiloop_is_clean(self) -> None:
        """A bulge chain that terminates in a multiloop (not a hairpin)."""
        chain = make_bulge_chain(10, n_before=2, n_after=1, loop_size=0)
        secstruct = chain[:-1] + "(...)(...)" + chain[-1]
        _assert_clean(secstruct)

    @pytest.mark.timeout(TIMEOUT)
    def test_reach_grows_linearly_not_exponentially(self) -> None:
        """`envelope.branch_reach` on a bulge chain of length `2N` must stay
        within a small constant factor of length `N`'s reach -- the OLD
        circular-envelope compounding grew ~3x PER LEVEL (so doubling chain
        length would blow the ratio up astronomically); the FIXED
        straight-continuation reach only adds a constant per level, so the
        ratio stays close to 2 (bounded generously at 3 to allow for the
        formula's other additive terms).
        """
        from rna_draw.parameters import DrawParameters

        params = DrawParameters()
        short_secstruct = make_bulge_chain(10)
        long_secstruct = make_bulge_chain(20)

        short_tree = build_structure_tree(get_pairmap_from_secstruct(short_secstruct))
        long_tree = build_structure_tree(get_pairmap_from_secstruct(long_secstruct))
        short_branch = short_tree.exterior.children[0]
        long_branch = long_tree.exterior.children[0]

        short_reach = envelope.branch_reach(
            short_tree, short_branch.closing_pair, params, envelope.ReachCache()
        )
        long_reach = envelope.branch_reach(
            long_tree, long_branch.closing_pair, params, envelope.ReachCache()
        )
        assert long_reach / short_reach < 3.0, (
            f"reach ratio {long_reach / short_reach:.2f} for doubled chain length "
            "-- suggests exponential compounding is back"
        )

    @pytest.mark.timeout(30)
    def test_deep_tail_length_lays_out_without_hanging(self) -> None:
        """A large (~3000nt) random well-nested structure -- representative
        of the deep-tail structures `_MAX_NUCLEOTIDES` now admits (raised
        1200 -> 4000 after this fix; see `engine.py`'s docstring) -- must
        complete well inside the benchmark harness's 30s per-structure kill,
        clean-or-flagged (S8's never-silent-overlap contract), never hang.
        """
        secstruct = random_structure(0, 3000)
        start = time.monotonic()
        try:
            x, y = ConstructiveEngine().layout(secstruct)
        except EngineError:
            return
        elapsed = time.monotonic() - start
        assert elapsed < 15.0, f"took {elapsed:.2f}s -- too slow for the deep tail"
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map, PARAMS_OVERLAP).passed


def make_degree2_chain(levels: int, side_loop: int = 3, main_loop: int = 3) -> str:
    """A chain of `levels` nested degree-2 multiloops (one small side
    hairpin + one continuing branch each), ending in a `main_loop`-nt
    terminal hairpin -- the pattern `envelope._degree2_packing`'s collinear
    straight-continuation targets: common in real rRNA, and what used to
    compound `envelope.branch_reach` ~3x PER NESTING LEVEL before this fix
    (see `engine.py`'s `_MAX_REACH` docstring).

    Args:
        levels: Number of nested degree-2 loops in the chain.
        side_loop: Unpaired nucleotide count in each level's side hairpin.
        main_loop: Unpaired nucleotide count in the terminal hairpin loop.

    Returns:
        A dot-bracket structure.
    """
    inner = "(" + "." * main_loop + ")"
    for _ in range(levels):
        side = "(" + "." * side_loop + ")"
        inner = "(" + side + inner + ")"
    return inner


class TestDegree2ChainStraightPlacement:
    """The degree-2-chain envelope-compounding fix
    (`envelope.lateral_reach`/`_degree2_packing`): once a run of nested
    degree-2 multiloops is long enough to matter
    (`envelope._DEGREE2_CHAIN_LENGTH_FLOOR`), the dominant child at every
    level is pinned as a collinear straight continuation, sized by its
    tighter, directional `lateral_reach` instead of its isotropic
    `branch_reach` -- collapsing the circular envelope's `~3x`-per-level
    compounding.
    """

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize("levels", [1, 3, 6, 10, 20, 30])
    def test_degree2_chain_is_clean(self, levels: int) -> None:
        _assert_clean(make_degree2_chain(levels))

    @pytest.mark.timeout(TIMEOUT)
    def test_lateral_reach_grows_linearly_not_exponentially(self) -> None:
        """`envelope.lateral_reach` (NOT `branch_reach` -- see
        `test_branch_reach_grows_boundedly_not_exponentially` below for why
        that one is a different story) on a degree-2 chain of length `2N`
        must stay within a small constant factor of length `N`'s reach:
        `lateral_reach`'s own recursion (`radius + branch_reach(side) +
        margin`, `envelope._degree2_lateral_reach`) adds a roughly CONSTANT
        increment per level once every level is pinned, so doubling the
        chain length should only roughly double it, not compound
        exponentially.
        """
        from rna_draw.parameters import DrawParameters

        params = DrawParameters()
        short_secstruct = make_degree2_chain(20)
        long_secstruct = make_degree2_chain(40)

        short_tree = build_structure_tree(get_pairmap_from_secstruct(short_secstruct))
        long_tree = build_structure_tree(get_pairmap_from_secstruct(long_secstruct))
        short_branch = short_tree.exterior.children[0]
        long_branch = long_tree.exterior.children[0]

        short_reach = envelope.lateral_reach(
            short_tree, short_branch.closing_pair, params, envelope.ReachCache()
        )
        long_reach = envelope.lateral_reach(
            long_tree, long_branch.closing_pair, params, envelope.ReachCache()
        )
        assert long_reach / short_reach < 3.0, (
            f"lateral_reach ratio {long_reach / short_reach:.2f} for doubled chain "
            "length -- suggests exponential compounding is back"
        )

    @pytest.mark.timeout(TIMEOUT)
    def test_branch_reach_grows_boundedly_not_exponentially(self) -> None:
        """`envelope.branch_reach` itself stays POLYNOMIAL (empirically
        quadratic), not exponential, down a degree-2 chain -- a dramatic
        improvement over the pre-fix `~3^N`, though not as tight as
        `lateral_reach`'s own linear growth (above). This is a genuine,
        understood limit of the fix, not a bug: `_two_seam_packing`'s own
        radius floor must include `dom`'s `lateral_reach` (needed for the
        packer's chord-based disjointness argument to stay sound -- see
        that function's docstring), so each level's own packing radius
        tracks the (linearly growing) `lateral_reach` below it; `branch_
        reach`'s `2 * radius + child_reach` formula (`_loop_branch_reach`,
        deliberately UNCHANGED) then accumulates that linearly-growing
        radius once per level, giving `O(N^2)` overall -- still comfortably
        bounded (see `test_deep_chain_stays_far_under_the_reach_cap`), just
        not asymptotically linear.
        """
        from rna_draw.parameters import DrawParameters

        params = DrawParameters()
        short_secstruct = make_degree2_chain(20)
        long_secstruct = make_degree2_chain(40)

        short_tree = build_structure_tree(get_pairmap_from_secstruct(short_secstruct))
        long_tree = build_structure_tree(get_pairmap_from_secstruct(long_secstruct))
        short_branch = short_tree.exterior.children[0]
        long_branch = long_tree.exterior.children[0]

        short_reach = envelope.branch_reach(
            short_tree, short_branch.closing_pair, params, envelope.ReachCache()
        )
        long_reach = envelope.branch_reach(
            long_tree, long_branch.closing_pair, params, envelope.ReachCache()
        )
        assert long_reach / short_reach < 6.0, (
            f"branch_reach ratio {long_reach / short_reach:.2f} for doubled chain "
            "length -- suggests exponential (not polynomial) compounding is back"
        )

    @pytest.mark.timeout(TIMEOUT)
    def test_deep_chain_stays_far_under_the_reach_cap(self) -> None:
        """A 60-level degree-2 chain (well past `_DEGREE2_CHAIN_LENGTH_
        FLOOR`, and deeper than any real hard-set structure measured for
        this fix) must still land far under `_MAX_REACH`, and its rendered
        bounding box must stay of a similarly bounded (not astronomical)
        size -- the pre-fix `~3^60` would be meaningless, unrepresentable
        magnitude; the post-fix quadratic growth keeps it small.
        """
        secstruct = make_degree2_chain(60)
        x, y = ConstructiveEngine().layout(secstruct)
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map, PARAMS_OVERLAP).passed
        width, height = max(x) - min(x), max(y) - min(y)
        assert width < 1e6 and height < 1e6, f"bbox {width} x {height} -- unexpectedly large"


class TestTwoLevelNesting:
    """2b: a multiloop whose child is itself a multiloop of hairpins."""

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize(
        "secstruct",
        [
            "(((...)(...))(...))",
            "((((...)(...)(...)))(...))",
            "((...)((...)(...)))",
            "(((...)(...))((...)(...)))",
            "((((...)(...))(...))(...))",
        ],
    )
    def test_nested_multiloop_is_clean(self, secstruct: str) -> None:
        _assert_clean(secstruct)


class TestDeeperMixedNesting:
    """2c: multiloops containing internal loops containing multiloops, etc."""

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize(
        "secstruct",
        [
            "((.((...)(...)).)(...))",
            "(((.(...)(...))(..(...)))(...))",
            "((((.(...).)(...))(...)))",
            "(.((...)(.(...)(...)))..)",
            "((((...)(...))((...)(...)))((...)(...)))",
        ],
    )
    def test_deep_mixed_nesting_is_clean(self, secstruct: str) -> None:
        _assert_clean(secstruct)

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize("seed", range(30))
    def test_random_well_nested_structures_are_clean_or_flagged(self, seed: int) -> None:
        """Random structures never come back silently dirty (S8 contract)."""
        secstruct = random_structure(seed, 40)
        try:
            x, y = ConstructiveEngine().layout(secstruct)
        except EngineError:
            return
        pair_map = get_pairmap_from_secstruct(secstruct)
        report = check_overlaps(x, y, pair_map, PARAMS_OVERLAP)
        assert report.passed


class TestExteriorLoop:
    """2d: the one OPEN-boundary loop -- tails and multiple top-level branches."""

    @pytest.mark.timeout(TIMEOUT)
    def test_five_and_three_prime_tails(self) -> None:
        _assert_clean("...((...))...")

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize("degree", range(2, 6))
    def test_multiple_top_level_stems(self, degree: int) -> None:
        _assert_clean("((...))" * degree)

    @pytest.mark.timeout(TIMEOUT)
    def test_all_unpaired_string(self) -> None:
        _assert_clean("." * 25)

    @pytest.mark.timeout(TIMEOUT)
    def test_single_hairpin_with_tails(self) -> None:
        _assert_clean(".....((...)).....")

    @pytest.mark.timeout(TIMEOUT)
    def test_bare_empty_loop_does_not_crash(self) -> None:
        """A degenerate `()` (no interior at all) is legal input; the
        engine must either lay it out clean or raise `EngineError` -- never
        crash or hang.
        """
        try:
            x, y = ConstructiveEngine().layout("().()")
        except EngineError:
            return
        pair_map = get_pairmap_from_secstruct("().()")
        assert check_overlaps(x, y, pair_map, PARAMS_OVERLAP).passed


class TestEnvelopeContainmentInvariant:
    """S5: every primitive in a branch's subtree lies within `branch_reach`
    of its attachment point -- checked directly against real coordinates,
    not just inferred from a clean checker run.
    """

    @pytest.mark.timeout(TIMEOUT)
    @pytest.mark.parametrize(
        "secstruct",
        [
            "((...)(...))",
            "(((...)(...))(...))",
            "((.(...).)(...)(...))",
            "((((...)(...))(...))((...)(...)))",
        ],
    )
    def test_subtree_stays_within_its_reach(self, secstruct: str) -> None:
        from rna_draw.parameters import DrawParameters

        pair_map = get_pairmap_from_secstruct(secstruct)
        tree = build_structure_tree(pair_map)
        params = DrawParameters()
        cache = envelope.ReachCache()
        x, y = ConstructiveEngine().layout(secstruct)

        for branch in tree.exterior.children:
            reach = envelope.branch_reach(tree, branch.closing_pair, params, cache)
            attachment_index = branch.closing_pair[0]
            ax, ay = x[attachment_index], y[attachment_index]
            for idx in range(branch.start, branch.end + 1):
                dist = math.hypot(x[idx] - ax, y[idx] - ay)
                assert dist <= reach + 1e-6, (
                    f"{secstruct!r} index {idx} at distance {dist} exceeds "
                    f"reach {reach} from attachment {attachment_index}"
                )


class TestNeverSilentOverlap:
    """S8: the engine either returns checker-clean coordinates or raises
    `EngineError` -- it must NEVER return a checker-dirty layout silently.
    """

    @pytest.mark.timeout(60)
    @pytest.mark.parametrize("seed", range(60))
    def test_fuzz_never_returns_silent_overlap(self, seed: int) -> None:
        secstruct = random_structure(seed, 50)
        try:
            x, y = ConstructiveEngine().layout(secstruct)
        except EngineError:
            return
        pair_map = get_pairmap_from_secstruct(secstruct)
        assert check_overlaps(x, y, pair_map, PARAMS_OVERLAP).passed

    @pytest.mark.timeout(TIMEOUT)
    def test_pseudoknot_raises(self) -> None:
        with pytest.raises(EngineError):
            ConstructiveEngine().layout("([)]")

    @pytest.mark.timeout(TIMEOUT)
    def test_over_max_nucleotides_raises_immediately(self) -> None:
        secstruct = "." * (_MAX_NUCLEOTIDES + 1)
        with pytest.raises(EngineError, match="declines structures over"):
            ConstructiveEngine().layout(secstruct)

    @pytest.mark.timeout(TIMEOUT)
    def test_a_genuinely_dirty_case_raises_not_returns(self) -> None:
        """A known-dirty structure (few, size-disparate loop slots -- the
        residual documented in the M2 diagnosis) raises `EngineError`
        rather than returning silently-dirty coordinates.
        """
        secstruct = ".(..(..)(..))(..)."
        with pytest.raises(EngineError, match="dirty layout"):
            ConstructiveEngine().layout(secstruct)

    @pytest.mark.timeout(TIMEOUT)
    def test_long_bulge_chain_is_clean_and_fast(self) -> None:
        """A long run of single-child bulge/interior loops -- common in real
        rRNA -- used to compound `envelope.branch_reach` ~3x PER LEVEL
        (circular envelope on every loop), reaching astronomical magnitudes
        well under `_MAX_NUCLEOTIDES` in length and forcing the `_MAX_REACH`
        guard to bail out (see the git history around this test for the
        prior "raises fast, not slow" contract). The straight-continuation
        bulge placement (`engine._place_bulge`) fixes the root cause: reach
        grows LINEARLY with chain length, so this now lays out fast AND
        checker-clean, no guard needed.
        """
        secstruct = "(" + ".(" * 40 + "." * 3 + ")." * 40 + ")"
        assert len(secstruct) < 200
        start = time.monotonic()
        x, y = ConstructiveEngine().layout(secstruct)
        elapsed = time.monotonic() - start
        pair_map = get_pairmap_from_secstruct(secstruct)
        report = check_overlaps(x, y, pair_map, PARAMS_OVERLAP)
        assert report.passed, f"left {report.num_overlaps} overlaps: {report.witnesses}"
        assert elapsed < 1.0, f"took {elapsed:.3f}s -- should be near-instant now"
