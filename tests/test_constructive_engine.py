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
    def test_long_bulge_chain_raises_fast_not_slow(self) -> None:
        """A long run of single-child bulge/interior loops -- common in real
        rRNA -- compounds `envelope.branch_reach` ~3x per level; well under
        `_MAX_NUCLEOTIDES` in length, this must still raise `EngineError`
        PROMPTLY (the `_MAX_REACH` guard) rather than spend minutes building
        astronomical coordinates the checker then grinds through (a real
        near-hang found while measuring on `benchmarks/hard_set.json`).
        """
        secstruct = "(" + ".(" * 40 + "." * 3 + ")." * 40 + ")"
        assert len(secstruct) < 200
        with pytest.raises(EngineError, match="envelope reach"):
            ConstructiveEngine().layout(secstruct)
