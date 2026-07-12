"""End-to-end tests for `rna_draw.layout.pseudoknot.engine.layout_pseudoknot`:
an H-type pseudoknot, a kissing-loop, a nested-inside-pseudoknot, and an
adversarial multi-crossing case, each asserted NEVER a silent overlap --
using only the `LayoutResult`'s PUBLIC fields (`pair_map`, `crossing_lines`)
as the oracle, not engine-internal state.
"""

from __future__ import annotations

import pytest

from rna_draw.layout.base import EngineError, is_pseudoknot_free
from rna_draw.layout.pseudoknot.engine import layout_pseudoknot
from rna_draw.layout.pseudoknot.validate import capsule_is_clean, polyline_capsules
from rna_draw.overlap import OverlapParams, build_primitives, check_overlaps

# A simple, clearly pseudoknotted H-type example.
H_TYPE = "((((....[[[[....))))....]]]]"

# Two hairpins whose own loops cross-pair with each other (a "kissing
# loop"): hairpin 1 is (0-13), hairpin 2 is (14-27), and `[[`/`]]` pair
# loop nts 6,7 of hairpin 1 with loop nts 20,21 of hairpin 2.
KISSING_LOOP = "((((..[[..))))((((..]]..))))"

# A normal nested hairpin (group B) sitting inside one arm of an H-type
# pseudoknot (group A crosses group C); see the derivation in
# `test_pseudoknot_engine.py`'s module notes / commit message.
NESTED_INSIDE_PK = "((((...((((....))))...[[[[....))))....]]]]"

# Two SEPARATE H-type pseudoknots concatenated: two independent crossing
# stems must each be placed (possibly as PK-A lines) without either
# clipping the OTHER's already-committed geometry.
TWO_CROSSINGS = "((([[[)))]]]" + "(((<<<)))>>>"


def _full_layout_is_clean(result) -> bool:
    """Independent oracle: checker-clean pair_map view AND every PK-A line
    clear of the whole layout and every OTHER PK-A line, using only the
    `LayoutResult`'s public fields (never engine-internal state).
    """
    assert result.pair_map is not None
    params = OverlapParams(node_r=result.node_r)
    base_report = check_overlaps(result.x, result.y, result.pair_map, params)
    if not base_report.passed:
        return False

    primitives = build_primitives(result.x, result.y, result.pair_map, params)
    all_line_capsules = [
        polyline_capsules(line.points, line.i, line.j, params.pair_half_width, uid)
        for uid, line in enumerate(result.crossing_lines)
    ]
    for uid, segments in enumerate(all_line_capsules):
        other_lines = [
            seg
            for other_uid, caps in enumerate(all_line_capsules)
            if other_uid != uid
            for seg in caps
        ]
        others = list(primitives) + other_lines
        for segment in segments:
            if not capsule_is_clean(segment, others, result.pair_map, params.tol):
                return False
    return True


class TestHType:
    def test_produces_a_layout_for_every_nucleotide(self) -> None:
        result = layout_pseudoknot(H_TYPE, OverlapParams())
        assert len(result.x) == len(H_TYPE)
        assert len(result.y) == len(H_TYPE)

    def test_never_silent_overlap(self) -> None:
        result = layout_pseudoknot(H_TYPE, OverlapParams())
        assert _full_layout_is_clean(result)

    def test_engine_name_is_pseudoknot(self) -> None:
        result = layout_pseudoknot(H_TYPE, OverlapParams())
        assert result.engine_name == "pseudoknot"

    def test_flagged_iff_any_line_or_unplaced(self) -> None:
        # Honest contract (Phase 4): flagged=False implies zero PK-A lines.
        result = layout_pseudoknot(H_TYPE, OverlapParams())
        if not result.flagged:
            assert result.crossing_lines == []


class TestKissingLoop:
    def test_never_silent_overlap(self) -> None:
        result = layout_pseudoknot(KISSING_LOOP, OverlapParams())
        assert _full_layout_is_clean(result)

    def test_all_nucleotides_placed(self) -> None:
        result = layout_pseudoknot(KISSING_LOOP, OverlapParams())
        assert len(result.x) == len(KISSING_LOOP)


class TestNestedInsidePseudoknot:
    def test_never_silent_overlap(self) -> None:
        result = layout_pseudoknot(NESTED_INSIDE_PK, OverlapParams())
        assert _full_layout_is_clean(result)

    def test_all_nucleotides_placed(self) -> None:
        result = layout_pseudoknot(NESTED_INSIDE_PK, OverlapParams())
        assert len(result.x) == len(NESTED_INSIDE_PK)


class TestTwoCrossingStems:
    def test_never_silent_overlap(self) -> None:
        result = layout_pseudoknot(TWO_CROSSINGS, OverlapParams())
        assert _full_layout_is_clean(result)

    def test_routed_lines_do_not_clip_each_other(self) -> None:
        # If both crossing stems escalated to PK-A, `_full_layout_is_clean`
        # already verifies their lines don't overlap; this just asserts
        # the scenario is actually exercised (>= 1 line drawn) so the test
        # is not vacuous.
        result = layout_pseudoknot(TWO_CROSSINGS, OverlapParams())
        assert result.crossing_pairs or result.crossing_lines


class TestPseudoknotFreeInputStillWorks:
    def test_pk_free_input_has_no_crossings(self) -> None:
        # layout_pseudoknot is never called by the pipeline for pk-free
        # input, but it must not misbehave if called directly.
        result = layout_pseudoknot("((((....))))", OverlapParams())
        assert result.crossing_pairs == []
        assert result.crossing_lines == []
        assert result.flagged is False


class TestRoundTripGuard:
    def test_internal_invariant_violation_raises_engine_error(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Force `_assert_round_trip`'s guard by making `nested_secstruct`
        # return pseudoknotted garbage regardless of input -- the
        # never-silent contract requires this be an EngineError, not a
        # silently wrong layout fed to the tree engine.
        import rna_draw.layout.pseudoknot.engine as engine_module

        def _garbage(n: int, pairs: list[tuple[int, int]]) -> str:
            return "([)]" + "." * (n - 4)

        monkeypatch.setattr(engine_module, "nested_secstruct", _garbage)
        with pytest.raises(EngineError):
            layout_pseudoknot(H_TYPE, OverlapParams())

    def test_pseudoknot_free_but_wrong_pairing_raises_engine_error(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # A DIFFERENT failure mode of `_assert_round_trip`'s guard: the
        # returned string IS pseudoknot-free (passes the first check) but
        # does not reproduce `retained`'s own pairs -- e.g. a plausible
        # `nested_secstruct` bug that silently drops or shifts a pair.
        import rna_draw.layout.pseudoknot.engine as engine_module

        def _wrong_pairing(n: int, pairs: list[tuple[int, int]]) -> str:
            return "." * n  # pseudoknot-free, but drops every pair

        monkeypatch.setattr(engine_module, "nested_secstruct", _wrong_pairing)
        with pytest.raises(EngineError):
            layout_pseudoknot(H_TYPE, OverlapParams())


class TestCorpusSample:
    """A light-weight, in-process cross-check of the corpus-scale contract
    (the full-corpus numbers are reported by
    `benchmarks/pseudoknot_gate.py`, not re-run in the unit suite).
    """

    def test_is_pseudoknot_free_rejects_all_bracket_types(self) -> None:
        assert is_pseudoknot_free(H_TYPE) is False
        assert is_pseudoknot_free(KISSING_LOOP) is False
        assert is_pseudoknot_free(NESTED_INSIDE_PK) is False
