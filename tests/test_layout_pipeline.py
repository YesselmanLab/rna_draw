"""Tests for `rna_draw.layout.base` helpers and the `layout_guaranteed`
pipeline: gating, fallback, engine selection, pseudoknots, and the
adaptive render-radius search.
"""

from __future__ import annotations

from collections.abc import Callable

import pytest
from conftest import random_structure

from rna_draw.layout.base import empty_report, has_empty_loop, is_pseudoknot_free
from rna_draw.layout.legacy import LegacyEngine
from rna_draw.layout.pipeline import (
    default_engine,
    layout_guaranteed,
    resolve_engine,
)
from rna_draw.layout.puzzler import PuzzlerEngine
from rna_draw.overlap import OverlapParams


class FakeEngine:
    """A `LayoutEngine` implementation returning caller-supplied coords.

    The mock seam for pipeline tests -- no subprocess, no tree recursion.
    """

    name = "fake"

    def __init__(
        self,
        coords_fn: Callable[[str], tuple[list[float], list[float]]] | None = None,
    ) -> None:
        """Store the coordinate-generating callable.

        Args:
            coords_fn: Maps a secstruct to `(x, y)`; defaults to placing
                every nucleotide at the origin (guaranteed overlapping for
                any structure with 2+ nucleotides).
        """
        self._coords_fn = coords_fn or (lambda ss: ([0.0] * len(ss), [0.0] * len(ss)))

    def layout(self, secstruct: str) -> tuple[list[float], list[float]]:
        """Return the stored coordinates for `secstruct`."""
        return self._coords_fn(secstruct)


class TestBaseHelpers:
    @pytest.mark.parametrize("secstruct", ["(())", "....", "", "(((.)))"])
    def test_is_pseudoknot_free_true_cases(self, secstruct: str) -> None:
        assert is_pseudoknot_free(secstruct) is True

    @pytest.mark.parametrize("secstruct", ["([)]", "((", "(]"])
    def test_is_pseudoknot_free_false_cases(self, secstruct: str) -> None:
        assert is_pseudoknot_free(secstruct) is False

    @pytest.mark.parametrize("secstruct", ["()", "(())"])
    def test_has_empty_loop_true_cases(self, secstruct: str) -> None:
        assert has_empty_loop(secstruct) is True

    @pytest.mark.parametrize("secstruct", ["(...)", "....", ""])
    def test_has_empty_loop_false_cases(self, secstruct: str) -> None:
        assert has_empty_loop(secstruct) is False

    def test_empty_report_passed(self) -> None:
        assert empty_report().passed is True


class TestLayoutGuaranteedLegacyClean:
    def test_legacy_engine_yields_clean_unflagged_result(self) -> None:
        result = layout_guaranteed("((((....))))", engine=LegacyEngine())
        assert result.flagged is False
        assert result.report.passed
        assert result.engine_name == "legacy"


class TestLayoutGuaranteedFallback:
    def test_bad_fake_engine_falls_back_and_stays_clean(self) -> None:
        result = layout_guaranteed("((((....))))", engine=FakeEngine())
        assert result.flagged is True
        assert result.engine_name == "fallback"
        assert result.report.passed is True

    @pytest.mark.parametrize("seed", range(15))
    @pytest.mark.parametrize("n", [2, 5, 20, 60])
    def test_never_silent_overlap_property(self, seed: int, n: int) -> None:
        secstruct = random_structure(seed, n)
        result = layout_guaranteed(secstruct, engine=FakeEngine())
        assert result.flagged or result.report.passed


class TestEngineSelection:
    def test_default_engine_is_puzzler_when_available(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(PuzzlerEngine, "is_available", staticmethod(lambda: True))
        assert isinstance(default_engine(), PuzzlerEngine)

    def test_default_engine_is_legacy_when_puzzler_unavailable(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(PuzzlerEngine, "is_available", staticmethod(lambda: False))
        assert isinstance(default_engine(), LegacyEngine)

    def test_pipeline_falls_back_to_legacy_or_fallback_when_puzzler_absent(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(PuzzlerEngine, "is_available", staticmethod(lambda: False))
        result = layout_guaranteed("((((....))))")
        assert result.engine_name in {"legacy", "fallback"}


class TestPseudoknotAndEmpty:
    def test_pseudoknot_is_flagged(self) -> None:
        result = layout_guaranteed("([)]")
        assert result.flagged is True

    def test_empty_structure(self) -> None:
        result = layout_guaranteed("")
        assert result.x == []
        assert result.y == []
        assert result.flagged is False


class TestResolveEngine:
    def test_auto_is_none(self) -> None:
        assert resolve_engine("auto") is None

    def test_legacy_returns_legacy_engine(self) -> None:
        assert isinstance(resolve_engine("legacy"), LegacyEngine)

    def test_bogus_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            resolve_engine("bogus")


class TestAdaptiveRenderRadius:
    """A layout that fails at `node_r=10` but clears a smaller radius must
    be accepted by `_largest_clean_node_r` rather than routed to the safe
    fallback.
    """

    # Three unpaired nucleotides: nt1 sits far away (perpendicular), so its
    # backbone capsules never approach the other two disks; nt0 and nt2 sit
    # 18 layout units apart -- closer than the disk-disk clearance at the
    # default `node_r=10` (required 2*10=20 > 18, overlap) but no closer
    # than the clearance at `node_r=9` (required 2*9=18 == 18, clean).
    NEAR_TOUCH_X = [0.0, 0.0, 18.0]
    NEAR_TOUCH_Y = [0.0, 200.0, 0.0]

    def test_shrinks_node_r_instead_of_falling_back(self) -> None:
        engine = FakeEngine(lambda ss: (self.NEAR_TOUCH_X, self.NEAR_TOUCH_Y))
        result = layout_guaranteed("...", engine=engine, params=OverlapParams())

        assert result.flagged is False
        assert result.engine_name == "fake"
        assert result.report.passed
        assert result.node_r < 10.0
        assert result.node_r == pytest.approx(9.0)
