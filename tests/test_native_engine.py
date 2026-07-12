"""Tests for `NativeTurtleEngine` (`rna_draw.layout.native`): the
`LayoutEngine` seam over the native `rna_layout` turtle-base port.

Coordinate parity against the vendored oracle lives in
`test_native_parity.py`; this file covers the engine contract itself
(rescale, malformed-input/pseudoknot guards, `LayoutEngine` Protocol
conformance) -- the same shape `test_vienna_binding.py` uses for
`ViennaTurtleEngine`.
"""

from __future__ import annotations

import math

import pytest

from rna_draw.layout.base import EngineError, EngineUnavailableError, LayoutEngine, has_empty_loop
from rna_draw.layout.native import NativeTurtleEngine
from rna_draw.overlap import check_overlaps
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct


class TestLayoutEngineProtocolConformance:
    def test_native_turtle_engine_satisfies_protocol(self) -> None:
        engine: LayoutEngine = NativeTurtleEngine()
        assert isinstance(engine, LayoutEngine)
        assert engine.name == "native_turtle"


class TestNativeTurtleEngineLayout:
    def test_empty_structure_returns_empty_lists(self) -> None:
        x, y = NativeTurtleEngine().layout("")
        assert x == []
        assert y == []

    def test_single_nucleotide_returns_origin(self) -> None:
        x, y = NativeTurtleEngine().layout(".")
        assert x == [0.0]
        assert y == [0.0]

    def test_output_length_matches_structure(self) -> None:
        structure = "((((....))))"
        x, y = NativeTurtleEngine().layout(structure)
        assert len(x) == len(structure)
        assert len(y) == len(structure)

    def test_rescales_to_primary_space(self) -> None:
        params = DrawParameters()
        structure = "((((....))))((((....))))"
        x, y = NativeTurtleEngine(params).layout(structure)
        steps = [
            math.hypot(x[i + 1] - x[i], y[i + 1] - y[i]) for i in range(len(x) - 1)
        ]
        steps.sort()
        median_step = steps[len(steps) // 2]
        assert math.isclose(median_step, params.PRIMARY_SPACE, rel_tol=0.15)

    def test_pseudoknot_raises_unavailable(self) -> None:
        with pytest.raises(EngineUnavailableError):
            NativeTurtleEngine().layout("(((...[[[...)))...]]]")

    def test_unbalanced_structure_raises_engine_error(self) -> None:
        with pytest.raises(EngineError):
            NativeTurtleEngine().layout("((")

    def test_empty_loop_raises_engine_error(self) -> None:
        with pytest.raises(EngineError):
            NativeTurtleEngine().layout("().()")

    @pytest.mark.parametrize("repeats", [1, 3, 7, 20, 50])
    def test_repeated_hairpin_structures_produce_finite_coords(self, repeats: int) -> None:
        # `random_structure` (`conftest.py`) reliably produces a bare "()"
        # empty loop for any n large enough to be interesting here (both
        # engines reject that shape uniformly -- see
        # `test_empty_loop_raises_engine_error` above, and
        # `test_vienna_binding.py`'s `_no_empty_loop_structure` for the
        # same known generator quirk), so this property test instead varies
        # STRUCTURE SIZE via repeated hairpin units for real, non-skipped
        # coverage across sizes.
        structure = "((((....))))" * repeats
        assert not has_empty_loop(structure)
        x, y = NativeTurtleEngine().layout(structure)
        assert len(x) == len(structure)
        assert all(math.isfinite(v) for v in x)
        assert all(math.isfinite(v) for v in y)

    def test_deterministic(self) -> None:
        structure = "((((..((((....))))..))))"
        first = NativeTurtleEngine().layout(structure)
        second = NativeTurtleEngine().layout(structure)
        assert first == second


class TestNativeTurtleEngineOverlapReport:
    """Sanity check the engine slots into `check_overlaps` honestly (it is
    not registered in the production selection chain this slice -- see
    `rna_draw.layout.native`'s module docstring -- but must still produce
    a checkable layout when called directly).
    """

    def test_hairpin_is_checkable(self) -> None:
        structure = "((((....))))"
        x, y = NativeTurtleEngine().layout(structure)
        pair_map = get_pairmap_from_secstruct(structure)
        report = check_overlaps(x, y, pair_map)
        assert report is not None
