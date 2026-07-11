"""Tests for `rna_draw.layout.production`: the escalation + monotone,
time-bounded post-pass composition wired into the real render pipeline.

Mirrors `tests/test_vienna_binding.py`'s convention of importing the
compiled `rna_draw._vienna_layout` extension unconditionally at module
level -- this module needs it to construct any engine here.
"""

from __future__ import annotations

import pytest
from conftest import random_structure

import rna_draw.layout.production as production_module
from rna_draw.layout.base import EngineError
from rna_draw.layout.postpass import POSTPASS_PARAMS
from rna_draw.layout.production import (
    PRODUCTION_CLEARANCE_LADDER,
    PRODUCTION_POSTPASS_BUDGET_S,
    EscalatingClearanceEngine,
    PostPassEngine,
    production_engine,
)
from rna_draw.overlap import check_overlaps
from rna_draw.render_rna import get_pairmap_from_secstruct

STRUCTURES = [
    "((((....))))",
    "((((....))))((((....))))((((....))))",
    "(((..((..((...))..))..)))",
]


class TestEscalatingClearanceEngineFinite:
    @pytest.mark.parametrize("structure", STRUCTURES)
    def test_finite_coords_matching_length(self, structure: str) -> None:
        x, y = EscalatingClearanceEngine().layout(structure)
        assert len(x) == len(structure)
        assert len(y) == len(structure)
        assert all(v == v and abs(v) != float("inf") for v in x)
        assert all(v == v and abs(v) != float("inf") for v in y)

    def test_empty_structure(self) -> None:
        x, y = EscalatingClearanceEngine().layout("")
        assert x == []
        assert y == []

    def test_single_nucleotide(self) -> None:
        x, y = EscalatingClearanceEngine().layout(".")
        assert x == [0.0]
        assert y == [0.0]


class TestEscalatingClearanceEngineGuards:
    @pytest.mark.timeout(20)
    @pytest.mark.parametrize("structure", ["()", "(().())", "().()"])
    def test_empty_loop_raises_engine_error_no_hang(self, structure: str) -> None:
        # The load-bearing hang-safety guard (see production.py's module
        # docstring): puzzler hangs unconditionally on a bare "()" loop, at
        # any clearance -- this must raise fast, never hang the caller.
        with pytest.raises(EngineError):
            EscalatingClearanceEngine().layout(structure)

    def test_pseudoknot_raises_engine_error(self) -> None:
        with pytest.raises(EngineError):
            EscalatingClearanceEngine().layout("([)]")

    def test_c_failure_becomes_engine_error(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # `pipeline._try_primary` catches ONLY `EngineError`; a raw
        # RuntimeError/ValueError escaping `layout()` would crash the
        # caller instead of routing to the safe fallback (Blocking Fix 2).
        def _boom(*args: object, **kwargs: object) -> tuple[list[float], list[float]]:
            raise RuntimeError("simulated C-level failure")

        monkeypatch.setattr(
            production_module._vienna_layout, "plot_coords_puzzler_opts", _boom
        )
        with pytest.raises(EngineError):
            EscalatingClearanceEngine().layout("((((....))))")


class TestEscalatingClearanceEngineABIGuard:
    def test_constructor_asserts_vienna_abi(self, monkeypatch: pytest.MonkeyPatch) -> None:
        def _drifted() -> tuple[int, int, int]:
            return (0, 0, 0)

        monkeypatch.setattr(production_module._vienna_layout, "abi_version", _drifted)
        with pytest.raises(Exception):
            EscalatingClearanceEngine()


class TestPostPassEngine:
    def test_default_name_is_postpass(self) -> None:
        assert PostPassEngine().name == "postpass"

    def test_custom_name_kwarg(self) -> None:
        assert PostPassEngine(name="custom").name == "custom"

    @pytest.mark.parametrize("structure", STRUCTURES)
    def test_never_worse_than_base(self, structure: str) -> None:
        base = EscalatingClearanceEngine()
        bx, by = base.layout(structure)
        pair_map = get_pairmap_from_secstruct(structure)
        before = check_overlaps(bx, by, pair_map, POSTPASS_PARAMS).num_overlaps

        px, py = PostPassEngine(base=base).layout(structure)
        after = check_overlaps(px, py, pair_map, POSTPASS_PARAMS).num_overlaps
        assert after <= before


class TestProductionEngineFactory:
    def test_name_is_production(self) -> None:
        assert production_engine().name == "production"

    def test_uses_capped_clearance_ladder(self) -> None:
        assert PRODUCTION_CLEARANCE_LADDER == (1.0, 1.25, 1.5)

    def test_uses_finite_postpass_budget(self) -> None:
        assert PRODUCTION_POSTPASS_BUDGET_S == 2.0

    @pytest.mark.timeout(20)
    @pytest.mark.parametrize("structure", STRUCTURES)
    def test_finite_coords(self, structure: str) -> None:
        engine = production_engine()
        x, y = engine.layout(structure)
        assert len(x) == len(structure)
        assert len(y) == len(structure)

    @pytest.mark.timeout(20)
    @pytest.mark.parametrize("seed", range(5))
    @pytest.mark.parametrize("n", [5, 20, 60])
    def test_no_hang_on_random_structures(self, seed: int, n: int) -> None:
        # `random_structure` can independently generate empty-loop or
        # pseudoknot-free-but-degenerate shapes; the guard(s) must route
        # these to a fast `EngineError`, never a hang.
        secstruct = random_structure(seed, n)
        try:
            production_engine().layout(secstruct)
        except EngineError:
            pass
