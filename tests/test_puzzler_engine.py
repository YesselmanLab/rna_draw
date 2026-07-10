"""Tests for `PuzzlerEngine`: EPS coordinate parsing (no binary needed),
graceful degradation when `RNAplot` is absent, and (skip-guarded) a real
ViennaRNA integration check that it clears the M2 checker where Legacy
does not.
"""

from __future__ import annotations

import math
import random
import shutil
import subprocess

import pytest

from rna_draw.layout.base import EngineError, EngineUnavailableError
from rna_draw.layout.legacy import LegacyEngine
from rna_draw.layout.puzzler import PuzzlerEngine, _parse_coor_block
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct

_STATIC_EPS = """
%!PS-Adobe
/coor [
[10.0 20.0]
[30.0 40.5]
[-5.25 0.0]
] def
/pairs [] def
"""


class TestParseCoorBlock:
    def test_parses_points_in_order(self) -> None:
        x, y = _parse_coor_block(_STATIC_EPS)
        assert x == [10.0, 30.0, -5.25]
        assert y == [20.0, 40.5, 0.0]

    def test_missing_coor_block_raises_engine_error(self) -> None:
        with pytest.raises(EngineError):
            _parse_coor_block("%!PS-Adobe\nno coordinates here\n")


class TestPuzzlerEngineUnavailablePath:
    def test_missing_binary_raises_engine_unavailable(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        def _raise_missing(*args: object, **kwargs: object) -> None:
            raise FileNotFoundError("RNAplot not found")

        monkeypatch.setattr("rna_draw.layout.puzzler.subprocess.run", _raise_missing)
        monkeypatch.setattr(shutil, "which", lambda _binary: None)
        with pytest.raises(EngineUnavailableError):
            PuzzlerEngine().layout("((..))")

    def test_pseudoknot_raises_engine_unavailable(self) -> None:
        with pytest.raises(EngineUnavailableError):
            PuzzlerEngine().layout("([)]")

    def test_empty_structure_returns_empty_lists(self) -> None:
        assert PuzzlerEngine().layout("") == ([], [])

    def test_single_nucleotide_returns_origin(self) -> None:
        assert PuzzlerEngine().layout(".") == ([0.0], [0.0])

    def test_is_available_reflects_which(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr(shutil, "which", lambda _binary: None)
        assert PuzzlerEngine.is_available() is False


HAVE_VIENNARNA = shutil.which("RNAfold") is not None and shutil.which("RNAplot") is not None
requires_viennarna = pytest.mark.skipif(
    not HAVE_VIENNARNA, reason="RNAfold/RNAplot (ViennaRNA) not found on PATH"
)


def _random_seq(n: int, seed: int) -> str:
    """A reproducible random RNA sequence of length `n`."""
    rng = random.Random(seed)
    return "".join(rng.choice("ACGU") for _ in range(n))


def _fold(seq: str) -> str:
    """Fold `seq` with `RNAfold --noPS`; return its dot-bracket structure."""
    result = subprocess.run(
        ["RNAfold", "--noPS"], input=f"{seq}\n", capture_output=True, text=True, check=True
    )
    return result.stdout.strip().split("\n")[1].split()[0]


@requires_viennarna
class TestPuzzlerEngineIntegration:
    # Same length/seeds as `tests/test_overlap_oracle.py`: hand-picked for a
    # comfortable >10x naive/puzzler overlap-count margin at this length.
    STRUCTURE_LEN = 800
    SEEDS = (3, 7, 9)

    def test_layout_length_matches_structure(self) -> None:
        secstruct = _fold(_random_seq(self.STRUCTURE_LEN, self.SEEDS[0]))
        x, y = PuzzlerEngine().layout(secstruct)
        assert len(x) == len(secstruct)
        assert len(y) == len(secstruct)

    def test_backbone_step_matches_primary_space(self) -> None:
        secstruct = _fold(_random_seq(self.STRUCTURE_LEN, self.SEEDS[0]))
        x, y = PuzzlerEngine().layout(secstruct)
        steps = [math.hypot(x[i + 1] - x[i], y[i + 1] - y[i]) for i in range(len(x) - 1)]
        steps.sort()
        median_step = steps[len(steps) // 2]
        target = DrawParameters().PRIMARY_SPACE
        assert math.isclose(median_step, target, rel_tol=0.01)

    def test_pair_distance_is_a_sane_finite_multiple_of_backbone_step(self) -> None:
        secstruct = _fold(_random_seq(self.STRUCTURE_LEN, self.SEEDS[0]))
        pair_map = get_pairmap_from_secstruct(secstruct)
        x, y = PuzzlerEngine().layout(secstruct)
        target = DrawParameters().PRIMARY_SPACE
        distances = [
            math.hypot(x[i] - x[j], y[i] - y[j])
            for i, j in enumerate(pair_map)
            if j != -1 and i < j
        ]
        assert distances
        for distance in distances:
            assert math.isfinite(distance)
            assert 0 < distance < 20 * target

    @pytest.mark.parametrize("seed", SEEDS)
    def test_puzzler_much_cleaner_than_legacy(self, seed: int) -> None:
        secstruct = _fold(_random_seq(self.STRUCTURE_LEN, seed))
        pair_map = get_pairmap_from_secstruct(secstruct)
        params = OverlapParams(node_r=DrawParameters().NODE_R)

        legacy_x, legacy_y = LegacyEngine().layout(secstruct)
        legacy_report = check_overlaps(legacy_x, legacy_y, pair_map, params)

        puzzler_x, puzzler_y = PuzzlerEngine().layout(secstruct)
        puzzler_report = check_overlaps(puzzler_x, puzzler_y, pair_map, params)

        assert puzzler_report.num_overlaps * 10 < legacy_report.num_overlaps

    # A folded 100-nt structure (`RNAfold` on `_random_seq(100, seed=2)`) hand-picked
    # because Puzzler lays it out fully checker-clean at the default `node_r` while
    # Legacy's tree-recursion layout does not.
    CURATED_PUZZLER_CLEAN_LEGACY_DIRTY = (
        ".(((((((....(((..........))))))))))(((((((.((((.....))))..(((((....."
        "..)))))(((....)))......)))))))..\n"
    ).strip()

    def test_curated_case_puzzler_passes_while_legacy_does_not(self) -> None:
        secstruct = self.CURATED_PUZZLER_CLEAN_LEGACY_DIRTY
        pair_map = get_pairmap_from_secstruct(secstruct)
        params = OverlapParams(node_r=DrawParameters().NODE_R)

        legacy_x, legacy_y = LegacyEngine().layout(secstruct)
        legacy_report = check_overlaps(legacy_x, legacy_y, pair_map, params)

        puzzler_x, puzzler_y = PuzzlerEngine().layout(secstruct)
        puzzler_report = check_overlaps(puzzler_x, puzzler_y, pair_map, params)

        assert puzzler_report.passed
        assert not legacy_report.passed
