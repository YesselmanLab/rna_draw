"""Tests for `rna_draw.layout.pseudoknot.extraction`: crossing detection and
the maximum-weight-nested-subset split (exact branch-and-bound + greedy
fallback), verified against real corpus pseudoknots.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from rna_draw.layout.base import is_pseudoknot_free
from rna_draw.layout.pseudoknot import extraction
from rna_draw.layout.pseudoknot.extraction import max_nested_subset, stems_cross
from rna_draw.layout.pseudoknot.parsing import (
    Stem,
    group_stems,
    nested_secstruct,
    parse_all_pairs,
    stem_pairs,
)

CORPUS = Path("/Users/jyesselman2/Downloads/dbnFiles")


def _parse_dbn(path: Path) -> tuple[str, str] | None:
    """Mirror `benchmarks.hard_gate.parse_dbn` (kept local to avoid an
    inter-test-module import): `(sequence, structure)` from a bpRNA `.dbn`.
    """
    seq = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if seq is None:
            seq = line
            continue
        if set(line) <= set(".()[]{}<>anAB"):
            return seq, line
    return None


def _real_pseudoknots(limit: int) -> list[str]:
    """Deterministically sample `limit` real pseudoknotted structures."""
    if not CORPUS.is_dir():
        return []
    files = sorted(CORPUS.glob("*.dbn"))
    rng = random.Random(3)
    rng.shuffle(files)
    found = []
    for path in files:
        parsed = _parse_dbn(path)
        if parsed is None:
            continue
        _seq, struct = parsed
        if is_pseudoknot_free(struct) or set(struct) - set(".()[]{}<>"):
            continue
        found.append(struct)
        if len(found) >= limit:
            break
    return found


requires_corpus = pytest.mark.skipif(not CORPUS.is_dir(), reason="corpus not available")


class TestStemsCross:
    def test_interleaved_stems_cross(self) -> None:
        a, b = Stem(i=0, j=8, length=1), Stem(i=3, j=11, length=1)
        assert stems_cross(a, b) is True
        assert stems_cross(b, a) is True

    def test_nested_stems_do_not_cross(self) -> None:
        outer, inner = Stem(i=0, j=10, length=1), Stem(i=2, j=5, length=1)
        assert stems_cross(outer, inner) is False

    def test_disjoint_stems_do_not_cross(self) -> None:
        a, b = Stem(i=0, j=2, length=1), Stem(i=4, j=6, length=1)
        assert stems_cross(a, b) is False


class TestMaxNestedSubset:
    def test_h_type_removes_one_of_two_equal_stems(self) -> None:
        stems = group_stems(parse_all_pairs("((([[[)))]]]"))
        retained, crossing = max_nested_subset(stems)
        assert len(retained) == 1
        assert len(crossing) == 1
        assert set(retained) | set(crossing) == set(stems)

    def test_no_crossings_retains_everything(self) -> None:
        stems = group_stems(parse_all_pairs("((()))"))
        retained, crossing = max_nested_subset(stems)
        assert retained == stems
        assert crossing == []

    def test_retained_subset_is_pseudoknot_free(self) -> None:
        stems = group_stems(parse_all_pairs("((([[[)))]]]"))
        retained, _crossing = max_nested_subset(stems)
        nested = nested_secstruct(12, stem_pairs(retained))
        assert is_pseudoknot_free(nested) is True

    def test_greedy_fallback_matches_exact_on_small_cover(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Force the greedy path (normally reserved for covers >
        # EXACT_COVER_LIMIT) on a small H-type case and confirm it finds
        # the SAME optimum the exact branch-and-bound does here (both
        # stems tie at weight 3, so removing either is optimal).
        stems = group_stems(parse_all_pairs("((([[[)))]]]"))
        monkeypatch.setattr(extraction, "EXACT_COVER_LIMIT", 0)
        retained, crossing = max_nested_subset(stems)
        assert len(retained) == 1
        assert len(crossing) == 1

    def test_three_way_crossing_keeps_max_weight_pair(self) -> None:
        # Three mutually-crossing stems of DIFFERENT weight: the MWIS
        # optimum for a triangle of crossings is exactly one stem (any
        # two would still conflict), so the heaviest single stem wins.
        heavy = Stem(i=0, j=9, length=5)
        b = Stem(i=2, j=11, length=2)
        c = Stem(i=4, j=13, length=1)
        # Make b and c cross too (both cross `heavy` already):
        # b=(2,11), c=(4,13): 2<4<11<13 -> crossing.
        retained, crossing = max_nested_subset([heavy, b, c])
        assert retained == [heavy]
        assert set(crossing) == {b, c}


class TestExtractionOnRealCorpus:
    @requires_corpus
    def test_retained_subset_always_pseudoknot_free(self) -> None:
        for struct in _real_pseudoknots(50):
            n = len(struct)
            stems = group_stems(parse_all_pairs(struct))
            retained, _crossing = max_nested_subset(stems)
            nested = nested_secstruct(n, stem_pairs(retained))
            assert is_pseudoknot_free(nested) is True, struct

    @requires_corpus
    def test_crossing_set_is_small(self) -> None:
        # Plan `Corpus facts`: minimum stems removed p90 == 4, max == 7
        # (316-structure measurement). A fresh 50-structure sample should
        # comfortably stay within a generous multiple of that ceiling.
        for struct in _real_pseudoknots(50):
            stems = group_stems(parse_all_pairs(struct))
            _retained, crossing = max_nested_subset(stems)
            assert len(crossing) <= 20, struct
