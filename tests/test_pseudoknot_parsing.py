"""Tests for `rna_draw.layout.pseudoknot.parsing`: multi-bracket parsing,
stem grouping, and the nested `()`-only reconstruction round-trip.
"""

from __future__ import annotations

import pytest

from rna_draw.layout.pseudoknot.parsing import (
    Stem,
    group_stems,
    nested_secstruct,
    parse_all_pairs,
    stem_pairs,
)
from rna_draw.render_rna import get_pairmap_from_secstruct


class TestParseAllPairs:
    def test_h_type_example_returns_six_pairs(self) -> None:
        # "((([[[)))]]]" from the plan's Phase 1 example.
        pairs = parse_all_pairs("((([[[)))]]]")
        assert pairs == [(0, 8), (1, 7), (2, 6), (3, 11), (4, 10), (5, 9)]

    def test_nested_only_structure(self) -> None:
        assert parse_all_pairs("((()))") == [(0, 5), (1, 4), (2, 3)]

    def test_all_four_bracket_types(self) -> None:
        pairs = parse_all_pairs("(<>)")
        assert (0, 3) in pairs
        assert (1, 2) in pairs

    def test_unpaired_only(self) -> None:
        assert parse_all_pairs("....") == []

    def test_unbalanced_bracket_left_unpaired(self) -> None:
        assert parse_all_pairs("(((.") == []

    def test_rejects_unsupported_character(self) -> None:
        with pytest.raises(ValueError):
            parse_all_pairs("(((xyz")


class TestGroupStems:
    def test_h_type_example_returns_two_stems(self) -> None:
        pairs = parse_all_pairs("((([[[)))]]]")
        stems = group_stems(pairs)
        assert stems == [Stem(i=0, j=8, length=3), Stem(i=3, j=11, length=3)]

    def test_nested_only_returns_one_stem_no_crossings(self) -> None:
        stems = group_stems(parse_all_pairs("((()))"))
        assert stems == [Stem(i=0, j=5, length=3)]

    def test_empty_pairs_returns_no_stems(self) -> None:
        assert group_stems([]) == []

    def test_two_disjoint_stacked_stems(self) -> None:
        # "(())" is ONE 2-rung stacked stem (0,3) then (1,2); two such
        # motifs separated by unpaired nts stay two SEPARATE stems.
        pairs = parse_all_pairs("(())..(())")
        stems = group_stems(pairs)
        assert stems == [Stem(i=0, j=3, length=2), Stem(i=6, j=9, length=2)]


class TestStemPairsAndDataclass:
    def test_stem_pairs_outermost_first(self) -> None:
        stem = Stem(i=0, j=8, length=3)
        assert stem.pairs() == [(0, 8), (1, 7), (2, 6)]

    def test_stem_pairs_flattens_multiple_stems(self) -> None:
        stems = [Stem(i=0, j=2, length=1), Stem(i=4, j=6, length=1)]
        assert stem_pairs(stems) == [(0, 2), (4, 6)]

    def test_length_one_stem(self) -> None:
        assert Stem(i=5, j=9, length=1).pairs() == [(5, 9)]


class TestNestedSecstructRoundTrip:
    @pytest.mark.parametrize(
        "secstruct",
        [
            "((()))",
            "((((....))))",
            "..((..))..",
            "....",
            "([{...}]).",  # nested through THREE different bracket types
        ],
    )
    def test_retained_subset_round_trips(self, secstruct: str) -> None:
        # Round-trip through the FULL multi-bracket parse -> group -> emit
        # pipeline: nested_secstruct(n, stem_pairs(all stems)) must
        # reproduce the exact retained pair_map when the input has no
        # crossings (every stem is independently retained). A genuine
        # crossing input is NOT expected to round-trip this way -- that is
        # exactly why extraction (tests/test_pseudoknot_extraction.py)
        # separates crossing stems out before calling `nested_secstruct`.
        stems = group_stems(parse_all_pairs(secstruct))
        pairs = stem_pairs(stems)
        n = len(secstruct)
        nested = nested_secstruct(n, pairs)

        expected = [-1] * n
        for i, j in pairs:
            expected[i], expected[j] = j, i
        assert get_pairmap_from_secstruct(nested) == expected

    def test_nested_secstruct_is_dots_and_parens_only(self) -> None:
        nested = nested_secstruct(6, [(0, 5), (1, 4)])
        assert set(nested) <= set("().")
        assert nested == "((..))"

    def test_nested_secstruct_all_unpaired(self) -> None:
        assert nested_secstruct(4, []) == "...."
