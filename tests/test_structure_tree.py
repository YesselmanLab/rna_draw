"""Tests for `rna_draw.layout.structure_tree`: tree builder + LCA correctness.

Sample structures follow the plan's fixtures:
  - `"((..))"`: a simple hairpin.
  - `"(((...)))"`: a stacked-helix hairpin.
  - `"((..)((..))..)"`: a two-hairpin multiloop nested under one outer pair.
  - `"...."`: a pure exterior loop (no pairs at all).
"""

from __future__ import annotations

import pytest

from rna_draw.layout.structure_tree import (
    Branch,
    Loop,
    StructureTree,
    build_structure_tree,
    lca_loop_and_branches,
    loop_center,
)
from rna_draw.render_rna import get_pairmap_from_secstruct

SIMPLE_HAIRPIN = "((..))"
STACKED_HAIRPIN = "(((...)))"
TWO_HAIRPIN_MULTILOOP = "((..)((..))..)"
PURE_EXTERIOR = "...."


def _tree(secstruct: str) -> StructureTree:
    """Build a `StructureTree` for a dot-bracket structure under test."""
    return build_structure_tree(get_pairmap_from_secstruct(secstruct))


class TestBuildStructureTree:
    """`build_structure_tree` must never raise/exit and must produce a tree
    whose branch ranges are disjoint, contiguous, and match their closing
    pair -- the load-bearing structural invariant rigid moves rely on.
    """

    def test_simple_hairpin_exterior_branch(self) -> None:
        tree = _tree(SIMPLE_HAIRPIN)
        assert [b.closing_pair for b in tree.exterior.children] == [(0, 5)]
        branch = tree.exterior.children[0]
        assert branch.start == 0
        assert branch.end == 5
        assert branch.closing_pair == (0, 5)

    def test_simple_hairpin_loop_members(self) -> None:
        tree = _tree(SIMPLE_HAIRPIN)
        hairpin_loop = tree.loop_by_closing_pair[(1, 4)]
        assert hairpin_loop.members == [1, 2, 3, 4]
        assert hairpin_loop.children == []

    def test_stacked_hairpin_branch_ranges_contiguous(self) -> None:
        tree = _tree(STACKED_HAIRPIN)
        # (0, 8) wraps (1, 7) wraps (2, 6): each inner loop has exactly the
        # next pair as its lone branch, and ranges strictly nest.
        outer = tree.loop_by_closing_pair[(0, 8)].children
        assert [b.closing_pair for b in outer] == [(1, 7)]
        middle = tree.loop_by_closing_pair[(1, 7)].children
        assert [b.closing_pair for b in middle] == [(2, 6)]
        assert tree.loop_by_closing_pair[(2, 6)].members == [2, 3, 4, 5, 6]

    def test_branch_ranges_disjoint_and_contiguous(self) -> None:
        tree = _tree(TWO_HAIRPIN_MULTILOOP)
        multiloop = tree.loop_by_closing_pair[(0, 13)]
        branches = sorted(multiloop.children, key=lambda b: b.start)
        assert [(b.start, b.end) for b in branches] == [(1, 4), (5, 10)]
        # disjoint: no overlap between consecutive branch ranges
        assert branches[0].end < branches[1].start

    def test_multiloop_members_include_closing_pair_and_branch_starts(self) -> None:
        tree = _tree(TWO_HAIRPIN_MULTILOOP)
        multiloop = tree.loop_by_closing_pair[(0, 13)]
        assert multiloop.members == [0, 1, 5, 11, 12, 13]

    def test_pure_exterior_has_no_branches(self) -> None:
        tree = _tree(PURE_EXTERIOR)
        assert tree.exterior.children == []
        assert tree.exterior.members == [0, 1, 2, 3]
        assert tree.loops == [tree.exterior]

    def test_does_not_raise_or_exit_on_valid_input(self) -> None:
        # A regression guard for the reuse-map warning: this module must be
        # its own pure recursion, never `render_rna.add_nodes_recursive`
        # (which calls `sys.exit(0)` on malformed input).
        tree = _tree(TWO_HAIRPIN_MULTILOOP)
        assert isinstance(tree.exterior, Loop)


class TestEnclosingPairs:
    """`enclosing_pairs[k]` must be the outer->inner nesting stack for `k`."""

    def test_hairpin_loop_member_stack(self) -> None:
        tree = _tree(TWO_HAIRPIN_MULTILOOP)
        assert tree.enclosing_pairs[2] == [(0, 13), (1, 4)]

    def test_deeply_nested_member_stack(self) -> None:
        tree = _tree(TWO_HAIRPIN_MULTILOOP)
        assert tree.enclosing_pairs[7] == [(0, 13), (5, 10), (6, 9)]

    def test_exterior_member_has_no_enclosing_pair(self) -> None:
        tree = _tree(PURE_EXTERIOR)
        assert tree.enclosing_pairs[0] == []

    def test_pair_endpoints_share_the_same_stack(self) -> None:
        tree = _tree(TWO_HAIRPIN_MULTILOOP)
        assert tree.enclosing_pairs[1] == tree.enclosing_pairs[4]


class TestLoopCenter:
    def test_centroid_of_members(self) -> None:
        x = [0.0, 10.0, 10.0, 0.0]
        y = [0.0, 0.0, 10.0, 10.0]
        assert loop_center([0, 1, 2, 3], x, y) == (5.0, 5.0)

    def test_empty_members_raises(self) -> None:
        with pytest.raises(ValueError, match="zero members"):
            loop_center([], [0.0], [0.0])


class TestLcaLoopAndBranches:
    """LCA correctness is what lets the post-pass find the right branches
    to move for a given pair of clashing nucleotide indices.
    """

    def test_cross_hairpin_witness_resolves_to_multiloop_and_both_branches(self) -> None:
        tree = _tree(TWO_HAIRPIN_MULTILOOP)
        # nt 2 is inside hairpin A's loop (1, 4); nt 7 is inside the
        # innermost hairpin (6, 9), itself nested inside branch (5, 10).
        loop, branch_a, branch_b = lca_loop_and_branches(tree, 2, 7)
        assert loop.closing_pair == (0, 13)
        assert branch_a == Branch(closing_pair=(1, 4), start=1, end=4)
        assert branch_b == Branch(closing_pair=(5, 10), start=5, end=10)
        assert branch_a != branch_b

    def test_top_level_exterior_clash_has_no_branches(self) -> None:
        tree = _tree(PURE_EXTERIOR)
        loop, branch_a, branch_b = lca_loop_and_branches(tree, 0, 3)
        assert loop.closing_pair is None
        assert branch_a is None
        assert branch_b is None

    def test_same_loop_members_have_no_branch(self) -> None:
        tree = _tree(TWO_HAIRPIN_MULTILOOP)
        loop, branch_a, branch_b = lca_loop_and_branches(tree, 2, 3)
        assert loop.closing_pair == (1, 4)
        assert branch_a is None
        assert branch_b is None

    def test_branch_start_index_resolves_into_its_own_branch(self) -> None:
        tree = _tree(TWO_HAIRPIN_MULTILOOP)
        # nt 1 is branch (1, 4)'s own opening index; against a nucleotide
        # in the sibling branch it should resolve as "inside branch (1, 4)".
        loop, branch_a, branch_b = lca_loop_and_branches(tree, 1, 7)
        assert loop.closing_pair == (0, 13)
        assert branch_a == Branch(closing_pair=(1, 4), start=1, end=4)
        assert branch_b == Branch(closing_pair=(5, 10), start=5, end=10)
