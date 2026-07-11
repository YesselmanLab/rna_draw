"""Pure loop/stem structure tree for a pseudoknot-free secondary structure.

This is a topological view only -- no coordinates are baked in. It answers
two questions the post-pass (`rna_draw.layout.postpass`) needs: "which
contiguous slice of nucleotide indices does this branch (stem + everything
it encloses) own?" and "which loop, and which of its child branches, do two
nucleotide indices belong to?" (their lowest common ancestor).

Because the input is assumed pseudoknot-free and well-nested (see
`rna_draw.layout.base.is_pseudoknot_free`; the gate already filters
pseudoknots before this module runs), every base pair `(i, j)` owns the
contiguous index range `[i, j]`, so a rigid transform on that range never
splits a stem or a loop.

This module deliberately does NOT reuse `render_rna.add_nodes_recursive`
(it calls `sys.exit(0)` on malformed input, which is unusable inside a
library function) -- `build_structure_tree` is its own pure recursive
descent over a pair map.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Branch:
    """A stem, and everything nested inside it, hanging off a loop.

    Args:
        closing_pair: The `(i, j)` base pair that opens this branch.
        start: First nucleotide index owned by this branch (equals `i`).
        end: Last nucleotide index owned by this branch (equals `j`).
    """

    closing_pair: tuple[int, int]
    start: int
    end: int


@dataclass
class Loop:
    """A loop: the flat ring of nucleotides directly inside a closing pair.

    Args:
        closing_pair: The base pair that closes this loop, or `None` for
            the exterior loop (the 5'/3' ends and any top-level branches).
        members: Nucleotide indices lying directly on this loop's
            boundary, in structure order: the closing pair's own two
            indices (first and last, if any), every unpaired nucleotide
            directly in the loop, and the first index of each child
            branch (the branch's own closing partner is not duplicated
            here -- it belongs to the branch's rigid slice, not the
            loop's boundary).
        children: Branches hanging directly off this loop.
    """

    closing_pair: tuple[int, int] | None
    members: list[int] = field(default_factory=list)
    children: list[Branch] = field(default_factory=list)


@dataclass
class StructureTree:
    """The full loop/branch decomposition of one pseudoknot-free structure.

    Args:
        exterior: The top-level loop (nucleotides not inside any pair).
        loops: Every loop, exterior included, in the order they were
            closed by `build_structure_tree`'s descent.
        enclosing_pairs: `enclosing_pairs[k]` is the stack of base pairs
            whose range contains index `k`, ordered outermost-first. If
            `k` is itself a pair endpoint, that pair is the innermost
            (last) entry.
        loop_by_closing_pair: Maps a closing pair (or `None` for the
            exterior loop) to its `Loop`.
    """

    exterior: Loop
    loops: list[Loop]
    enclosing_pairs: dict[int, list[tuple[int, int]]]
    loop_by_closing_pair: dict[tuple[int, int] | None, Loop]


@dataclass
class _BuildState:
    """Accumulators threaded through `_descend`'s recursion.

    Grouped into one object so `_descend` stays within the standard's
    4-parameter guideline instead of passing three accumulators loose.

    Args:
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        loops: Every completed `Loop` is appended here.
        loop_by_closing_pair: Closing-pair-to-`Loop` lookup being built.
        enclosing_pairs: Per-index nesting stacks being built.
    """

    pair_map: Sequence[int]
    loops: list[Loop]
    loop_by_closing_pair: dict[tuple[int, int] | None, Loop]
    enclosing_pairs: dict[int, list[tuple[int, int]]]


def _descend(
    state: _BuildState,
    span: tuple[int, int],
    closing_pair: tuple[int, int] | None,
    stack: list[tuple[int, int]],
) -> Loop:
    """Build one loop by scanning `span`, recursing into each branch.

    Args:
        state: Shared accumulators (see `_BuildState`).
        span: `(lo, hi)`, this loop's flat scan range (inclusive).
        closing_pair: This loop's own closing pair, or `None` for the
            exterior loop.
        stack: The enclosing-pair stack in effect for `span` (does NOT
            include `closing_pair` itself; callers pass the parent's
            stack, and a child branch's stack gets `closing_pair`
            appended when it recurses).

    Returns:
        The completed `Loop` for `span`.
    """
    lo, hi = span
    loop = Loop(closing_pair=closing_pair)
    if closing_pair is not None:
        loop.members.append(closing_pair[0])

    i = lo
    while i <= hi:
        j = state.pair_map[i]
        if j == -1:
            state.enclosing_pairs[i] = list(stack)
            loop.members.append(i)
            i += 1
            continue
        loop.members.append(i)
        loop.children.append(Branch(closing_pair=(i, j), start=i, end=j))
        child_stack = [*stack, (i, j)]
        state.enclosing_pairs[i] = child_stack
        state.enclosing_pairs[j] = child_stack
        _descend(state, (i + 1, j - 1), (i, j), child_stack)
        i = j + 1

    if closing_pair is not None:
        loop.members.append(closing_pair[1])
    state.loops.append(loop)
    state.loop_by_closing_pair[closing_pair] = loop
    return loop


def build_structure_tree(pair_map: Sequence[int]) -> StructureTree:
    """Build the loop/branch tree of a pseudoknot-free structure.

    Args:
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired (see `render_rna.get_pairmap_from_secstruct`).
            Must be well-nested (pseudoknot-free); see module docstring.

    Returns:
        The `StructureTree` describing every loop and branch.
    """
    state = _BuildState(pair_map=pair_map, loops=[], loop_by_closing_pair={}, enclosing_pairs={})
    exterior = _descend(state, (0, len(pair_map) - 1), None, [])
    return StructureTree(
        exterior=exterior,
        loops=state.loops,
        enclosing_pairs=state.enclosing_pairs,
        loop_by_closing_pair=state.loop_by_closing_pair,
    )


def loop_center(
    members: Sequence[int], x: Sequence[float], y: Sequence[float]
) -> tuple[float, float]:
    """Centroid of a set of nucleotides at their current coordinates.

    Computed on demand rather than stored on `Loop`, because coordinates
    mutate as the post-pass applies rigid moves.

    Args:
        members: Nucleotide indices to average (typically `Loop.members`,
            but any index sequence works -- e.g. a branch's own range).
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.

    Returns:
        `(cx, cy)`, the mean position of `members`.

    Raises:
        ValueError: If `members` is empty.
    """
    if not members:
        raise ValueError("cannot compute a center of zero members")
    cx = sum(x[m] for m in members) / len(members)
    cy = sum(y[m] for m in members) / len(members)
    return cx, cy


def _child_branch(loop: Loop, stack: Sequence[tuple[int, int]], depth: int) -> Branch | None:
    """Look up the child branch a nucleotide's stack diverges into.

    Args:
        loop: The candidate LCA loop.
        stack: The nucleotide's full enclosing-pair stack.
        depth: Index into `stack` where it diverges from the other side's
            stack (equivalently, `len` of the common prefix).

    Returns:
        The `Branch` of `loop` whose closing pair is `stack[depth]`, or
        `None` if `stack` ends exactly at `depth` (the nucleotide is a
        bare member of `loop`, not inside any of its child branches).
    """
    if depth >= len(stack):
        return None
    next_pair = stack[depth]
    for branch in loop.children:
        if branch.closing_pair == next_pair:
            return branch
    return None


def lca_loop_and_branches(
    tree: StructureTree, ka: int, kb: int
) -> tuple[Loop, Branch | None, Branch | None]:
    """Find the lowest common ancestor loop of two nucleotide indices.

    Walks `ka` and `kb`'s `enclosing_pairs` stacks (outermost first) to
    their deepest common pair -- the LCA loop -- then, for each side,
    identifies the direct child branch of that loop the index diverges
    into (or `None` if the index is itself a bare member of the LCA loop).

    Args:
        tree: The structure tree built by `build_structure_tree`.
        ka: First nucleotide index.
        kb: Second nucleotide index.

    Returns:
        `(loop, branch_a, branch_b)`: the LCA `Loop`, and each side's
        child `Branch` (or `None` when that side is a direct member of
        `loop` rather than nested inside one of its branches).
    """
    stack_a = tree.enclosing_pairs.get(ka, [])
    stack_b = tree.enclosing_pairs.get(kb, [])
    depth = 0
    while depth < len(stack_a) and depth < len(stack_b) and stack_a[depth] == stack_b[depth]:
        depth += 1
    lca_pair = stack_a[depth - 1] if depth > 0 else None
    loop = tree.loop_by_closing_pair[lca_pair]
    return loop, _child_branch(loop, stack_a, depth), _child_branch(loop, stack_b, depth)


__all__ = [
    "Branch",
    "Loop",
    "StructureTree",
    "build_structure_tree",
    "loop_center",
    "lca_loop_and_branches",
]
