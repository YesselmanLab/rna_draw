"""Phase 1b: split parsed stems into a maximum-weight nested subset (the
retained tree) and the crossing stems (the pseudoknot itself), via Maximum
Weight Independent Set on the stem-crossing conflict graph.

Exact by branch-and-bound over only the stems actually involved in a
crossing -- the corpus's crossing-involved-stem count is tiny (plan
`Corpus facts`: minimum-stems-to-remove-to-reach-nested, measured on 316
real pseudoknots, median 3, p90 4, max 7), so a full polynomial circle-graph
algorithm is unnecessary at this scale. Falls back to a greedy max-degree
removal only for the never-observed pathological case above
`EXACT_COVER_LIMIT`. Cite: Gavril, "Algorithms for a maximum clique and a
maximum independent set of a circle graph," Networks 3 (1973), for the
polynomial-time guarantee this shortcuts around.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .parsing import Stem

# Branch-and-bound is exponential worst case; this caps how many
# crossing-involved stems it is ever attempted on (corpus max measured: 7,
# see module docstring) before falling back to the greedy heuristic.
EXACT_COVER_LIMIT = 22


def stems_cross(a: Stem, b: Stem) -> bool:
    """Whether two stems' outer pairs interleave (a genuine crossing).

    Args:
        a: First stem.
        b: Second stem.

    Returns:
        True iff `a`'s and `b`'s outer-pair index intervals properly
        interleave (`a.i < b.i < a.j < b.j`, or the symmetric case) rather
        than nest or sit disjoint -- the standard chord-crossing test. This
        is sufficient at the whole-stem level because every stem's own
        indices are disjoint from every other stem's (each index appears
        in at most one base pair), so a stem's two arms and everything it
        encloses share one of exactly these three interval relationships
        with another stem's.
    """
    return (a.i < b.i < a.j < b.j) or (b.i < a.i < b.j < a.j)


def max_nested_subset(stems: list[Stem]) -> tuple[list[Stem], list[Stem]]:
    """Split `stems` into a maximum-weight nested subset + the crossing rest.

    Args:
        stems: Every stem parsed from a (possibly pseudoknotted) structure.

    Returns:
        `(retained, crossing)`: `retained` stems are pairwise non-crossing
        (their reconstructed dot-bracket string is pseudoknot-free) and
        have maximum total base-pair count among such subsets; `crossing`
        is every other stem.
    """
    conflict = _conflict_graph(stems)
    involved = [idx for idx, neighbors in enumerate(conflict) if neighbors]
    free = [idx for idx, neighbors in enumerate(conflict) if not neighbors]

    if len(involved) <= EXACT_COVER_LIMIT:
        kept = _mwis_exact(stems, conflict, involved)
    else:
        kept = _mwis_greedy(stems, conflict, involved)

    retained_idx = set(free) | kept
    retained = [stems[idx] for idx in sorted(retained_idx)]
    crossing = [stems[idx] for idx in range(len(stems)) if idx not in retained_idx]
    return retained, crossing


def _conflict_graph(stems: list[Stem]) -> list[set[int]]:
    """Adjacency list: `conflict[k]` is every stem index crossing stem `k`."""
    conflict: list[set[int]] = [set() for _ in stems]
    for a in range(len(stems)):
        for b in range(a + 1, len(stems)):
            if stems_cross(stems[a], stems[b]):
                conflict[a].add(b)
                conflict[b].add(a)
    return conflict


@dataclass
class _BranchState:
    """Shared accumulators for `_branch`'s recursive search.

    Args:
        weight: Stem index -> base-pair count (the MWIS objective).
        conflict: Stem index -> set of crossing stem indices.
        best_value: Highest total weight found so far.
        best_chosen: The stem-index set achieving `best_value`.
    """

    weight: dict[int, int]
    conflict: list[set[int]]
    best_value: int = 0
    best_chosen: set[int] = field(default_factory=set)


def _mwis_exact(stems: list[Stem], conflict: list[set[int]], involved: list[int]) -> set[int]:
    """Exact Maximum Weight Independent Set via branch-and-bound.

    Args:
        stems: Every stem (indexed as `conflict` and `involved` refer to).
        conflict: Crossing adjacency, see `_conflict_graph`.
        involved: Stem indices with at least one crossing (only these are
            searched; crossing-free stems are always retained by the
            caller, `max_nested_subset`).

    Returns:
        The subset of `involved` achieving maximum total `Stem.length`,
        with no two chosen stems crossing.
    """
    weight = {idx: stems[idx].length for idx in involved}
    state = _BranchState(weight=weight, conflict=conflict)
    ordered = sorted(involved, key=lambda idx: -weight[idx])
    _branch(ordered, set(), 0, state)
    return state.best_chosen


def _branch(remaining: list[int], chosen: set[int], value: int, state: _BranchState) -> None:
    """Recursive include/exclude search with a weight-sum upper-bound prune."""
    if value > state.best_value:
        state.best_value = value
        state.best_chosen = set(chosen)
    if not remaining:
        return
    if value + sum(state.weight[idx] for idx in remaining) <= state.best_value:
        return

    node, *rest = remaining
    _branch(rest, chosen, value, state)

    non_conflicting = [idx for idx in rest if idx not in state.conflict[node]]
    chosen.add(node)
    _branch(non_conflicting, chosen, value + state.weight[node], state)
    chosen.discard(node)


def _mwis_greedy(stems: list[Stem], conflict: list[set[int]], involved: list[int]) -> set[int]:
    """Greedy fallback: repeatedly drop the max-degree (tie: shortest) stem.

    Args:
        stems: Every stem.
        conflict: Crossing adjacency, see `_conflict_graph`.
        involved: Stem indices with at least one crossing.

    Returns:
        A (not necessarily maximum) independent subset of `involved`.
    """
    remaining = set(involved)
    while True:
        degrees = {idx: len(conflict[idx] & remaining) for idx in remaining}
        crossing_now = {idx: deg for idx, deg in degrees.items() if deg > 0}
        if not crossing_now:
            return remaining
        drop = max(crossing_now, key=lambda idx: (crossing_now[idx], -stems[idx].length))
        remaining.discard(drop)


__all__ = ["stems_cross", "max_nested_subset", "EXACT_COVER_LIMIT"]
