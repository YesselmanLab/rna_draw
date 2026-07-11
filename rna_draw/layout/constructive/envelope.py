"""Bottom-up subtree envelope ("reach") computation (runbook S5).

A pure topology pass: for every branch in a `StructureTree`, compute a
radius bounding its ENTIRE RHO-inflated subtree in a disk centered at the
branch's own attachment point -- the point later fed to `place_stem` as
`base_point` (see `rna_draw.layout.constructive.engine`). A parent's packer
(`pack_loop_angles` for a loop's own circle, `pack_line_positions` for the
exterior's open boundary) then sizes each child's slot from this
provably-sufficient bound instead of a fixed guess -- M1's gap, and the
reason the first M2 attempt's composition came back dirty (see the
engine's module docstring).

SOUNDNESS ARGUMENT: `pack_loop_angles` guarantees, for ANY two packed slots
i, k (adjacent or not, wrapping through the reserved seam or not), that the
straight-line (chord) distance between their anchor points is at least
`half_width_i + half_width_k`. Proof: each slot's half-angle is
`asin(half_width / radius) <= pi / 2`, so any two slots' combined
half-angle `a + b <= pi`; by concavity of `sin` on `[0, pi]`,
`sin(a) + sin(b) <= 2 * sin((a + b) / 2)`, i.e. `chord(a + b) >=
half_width_i + half_width_k` (the boundary case, when the slots are
adjacent). For NON-adjacent slots, the true angular gap `theta` (summed
through every intervening slot) satisfies `theta + a + b <= 2 * pi` --
because `theta + a + b` is exactly the sum of the FULL widths of every
slot from `i` to `k` inclusive, a subset of the sum over ALL packed slots,
which `pack_loop_angles`' own fitting condition already bounds by `2 * pi`.
Since chord only depends on `min(theta, 2*pi - theta)` (`sin(theta / 2) ==
sin(pi - theta / 2)`), `theta + a + b <= 2 * pi` is exactly what's needed
for `chord(theta) >= chord(a + b) >= half_width_i + half_width_k` to hold
even when `theta > pi`. So feeding `branch_reach(child)` as a slot's
`half_width` keeps two sibling subtrees' RHO-inflated disks disjoint,
which the runbook's geometry contract states is enough clearance for the
checker's disk/capsule primitives. `pack_line_positions` gets the same
guarantee directly (monotone positions on an open line put every pair,
not just adjacent ones, at least as far apart as their combined extents).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rna_draw.layout.structure_tree import Branch, Loop, StructureTree, collapse_stem
from rna_draw.parameters import DrawParameters

from .geometry_helpers import LoopPacking, pack_bulge_linear, pack_loop_angles

# node_r + max(backbone_half_width, pair_half_width) at the checker's
# target geometry (`OverlapParams()`'s defaults: 10 + 7.5) -- the runbook's
# Minkowski-inflation radius. Added once per branch, so nesting depth `d`
# compounds it `d` times -- deliberately generous (correctness first; M3
# owns tightening this for compactness).
RHO_MARGIN = 17.5

# Same clearance-margin constants M1 used for a loop's own reserved seam
# and a bare unpaired member's slot (see `engine.py`'s original comment):
# just enough to clear float noise and a stem rung's own disk offset from
# its anchor point.
_UNPAIRED_MARGIN = 0.5
_STEM_RUNG_MARGIN = 9.0

_RADIUS_STEP = 1.0
_RADIUS_SEARCH_SLACK = 1.0


@dataclass
class ReachCache:
    """Memoized S5 results, keyed by closing pair.

    Populated once, bottom-up, by the first call to `branch_reach` on a
    structure's top-level branches (see `engine._place_exterior`); every
    descendant loop's packing and every descendant branch's reach get
    filled in as a side effect of that single recursive walk, so the
    placement pass (`engine._place_loop_members`) never recomputes them --
    it just looks them up, guaranteeing pass 1 and pass 2 agree exactly.

    Args:
        branch_reach: `reach(branch)` for every branch closing pair seen so
            far -- a disk radius, centered at the branch's own attachment
            point, bounding its entire RHO-inflated subtree.
        loop_packing: Each interior loop's precomputed `LoopPacking` (its
            closing pair -> `(radius, angles)`).
    """

    branch_reach: dict[tuple[int, int], float] = field(default_factory=dict)
    loop_packing: dict[tuple[int, int], LoopPacking] = field(default_factory=dict)


def stem_rung_half_width(params: DrawParameters, margin_scale: float = 1.0) -> float:
    """Chord-clearance half-width a stem's outermost rung needs on a loop circle.

    Args:
        params: Target geometry (`NODE_R`, `PAIR_SPACE`).
        margin_scale: Multiplier on `_STEM_RUNG_MARGIN` -- S8's bounded
            re-inflation schedule (see `engine._build_and_place`).

    Returns:
        Half the rung's physical `PAIR_SPACE` separation, plus a disk
        radius and a clearance margin (`_STEM_RUNG_MARGIN`).
    """
    return params.PAIR_SPACE / 2.0 + params.NODE_R + _STEM_RUNG_MARGIN * margin_scale


def unpaired_half_width(params: DrawParameters, margin_scale: float = 1.0) -> float:
    """Chord-clearance half-width a bare unpaired member's own slot needs.

    Args:
        params: Target geometry (`NODE_R`).
        margin_scale: Multiplier on `_UNPAIRED_MARGIN` (see `stem_rung_half_width`).

    Returns:
        The nucleotide disk radius plus a clearance margin (`_UNPAIRED_MARGIN`).
    """
    return params.NODE_R + _UNPAIRED_MARGIN * margin_scale


def branch_reach(
    tree: StructureTree,
    closing_pair: tuple[int, int],
    params: DrawParameters,
    cache: ReachCache,
    margin_scale: float = 1.0,
) -> float:
    """Bound one branch's whole subtree in a disk centered at its attachment point.

    The attachment point is `closing_pair[0]`'s own coordinate -- what a
    sibling's backbone bond actually connects to, and (via
    `geometry_helpers.stem_base_for_attachment`) exactly where the
    placement pass puts it, not the ladder's rung-center axis, which sits
    `PAIR_SPACE / 2` off to one side of it (see that function's docstring
    for why the distinction is load-bearing).

    Dispatches on the branch point at the bottom of `closing_pair`'s
    collapsed stem run: a genuine multiloop (2+ children) or a terminal
    hairpin loop (0 children) uses the circular bounding-disk envelope
    (`_loop_branch_reach`); a bulge/interior loop (exactly 1 child) uses a
    tighter, LINEAR straight-continuation bound (`_bulge_branch_reach`) --
    the circular formula compounds `~3x` per nesting level on a long bulge
    CHAIN (see `engine.py`'s `_MAX_REACH` docstring for the blow-up this
    fixes), while the straight bound only adds a constant per level.

    Args:
        tree: The structure tree containing `closing_pair`.
        closing_pair: The branch's outermost pair.
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin
            (`RHO_MARGIN` included) -- S8's bounded re-inflation schedule.

    Returns:
        `reach(branch)`: every RHO-inflated primitive in `branch`'s subtree
        lies within this radius of its attachment point (see the module
        docstring's soundness argument for why this radius is what a
        parent packer needs as the branch's slot half-width).
    """
    cached = cache.branch_reach.get(closing_pair)
    if cached is not None:
        return cached
    depth, loop = collapse_stem(tree, closing_pair)
    if len(loop.children) == 1:
        reach = _bulge_branch_reach(tree, loop, depth, params, cache, margin_scale)
    else:
        reach = _loop_branch_reach(tree, loop, depth, params, cache, margin_scale)
    cache.branch_reach[closing_pair] = reach
    return reach


def _loop_branch_reach(
    tree: StructureTree,
    loop: Loop,
    depth: int,
    params: DrawParameters,
    cache: ReachCache,
    margin_scale: float,
) -> float:
    """`branch_reach`'s circular case: a multiloop (2+ children) or a
    terminal hairpin loop (0 children), packed on `loop_packing`'s circle.

    By triangle inequality: any point on the collapsed stem's far strand is
    within `PAIR_SPACE` of the near strand at the same rung, so within
    `PAIR_SPACE / 2` of the attachment point BEYOND the near strand's own
    on-axis distance; the loop this stem closes sits `packing.radius`
    beyond the stem's tip along that same axis, and its own members (or a
    child's own attachment point) lie within another `packing.radius` of
    that loop's center; and a child's own subtree is bounded, recursively,
    within `branch_reach(child)` of ITS attachment point. Summing every
    term (a safe, deliberately non-tight over-approximation) gives the
    formula below.

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop at the bottom of the collapsed stem run.
        depth: The collapsed stem's stacked-pair count.
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        The branch's reach.
    """
    packing = loop_packing(tree, loop, params, cache, margin_scale)
    child_reach = max(
        (
            branch_reach(tree, child.closing_pair, params, cache, margin_scale)
            for child in loop.children
        ),
        default=0.0,
    )
    return (
        (depth - 1) * params.PRIMARY_SPACE
        + params.PAIR_SPACE / 2.0
        + 2.0 * packing.radius
        + child_reach
        + RHO_MARGIN * margin_scale
    )


def bulge_split(loop: Loop, child: Branch) -> tuple[int, int]:
    """Split a single-child loop's interior into `(n_before, n_after)`.

    Args:
        loop: A loop with exactly one child branch.
        child: That one child branch.

    Returns:
        `(n_before, n_after)`: unpaired member counts before and after
        `child.start` in `loop.members[1:-1]` (structure order).
    """
    interior = loop.members[1:-1]
    index = interior.index(child.start)
    return index, len(interior) - index - 1


def _bulge_branch_reach(
    tree: StructureTree,
    loop: Loop,
    depth: int,
    params: DrawParameters,
    cache: ReachCache,
    margin_scale: float,
) -> float:
    """`branch_reach`'s straight-continuation case: a bulge/interior loop
    (exactly 1 child), placed by `geometry_helpers.place_bulge_geometry`
    instead of a circular envelope.

    Same triangle-inequality style as `_loop_branch_reach`, but the "beyond
    the tip" term is `packing.steps * PRIMARY_SPACE + PAIR_SPACE / 2`
    (the straight-line axial reach of the bulge's own unpaired content,
    from `geometry_helpers.pack_bulge_linear`) instead of `2 *
    packing.radius` -- LINEAR in the bulge's own size, not compounded by a
    circular `radius_floor`. Composed over a chain of `N` such loops, total
    reach is `O(N)` instead of the circular formula's `O(3^N)`.

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop at the bottom of the collapsed stem run (exactly 1 child).
        depth: The collapsed stem's stacked-pair count.
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        The branch's reach.
    """
    child = loop.children[0]
    n_before, n_after = bulge_split(loop, child)
    packing = pack_bulge_linear(n_before, n_after)
    child_reach = branch_reach(tree, child.closing_pair, params, cache, margin_scale)
    return (
        (depth - 1) * params.PRIMARY_SPACE
        + params.PAIR_SPACE / 2.0
        + packing.steps * params.PRIMARY_SPACE
        + params.PAIR_SPACE / 2.0
        + child_reach
        + RHO_MARGIN * margin_scale
    )


def loop_packing(
    tree: StructureTree,
    loop: Loop,
    params: DrawParameters,
    cache: ReachCache,
    margin_scale: float = 1.0,
) -> LoopPacking:
    """Pack one interior loop's members (unpaired nts + child branches).

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop to pack; must have a closing pair (the exterior loop
            has no circle to pack -- see `engine._place_exterior`, which
            uses `pack_line_positions` instead).
        params: Target geometry.
        cache: Memoization; mutated in place (also fills in every child
            branch's `branch_reach` as a side effect, recursively).
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        The loop's `LoopPacking` (radius + one angle per interior member).

    Raises:
        ValueError: If `loop.closing_pair is None` (the exterior loop).
    """
    if loop.closing_pair is None:
        raise ValueError("loop_packing needs a closing pair; pack the exterior loop separately")
    cached = cache.loop_packing.get(loop.closing_pair)
    if cached is not None:
        return cached
    reserved = stem_rung_half_width(params, margin_scale)
    slots = _loop_slots(tree, loop, params, cache, margin_scale)
    radius_floor = max(reserved, params.NODE_R, *slots) + _RADIUS_SEARCH_SLACK
    packing = pack_loop_angles(slots, reserved, radius_floor, _RADIUS_STEP)
    cache.loop_packing[loop.closing_pair] = packing
    return packing


def _loop_slots(
    tree: StructureTree,
    loop: Loop,
    params: DrawParameters,
    cache: ReachCache,
    margin_scale: float,
) -> list[float]:
    """Each interior member's required chord-clearance half-width.

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop being packed.
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        One half-width per `loop.members[1:-1]` entry, in order: a child
        branch's own `branch_reach`, or a bare unpaired member's fixed
        `unpaired_half_width`.
    """
    branch_by_start = {branch.start: branch for branch in loop.children}
    slots = []
    for member in loop.members[1:-1]:
        branch = branch_by_start.get(member)
        if branch is None:
            slots.append(unpaired_half_width(params, margin_scale))
        else:
            slots.append(branch_reach(tree, branch.closing_pair, params, cache, margin_scale))
    return slots


__all__ = [
    "RHO_MARGIN",
    "ReachCache",
    "stem_rung_half_width",
    "unpaired_half_width",
    "branch_reach",
    "bulge_split",
    "loop_packing",
]
