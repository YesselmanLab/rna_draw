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

from .geometry_helpers import LoopPacking, pack_bulge_linear, pack_loop_angles, pack_two_seam_angles

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

# `_degree2_packing` only ATTEMPTS the two-seam (collinear) candidate for a
# loop that sits in a long enough run of consecutive degree-2 (and
# straight-through bulge/interior, which don't themselves compound but
# don't break the run either) loops -- otherwise it returns the full-circle
# packing untouched, so a loop nowhere near the exponential-blowup regime
# this fix targets gets BYTE-IDENTICAL geometry to before the fix, not just
# an equally-sized one. This matters because pinning changes WHICH ANGLE
# every sibling on the circle lands at (not just the radius): even a
# coincidental tie (`_degree2_packing`'s `<=`) reshuffles a loop's own
# children's angles, and that reshuffle propagates to every ancestor's own
# packing (a smaller `branch_reach` changes that ancestor's own slot size
# too) -- harmless for a loop actually escaping exponential blowup, but for
# an UNRELATED loop elsewhere in the same structure, a reshuffled angle can
# coincidentally swing a subtree closer to a totally different, distant
# part of the tree that was never a collision risk before (measured: 19 of
# the 436 previously-clean hard-set structures newly failed the checker
# with unconditional pinning; none of them had ANY degree-2 run longer than
# 3 anywhere in the whole structure).
#
# A RADIUS-based floor (tried first) does not work here: the compounding is
# ~3x PER LEVEL, so a short run's own radius is still small at every level
# except (if any) the outermost one or two -- gating on radius alone either
# never fires early enough in a short-but-genuinely-exponential run to
# provide real benefit, or has to be set so low it re-triggers on isolated,
# harmless degree-2 loops. Run LENGTH is a direct, purely topological
# (reach-independent) proxy for "will this compound" and does not have that
# problem: it is the same whether or not any pinning has happened yet.
#
# THRESHOLD, measured directly on `benchmarks/hard_set.json`: the run
# length (`_degree2_chain_length`, counted through straight-continuation
# bulges too) of every currently-clean structure's longest degree-2 run
# tops out at 5, while the shortest of the three SMALL exponential-blowup
# rejects (`bpRNA_RFAM_35151`, 754nt) reaches 6 -- `6` is the smallest
# integer that fixes those three (`bpRNA_RFAM_35151/40000/40002`) without
# ever firing on any currently-clean structure. The eleven LARGER rejects
# (~2860-3024nt) top out at run length 3 -- their blowup is dominated by a
# few 3+/4+/7-way junctions aggregating several only-modestly-deep degree-2
# runs, not one single long run (see `engine.py`'s `_MAX_REACH` docstring);
# fixing THAT is out of this fix's scope (multi-branch aggregation, not a
# straight-continuation case), so those eleven remain errors-not-dirty,
# unchanged by this fix, same as before it.
_DEGREE2_CHAIN_LENGTH_FLOOR = 6


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
        lateral_reach: `lateral_reach(branch)` for every branch closing pair
            seen so far -- an ANISOTROPIC bound (the subtree's perpendicular
            half-width about its own attachment axis, not an isotropic disk
            radius; see `lateral_reach`'s docstring). Separate from
            `branch_reach` because the two bound different things and a
            degree-2 loop's own `lateral_reach` needs both.
        degree2_pinned: Whether a degree-2 loop's `loop_packing` actually
            pinned its dominant child collinear at angle `pi`
            (`_degree2_packing`) rather than falling back to the ordinary
            full-circle packing (see that function's docstring for why it
            sometimes doesn't pay off) -- `lateral_reach` reads this to know
            which formula is sound for that loop.
        degree2_chain_length: `_degree2_chain_length(loop)` for every loop
            seen so far -- a purely topological (reach-independent) count
            of how long a run of consecutive degree-2 loops `loop` sits in,
            used to scope pinning to loops that can actually benefit from
            it (see `_DEGREE2_CHAIN_LENGTH_FLOOR`).
    """

    branch_reach: dict[tuple[int, int], float] = field(default_factory=dict)
    loop_packing: dict[tuple[int, int], LoopPacking] = field(default_factory=dict)
    lateral_reach: dict[tuple[int, int], float] = field(default_factory=dict)
    degree2_pinned: dict[tuple[int, int], bool] = field(default_factory=dict)
    degree2_chain_length: dict[tuple[int, int], int] = field(default_factory=dict)


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


def lateral_reach(
    tree: StructureTree,
    closing_pair: tuple[int, int],
    params: DrawParameters,
    cache: ReachCache,
    margin_scale: float = 1.0,
) -> float:
    """Bound one branch's subtree perpendicular to its OWN attachment axis.

    Unlike `branch_reach` (an isotropic disk radius -- the subtree can
    extend that far in ANY direction from its attachment point), this is a
    directional (anisotropic) bound: the subtree stays within
    `lateral_reach` of the RAY from the attachment point along the
    branch's own outward `axis_dir`, but may extend arbitrarily far ALONG
    that ray (up to `branch_reach`, which remains the bound on axial
    extent). Tighter than `branch_reach` only for a genuinely straight
    (chain-shaped) subtree -- a bulge/interior loop (1 child) or a degree-2
    multiloop (2 children) whose dominant child is placed as a collinear
    straight continuation (see `_degree2_packing`) -- because only those
    stay confined to a corridor around their own axis; a hairpin (0
    children) or a genuine (3+-way) multiloop spreads isotropically, so its
    `lateral_reach` is just its `branch_reach`.

    This is what lets a long CHAIN of degree-2 multiloops (one dominant
    "continuing" branch plus one small side branch, repeated many nesting
    levels deep, common in real rRNA) size its own loop-circle radius by
    the dominant child's `lateral_reach` instead of its full isotropic
    `branch_reach` (`_degree2_packing`) -- collapsing the circular
    envelope's `~3x`-per-level compounding (`_loop_branch_reach`'s
    docstring) to linear, the same way `_bulge_branch_reach` already does
    for single-child chains.

    Args:
        tree: The structure tree containing `closing_pair`.
        closing_pair: The branch's outermost pair.
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        A radius such that the branch's entire subtree lies within it of
        the ray from its attachment point along its own outward axis.
    """
    cached = cache.lateral_reach.get(closing_pair)
    if cached is not None:
        return cached
    _depth, loop = collapse_stem(tree, closing_pair)
    nch = len(loop.children)
    if nch == 1:
        result = _bulge_lateral_reach(tree, loop, params, cache, margin_scale)
    elif nch == 2:
        result = _degree2_lateral_reach(tree, closing_pair, loop, params, cache, margin_scale)
    else:
        result = branch_reach(tree, closing_pair, params, cache, margin_scale)
    cache.lateral_reach[closing_pair] = result
    return result


def _bulge_lateral_reach(
    tree: StructureTree, loop: Loop, params: DrawParameters, cache: ReachCache, margin_scale: float
) -> float:
    """`lateral_reach`'s single-child case: already a straight rail.

    A bulge/interior loop's unpaired members sit a fixed, small offset
    (`stem_rung_half_width`'s own rail-clearance formula) to either side of
    the SAME axis its child continues along (`_place_bulge`); the child's
    own subtree can only widen this further via its own `lateral_reach`.

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop at the bottom of the collapsed stem run (1 child).
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        The branch's `lateral_reach`.
    """
    child = loop.children[0]
    rail_half_width = stem_rung_half_width(params, margin_scale)
    child_lateral = lateral_reach(tree, child.closing_pair, params, cache, margin_scale)
    return max(rail_half_width, child_lateral)


def _dominant_and_side(
    tree: StructureTree,
    loop: Loop,
    params: DrawParameters,
    cache: ReachCache,
    margin_scale: float,
) -> tuple[Branch, Branch]:
    """A degree-2 loop's dominant (max-`branch_reach`) child and its sibling.

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop being packed (exactly 2 children).
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        `(dom, side)`.
    """
    first, second = loop.children
    first_reach = branch_reach(tree, first.closing_pair, params, cache, margin_scale)
    second_reach = branch_reach(tree, second.closing_pair, params, cache, margin_scale)
    return (first, second) if first_reach >= second_reach else (second, first)


def _degree2_lateral_reach(
    tree: StructureTree,
    closing_pair: tuple[int, int],
    loop: Loop,
    params: DrawParameters,
    cache: ReachCache,
    margin_scale: float,
) -> float:
    """`lateral_reach`'s degree-2 case: sound only if `dom` was pinned.

    `loop_packing` (via `_degree2_packing`) only pins the dominant child
    `dom` collinear at angle `pi` when doing so does not INFLATE the loop's
    own radius versus the ordinary full-circle packing (see that
    function's docstring). When it didn't pin, this loop packed exactly
    like a genuine multiloop, so only the isotropic `branch_reach` bound is
    sound here too. When it did, the loop's own subtree stays within
    `dom`'s own `lateral_reach` of the (now-collinear) axis beyond the
    loop, OR within the loop's circle (`packing.radius`) plus the side
    child's own isotropic `branch_reach`, whichever is larger.

    Args:
        tree: The structure tree containing `loop`.
        closing_pair: `loop`'s own closing pair (the branch's attachment).
        loop: The loop at the bottom of the collapsed stem run (2 children).
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        The branch's `lateral_reach`.
    """
    packing = loop_packing(tree, loop, params, cache, margin_scale)
    assert loop.closing_pair is not None  # `loop_packing` above already enforced this
    if not cache.degree2_pinned.get(loop.closing_pair, False):
        return branch_reach(tree, closing_pair, params, cache, margin_scale)
    dom, side = _dominant_and_side(tree, loop, params, cache, margin_scale)
    dom_lateral = lateral_reach(tree, dom.closing_pair, params, cache, margin_scale)
    side_reach = branch_reach(tree, side.closing_pair, params, cache, margin_scale)
    return max(dom_lateral, packing.radius + side_reach + RHO_MARGIN * margin_scale)


def loop_packing(
    tree: StructureTree,
    loop: Loop,
    params: DrawParameters,
    cache: ReachCache,
    margin_scale: float = 1.0,
) -> LoopPacking:
    """Pack one interior loop's members (unpaired nts + child branches).

    A degree-2 loop (exactly 2 children) additionally tries pinning its
    dominant child collinear at angle `pi` (`_degree2_packing`), keeping
    whichever candidate packing gives the smaller radius -- so this never
    packs WORSE than the ordinary full-circle packer used here for every
    other loop shape.

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
    if len(loop.children) == 2:
        packing = _degree2_packing(tree, loop, params, cache, margin_scale)
    else:
        packing = _full_circle_packing(tree, loop, params, cache, margin_scale)
    cache.loop_packing[loop.closing_pair] = packing
    return packing


def _full_circle_packing(
    tree: StructureTree, loop: Loop, params: DrawParameters, cache: ReachCache, margin_scale: float
) -> LoopPacking:
    """The ordinary full-circle candidate: every interior slot isotropic.

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop to pack.
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        The full-circle `LoopPacking`.
    """
    reserved = stem_rung_half_width(params, margin_scale)
    slots = _loop_slots(tree, loop, params, cache, margin_scale)
    radius_floor = max(reserved, params.NODE_R, *slots) + _RADIUS_SEARCH_SLACK
    return pack_loop_angles(slots, reserved, radius_floor, _RADIUS_STEP)


def _split_slots_around(
    loop: Loop, dom: Branch, slots: list[float]
) -> tuple[list[float], list[float]]:
    """Split a degree-2 loop's interior slots into before/after `dom`.

    Args:
        loop: The loop being packed (exactly 2 children).
        dom: The dominant child (`_dominant_and_side`); excluded from both
            returned lists (its own slot is sized by `lateral_reach`, not
            `slots`' isotropic `branch_reach` entry).
        slots: `_loop_slots`'s per-interior-member half-widths, in
            structure order (same order as `loop.members[1:-1]`).

    Returns:
        `(before, after)`: `slots` split at `dom`'s own index, `dom`'s own
        entry dropped.
    """
    interior = loop.members[1:-1]
    index = interior.index(dom.start)
    return slots[:index], slots[index + 1 :]


def _two_seam_packing(
    tree: StructureTree,
    loop: Loop,
    dom: Branch,
    params: DrawParameters,
    cache: ReachCache,
    margin_scale: float,
) -> LoopPacking:
    """The degree-2 collinear-continuation candidate: `dom` pinned at `pi`.

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop being packed (exactly 2 children).
        dom: The dominant child, pinned collinear (`_dominant_and_side`).
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        The two-seam `LoopPacking` (see `geometry_helpers.pack_two_seam_angles`).
    """
    reserved = stem_rung_half_width(params, margin_scale)
    slots = _loop_slots(tree, loop, params, cache, margin_scale)
    before, after = _split_slots_around(loop, dom, slots)
    dom_half_width = lateral_reach(tree, dom.closing_pair, params, cache, margin_scale)
    floor_terms = (reserved, params.NODE_R, dom_half_width, *before, *after)
    radius_floor = max(floor_terms) + _RADIUS_SEARCH_SLACK
    return pack_two_seam_angles(before, after, dom_half_width, reserved, radius_floor, _RADIUS_STEP)


def _degree2_chain_length(tree: StructureTree, loop: Loop, cache: ReachCache) -> int:
    """Length of the longest run of consecutive degree-2 loops through `loop`.

    Purely topological (no `params`/reach values involved, so it is stable
    regardless of any pinning decision made elsewhere) -- see
    `_DEGREE2_CHAIN_LENGTH_FLOOR`'s comment for why this, not a reach or
    radius value, is the sound way to scope pinning to loops that actually
    need it. A straight-continuation bulge/interior loop (1 child) does not
    itself compound (`_bulge_branch_reach` is already linear), but does not
    BREAK a run either -- the run continues through it to whatever is on
    the other side. A hairpin (0 children) or a genuine (3+-way) multiloop
    ends a run (its own children spread isotropically, not collinearly).

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop to measure.
        cache: Memoization; mutated in place.

    Returns:
        `0` if `loop` is not itself degree-2; otherwise `1 +` the longer of
        its two children's own runs.
    """
    if loop.closing_pair is None:
        return 0
    cached = cache.degree2_chain_length.get(loop.closing_pair)
    if cached is not None:
        return cached
    nch = len(loop.children)
    if nch == 1:
        _depth, child_loop = collapse_stem(tree, loop.children[0].closing_pair)
        result = _degree2_chain_length(tree, child_loop, cache)
    elif nch == 2:
        runs = []
        for child in loop.children:
            _depth, child_loop = collapse_stem(tree, child.closing_pair)
            runs.append(_degree2_chain_length(tree, child_loop, cache))
        result = 1 + max(runs)
    else:
        result = 0
    cache.degree2_chain_length[loop.closing_pair] = result
    return result


def _degree2_run_top(tree: StructureTree, loop: Loop) -> Loop:
    """Walk up to the outermost loop in `loop`'s own degree-2/bulge run.

    `_degree2_chain_length` only counts DOWNWARD (how much run is left
    BENEATH a loop) -- a loop near the BOTTOM of a long run sees a short
    count there (as little as `1`), even though it is part of a run long
    enough to be worth pinning. This walks up to the run's own top FIRST,
    so `_degree2_chain_length` at THAT loop gives the run's TOTAL length,
    the same for every loop within it (see `_degree2_packing`'s gate).

    Purely topological, via `tree.enclosing_pairs` (no `params`/reach
    values, no recursion into `branch_reach`): a parent loop with 1 or 2
    children continues the same run `loop` is part of (a bulge parent has
    no other child to have arrived from; a degree-2 parent might reach
    `loop` via either child -- either way `loop` is still
    collinear-pinnable content of that parent's own run). Stops at a
    hairpin, a genuine (3+-way) multiloop, or the exterior, whichever ends
    the run first.

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop to walk up from.

    Returns:
        The outermost `Loop` in `loop`'s own run (`loop` itself if it has
        no degree-2/bulge parent).
    """
    if loop.closing_pair is None:
        return loop
    stack = tree.enclosing_pairs[loop.closing_pair[0]]
    if len(stack) < 2:
        return loop
    parent = tree.loop_by_closing_pair[stack[-2]]
    if len(parent.children) in (1, 2):
        return _degree2_run_top(tree, parent)
    return loop


def _degree2_total_run_length(tree: StructureTree, loop: Loop, cache: ReachCache) -> int:
    """The TOTAL length of the degree-2/bulge run `loop` belongs to.

    Args:
        tree: The structure tree containing `loop`.
        loop: Any loop within the run (need not be its top or bottom).
        cache: Memoization; mutated in place.

    Returns:
        `_degree2_chain_length` measured from the run's own outermost loop
        (`_degree2_run_top`) -- the same value for every loop in the run.
    """
    return _degree2_chain_length(tree, _degree2_run_top(tree, loop), cache)


def _degree2_packing(
    tree: StructureTree, loop: Loop, params: DrawParameters, cache: ReachCache, margin_scale: float
) -> LoopPacking:
    """A degree-2 loop's packing: pin unless it is unneeded or would enlarge the loop.

    Tries pinning the dominant child `dom` collinear at angle `pi`
    (`_two_seam_packing`), sized by `dom`'s tighter, directional
    `lateral_reach` instead of its isotropic `branch_reach` -- this is what
    makes a long CHAIN of degree-2 loops linear instead of exponential
    (`lateral_reach`'s docstring).

    Two guards keep this SCOPED to loops that actually need it, both
    leaving a loop's packing BYTE-IDENTICAL to the plain full-circle
    packer (used everywhere else) when they don't apply:
    - `_DEGREE2_CHAIN_LENGTH_FLOOR`: below this run length, always use the
      full-circle packing untouched -- pinning changes WHICH ANGLE every
      sibling on the circle lands at (not just `dom`'s), and that reshuffle
      propagates to every ancestor's own packing (a smaller `branch_reach`
      changes its slot size there too); harmless for a loop actually
      escaping exponential blowup, but for an unrelated, already-short run,
      an unnecessary reshuffle risks coincidentally swinging some subtree
      closer to a totally different, distant part of the tree that was
      never a collision risk before (see the constant's own comment for
      the measured regression this closes).
    - Even above the floor, pinning `dom` to an EXACT angle, rather than
      letting the ordinary packer place it wherever the additive slot
      order is cheapest, can lose on a lopsided loop (e.g. most of its
      OTHER interior content on one side of `dom`, forcing one of the two
      `<= pi` seam-to-seam arcs over budget even though the two-seam
      packing's own `dom` slot is tighter) -- so this still declines to
      pin whenever doing so would strictly ENLARGE the radius.

    `lateral_reach` (`_degree2_lateral_reach`) reads back which candidate
    won (`cache.degree2_pinned`) so it stays consistent with whatever this
    function actually chose.

    Args:
        tree: The structure tree containing `loop`.
        loop: The loop to pack (exactly 2 children).
        params: Target geometry.
        cache: Memoization; mutated in place.
        margin_scale: Multiplier on every additive clearance margin.

    Returns:
        The two-seam `LoopPacking` if it was tried and is no larger than
        the full-circle one, else the full-circle one.
    """
    assert loop.closing_pair is not None  # `loop_packing`'s caller already enforced this
    full_packing = _full_circle_packing(tree, loop, params, cache, margin_scale)
    if _degree2_total_run_length(tree, loop, cache) < _DEGREE2_CHAIN_LENGTH_FLOOR:
        cache.degree2_pinned[loop.closing_pair] = False
        return full_packing
    dom, _side = _dominant_and_side(tree, loop, params, cache, margin_scale)
    two_seam_packing = _two_seam_packing(tree, loop, dom, params, cache, margin_scale)
    pinned = two_seam_packing.radius <= full_packing.radius
    cache.degree2_pinned[loop.closing_pair] = pinned
    return two_seam_packing if pinned else full_packing


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
    "lateral_reach",
    "bulge_split",
    "loop_packing",
]
