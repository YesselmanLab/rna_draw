#pragma once

/**
 * @file resolve_exterior.hpp
 * @brief The exterior-loop child-spacing pass, ported from
 *        `resolveExteriorChildrenIntersectionXY`
 *        (`resolveExteriorChildIntersections.inc:86`) -- Milestone A step 6
 *        ("finalization with resolver OFF").
 *
 * Despite its name and its use of `intersect_trees` (an Milestone-A-step-5
 * DETECTION predicate), this is NOT part of "the resolver"
 * (`PuzzlerOptions::check_sibling`/`check_ancestor`/`optimize`, Milestone A
 * steps 7-9, which mutate `Config` angles): `RNApuzzler.c:498-507` calls it
 * UNCONDITIONALLY, regardless of those three options -- it is the last
 * stage of coordinate FINALIZATION, spacing (and, if `allow_flipping`,
 * flipping) the exterior loop's direct children left-to-right along the x
 * axis so they do not overlap each other, using simple rigid translation
 * (never touching any `Config`). That is exactly why the plan groups it
 * with `determine_nucleotide_coordinates` under "finalization with resolver
 * OFF": both run whether or not the sibling/ancestor/optimize resolver
 * passes ever do.
 *
 * SCOPE NOTE -- dead code NOT ported: `resolveExteriorChildIntersections.inc`
 * defines two more variants of this pass
 * (`resolveExteriorChildrenIntersectionAffin`,
 * `resolveExteriorChildIntersections`) plus a `getSimpleBoundingBox` helper;
 * none of the three is ever called anywhere in the vendored codebase
 * (verified by `grep` over the whole tree) -- only
 * `resolveExteriorChildrenIntersectionXY` is (from `RNApuzzler.c`). Not
 * ported.
 */

#include "rna_layout/tree.hpp"
#include "rna_layout/turtle.hpp"

namespace rna_layout {

/**
 * Space (and, if @p allow_flipping, flip) @p exterior_root's direct
 * children along the x axis so no two intersect, updating both their
 * `StemBox`/`LoopBox`/`Aabb` (via `translate_bounding_boxes`) and every
 * nucleotide coordinate @p coords already holds for them (via direct
 * `coords.x`/`coords.y` writes) -- see this file's header for why boxes
 * AND coordinates both need updating (the boxes drive this pass's own
 * `intersect_trees` checks; `coords` is the actual output). A no-op if
 * @p exterior_root has fewer than 2 children.
 *
 * Ported verbatim from `resolveExteriorChildrenIntersectionXY`
 * (`resolveExteriorChildIntersections.inc:86-298`), INCLUDING its `coords.x`/
 * `coords.y` index arithmetic (direct `[base]`, not `[base - 1]`, throughout
 * this function -- the SAME one-past-the-reference-nucleotide convention
 * `nucleotide_coords.cpp`'s `handle_loop` documents; preserved here for the
 * same reason: not a bug, a deliberate reference-implementation convention
 * this port reproduces exactly rather than "fixes").
 *
 * @param exterior_root The tree's root; its children's boxes are translated
 *     in place.
 * @param pair_table 1-indexed pair table.
 * @param unpaired Default backbone-step distance (the per-collision spacing
 *     increment).
 * @param allow_flipping Whether a colliding child may be flipped below the
 *     baseline instead of pushed further right.
 * @param clearance `PuzzlerOptions::clearance`, threaded into every
 *     `intersect_trees` call this pass makes.
 * @param coords The nucleotide coordinates to update in place (must already
 *     hold `determine_nucleotide_coordinates`'s output).
 */
void resolve_exterior_children_intersection(TreeNode& exterior_root,
                                            const std::vector<int>& pair_table, double unpaired,
                                            bool allow_flipping, double clearance, Coords& coords);

}  // namespace rna_layout
