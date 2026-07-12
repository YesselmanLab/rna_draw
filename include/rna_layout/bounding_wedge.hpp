#pragma once

/**
 * @file bounding_wedge.hpp
 * @brief The bounding-wedge predicate, ported from `boundingWedge.inc`.
 *
 * A bounding wedge is the angular range (about @p root's loop center) that
 * fully contains one of @p root's subtrees -- the resolver's ancestor-
 * intersection check (Milestone A step 8, `handleAncestorIntersections.inc`)
 * uses it to bound how far a descendant subtree reaches around its
 * ancestor. Ported now (Milestone A step 4, alongside the config tree/boxes
 * it reads) as a pure, standalone predicate over an already-built tree --
 * it is not yet wired into any pipeline; see the plan's port order.
 *
 * `root` here is a NAMING artifact of the vendored function, not
 * necessarily the config tree's actual root: `getBoundingWedge` is called
 * by the resolver on the ancestor node currently under test, which is
 * usually an interior node, not the exterior loop. Preserves the
 * reference's arithmetic and operation order (Task 1's fidelity rule).
 */

#include "rna_layout/tree.hpp"

namespace rna_layout {

/// The angular range `[min_angle, max_angle]` (radians, about `root`'s loop
/// center) a bounding wedge spans; may extend outside `[0, 2*pi)`.
struct AngleRange {
  double min_angle = 0.0;
  double max_angle = 0.0;
};

/**
 * The bounding wedge of @p root's @p child_index-th child's whole subtree.
 * Ported from `getBoundingWedge` (`boundingWedge.inc:249`).
 *
 * @param root The node the wedge is measured about; must have a `sbox`
 *     (i.e. must not be the config tree's true exterior root -- see this
 *     header's `root` naming note).
 * @param child_index Index into `root.children`.
 * @param clearance `PuzzlerOptions::clearance`; threads into the vendored
 *     `epsilonFix` tolerance (`geom::epsilon_fix`).
 */
[[nodiscard]] AngleRange bounding_wedge(const TreeNode& root, int child_index, double clearance);

}  // namespace rna_layout
