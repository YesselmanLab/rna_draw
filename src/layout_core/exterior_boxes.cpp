/**
 * @file exterior_boxes.cpp
 * @brief `setup_exterior_bounding_boxes`, ported from the LIVE variant,
 *        `TENTATIVE3_setupExteriorBoundingBoxes`
 *        (`handleAncestorIntersections.inc:477`) -- see
 *        `resolve_internal.hpp`'s scope note on the dead
 *        `TENTATIVE`/`TENTATIVE2` variants this does NOT port.
 *
 * Rebuilds the exterior root's `LoopBox`/`StemBox`/`Aabb` as a SYNTHETIC
 * box spanning from a top-level ancestor's loop center out towards an
 * intersector, so `intersect_node_node` can test it like any other
 * box-bearing node (`resolve_ancestors.cpp`'s `check_node_against_ancestors`,
 * Milestone A step 8).
 */

#include "resolve_internal.hpp"

#include "rna_layout/bounding_boxes.hpp"
#include "rna_layout/config_tree.hpp"
#include "rna_layout/geometry.hpp"
#include "rna_layout/intersect_tree.hpp"

namespace rna_layout {

namespace {

/// Overwrite @p exterior's `LoopBox`/`StemBox`/`Aabb` with @p loop and a
/// stem box built from the given rectangle corners (south/west, north/west,
/// south/east). Ported from the LIVE variant,
/// `TENTATIVE2_updateExteriorBoundingBoxes`
/// (`handleAncestorIntersections.inc:443`) -- see this file's header for why
/// only this "2" variant (of three near-duplicates in the vendored file) is
/// ported.
void update_exterior_bounding_boxes(TreeNode& exterior, LoopBox loop, double stem_north_x,
                                    double stem_south_x, double stem_west_y, double stem_east_y) {
  const Vec2 s{stem_south_x, stem_west_y};
  const Vec2 e{stem_north_x, stem_west_y};
  const Vec2 sp{stem_south_x, stem_east_y};
  const StemBox stem = create_stem_box(s, e, sp);

  exterior.lbox = loop;
  exterior.sbox = stem;
  exterior.aabb = compute_aabb(stem, loop);
}

}  // namespace

void setup_exterior_bounding_boxes(TreeNode& exterior, const TreeNode& top_level_ancestor,
                                   const TreeNode& intersector, const PuzzlerOptions& opts) {
  const double upper_y = geom::kExteriorY;
  const double lower_y = upper_y - opts.paired;

  const double loop_x = top_level_ancestor.lbox->center.x;  // NOLINT(bugprone-unchecked-optional-access)

  const double radius = 0.5 * (upper_y - lower_y);
  const Vec2 center{loop_x, upper_y - radius};

  const Aabb& intersector_aabb = intersector.aabb;

  if (intersector_aabb.max.x < loop_x) {
    // Intersector is left of the top-level ancestor; use the distance
    // aabb->min.x .. loop_x for the stem setup.
    const double stem_north_x = loop_x;
    const double stem_south_x = intersector_aabb.min.x;
    const double stem_west_y = upper_y;
    const double stem_east_y = lower_y;
    update_exterior_bounding_boxes(exterior, create_loop_box(center, radius), stem_north_x,
                                   stem_south_x, stem_west_y, stem_east_y);
  } else if (loop_x < intersector_aabb.min.x) {
    // Intersector is right of the top-level ancestor; use the distance
    // loop_x .. aabb->max.x for the stem setup.
    const double stem_north_x = loop_x;
    const double stem_south_x = intersector_aabb.max.x;
    const double stem_west_y = lower_y;
    const double stem_east_y = upper_y;
    update_exterior_bounding_boxes(exterior, create_loop_box(center, radius), stem_north_x,
                                   stem_south_x, stem_west_y, stem_east_y);
  } else {
    // Intersector shares some space in the x direction with the top-level
    // ancestor: try the left-side setup first, then (if that does not
    // actually intersect) the right-side setup.
    const double stem_north_x = loop_x;
    const double stem_south_x = intersector_aabb.min.x;
    const double stem_west_y = upper_y;
    const double stem_east_y = lower_y;
    update_exterior_bounding_boxes(exterior, create_loop_box(center, radius), stem_north_x,
                                   stem_south_x, stem_west_y, stem_east_y);

    if (intersect_node_node(intersector, exterior, opts.clearance).type == IntersectionType::none) {
      const double stem_north_x2 = loop_x;
      const double stem_south_x2 = intersector_aabb.max.x;
      const double stem_west_y2 = lower_y;
      const double stem_east_y2 = upper_y;
      update_exterior_bounding_boxes(exterior, create_loop_box(center, radius), stem_north_x2,
                                     stem_south_x2, stem_west_y2, stem_east_y2);
    }
  }
}

}  // namespace rna_layout
