/**
 * @file rotation_angle2.cpp
 * @brief The nine per-intersection-type rotation-angle solvers plus
 *        `get_rotation_angle`'s dispatcher, ported from `rotationAngle.inc`.
 *        Each solver reduces its specific box-pair shape to one call of
 *        `fix_intersection_of_rectangle_and_circle` or
 *        `fix_intersection_of_circles` (`rotation_angle.cpp`).
 *
 * FIDELITY NOTE: see `rotation_angle.cpp`'s file header -- the same
 * exact-expression-tree rule applies throughout.
 */

#include "resolve_internal.hpp"

#include "rna_layout/bounding_boxes.hpp"
#include "rna_layout/geometry.hpp"
#include "rna_layout/intersect_tree.hpp"

namespace rna_layout {

namespace {

double get_rotation_angle_lxl(const TreeNode& ancestor, const TreeNode& rotation_node,
                              const TreeNode& intersector, short rotation_sign, double clearance) {
  // Ported from `getRotationAngleLxL` (`rotationAngle.inc:304`).
  const LoopBox& static_loop = *ancestor.lbox;
  const LoopBox& rotation_loop = *rotation_node.lbox;
  const LoopBox& mobile_loop = *intersector.lbox;

  return fix_intersection_of_circles(static_loop.center, static_loop.radius, mobile_loop.center,
                                     mobile_loop.radius, rotation_loop.center, rotation_sign,
                                     clearance);
}

double get_rotation_angle_lxs(const TreeNode& ancestor, const TreeNode& rotation_node,
                              const TreeNode& intersector, short rotation_sign, double clearance) {
  // Ported from `getRotationAngleLxS` (`rotationAngle.inc:343`).
  const StemBox& static_rect = *intersector.sbox;
  const LoopBox& mobile_circ = *ancestor.lbox;
  const LoopBox& rotation_loop = *rotation_node.lbox;

  const short inverse_rotation_sign = static_cast<short>(-1 * rotation_sign);
  const double inverse_rotation_angle = fix_intersection_of_rectangle_and_circle(
      static_rect.c, static_rect.a, static_rect.b, static_rect.e.x, static_rect.e.y,
      mobile_circ.center, mobile_circ.radius, rotation_loop.center, inverse_rotation_sign,
      clearance);
  return -1 * inverse_rotation_angle;
}

double get_rotation_angle_sxl(const TreeNode& ancestor, const TreeNode& rotation_node,
                              const TreeNode& intersector, short rotation_sign, double clearance) {
  // Ported from `getRotationAngleSxL` (`rotationAngle.inc:376`).
  const StemBox& static_rect = *ancestor.sbox;
  const LoopBox& mobile_circ = *intersector.lbox;
  const LoopBox& rotation_loop = *rotation_node.lbox;

  return fix_intersection_of_rectangle_and_circle(static_rect.c, static_rect.a, static_rect.b,
                                                  static_rect.e.x, static_rect.e.y,
                                                  mobile_circ.center, mobile_circ.radius,
                                                  rotation_loop.center, rotation_sign, clearance);
}

double get_rotation_angle_lxb(const TreeNode& ancestor, const TreeNode& rotation_node,
                              const TreeNode& intersector, short rotation_sign, double clearance) {
  // Ported from `getRotationAngleLxB` (`rotationAngle.inc:407`): construct
  // circles around the intersecting loop and bulge and resolve their
  // intersection.
  const LoopBox& static_loop = *ancestor.lbox;
  const StemBox& mobile_stem = *intersector.sbox;

  const BulgeHit hit = intersect_loop_bulges(static_loop, mobile_stem, clearance);
  const BulgePoints mobile_bulge = bulge_coordinates(mobile_stem, hit.bulge);
  const geom::Circle mobile_circle =
      geom::circumcircle(mobile_bulge.prev, mobile_bulge.at, mobile_bulge.next);

  const LoopBox& rotation_loop = *rotation_node.lbox;

  return fix_intersection_of_circles(static_loop.center, static_loop.radius, mobile_circle.center,
                                     mobile_circle.radius, rotation_loop.center, rotation_sign,
                                     clearance);
}

double get_rotation_angle_bxl(const TreeNode& ancestor, const TreeNode& rotation_node,
                              const TreeNode& intersector, short rotation_sign, double clearance) {
  // Ported from `getRotationAngleBxL` (`rotationAngle.inc:471`): construct
  // circles around the intersecting bulge and loop and resolve their
  // intersection.
  const StemBox& static_stem = *ancestor.sbox;
  const LoopBox& mobile_loop = *intersector.lbox;

  const BulgeHit hit = intersect_loop_bulges(mobile_loop, static_stem, clearance);
  const BulgePoints static_bulge = bulge_coordinates(static_stem, hit.bulge);
  const geom::Circle static_circle =
      geom::circumcircle(static_bulge.prev, static_bulge.at, static_bulge.next);

  const LoopBox& rotation_loop = *rotation_node.lbox;

  return fix_intersection_of_circles(static_circle.center, static_circle.radius, mobile_loop.center,
                                     mobile_loop.radius, rotation_loop.center, rotation_sign,
                                     clearance);
}

double get_rotation_angle_sxs(const TreeNode& ancestor, const TreeNode& rotation_node,
                              const TreeNode& intersector, short rotation_sign, double clearance) {
  // Ported from `getRotationAngleSxS` (`rotationAngle.inc:535`): delegates
  // to the SxL solver unchanged.
  return get_rotation_angle_sxl(ancestor, rotation_node, intersector, rotation_sign, clearance);
}

double get_rotation_angle_sxb(const TreeNode& ancestor, const TreeNode& rotation_node,
                              const TreeNode& intersector, short rotation_sign, double clearance) {
  // Ported from `getRotationAngleSxB` (`rotationAngle.inc:554`).
  const StemBox& static_stem = *ancestor.sbox;
  const StemBox& mobile_stem = *intersector.sbox;
  const LoopBox& rotation_loop = *rotation_node.lbox;

  const BulgeHit hit = intersect_stem_bulges(static_stem, mobile_stem, clearance);
  const BulgePoints mobile_bulge = bulge_coordinates(mobile_stem, hit.bulge);
  const geom::Circle mobile_circle =
      geom::circumcircle(mobile_bulge.prev, mobile_bulge.at, mobile_bulge.next);

  return fix_intersection_of_rectangle_and_circle(
      static_stem.c, static_stem.a, static_stem.b, static_stem.e.x, static_stem.e.y,
      mobile_circle.center, mobile_circle.radius, rotation_loop.center, rotation_sign, clearance);
}

double get_rotation_angle_bxs(const TreeNode& ancestor, const TreeNode& rotation_node,
                              const TreeNode& intersector, short rotation_sign, double clearance) {
  // Ported from `getRotationAngleBxS` (`rotationAngle.inc:599`). NOTE the
  // vendored function's own swap: `staticStem` is the INTERSECTOR's box and
  // `mobileStem` is the ANCESTOR's box here (opposite of every other
  // XxS/SxX solver) -- preserved as written.
  const StemBox& static_stem = *intersector.sbox;
  const StemBox& mobile_stem = *ancestor.sbox;
  const LoopBox& rotation_loop = *rotation_node.lbox;

  const BulgeHit hit = intersect_stem_bulges(static_stem, mobile_stem, clearance);
  const BulgePoints mobile_bulge = bulge_coordinates(mobile_stem, hit.bulge);
  const geom::Circle mobile_circle =
      geom::circumcircle(mobile_bulge.prev, mobile_bulge.at, mobile_bulge.next);

  return fix_intersection_of_rectangle_and_circle(
      static_stem.c, static_stem.a, static_stem.b, static_stem.e.x, static_stem.e.y,
      mobile_circle.center, mobile_circle.radius, rotation_loop.center, rotation_sign, clearance);
}

double get_rotation_angle_bxb(const TreeNode& ancestor, const TreeNode& rotation_node,
                              const TreeNode& intersector, short rotation_sign, double clearance) {
  // Ported from `getRotationAngleBxB` (`rotationAngle.inc:644`): construct
  // circles around both bulges and resolve their intersection.
  const StemBox& static_stem = *ancestor.sbox;
  const StemBox& mobile_stem = *intersector.sbox;

  const BulgeBulgeHit hit = intersect_bulges_bulges(static_stem, mobile_stem, clearance);
  const BulgePoints static_bulge = bulge_coordinates(static_stem, hit.bulge1);
  const geom::Circle static_circle =
      geom::circumcircle(static_bulge.prev, static_bulge.at, static_bulge.next);

  const BulgePoints mobile_bulge = bulge_coordinates(mobile_stem, hit.bulge2);
  const geom::Circle mobile_circle =
      geom::circumcircle(mobile_bulge.prev, mobile_bulge.at, mobile_bulge.next);

  const LoopBox& rotation_loop = *rotation_node.lbox;

  return fix_intersection_of_circles(static_circle.center, static_circle.radius,
                                     mobile_circle.center, mobile_circle.radius,
                                     rotation_loop.center, rotation_sign, clearance);
}

}  // namespace

double get_rotation_angle(const TreeNode& ancestor, const TreeNode& rotation_node,
                          const TreeNode& intersector, IntersectionType it, short rotation_sign,
                          double clearance) {
  // Ported from `getRotationAngle` (`rotationAngle.inc:715`): dispatch on
  // intersection type; any other type (`none`/`siblings`/`exterior`, never
  // actually passed here) returns `0.0` (mirrors the vendored `default:`
  // case, which prints an error and falls through with `rotationAngle`
  // still at its `0.0` initializer).
  switch (it) {
    case IntersectionType::loop_loop:
      return get_rotation_angle_lxl(ancestor, rotation_node, intersector, rotation_sign, clearance);
    case IntersectionType::loop_stem:
      return get_rotation_angle_lxs(ancestor, rotation_node, intersector, rotation_sign, clearance);
    case IntersectionType::loop_bulge:
      return get_rotation_angle_lxb(ancestor, rotation_node, intersector, rotation_sign, clearance);
    case IntersectionType::stem_loop:
      return get_rotation_angle_sxl(ancestor, rotation_node, intersector, rotation_sign, clearance);
    case IntersectionType::stem_stem:
      return get_rotation_angle_sxs(ancestor, rotation_node, intersector, rotation_sign, clearance);
    case IntersectionType::stem_bulge:
      return get_rotation_angle_sxb(ancestor, rotation_node, intersector, rotation_sign, clearance);
    case IntersectionType::bulge_loop:
      return get_rotation_angle_bxl(ancestor, rotation_node, intersector, rotation_sign, clearance);
    case IntersectionType::bulge_stem:
      return get_rotation_angle_bxs(ancestor, rotation_node, intersector, rotation_sign, clearance);
    case IntersectionType::bulge_bulge:
      return get_rotation_angle_bxb(ancestor, rotation_node, intersector, rotation_sign, clearance);
    case IntersectionType::none:
    case IntersectionType::siblings:
    case IntersectionType::exterior:
    default:
      return 0.0;
  }
}

}  // namespace rna_layout
