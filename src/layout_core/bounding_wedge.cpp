/**
 * @file bounding_wedge.cpp
 * @brief Implements `bounding_wedge.hpp`, ported from `boundingWedge.inc`.
 */

#include "rna_layout/bounding_wedge.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "rna_layout/bounding_boxes.hpp"
#include "rna_layout/config_tree.hpp"
#include "rna_layout/geometry.hpp"

namespace rna_layout {

namespace {

/// The node-relative angle of @p node about @p root, following the same
/// `getChildAngle`-equivalent math `getBoundingWedgeRec` inlines for
/// performance (`boundingWedge.inc:69-109`) rather than calling
/// `get_child_angle` on every recursive step.
double node_angle_about_root(const TreeNode& root, const TreeNode& node, double parent_angle,
                             Vec2 center_root, Vec2 v_root_node, double& min_angle,
                             double& max_angle) {
  const TreeNode* parent = node.parent;
  if (parent == &root) {
    // Only for the initial call, not the recursive ones: initialize
    // min/max with the direct child's angle.
    const double node_angle = get_child_angle(root, node);
    min_angle = node_angle;
    max_angle = node_angle;
    return node_angle;
  }

  const Vec2 center_parent = get_loop_center(*parent);
  const Vec2 v_root_parent = geom::vector_from_to(center_root, center_parent);
  double diff_parent = geom::angle_between(v_root_parent, v_root_node);
  if (!geom::is_to_the_right_point_vector(center_root, v_root_parent, get_loop_center(node))) {
    diff_parent *= -1;
  }
  return parent_angle + diff_parent;
}

/// The points of interest a bounding wedge must contain: every bulge peak
/// of `node`'s stem (each offset an extra `distance` further out, matching
/// the resolver's clearance), plus -- only for a direct child of `root` --
/// the two bottom corners of that child's stem where it meets `root`'s
/// loop. Ported from `boundingWedge.inc:111-162`.
std::vector<Vec2> wedge_points_of_interest(const TreeNode& root, const TreeNode& node,
                                           double distance) {
  // `node` is always a descendant reached via `root.children` (never the
  // config tree's true exterior root -- see `bounding_wedge.hpp`'s `root`
  // naming note), so it always has a `sbox` on an already-built (T1) tree.
  const StemBox& stem_node = *node.sbox;  // NOLINT(bugprone-unchecked-optional-access)

  std::vector<Vec2> points;
  points.reserve(stem_node.bulges.size() + 2);
  for (std::size_t i = 0; i < stem_node.bulges.size(); ++i) {
    points.push_back(bulge_coordinates_extra_distance(stem_node, static_cast<int>(i), distance).at);
  }

  if (node.parent == &root) {
    points.push_back(Vec2{
        stem_node.c.x - stem_node.e.x * stem_node.a.x + stem_node.e.y * stem_node.b.x,
        stem_node.c.y - stem_node.e.x * stem_node.a.y + stem_node.e.y * stem_node.b.y,
    });
    points.push_back(Vec2{
        stem_node.c.x - stem_node.e.x * stem_node.a.x - stem_node.e.y * stem_node.b.x,
        stem_node.c.y - stem_node.e.x * stem_node.a.y - stem_node.e.y * stem_node.b.y,
    });
  }
  return points;
}

/// Widen `[min_angle, max_angle]` to include `node_angle + diff_angle` for
/// every `diff_angle` in @p diffs. Shared tail of the two update passes
/// `getBoundingWedgeRec` runs (tangent-touch angles, then points of
/// interest) -- both do exactly this, so it is factored out once (not a
/// fidelity concern: both call sites are the same min/max compare).
template <typename Range>
void widen_angle_range(double node_angle, const Range& diffs, double& min_angle,
                       double& max_angle) {
  for (double diff_angle : diffs) {
    const double point_angle = node_angle + diff_angle;
    if (point_angle < min_angle) {
      min_angle = point_angle;
    }
    if (point_angle > max_angle) {
      max_angle = point_angle;
    }
  }
}

/// Ported from `getBoundingWedgeRec` (`boundingWedge.inc:40`).
void bounding_wedge_recursive(const TreeNode& root, const TreeNode& node, double parent_angle,
                              double clearance, double& min_angle, double& max_angle) {
  const double distance = geom::epsilon_fix(clearance);

  const Vec2 center_root = get_loop_center(root);
  const Vec2 center_node = get_loop_center(node);
  const Vec2 v_root_node = geom::vector_from_to(center_root, center_node);

  const double node_angle = node_angle_about_root(root, node, parent_angle, center_root,
                                                  v_root_node, min_angle, max_angle);

  // Tangent-touch angles: the two points where lines from `center_root`
  // touch `node`'s loop circle (grown by `distance`) guarantee the whole
  // loop is contained in the wedge.
  const LoopBox& loop_node = *node.lbox;  // NOLINT(bugprone-unchecked-optional-access) -- see above
  const double radius_node = loop_node.radius + distance;
  const double distance_root_node = geom::length(v_root_node);
  const double angle1 = std::asin(radius_node / distance_root_node);
  const double angle2 = -angle1;
  widen_angle_range(node_angle, std::array<double, 2>{angle1, angle2}, min_angle, max_angle);

  // Bulge points (and, for a direct child of root, its stem's bottom
  // corners).
  const std::vector<Vec2> points = wedge_points_of_interest(root, node, distance);
  std::vector<double> point_diffs;
  point_diffs.reserve(points.size());
  for (const Vec2& point : points) {
    const Vec2 v_center_point = geom::vector_from_to(center_root, point);
    double diff_angle = geom::angle_between(v_root_node, v_center_point);
    const double sign =
        geom::is_to_the_right_point_vector(center_root, v_root_node, point) ? 1.0 : -1.0;
    diff_angle *= sign;
    point_diffs.push_back(diff_angle);
  }
  widen_angle_range(node_angle, point_diffs, min_angle, max_angle);

  for (const auto& child : node.children) {
    bounding_wedge_recursive(root, *child, node_angle, clearance, min_angle, max_angle);
  }
}

}  // namespace

AngleRange bounding_wedge(const TreeNode& root, int child_index, double clearance) {
  const TreeNode& child = *root.children.at(static_cast<std::size_t>(child_index));
  AngleRange range;
  bounding_wedge_recursive(root, child, /*parent_angle=*/0.0, clearance, range.min_angle,
                           range.max_angle);
  return range;
}

}  // namespace rna_layout
