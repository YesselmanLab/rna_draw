/**
 * @file resolve_ancestors.cpp
 * @brief The ancestor-intersection resolver, ported from the LIVE half of
 *        `handleAncestorIntersections.inc` (989 LOC; see
 *        `resolve_internal.hpp`'s scope note on the six dead functions this
 *        does NOT port): `is_straight_interior_loop`,
 *        `construct_reduced_intersection_path`, `get_rotation_sign`,
 *        `fix_intersection_with_ancestor`, `handle_intersection_with_ancestor`,
 *        and `check_node_against_ancestors` -- `resolve.cpp`'s ancestor
 *        branch entry point (Milestone A step 8).
 *
 * FIDELITY NOTE (`.claude/plans/current-plan.md`'s "CRITICAL" section): this
 * is the highest-risk module in the whole port; every expression below
 * preserves the reference's exact operation order.
 */

#include "resolve_internal.hpp"

#include <cmath>

#include "rna_layout/config_tree.hpp"
#include "rna_layout/geometry.hpp"
#include "rna_layout/intersect_tree.hpp"

namespace rna_layout {

bool is_straight_interior_loop(const TreeNode& node) {
  return is_interior_loop(node) && get_child_angle_by_index(node, 0) == geom::kPi;
}

std::vector<TreeNode*> construct_reduced_intersection_path(TreeNode& ancestor,
                                                            TreeNode& intersector,
                                                            IntersectionType it) {
  // Ported from `constructReducedIntersectionPath`
  // (`handleAncestorIntersections.inc:311`).
  int path_length = 1;
  TreeNode* node = &intersector;
  while (node != &ancestor) {
    node = node->parent;
    if (!is_straight_interior_loop(*node)) {
      ++path_length;
    }
  }

  switch (it) {
    case IntersectionType::loop_loop:
    case IntersectionType::loop_stem:
    case IntersectionType::loop_bulge:
      // Start at the child of the ancestor node for Lx? intersections.
      if (!is_straight_interior_loop(ancestor)) {
        --path_length;
      }
      break;
    default:
      // Include the ancestor node for all other intersections (unless it is
      // itself a straight interior loop).
      ;
  }

  std::vector<TreeNode*> path(static_cast<std::size_t>(path_length));
  node = &intersector;
  for (int i = path_length - 1; i >= 0; node = node->parent) {
    if (i == path_length - 1 || !is_straight_interior_loop(*node)) {
      path[static_cast<std::size_t>(i)] = node;
      --i;
    }
  }
  return path;
}

short get_rotation_sign(const std::vector<TreeNode*>& path) {
  // Ported from the LIVE variant, `TENTATIVE2_getRotationSign`
  // (`handleAncestorIntersections.inc:101`) -- see `resolve_internal.hpp`'s
  // scope note.
  if (path.size() < 2) {
    return 0;
  }

  double angle = 0.0;
  TreeNode* current_node = path[0];
  for (std::size_t i = 1; i < path.size(); ++i) {
    TreeNode* child_node = path[i];
    angle += get_child_angle(*current_node, *child_node);
    angle -= geom::kPi;
    current_node = child_node;
  }

  if (angle < 0.0) {
    return 1;
  }
  if (angle > 0.0) {
    return -1;
  }
  return 0;
}

TreeNode* fix_intersection_with_ancestor(TreeNode& ancestor, TreeNode& rotation_node,
                                         TreeNode& intersector, int rotation_index,
                                         short rotation_sign, IntersectionType it,
                                         const PuzzlerOptions& opts, ResolverState& state) {
  // Ported from `fixIntersectionWithAncestor`
  // (`handleAncestorIntersections.inc:214`).
  if (&rotation_node == &ancestor &&
      (it == IntersectionType::loop_loop || it == IntersectionType::loop_stem ||
       it == IntersectionType::loop_bulge)) {
    return nullptr;
  }

  double internal_child_angle = 0.0;
  if (is_interior_loop(rotation_node)) {
    // Prevent interior loops from increasing the distance to their
    // "straight" state -- only allow rotations towards straight.
    internal_child_angle = get_child_angle_by_index(rotation_node, 0);
    short allowed_rotation_sign = 0;
    if (internal_child_angle > geom::kPi) {
      allowed_rotation_sign = -1;
    } else if (internal_child_angle < geom::kPi) {
      allowed_rotation_sign = 1;
    }
    if (rotation_sign != allowed_rotation_sign) {
      return nullptr;
    }
  }

  double rotation_angle =
      get_rotation_angle(ancestor, rotation_node, intersector, it, rotation_sign, opts.clearance);

  if (is_interior_loop(rotation_node)) {
    // Prevent interior loops from rotating over their "straight" state --
    // limit the rotation so this interior loop becomes straight, at most.
    const double diff_to_straight = geom::kPi - internal_child_angle;
    if (std::fabs(rotation_angle) > std::fabs(diff_to_straight)) {
      rotation_angle = diff_to_straight;
    }
  }

  bool changed = false;
  if (rotation_angle != 0.0) {
    const double delta_angle = std::fabs(rotation_angle);
    int index_left = -2;
    int index_right = -2;
    if (rotation_angle > 0.0) {
      index_left = -1;
      index_right = rotation_index;
    } else {
      index_left = rotation_index;
      index_right = -1;
    }

    // The return value (how much of `delta_angle` was actually achieved) is
    // discarded here, unlike `fix_intersection_of_siblings`'s use of the
    // same primitive -- matches the vendored `fixIntersectionWithAncestor`
    // (`handleAncestorIntersections.inc:278`), which never captures it
    // either; `check_and_apply_config_changes` below is the real gate.
    std::vector<double> deltas;
    (void)calc_deltas(rotation_node, &ancestor, index_left, index_right, delta_angle, opts.paired,
                      opts.clearance, deltas);

    const IntersectionType it_log = is_exterior(ancestor) ? IntersectionType::exterior : it;
    changed = check_and_apply_config_changes(rotation_node, deltas, it_log, opts.unpaired,
                                             opts.paired, state);
  }

  return changed ? &rotation_node : nullptr;
}

TreeNode* handle_intersection_with_ancestor(TreeNode& ancestor, TreeNode& intersector,
                                            const PuzzlerOptions& opts, ResolverState& state) {
  // Ported from `handleIntersectionWithAncestor`
  // (`handleAncestorIntersections.inc:361`).
  const NodeIntersection ni = intersect_node_node(ancestor, intersector, opts.clearance);
  const IntersectionType it = ni.type;
  if (it == IntersectionType::none) {
    return nullptr;
  }

  std::vector<TreeNode*> path = construct_reduced_intersection_path(ancestor, intersector, it);
  const int path_length = static_cast<int>(path.size());

  std::vector<int> child_index(static_cast<std::size_t>(path_length - 1));
  for (int i = 0; i < path_length - 1; ++i) {
    child_index[static_cast<std::size_t>(i)] =
        get_child_index(*path[static_cast<std::size_t>(i)], path[static_cast<std::size_t>(i + 1)]->id);
  }

  TreeNode* changed_node = nullptr;
  const short rotation_sign = get_rotation_sign(path);

  if (rotation_sign != 0) {
    // Run from intersector to ancestor twice: first pass only considers
    // interior loops as rotation candidates.
    int node_number = path_length - 2;  // skip the intersector itself
    while (changed_node == nullptr && node_number >= 0) {
      TreeNode& candidate = *path[static_cast<std::size_t>(node_number)];
      if (is_interior_loop(candidate)) {
        changed_node =
            fix_intersection_with_ancestor(ancestor, candidate, intersector,
                                           child_index[static_cast<std::size_t>(node_number)],
                                           rotation_sign, it, opts, state);
      }
      --node_number;
    }

    // Second pass only considers multi loops as rotation candidates.
    node_number = path_length - 2;
    while (changed_node == nullptr && node_number >= 0) {
      TreeNode& candidate = *path[static_cast<std::size_t>(node_number)];
      if (is_multi_loop(candidate)) {
        changed_node =
            fix_intersection_with_ancestor(ancestor, candidate, intersector,
                                           child_index[static_cast<std::size_t>(node_number)],
                                           rotation_sign, it, opts, state);
      }
      --node_number;
    }
  }

  return changed_node;
}

TreeNode* check_node_against_ancestors(TreeNode& node, const PuzzlerOptions& opts,
                                       ResolverState& state) {
  // Ported from `checkNodeAgainstAncestors`
  // (`handleAncestorIntersections.inc:947`).
  TreeNode* changed_node = nullptr;
  TreeNode* ancestor = node.parent;
  TreeNode* top_level_ancestor = &node;

  // Move towards the root, checking and fixing ancestor intersections.
  while (!is_exterior(*ancestor)) {
    top_level_ancestor = ancestor;
    const NodeIntersection ni = intersect_node_node(node, *ancestor, opts.clearance);
    if (ni.type != IntersectionType::none) {
      changed_node = handle_intersection_with_ancestor(*ancestor, node, opts, state);
      if (changed_node != nullptr) {
        return changed_node;
      }
    }
    ancestor = ancestor->parent;
  }

  // Check and fix an ancestor intersection against the exterior.
  if (opts.check_exterior) {
    if (intersect_node_exterior(node, opts.check_exterior, opts.clearance)) {
      TreeNode* exterior = top_level_ancestor->parent;
      setup_exterior_bounding_boxes(*exterior, *top_level_ancestor, node, opts);
      changed_node = handle_intersection_with_ancestor(*exterior, node, opts, state);
    }
  }

  return changed_node;
}

}  // namespace rna_layout
