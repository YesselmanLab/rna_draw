/**
 * @file resolve_exterior.cpp
 * @brief Implements `resolve_exterior.hpp`.
 */

#include "rna_layout/resolve_exterior.hpp"

#include <vector>

#include "rna_layout/config_tree.hpp"
#include "rna_layout/intersect_tree.hpp"

namespace rna_layout {

namespace {

/// Ported from the "for each subtree: compute first base / backbone count"
/// block (`resolveExteriorChildIntersections.inc:117-141`).
struct SubtreeLayout {
  std::vector<int> first_base;
  std::vector<int> backbone;
};

SubtreeLayout compute_subtree_layout(const std::vector<int>& pair_table, int subtree_count) {
  SubtreeLayout layout;
  layout.first_base.assign(static_cast<std::size_t>(subtree_count), 0);
  layout.backbone.assign(static_cast<std::size_t>(subtree_count), 0);

  int subtree = 0;
  int base = 1;
  const int length = pair_table[0];
  while (base < length && subtree < subtree_count) {
    if (pair_table[base] > base) {
      layout.first_base[static_cast<std::size_t>(subtree)] = base;
      ++subtree;
      base = pair_table[base];
    } else {
      ++base;
      ++layout.backbone[static_cast<std::size_t>(subtree)];
    }
  }
  return layout;
}

/// Ported from the `while (changed)` collision-resolution loop
/// (`resolveExteriorChildIntersections.inc:170-245`): repeatedly translate
/// `child_tree_node[subtree]` right by `unpaired * backbone[subtree]` until
/// it no longer intersects any "upper" (and, if `allow_flipping`, "lower")
/// sibling already placed, then file it into whichever list it settled on.
/// A pure extraction (not a reordering) of the vendored inner loop -- kept
/// as its own function only because `resolve_exterior_children_intersection`
/// would otherwise exceed a reasonable single-function length, per the
/// plan's style notes on real-seam decomposition.
void resolve_one_subtree_collision(int subtree, const std::vector<TreeNode*>& child_tree_node,
                                   const std::vector<int>& backbone, bool allow_flipping,
                                   double clearance, double unpaired, std::vector<double>& distance,
                                   double& offset, std::vector<int>& upper,
                                   std::vector<int>& lower) {
  bool changed = true;
  while (changed) {
    changed = false;
    bool intersect_upper = false;
    bool intersect_lower = false;

    for (int upper_stem : upper) {
      intersect_upper =
          intersect_trees(*child_tree_node[static_cast<std::size_t>(subtree)],
                          *child_tree_node[static_cast<std::size_t>(upper_stem)], clearance);
      if (intersect_upper) {
        break;
      }
    }

    if (allow_flipping) {
      for (int lower_stem : lower) {
        intersect_lower =
            intersect_trees(*child_tree_node[static_cast<std::size_t>(subtree)],
                            *child_tree_node[static_cast<std::size_t>(lower_stem)], clearance);
        if (intersect_lower) {
          break;
        }
      }
    }

    if ((!allow_flipping && intersect_upper) ||
        (allow_flipping && intersect_upper && intersect_lower)) {
      distance[static_cast<std::size_t>(subtree)] += unpaired;
      const double fix_overlap =
          unpaired * static_cast<double>(backbone[static_cast<std::size_t>(subtree)]);
      translate_bounding_boxes(*child_tree_node[static_cast<std::size_t>(subtree)],
                               Vec2{fix_overlap, 0.0});
      offset += fix_overlap;
      changed = true;
    } else {
      if (allow_flipping && intersect_upper) {
        lower.push_back(subtree);
      } else {
        upper.push_back(subtree);
      }
    }
  }
}

/// Ported from the "for all subtrees: translate exterior bases between
/// previous and current subtree" tail of the main loop
/// (`resolveExteriorChildIntersections.inc:247-254`). Direct `[base]`
/// indexing (not `[base - 1]`) is INTENTIONAL -- see `resolve_exterior.hpp`.
void translate_gap_coordinates(const std::vector<int>& pair_table, int subtree,
                               const SubtreeLayout& layout, const std::vector<double>& distance,
                               double accumulated_translation, Coords& coords) {
  int current_base = 1;
  const int gap_start = pair_table[static_cast<std::size_t>(
      layout.first_base[static_cast<std::size_t>(subtree - 1)])];
  const int gap_end = layout.first_base[static_cast<std::size_t>(subtree)];
  for (int base = gap_start; base < gap_end; ++base, ++current_base) {
    coords.x[static_cast<std::size_t>(base)] +=
        current_base * distance[static_cast<std::size_t>(subtree)] + accumulated_translation;
  }
}

/// Ported from the "modify x- and y-coordinates for all subtrees" block
/// (`resolveExteriorChildIntersections.inc:262-283`), the loop AFTER the
/// main per-subtree loop in the vendored source (unchanged position here).
void translate_and_flip_subtree_coordinates(const std::vector<int>& pair_table, int subtree_count,
                                            const SubtreeLayout& layout,
                                            const std::vector<double>& distance,
                                            const std::vector<int>& lower, Coords& coords) {
  std::size_t current_lower = 0;
  double translation = 0.0;
  for (int subtree = 1; subtree < subtree_count; ++subtree) {
    translation += distance[static_cast<std::size_t>(subtree)] *
                   static_cast<double>(layout.backbone[static_cast<std::size_t>(subtree)]);
    const int base_start = layout.first_base[static_cast<std::size_t>(subtree)];
    const int base_end = pair_table[static_cast<std::size_t>(base_start)];
    for (int base = base_start; base < base_end; ++base) {
      coords.x[static_cast<std::size_t>(base)] += translation;
    }

    if (current_lower < lower.size() && subtree == lower[current_lower]) {
      const double exterior_y = coords.y[1];
      for (int base = base_start; base < base_end; ++base) {
        coords.y[static_cast<std::size_t>(base)] =
            2 * exterior_y - coords.y[static_cast<std::size_t>(base)];
      }
      ++current_lower;
    }
  }
}

}  // namespace

void resolve_exterior_children_intersection(TreeNode& exterior_root,
                                            const std::vector<int>& pair_table, double unpaired,
                                            bool allow_flipping, double clearance, Coords& coords) {
  const int subtree_count = static_cast<int>(exterior_root.children.size());
  if (subtree_count < 2) {
    return;
  }

  std::vector<TreeNode*> child_tree_node;
  child_tree_node.reserve(static_cast<std::size_t>(subtree_count));
  for (auto& child : exterior_root.children) {
    child_tree_node.push_back(child.get());
  }

  const SubtreeLayout layout = compute_subtree_layout(pair_table, subtree_count);
  std::vector<double> distance(static_cast<std::size_t>(subtree_count), 0.0);

  std::vector<int> upper{0};  // The first subtree starts on the upper side.
  std::vector<int> lower;

  double offset = 0.0;
  double accumulated_translation = 0.0;

  // Ported from the main per-subtree loop
  // (`resolveExteriorChildIntersections.inc:159-254`): translate by the
  // running offset, resolve this subtree's own collisions, THEN patch the
  // gap coordinates before it and roll `accumulated_translation` forward --
  // in that exact order, one subtree at a time (not split into separate
  // passes: the gap patch for subtree K reads `distance[K]`, which is only
  // final once subtree K's own collision resolution above it has run).
  for (int subtree = 1; subtree < subtree_count; ++subtree) {
    if (offset > 0.0) {
      translate_bounding_boxes(*child_tree_node[static_cast<std::size_t>(subtree)],
                               Vec2{offset, 0.0});
    }

    resolve_one_subtree_collision(subtree, child_tree_node, layout.backbone, allow_flipping,
                                  clearance, unpaired, distance, offset, upper, lower);

    translate_gap_coordinates(pair_table, subtree, layout, distance, accumulated_translation,
                              coords);
    accumulated_translation +=
        distance[static_cast<std::size_t>(subtree)] *
        static_cast<double>(layout.backbone[static_cast<std::size_t>(subtree)]);
  }

  // Last part of the exterior loop, after the final subtree.
  const int tail_start = pair_table[static_cast<std::size_t>(
      layout.first_base[static_cast<std::size_t>(subtree_count - 1)])];
  for (int base = tail_start; base < pair_table[0]; ++base) {
    coords.x[static_cast<std::size_t>(base)] += accumulated_translation;
  }

  translate_and_flip_subtree_coordinates(pair_table, subtree_count, layout, distance, lower,
                                         coords);
}

}  // namespace rna_layout
