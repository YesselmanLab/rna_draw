/**
 * @file intersect_tree.cpp
 * @brief Implements `intersect_tree.hpp`'s node-level predicates + the
 *        `detect_intersections` parity seam, ported from
 *        `intersectLevelTreeNodes.inc` + `intersectionType.inc`.
 */

#include "rna_layout/intersect_tree.hpp"

#include "rna_layout/config_tree.hpp"
#include "rna_layout/geometry.hpp"

namespace rna_layout {

namespace {

/// Ported from `intersectNodesBoundingBoxes` (`intersectLevelTreeNodes.inc:136-172`):
/// a cheap AABB-vs-AABB reject test `intersect_node_node` runs before any of
/// the more expensive box-level predicates.
bool intersect_nodes_bounding_boxes(const Aabb& aabb1, const Aabb& aabb2, const StemBox& stem1,
                                    const StemBox& stem2, double clearance) {
  double extra_distance = geom::epsilon_recognize(clearance);

  int count = 0;
  if (stem1.bulge_dist > 0.0) {
    ++count;
  }
  if (stem2.bulge_dist > 0.0) {
    ++count;
  }
  if (count > 0) {
    extra_distance += (1.0 / count) * (stem1.bulge_dist + stem2.bulge_dist);
  }

  if (aabb1.max.x < aabb2.min.x - extra_distance || aabb2.max.x < aabb1.min.x - extra_distance ||
      aabb1.max.y < aabb2.min.y - extra_distance || aabb2.max.y < aabb1.min.y - extra_distance) {
    return false;
  }
  return true;
}

/// DFS pre-order flatten, matching `debug_dump.cpp`'s `dump_tree_recursive`
/// numbering exactly (index == id) -- see `detect_intersections`'s doc
/// comment for why this must line up with the oracle's own dump.
void flatten_tree(const TreeNode& node, std::vector<const TreeNode*>& out) {
  out.push_back(&node);
  for (const auto& child : node.children) {
    flatten_tree(*child, out);
  }
}

}  // namespace

const char* intersection_type_to_string(IntersectionType type) {
  // Ported from `intersectionTypeToString` (`intersectionType.inc:27-56`).
  switch (type) {
    case IntersectionType::loop_loop:
      return "LxL";
    case IntersectionType::loop_stem:
      return "LxS";
    case IntersectionType::loop_bulge:
      return "LxB";
    case IntersectionType::stem_loop:
      return "SxL";
    case IntersectionType::stem_stem:
      return "SxS";
    case IntersectionType::stem_bulge:
      return "SxB";
    case IntersectionType::bulge_loop:
      return "BxL";
    case IntersectionType::bulge_stem:
      return "BxS";
    case IntersectionType::bulge_bulge:
      return "BxB";
    case IntersectionType::siblings:
      return "BRA";
    case IntersectionType::exterior:
      return "EXT";
    case IntersectionType::none:
    default:
      return "UNK";
  }
}

NodeIntersection intersect_node_node(const TreeNode& node1, const TreeNode& node2,
                                     double clearance) {
  // Ported verbatim from `intersectNodeNode` (`intersectLevelTreeNodes.inc:175-266`):
  // AABB reject, then SxS/LxL/SxL/LxS/LxB/BxL/SxB/BxS/BxB in that exact
  // order, each gated by the "is one the other's parent" relationship.
  if (&node1 == &node2) {
    return NodeIntersection{};
  }

  // NOLINTBEGIN(bugprone-unchecked-optional-access) -- see
  // `intersect_node_node`'s doc comment: callers never pass the true
  // exterior root here.
  const StemBox& sbox1 = *node1.sbox;
  const LoopBox& lbox1 = *node1.lbox;
  const StemBox& sbox2 = *node2.sbox;
  const LoopBox& lbox2 = *node2.lbox;
  // NOLINTEND(bugprone-unchecked-optional-access)

  const bool intersect =
      intersect_nodes_bounding_boxes(node1.aabb, node2.aabb, sbox1, sbox2, clearance);
  if (!intersect) {
    return NodeIntersection{};
  }

  const TreeNode* parent_of_node1 = node1.parent;
  const TreeNode* parent_of_node2 = node2.parent;
  const bool node1_is_parent_of_node2 = (&node1 == parent_of_node2);
  const bool node2_is_parent_of_node1 = (&node2 == parent_of_node1);
  const bool nodes_have_common_parent = (parent_of_node1 == parent_of_node2);

  if (!node1_is_parent_of_node2 && !node2_is_parent_of_node1 && !nodes_have_common_parent &&
      intersect_stem_stem(sbox1, sbox2)) {
    return NodeIntersection{IntersectionType::stem_stem, -1, -1};
  }

  if (!node1_is_parent_of_node2 && !node2_is_parent_of_node1 &&
      intersect_loop_loop(lbox1, lbox2, clearance)) {
    return NodeIntersection{IntersectionType::loop_loop, -1, -1};
  }

  if (!node2_is_parent_of_node1 && intersect_stem_loop(sbox1, lbox2, clearance)) {
    return NodeIntersection{IntersectionType::stem_loop, -1, -1};
  }

  if (!node1_is_parent_of_node2 && intersect_stem_loop(sbox2, lbox1, clearance)) {
    return NodeIntersection{IntersectionType::loop_stem, -1, -1};
  }

  if (!node1_is_parent_of_node2) {
    const BulgeHit hit = intersect_loop_bulges(lbox1, sbox2, clearance);
    if (hit.intersects) {
      return NodeIntersection{IntersectionType::loop_bulge, -1, hit.bulge};
    }
  }

  if (!node2_is_parent_of_node1) {
    const BulgeHit hit = intersect_loop_bulges(lbox2, sbox1, clearance);
    if (hit.intersects) {
      return NodeIntersection{IntersectionType::bulge_loop, hit.bulge, -1};
    }
  }

  {
    const BulgeHit hit = intersect_stem_bulges(sbox1, sbox2, clearance);
    if (hit.intersects) {
      return NodeIntersection{IntersectionType::stem_bulge, -1, hit.bulge};
    }
  }

  {
    const BulgeHit hit = intersect_stem_bulges(sbox2, sbox1, clearance);
    if (hit.intersects) {
      return NodeIntersection{IntersectionType::bulge_stem, hit.bulge, -1};
    }
  }

  {
    const BulgeBulgeHit hit = intersect_bulges_bulges(sbox1, sbox2, clearance);
    if (hit.intersects) {
      return NodeIntersection{IntersectionType::bulge_bulge, hit.bulge1, hit.bulge2};
    }
  }

  return NodeIntersection{};
}

bool intersect_node_exterior(const TreeNode& node, bool check_exterior_intersections,
                             double clearance) {
  if (is_exterior(node)) {
    return false;
  }
  if (is_exterior(*node.parent)) {
    return false;
  }

  // NOLINTBEGIN(bugprone-unchecked-optional-access) -- `node` is a direct
  // child of the exterior root (guaranteed by the two guards above), which
  // always has an `lbox` (see `config_tree.cpp`'s invariant note).
  const double cy = node.lbox->center.y;
  const double r = node.lbox->radius + geom::epsilon_recognize(clearance);
  // NOLINTEND(bugprone-unchecked-optional-access)

  if (check_exterior_intersections) {
    return (cy - r) <= geom::kExteriorY;
  }
  return false;
}

const TreeNode* intersect_node_tree(const TreeNode& node, const TreeNode& tree, double clearance) {
  const NodeIntersection it = intersect_node_node(node, tree, clearance);
  if (it.type != IntersectionType::none) {
    return &tree;
  }
  for (const auto& child : tree.children) {
    const TreeNode* found = intersect_node_tree(node, *child, clearance);
    if (found != nullptr) {
      return found;
    }
  }
  return nullptr;
}

namespace {

bool intersect_iterate_tree(const TreeNode& tree1, const TreeNode& tree2, double clearance) {
  if (intersect_node_tree(tree1, tree2, clearance) != nullptr) {
    return true;
  }
  for (const auto& child : tree1.children) {
    if (intersect_iterate_tree(*child, tree2, clearance)) {
      return true;
    }
  }
  return false;
}

}  // namespace

bool intersect_trees(const TreeNode& tree1, const TreeNode& tree2, double clearance) {
  return intersect_iterate_tree(tree1, tree2, clearance);
}

bool intersect_node_lists(const std::vector<const TreeNode*>& list1,
                          const std::vector<const TreeNode*>& list2,
                          bool check_exterior_intersections, double clearance) {
  for (const TreeNode* node1 : list1) {
    const bool is_exterior1 = is_exterior(*node1);
    for (const TreeNode* node2 : list2) {
      if (is_exterior1) {
        if (intersect_node_exterior(*node2, check_exterior_intersections, clearance)) {
          return true;
        }
      } else if (is_exterior(*node2)) {
        if (intersect_node_exterior(*node1, check_exterior_intersections, clearance)) {
          return true;
        }
      } else if (intersect_node_node(*node1, *node2, clearance).type != IntersectionType::none) {
        return true;
      }
    }
  }
  return false;
}

std::vector<Detection> detect_intersections(const TreeNode& root, double clearance) {
  std::vector<const TreeNode*> nodes;
  flatten_tree(root, nodes);
  const int n = static_cast<int>(nodes.size());

  std::vector<Detection> detections;
  for (int i = 1; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const NodeIntersection it = intersect_node_node(*nodes[i], *nodes[j], clearance);
      if (it.type != IntersectionType::none) {
        detections.push_back(Detection{i, j, it.type});
      }
    }
  }

  for (int i = 1; i < n; ++i) {
    if (nodes[i]->parent == &root &&
        intersect_node_exterior(*nodes[i], /*check_exterior_intersections=*/true, clearance)) {
      detections.push_back(Detection{i, 0, IntersectionType::exterior});
    }
  }

  return detections;
}

}  // namespace rna_layout
