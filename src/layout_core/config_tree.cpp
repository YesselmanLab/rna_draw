/**
 * @file config_tree.cpp
 * @brief Implements `config_tree.hpp`, ported from `configtree.inc`.
 */

#include "rna_layout/config_tree.hpp"

#include <cmath>

#include "rna_layout/bounding_boxes.hpp"
#include "rna_layout/geometry.hpp"

namespace rna_layout {

namespace {

// Mutually recursive: a loop's children are stems, each of which opens
// another loop.
void build_loop_children(TreeNode& node, int loop_start, int& node_id,
                         const std::vector<int>& pair_table, const std::vector<BaseInfo>& base_info,
                         const std::vector<Config>& configs);

/// Ported from `treeHandleStem` (`configtree.inc:807`): find the stem's
/// enclosed loop (skipping any bulge-folded loop-opening bases, which have
/// no `Config`/`loop_id` of their own -- `config.cpp`'s bulge branch),
/// append a new child node for it, and recurse into that loop's own
/// children.
void build_tree_stem(TreeNode& parent, int& node_id, int stem_start,
                     const std::vector<int>& pair_table, const std::vector<BaseInfo>& base_info,
                     const std::vector<Config>& configs) {
  ++node_id;  // Assigned to the child about to be created; not stored (see
              // `config_tree.hpp`'s DFS pre-order note in `build_config_tree`).
  int i = stem_start;
  while (!base_info[i].loop_id.has_value()) {
    ++i;
  }

  TreeNode* child = parent.add_child();
  child->loop_start = i;
  child->stem_start = stem_start;
  // Same invariant/tradeoff as `turtle.cpp`'s `configs[base_info[start]
  // .loop_id.value()]`: by construction the `while` loop above only stops
  // at a base whose `loop_id` is set.
  child->cfg = configs[base_info[i].loop_id.value()];  // NOLINT(bugprone-unchecked-optional-access)

  build_loop_children(*child, i, node_id, pair_table, base_info, configs);
}

void build_loop_children(TreeNode& node, int loop_start, int& node_id,
                         const std::vector<int>& pair_table, const std::vector<BaseInfo>& base_info,
                         const std::vector<Config>& configs) {
  const int end = pair_table[loop_start];
  for (int i = loop_start + 1; i < end; ++i) {
    if (pair_table[i] > i) {
      build_tree_stem(node, node_id, i, pair_table, base_info, configs);
      i = pair_table[i];
    }
  }
}

/// Attach each node's INITIAL (turtle-coordinate-derived) boxes. Ported
/// from `buildBoundingBoxes` (`configtree.inc:835`).
void attach_initial_bounding_boxes(TreeNode& node, const std::vector<int>& pair_table,
                                   const std::vector<BaseInfo>& base_info,
                                   const std::vector<Config>& configs, const Coords& coords,
                                   double bulge_dist) {
  const bool is_root = (node.parent == nullptr);
  if (!is_root) {
    node.lbox = build_loop_box(node.loop_start, pair_table, base_info, configs, coords);
    node.sbox = build_stem_box(node.stem_start, node.loop_start, pair_table, coords, bulge_dist);
    node.aabb = compute_aabb(*node.sbox, *node.lbox);
  }
  for (auto& child : node.children) {
    attach_initial_bounding_boxes(*child, pair_table, base_info, configs, coords, bulge_dist);
  }
}

}  // namespace

double stem_bulge_distance(double unpaired) {
  return std::sqrt(unpaired * unpaired - 0.25 * unpaired * unpaired);
}

std::unique_ptr<TreeNode> build_config_tree(const std::vector<int>& pair_table,
                                            const std::vector<BaseInfo>& base_info,
                                            const std::vector<Config>& configs,
                                            const Coords& coords, double bulge_dist) {
  auto root = std::make_unique<TreeNode>();
  root->loop_start = 1;
  root->stem_start = -1;
  // `root->cfg` stays `std::nullopt`: the exterior loop has no Config
  // (`createTreeNode(nodeID, NULL, 1, -1, pair_table, NULL)`,
  // `configtree.inc:875`).

  // The vendored `buildConfigtree` scans the exterior loop's direct stems
  // with `loopStart = 0` (NOT the root's own `loop_start = 1`), exploiting
  // `pair_table[0] == length` so the scan range covers the whole sequence
  // (`configtree.inc:880`, and see `createTreeNode`'s matching `cfg == NULL`
  // branch, `configtree.inc:703-706`) -- preserved here via the same `0`.
  int node_id = 0;
  build_loop_children(*root, /*loop_start=*/0, node_id, pair_table, base_info, configs);

  attach_initial_bounding_boxes(*root, pair_table, base_info, configs, coords, bulge_dist);
  return root;
}

// NOLINTBEGIN(bugprone-unchecked-optional-access) -- `cfg`/`lbox`/`sbox` are
// `std::optional` only because `std::nullopt` is how the ROOT (the exterior
// loop, which has neither) is represented; every non-root `TreeNode` is
// GUARANTEED to have all three set by the time this function runs on it
// (`build_config_tree`'s `attach_initial_bounding_boxes` sets `lbox`/`sbox`
// on every non-root node it visits, and `build_tree_stem` sets `cfg` on
// every node it creates -- there is no code path that creates a non-root
// node without them). Every access below is guarded by an `is_exterior`
// check (this function's own, or its caller's, since `build_config_tree`
// never calls this on a node without first establishing the invariant for
// that node's PARENT too) -- same "provable by construction, not by the
// type system" tradeoff `turtle.cpp` documents for its own `.value()` use.
void update_bounding_boxes(TreeNode& node, double paired, double unpaired) {
  if (!is_exterior(node)) {
    const long num_stem_backbones = std::lround((2.0 * node.sbox->e.x) / unpaired);
    const double stem_length = unpaired * static_cast<double>(num_stem_backbones);
    const double distance_stem_end_to_loop_center =
        std::sqrt(node.cfg->radius * node.cfg->radius - 0.25 * paired * paired);
    const double distance_stem_center_to_loop_center =
        0.5 * stem_length + distance_stem_end_to_loop_center;
    node.lbox->center = Vec2{
        node.sbox->c.x + distance_stem_center_to_loop_center * node.sbox->a.x,
        node.sbox->c.y + distance_stem_center_to_loop_center * node.sbox->a.y,
    };
    node.lbox->radius = node.cfg->radius;

    node.aabb = compute_aabb(*node.sbox, *node.lbox);
  }

  double child_angle_rad = 0.0;
  for (std::size_t i = 0; i < node.children.size(); ++i) {
    TreeNode& child = *node.children[i];
    StemBox& sbox = *child.sbox;
    // `lbox` here is the CHILD's own (not-yet-fixed-this-pass) loop box --
    // read BEFORE anything below writes to `sbox`, matching the vendored
    // read-then-write order (`configtree.inc:355-361`) exactly: for an
    // exterior child this reuses the child's already-built x as the new
    // stem's x, since the exterior "loop" is the straight EXTERIOR_Y line,
    // not a circle with one shared center.
    const LoopBox& lbox = *child.lbox;

    Vec2 parent_loop_center;
    if (is_exterior(node)) {
      parent_loop_center = Vec2{lbox.center.x, geom::kExteriorY};
    } else {
      parent_loop_center = get_loop_center(node);
    }

    const long num_stem_backbones = std::lround((2.0 * sbox.e.x) / unpaired);
    const double stem_length = unpaired * static_cast<double>(num_stem_backbones);
    sbox.e.x = 0.5 * stem_length;
    sbox.e.y = 0.5 * paired;

    if (is_exterior(node)) {
      child_angle_rad = geom::kPi;
    } else {
      child_angle_rad += node.cfg->arcs[i].angle;  // getArcAngle(node->cfg, i)
    }

    Vec2 a_fixed;
    if (is_exterior(node)) {
      a_fixed = Vec2{0.0, 1.0};
    } else {
      const double gamma = child_angle_rad - geom::kPi;
      a_fixed = geom::rotate_vector_by_angle(node.sbox->a, gamma);
    }
    sbox.a = a_fixed;

    Vec2 b_fixed = geom::normal(a_fixed);
    b_fixed.x *= -1;
    b_fixed.y *= -1;
    sbox.b = b_fixed;

    double s0 = 0.0;
    if (!is_exterior(node)) {
      s0 = std::sqrt(node.cfg->radius * node.cfg->radius - 0.25 * paired * paired);
    }
    const double distance_stem_center = s0 + 0.5 * stem_length;

    sbox.c = Vec2{
        parent_loop_center.x + distance_stem_center * a_fixed.x,
        parent_loop_center.y + distance_stem_center * a_fixed.y,
    };

    if (stem_length == 0) {
      sbox.e.x = geom::kEpsilon7;
    }
  }

  for (auto& child : node.children) {
    update_bounding_boxes(*child, paired, unpaired);
  }
}
// NOLINTEND(bugprone-unchecked-optional-access)

bool is_exterior(const TreeNode& node) { return node.parent == nullptr; }

// NOLINTBEGIN(bugprone-unchecked-optional-access) -- same invariant as
// `update_bounding_boxes`'s block comment above: `lbox`/`sbox` are set on
// every non-root node these getters are documented (`config_tree.hpp`) to
// require.
Vec2 get_loop_center(const TreeNode& node) { return node.lbox->center; }

Vec2 get_stem_center(const TreeNode& node) { return node.sbox->c; }

double get_child_angle(const TreeNode& parent, const TreeNode& child) {
  const Vec2 parent_loop_center = parent.lbox->center;
  const Vec2 parent_stem_center = parent.sbox->c;
  const Vec2 parent_loop_stem_vector = geom::vector_from_to(parent_loop_center, parent_stem_center);

  const Vec2 child_loop_center = child.lbox->center;
  double angle = geom::angle_pt_pt_pt(parent_stem_center, parent_loop_center, child_loop_center);

  if (!geom::is_to_the_right_point_vector(parent_loop_center, parent_loop_stem_vector,
                                          child_loop_center)) {
    angle = geom::kTwoPi - angle;
  }
  return angle;
}

void translate_bounding_boxes(TreeNode& node, Vec2 vector) {
  translate_stem_box(*node.sbox, vector);
  translate_loop_box(*node.lbox, vector);
  node.aabb = compute_aabb(*node.sbox, *node.lbox);

  for (auto& child : node.children) {
    translate_bounding_boxes(*child, vector);
  }
}
// NOLINTEND(bugprone-unchecked-optional-access)

}  // namespace rna_layout
