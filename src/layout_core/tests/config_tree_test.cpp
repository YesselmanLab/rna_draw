// ctest: `rna_layout::build_config_tree`/`update_bounding_boxes` sanity
// checks (`config_tree.hpp`). Parity against the vendored oracle (topology
// exact, config/box tight numeric) is the Python-level gate
// (`tests/test_native_parity.py::TestTreeParity`); this file only proves
// the port's own internal shape/invariants.

#include "rna_layout/config_tree.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "test_tree_helper.hpp"

using rna_layout::TreeNode;
using rna_layout::Vec2;

namespace {

int g_failures = 0;

void expect(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

void expect_near(double actual, double expected, double tol, const char* what) {
  if (std::fabs(actual - expected) > tol) {
    std::cerr << "FAIL " << what << ": expected " << expected << ", got " << actual << "\n";
    ++g_failures;
  }
}

void expect_finite(double value, const char* what) {
  if (!std::isfinite(value)) {
    std::cerr << "FAIL " << what << ": non-finite value " << value << "\n";
    ++g_failures;
  }
}

constexpr double kPaired = 35.0;
constexpr double kUnpaired = 25.0;

void test_hairpin_topology() {
  const auto tree = rna_layout::test::build_updated_tree("((((....))))", kPaired, kUnpaired);

  expect(rna_layout::is_exterior(*tree), "hairpin: root is exterior");
  expect(!tree->cfg.has_value(), "hairpin: root has no cfg");
  expect(tree->children.size() == 1, "hairpin: root has exactly one child (the outer stem)");

  const TreeNode& hairpin_loop = *tree->children[0];
  expect(!rna_layout::is_exterior(hairpin_loop), "hairpin: the loop node is not exterior");
  expect(hairpin_loop.cfg.has_value(), "hairpin: the loop node has a cfg");
  expect(hairpin_loop.children.empty(), "hairpin: the hairpin loop has no children");
  expect(hairpin_loop.lbox.has_value(), "hairpin: the loop node has a lbox after update");
  expect(hairpin_loop.sbox.has_value(), "hairpin: the loop node has a sbox after update");
}

void test_two_way_multiloop_topology() {
  // Outer stem opens a multiloop with two child stems -- root -> multiloop
  // node -> two hairpin nodes.
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...))))", kPaired, kUnpaired);

  expect(tree->children.size() == 1, "multiloop: root has one child (the outer stem)");
  const TreeNode& multiloop = *tree->children[0];
  expect(multiloop.children.size() == 2, "multiloop: the multiloop node has two children");
}

void test_multi_branch_exterior_topology() {
  // Three sibling stems directly off the exterior loop -- root has three
  // children, none nested inside each other.
  const auto tree =
      rna_layout::test::build_updated_tree("((....))..((....))..((....))", kPaired, kUnpaired);

  expect(tree->children.size() == 3, "multi-branch exterior: root has three direct children");
  for (const auto& child : tree->children) {
    expect(child->children.empty(), "multi-branch exterior: each hairpin child is a leaf");
  }
}

// NOLINTBEGIN(bugprone-unchecked-optional-access) -- each `->`/`.value()`
// below is preceded by (or, one line down from) an `expect(...has_value()..)`
// assertion; clang-tidy's flow analysis does not understand this file's
// plain `expect()` helper as a gate the way it understands `assert()`.
void test_updated_boxes_are_geometrically_sane() {
  const auto tree =
      rna_layout::test::build_updated_tree("((((..((((....))))..))))", kPaired, kUnpaired);

  const TreeNode& outer = *tree->children[0];
  const TreeNode& inner = *outer.children[0];

  for (const TreeNode* node : {&outer, &inner}) {
    expect(node->lbox.has_value() && node->sbox.has_value(), "internal loop: boxes are set");
    expect(node->lbox->radius > 0.0, "internal loop: loop radius is positive");

    // `a`/`b` are unit directions by construction (`update_bounding_boxes`
    // always writes a freshly rotated/normal-derived pair).
    const double len_a = std::hypot(node->sbox->a.x, node->sbox->a.y);
    const double len_b = std::hypot(node->sbox->b.x, node->sbox->b.y);
    expect_near(len_a, 1.0, 1e-9, "internal loop: stem 'a' direction is unit length");
    expect_near(len_b, 1.0, 1e-9, "internal loop: stem 'b' direction is unit length");

    expect(node->aabb.min.x <= node->aabb.max.x, "internal loop: aabb.min.x <= aabb.max.x");
    expect(node->aabb.min.y <= node->aabb.max.y, "internal loop: aabb.min.y <= aabb.max.y");

    expect_finite(node->lbox->center.x, "internal loop: lbox center x finite");
    expect_finite(node->lbox->center.y, "internal loop: lbox center y finite");
  }
}

void test_navigation_getters() {
  const auto tree = rna_layout::test::build_updated_tree("((((....))))", kPaired, kUnpaired);
  const TreeNode& loop = *tree->children[0];

  const Vec2 loop_center = rna_layout::get_loop_center(loop);
  const Vec2 stem_center = rna_layout::get_stem_center(loop);
  expect_near(loop_center.x, loop.lbox->center.x, 1e-12, "get_loop_center matches lbox.center.x");
  expect_near(stem_center.x, loop.sbox->c.x, 1e-12, "get_stem_center matches sbox.c.x");
}
// NOLINTEND(bugprone-unchecked-optional-access)

void test_get_child_angle_is_finite() {
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...))))", kPaired, kUnpaired);
  const TreeNode& multiloop = *tree->children[0];

  for (const auto& child : multiloop.children) {
    const double angle = rna_layout::get_child_angle(multiloop, *child);
    expect_finite(angle, "get_child_angle is finite");
  }
}

void test_bulge_folds_into_enclosing_stem() {
  // A single-unpaired-base bulge on the strand entering the loop has no
  // config tree node of its own (`config.cpp`'s bulge branch); the tree
  // should skip straight from the outer stem to the hairpin loop.
  const auto tree = rna_layout::test::build_updated_tree("(.((....)))", kPaired, kUnpaired);

  expect(tree->children.size() == 1, "bulge: root has one child");
  const TreeNode& hairpin = *tree->children[0];
  expect(hairpin.children.empty(), "bulge: the bulge is folded away -- one node, no children");
  expect(hairpin.sbox.has_value() && !hairpin.sbox->bulges.empty(),
         "bulge: the stem box records at least one bulge notch");
}

}  // namespace

int main() {
  test_hairpin_topology();
  test_two_way_multiloop_topology();
  test_multi_branch_exterior_topology();
  test_updated_boxes_are_geometrically_sane();
  test_navigation_getters();
  test_get_child_angle_is_finite();
  test_bulge_folds_into_enclosing_stem();

  if (g_failures > 0) {
    std::cerr << "config_tree_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
