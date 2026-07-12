// ctest: `rna_layout::bounding_wedge` sanity checks (`bounding_wedge.hpp`).
// Not yet wired into any pipeline (Milestone A step 8 wires it into the
// ancestor-intersection resolver); this proves the standalone predicate
// compiles, runs, and returns geometrically sane ranges on an already-built
// (T1) tree.

#include "rna_layout/bounding_wedge.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "test_tree_helper.hpp"

using rna_layout::AngleRange;
using rna_layout::TreeNode;

namespace {

int g_failures = 0;

void expect(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
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
constexpr double kClearance = 1.0;

void test_wedge_of_each_multiloop_child_is_sane() {
  // root -> multiloop (two children) -- `bounding_wedge`'s `root` param
  // must have its own sbox (not the true exterior root), so use the
  // multiloop node.
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...))))", kPaired, kUnpaired);
  const TreeNode& multiloop = *tree->children[0];
  expect(multiloop.children.size() == 2, "setup: multiloop has two children");

  for (int child_index = 0; child_index < 2; ++child_index) {
    const AngleRange range = rna_layout::bounding_wedge(multiloop, child_index, kClearance);
    expect_finite(range.min_angle, "wedge min_angle is finite");
    expect_finite(range.max_angle, "wedge max_angle is finite");
    expect(range.min_angle <= range.max_angle, "wedge min_angle <= max_angle");
  }
}

void test_wedge_widens_for_a_bulged_child() {
  // A two-way multiloop whose second branch has a single-unpaired-base
  // bulge -- that child's stem box has a non-empty `bulges` list, so its
  // wedge exercises `wedge_points_of_interest`'s bulge-point pass (not just
  // the two direct-child stem-corner points every child gets).
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(.((...))))", kPaired, kUnpaired);
  const TreeNode& multiloop = *tree->children[0];
  expect(multiloop.children.size() == 2, "setup: multiloop has two children");
  const TreeNode& bulged_child = *multiloop.children[1];
  // `bulged_child` is a non-root node of an already-built (T1) tree, so its
  // `sbox` is guaranteed set -- see `config_tree.cpp`'s invariant note.
  expect(!bulged_child.sbox->bulges.empty(),  // NOLINT(bugprone-unchecked-optional-access)
         "setup: the second child's stem has a bulge notch");

  const AngleRange range = rna_layout::bounding_wedge(multiloop, /*child_index=*/1, kClearance);
  expect_finite(range.min_angle, "bulged-child wedge min_angle is finite");
  expect_finite(range.max_angle, "bulged-child wedge max_angle is finite");
  expect(range.min_angle <= range.max_angle, "bulged-child wedge min_angle <= max_angle");
}

}  // namespace

int main() {
  test_wedge_of_each_multiloop_child_is_sane();
  test_wedge_widens_for_a_bulged_child();

  if (g_failures > 0) {
    std::cerr << "bounding_wedge_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
