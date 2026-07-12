// ctest: `rna_layout::check_siblings` (`../resolve_internal.hpp`, ported
// from `checkSiblings`, `handleSiblingIntersections.inc`) and the driver
// `rna_layout::check_and_fix_intersections` (`resolve.hpp`), Milestone A
// step 7 (SIBLING intersection resolution).

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "../resolve_internal.hpp"
#include "rna_layout/config_tree.hpp"
#include "rna_layout/intersect_tree.hpp"
#include "rna_layout/resolve.hpp"
#include "test_tree_helper.hpp"

using rna_layout::check_and_fix_intersections;
using rna_layout::check_siblings;
using rna_layout::IntersectionType;
using rna_layout::PuzzlerOptions;
using rna_layout::ResolverState;
using rna_layout::SiblingCheckResult;
using rna_layout::TreeNode;

namespace {

int g_failures = 0;

void expect(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

constexpr double kPaired = 35.0;
constexpr double kUnpaired = 25.0;

PuzzlerOptions sibling_only_options() {
  PuzzlerOptions opts;
  opts.paired = kPaired;
  opts.unpaired = kUnpaired;
  opts.check_sibling = true;
  opts.check_ancestor = false;
  opts.optimize = false;
  return opts;
}

/// Force `multiloop.children[1]`'s whole subtree to exactly overlap
/// `multiloop.children[0]`'s (same loop center) -- a deliberately
/// artificial, geometry-inconsistent-with-`Config` state, but a
/// deterministic way to guarantee `intersect_trees` (and a wide-open
/// bounding wedge overlap) without depending on turtle-layout's actual
/// initial placement happening to collide.
void collide_children_zero_and_one(TreeNode& multiloop) {
  TreeNode& first = *multiloop.children[0];
  TreeNode& second = *multiloop.children[1];
  const rna_layout::Vec2 onto_first = rna_layout::get_loop_center(first);
  const rna_layout::Vec2 second_center = rna_layout::get_loop_center(second);
  const rna_layout::Vec2 shift{onto_first.x - second_center.x, onto_first.y - second_center.y};
  rna_layout::translate_bounding_boxes(second, shift);
}

void test_check_siblings_none_when_no_intersections() {
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...)))(((...))))", kPaired, kUnpaired);
  TreeNode& multiloop = *tree->children[0];
  expect(multiloop.children.size() == 3, "setup: three-way multiloop has three children");

  ResolverState state;
  const SiblingCheckResult result = check_siblings(multiloop, sibling_only_options(), state);
  expect(result == SiblingCheckResult::none,
         "check_siblings: a freshly built (non-colliding) multiloop has no sibling intersections");
  expect(state.changes_applied == 0, "check_siblings: no changes applied when nothing intersects");
  expect(state.trace.empty(), "check_siblings: no trace entries when nothing intersects");
}

void test_check_siblings_restarts_and_records_a_trace_on_collision() {
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...)))(((...))))", kPaired, kUnpaired);
  TreeNode& multiloop = *tree->children[0];

  collide_children_zero_and_one(multiloop);
  const bool collided_before = rna_layout::intersect_trees(
      *multiloop.children[0], *multiloop.children[1], /*clearance=*/1.0);
  expect(collided_before, "setup: children 0 and 1 collide before check_siblings runs");

  ResolverState state;
  const SiblingCheckResult result = check_siblings(multiloop, sibling_only_options(), state);
  expect(result == SiblingCheckResult::restart,
         "check_siblings: an intersecting pair (fixed or not) always requests a restart");
  expect(state.changes_applied >= 1,
         "check_siblings: at least one check_and_apply_config_changes call was made");
  expect(!state.trace.empty(),
         "check_siblings: a colliding pair produces at least one trace entry");
  for (const auto& entry : state.trace) {
    expect(entry.type == IntersectionType::siblings,
           "check_siblings: every trace entry's type is siblings ('BRA')");
    expect(entry.node_id == multiloop.id,
           "check_siblings: every trace entry targets the common parent");
  }
}

void test_max_config_changes_budget_aborts() {
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...)))(((...))))", kPaired, kUnpaired);
  TreeNode& multiloop = *tree->children[0];
  collide_children_zero_and_one(multiloop);

  ResolverState state;
  state.max_config_changes = -1;  // already exceeded before the first call
  const SiblingCheckResult result = check_siblings(multiloop, sibling_only_options(), state);
  expect(result == SiblingCheckResult::aborted,
         "check_siblings: an already-exhausted change budget aborts immediately");
  expect(state.trace.empty(), "check_siblings: an aborted call makes no config changes");
}

// `check_ancestor == true` is exercised end-to-end in
// `resolve_ancestors_test.cpp` (Milestone A step 8) -- it no longer throws.

// `optimize == true` is exercised end-to-end in `optimize_test.cpp`
// (Milestone A step 9) -- it no longer throws.

void test_driver_resolves_a_colliding_multiloop_end_to_end() {
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...)))(((...))))", kPaired, kUnpaired);
  TreeNode& multiloop = *tree->children[0];
  collide_children_zero_and_one(multiloop);

  ResolverState state;
  TreeNode* result = check_and_fix_intersections(tree.get(), sibling_only_options(), state);
  expect(result == nullptr,
         "check_and_fix_intersections: returns nullptr (Milestone A step 7 scope)");
  expect(state.changes_applied >= 1,
         "check_and_fix_intersections: resolved the artificial collision");
}

void test_driver_is_a_no_op_on_an_already_clean_tree() {
  const auto tree = rna_layout::test::build_updated_tree("((((...))))", kPaired, kUnpaired);
  ResolverState state;
  TreeNode* result = check_and_fix_intersections(tree.get(), sibling_only_options(), state);
  expect(result == nullptr, "check_and_fix_intersections: returns nullptr on a clean tree");
  expect(state.changes_applied == 0,
         "check_and_fix_intersections: no changes on an already-clean tree");
}

}  // namespace

int main() {
  test_check_siblings_none_when_no_intersections();
  test_check_siblings_restarts_and_records_a_trace_on_collision();
  test_max_config_changes_budget_aborts();
  test_driver_resolves_a_colliding_multiloop_end_to_end();
  test_driver_is_a_no_op_on_an_already_clean_tree();

  if (g_failures > 0) {
    std::cerr << "resolve_siblings_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
