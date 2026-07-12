// ctest: `rna_layout::check_node_against_ancestors` (`../resolve_internal.hpp`,
// ported from `checkNodeAgainstAncestors`, `handleAncestorIntersections.inc`)
// and the driver `rna_layout::check_and_fix_intersections` (`resolve.hpp`)
// with `check_ancestor` on, Milestone A step 8 (ANCESTOR intersection
// resolution).

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "../resolve_internal.hpp"
#include "rna_layout/config_tree.hpp"
#include "rna_layout/intersect_tree.hpp"
#include "rna_layout/resolve.hpp"
#include "test_tree_helper.hpp"

using rna_layout::check_and_fix_intersections;
using rna_layout::check_node_against_ancestors;
using rna_layout::IntersectionType;
using rna_layout::NodeIntersection;
using rna_layout::PuzzlerOptions;
using rna_layout::ResolverState;
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

PuzzlerOptions ancestor_only_options() {
  PuzzlerOptions opts;
  opts.paired = kPaired;
  opts.unpaired = kUnpaired;
  opts.check_sibling = false;
  opts.check_ancestor = true;
  opts.optimize = false;
  return opts;
}

// A three-way multiloop (branch index 1 has an ASYMMETRIC interior loop --
// 4 unpaired bases on one side, 1 on the other -- enclosing a deep hairpin)
// wrapped in one enclosing pair, giving a tree shaped exactly like
// `check_node_against_ancestors`'s intended use: `outer` (the multiloop-
// bearing node) is TWO tree levels above `deep_hairpin`, with a real
// (non-straight) interior loop, `branch1`, in between -- the minimum depth
// at which an LxL ancestor fix has an eligible rotation candidate (a DIRECT
// parent/child ancestor pair can never fix an LxL intersection: the direct
// ancestor itself is excluded from `construct_reduced_intersection_path`
// for Lx? types, and there is nothing else on the path to rotate -- see
// `resolve_ancestors.cpp`'s doc comments).
constexpr const char* kStructure = "((((...)))((....(((.....))).))(((...))))";

/// Translate @p node's whole subtree onto @p onto's loop center -- the same
/// "deliberately artificial, geometry-inconsistent-with-Config" collision
/// trick `resolve_siblings_test.cpp` uses, applied to an ancestor pair
/// instead of a sibling pair.
void collide_onto(TreeNode& node, const TreeNode& onto) {
  const rna_layout::Vec2 target = rna_layout::get_loop_center(onto);
  const rna_layout::Vec2 current = rna_layout::get_loop_center(node);
  const rna_layout::Vec2 shift{target.x - current.x, target.y - current.y};
  rna_layout::translate_bounding_boxes(node, shift);
}

void test_check_node_against_ancestors_none_on_a_clean_tree() {
  const auto tree = rna_layout::test::build_updated_tree(kStructure, kPaired, kUnpaired);
  TreeNode& outer = *tree->children[0];
  TreeNode& branch1 = *outer.children[1];
  TreeNode& deep_hairpin = *branch1.children[0];
  expect(branch1.children.size() == 1, "setup: branch1 is an interior loop (one child)");

  ResolverState state;
  TreeNode* result = check_node_against_ancestors(deep_hairpin, ancestor_only_options(), state);
  expect(result == nullptr,
         "check_node_against_ancestors: a freshly built (non-colliding) tree has no ancestor "
         "intersection");
  expect(state.changes_applied == 0,
         "check_node_against_ancestors: no changes applied when nothing intersects");
}

void test_check_node_against_ancestors_fixes_a_forced_collision() {
  const auto tree = rna_layout::test::build_updated_tree(kStructure, kPaired, kUnpaired);
  TreeNode& outer = *tree->children[0];
  TreeNode& branch1 = *outer.children[1];
  TreeNode& deep_hairpin = *branch1.children[0];

  collide_onto(deep_hairpin, outer);
  const NodeIntersection forced =
      rna_layout::intersect_node_node(deep_hairpin, outer, /*clearance=*/1.0);
  expect(forced.type != IntersectionType::none,
         "setup: deep_hairpin and outer collide before check_node_against_ancestors runs");

  ResolverState state;
  TreeNode* result = check_node_against_ancestors(deep_hairpin, ancestor_only_options(), state);
  expect(result == &branch1,
         "check_node_against_ancestors: the fix rotates branch1 (the only eligible node on the "
         "reduced path between deep_hairpin and outer)");
  expect(state.changes_applied >= 1,
         "check_node_against_ancestors: at least one check_and_apply_config_changes call was made");
  expect(!state.trace.empty(), "check_node_against_ancestors: a fixed collision leaves a trace");
  if (!state.trace.empty()) {
    expect(state.trace.front().node_id == branch1.id,
           "check_node_against_ancestors: the trace entry targets branch1");
    expect(state.trace.front().accepted,
           "check_node_against_ancestors: the recorded change was accepted");
  }
}

void test_driver_resolves_a_colliding_ancestor_pair_end_to_end() {
  const auto tree = rna_layout::test::build_updated_tree(kStructure, kPaired, kUnpaired);
  TreeNode& outer = *tree->children[0];
  TreeNode& branch1 = *outer.children[1];
  TreeNode& deep_hairpin = *branch1.children[0];
  collide_onto(deep_hairpin, outer);

  ResolverState state;
  TreeNode* result = check_and_fix_intersections(tree.get(), ancestor_only_options(), state);
  expect(result == nullptr,
         "check_and_fix_intersections: the outermost call always returns nullptr (see "
         "resolve.hpp)");
  expect(state.changes_applied >= 1,
         "check_and_fix_intersections: resolved the artificial ancestor collision");

  const NodeIntersection after =
      rna_layout::intersect_node_node(deep_hairpin, outer, /*clearance=*/1.0);
  expect(after.type == IntersectionType::none,
         "check_and_fix_intersections: the forced collision is actually gone afterwards");
}

void test_driver_is_a_no_op_on_an_already_clean_tree() {
  // A single hairpin, not `kStructure`: `kStructure`'s 3-way multiloop has a
  // GENUINE (unforced) exterior-baseline dip on one of its deeper branches
  // -- the ancestor resolver correctly fixes it (exercised implicitly by
  // `test_driver_resolves_a_colliding_ancestor_pair_end_to_end`, which
  // starts from `kStructure` too), so it is not a clean-tree fixture. A
  // trivial one-node tree, whose only node is a DIRECT child of the root,
  // is: `intersect_node_exterior`'s own guard (`is_exterior(*node.parent)`,
  // mirroring the vendored `isExterior(getParent(node))`) always excludes
  // direct root children from the exterior-baseline check entirely.
  const auto tree = rna_layout::test::build_updated_tree("((((...))))", kPaired, kUnpaired);
  ResolverState state;
  TreeNode* result = check_and_fix_intersections(tree.get(), ancestor_only_options(), state);
  expect(result == nullptr, "check_and_fix_intersections: returns nullptr on a clean tree");
  expect(state.changes_applied == 0,
         "check_and_fix_intersections: no changes on an already-clean tree");
}

void test_driver_combines_sibling_and_ancestor_checks() {
  const auto tree = rna_layout::test::build_updated_tree(kStructure, kPaired, kUnpaired);
  TreeNode& outer = *tree->children[0];
  TreeNode& branch1 = *outer.children[1];
  TreeNode& deep_hairpin = *branch1.children[0];
  collide_onto(deep_hairpin, outer);

  PuzzlerOptions opts = ancestor_only_options();
  opts.check_sibling = true;
  ResolverState state;
  TreeNode* result = check_and_fix_intersections(tree.get(), opts, state);
  expect(result == nullptr, "check_and_fix_intersections: sibling+ancestor still returns nullptr");
  expect(state.changes_applied >= 1,
         "check_and_fix_intersections: sibling+ancestor resolved the forced collision");
}

void test_driver_throws_when_optimize_requested() {
  const auto tree = rna_layout::test::build_updated_tree("((((...))))", kPaired, kUnpaired);
  PuzzlerOptions opts = ancestor_only_options();
  opts.optimize = true;
  ResolverState state;

  bool threw = false;
  try {
    check_and_fix_intersections(tree.get(), opts, state);
  } catch (const std::logic_error&) {
    threw = true;
  }
  expect(threw, "check_and_fix_intersections: optimize == true throws (Milestone A step 9)");
}

}  // namespace

int main() {
  test_check_node_against_ancestors_none_on_a_clean_tree();
  test_check_node_against_ancestors_fixes_a_forced_collision();
  test_driver_resolves_a_colliding_ancestor_pair_end_to_end();
  test_driver_is_a_no_op_on_an_already_clean_tree();
  test_driver_combines_sibling_and_ancestor_checks();
  test_driver_throws_when_optimize_requested();

  if (g_failures > 0) {
    std::cerr << "resolve_ancestors_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
