// ctest: `rna_layout::optimize_node`/`optimize_tree`
// (`../resolve_internal.hpp`, ported from `optimizeNode`/`optimizeTree`,
// `optimize.inc`) and the driver's own optimize gate (`resolve.cpp`),
// Milestone A step 9 (the OPTIMIZE resolver pass).

#include <cstdlib>
#include <iostream>

#include "../resolve_internal.hpp"
#include "rna_layout/config_tree.hpp"
#include "rna_layout/intersect_tree.hpp"
#include "rna_layout/resolve.hpp"
#include "test_tree_helper.hpp"

using rna_layout::check_and_fix_intersections;
using rna_layout::Config;
using rna_layout::optimize_node;
using rna_layout::optimize_tree;
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

PuzzlerOptions optimize_only_options() {
  PuzzlerOptions opts;
  opts.paired = kPaired;
  opts.unpaired = kUnpaired;
  opts.check_sibling = false;
  opts.check_ancestor = false;
  opts.check_exterior = false;
  opts.optimize = true;
  return opts;
}

// A three-way multiloop with an empty `subtree`/`ancestor_list` context (no
// other node can ever intersect it in this fixture): the simplest setting
// in which `optimize_node`'s own radius-shrinking logic is directly
// observable in isolation from the rest of the resolver.
constexpr const char* kMultiloopStructure = "((((...)))(((...)))(((...))))";

void test_optimize_node_noop_on_hairpin() {
  // "((((...))))" is a single stacked stem leading straight to one hairpin
  // loop -- the config tree has exactly one non-root node (the hairpin
  // itself, `tree->children[0]`, with no children of its own); see
  // `resolve_ancestors_test.cpp`'s `test_driver_is_a_no_op_on_an_already_
  // clean_tree` for the same one-node-tree fixture.
  const auto tree = rna_layout::test::build_updated_tree("((((...))))", kPaired, kUnpaired);
  TreeNode& hairpin = *tree->children[0];
  expect(hairpin.children.empty(), "setup: hairpin has no children");

  ResolverState state;
  const std::vector<const TreeNode*> empty_list;
  const double ratio =
      optimize_node(hairpin, empty_list, empty_list, optimize_only_options(), state);
  expect(ratio == 1.0, "optimize_node: hairpin loops (no children) are a no-op (ratio 1.0)");
  expect(state.changes_applied == 0, "optimize_node: hairpin loops apply no changes");
}

void test_optimize_node_noop_when_radius_not_inflated() {
  // A freshly built loop's radius equals its own default_radius (untouched
  // by any resolver) -- `cfg.radius - cfg.default_radius < 5.0` short-
  // circuits before any shrink attempt.
  const auto tree = rna_layout::test::build_updated_tree(kMultiloopStructure, kPaired, kUnpaired);
  TreeNode& multiloop = *tree->children[0];
  const Config& cfg = *multiloop.cfg;  // NOLINT(bugprone-unchecked-optional-access)
  expect(cfg.radius == cfg.default_radius,
         "setup: a freshly built loop's radius equals its default_radius");

  ResolverState state;
  const std::vector<const TreeNode*> empty_list;
  const double ratio =
      optimize_node(multiloop, empty_list, empty_list, optimize_only_options(), state);
  expect(ratio == 1.0, "optimize_node: an un-inflated radius is a no-op (ratio 1.0)");
  expect(state.changes_applied == 0,
         "optimize_node: no changes applied when radius is not inflated");
}

void test_optimize_node_shrinks_an_artificially_inflated_radius() {
  // Simulate what the sibling/ancestor resolver leaves behind on a node it
  // widened: apply a much larger radius directly (`radius_new > 0`, the
  // same sentinel `check_and_apply_config_changes` uses) with no angle
  // deltas, then let `optimize_node` try to shrink it back down. With an
  // EMPTY subtree/ancestor_list (nothing else in the tree to collide
  // against), `shrink_loop_radius`'s intersection checks always pass, so
  // the loop shrinks all the way back to its own geometric minimum.
  const auto tree = rna_layout::test::build_updated_tree(kMultiloopStructure, kPaired, kUnpaired);
  TreeNode& multiloop = *tree->children[0];
  // A live reference (not a copy): stays valid and up to date across the
  // in-place mutation below (`apply_changes_to_config_and_bounding_boxes`
  // never reallocates `multiloop.cfg`'s contained `Config`).
  const Config& cfg = *multiloop.cfg;  // NOLINT(bugprone-unchecked-optional-access)
  const double default_radius = cfg.default_radius;
  const double inflated_radius = default_radius + 50.0;
  rna_layout::apply_changes_to_config_and_bounding_boxes(multiloop, {}, inflated_radius, kPaired,
                                                         kUnpaired);
  expect(cfg.radius > default_radius + 5.0,
         "setup: the loop's radius is now well above default_radius + 5.0");

  ResolverState state;
  const std::vector<const TreeNode*> empty_list;
  const double ratio =
      optimize_node(multiloop, empty_list, empty_list, optimize_only_options(), state);

  expect(ratio < 1.0, "optimize_node: shrinks an artificially inflated, uncontested radius");
  expect(cfg.radius < inflated_radius, "optimize_node: the applied radius actually decreased");
  expect(state.changes_applied == 1,
         "optimize_node: exactly one improvement was logged (numberOfChangesAppliedToConfig++)");
}

void test_optimize_tree_leaves_an_uninflated_tree_unchanged() {
  const auto tree = rna_layout::test::build_updated_tree(kMultiloopStructure, kPaired, kUnpaired);
  TreeNode& multiloop = *tree->children[0];
  const Config& cfg = *multiloop.cfg;  // NOLINT(bugprone-unchecked-optional-access)
  const double radius_before = cfg.radius;

  ResolverState state;
  const double ratio = optimize_tree(multiloop, optimize_only_options(), state);
  expect(ratio == 1.0, "optimize_tree: an already-minimal tree is a no-op (ratio 1.0)");
  expect(cfg.radius == radius_before,
         "optimize_tree: an already-minimal loop's radius is unchanged");
  expect(state.changes_applied == 0,
         "optimize_tree: no changes applied to an already-minimal tree");
}

void test_optimize_tree_returns_one_when_optimize_option_is_false() {
  const auto tree = rna_layout::test::build_updated_tree(kMultiloopStructure, kPaired, kUnpaired);
  TreeNode& multiloop = *tree->children[0];

  PuzzlerOptions opts = optimize_only_options();
  opts.optimize = false;
  ResolverState state;
  const double ratio = optimize_tree(multiloop, opts, state);
  expect(ratio == 1.0, "optimize_tree: opts.optimize == false is a no-op (vendored early return)");
}

void test_driver_runs_the_full_pipeline_and_produces_a_clean_tree() {
  const std::string structure = "((((...)))((....(((.....))).))(((...))))";
  const auto tree = rna_layout::test::build_updated_tree(structure, kPaired, kUnpaired);

  PuzzlerOptions opts;
  opts.paired = kPaired;
  opts.unpaired = kUnpaired;
  ResolverState state;
  TreeNode* result = check_and_fix_intersections(tree.get(), opts, state);
  expect(result == nullptr,
         "check_and_fix_intersections: the outermost call always returns nullptr");

  const auto remaining = rna_layout::detect_intersections(*tree, opts.clearance);
  expect(remaining.empty(),
         "check_and_fix_intersections (full pipeline): the resolved tree is intersection-free");
}

void test_driver_is_deterministic_with_optimize_on() {
  const std::string structure = "((((..((((....))))..))))(((...)))(((...)))(((...)))";
  const auto tree1 = rna_layout::test::build_updated_tree(structure, kPaired, kUnpaired);
  const auto tree2 = rna_layout::test::build_updated_tree(structure, kPaired, kUnpaired);

  PuzzlerOptions opts;
  opts.paired = kPaired;
  opts.unpaired = kUnpaired;
  ResolverState state1;
  ResolverState state2;
  check_and_fix_intersections(tree1.get(), opts, state1);
  check_and_fix_intersections(tree2.get(), opts, state2);

  const Config& cfg1 = *tree1->children[0]->cfg;  // NOLINT(bugprone-unchecked-optional-access)
  const Config& cfg2 = *tree2->children[0]->cfg;  // NOLINT(bugprone-unchecked-optional-access)
  expect(cfg1.radius == cfg2.radius,
         "check_and_fix_intersections (full pipeline): repeated runs converge to the identical "
         "radius");
  expect(state1.changes_applied == state2.changes_applied,
         "check_and_fix_intersections (full pipeline): repeated runs apply the identical number "
         "of changes");
}

}  // namespace

int main() {
  test_optimize_node_noop_on_hairpin();
  test_optimize_node_noop_when_radius_not_inflated();
  test_optimize_node_shrinks_an_artificially_inflated_radius();
  test_optimize_tree_leaves_an_uninflated_tree_unchanged();
  test_optimize_tree_returns_one_when_optimize_option_is_false();
  test_driver_runs_the_full_pipeline_and_produces_a_clean_tree();
  test_driver_is_deterministic_with_optimize_on();

  if (g_failures > 0) {
    std::cerr << "optimize_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
