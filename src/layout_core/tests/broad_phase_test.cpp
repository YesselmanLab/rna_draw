// ctest: `any_intersection`'s (`intersect_tree.hpp`) SUPERSET property --
// SPEED lever A1 (`.claude/plans/current-plan-speed.md`). For many random
// trees, many (subtree, ancestor_list) node selections drawn from them, and
// a sweep of `clearance` values, `any_intersection` must return EXACTLY the
// same boolean as the brute-force formula it replaces
// (`check_optimize_intersections`, `optimize.cpp:54`):
//   intersect_node_lists(subtree, subtree, ...) ||
//   intersect_node_lists(subtree, ancestor_list, ...)
// This is the correctness argument's empirical half (the doc-comment
// argument in `broad_phase.cpp` is the proof; this is the check that the
// implementation actually matches it) -- exercises both the brute-force
// fallback (small `m`) and the real grid path (`m` above the threshold, via
// deliberately large generated structures) on the exact geometry real
// `check_optimize_intersections` calls see.

#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "rna_layout/config_tree.hpp"
#include "rna_layout/intersect_tree.hpp"
#include "rna_layout/tree.hpp"
#include "test_tree_helper.hpp"

using rna_layout::TreeNode;

namespace {

int g_failures = 0;

void expect(bool condition, const std::string& what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

/// A random, well-nested, empty-loop-free dot-bracket string of exactly
/// @p len characters -- deep/branchy enough that some of its nodes' subtrees
/// exceed `broad_phase.cpp`'s brute-force threshold, so the sweep below
/// actually exercises the grid path, not just its small-`m` fallback.
std::string gen_region(std::mt19937& rng, int len) {
  if (len <= 0) {
    return "";
  }
  if (len < 5) {
    return std::string(static_cast<std::size_t>(len), '.');
  }
  std::uniform_int_distribution<int> pick(0, 2);
  if (pick(rng) == 0) {
    const int lead = 1 + static_cast<int>(rng() % 3);
    return std::string(static_cast<std::size_t>(lead), '.') + gen_region(rng, len - lead);
  }
  const int inner_len = 1 + static_cast<int>(rng() % static_cast<unsigned>(len - 1));
  const std::string inner = gen_region(rng, inner_len);
  const int rest_len = len - 2 - inner_len;
  return "(" + inner + ")" + gen_region(rng, rest_len);
}

std::string gen_structure(std::mt19937& rng, int len) { return gen_region(rng, len); }

/// @p node itself, then every descendant, DFS pre-order -- a standalone
/// copy of `optimize2.cpp`'s private `collect_subtree_nodes` (not exposed
/// outside that TU), since this test only needs the public `TreeNode` shape.
void collect_subtree(const TreeNode& node, std::vector<const TreeNode*>& out) {
  out.push_back(&node);
  for (const auto& child : node.children) {
    collect_subtree(*child, out);
  }
}

/// @p node's ancestor chain, nearest first, up to and including the root --
/// a standalone copy of `optimize2.cpp`'s private `collect_ancestor_nodes`.
std::vector<const TreeNode*> collect_ancestors(const TreeNode& node) {
  std::vector<const TreeNode*> ancestors;
  const TreeNode* ancestor = node.parent;
  while (ancestor != nullptr) {
    ancestors.push_back(ancestor);
    ancestor = ancestor->parent;
  }
  return ancestors;
}

/// The brute-force formula `check_optimize_intersections` used before A1 --
/// the reference `any_intersection` must reproduce bit-for-bit.
bool brute_force_reference(const std::vector<const TreeNode*>& subtree,
                           const std::vector<const TreeNode*>& ancestor_list, double clearance) {
  return rna_layout::intersect_node_lists(subtree, subtree, /*check_exterior_intersections=*/true,
                                          clearance) ||
         rna_layout::intersect_node_lists(subtree, ancestor_list,
                                          /*check_exterior_intersections=*/true, clearance);
}

/// For every non-root node of @p root, compares `any_intersection` against
/// `brute_force_reference` at @p clearance; records any mismatch.
void check_tree_at_clearance(const TreeNode& root, double clearance, const std::string& label) {
  std::vector<const TreeNode*> all_nodes;
  collect_subtree(root, all_nodes);
  for (const TreeNode* node : all_nodes) {
    if (node->parent == nullptr) {
      continue;  // the root itself is never optimize_tree's `node` argument.
    }
    std::vector<const TreeNode*> subtree;
    collect_subtree(*node, subtree);
    const std::vector<const TreeNode*> ancestor_list = collect_ancestors(*node);

    const bool expected = brute_force_reference(subtree, ancestor_list, clearance);
    const bool actual = rna_layout::any_intersection(
        subtree, ancestor_list, /*check_exterior_intersections=*/true, clearance);
    expect(expected == actual, label + " node id=" + std::to_string(node->id) +
                                   " subtree=" + std::to_string(subtree.size()) +
                                   " ancestors=" + std::to_string(ancestor_list.size()) +
                                   " clearance=" + std::to_string(clearance));
  }
}

void test_random_trees_superset_property() {
  std::mt19937 rng(12345);
  constexpr double kClearanceSweep[] = {0.0, 0.5, 1.0, 2.0};
  constexpr int kNumTrees = 30;

  int trees_built = 0;
  for (int i = 0; i < kNumTrees; ++i) {
    const int len = 100 + static_cast<int>(rng() % 500);  // 100..599 nt
    const std::string structure = gen_structure(rng, len);
    if (structure.find('(') == std::string::npos) {
      continue;  // an all-unpaired draw has no tree to test.
    }
    const auto tree = rna_layout::test::build_updated_tree(structure, /*paired=*/35.0,
                                                           /*unpaired=*/25.0);
    ++trees_built;
    for (double clearance : kClearanceSweep) {
      check_tree_at_clearance(*tree, clearance, "tree#" + std::to_string(i));
    }
  }
  expect(trees_built > 0, "test setup: at least one random tree was actually built");
}

}  // namespace

int main() {
  test_random_trees_superset_property();

  if (g_failures > 0) {
    std::cerr << "broad_phase_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
