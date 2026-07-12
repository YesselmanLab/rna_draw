// ctest: `rna_layout::TreeNode` RAII sanity checks (`tree.hpp`). Not yet
// exercised by any layout pass (see that header's comment) -- this proves
// the ownership shape (unique_ptr children, non-owning parent) compiles
// and behaves before the config-tree step (Milestone A step 4) builds on
// it.

#include "rna_layout/tree.hpp"

#include <cstdlib>
#include <iostream>

using rna_layout::TreeNode;

namespace {

int g_failures = 0;

void expect(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

void test_add_child_wires_parent() {
  TreeNode root;
  TreeNode* child = root.add_child();
  expect(root.children.size() == 1, "add_child appends exactly one child");
  expect(child->parent == &root, "add_child wires the child's parent back to the root");
}

void test_multiple_children_and_grandchildren() {
  TreeNode root;
  TreeNode* first = root.add_child();
  TreeNode* second = root.add_child();
  TreeNode* grandchild = first->add_child();

  expect(root.children.size() == 2, "root has both children");
  expect(second->parent == &root, "second child's parent is the root");
  expect(first->children.size() == 1, "first child has one grandchild");
  expect(grandchild->parent == first, "grandchild's parent is the first child, not the root");
}

void test_fields_default_construct() {
  TreeNode node;
  expect(node.parent == nullptr, "a fresh TreeNode has no parent");
  expect(node.children.empty(), "a fresh TreeNode has no children");
  expect(!node.cfg.has_value(), "a fresh TreeNode has no config yet");
  expect(node.loop_start == -1, "a fresh TreeNode's loop_start defaults to -1");
  expect(node.stem_start == -1, "a fresh TreeNode's stem_start defaults to -1");
}

}  // namespace

int main() {
  test_add_child_wires_parent();
  test_multiple_children_and_grandchildren();
  test_fields_default_construct();
  // `root`'s destructor recursively frees every descendant via
  // unique_ptr, replacing the vendored `freeTree` -- there is nothing to
  // assert here directly, but running under `ctest` (no leak sanitizer
  // yet -- see plan step 4's asan note for the config-tree step) at least
  // proves it compiles and does not crash on teardown.

  if (g_failures > 0) {
    std::cerr << "tree_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
