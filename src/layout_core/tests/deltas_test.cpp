// ctest: `rna_layout::calc_deltas` (`../resolve_internal.hpp`, ported from
// `calcDeltas.inc`). Property-based (the redistribution algorithm has no
// simple closed form to golden-check against): every returned `deltas`
// vector either satisfies `cfg_is_valid` or is all-zero, always sums to
// (near) zero, and the function is deterministic.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "../resolve_internal.hpp"
#include "rna_layout/config.hpp"
#include "test_tree_helper.hpp"

using rna_layout::calc_deltas;
using rna_layout::cfg_is_valid;
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

/// Every `deltas` entry `calc_deltas` returns is EITHER all-zero (the
/// change was rejected/couldn't sum to zero) OR passes `cfg_is_valid` --
/// `calc_deltas`'s own final guard (`calcDeltas.inc:506-525`).
bool deltas_are_valid_or_all_zero(const TreeNode& node, const std::vector<double>& deltas) {
  bool all_zero = true;
  for (double delta : deltas) {
    if (delta != 0.0) {
      all_zero = false;
      break;
    }
  }
  return all_zero || cfg_is_valid(*node.cfg, deltas);
}

double sum(const std::vector<double>& deltas) {
  double total = 0.0;
  for (double delta : deltas) {
    total += delta;
  }
  return total;
}

void test_negative_delta_angle_returns_zero_and_untouched_deltas() {
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...)))(((...))))", kPaired, kUnpaired);
  const TreeNode& multiloop = *tree->children[0];

  std::vector<double> deltas(multiloop.children.size() + 1, -7.0);
  const double changed =
      calc_deltas(multiloop, multiloop.parent, /*index_left=*/0,
                  /*index_right=*/1, /*delta_angle=*/-0.1, kPaired, kClearance, deltas);
  expect(changed == 0.0, "calc_deltas: negative delta_angle returns 0.0");
  for (double delta : deltas) {
    expect(delta == -7.0, "calc_deltas: negative delta_angle leaves the out-param untouched");
  }
}

void test_three_way_multiloop_redistribution_is_sane() {
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...)))(((...))))", kPaired, kUnpaired);
  const TreeNode& multiloop = *tree->children[0];
  expect(multiloop.children.size() == 3, "setup: three-way multiloop has three children");

  std::vector<double> deltas;
  const double changed_angle =
      calc_deltas(multiloop, multiloop.parent, /*index_left=*/0,
                  /*index_right=*/1, /*delta_angle=*/0.05, kPaired, kClearance, deltas);

  expect_finite(changed_angle, "calc_deltas: changed_angle is finite");
  expect(changed_angle >= 0.0, "calc_deltas: changed_angle is non-negative");
  expect(changed_angle <= 0.05 + 1e-9, "calc_deltas: changed_angle does not exceed the request");
  expect(deltas.size() == multiloop.children.size() + 1,
         "calc_deltas: deltas sized children.size() + 1");
  for (double delta : deltas) {
    expect_finite(delta, "calc_deltas: every delta is finite");
  }
  expect(std::fabs(sum(deltas)) < 1e-6, "calc_deltas: deltas sum to ~0");
  expect(deltas_are_valid_or_all_zero(multiloop, deltas),
         "calc_deltas: deltas are cfg_is_valid or all-zero");
}

void test_is_deterministic() {
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...)))(((...))))", kPaired, kUnpaired);
  const TreeNode& multiloop = *tree->children[0];

  std::vector<double> first;
  std::vector<double> second;
  const double changed_first =
      calc_deltas(multiloop, multiloop.parent, 0, 1, 0.05, kPaired, kClearance, first);
  const double changed_second =
      calc_deltas(multiloop, multiloop.parent, 0, 1, 0.05, kPaired, kClearance, second);
  expect(changed_first == changed_second,
         "calc_deltas: repeated calls return the same changed_angle");
  expect(first == second, "calc_deltas: repeated calls return identical deltas");
}

void test_zero_delta_angle_is_a_no_op() {
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...)))(((...))))", kPaired, kUnpaired);
  const TreeNode& multiloop = *tree->children[0];

  std::vector<double> deltas;
  const double changed_angle = calc_deltas(multiloop, multiloop.parent, 0, 1, /*delta_angle=*/0.0,
                                           kPaired, kClearance, deltas);
  expect(changed_angle == 0.0, "calc_deltas: delta_angle == 0.0 accomplishes nothing");
  for (double delta : deltas) {
    expect(delta == 0.0, "calc_deltas: delta_angle == 0.0 leaves every delta at 0");
  }
}

}  // namespace

int main() {
  test_negative_delta_angle_returns_zero_and_untouched_deltas();
  test_three_way_multiloop_redistribution_is_sane();
  test_is_deterministic();
  test_zero_delta_angle_is_a_no_op();

  if (g_failures > 0) {
    std::cerr << "deltas_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
