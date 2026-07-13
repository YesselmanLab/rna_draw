// ctest: `rna_layout::overlap::check_overlaps_batch` (`overlap_batch.hpp`).
// Covers the CPP-REVIEWER WARNING fix: mismatched outer-list lengths must
// raise `std::invalid_argument` up front rather than reading past the end
// of a shorter list inside a worker thread (UB, not catchable by the
// per-element try/catch). Plain assert-based `main()`, matching
// `overlap_check_test.cpp`'s style.

#include "rna_layout/overlap_batch.hpp"

#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace rna_layout::overlap;  // NOLINT(google-build-using-namespace)

namespace {

int g_failures = 0;

void expect_true(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

void expect_throws_invalid_argument(const std::function<void()>& fn, const char* what) {
  try {
    fn();
  } catch (const std::invalid_argument&) {
    return;
  } catch (...) {
    std::cerr << "FAIL " << what << ": threw the wrong exception type\n";
    ++g_failures;
    return;
  }
  std::cerr << "FAIL " << what << ": did not throw\n";
  ++g_failures;
}

// -----------------------------------------------------------------------
// Well-formed batch: unaffected by the new guard.
// -----------------------------------------------------------------------
void test_well_formed_batch_returns_one_count_per_structure() {
  const std::vector<std::vector<double>> xs = {{0.0, 30.0}, {0.0, 60.0}};
  const std::vector<std::vector<double>> ys = {{0.0, 0.0}, {0.0, 0.0}};
  const std::vector<std::vector<int>> pair_maps = {{-1, -1}, {-1, -1}};
  const std::vector<double> node_rs = {10.0, 10.0};
  const std::vector<int> counts = check_overlaps_batch(xs, ys, pair_maps, node_rs, 0.75, 1e-6, 1);
  expect_true(counts.size() == 2, "well-formed batch: one count per structure");
}

// -----------------------------------------------------------------------
// CPP-REVIEWER WARNING: mismatched outer-list lengths must raise, not
// crash. Each case below shortens exactly one of ys/pair_maps/node_rs
// relative to xs.
// -----------------------------------------------------------------------
void test_mismatched_ys_length_raises() {
  const std::vector<std::vector<double>> xs = {{0.0, 30.0}, {0.0, 60.0}};
  const std::vector<std::vector<double>> ys = {{0.0, 0.0}};
  const std::vector<std::vector<int>> pair_maps = {{-1, -1}, {-1, -1}};
  const std::vector<double> node_rs = {10.0, 10.0};
  expect_throws_invalid_argument(
      [&]() { (void)check_overlaps_batch(xs, ys, pair_maps, node_rs, 0.75, 1e-6, 1); },
      "shorter ys than xs raises invalid_argument");
}

void test_mismatched_pair_maps_length_raises() {
  const std::vector<std::vector<double>> xs = {{0.0, 30.0}, {0.0, 60.0}};
  const std::vector<std::vector<double>> ys = {{0.0, 0.0}, {0.0, 0.0}};
  const std::vector<std::vector<int>> pair_maps = {{-1, -1}};
  const std::vector<double> node_rs = {10.0, 10.0};
  expect_throws_invalid_argument(
      [&]() { (void)check_overlaps_batch(xs, ys, pair_maps, node_rs, 0.75, 1e-6, 1); },
      "shorter pair_maps than xs raises invalid_argument");
}

void test_mismatched_node_rs_length_raises() {
  const std::vector<std::vector<double>> xs = {{0.0, 30.0}, {0.0, 60.0}};
  const std::vector<std::vector<double>> ys = {{0.0, 0.0}, {0.0, 0.0}};
  const std::vector<std::vector<int>> pair_maps = {{-1, -1}, {-1, -1}};
  const std::vector<double> node_rs = {10.0};
  expect_throws_invalid_argument(
      [&]() { (void)check_overlaps_batch(xs, ys, pair_maps, node_rs, 0.75, 1e-6, 1); },
      "shorter node_rs than xs raises invalid_argument");
}

/// Same mismatch as `test_mismatched_ys_length_raises`, but forced onto the
/// multi-threaded path (the fix must guard before any worker is spawned,
/// not just the `thread_count <= 1` shortcut).
void test_mismatched_length_raises_on_multithreaded_path() {
  const std::vector<std::vector<double>> xs = {{0.0, 30.0}, {0.0, 60.0}, {0.0, 90.0}, {0.0, 120.0}};
  const std::vector<std::vector<double>> ys = {{0.0, 0.0}, {0.0, 0.0}};
  const std::vector<std::vector<int>> pair_maps = {{-1, -1}, {-1, -1}, {-1, -1}, {-1, -1}};
  const std::vector<double> node_rs = {10.0, 10.0, 10.0, 10.0};
  expect_throws_invalid_argument(
      [&]() { (void)check_overlaps_batch(xs, ys, pair_maps, node_rs, 0.75, 1e-6, 4); },
      "shorter ys than xs raises invalid_argument on the multi-threaded path");
}

}  // namespace

int main() {
  test_well_formed_batch_returns_one_count_per_structure();
  test_mismatched_ys_length_raises();
  test_mismatched_pair_maps_length_raises();
  test_mismatched_node_rs_length_raises();
  test_mismatched_length_raises_on_multithreaded_path();

  if (g_failures > 0) {
    std::cerr << g_failures << " failure(s)\n";
    return 1;
  }
  std::cout << "overlap_batch_test: all checks passed\n";
  return 0;
}
