// ctest: `rna_layout::overlap::SpatialHash` (`spatial_hash.hpp`), port of
// `rna_draw/spatial_hash.py:33-113`. Covers the superset property (every
// AABB-overlapping pair is a candidate), dedupe (Python uses a `set`),
// `insert_segment` tiling of a long diagonal, and the zero-length
// fallback. Plain assert-based `main()`, matching `geometry_test.cpp`'s
// style.

#include "rna_layout/spatial_hash.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

#include "rna_layout/types.hpp"

using rna_layout::Aabb;
using rna_layout::Vec2;
using rna_layout::overlap::SpatialHash;

namespace {

int g_failures = 0;

void expect_true(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

bool aabbs_overlap(const Aabb& a, const Aabb& b) {
  return a.min.x <= b.max.x && b.min.x <= a.max.x && a.min.y <= b.max.y && b.min.y <= a.max.y;
}

/// Every AABB-overlapping pair (brute-force ground truth) must appear in
/// `candidate_pairs()` -- the superset argument `spatial_hash.hpp`'s file
/// header states. Cell size is deliberately small relative to the AABBs so
/// the superset property is exercised across several cells, not just one.
void test_superset_property_over_random_boxes() {
  const std::vector<Aabb> boxes = {
      {Vec2{0, 0}, Vec2{5, 5}},   {Vec2{3, 3}, Vec2{8, 8}},    {Vec2{100, 100}, Vec2{105, 105}},
      {Vec2{-2, -2}, Vec2{2, 2}}, {Vec2{50, 0}, Vec2{60, 10}}, {Vec2{52, 2}, Vec2{58, 8}},
  };
  SpatialHash grid(3.0);
  for (std::size_t i = 0; i < boxes.size(); ++i) {
    grid.insert(static_cast<int>(i), boxes[i]);
  }
  const std::vector<std::pair<int, int>> candidates = grid.candidate_pairs();
  const std::set<std::pair<int, int>> candidate_set(candidates.begin(), candidates.end());

  for (std::size_t i = 0; i < boxes.size(); ++i) {
    for (std::size_t j = i + 1; j < boxes.size(); ++j) {
      if (!aabbs_overlap(boxes[i], boxes[j])) {
        continue;
      }
      const bool present = candidate_set.count({static_cast<int>(i), static_cast<int>(j)}) > 0;
      expect_true(present, "superset property: overlapping AABB pair present in candidate_pairs");
    }
  }
}

/// `candidate_pairs()` must never emit the same `(i, j)` twice, even when
/// two items share MANY cells (Python's `set`-based dedup, `spatial_hash.py:
/// 99-113` -- a double-counted pair would double-count a checker witness).
void test_candidate_pairs_deduped() {
  SpatialHash grid(1.0);
  // Two large, fully-overlapping boxes span many shared cells at this
  // cell size.
  grid.insert(0, Aabb{Vec2{0, 0}, Vec2{10, 10}});
  grid.insert(1, Aabb{Vec2{0, 0}, Vec2{10, 10}});
  const std::vector<std::pair<int, int>> candidates = grid.candidate_pairs();
  int count01 = 0;
  for (const auto& [i, j] : candidates) {
    if (i == 0 && j == 1) {
      ++count01;
    }
  }
  expect_true(count01 == 1, "candidate_pairs dedupes a pair sharing many cells");
}

/// `candidate_pairs()` returns `(i, j)` with `i < j` only.
void test_candidate_pairs_ascending() {
  SpatialHash grid(5.0);
  grid.insert(3, Aabb{Vec2{0, 0}, Vec2{1, 1}});
  grid.insert(1, Aabb{Vec2{0, 0}, Vec2{1, 1}});
  for (const auto& [i, j] : grid.candidate_pairs()) {
    expect_true(i < j, "candidate_pairs emits i < j");
  }
}

/// `insert_segment` tiles a long diagonal so a small, unrelated item near
/// the segment's midpoint still becomes a candidate -- the same superset
/// argument as `insert`, extended to segment pieces (`spatial_hash.py:
/// 60-97`'s docstring).
void test_insert_segment_tiling_long_diagonal() {
  SpatialHash grid(10.0);
  grid.insert_segment(/*item_index=*/0, Vec2{0, 0}, Vec2{1000, 1000}, /*pad=*/2.0);
  grid.insert(/*item_index=*/1, Aabb{Vec2{499, 499}, Vec2{501, 501}});  // sits on the midpoint
  const std::vector<std::pair<int, int>> candidates = grid.candidate_pairs();
  const bool present =
      std::find(candidates.begin(), candidates.end(), std::make_pair(0, 1)) != candidates.end();
  expect_true(present, "insert_segment tiling: midpoint-adjacent item is a candidate");
}

/// A zero-length segment (`x0==x1, y0==y1`) falls back to a single padded
/// point AABB (`spatial_hash.py:83-85`), not a division by a zero length.
void test_insert_segment_zero_length_fallback() {
  SpatialHash grid(5.0);
  grid.insert_segment(/*item_index=*/0, Vec2{10, 10}, Vec2{10, 10}, /*pad=*/3.0);
  grid.insert(/*item_index=*/1, Aabb{Vec2{9, 9}, Vec2{11, 11}});
  const std::vector<std::pair<int, int>> candidates = grid.candidate_pairs();
  const bool present =
      std::find(candidates.begin(), candidates.end(), std::make_pair(0, 1)) != candidates.end();
  expect_true(present, "insert_segment zero-length fallback still yields a candidate");
}

}  // namespace

int main() {
  test_superset_property_over_random_boxes();
  test_candidate_pairs_deduped();
  test_candidate_pairs_ascending();
  test_insert_segment_tiling_long_diagonal();
  test_insert_segment_zero_length_fallback();

  if (g_failures > 0) {
    std::cerr << g_failures << " failure(s)\n";
    return 1;
  }
  std::cout << "spatial_hash_test: all checks passed\n";
  return 0;
}
