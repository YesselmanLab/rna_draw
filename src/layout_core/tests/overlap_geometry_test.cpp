// ctest: predicate parity for `rna_layout::overlap` geometry
// (`overlap_geometry.hpp`) against hardcoded outputs computed once from the
// frozen Python arbiter (`rna_draw/geometry.py`) -- see that script's
// invocation in this port's plan (`.claude/plans/current-plan-checker.md`,
// step 1). Plain assert-based `main()`, matching `geometry_test.cpp`'s style.
//
// Tolerance: 1e-9, not the tighter 1e-12 a pure-arithmetic check could use --
// several of these values round-trip through `hypot`, and `std::hypot`
// (not guaranteed correctly-rounded) vs Python's `math.hypot` (correctly
// rounded, Py3.8+) can differ by ~1 ULP (~1e-13 at these magnitudes). 1e-9
// stays far tighter than any real port bug would produce while not
// false-failing on that last-ULP libm noise (see overlap_geometry.hpp's
// file header and the plan's Risk R1).

#include "rna_layout/overlap_geometry.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>

#include "rna_layout/types.hpp"

using rna_layout::Vec2;
namespace ov = rna_layout::overlap;

namespace {

int g_failures = 0;

void expect_near(double actual, double expected, const char* what) {
  constexpr double kTol = 1e-9;
  if (std::fabs(actual - expected) > kTol) {
    std::cerr << "FAIL " << what << ": expected " << expected << ", got " << actual << "\n";
    ++g_failures;
  }
}

void expect_none(const std::optional<double>& actual, const char* what) {
  if (actual.has_value()) {
    std::cerr << "FAIL " << what << ": expected std::nullopt, got " << *actual << "\n";
    ++g_failures;
  }
}

void expect_some(const std::optional<double>& actual, double expected, const char* what) {
  if (!actual.has_value()) {
    std::cerr << "FAIL " << what << ": expected " << expected << ", got std::nullopt\n";
    ++g_failures;
    return;
  }
  expect_near(*actual, expected, what);
}

void test_clamp() {
  expect_near(ov::clamp(5.0, 0.0, 10.0), 5.0, "clamp(5,0,10) == 5");
  expect_near(ov::clamp(-5.0, 0.0, 10.0), 0.0, "clamp(-5,0,10) == 0");
  expect_near(ov::clamp(15.0, 0.0, 10.0), 10.0, "clamp(15,0,10) == 10");
}

void test_point_segment_distance() {
  expect_near(ov::point_segment_distance(Vec2{5, 5}, Vec2{0, 0}, Vec2{10, 0}), 5.0,
              "point_segment_distance interior projection");
  expect_near(ov::point_segment_distance(Vec2{-5, 5}, Vec2{0, 0}, Vec2{10, 0}), 7.0710678118654755,
              "point_segment_distance clamp t=0");
  expect_near(ov::point_segment_distance(Vec2{15, 5}, Vec2{0, 0}, Vec2{10, 0}), 7.0710678118654755,
              "point_segment_distance clamp t=1");
  expect_near(ov::point_segment_distance(Vec2{3, 4}, Vec2{1, 1}, Vec2{1, 1}), 3.605551275463989,
              "point_segment_distance degenerate segment");
}

void test_closest_params() {
  auto [s1, t1] = ov::closest_params(Vec2{0, 0}, Vec2{10, 0}, Vec2{5, -5}, Vec2{0, 10});
  expect_near(s1, 0.5, "closest_params crossing s");
  expect_near(t1, 0.5, "closest_params crossing t");

  // denom == 0.0: two non-degenerate, parallel segments.
  auto [s2, t2] = ov::closest_params(Vec2{0, 0}, Vec2{10, 0}, Vec2{0, 5}, Vec2{10, 0});
  expect_near(s2, 0.0, "closest_params parallel (denom==0) s");
  expect_near(t2, 0.0, "closest_params parallel (denom==0) t");

  // c_dot_c == 0.0 guard, exercised directly (segment 2's direction is the
  // zero vector -- segment_segment_distance filters this case out before
  // ever reaching closest_params, but the guard itself is unit-testable).
  auto [s3, t3] = ov::closest_params(Vec2{0, 0}, Vec2{10, 0}, Vec2{5, 5}, Vec2{0, 0});
  expect_near(s3, 0.0, "closest_params c_dot_c==0 guard s");
  expect_near(t3, 0.0, "closest_params c_dot_c==0 guard t");

  // a_dot_a == 0.0 guard, exercised directly (segment 1 degenerate).
  auto [s4, t4] = ov::closest_params(Vec2{5, 5}, Vec2{0, 0}, Vec2{0, 0}, Vec2{10, 0});
  expect_near(s4, 0.0, "closest_params a_dot_a==0 guard s");
  expect_near(t4, 0.5, "closest_params a_dot_a==0 guard t");
}

void test_segment_segment_distance() {
  expect_near(ov::segment_segment_distance(Vec2{0, 0}, Vec2{10, 10}, Vec2{0, 10}, Vec2{10, 0}), 0.0,
              "segment_segment_distance crossing segments");
  expect_near(ov::segment_segment_distance(Vec2{0, 0}, Vec2{10, 0}, Vec2{0, 5}, Vec2{10, 5}), 5.0,
              "segment_segment_distance parallel segments");
  expect_near(ov::segment_segment_distance(Vec2{0, 0}, Vec2{1, 0}, Vec2{5, 5}, Vec2{6, 6}),
              6.4031242374328485, "segment_segment_distance endpoint-to-endpoint");
  expect_near(ov::segment_segment_distance(Vec2{3, 3}, Vec2{3, 3}, Vec2{0, 0}, Vec2{10, 0}), 3.0,
              "segment_segment_distance degenerate seg1");
  expect_near(ov::segment_segment_distance(Vec2{0, 0}, Vec2{10, 0}, Vec2{3, 3}, Vec2{3, 3}), 3.0,
              "segment_segment_distance degenerate seg2");
  expect_near(ov::segment_segment_distance(Vec2{0, 0}, Vec2{0, 0}, Vec2{3, 4}, Vec2{3, 4}), 5.0,
              "segment_segment_distance both degenerate");
}

void test_disks_overlap() {
  expect_some(ov::disks_overlap(ov::DiskGeom{0, 0, 10}, ov::DiskGeom{15, 0, 10}, 1e-6), 5.0,
              "disks_overlap overlapping");
  expect_none(ov::disks_overlap(ov::DiskGeom{0, 0, 10}, ov::DiskGeom{25, 0, 10}, 1e-6),
              "disks_overlap far apart");
  expect_none(ov::disks_overlap(ov::DiskGeom{0, 0, 10}, ov::DiskGeom{20, 0, 10}, 1e-6),
              "disks_overlap exactly touching, not flagged");
}

void test_disk_capsule_overlap() {
  const ov::SegmentGeom capsule{0, -20, 0, 20, 5};
  expect_some(ov::disk_capsule_overlap(ov::DiskGeom{3, 0, 10}, capsule, 1e-6), 12.0,
              "disk_capsule_overlap overlapping");
  expect_none(ov::disk_capsule_overlap(ov::DiskGeom{20, 0, 10}, capsule, 1e-6),
              "disk_capsule_overlap far apart");
}

void test_capsule_capsule_overlap() {
  const ov::SegmentGeom a{0, 0, 10, 10, 2};
  const ov::SegmentGeom crossing{0, 10, 10, 0, 2};
  expect_some(ov::capsule_capsule_overlap(a, crossing, 1e-6), 4.0,
              "capsule_capsule_overlap crossing");
  const ov::SegmentGeom far{0, 50, 10, 60, 2};
  expect_none(ov::capsule_capsule_overlap(a, far, 1e-6), "capsule_capsule_overlap far apart");
}

}  // namespace

int main() {
  test_clamp();
  test_point_segment_distance();
  test_closest_params();
  test_segment_segment_distance();
  test_disks_overlap();
  test_disk_capsule_overlap();
  test_capsule_capsule_overlap();

  if (g_failures > 0) {
    std::cerr << g_failures << " failure(s)\n";
    return 1;
  }
  std::cout << "overlap_geometry_test: all checks passed\n";
  return 0;
}
