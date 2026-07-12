// ctest: analytic checks for `rna_layout::geom` primitives (`geometry.hpp`).
// Plain assert-based `main()`, matching
// `src/vienna_layout/tests/test_puzzler_smoke.cpp`'s style (no gtest
// dependency in this tree yet).

#include "rna_layout/geometry.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "rna_layout/types.hpp"

using rna_layout::Vec2;
namespace geom = rna_layout::geom;

namespace {

int g_failures = 0;

void expect_near(double actual, double expected, double tol, const char* what) {
  if (std::fabs(actual - expected) > tol) {
    std::cerr << "FAIL " << what << ": expected " << expected << ", got " << actual << "\n";
    ++g_failures;
  }
}

void test_length() {
  expect_near(geom::length(Vec2{3.0, 4.0}), 5.0, 1e-12, "length({3,4}) == 5");
  expect_near(geom::length(Vec2{0.0, 0.0}), 0.0, 1e-12, "length({0,0}) == 0");
}

void test_dot_and_cross() {
  expect_near(geom::dot(Vec2{1, 0}, Vec2{0, 1}), 0.0, 1e-12, "dot of perpendicular unit vectors");
  expect_near(geom::dot(Vec2{2, 3}, Vec2{4, 5}), 23.0, 1e-12, "dot({2,3},{4,5}) == 23");
  expect_near(geom::cross(Vec2{1, 0}, Vec2{0, 1}), 1.0, 1e-12, "cross of the standard basis");
}

void test_angle_between() {
  expect_near(geom::angle_between(Vec2{1, 0}, Vec2{0, 1}), geom::kPiHalf, 1e-9,
              "angle between perpendicular vectors == pi/2");
  expect_near(geom::angle_between(Vec2{1, 0}, Vec2{1, 0}), 0.0, 1e-9,
              "angle between identical vectors == 0");
  expect_near(geom::angle_between(Vec2{1, 0}, Vec2{-1, 0}), geom::kPi, 1e-9,
              "angle between opposite vectors == pi");
  // A non-axis-aligned pair: order-only divergences (e.g. accidentally
  // swapping the normalize-then-dot expression tree for dot-then-divide)
  // are invisible on axis-aligned inputs but show up here.
  expect_near(geom::angle_between(Vec2{2, 1}, Vec2{-1, 3}), 1.4288992721907325, 1e-9,
              "angle between non-axis-aligned vectors {2,1} vs {-1,3}");
}

void test_degree_radian_round_trip() {
  expect_near(geom::to_radian(geom::to_degree(1.2345)), 1.2345, 1e-12,
              "to_radian(to_degree(x)) == x");
  expect_near(geom::to_degree(geom::kPi), 180.0, 1e-9, "to_degree(pi) == 180");
}

void test_rotate_point_around_point() {
  // Rotating (1, 0) by +pi/2 (clockwise, per the vendored sign convention)
  // around the origin lands at (0, -1).
  const Vec2 rotated = geom::rotate_point_around_point(Vec2{1, 0}, Vec2{0, 0}, geom::kPiHalf);
  expect_near(rotated.x, 0.0, 1e-9, "rotate_point_around_point x");
  expect_near(rotated.y, -1.0, 1e-9, "rotate_point_around_point y");

  // A full 2*pi rotation returns to the start.
  const Vec2 full_turn =
      geom::rotate_point_around_point(Vec2{3.0, -2.0}, Vec2{1.0, 1.0}, geom::kTwoPi);
  expect_near(full_turn.x, 3.0, 1e-9, "full-turn rotation x");
  expect_near(full_turn.y, -2.0, 1e-9, "full-turn rotation y");
}

void test_distance_angle_round_trip() {
  const double radius = 10.0;
  const double distance = 6.0;
  const double angle = geom::distance_to_angle(radius, distance);
  expect_near(geom::angle_to_distance(radius, angle), distance, 1e-9,
              "angle_to_distance(distance_to_angle(d)) == d");
}

void test_epsilon_clearance_scaling() {
  expect_near(geom::epsilon_recognize(1.0), 14.0, 1e-12, "epsilon_recognize(1.0) == 14 (stock)");
  expect_near(geom::epsilon_fix(1.0), 19.0, 1e-12, "epsilon_fix(1.0) == 19 (stock)");
  expect_near(geom::epsilon_recognize(2.0), 28.0, 1e-12, "epsilon_recognize scales linearly");
}

}  // namespace

int main() {
  test_length();
  test_dot_and_cross();
  test_angle_between();
  test_degree_radian_round_trip();
  test_rotate_point_around_point();
  test_distance_angle_round_trip();
  test_epsilon_clearance_scaling();

  if (g_failures > 0) {
    std::cerr << "geometry_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
