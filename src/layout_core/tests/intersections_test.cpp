// ctest: `rna_layout::geom`'s line/segment/circle intersection primitives
// (`intersections.hpp`).
//
// `solve_square_equation`/`get_cut_points_of_circles`/
// `get_cut_points_of_circle_and_line` are NOT exercised by any parity dump
// yet (see `intersections.hpp`'s scope note: their only vendored callers
// are out-of-scope arc/PostScript code or dead code) -- per the plan's
// fidelity rule, these get TARGETED tests against HAND-COMPUTED reference
// values (traced through the exact ported expression tree, not just
// "geometrically plausible" checks), not just parity comparison.

#include "rna_layout/intersections.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "rna_layout/types.hpp"

using rna_layout::StemBox;
using rna_layout::Vec2;
namespace geom = rna_layout::geom;

namespace {

int g_failures = 0;

void expect(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

void expect_near(double actual, double expected, double tol, const char* what) {
  if (std::fabs(actual - expected) > tol) {
    std::cerr << "FAIL " << what << ": expected " << expected << ", got " << actual << "\n";
    ++g_failures;
  }
}

constexpr double kTol = 1e-9;

/*=============================================================================
 *  solve_square_equation -- hand-computed reference values
 *============================================================================*/

void test_solve_square_equation_two_roots() {
  // x^2 - 5x + 6 = 0 -> roots 2, 3. Reference computes solution1 via
  // `(-b + sqrt(discr)) / (2a)` (the LARGER root here) and solution2 via
  // `(-b - sqrt(discr)) / (2a)` (the smaller) -- traced by hand:
  // discr = 25 - 24 = 1, sqrt(discr) = 1, sol1 = (5+1)/2 = 3, sol2 = (5-1)/2 = 2.
  const geom::QuadraticSolutions sol = geom::solve_square_equation(1.0, -5.0, 6.0);
  expect(sol.count == 2, "solve_square_equation: two real roots");
  expect_near(sol.solution1, 3.0, kTol, "solve_square_equation: solution1 (the +sqrt root)");
  expect_near(sol.solution2, 2.0, kTol, "solve_square_equation: solution2 (the -sqrt root)");
}

void test_solve_square_equation_double_root() {
  // x^2 - 4x + 4 = 0 -> double root 2 (discr == 0 exactly).
  const geom::QuadraticSolutions sol = geom::solve_square_equation(1.0, -4.0, 4.0);
  expect(sol.count == 1, "solve_square_equation: double root reports count == 1");
  expect_near(sol.solution1, 2.0, kTol, "solve_square_equation: double root solution1");
  expect_near(sol.solution2, 2.0, kTol, "solve_square_equation: double root solution2");
}

void test_solve_square_equation_no_real_roots() {
  // x^2 + 1 = 0 -> discr = -4 < 0.
  const geom::QuadraticSolutions sol = geom::solve_square_equation(1.0, 0.0, 1.0);
  expect(sol.count == 0, "solve_square_equation: no real roots reports count == 0");
}

/*=============================================================================
 *  get_cut_points_of_circles -- hand-computed reference values
 *============================================================================*/

void test_get_cut_points_of_circles_classic_symmetric_pair() {
  // c1=(0,0) r=5, c2=(8,0) r=5: classic 3-4-5 construction -> cut points at
  // (4, 3) and (4, -3). Traced by hand through the "smallDY && !smallDX"
  // branch (dy == 0 here): k=0, l=-64, m=-16, p=-4, q=0, a=1, b=0, c=-9;
  // solve_square_equation(1, 0, -9) -> discr=36, sqrt=6, sol1=3, sol2=-3;
  // point1 = ((sol1*k+l)/m, sol1) = (4, 3), point2 = (4, -3).
  const geom::CutPoints cut =
      geom::get_cut_points_of_circles(Vec2{0.0, 0.0}, 5.0, Vec2{8.0, 0.0}, 5.0);
  expect(cut.count == 2, "get_cut_points_of_circles: two cut points");
  expect_near(cut.point1.x, 4.0, kTol, "get_cut_points_of_circles: point1.x");
  expect_near(cut.point1.y, 3.0, kTol, "get_cut_points_of_circles: point1.y");
  expect_near(cut.point2.x, 4.0, kTol, "get_cut_points_of_circles: point2.x");
  expect_near(cut.point2.y, -3.0, kTol, "get_cut_points_of_circles: point2.y");
}

void test_get_cut_points_of_circles_coincident() {
  const geom::CutPoints cut =
      geom::get_cut_points_of_circles(Vec2{1.0, 1.0}, 3.0, Vec2{1.0, 1.0}, 3.0);
  expect(cut.count == -1, "get_cut_points_of_circles: coincident circles report count == -1");
}

void test_get_cut_points_of_circles_concentric_different_radius() {
  const geom::CutPoints cut =
      geom::get_cut_points_of_circles(Vec2{2.0, 2.0}, 3.0, Vec2{2.0, 2.0}, 5.0);
  expect(cut.count == 0, "get_cut_points_of_circles: concentric circles report count == 0");
}

/*=============================================================================
 *  get_cut_points_of_circle_and_line -- hand-computed reference values
 *============================================================================*/

void test_get_cut_points_of_circle_and_line() {
  // Circle center=(0,0) r=5; horizontal line through y=0, anchored at
  // (-10, 0) with direction (1, 0) -> hits the circle at t=5 (x=-5) and
  // t=15 (x=5). Traced by hand: a=1, b=2*1*(-10)=-20, c=100-25=75;
  // discr=400-300=100, sqrt=10; sol1=(20+10)/2=15, sol2=(20-10)/2=5;
  // point1 = anchor + 15*dir = (5, 0), point2 = anchor + 5*dir = (-5, 0).
  const geom::CutPoints cut = geom::get_cut_points_of_circle_and_line(
      Vec2{0.0, 0.0}, 5.0, Vec2{-10.0, 0.0}, Vec2{1.0, 0.0});
  expect(cut.count == 2, "get_cut_points_of_circle_and_line: two cut points");
  expect_near(cut.point1.x, 5.0, kTol, "get_cut_points_of_circle_and_line: point1.x");
  expect_near(cut.point1.y, 0.0, kTol, "get_cut_points_of_circle_and_line: point1.y");
  expect_near(cut.point2.x, -5.0, kTol, "get_cut_points_of_circle_and_line: point2.x");
  expect_near(cut.point2.y, 0.0, kTol, "get_cut_points_of_circle_and_line: point2.y");
}

/*=============================================================================
 *  The remaining primitives: sanity checks (simple, hand-verifiable cases)
 *============================================================================*/

void test_match_line_point() {
  expect(geom::match_line_point(Vec2{0.0, 0.0}, Vec2{4.0, 4.0}, Vec2{2.0, 2.0}),
         "match_line_point: midpoint is on the segment");
  expect(!geom::match_line_point(Vec2{0.0, 0.0}, Vec2{4.0, 4.0}, Vec2{10.0, 10.0}),
         "match_line_point: point past the segment's end is not on it");
}

void test_intersect_circle_circle() {
  expect(geom::intersect_circle_circle(Vec2{0.0, 0.0}, 3.0, Vec2{5.0, 0.0}, 3.0),
         "intersect_circle_circle: overlapping circles");
  expect(!geom::intersect_circle_circle(Vec2{0.0, 0.0}, 1.0, Vec2{10.0, 0.0}, 1.0),
         "intersect_circle_circle: far-apart circles do not overlap");
}

void test_intersect_line_segments() {
  expect(
      geom::intersect_line_segments(Vec2{0.0, 0.0}, Vec2{4.0, 4.0}, Vec2{0.0, 4.0}, Vec2{4.0, 0.0}),
      "intersect_line_segments: crossing diagonals");
  expect(!geom::intersect_line_segments(Vec2{0.0, 0.0}, Vec2{1.0, 1.0}, Vec2{5.0, 5.0},
                                        Vec2{6.0, 6.0}),
         "intersect_line_segments: parallel, non-overlapping segments");
}

void test_project_point_onto_line() {
  // A=(0,0), B=(4,4) (45-degree line, non-axis-aligned -- avoids the
  // reference's own w.x==0 degeneracy for exactly horizontal/vertical
  // segments), p=(0,4). Traced by hand: u=(0,4), v=(4,4), w=(-4,4);
  // r = (4 - 0*4/-4) / (4 - 4*4/-4) = 4/8 = 0.5; point = (0,0)+0.5*(4,4) = (2,2).
  const Vec2 projected =
      geom::project_point_onto_line(Vec2{0.0, 0.0}, Vec2{4.0, 4.0}, Vec2{0.0, 4.0});
  expect_near(projected.x, 2.0, kTol, "project_point_onto_line: x");
  expect_near(projected.y, 2.0, kTol, "project_point_onto_line: y");
}

void test_closest_pt_point_obb() {
  // Axis-aligned stem rectangle: a=(1,0), b=(0,1), c=(0,0), e=(2,1) ->
  // spans x in [-2,2], y in [-1,1].
  StemBox stem;
  stem.a = Vec2{1.0, 0.0};
  stem.b = Vec2{0.0, 1.0};
  stem.c = Vec2{0.0, 0.0};
  stem.e = Vec2{2.0, 1.0};

  const Vec2 outside = geom::closest_pt_point_obb(stem, Vec2{5.0, 5.0});
  expect_near(outside.x, 2.0, kTol, "closest_pt_point_obb: clamped x for a far point");
  expect_near(outside.y, 1.0, kTol, "closest_pt_point_obb: clamped y for a far point");

  const Vec2 inside = geom::closest_pt_point_obb(stem, Vec2{1.0, 0.5});
  expect_near(inside.x, 1.0, kTol, "closest_pt_point_obb: unclamped x for an interior point");
  expect_near(inside.y, 0.5, kTol, "closest_pt_point_obb: unclamped y for an interior point");
}

void test_test_circle_triangle() {
  // Non-axis-aligned triangle (avoids the same w.x==0 degeneracy noted
  // above for `project_point_onto_line`/`closest_pt_point_bulge`'s
  // edge-projection fallback).
  const Vec2 a{0.0, 0.0};
  const Vec2 b{3.0, 0.5};
  const Vec2 c{1.0, 3.0};
  const Vec2 centroid{(a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0};

  expect(geom::test_circle_triangle(centroid, 0.001, a, b, c),
         "test_circle_triangle: a tiny circle centered inside the triangle intersects it");
  expect(geom::test_circle_triangle(a, 0.001, a, b, c),
         "test_circle_triangle: a tiny circle centered ON a vertex intersects it");
  expect(!geom::test_circle_triangle(Vec2{1000.0, 1000.0}, 1.0, a, b, c),
         "test_circle_triangle: a far-away circle does not intersect");
}

void test_closest_pt_point_bulge_interior_point_is_itself() {
  const Vec2 a{0.0, 0.0};
  const Vec2 b{3.0, 0.5};
  const Vec2 c{1.0, 3.0};
  const Vec2 centroid{(a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0};

  const Vec2 closest = geom::closest_pt_point_bulge(centroid, a, b, c);
  expect_near(closest.x, centroid.x, kTol,
              "closest_pt_point_bulge: an interior point's closest point is itself (x)");
  expect_near(closest.y, centroid.y, kTol,
              "closest_pt_point_bulge: an interior point's closest point is itself (y)");
}

}  // namespace

int main() {
  test_solve_square_equation_two_roots();
  test_solve_square_equation_double_root();
  test_solve_square_equation_no_real_roots();

  test_get_cut_points_of_circles_classic_symmetric_pair();
  test_get_cut_points_of_circles_coincident();
  test_get_cut_points_of_circles_concentric_different_radius();

  test_get_cut_points_of_circle_and_line();

  test_match_line_point();
  test_intersect_circle_circle();
  test_intersect_line_segments();
  test_project_point_onto_line();
  test_closest_pt_point_obb();
  test_test_circle_triangle();
  test_closest_pt_point_bulge_interior_point_is_itself();

  if (g_failures > 0) {
    std::cerr << "intersections_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
