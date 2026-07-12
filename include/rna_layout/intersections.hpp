#pragma once

/**
 * @file intersections.hpp
 * @brief Line/segment/circle/point-vs-shape intersection primitives, ported
 *        from `intersectLevelLines.inc` plus the generic point-projection
 *        helpers `intersectLevelBoundingBoxes.inc` defines ahead of its own
 *        stem/loop/bulge-specific predicates (`intersect_tree.hpp`).
 *
 * This is the "deferred from module 1" half of `vector_math.inc` +
 * `intersectLevelLines.inc`: `geometry.hpp`'s own docstring notes these
 * (circle/line intersection solvers, `solveSquareEquation`,
 * `getCutPointsOf*`) are ported here, at the detection-predicate step
 * (Milestone A step 5).
 *
 * FIDELITY NOTE (`.claude/plans/current-plan.md`'s "CRITICAL" section):
 * every function here preserves the reference's EXACT expression tree
 * (operation order), not just its math -- these are exactly the
 * order-sensitive primitives (circle/line cut points, `isToTheRight`) the
 * module-1 review flagged as the next trap after `angle_between`.
 *
 * SCOPE NOTE -- dead code NOT ported (verified via `grep` over the whole
 * vendored call graph reachable from `vrna_plot_coords_puzzler_pt`):
 * `intersectLineArc`/`intersectArcArc`/`matchPointArc` (`intersectLevelLines.inc`)
 * are called only from `RNApuzzler.c`'s `checkRemainingIntersections`, which
 * is itself only reachable when `arc_coords != NULL` -- rna_draw always
 * calls with `arc_coords = NULL` (`bindings.cpp:146`), so this whole
 * arc/PostScript path is OUT OF SCOPE per the plan's reference map.
 * `TestLoopBulge`/`ClosestPtPointTriangle`/`checkBounds`
 * (`intersectLevelBoundingBoxes.inc`, `intersectLevelTreeNodes.inc`) are
 * DEFINED but never CALLED anywhere in the vendored sources -- genuine dead
 * code, not ported.
 */

#include "rna_layout/types.hpp"

namespace rna_layout::geom {

/// The number of solutions `solve_square_equation` found (0, 1, or 2) plus
/// their values (only as many of `solution1`/`solution2` are meaningful as
/// `count` indicates). Mirrors `solveSquareEquation`'s `short` return +
/// `double*` out-params (`vector_math.inc:577`).
struct QuadraticSolutions {
  int count = 0;
  double solution1 = 0.0;
  double solution2 = 0.0;
};

/**
 * Solve `a*x^2 + b*x + c = 0` for real roots. Mirrors `solveSquareEquation`
 * (`vector_math.inc:577`) expression-for-expression (discriminant, then the
 * `(-b +- sqrt(discr)) / (2*a)` pair computed as two SEPARATE expressions,
 * not one shared `sqrt` reused with a sign flip).
 */
[[nodiscard]] QuadraticSolutions solve_square_equation(double a, double b, double c);

/// The number of cut points `get_cut_points_of_circles`/
/// `get_cut_points_of_circle_and_line` found, plus their coordinates (as
/// many of `point1`/`point2` as `count` indicates meaningful). `count == -1`
/// is `get_cut_points_of_circles`'s "circles coincide" case (infinite cut
/// points); `point1`/`point2` are then unset, matching the vendored
/// function's own contract (`ret1`/`ret2` left untouched by that branch).
struct CutPoints {
  int count = 0;
  Vec2 point1{};
  Vec2 point2{};
};

/**
 * The common point(s) of two circles. Mirrors `getCutPointsOfCircles`
 * (`vector_math.inc:607`) expression-for-expression, INCLUDING its
 * small-delta classification (`eps = 1.0`, a hardcoded literal distinct from
 * `EPSILON_7` -- preserved as written) and its two near-mirror-image
 * branches (`!smallDY` vs `smallDY && !smallDX`), each solving for the OTHER
 * axis first -- not simplified to one shared branch.
 *
 * NOT on the Step 5 detection call graph (see this file's scope note: the
 * only vendored caller, `intersectArcArc`, is itself out-of-scope
 * arc/PostScript code) -- ported now per the module map anyway, since
 * `rotationAngle.inc` (Milestone A step 9) needs it unchanged, and it is
 * exactly the kind of order-sensitive primitive the fidelity rule flags.
 * NOT parity-exercised by any dump yet; see `intersections_test.cpp` for the
 * targeted hand-computed-reference tests this fidelity rule requires for
 * primitives in that position.
 */
[[nodiscard]] CutPoints get_cut_points_of_circles(Vec2 c1, double r1, Vec2 c2, double r2);

/**
 * The common point(s) of a circle and an infinite line (given by an anchor
 * point and a direction vector). Mirrors `getCutPointsOfCircleAndLine`
 * (`vector_math.inc:722`); a direct application of `solve_square_equation`
 * to the line's parametric form substituted into the circle equation.
 *
 * NOT on the Step 5 detection call graph -- see `get_cut_points_of_circles`'s
 * doc comment (same reasoning; the only vendored caller of a function that
 * calls this, `intersectLineArc`, is out-of-scope arc/PostScript code, and
 * `TestLoopBulge` -- `intersectLevelBoundingBoxes.inc`'s only other caller --
 * is itself dead code, never invoked). Ported now for the same reason.
 */
[[nodiscard]] CutPoints get_cut_points_of_circle_and_line(Vec2 center, double radius, Vec2 anchor,
                                                          Vec2 direction);

/**
 * Whether @p p lies within the `[0, 1]` parametric range of the line through
 * @p p_line with direction @p dir_line (i.e. on the segment from `p_line` to
 * `p_line + dir_line`). Mirrors `matchLinePoint` (`intersectLevelLines.inc:134`),
 * including its `0.0001` (not `EPSILON_7`) axis-degeneracy threshold and its
 * "neither axis is non-degenerate" `return 0` fallback.
 */
[[nodiscard]] bool match_line_point(Vec2 p_line, Vec2 dir_line, Vec2 p);

/**
 * Whether two circles overlap (their boundaries are closer than the sum of
 * their radii). Mirrors `intersectCircleCircle` (`intersectLevelLines.inc:157`).
 */
[[nodiscard]] bool intersect_circle_circle(Vec2 c1, double r1, Vec2 c2, double r2);

/**
 * Whether segment `AB` intersects segment `XY`. Mirrors `intersectLineSegments`
 * (`intersectLevelLines.inc:174`): the AABB-reject fast path, the parallel
 * (`fabs(denominator) < EPSILON_7`) collinear-overlap case, and the general
 * parametric-intersection case, in that order.
 *
 * The vendored signature also returns the cut point via an out-param `P`
 * (`double P[2]`); EVERY call site in the whole vendored codebase reachable
 * from the puzzler entry point passes `NULL` for it (verified by `grep` --
 * the one non-`NULL`-eligible caller, `RNApuzzler.c`'s
 * `checkRemainingIntersections`, is itself out-of-scope debug/arc code), so
 * this port drops the unused out-param rather than declaring-and-ignoring it.
 */
[[nodiscard]] bool intersect_line_segments(Vec2 a, Vec2 b, Vec2 x, Vec2 y);

/**
 * The point on the segment `AB` closest to @p p (clamped projection).
 * Mirrors `projectPointOntoLine` (`intersectLevelBoundingBoxes.inc:111`).
 */
[[nodiscard]] Vec2 project_point_onto_line(Vec2 a, Vec2 b, Vec2 p);

/**
 * The point on triangle `ABC`'s boundary closest to @p p, for the
 * RNApuzzler-specific case where `ABC` is a bulge's near-equilateral,
 * no-obtuse-angle triangle (see the vendored comment this preserves).
 * Mirrors `ClosestPtPointBulge` (`intersectLevelBoundingBoxes.inc:144`): the
 * three "is @p p on the outer side of this edge" tests via
 * `is_to_the_right_point_point`, each short-circuiting to
 * `project_point_onto_line` on that edge, with the "inside the triangle"
 * fallback returning @p p unchanged.
 */
[[nodiscard]] Vec2 closest_pt_point_bulge(Vec2 p, Vec2 a, Vec2 b, Vec2 c);

/**
 * The point on stem @p stem's oriented rectangle (its `a`/`b`/`c`/`e`
 * fields) closest to @p p. Mirrors `ClosestPtPointOBB`
 * (`intersectLevelBoundingBoxes.inc:317`), including its sign/abs-based
 * clamp (not `std::clamp`) -- preserved for expression-tree parity.
 */
[[nodiscard]] Vec2 closest_pt_point_obb(const StemBox& stem, Vec2 p);

/**
 * Whether the circle at @p circle_center with radius @p circle_radius
 * intersects (touches or overlaps) triangle `ABC`. Mirrors
 * `TestCircleTriangle` (`intersectLevelBoundingBoxes.inc:446`):
 * `closest_pt_point_bulge` then a squared-distance compare against the
 * squared radius.
 */
[[nodiscard]] bool test_circle_triangle(Vec2 circle_center, double circle_radius, Vec2 a, Vec2 b,
                                        Vec2 c);

}  // namespace rna_layout::geom
