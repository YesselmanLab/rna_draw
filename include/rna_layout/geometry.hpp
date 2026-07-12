#pragma once

/**
 * @file geometry.hpp
 * @brief 2D math primitives, ported from `vector_math.inc` and the constants
 *        in `definitions.inc`.
 *
 * Everything here is a pure function of its arguments (no globals, no
 * mutable state) -- a deliberate departure from the vendored file, whose
 * `epsilonRecognize`/`epsilonFix` macros read a process-global `clearance`.
 * Those two constants become plain functions of a `clearance` parameter
 * here (`epsilon_recognize`/`epsilon_fix`), matching `PuzzlerOptions::clearance`
 * (`types.hpp`).
 *
 * Scope note: the turtle-base pass's primitives (angles, rotation,
 * radius/angle conversion) plus the config-tree/bounding-box/wedge
 * primitives (`vector`, `normal`, `rotateVectorByAngle`,
 * `isToTheRightPoint{Point,Vector}`, `anglePtPtPt2D`) needed by Milestone A
 * step 4 are ported here. The remaining `vector_math.inc` primitives, used
 * only by the resolver (circle/line intersection solvers, `solveSquare
 * Equation`, `getCutPointsOf*`), are ported alongside the detection-predicate
 * step (Milestone A step 5), which is where they are first exercised.
 */

#include "rna_layout/types.hpp"

namespace rna_layout::geom {

/// pi, exactly as the vendored `MATH_PI` (not `M_PI`, to keep the ported
/// arithmetic bit-for-bit identical to the reference in the deterministic
/// passes -- see `.claude/plans/current-plan.md`'s style notes).
inline constexpr double kPi = 3.141592653589793;
inline constexpr double kPiHalf = kPi / 2.0;
inline constexpr double kTwoPi = 2.0 * kPi;

/// Newton-iteration convergence tolerance for `approximate_config_arc_radius`
/// (`config.cpp`), matching the vendored `EPSILON_3` (`definitions.inc`).
inline constexpr double kEpsilon3 = 1e-3;

/// The vendored `EPSILON_7` (`definitions.inc:63`): the general-purpose
/// floating-point tie-break tolerance used at circle/line-tangency and
/// degenerate-geometry boundaries throughout the config-tree/bounding-box
/// code (`updateBoundingBoxes`'s zero-length-stem guard,
/// `config_tree.cpp`), in addition to `angle_between`'s `acos` domain guard.
inline constexpr double kEpsilon7 = 1e-7;

/// The exterior loop's fixed y-coordinate, matching the vendored
/// `EXTERIOR_Y` (`definitions.inc:18`): the turtle-base pass's affine walk
/// starts every base there (`turtle.cpp`), and `updateBoundingBoxes`
/// (`config_tree.cpp`) re-anchors every exterior-loop child stem to the same
/// line rather than a circular loop center.
inline constexpr double kExteriorY = 100.0;

/**
 * Length of @p v.
 */
[[nodiscard]] double length(Vec2 v);

/**
 * Dot product of @p a and @p b.
 */
[[nodiscard]] double dot(Vec2 a, Vec2 b);

/**
 * 2D cross product (the scalar z-component of the 3D cross product of
 * @p a and @p b extended into the xy-plane).
 */
[[nodiscard]] double cross(Vec2 a, Vec2 b);

/**
 * The angle between two vectors, in `[0, pi]` radians. Mirrors
 * `angleBetweenVectors2D` (`vector_math.inc:459`) with the same
 * normalize-then-dot expression tree (not just the equivalent math) and its
 * `EPSILON_7` guard against `acos` domain error at +-1 from floating-point
 * round-off -- see the fidelity note in `geometry.cpp`.
 */
[[nodiscard]] double angle_between(Vec2 a, Vec2 b);

/**
 * Convert an angle from radians to degrees.
 */
[[nodiscard]] double to_degree(double angle_rad);

/**
 * Convert an angle from degrees to radians.
 */
[[nodiscard]] double to_radian(double angle_deg);

/**
 * Rotate @p point clockwise (for positive @p angle_rad) around
 * @p center, matching the vendored `rotatePointAroundPoint`'s sign
 * convention (`vector_math.inc:538`).
 *
 * @param point Point to rotate.
 * @param center Center of rotation.
 * @param angle_rad Rotation angle, radians; positive is clockwise.
 * @return The rotated point.
 */
[[nodiscard]] Vec2 rotate_point_around_point(Vec2 point, Vec2 center, double angle_rad);

/**
 * The angle between two points on a circle of the given radius, connected
 * by a chord of the given distance. Mirrors `definitions.inc`'s
 * `distanceToAngle`.
 *
 * @param radius Circle radius.
 * @param distance Chord length between the two points.
 * @return Angle in radians, in `[0, pi]`.
 */
[[nodiscard]] double distance_to_angle(double radius, double distance);

/**
 * The chord length subtended by @p angle_rad on a circle of the given
 * radius. Mirrors `definitions.inc`'s `angleToDistance`; the inverse of
 * `distance_to_angle`.
 */
[[nodiscard]] double angle_to_distance(double radius, double angle_rad);

/**
 * The rna_draw-only intersection-recognition tolerance (vendored
 * `epsilonRecognize`, `definitions.inc:59`), as a pure function of
 * @p clearance rather than a process global.
 *
 * @param clearance `PuzzlerOptions::clearance`; `1.0` reproduces the
 *     vendored stock constant (14.0).
 */
[[nodiscard]] double epsilon_recognize(double clearance);

/**
 * The rna_draw-only intersection-resolution target gap (vendored
 * `epsilonFix`, `definitions.inc:60`), as a pure function of @p clearance.
 *
 * @param clearance `PuzzlerOptions::clearance`; `1.0` reproduces the
 *     vendored stock constant (19.0).
 */
[[nodiscard]] double epsilon_fix(double clearance);

/**
 * The vector from @p from to @p to. Mirrors `vector` (`vector_math.inc:753`);
 * named `vector_from_to` here since `vector` collides with `std::vector`.
 */
[[nodiscard]] Vec2 vector_from_to(Vec2 from, Vec2 to);

/**
 * The unit vector perpendicular to @p v, rotated -90 degrees (i.e.
 * `(v.y, -v.x)`, normalized). Mirrors `normal` (`vector_math.inc:763`).
 */
[[nodiscard]] Vec2 normal(Vec2 v);

/**
 * Rotate @p v clockwise (for positive @p angle_rad) about the origin.
 * Mirrors `rotateVectorByAngle` (`vector_math.inc:554`), itself
 * `rotatePointAroundPoint(v, {0,0}, angle)`.
 */
[[nodiscard]] Vec2 rotate_vector_by_angle(Vec2 v, double angle_rad);

/**
 * Whether @p point lies to the right of the directed line from
 * @p line_start to @p line_end. Mirrors `isToTheRightPointPoint`
 * (`vector_math.inc:392`), including its squared-distance-comparison
 * implementation (no `sqrt`) -- preserved for expression-tree parity, not
 * just for speed.
 */
[[nodiscard]] bool is_to_the_right_point_point(Vec2 line_start, Vec2 line_end, Vec2 point);

/**
 * Whether @p point lies to the right of the directed line starting at
 * @p line_start with direction @p line_vector. Mirrors
 * `isToTheRightPointVector` (`vector_math.inc:445`).
 */
[[nodiscard]] bool is_to_the_right_point_vector(Vec2 line_start, Vec2 line_vector, Vec2 point);

/**
 * The angle at @p center between the rays to @p p1 and @p p3. Mirrors
 * `anglePtPtPt2D` (`vector_math.inc:489`): `angle_between(p1 - center,
 * p3 - center)`.
 */
[[nodiscard]] double angle_pt_pt_pt(Vec2 p1, Vec2 center, Vec2 p3);

}  // namespace rna_layout::geom
