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
 * Scope note: only the primitives the turtle-base pass and its config
 * generation need are ported this slice (angles, rotation, radius/angle
 * conversion). The remaining `vector_math.inc` primitives used only by the
 * resolver (circle/line intersection solvers) are ported alongside the
 * detection-predicate step (Milestone A step 5), which is where they are
 * first exercised.
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
 * `angleBetweenVectors2D` (`vector_math.inc:459`), including its
 * `EPSILON_7` guard against `acos` domain error at +-1 from floating-point
 * round-off.
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

}  // namespace rna_layout::geom
