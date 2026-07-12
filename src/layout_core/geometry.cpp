/**
 * @file geometry.cpp
 * @brief Implements the 2D math primitives declared in `geometry.hpp`.
 */

#include "rna_layout/geometry.hpp"

#include <cmath>

namespace rna_layout::geom {

double length(Vec2 v) { return std::sqrt(v.x * v.x + v.y * v.y); }

double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }

double cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }

double angle_between(Vec2 a, Vec2 b) {
  // Preserve the reference's exact expression tree (normalize-then-dot),
  // not just its math: `angleBetweenVectors2D` (vector_math.inc:458) calls
  // `normalize()` on each vector (which divides each COMPONENT by its own
  // length in place, vector_math.inc:382-389) and only then dots the unit
  // vectors via `scalarProduct2D` (vector_math.inc:368-379). That is NOT
  // bit-identical to `dot(a, b) / (len_a * len_b)` -- the two forms differ
  // by ~1e-16, which is enough to break the resolver's discrete tie-breaks
  // (EPSILON_7 boundary) under the parity gate.
  const double len_a = length(a);
  const double len_b = length(b);
  const double cos_angle = (a.x / len_a) * (b.x / len_b) + (a.y / len_a) * (b.y / len_b);

  if (std::fabs(cos_angle - -1.0) < kEpsilon7) {
    return kPi;
  }
  if (std::fabs(cos_angle - 1.0) < kEpsilon7) {
    return 0.0;
  }
  return std::acos(cos_angle);
}

double to_degree(double angle_rad) { return angle_rad * (180.0 / kPi); }

double to_radian(double angle_deg) { return angle_deg * (kPi / 180.0); }

Vec2 rotate_point_around_point(Vec2 point, Vec2 center, double angle_rad) {
  // Negated because a positive input angle means "rotate clockwise" here,
  // matching the vendored sign convention (`vector_math.inc:547`).
  const double phi = -angle_rad;
  const double dx = point.x - center.x;
  const double dy = point.y - center.y;
  return Vec2{
      center.x + dx * std::cos(phi) - dy * std::sin(phi),
      center.y + dx * std::sin(phi) + dy * std::cos(phi),
  };
}

double distance_to_angle(double radius, double distance) {
  return 2.0 * std::asin(distance / (2.0 * radius));
}

double angle_to_distance(double radius, double angle_rad) {
  return 2.0 * radius * std::sin(angle_rad / 2.0);
}

double epsilon_recognize(double clearance) { return 14.0 * clearance; }

double epsilon_fix(double clearance) { return 19.0 * clearance; }

Vec2 vector_from_to(Vec2 from, Vec2 to) { return Vec2{to.x - from.x, to.y - from.y}; }

Vec2 normal(Vec2 v) {
  // Reference builds the perpendicular first (`vNormal = {v.y, -v.x}`), THEN
  // normalizes it via its own `vectorLength2D` (`vector_math.inc:763-777`) --
  // not `length(v)` -- kept as two separate steps for expression-tree parity.
  const Vec2 v_normal{v.y, -v.x};
  const double normal_length = length(v_normal);
  return Vec2{v_normal.x / normal_length, v_normal.y / normal_length};
}

Vec2 rotate_vector_by_angle(Vec2 v, double angle_rad) {
  return rotate_point_around_point(v, Vec2{0.0, 0.0}, angle_rad);
}

bool is_to_the_right_point_point(Vec2 line_start, Vec2 line_end, Vec2 point) {
  // Ported verbatim from `isToTheRightPointPoint` (`vector_math.inc:392`):
  // compares SQUARED distances from `point` to the line's left/right normal
  // offsets rather than normalizing/using `sqrt`, by design (see that
  // function's own comment) -- preserved here, not just for speed, but
  // because it is a different (parity-relevant) expression tree than a
  // sign-of-cross-product test would be.
  const Vec2 line_vector{line_end.x - line_start.x, line_end.y - line_start.y};
  const Vec2 normal_vector{line_vector.y, -line_vector.x};
  const Vec2 right{line_end.x + normal_vector.x, line_end.y + normal_vector.y};
  const Vec2 left{line_end.x - normal_vector.x, line_end.y - normal_vector.y};
  const Vec2 v_right{point.x - right.x, point.y - right.y};
  const Vec2 v_left{point.x - left.x, point.y - left.y};
  const double squared_distance_right = dot(v_right, v_right);
  const double squared_distance_left = dot(v_left, v_left);
  return squared_distance_right < squared_distance_left;
}

bool is_to_the_right_point_vector(Vec2 line_start, Vec2 line_vector, Vec2 point) {
  const Vec2 line_end{line_start.x + line_vector.x, line_start.y + line_vector.y};
  return is_to_the_right_point_point(line_start, line_end, point);
}

double angle_pt_pt_pt(Vec2 p1, Vec2 center, Vec2 p3) {
  const Vec2 v1{p1.x - center.x, p1.y - center.y};
  const Vec2 v2{p3.x - center.x, p3.y - center.y};
  return angle_between(v1, v2);
}

}  // namespace rna_layout::geom
