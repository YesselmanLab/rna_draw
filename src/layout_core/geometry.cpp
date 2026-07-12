/**
 * @file geometry.cpp
 * @brief Implements the 2D math primitives declared in `geometry.hpp`.
 */

#include "rna_layout/geometry.hpp"

#include <cmath>

namespace rna_layout::geom {

namespace {

/// Matches the vendored `EPSILON_7` (`definitions.inc:63`): the tolerance
/// `angle_between` uses to detect a dot product at exactly +-1 before
/// calling `acos`, which is otherwise +-0 derivative there and amplifies
/// round-off into a visible angle error.
constexpr double kEpsilon7 = 1e-7;

}  // namespace

double length(Vec2 v) { return std::sqrt(v.x * v.x + v.y * v.y); }

double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }

double cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }

double angle_between(Vec2 a, Vec2 b) {
  const double len_a = length(a);
  const double len_b = length(b);
  const double cos_angle = dot(a, b) / (len_a * len_b);

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

}  // namespace rna_layout::geom
