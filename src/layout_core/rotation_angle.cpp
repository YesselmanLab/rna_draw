/**
 * @file rotation_angle.cpp
 * @brief The shared circle/rectangle rotation-angle solvers, ported from
 *        `rotationAngle.inc`: `point_to_angle` (a private helper local to
 *        this file), `fix_intersection_of_rectangle_and_circle`, and
 *        `fix_intersection_of_circles` -- the two geometry solvers every
 *        `get_rotation_angle_*` type-specific function in
 *        `rotation_angle2.cpp` reduces to.
 *
 * FIDELITY NOTE (`.claude/plans/current-plan.md`'s "CRITICAL" section): this
 * is the densest, most order-sensitive geometry in the whole port -- every
 * expression here preserves the reference's exact operation order, not just
 * its math.
 */

#include "resolve_internal.hpp"

#include <cmath>

#include "rna_layout/geometry.hpp"
#include "rna_layout/intersections.hpp"

namespace rna_layout {

namespace {

/// The signed angle (about @p center, relative to reference direction
/// @p v_ref) to @p point, oriented so a rotation of @p rotation_sign moves
/// monotonically from `0` outward. Ported from `pointToAngle`
/// (`rotationAngle.inc:25`).
double point_to_angle(Vec2 center, Vec2 v_ref, short rotation_sign, Vec2 point) {
  const Vec2 v_center_to_point = geom::vector_from_to(center, point);
  double angle = geom::angle_between(v_ref, v_center_to_point);
  const bool cw = geom::is_to_the_right_point_vector(center, v_ref, point);

  if (rotation_sign > 0 && cw) {
    // angle = angle;
  } else if (rotation_sign > 0 && !cw) {
    angle = geom::kTwoPi - angle;
  } else if (rotation_sign < 0 && cw) {
    angle = -geom::kTwoPi + angle;
  } else if (rotation_sign < 0 && !cw) {
    angle = -angle;
  }
  return angle;
}

}  // namespace

double fix_intersection_of_rectangle_and_circle(Vec2 static_rect_center, Vec2 static_rect_vec_a,
                                                 Vec2 static_rect_vec_b,
                                                 double /*static_rect_length_a*/,
                                                 double static_rect_length_b,
                                                 Vec2 mobile_circ_center, double mobile_circ_radius,
                                                 Vec2 rotation_center, short rotation_sign,
                                                 double clearance) {
  // Ported from `fixIntersectionOfRectangleAndCircle` (`rotationAngle.inc:64`).
  // `static_rect_length_a` is accepted (matching the vendored signature,
  // "used for debug only" per its own comment) but never read -- the
  // vendored function itself never reads it either.
  if (rotation_sign == 0) {
    return 0.0;
  }

  const double distance = geom::epsilon_fix(clearance) + mobile_circ_radius;

  const Vec2 v_rotation_center_to_in_point = geom::vector_from_to(rotation_center, mobile_circ_center);
  const double rotation_radius = geom::length(v_rotation_center_to_in_point);

  const double axis_offset = static_rect_length_b + distance;
  const Vec2 axis_direction = static_rect_vec_a;
  const Vec2 axis_anchor_positive{
      static_rect_center.x + axis_offset * static_rect_vec_b.x,
      static_rect_center.y + axis_offset * static_rect_vec_b.y,
  };
  const Vec2 axis_anchor_negative{
      static_rect_center.x - axis_offset * static_rect_vec_b.x,
      static_rect_center.y - axis_offset * static_rect_vec_b.y,
  };

  std::vector<Vec2> cut;
  const geom::CutPoints positive = geom::get_cut_points_of_circle_and_line(
      rotation_center, rotation_radius, axis_anchor_positive, axis_direction);
  if (positive.count > 0) {
    cut.push_back(positive.point1);
  }
  if (positive.count > 1) {
    cut.push_back(positive.point2);
  }
  const geom::CutPoints negative = geom::get_cut_points_of_circle_and_line(
      rotation_center, rotation_radius, axis_anchor_negative, axis_direction);
  if (negative.count > 0) {
    cut.push_back(negative.point1);
  }
  if (negative.count > 1) {
    cut.push_back(negative.point2);
  }

  if (cut.empty()) {
    // No cut points found (a known rare scenario the vendored comment
    // documents by example) -- fall back to the two points on the circle
    // closest to the axis lines.
    const Vec2 axis_normal = geom::normal(axis_direction);
    cut.push_back(Vec2{rotation_center.x + rotation_radius * axis_normal.x,
                       rotation_center.y + rotation_radius * axis_normal.y});
    cut.push_back(Vec2{rotation_center.x - rotation_radius * axis_normal.x,
                       rotation_center.y - rotation_radius * axis_normal.y});
  }

  std::vector<double> angles(cut.size());
  for (std::size_t i = 0; i < cut.size(); ++i) {
    angles[i] = point_to_angle(rotation_center, v_rotation_center_to_in_point, rotation_sign, cut[i]);
  }
  for (double& a : angles) {
    if (a == 0.0) {
      a = std::signbit(a) ? geom::kMinNegativeAngle : geom::kMinPositiveAngle;
    }
  }

  double angle = rotation_sign * geom::kTwoPi;
  for (double a : angles) {
    if (rotation_sign > 0.0 && a > 0.0) {
      angle = std::fmin(angle, a);
    }
    if (rotation_sign < 0.0 && a < 0.0) {
      angle = std::fmax(angle, a);
    }
  }

  if (std::fabs(angle) == 0.0 || std::fabs(angle) == geom::kTwoPi) {
    angle = 0.0;
  }
  return angle;
}

double fix_intersection_of_circles(Vec2 static_circle_center, double static_circle_radius,
                                   Vec2 mobile_circle_center, double mobile_circle_radius,
                                   Vec2 rotation_center, short rotation_sign, double clearance) {
  // Ported from `fixIntersectionOfCircles` (`rotationAngle.inc:196`).
  if (rotation_sign == 0) {
    return 0.0;
  }

  const double distance = geom::epsilon_fix(clearance);

  const Vec2 v_rotation_center_to_circle_loop_center =
      geom::vector_from_to(rotation_center, mobile_circle_center);
  const double rotation_radius = geom::length(v_rotation_center_to_circle_loop_center);

  const double extended_static_circle_radius = static_circle_radius + mobile_circle_radius + distance;

  const geom::CutPoints cuts = geom::get_cut_points_of_circles(
      rotation_center, rotation_radius, static_circle_center, extended_static_circle_radius);

  // HARDENED (see this function's declaration doc comment in
  // resolve_internal.hpp): `count <= 0` catches both the vendored
  // `numCutPoints == 0` early-return AND the "circles coincide" `count ==
  // -1` case, whose vendored counterpart reads uninitialized cut-point
  // buffers -- provably never reached on real RNA-structure geometry, but
  // guarded here rather than replicated.
  if (cuts.count <= 0) {
    return 0.0;
  }

  double angle1 = 0.0;
  double angle2 = 0.0;
  {
    const Vec2 v_circle_center_to_cut1 = geom::vector_from_to(rotation_center, cuts.point1);
    angle1 = geom::angle_between(v_rotation_center_to_circle_loop_center, v_circle_center_to_cut1);
    const bool is_cw1 =
        geom::is_to_the_right_point_vector(rotation_center, v_rotation_center_to_circle_loop_center,
                                           cuts.point1);
    if (!is_cw1) {
      angle1 *= -1;
    }
    if (angle1 == 0.0) {
      angle1 = std::signbit(angle1) ? geom::kMinNegativeAngle : geom::kMinPositiveAngle;
    }

    const Vec2 v_circle_center_to_cut2 = geom::vector_from_to(rotation_center, cuts.point2);
    angle2 = geom::angle_between(v_rotation_center_to_circle_loop_center, v_circle_center_to_cut2);
    const bool is_cw2 =
        geom::is_to_the_right_point_vector(rotation_center, v_rotation_center_to_circle_loop_center,
                                           cuts.point2);
    if (!is_cw2) {
      angle2 *= -1;
    }
    if (angle2 == 0.0) {
      angle2 = std::signbit(angle2) ? geom::kMinNegativeAngle : geom::kMinPositiveAngle;
    }

    if (is_cw1 == is_cw2) {
      if (std::fabs(angle1) < std::fabs(angle2)) {
        angle2 = is_cw2 ? angle2 - geom::kTwoPi : geom::kTwoPi - angle2;
      } else {
        angle1 = is_cw1 ? angle1 - geom::kTwoPi : geom::kTwoPi - angle1;
      }
    }
  }

  double rotation_angle = 0.0;
  if (rotation_sign == 1) {
    rotation_angle = std::fmax(angle1, angle2);
  } else if (rotation_sign == -1) {
    rotation_angle = std::fmin(angle1, angle2);
  }
  return rotation_angle;
}

}  // namespace rna_layout
