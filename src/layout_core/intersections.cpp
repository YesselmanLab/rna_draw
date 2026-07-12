/**
 * @file intersections.cpp
 * @brief Implements `intersections.hpp`.
 */

#include "rna_layout/intersections.hpp"

#include <cmath>

#include "rna_layout/geometry.hpp"

namespace rna_layout::geom {

QuadraticSolutions solve_square_equation(double a, double b, double c) {
  // Ported verbatim from `solveSquareEquation` (`vector_math.inc:577-604`):
  // discriminant sign gate, then TWO separate `sqrt(discr)` evaluations (not
  // one shared value with a sign flip) for `answer1`/`answer2`.
  const double discr = b * b - 4 * a * c;
  if (discr < 0.0) {
    return QuadraticSolutions{};
  }

  QuadraticSolutions result;
  result.count = (discr == 0.0) ? 1 : 2;

  const double answer1 = (-b + std::sqrt(discr)) / (2 * a);
  const double answer2 = (-b - std::sqrt(discr)) / (2 * a);
  result.solution1 = answer1;
  result.solution2 = answer2;
  return result;
}

CutPoints get_cut_points_of_circles(Vec2 c1, double r1, Vec2 c2, double r2) {
  // Ported verbatim from `getCutPointsOfCircles` (`vector_math.inc:607-719`).
  const double c1x = c1.x;
  const double c1y = c1.y;
  const double c2x = c2.x;
  const double c2y = c2.y;

  double dx = c1x - c2x;
  dx = dx < 0 ? -dx : dx;
  double dy = c1y - c2y;
  dy = dy < 0 ? -dy : dy;
  double dr = r1 - r2;
  dr = dr < 0 ? -dr : dr;

  // Hardcoded literal, distinct from `kEpsilon7` -- preserved as written.
  const double eps = 1.0;

  const bool small_dx = dx < eps;
  const bool small_dy = dy < eps;
  const bool small_dr = dr < eps;

  CutPoints result;

  if (small_dx && small_dy) {
    result.count = small_dr ? -1 : 0;
    return result;
  }

  if (!small_dy) {
    // (smallDX || !smallDX) && !smallDY: solve for y in terms of x, then
    // substitute into circle1's equation.
    const double k = -2 * c1x + 2 * c2x;
    const double l = c1x * c1x - c2x * c2x + c1y * c1y - c2y * c2y - r1 * r1 + r2 * r2;
    const double m = (-1) * (-2 * c1y + 2 * c2y);

    const double p = c1y - (l / m);
    const double q = k / m;
    const double a = q * q + 1;
    const double b = -2 * c1x - 2 * p * q;
    const double c = c1x * c1x + p * p - r1 * r1;

    const QuadraticSolutions sol = solve_square_equation(a, b, c);
    result.count = sol.count;
    if (sol.count > 0) {
      result.point1 = Vec2{sol.solution1, (sol.solution1 * k + l) / m};
    }
    if (sol.count > 1) {
      result.point2 = Vec2{sol.solution2, (sol.solution2 * k + l) / m};
    }
  } else {
    // smallDY && !smallDX: mirror-image branch, solving for x in terms of y.
    const double k = -2 * c1y + 2 * c2y;
    const double l = (c1x * c1x - c2x * c2x) + (c1y * c1y - c2y * c2y) + (r2 * r2 - r1 * r1);
    const double m = (-1) * (-2 * c1x + 2 * c2x);

    const double p = c1x - (l / m);
    const double q = k / m;
    const double a = q * q + 1;
    const double b = -2 * c1y - 2 * p * q;
    const double c = c1y * c1y + p * p - r1 * r1;

    const QuadraticSolutions sol = solve_square_equation(a, b, c);
    result.count = sol.count;
    if (sol.count > 0) {
      result.point1 = Vec2{(sol.solution1 * k + l) / m, sol.solution1};
    }
    if (sol.count > 1) {
      result.point2 = Vec2{(sol.solution2 * k + l) / m, sol.solution2};
    }
  }

  return result;
}

CutPoints get_cut_points_of_circle_and_line(Vec2 center, double radius, Vec2 anchor,
                                            Vec2 direction) {
  // Ported verbatim from `getCutPointsOfCircleAndLine` (`vector_math.inc:722-750`).
  const double a = direction.x * direction.x + direction.y * direction.y;
  const double b =
      2 * direction.x * (anchor.x - center.x) + 2 * direction.y * (anchor.y - center.y);
  const double c = (anchor.x - center.x) * (anchor.x - center.x) +
                   (anchor.y - center.y) * (anchor.y - center.y) - radius * radius;

  const QuadraticSolutions sol = solve_square_equation(a, b, c);
  CutPoints result;
  result.count = sol.count;
  if (sol.count > 0) {
    result.point1 =
        Vec2{anchor.x + sol.solution1 * direction.x, anchor.y + sol.solution1 * direction.y};
  }
  if (sol.count > 1) {
    result.point2 =
        Vec2{anchor.x + sol.solution2 * direction.x, anchor.y + sol.solution2 * direction.y};
  }
  return result;
}

bool match_line_point(Vec2 p_line, Vec2 dir_line, Vec2 p) {
  // Ported verbatim from `matchLinePoint` (`intersectLevelLines.inc:134-154`):
  // `0.0001`, not `kEpsilon7`, is the reference's own literal here.
  double t = -1.0;
  if (std::fabs(dir_line.x) > 0.0001) {
    t = (p.x - p_line.x) / dir_line.x;
  } else if (std::fabs(dir_line.y) > 0.0001) {
    t = (p.y - p_line.y) / dir_line.y;
  } else {
    return false;
  }
  return 0.0 <= t && t <= 1.0;
}

bool intersect_circle_circle(Vec2 c1, double r1, Vec2 c2, double r2) {
  const Vec2 v_c1_c2 = vector_from_to(c1, c2);
  const double distance = length(v_c1_c2);
  return distance < (r1 + r2);
}

bool intersect_line_segments(Vec2 a, Vec2 b, Vec2 x, Vec2 y) {
  // Ported verbatim from `intersectLineSegments` (`intersectLevelLines.inc:174-277`).
  if ((x.x < a.x - kEpsilon7 && x.x < b.x - kEpsilon7 && y.x < a.x - kEpsilon7 &&
       y.x < b.x - kEpsilon7) ||
      (x.x > a.x + kEpsilon7 && x.x > b.x + kEpsilon7 && y.x > a.x + kEpsilon7 &&
       y.x > b.x + kEpsilon7)) {
    return false;
  }

  if ((x.y < a.y - kEpsilon7 && x.y < b.y - kEpsilon7 && y.y < a.y - kEpsilon7 &&
       y.y < b.y - kEpsilon7) ||
      (x.y > a.y + kEpsilon7 && x.y > b.y + kEpsilon7 && y.y > a.y + kEpsilon7 &&
       y.y > b.y + kEpsilon7)) {
    return false;
  }

  const double denominator = (b.x - a.x) * (x.y - y.y) - (b.y - a.y) * (x.x - y.x);

  if (std::fabs(denominator) < kEpsilon7) {
    // Lines are parallel: check whether X is situated on the AB line.
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    double s_x = 0.0;
    double s_y = 0.0;

    if (std::fabs(dx) > kEpsilon7) {
      s_x = (x.x - a.x) / dx;
      const double ref_xy = a.y + s_x * dy;
      if (std::fabs(ref_xy - x.y) > kEpsilon7) {
        return false;  // AB and XY are not part of the same line.
      }
      s_y = (y.x - a.x) / dx;
    } else {
      s_x = (x.y - a.y) / dy;
      const double ref_xx = a.x + s_x * dx;
      if (std::fabs(ref_xx - x.x) > kEpsilon7) {
        return false;  // AB and XY are not part of the same line.
      }
      s_y = (y.y - a.y) / dy;
    }

    // X or Y situated directly on AB.
    if ((0.0 <= s_x && s_x <= 1.0) || (0.0 <= s_y && s_y <= 1.0)) {
      return true;
    }

    // XY encloses AB.
    if ((s_x < 0.0 && 1.0 < s_y) || (s_y < 0.0 && 1.0 < s_x)) {
      return true;
    }
    return false;
  }

  // Lines are not parallel and might intersect (default case).
  const double nominator_s = (x.x - y.x) * (a.y - x.y) - (x.y - y.y) * (a.x - x.x);
  const double nominator_t = (a.x - x.x) * (b.y - a.y) - (a.y - x.y) * (b.x - a.x);
  const double s = nominator_s / denominator;
  const double t = nominator_t / denominator;

  if (0.0 <= s && s <= 1.0 && 0.0 <= t && t <= 1.0) {
    const Vec2 p_s{a.x + s * (b.x - a.x), a.y + s * (b.y - a.y)};
    const Vec2 p_t{x.x + t * (y.x - x.x), x.y + t * (y.y - x.y)};
    if (std::fabs(p_s.x - p_t.x) < kEpsilon7 && std::fabs(p_s.y - p_t.y) < kEpsilon7) {
      return true;
    }
  }

  return false;
}

Vec2 project_point_onto_line(Vec2 a, Vec2 b, Vec2 p) {
  // Ported verbatim from `projectPointOntoLine` (`intersectLevelBoundingBoxes.inc:111-141`).
  const Vec2 u = vector_from_to(a, p);
  const Vec2 v = vector_from_to(a, b);
  const Vec2 w{-v.y, v.x};

  const double r = (u.y - u.x * w.y / w.x) / (v.y - v.x * w.y / w.x);

  if (r < 0.0) {
    return a;
  }
  if (r > 1.0) {
    return b;
  }
  return Vec2{a.x + r * v.x, a.y + r * v.y};
}

Vec2 closest_pt_point_bulge(Vec2 p, Vec2 a, Vec2 b, Vec2 c) {
  // Ported verbatim from `ClosestPtPointBulge` (`intersectLevelBoundingBoxes.inc:144-203`):
  // three "outer side of this edge" tests, in order AB, BC, CA, each
  // short-circuiting to a projection onto that edge.
  const bool orient_abc = is_to_the_right_point_point(a, b, c);
  const bool orient_abp = is_to_the_right_point_point(a, b, p);
  if (orient_abc != orient_abp) {
    return project_point_onto_line(a, b, p);
  }

  const bool orient_bca = is_to_the_right_point_point(b, c, a);
  const bool orient_bcp = is_to_the_right_point_point(b, c, p);
  if (orient_bca != orient_bcp) {
    return project_point_onto_line(b, c, p);
  }

  const bool orient_cab = is_to_the_right_point_point(c, a, b);
  const bool orient_cap = is_to_the_right_point_point(c, a, p);
  if (orient_cab != orient_cap) {
    return project_point_onto_line(c, a, p);
  }

  // p is inside ABC.
  return p;
}

Vec2 closest_pt_point_obb(const StemBox& stem, Vec2 p) {
  // Ported verbatim from `ClosestPtPointOBB` (`intersectLevelBoundingBoxes.inc:317-350`):
  // the sign/abs clamp (not `std::clamp`) is the reference's own idiom, kept
  // as written for expression-tree parity.
  const Vec2 u0 = stem.a;
  const Vec2 u1 = stem.b;
  const Vec2 dv = vector_from_to(stem.c, p);

  const double dist_0 = dot(dv, u0);
  const double dist_1 = dot(dv, u1);

  const double sign_d0 = dist_0 < 0 ? -1 : 1;
  const double sign_d1 = dist_1 < 0 ? -1 : 1;
  const double sign_e0 = stem.e.x < 0 ? -1 : 1;
  const double sign_e1 = stem.e.y < 0 ? -1 : 1;
  const double abs_d0 = sign_d0 * dist_0;
  const double abs_d1 = sign_d1 * dist_1;
  const double abs_e0 = sign_e0 * stem.e.x;
  const double abs_e1 = sign_e1 * stem.e.y;

  const double clamped_0 = abs_d0 > abs_e0 ? sign_d0 * abs_e0 : sign_d0 * abs_d0;
  const double clamped_1 = abs_d1 > abs_e1 ? sign_d1 * abs_e1 : sign_d1 * abs_d1;

  return Vec2{
      stem.c.x + clamped_0 * stem.a.x + clamped_1 * stem.b.x,
      stem.c.y + clamped_0 * stem.a.y + clamped_1 * stem.b.y,
  };
}

bool test_circle_triangle(Vec2 circle_center, double circle_radius, Vec2 a, Vec2 b, Vec2 c) {
  const Vec2 p = closest_pt_point_bulge(circle_center, a, b, c);
  const Vec2 v = vector_from_to(p, circle_center);
  return dot(v, v) <= (circle_radius * circle_radius);
}

}  // namespace rna_layout::geom
