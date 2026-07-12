/**
 * @file turtle.cpp
 * @brief Implements `layout_turtle` (`turtle.hpp`), ported from
 *        `coordinates.inc`'s `computeAffineCoordinates` /
 *        `affineToCartesianCoordinates` and their `handle*`/`detect*` helpers.
 *
 * Preserves the reference's arithmetic and operation order throughout (per
 * `.claude/plans/current-plan.md`'s style notes): this is a deterministic
 * pass, and the parity gate (`tests/test_native_parity.py`) compares its
 * output to the vendored C bit-for-bit (<=1e-6 after alignment).
 */

#include "rna_layout/turtle.hpp"

#include <cmath>
#include <stdexcept>

#include "rna_layout/config.hpp"
#include "rna_layout/geometry.hpp"
#include "rna_layout/pair_table.hpp"

namespace rna_layout {

namespace {

// Mutually recursive: a stem's loop can contain child stems.
void handle_stem(const std::vector<int>& pair_table, int i, double paired, double unpaired,
                 std::vector<BaseInfo>& base_info, const std::vector<Config>& configs,
                 int direction);

/*===========================================================================
 *  Exterior-loop / dangle handling
 *==========================================================================*/

/// Assign angles/types to a run of unpaired exterior-loop bases starting at
/// `current_base`, up to (and including the entry angle of) the next
/// paired base. Ported from `handleExteriorBases` (`coordinates.inc:42`).
///
/// @return The first paired base reached (or `length` if none remain).
int handle_exterior_bases(const std::vector<int>& pair_table, int current_base,
                          std::vector<BaseInfo>& base_info, int direction) {
  const int length = pair_table[0];

  if (current_base > 1) {
    base_info[current_base].angle += direction * geom::kPiHalf;
    base_info[current_base].type = BaseType::exterior;
  }

  while (current_base < length && pair_table[current_base] <= 0) {
    base_info[current_base + 1].angle = 0.0;
    base_info[current_base].type = BaseType::exterior;
    ++current_base;
  }

  // The following stem cannot alter the first two bases of a stem for its
  // own exit angle, so that angle is written here instead.
  if (current_base + 1 <= length) {
    base_info[current_base + 1].angle = direction * geom::kPiHalf;
    base_info[current_base].type = BaseType::exterior;
  } else {
    base_info[current_base].type = BaseType::exterior;
  }
  return current_base;
}

/*===========================================================================
 *  Bulge detection (loops with exactly one unpaired base)
 *==========================================================================*/

/// Count `m` (base pairs) and `n` (consecutive/backbone steps) around the
/// loop opened by `pair_table[i]`, used to detect the single-unpaired-base
/// bulge special case. Ported from `countLoopPairs` (`coordinates.inc:77`).
void count_loop_pairs(const std::vector<int>& pair_table, int i, int& m, int& n) {
  const int end = pair_table[i];
  ++i;
  m = n = 1;
  while (i < end) {
    if (pair_table[i] <= 0 || pair_table[i] < i) {
      ++n;
      ++i;
    } else {
      ++m;
      i = pair_table[i];
    }
  }
}

/// Detect whether the loop starting at `start` is a bulge (an internal
/// loop with a single unpaired base on one strand), returning the position
/// of the base pair right after the bulge, or 0 if this is not a bulge.
///
/// Ported verbatim from `detectBulge` (`coordinates.inc:117`): a small
/// state machine walking the loop and giving up as soon as the shape
/// stops looking like a bulge. Preserved as one function -- it has no
/// natural sub-seams and is a direct, well-tested port of the reference.
int detect_bulge(int start, const std::vector<int>& pair_table) {
  int bulge = 0;
  const int end = pair_table[start];
  int iterate = 1;
  int old = 0;
  int i = start + 1;
  int side = 0;

  // Reference is a `do { ... } while (i > start)` (checks the loop-around
  // condition at the END, since `i` always starts `> start`); written as
  // `while (true)` + a trailing break for the same "check at the end"
  // shape without the do-while form clang-tidy flags
  // (cppcoreguidelines-avoid-do-while).
  while (true) {
    if (pair_table[i] > 0) {
      if (iterate > 0) {
        if (pair_table[i] == old) {
          ++side;
          ++i;
        } else if (start == pair_table[i] || pair_table[i] == end - 2) {
          bulge = pair_table[i];
          break;
        } else {
          break;
        }
      } else {
        ++iterate;
        old = i;
        i = pair_table[i];
      }
    } else if (iterate > 0) {
      iterate = 0;
      ++i;
    } else {
      ++i;
    }
    if (i <= start) {
      break;
    }
  }

  (void)side;  // computed for parity with the reference; unused thereafter
  return bulge;
}

/*===========================================================================
 *  Loop drawing: bulge case
 *==========================================================================*/

/// The bulge sub-case where the single unpaired base sits on the strand
/// entering the loop (`pair_table[i + 1] == 0`): angle it directly, then
/// recurse into the far stem. Ported from `handleLoop`'s bulge branch,
/// `pair_table[i + 1] == 0` arm (`coordinates.inc:192-210`).
void handle_bulge_unpaired_first(int i, const std::vector<int>& pair_table, double paired,
                                 double unpaired, double alpha, int direction,
                                 std::vector<BaseInfo>& base_info,
                                 const std::vector<Config>& configs) {
  base_info[i + 1].angle += direction * alpha;
  base_info[i].type = BaseType::bulge;
  base_info[pair_table[i]].type = BaseType::bulge;
  ++i;

  base_info[i + 1].angle = -direction * alpha * 2;
  base_info[i].type = BaseType::bulge;
  ++i;
  if (i + 1 <= pair_table[0]) {
    base_info[i + 1].angle = direction * alpha;
  }
  base_info[i].type = BaseType::bulge;
  base_info[pair_table[i]].type = BaseType::bulge;

  handle_stem(pair_table, i, paired, unpaired, base_info, configs, direction);
}

/// The bulge sub-case where the single unpaired base sits on the far
/// strand: recurse into the near stem first, then angle the unpaired base
/// on the way back out. Ported from `handleLoop`'s bulge branch, `else`
/// arm (`coordinates.inc:211-239`).
void handle_bulge_unpaired_second(int i, const std::vector<int>& pair_table, double paired,
                                  double unpaired, double alpha, int direction,
                                  std::vector<BaseInfo>& base_info,
                                  const std::vector<Config>& configs) {
  base_info[i].type = BaseType::bulge;
  ++i;

  // Reference writes `baseInformation[i + 1].baseType = TYPE_BULGE` twice in
  // a row here (interleaved with two angle `+= 0.0` no-ops); collapsed to
  // one assignment -- same observable state, one fewer redundant write.
  base_info[i + 1].type = BaseType::bulge;

  handle_stem(pair_table, i, paired, unpaired, base_info, configs, direction);
  i = pair_table[i];

  base_info[i + 1].angle += direction * alpha;
  base_info[i].type = BaseType::bulge;
  ++i;

  base_info[i + 1].angle = -direction * alpha * 2;
  base_info[i].type = BaseType::bulge;
  ++i;
  if (i + 1 <= pair_table[0]) {
    base_info[i + 1].angle = direction * alpha;
  }
  base_info[i].type = BaseType::bulge;
}

/// Draw a bulge loop (a loop with exactly one unpaired base between its two
/// enclosed stems): the bulge's unpaired base gets a computed opening
/// angle `alpha` instead of a full loop `Config`. Ported from `handleLoop`'s
/// `bulge > 0 && n - m == 1` branch (`coordinates.inc:186-239`).
void handle_loop_bulge_case(int i, const std::vector<int>& pair_table, double paired,
                            double unpaired, int direction, std::vector<BaseInfo>& base_info,
                            const std::vector<Config>& configs) {
  // Reference computes this as `(unpaired * (n - m + 1)) / 2` truncated to
  // int; since this branch only fires when `n - m == 1`, that reduces
  // exactly to `unpaired` (kept as the literal expression for parity).
  const int chord_length = static_cast<int>((unpaired * 2) / 2);
  const double alpha = std::acos(unpaired / (2 * chord_length));

  if (pair_table[i + 1] == 0) {
    handle_bulge_unpaired_first(i, pair_table, paired, unpaired, alpha, direction, base_info,
                                configs);
  } else {
    handle_bulge_unpaired_second(i, pair_table, paired, unpaired, alpha, direction, base_info,
                                 configs);
  }
}

/*===========================================================================
 *  Loop drawing: general (non-bulge) case
 *==========================================================================*/

/// The per-arc geometry a loop's `Config` implies: how much of the loop's
/// full circle a single backbone step spans, and the resulting chord
/// distance/turn angles. Ported from the setup block that opens
/// `handleLoop`'s general-case branch and repeats each time a new arc
/// starts (`coordinates.inc:261-270` and `:318-324`).
struct ArcGeometry {
  double angle_over_paired = 0.0;
  double bb_angle = 0.0;
  double distance = 0.0;
  double delta_ab = 0.0;
  double delta_bb = 0.0;
};

ArcGeometry compute_arc_geometry(const Config& cfg, int arc_index, double paired) {
  const double r = cfg.radius;
  ArcGeometry geo;
  geo.angle_over_paired = 2 * std::asin(paired / (2 * r));
  const double current_angle = cfg.arcs[arc_index].angle;
  geo.bb_angle = (current_angle - geo.angle_over_paired) / cfg.arcs[arc_index].segments;
  geo.distance = std::sqrt(2 * r * r * (1 - std::cos(geo.bb_angle)));
  geo.delta_ab = 0.5 * (geom::kPi + geo.angle_over_paired + geo.bb_angle);
  geo.delta_bb = geom::kPi + geo.bb_angle;
  return geo;
}

/// Draw a general (non-bulge) loop: walk its bases, alternating between
/// unpaired backbone steps, stems opening (recurse via `handle_stem`), and
/// stems closing back into the loop (advance to the next `Config` arc).
/// Ported from `handleLoop`'s `else` branch, the "Loop Drawing Algorithm"
/// block (`coordinates.inc:243-347`). Kept as one function -- its three
/// per-base cases share too much loop-local state (`geo`, `current_arc`,
/// `current_stem_count`) to split into free helpers without turning that
/// state into a pile of out-parameters, which would read worse.
void handle_loop_general_case(int i, const std::vector<int>& pair_table, double paired,
                              double unpaired, std::vector<BaseInfo>& base_info,
                              const std::vector<Config>& configs, int direction) {
  const int start = i;
  const int end = pair_table[i];
  // `.value()` (throws, rather than `*` which is UB on empty) is a cheap
  // safety net: by construction every non-bulge loop-opening base was
  // given a `loop_id` by `generate_config`, but nothing in the type system
  // proves that at this call site (hence the NOLINT: clang-tidy cannot see
  // that invariant either, and `.value()` is already its own recommended
  // fix for the unchecked-`*`-access case).
  const Config& cfg =
      configs[base_info[start].loop_id.value()];  // NOLINT(bugprone-unchecked-optional-access)

  int current_arc = 0;
  ArcGeometry geo = compute_arc_geometry(cfg, current_arc, paired);
  ++current_arc;

  base_info[i + 1].angle += direction * (geom::kPi - geo.delta_ab);
  base_info[i].distance = geo.distance;

  base_info[i].type = (base_info[i].type == BaseType::loop1) ? BaseType::loop2 : BaseType::loop1;
  ++i;

  int current_stem_count = 0;
  while (i < end) {
    if (pair_table[i] <= 0) {
      // Unpaired backbone step within the loop's current arc.
      base_info[i + 1].angle = -direction * (geo.delta_bb - geom::kPi);
      base_info[i].distance = geo.distance;
      base_info[i].type = BaseType::loop1;
      ++i;
    } else if (pair_table[i] > i) {
      // A stem opens here: it ends this arc; recurse to draw it, then jump
      // past it to resume this loop at its far side.
      base_info[i + 1].angle = direction * (geom::kPi - geo.delta_ab);
      ++current_stem_count;
      base_info[i].type = BaseType::loop1;
      handle_stem(pair_table, i, paired, unpaired, base_info, configs, direction);
      i = pair_table[i];
    } else {
      // A stem returns to this loop: if it was a simple pass-through
      // (`current_stem_count == 1`), the loop has reached its next arc.
      if (current_stem_count == 1) {
        current_stem_count = 0;
        geo = compute_arc_geometry(cfg, current_arc, paired);
        ++current_arc;
      }
      base_info[i + 1].angle += direction * (geom::kPi - geo.delta_ab);
      base_info[i].distance = geo.distance;
      base_info[i].type = BaseType::loop1;
      ++i;
    }
  }

  if (i + 1 <= pair_table[0]) {
    base_info[i + 1].angle = direction * (geom::kPi - geo.delta_ab);
  }
  base_info[i].type = BaseType::loop1;
}

/*===========================================================================
 *  Loop / stem dispatch
 *==========================================================================*/

/// Ported from `handleLoop` (`coordinates.inc:169`): dispatch to the bulge
/// or general-case drawing routine based on `detect_bulge`.
void handle_loop(int i, const std::vector<int>& pair_table, double paired, double unpaired,
                 std::vector<BaseInfo>& base_info, const std::vector<Config>& configs,
                 int direction) {
  int m = 0;
  int n = 0;
  count_loop_pairs(pair_table, i, m, n);
  const int bulge = detect_bulge(i, pair_table);

  if (bulge > 0 && n - m == 1) {
    handle_loop_bulge_case(i, pair_table, paired, unpaired, direction, base_info, configs);
  } else {
    handle_loop_general_case(i, pair_table, paired, unpaired, base_info, configs, direction);
  }
}

/// Ported from `handleStem` (`coordinates.inc:356`): walk a stem's paired
/// bases (which need no angle work beyond what the loop/exit already
/// wrote), find its enclosed loop (if any) via `handle_loop`, then mark
/// the stem's second half.
void handle_stem(const std::vector<int>& pair_table, int i, double paired, double unpaired,
                 std::vector<BaseInfo>& base_info, const std::vector<Config>& configs,
                 int direction) {
  const int end = pair_table[i] + 1;

  base_info[i].type = BaseType::stem;
  ++i;

  while (pair_table[i] > 0 &&
         (pair_table[i] == end - 1 || pair_table[i] + 1 == pair_table[i - 1])) {
    base_info[i + 1].angle = 0.0;
    base_info[i].type = BaseType::stem;
    ++i;
  }
  if (pair_table[i] != end - 1) {
    // Reference does `handleLoop(--i, ...)`: the pre-decrement permanently
    // moves `i` back one base (to the last base of the stem's first half)
    // before recursing, and that moved `i` is still live below.
    --i;
    handle_loop(i, pair_table, paired, unpaired, base_info, configs, direction);
  }

  i = pair_table[i];
  base_info[i].type = BaseType::stem;
  ++i;

  while (i < end && i < pair_table[0]) {
    base_info[i].type = BaseType::stem;
    ++i;
  }
}

/*===========================================================================
 *  Top-level affine walk + Cartesian integration
 *==========================================================================*/

/// The singular case where the walk's very first stem starts at base 1:
/// `handle_stem` cannot look one base further back for an exit angle
/// there, so it is seeded directly. Ported from the `currentBase == 1`
/// branch of `computeAffineCoordinates` (`coordinates.inc:446-464`).
///
/// @return True if `current_base == 1` (the caller must treat this like
///     the reference's `continue`: skip the rest of this loop iteration).
bool handle_first_stem(int& current_base, int dangle_count, int length,
                       const std::vector<int>& pair_table, double paired, double unpaired,
                       std::vector<BaseInfo>& base_info, const std::vector<Config>& configs,
                       int direction) {
  if (current_base != 1) {
    return false;
  }
  if (dangle_count < 1) {
    base_info[0].angle = base_info[1].angle = base_info[2].angle = -geom::kPiHalf;
    base_info[current_base].type = BaseType::exterior;
  }
  handle_stem(pair_table, current_base, paired, unpaired, base_info, configs, direction);
  current_base = pair_table[current_base] + 1;
  if (current_base == length) {
    base_info[current_base - 1].type = BaseType::exterior;
    base_info[current_base].type = BaseType::exterior;
    base_info[current_base].angle = -geom::kPiHalf;
  }
  return true;
}

/// Mark the pre-stem exterior dangle when a stem doesn't immediately
/// follow its predecessor. Ported from the `else` arm alongside
/// `handle_first_stem` (`coordinates.inc:465-475`).
void mark_pre_stem_dangle(int current_base, double unpaired, std::vector<BaseInfo>& base_info,
                          int direction, int& dangle_count) {
  base_info[current_base].angle += direction * geom::kPiHalf;
  base_info[current_base + 1].distance = unpaired;
  base_info[current_base - 1].type = BaseType::exterior;
  base_info[current_base + 1].angle += direction * geom::kPiHalf;
  base_info[current_base].type = BaseType::exterior;
  ++dangle_count;
}

/// Ported from `computeAffineCoordinates` (`coordinates.inc:407`): the
/// top-level walk that assigns every base an affine angle/distance by
/// alternating exterior-loop runs with `handle_stem` recursions.
void compute_affine_coordinates(const std::vector<int>& pair_table, double paired, double unpaired,
                                std::vector<BaseInfo>& base_info,
                                const std::vector<Config>& configs) {
  const int length = pair_table[0];
  int current_base = 1;
  constexpr int direction = -1;

  base_info[0].angle = 0.0;
  if (2 <= length) {
    base_info[1].angle = base_info[0].angle;
    base_info[2].angle = base_info[1].angle;
  }

  int dangle_count = 0;
  while (current_base < length) {
    if (pair_table[current_base] <= 0) {
      if (current_base > 1) {
        base_info[current_base - 1].type = BaseType::exterior;
      }
      current_base = handle_exterior_bases(pair_table, current_base, base_info, direction);
      ++dangle_count;
    }
    // Reference re-checks `currentBase < length` here; equivalent to
    // ending this iteration early, since nothing changes `current_base`
    // again before the `while` condition is re-tested.
    if (current_base >= length) {
      break;
    }

    const bool gap_before_stem = pair_table[current_base] - pair_table[current_base - 1] != 1 &&
                                 pair_table[current_base] != 0 && pair_table[current_base - 1] != 0;
    if (gap_before_stem) {
      if (handle_first_stem(current_base, dangle_count, length, pair_table, paired, unpaired,
                            base_info, configs, direction)) {
        continue;
      }
      mark_pre_stem_dangle(current_base, unpaired, base_info, direction, dangle_count);
    }

    handle_stem(pair_table, current_base, paired, unpaired, base_info, configs, direction);
    current_base = pair_table[current_base] + 1;
    if (current_base == length) {
      base_info[current_base - 1].type = BaseType::exterior;
      current_base = handle_exterior_bases(pair_table, current_base, base_info, direction);
    }
  }
  base_info[length].type = BaseType::exterior;
}

/// Ported from `affineToCartesianCoordinates` (`coordinates.inc:500`):
/// integrate each base's affine angle/distance into a cumulative-angle
/// Cartesian walk.
Coords affine_to_cartesian(const std::vector<BaseInfo>& base_info, int length) {
  Coords coords;
  if (length < 1) {
    return coords;
  }
  coords.x.resize(length);
  coords.y.resize(length);

  double angle = 0.0;
  coords.x[0] = coords.y[0] = geom::kExteriorY;
  for (int i = 1; i < length; ++i) {
    angle -= base_info[i + 1].angle;
    coords.x[i] = coords.x[i - 1] + base_info[i].distance * std::cos(angle);
    coords.y[i] = coords.y[i - 1] + base_info[i].distance * std::sin(angle);
  }
  return coords;
}

/*===========================================================================
 *  String-input guards (mirrors src/vienna_layout/bindings.cpp)
 *==========================================================================*/

void validate_nonempty(const std::string& structure) {
  if (structure.empty()) {
    throw std::invalid_argument("layout_turtle requires a non-empty structure");
  }
}

/// A bare `"()"` empty loop segfaults the vendored turtle
/// (`src/vienna_layout/bindings.cpp`'s `validate_no_empty_loop`); reject it
/// uniformly here too so the native engine has the same failure contract.
void validate_no_empty_loop(const std::string& structure) {
  if (structure.find("()") != std::string::npos) {
    throw std::invalid_argument(
        "structure contains an empty loop \"()\": not supported by layout_turtle");
  }
}

}  // namespace

TurtleLayout run_turtle_layout(const std::vector<int>& pair_table, double paired, double unpaired) {
  const int length = pair_table[0];

  TurtleLayout layout;
  layout.base_info.assign(length + 1, BaseInfo{});
  for (BaseInfo& base : layout.base_info) {
    base.distance = unpaired;
  }

  layout.configs = generate_config(pair_table, unpaired, paired, layout.base_info);
  compute_affine_coordinates(pair_table, paired, unpaired, layout.base_info, layout.configs);
  layout.coords = affine_to_cartesian(layout.base_info, length);
  return layout;
}

Coords layout_turtle(const std::vector<int>& pair_table) {
  // Turtle's own hardcoded geometry constants -- see this function's doc
  // comment in `turtle.hpp`.
  constexpr double paired = 35.0;
  constexpr double unpaired = 25.0;
  return run_turtle_layout(pair_table, paired, unpaired).coords;
}

Coords layout_turtle(const std::string& dot_bracket) {
  validate_nonempty(dot_bracket);
  validate_no_empty_loop(dot_bracket);
  const std::vector<int> pair_table = make_pair_table(dot_bracket);  // validates well-nestedness
  return layout_turtle(pair_table);
}

}  // namespace rna_layout
