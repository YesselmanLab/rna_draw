#pragma once

/**
 * @file types.hpp
 * @brief Value types for the owned `rna_layout` port of RNApuzzler/RNAturtle.
 *
 * Replaces the vendored C structs (`src/vienna_layout/vendor/RNApuzzler/headers`)
 * with plain, RAII-friendly value types: no raw pointers between sibling structs,
 * no process-global mutable state (`clearance` is a `PuzzlerOptions` field), and
 * ownership expressed by `std::vector`/`std::optional` rather than `malloc`/`free`.
 *
 * See `.claude/plans/current-plan.md` ("Value types") for the C ->  C++ mapping
 * this file implements.
 */

#include <cstdint>
#include <optional>
#include <vector>

namespace rna_layout {

/**
 * A 2D point or vector, used throughout the layout core in place of the
 * vendored `double[2]` idiom (`vector_math.inc`).
 */
struct Vec2 {
  double x = 0.0;
  double y = 0.0;
};

/**
 * The role a nucleotide plays in the turtle/puzzler layout, mirroring the
 * vendored `BASE_TYPE` enum (`tBaseInformation_struct.h`).
 *
 * Values MUST match the C enum's integer order: downstream code (ported
 * later) compares raw ints against `TYPE_LOOP1`/`TYPE_LOOP2` boundaries.
 */
enum class BaseType : std::uint8_t {
  none = 0,
  exterior = 1,
  stem = 2,
  bulge = 3,
  loop1 = 4,
  loop2 = 5,
};

/**
 * One arc of a loop's `Config`: the angle spanned between two consecutive
 * stems (or a stem and itself, for a hairpin) and how many backbone
 * segments subdivide it. Mirrors `configArc_struct.h`.
 */
struct ConfigArc {
  int segments = 0;
  double angle = 0.0;
};

/**
 * The per-loop drawing configuration: the loop's radius plus one `ConfigArc`
 * per gap between stems around the loop. Mirrors `config_struct.h`.
 */
struct Config {
  double radius = 0.0;
  double min_radius = 0.0;
  double default_radius = 0.0;
  std::vector<ConfigArc> arcs;
};

/**
 * Per-nucleotide layout state, mirroring `tBaseInformation_struct.h`.
 *
 * `loop_id` replaces the vendored `config*` pointer: it indexes into a
 * `std::vector<Config>` owned by the caller (see `turtle.hpp`), removing
 * the aliasing between `tBaseInformation` instances that share a loop.
 */
struct BaseInfo {
  BaseType type = BaseType::none;
  double angle = 0.0;
  double distance = 0.0;
  std::optional<int> loop_id;
};

/**
 * An axis-aligned bounding circle around a loop. Mirrors `boundingBoxes_struct.h`'s
 * `loopBox`. Built by `bounding_boxes.hpp`'s `build_loop_box` (Milestone A step 4);
 * unused by the turtle-base pass itself.
 */
struct LoopBox {
  Vec2 center;
  double radius = 0.0;
};

/**
 * One "notch" a single unpaired base cuts into a `StemBox`'s rectangle,
 * mirroring the vendored `stemBox::bulges[i]` (`double[4]`,
 * `boundingBoxes.inc:452-465`): NOT a 2D point (an earlier plan draft's
 * `std::vector<Vec2>` would silently drop two of the four fields every
 * bulge geometry helper needs) but the tuple `createBulge` actually builds --
 * a strand @ref sign plus three positions along the stem's `a`-axis
 * (projections of the base before/at/after the bulge, via `getA`,
 * `boundingBoxes.inc:412`).
 */
struct Bulge {
  /// +1.0 for a bulge on the stem's "start" strand, -1.0 for its "end"
  /// strand (`boundingBoxes.inc`'s `setBulges`, the two `bSign` call sites).
  double sign = 0.0;
  double a_prev = 0.0;
  double a_this = 0.0;
  double a_next = 0.0;
};

/**
 * An oriented bounding box around a stem, plus any bulge notches cut into it.
 * Mirrors `boundingBoxes_struct.h`'s `stemBox` (`a`/`b` unit directions,
 * `c` center, `e` half-extents). Built by `bounding_boxes.hpp`'s
 * `build_stem_box`; see `LoopBox`.
 */
struct StemBox {
  Vec2 a;
  Vec2 b;
  Vec2 c;
  Vec2 e;
  std::vector<Bulge> bulges;
  double bulge_dist = 0.0;
};

/**
 * An axis-aligned bounding box, mirroring `AABB_struct.h`. Unused by the
 * turtle-base pass; see `LoopBox`.
 */
struct Aabb {
  Vec2 min;
  Vec2 max;
};

/**
 * Levers for the (not-yet-ported) RNApuzzler resolver, mirroring the
 * vendored `vrna_plot_options_puzzler_t` plus the rna_draw-only `clearance`
 * fork lever. `clearance` is a plain field here -- never a process global
 * like the vendored `rnadraw_clearance_value` -- so intersection tolerances
 * (`epsilon_recognize`/`epsilon_fix`, see `geometry.hpp`) are threaded
 * explicitly through `const PuzzlerOptions&` rather than mutated at a
 * distance.
 *
 * Unused by the turtle-base pass (`layout_turtle` takes no options); kept
 * here now so the resolver steps do not need a second foundation pass.
 */
struct PuzzlerOptions {
  double paired = 35.0;
  double unpaired = 25.0;
  bool check_ancestor = true;
  bool check_sibling = true;
  bool check_exterior = true;
  bool optimize = true;
  bool allow_flipping = false;
  int max_config_changes = 25000;
  double clearance = 1.0;
};

}  // namespace rna_layout
