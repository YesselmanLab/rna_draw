#pragma once

/**
 * @file config.hpp
 * @brief Per-loop drawing config generation, ported from `drawingconfig.inc`.
 *
 * Scope note: only the "generate a default config for every loop" half of
 * `drawingconfig.inc` is ported this slice (`cfgGenerateConfig` and the
 * radius helpers it needs) -- the turtle-base pass never calls the
 * mutate-in-place half (`cfgApplyChanges`/`cfgIsValid`), which belongs to
 * the resolver (Milestone A step 7+).
 */

#include <vector>

#include "rna_layout/types.hpp"

namespace rna_layout {

/**
 * Generate a default `Config` for every loop in the structure (including
 * the exterior loop's stems, but not the exterior loop itself -- it has no
 * enclosing radius), and annotate each loop's opening base with a
 * `loop_id` that indexes into the returned vector.
 *
 * Ported from `cfgGenerateConfig` (`drawingconfig.inc:578`) plus its
 * recursive helpers `cfgGenHandleStem`/`cfgGenHandleLoop`/
 * `cfgGenerateDefaultConfig`; bulge loops (exactly one unpaired base
 * between exactly two stems) are recognized and skipped, matching the
 * vendored code's choice to fold a bulge's single unpaired base into its
 * enclosing stem rather than giving it its own loop config.
 *
 * @param pair_table 1-indexed pair table (`pair_table[0]` is the sequence
 *     length; `pair_table[i]` is `i`'s partner, 1-indexed, or 0 if
 *     unpaired), matching `make_pair_table`'s output.
 * @param unpaired Default backbone-step distance between consecutive
 *     unpaired (or loop-adjacent) bases.
 * @param paired Default distance between the two bases of a base pair.
 * @param base_info In/out per-base array, size `pair_table[0] + 1`; every
 *     loop-opening base's `loop_id` is set to its `Config`'s index in the
 *     returned vector. Every other field is left untouched.
 * @return One `Config` per loop, in the order `cfgGenHandleLoop` visits
 *     them (parent loops before their children's siblings, matching the
 *     reference's traversal -- later steps rely on this order for parity).
 */
[[nodiscard]] std::vector<Config> generate_config(const std::vector<int>& pair_table,
                                                  double unpaired, double paired,
                                                  std::vector<BaseInfo>& base_info);

/**
 * Approximate the radius of a circle required to draw @p stems base-pair
 * chords (length @p paired) and @p backbones unpaired chords (length
 * @p unpaired) spread evenly over @p angle radians, via Newton iteration.
 *
 * Ported from `approximateConfigArcRadius` (`drawingconfig.inc:270`);
 * preserves the reference's bracketing (`[0.5*paired, 0.5*unpaired]`-based
 * lower/upper bounds) and iteration/clamp order for parity.
 *
 * @param paired Distance spanned by a paired-base chord.
 * @param unpaired Distance spanned by an unpaired-base chord.
 * @param stems Number of paired-base chords in the arc.
 * @param backbones Number of unpaired chords in the arc.
 * @param angle Total angle, radians, the arc must span.
 * @return The fitted radius.
 */
[[nodiscard]] double approximate_config_arc_radius(double paired, double unpaired, int stems,
                                                   int backbones, double angle);

/**
 * The radius that best fits every arc of @p cfg without compressing or
 * stretching any of them: the max of `approximate_config_arc_radius` (with
 * `stems = 1`) over each of the config's arcs.
 *
 * Ported from `approximateConfigRadius` (`drawingconfig.inc:333`). Not
 * called by `layout_turtle` (only by the resolver's `cfgUpdateMinRadius`,
 * Milestone A step 7+); ported alongside `generate_config` per the plan's
 * Step 2 parity gate, which validates both together.
 *
 * @param cfg The loop config whose arcs to fit.
 * @param unpaired Default backbone-step distance.
 * @param paired Default paired-base distance.
 * @return The fitted radius, `>= 0`.
 */
[[nodiscard]] double approximate_config_radius(const Config& cfg, double unpaired, double paired);

}  // namespace rna_layout
