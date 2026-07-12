#pragma once

/**
 * @file turtle.hpp
 * @brief The RNAturtle base layout: the deterministic, resolver-free half of
 *        the RNApuzzler algorithm.
 *
 * Ported from `RNAturtle.c`'s `vrna_plot_coords_turtle_pt`
 * (`cfgGenerateConfig` -> `computeAffineCoordinates` ->
 * `affineToCartesianCoordinates`, `coordinates.inc`). This is the first
 * SHIPPABLE native engine in the `rna_layout` port: turtle terminates
 * before `buildConfigtree`/the intersection resolver ever run, so it is a
 * complete, independently parity-testable vertical slice.
 */

#include <string>
#include <vector>

namespace rna_layout {

/// Per-nucleotide layout output, 0-indexed, one entry per sequence position.
struct Coords {
  std::vector<double> x{};
  std::vector<double> y{};
};

/**
 * Lay out @p pair_table with RNAturtle: generate a default per-loop
 * `Config` (`generate_config`), walk the structure to assign each base an
 * affine angle/distance (`computeAffineCoordinates`), then integrate those
 * into Cartesian coordinates (`affineToCartesianCoordinates`).
 *
 * Uses the vendored turtle's hardcoded geometry constants (`paired = 35`,
 * `unpaired = 25`), matching `vrna_plot_coords_turtle_pt`
 * (`RNAturtle.c:62-63`) -- turtle (unlike puzzler) does not expose these as
 * options.
 *
 * @param pair_table 1-indexed pair table (`make_pair_table`'s output).
 * @return `Coords` of size `pair_table[0]` (the sequence length); empty if
 *     the structure is empty.
 */
[[nodiscard]] Coords layout_turtle(const std::vector<int>& pair_table);

/**
 * Convenience overload: validate @p dot_bracket the same way
 * `src/vienna_layout/bindings.cpp` validates input to the vendored engines
 * (non-empty, well-nested `"()."`-only, no bare `"()"` empty loop -- see
 * that file's guards, which this reproduces verbatim for a matching
 * failure contract), build its pair table, and delegate to the pair-table
 * overload.
 *
 * @param dot_bracket A `"()."`-only secondary structure string.
 * @return `Coords` of size `dot_bracket.size()`.
 * @throws std::invalid_argument If @p dot_bracket is empty, not
 *     well-nested, or contains a bare `"()"` empty loop.
 */
[[nodiscard]] Coords layout_turtle(const std::string& dot_bracket);

}  // namespace rna_layout
