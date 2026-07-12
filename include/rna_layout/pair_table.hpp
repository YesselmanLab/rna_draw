#pragma once

/**
 * @file pair_table.hpp
 * @brief Owns the dot-bracket -> pair-table conversion, replacing the
 *        vendored `vrna_ptable` dependency (`vrna_compat.c`).
 */

#include <string>
#include <vector>

namespace rna_layout {

/**
 * Convert a well-nested `"()."`-only dot-bracket structure into a
 * 1-indexed pair table: `result[0]` is the sequence length, and
 * `result[i]` (`1 <= i <= length`) is `i`'s 1-indexed partner, or `0` if
 * `i` is unpaired.
 *
 * Ported from `vrna_ptable`/`extract_pairs`
 * (`vendor/vrna_compat.c`, itself a faithful copy of ViennaRNA 2.7.0's
 * `structure_pairtable.c`), restricted to the round-bracket alphabet
 * `"()."` this port supports (no `vrna_pt_pk_get`/angle/curly/square
 * variants -- pseudoknots are out of scope for `rna_layout`, matching
 * `rna_draw.layout.base.is_pseudoknot_free`).
 *
 * @param dot_bracket A `"()."`-only secondary structure string.
 * @return The 1-indexed pair table described above, size
 *     `dot_bracket.size() + 1`.
 * @throws std::invalid_argument If @p dot_bracket contains a character
 *     outside `"()."`, or its parentheses are unbalanced.
 */
[[nodiscard]] std::vector<int> make_pair_table(const std::string& dot_bracket);

}  // namespace rna_layout
