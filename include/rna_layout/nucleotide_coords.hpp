#pragma once

/**
 * @file nucleotide_coords.hpp
 * @brief Final per-nucleotide coordinate derivation from an already-built
 *        config tree, ported from `determineNucleotideCoordinates`
 *        (`RNApuzzler.c:188`) -- Milestone A step 6 ("finalization with
 *        resolver OFF").
 *
 * This turns a `TreeNode` tree's `Config`/`StemBox`/`LoopBox` geometry into
 * one `(x, y)` per nucleotide: stem bases are placed along the stem
 * rectangle's two long edges (with bulges cut in at their own triangle
 * peaks), loop-arc bases are placed on the loop's circle, and exterior
 * (unpaired, tree-root-level) bases walk left to right along the fixed
 * `geom::kExteriorY` baseline. It reads the tree; it never mutates it -- the
 * (not-yet-ported) resolver is what changes the tree's `Config`/boxes
 * BEFORE this runs (Milestone A steps 7-9).
 *
 * HOISTED REDUNDANCY (a deliberate, PROVEN-safe deviation from the
 * vendored control flow -- not a math simplification the fidelity rule
 * would forbid): the vendored `determineNucleotideCoordinates` is
 * recursive over the WHOLE tree, and its final "exterior" block (the
 * baseline walk) is the LAST statement in that SAME recursive function --
 * so it re-executes once per TREE NODE (root's call, then again for every
 * child's call, grandchild's call, ...), each time recomputing the exact
 * same values from the exact same inputs (`pair_table`/`unpairedDistance`/
 * `EXTERIOR_Y` alone -- it reads no node-specific state). This port keeps
 * that block's own internal arithmetic expression-for-expression, but runs
 * it exactly ONCE, after the (also-preserved) recursive stem/loop walk
 * finishes -- provably identical final coordinates, without the O(n^2)
 * (n = tree node count) redundant recomputation.
 */

#include "rna_layout/tree.hpp"
#include "rna_layout/turtle.hpp"

namespace rna_layout {

/**
 * Compute every nucleotide's final `(x, y)` from @p root's tree, writing
 * into @p coords (which must already be sized to `pair_table[0]`, e.g. via
 * `Coords{std::vector<double>(length), std::vector<double>(length)}`).
 *
 * @param root The tree's root (exterior loop).
 * @param pair_table 1-indexed pair table.
 * @param unpaired Default backbone-step distance (exterior baseline walk).
 * @param paired Distance between the two bases of a base pair (loop-arc
 *     start-angle offset).
 * @param coords Output; overwritten in place.
 */
void determine_nucleotide_coordinates(const TreeNode& root, const std::vector<int>& pair_table,
                                      double unpaired, double paired, Coords& coords);

}  // namespace rna_layout
