/**
 * @file bindings.cpp
 * @brief pybind11 bindings for the owned `rna_layout` core, exposed to
 *        Python as `rna_draw._layout_core`.
 *
 * Mirrors `src/vienna_layout/bindings.cpp`'s shape (guard -> call -> return
 * `(x, y)`) so the two engines are drop-in comparable from Python, but this
 * module links no ViennaRNA header or symbol at all -- `rna_layout` owns
 * its own pair-table conversion (`pair_table.hpp`) and value types
 * (`types.hpp`).
 *
 * SCOPE (Milestone A step 3): only `layout_turtle` is wired up. The
 * resolver-backed `layout_puzzler` is a later step (`.claude/plans/
 * current-plan.md`'s port order, steps 6-10) and is deliberately absent
 * here rather than declared-and-stubbed.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <utility>
#include <vector>

#include "rna_layout/turtle.hpp"

namespace py = pybind11;

namespace {

using CoordVectors = std::pair<std::vector<double>, std::vector<double>>;

/// Lay out `structure` with the native RNAturtle port; returns `(x, y)`.
/// `rna_layout::layout_turtle` throws `std::invalid_argument` on malformed
/// input, which pybind11 translates to a Python `ValueError` -- the same
/// failure contract `src/vienna_layout/bindings.cpp`'s `plot_coords_turtle`
/// gives (there it is `std::invalid_argument`/`std::runtime_error`, mapped
/// to `ValueError`/`RuntimeError`).
CoordVectors plot_coords_turtle(const std::string& structure) {
  rna_layout::Coords coords = rna_layout::layout_turtle(structure);
  return {std::move(coords.x), std::move(coords.y)};
}

}  // namespace

PYBIND11_MODULE(_layout_core, m) {
  m.doc() =
      "Owned modern-C++ port of RNApuzzler/RNAturtle's layout core "
      "(rna_layout namespace). This build exposes the turtle-base engine "
      "only; the resolver-backed puzzler engine lands in a later "
      "Milestone A step. No ViennaRNA header/runtime dependency.";

  m.def("plot_coords_turtle", &plot_coords_turtle, py::arg("structure"),
        "Lay out a dot-bracket structure with the native RNAturtle port; "
        "returns (x, y).");

  // Named to match the parity-oracle instrumentation seam
  // (`.claude/plans/current-plan.md`'s "Parity oracle harness"): an alias
  // of `plot_coords_turtle`, since turtle has no separate pre-tree dump
  // point -- it terminates before any tree is ever built.
  m.def("dump_turtle", &plot_coords_turtle, py::arg("structure"),
        "T0 turtle-coordinate dump for parity testing against the vendored "
        "oracle's dump_turtle; identical to plot_coords_turtle.");
}
