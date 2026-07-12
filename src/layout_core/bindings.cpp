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
 *
 * `dump_tree` (Milestone A step 4) is PARITY-ONLY instrumentation: it runs
 * the turtle pass + config-tree build + `update_bounding_boxes` and returns
 * the T1 tree/box dump `tests/test_native_parity.py` compares against the
 * vendored oracle's `dump_tree` (`vendor_instrument.c`); no production path
 * calls it.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "rna_layout/config_tree.hpp"
#include "rna_layout/debug_dump.hpp"
#include "rna_layout/intersect_tree.hpp"
#include "rna_layout/pair_table.hpp"
#include "rna_layout/puzzler.hpp"
#include "rna_layout/resolve.hpp"
#include "rna_layout/turtle.hpp"

namespace py = pybind11;

namespace {

using CoordVectors = std::pair<std::vector<double>, std::vector<double>>;

/// `dump_tree` calls `make_pair_table` directly (the pair-table overload of
/// the turtle pass, not `layout_turtle(const std::string&)`), so it needs
/// its own copy of the same two guards that overload applies
/// (`turtle.cpp`'s `validate_nonempty`/`validate_no_empty_loop`) -- kept
/// small and duplicated here rather than exposed from `turtle.cpp`, mirroring
/// `src/vienna_layout/bindings.cpp`'s own precedent of duplicating guards on
/// a directly-callable entry point.
void validate_dump_tree_input(const std::string& structure) {
  if (structure.empty()) {
    throw std::invalid_argument("dump_tree requires a non-empty structure");
  }
  if (structure.find("()") != std::string::npos) {
    throw std::invalid_argument(
        "structure contains an empty loop \"()\": not supported by dump_tree");
  }
}

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

/// `rna_layout::DumpConfigArc`/`DumpConfig`/`DumpLoopBox`/`DumpStemBox` ->
/// the same field-name `py::dict` shape the vendored oracle's
/// `vendor_instrument.c` JSON dump uses, so `test_native_parity.py`
/// compares both sides without a native/vendored-specific code path.
py::object to_python(const rna_layout::DumpTreeNode& node) {
  py::dict entry;
  entry["id"] = node.id;
  entry["parent_id"] = node.parent_id;
  entry["loop_start"] = node.loop_start;
  entry["stem_start"] = node.stem_start;

  entry["cfg"] = py::none();
  if (node.cfg.has_value()) {
    py::list arcs;
    for (const rna_layout::DumpConfigArc& arc : node.cfg->arcs) {
      py::dict arc_dict;
      arc_dict["segments"] = arc.segments;
      arc_dict["angle"] = arc.angle;
      arcs.append(std::move(arc_dict));
    }
    py::dict cfg;
    cfg["radius"] = node.cfg->radius;
    cfg["min_radius"] = node.cfg->min_radius;
    cfg["default_radius"] = node.cfg->default_radius;
    cfg["arcs"] = std::move(arcs);
    entry["cfg"] = std::move(cfg);
  }

  entry["lbox"] = py::none();
  if (node.lbox.has_value()) {
    py::dict lbox;
    lbox["cx"] = node.lbox->cx;
    lbox["cy"] = node.lbox->cy;
    lbox["r"] = node.lbox->r;
    entry["lbox"] = std::move(lbox);
  }

  entry["sbox"] = py::none();
  if (node.sbox.has_value()) {
    py::dict sbox;
    sbox["ax"] = node.sbox->ax;
    sbox["ay"] = node.sbox->ay;
    sbox["bx"] = node.sbox->bx;
    sbox["by"] = node.sbox->by;
    sbox["cx"] = node.sbox->cx;
    sbox["cy"] = node.sbox->cy;
    sbox["ex"] = node.sbox->ex;
    sbox["ey"] = node.sbox->ey;
    sbox["bulge_count"] = node.sbox->bulge_count;
    sbox["bulge_dist"] = node.sbox->bulge_dist;
    entry["sbox"] = std::move(sbox);
  }

  return std::move(entry);
}

/// Run the turtle pass + config-tree build + `update_bounding_boxes` on
/// `structure` and return the T1 tree/box dump as a `list[dict]` (see
/// `to_python`). Mirrors `RNApuzzler.c:421-476`'s setup through
/// `updateBoundingBoxes` -- everything up to (not including) the resolver.
/// Named distinctly from `rna_layout::dump_tree` (which this calls,
/// qualified) -- this one is the Python-facing entry point, registered
/// below as `_layout_core.dump_tree`.
py::list dump_config_tree_binding(const std::string& structure, double paired, double unpaired) {
  validate_dump_tree_input(structure);
  const std::vector<int> pair_table = rna_layout::make_pair_table(structure);
  const rna_layout::TurtleLayout turtle =
      rna_layout::run_turtle_layout(pair_table, paired, unpaired);
  const double bulge_dist = rna_layout::stem_bulge_distance(unpaired);

  std::unique_ptr<rna_layout::TreeNode> tree = rna_layout::build_config_tree(
      pair_table, turtle.base_info, turtle.configs, turtle.coords, bulge_dist);
  rna_layout::update_bounding_boxes(*tree, paired, unpaired);

  py::list result;
  for (const rna_layout::DumpTreeNode& node : rna_layout::dump_tree(*tree)) {
    result.append(to_python(node));
  }
  return result;
}

/// Run the turtle pass + config-tree build + `update_bounding_boxes` on
/// `structure` and return the T0 tree's full intersection detection set
/// (Milestone A step 5) as a `list[tuple[int, int, str]]`, matching the
/// vendored oracle's `dump_detections` JSON shape (`vendor_instrument.c`'s
/// `rnadraw_oracle_dump_detections`) field-for-field so
/// `tests/test_native_parity.py` can compare both sides directly.
py::list dump_detections_binding(const std::string& structure, double paired, double unpaired,
                                 double clearance) {
  validate_dump_tree_input(structure);
  const std::vector<int> pair_table = rna_layout::make_pair_table(structure);
  const rna_layout::TurtleLayout turtle =
      rna_layout::run_turtle_layout(pair_table, paired, unpaired);
  const double bulge_dist = rna_layout::stem_bulge_distance(unpaired);

  std::unique_ptr<rna_layout::TreeNode> tree = rna_layout::build_config_tree(
      pair_table, turtle.base_info, turtle.configs, turtle.coords, bulge_dist);
  rna_layout::update_bounding_boxes(*tree, paired, unpaired);

  py::list result;
  for (const rna_layout::Detection& detection :
       rna_layout::detect_intersections(*tree, clearance)) {
    result.append(py::make_tuple(detection.node1_id, detection.node2_id,
                                 rna_layout::intersection_type_to_string(detection.type)));
  }
  return result;
}

/// `layout_puzzler` with `check_sibling`/`check_ancestor`/`optimize` forced
/// false -- the Python-facing entry point for the resolver-off
/// finalization path (Milestone A step 6), so
/// `tests/test_native_parity.py` can call it without constructing a
/// `PuzzlerOptions` binding (not yet exposed; the resolver-on path is a
/// later Milestone A step). `check_exterior`/`allow_flipping`/`clearance`
/// stay at `rna_layout::PuzzlerOptions`'s own defaults, matching the
/// vendored oracle's `plot_coords_puzzler_resolver_off` counterpart
/// (`src/vienna_layout/bindings.cpp`).
CoordVectors plot_coords_puzzler_resolver_off(const std::string& structure) {
  rna_layout::PuzzlerOptions opts;
  opts.check_sibling = false;
  opts.check_ancestor = false;
  opts.optimize = false;
  rna_layout::Coords coords = rna_layout::layout_puzzler(structure, opts);
  return {std::move(coords.x), std::move(coords.y)};
}

/// `layout_puzzler` with `check_sibling` forced true and
/// `check_ancestor`/`optimize` forced false -- the Python-facing entry
/// point for the SIBLING-only resolver path (Milestone A step 7), so
/// `tests/test_native_parity.py` can call it without constructing a
/// `PuzzlerOptions` binding. `check_exterior`/`allow_flipping`/`clearance`/
/// `max_config_changes` stay at `rna_layout::PuzzlerOptions`'s own
/// defaults, matching the vendored oracle's
/// `plot_coords_puzzler_sibling_only` counterpart
/// (`src/vienna_layout/bindings.cpp`).
CoordVectors plot_coords_puzzler_sibling_only(const std::string& structure) {
  rna_layout::PuzzlerOptions opts;
  opts.check_sibling = true;
  opts.check_ancestor = false;
  opts.optimize = false;
  rna_layout::Coords coords = rna_layout::layout_puzzler(structure, opts);
  return {std::move(coords.x), std::move(coords.y)};
}

/// `rna_layout::ChangeTraceEntry` -> the same field-name `py::dict` shape
/// the vendored oracle's `vendor_instrument.c` change-trace JSON dump uses
/// (Milestone A step 7's change-trace parity gate), so
/// `tests/test_native_parity.py` compares both sides without a
/// native/vendored-specific code path.
py::object to_python(const rna_layout::ChangeTraceEntry& entry) {
  py::dict result;
  result["node_id"] = entry.node_id;
  result["type"] = rna_layout::intersection_type_to_string(entry.type);
  result["deltas"] = entry.deltas;
  result["accepted"] = entry.accepted;
  return std::move(result);
}

/// Run the turtle pass + config-tree build + `update_bounding_boxes` +
/// SIBLING-only `check_and_fix_intersections` on `structure`, and return the
/// ORDERED sequence of config-change decisions the resolver made (Milestone
/// A step 7's `dump_change_trace` parity seam) as a `list[dict]`. Mirrors
/// `plot_coords_puzzler_sibling_only`'s pipeline but exposes the resolver's
/// internal decisions rather than only the final coordinates.
py::list dump_change_trace_binding(const std::string& structure, double paired, double unpaired,
                                   double clearance, int max_config_changes) {
  validate_dump_tree_input(structure);
  const std::vector<int> pair_table = rna_layout::make_pair_table(structure);
  const rna_layout::TurtleLayout turtle =
      rna_layout::run_turtle_layout(pair_table, paired, unpaired);
  const double bulge_dist = rna_layout::stem_bulge_distance(unpaired);

  std::unique_ptr<rna_layout::TreeNode> tree = rna_layout::build_config_tree(
      pair_table, turtle.base_info, turtle.configs, turtle.coords, bulge_dist);
  rna_layout::update_bounding_boxes(*tree, paired, unpaired);

  rna_layout::PuzzlerOptions opts;
  opts.paired = paired;
  opts.unpaired = unpaired;
  opts.clearance = clearance;
  opts.check_sibling = true;
  opts.check_ancestor = false;
  opts.optimize = false;

  rna_layout::ResolverState state;
  state.max_config_changes = max_config_changes <= 0 ? 25000 : max_config_changes;
  rna_layout::check_and_fix_intersections(tree.get(), opts, state);

  py::list result;
  for (const rna_layout::ChangeTraceEntry& entry : state.trace) {
    result.append(to_python(entry));
  }
  return result;
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

  m.def("dump_tree", &dump_config_tree_binding, py::arg("structure"), py::arg("paired") = 35.0,
        py::arg("unpaired") = 25.0,
        "T1 config-tree/bounding-box dump (post update_bounding_boxes, "
        "pre-resolver) for parity testing against the vendored oracle's "
        "dump_tree; list[dict], one entry per tree node in DFS pre-order "
        "(id == its index; parent_id == -1 for the root).");

  m.def("dump_detections", &dump_detections_binding, py::arg("structure"), py::arg("paired") = 35.0,
        py::arg("unpaired") = 25.0, py::arg("clearance") = 1.0,
        "The full intersection detection set over the T0 tree (Milestone A "
        "step 5) for parity testing against the vendored oracle's "
        "dump_detections; list[tuple[int, int, str]] of (node1_id, "
        "node2_id, type).");

  m.def("plot_coords_puzzler_resolver_off", &plot_coords_puzzler_resolver_off, py::arg("structure"),
        "Lay out a dot-bracket structure with the native RNApuzzler port, "
        "resolver disabled (check_sibling/check_ancestor/optimize all "
        "false; Milestone A step 6) -- for parity testing against the "
        "vendored oracle's plot_coords_puzzler_resolver_off. Returns (x, y).");

  m.def("plot_coords_puzzler_sibling_only", &plot_coords_puzzler_sibling_only, py::arg("structure"),
        "Lay out a dot-bracket structure with the native RNApuzzler port, "
        "SIBLING intersection resolver only (check_sibling true, "
        "check_ancestor/optimize false; Milestone A step 7) -- for parity "
        "testing against the vendored oracle's "
        "plot_coords_puzzler_sibling_only. Returns (x, y).");

  m.def("dump_change_trace", &dump_change_trace_binding, py::arg("structure"),
        py::arg("paired") = 35.0, py::arg("unpaired") = 25.0, py::arg("clearance") = 1.0,
        py::arg("max_config_changes") = 25000,
        "The ordered sequence of config-change decisions the SIBLING "
        "resolver made (Milestone A step 7) for parity testing against the "
        "vendored oracle's dump_change_trace; list[dict] of (node_id, type, "
        "deltas, accepted).");
}
