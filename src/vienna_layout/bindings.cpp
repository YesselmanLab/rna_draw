// pybind11 bindings for ViennaRNA's in-process secondary-structure layout
// algorithms (puzzler, turtle), exposed to Python as
// `rna_draw._vienna_layout`.
//
// STANDALONE BUILD: the two vendored layout translation units
// (`vendor/RNApuzzler/{RNApuzzler,RNAturtle}.c`) are compiled directly into
// this extension against a ~120 LOC compat shim
// (`vendor/vrna_compat.c`, providing `vrna_alloc` + `vrna_ptable`), so this
// module links NO `libRNA.a` -- rna_draw has no ViennaRNA runtime
// dependency. (The naview binding was dropped: it lived only in
// `libRNA.a`'s non-reentrant `naview.o`, and per-structure benchmarking
// showed it never rescues a puzzler-dirty structure -- a documented
// dead-end.)
//
// REENTRANCY: both RNApuzzler and RNAturtle are reentrant (no mutable
// file-scope state), so both bindings are safe to call concurrently.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// ViennaRNA's headers do not wrap their declarations in `extern "C"`
// (verified: only 2 of 208 headers in this install do), so a C++
// translation unit including them gets C++-mangled declarations that do
// not match the plain C symbols the vendored layout objects export -- wrap
// the include ourselves. The vendored `RNApuzzler.c`/`RNAturtle.c` define
// these entry points; we only need the installed headers here for their
// declarations (the include paths are still on the build's search path).
// clang-format off
extern "C" {
#include <ViennaRNA/vrna_config.h>
#include <ViennaRNA/plotting/layouts.h>
}
// clang-format on

namespace py = pybind11;

namespace {

using CoordVectors = std::pair<std::vector<double>, std::vector<double>>;

/// RAII owner for a malloc'd buffer returned by a `vrna_plot_coords_*` call.
/// The ViennaRNA C API hands back ownership of `float*`/`double*` buffers
/// allocated with `malloc`; wrapping each in this guard makes the
/// copy-into-`std::vector` path exception-safe and leak-free without
/// hand-rolled `free()` calls at every return.
template <typename T>
struct MallocBuffer {
  T* ptr = nullptr;

  ~MallocBuffer() { std::free(ptr); }
  MallocBuffer() = default;
  MallocBuffer(const MallocBuffer&) = delete;
  MallocBuffer& operator=(const MallocBuffer&) = delete;
  // Scope-local RAII guard only (never returned/stored elsewhere) -- moving
  // it is never needed, so the Rule of 5 is satisfied by deleting it too.
  MallocBuffer(MallocBuffer&&) = delete;
  MallocBuffer& operator=(MallocBuffer&&) = delete;
};

/// Copy `n` entries from two malloc'd coordinate buffers into `std::vector`s.
CoordVectors to_vectors(const MallocBuffer<float>& x, const MallocBuffer<float>& y, int n) {
  std::vector<double> xs(x.ptr, x.ptr + n);
  std::vector<double> ys(y.ptr, y.ptr + n);
  return {std::move(xs), std::move(ys)};
}

/// Reject anything `vrna_plot_coords_*` cannot safely handle.
///
/// Verified (M5.1 spike): `vrna_plot_coords_puzzler` SEGFAULTS -- does not
/// return 0 -- on an unbalanced-parenthesis structure (e.g. `"("`, `"))"`)
/// because it builds a pair table assuming balance and indexes past it.
/// The Python `LayoutEngine` wrappers already guard with
/// `is_pseudoknot_free` before reaching this module, but this module is
/// itself importable and callable directly, so it must not be a segfault
/// oracle: validate well-nestedness (same rule as
/// `rna_draw.layout.base.is_pseudoknot_free`, duplicated here since this
/// translation unit has no Python dependency) before ever calling into the
/// C API.
void validate_well_nested(const std::string& structure) {
  int depth = 0;
  for (char c : structure) {
    if (c != '(' && c != ')' && c != '.') {
      throw std::invalid_argument(
          "structure contains a character outside \"().\": not "
          "well-nested dot-bracket input");
    }
    depth += (c == '(') ? 1 : (c == ')') ? -1 : 0;
    if (depth < 0) {
      throw std::invalid_argument(
          "unbalanced structure: an unmatched ')' would "
          "crash vrna_plot_coords_*");
    }
  }
  if (depth != 0) {
    throw std::invalid_argument(
        "unbalanced structure: unmatched '(' would crash "
        "vrna_plot_coords_*");
  }
}

/// Reject the empty structure before it reaches the C API.
///
/// puzzler/turtle return 0 cleanly for this case (caught by the `n == 0`
/// check after the call); rejecting it up front keeps the guard uniform
/// and self-documenting.
void validate_nonempty(const std::string& structure) {
  if (structure.empty()) {
    throw std::runtime_error("vrna_plot_coords_* requires a non-empty structure");
  }
}

/// Reject a structure containing a bare empty loop (`"()"` with nothing
/// between the pair).
///
/// Verified (M5.1 spike): `vrna_plot_coords_puzzler("().()"))` does not
/// return -- it loops effectively forever (its iterative
/// intersection-resolution never converges on a degenerate zero-nucleotide
/// loop) -- and `vrna_plot_coords_turtle("().()")` SEGFAULTS on the same
/// input. This is the same degenerate shape that already crashes
/// `render_rna.py`'s tree recursion (see
/// `rna_draw.layout.base.has_empty_loop`); guard uniformly across both
/// engines rather than special-case per algorithm.
void validate_no_empty_loop(const std::string& structure) {
  if (structure.find("()") != std::string::npos) {
    throw std::invalid_argument(
        "structure contains an empty loop \"()\": vrna_plot_coords_puzzler "
        "hangs and vrna_plot_coords_turtle segfaults on this shape");
  }
}

/// Call `vrna_plot_coords_puzzler(structure, &x, &y, NULL, NULL)`, copy the
/// result into `std::vector`s, and free the malloc'd buffers.
CoordVectors plot_coords_puzzler(const std::string& structure) {
  validate_nonempty(structure);
  validate_well_nested(structure);
  validate_no_empty_loop(structure);
  MallocBuffer<float> x;
  MallocBuffer<float> y;
  int n = vrna_plot_coords_puzzler(structure.c_str(), &x.ptr, &y.ptr, nullptr, nullptr);
  if (n == 0 || static_cast<size_t>(n) != structure.size()) {
    throw std::runtime_error("vrna_plot_coords_puzzler failed on structure of length " +
                             std::to_string(structure.size()));
  }
  return to_vectors(x, y, n);
}

/// Call `vrna_plot_coords_turtle(structure, &x, &y, &arc_coords)`; turtle
/// (unlike puzzler) always wants a real `arc_coords` pointer, which is
/// freed via the same RAII guard and otherwise discarded.
CoordVectors plot_coords_turtle(const std::string& structure) {
  validate_nonempty(structure);
  validate_well_nested(structure);
  validate_no_empty_loop(structure);
  MallocBuffer<float> x;
  MallocBuffer<float> y;
  MallocBuffer<double> arcs;
  int n = vrna_plot_coords_turtle(structure.c_str(), &x.ptr, &y.ptr, &arcs.ptr);
  if (n == 0 || static_cast<size_t>(n) != structure.size()) {
    throw std::runtime_error("vrna_plot_coords_turtle failed on structure of length " +
                             std::to_string(structure.size()));
  }
  return to_vectors(x, y, n);
}

/// RAII owner for a `vrna_plot_options_puzzler_t*` allocated by
/// `vrna_plot_options_puzzler()`, freed via `vrna_plot_options_puzzler_free`.
struct PuzzlerOptions {
  vrna_plot_options_puzzler_t* ptr = vrna_plot_options_puzzler();

  ~PuzzlerOptions() { vrna_plot_options_puzzler_free(ptr); }
  PuzzlerOptions() = default;
  PuzzlerOptions(const PuzzlerOptions&) = delete;
  PuzzlerOptions& operator=(const PuzzlerOptions&) = delete;
  // Scope-local RAII guard only -- see MallocBuffer's rationale above.
  PuzzlerOptions(PuzzlerOptions&&) = delete;
  PuzzlerOptions& operator=(PuzzlerOptions&&) = delete;
};

extern "C" void rnadraw_set_clearance(double factor);

/// Call `vrna_plot_coords_puzzler(structure, &x, &y, NULL, options)` with an
/// explicit options struct exposing the resolver levers:
/// `allow_flipping` (RNApuzzler's exterior-branch flip heuristic),
/// `max_config_changes` (the config-change search budget; <= 0 falls back
/// to the engine's 25000 default -- see the vendored `RNApuzzler.c` edit
/// that makes this budget caller-respectable instead of hardcoded), and
/// `clearance` (>1 scales puzzler's intersection clearance so it resolves
/// near-touches that rna_draw's stricter M2 checker flags; <= 0 or 1.0 ==
/// stock). Clearance is a process-global set for the duration of this call
/// and reset after (single-threaded per process -- see definitions.inc).
CoordVectors plot_coords_puzzler_opts(const std::string& structure, bool allow_flipping,
                                      int max_config_changes, double clearance) {
  validate_nonempty(structure);
  validate_well_nested(structure);
  validate_no_empty_loop(structure);
  PuzzlerOptions options;
  options.ptr->allowFlipping = allow_flipping ? 1 : 0;
  options.ptr->maximumNumberOfConfigChangesAllowed = max_config_changes;
  MallocBuffer<float> x;
  MallocBuffer<float> y;
  rnadraw_set_clearance(clearance);
  int n = vrna_plot_coords_puzzler(structure.c_str(), &x.ptr, &y.ptr, nullptr, options.ptr);
  rnadraw_set_clearance(0.0);
  if (n == 0 || static_cast<size_t>(n) != structure.size()) {
    throw std::runtime_error("vrna_plot_coords_puzzler_opts failed on structure of length " +
                             std::to_string(structure.size()));
  }
  return to_vectors(x, y, n);
}

/// Call `vrna_plot_coords_puzzler` with `checkSiblingIntersections = 0`,
/// `checkAncestorIntersections = 0`, `optimize = 0` -- the vendored-side
/// counterpart of `rna_layout::layout_puzzler`'s resolver-off finalization
/// path (Milestone A step 6). `checkExteriorIntersections` and
/// `allowFlipping` are left at `vrna_plot_options_puzzler()`'s own defaults
/// (`1` and `0` respectively), matching `rna_layout::PuzzlerOptions`'s
/// defaults for those two fields.
CoordVectors plot_coords_puzzler_resolver_off(const std::string& structure) {
  validate_nonempty(structure);
  validate_well_nested(structure);
  validate_no_empty_loop(structure);
  PuzzlerOptions options;
  options.ptr->checkSiblingIntersections = 0;
  options.ptr->checkAncestorIntersections = 0;
  options.ptr->optimize = 0;
  MallocBuffer<float> x;
  MallocBuffer<float> y;
  int n = vrna_plot_coords_puzzler(structure.c_str(), &x.ptr, &y.ptr, nullptr, options.ptr);
  if (n == 0 || static_cast<size_t>(n) != structure.size()) {
    throw std::runtime_error(
        "vrna_plot_coords_puzzler (resolver off) failed on structure of length " +
        std::to_string(structure.size()));
  }
  return to_vectors(x, y, n);
}

/// Call `vrna_plot_coords_puzzler` with `checkSiblingIntersections = 1`,
/// `checkAncestorIntersections = 0`, `optimize = 0` -- the vendored-side
/// counterpart of `rna_layout::layout_puzzler`'s SIBLING-only resolver path
/// (Milestone A step 7). `checkExteriorIntersections` and `allowFlipping`
/// are left at `vrna_plot_options_puzzler()`'s own defaults, matching
/// `rna_layout::PuzzlerOptions`'s defaults for those two fields.
CoordVectors plot_coords_puzzler_sibling_only(const std::string& structure) {
  validate_nonempty(structure);
  validate_well_nested(structure);
  validate_no_empty_loop(structure);
  PuzzlerOptions options;
  options.ptr->checkSiblingIntersections = 1;
  options.ptr->checkAncestorIntersections = 0;
  options.ptr->optimize = 0;
  MallocBuffer<float> x;
  MallocBuffer<float> y;
  int n = vrna_plot_coords_puzzler(structure.c_str(), &x.ptr, &y.ptr, nullptr, options.ptr);
  if (n == 0 || static_cast<size_t>(n) != structure.size()) {
    throw std::runtime_error(
        "vrna_plot_coords_puzzler (sibling only) failed on structure of length " +
        std::to_string(structure.size()));
  }
  return to_vectors(x, y, n);
}

// Defined in vendor_instrument.c, compiled in only when RNA_DRAW_BUILD_ORACLE
// is on (see CMakeLists.txt); that file documents the (macro-interposition-
// adjacent) mechanism used to reach the vendored tree/box internals, and the
// dump_detections/dump_change_trace entry points still to come. Guarded by
// the same preprocessor define CMake sets for the oracle build, so this
// binding module still links when the option is off.
#ifdef RNA_DRAW_BUILD_ORACLE
extern "C" const char* rnadraw_oracle_instrumentation_version(void);
extern "C" char* rnadraw_oracle_dump_tree(const char* structure, double paired, double unpaired);

/// T1 config-tree/bounding-box dump: calls `rnadraw_oracle_dump_tree`
/// (`vendor_instrument.c`), which returns a malloc'd JSON string; copy it
/// into a `std::string` and free the buffer (same `MallocBuffer` RAII
/// pattern as `plot_coords_*`'s float/double output buffers, specialized
/// for `char*` here since the buffer is NUL-terminated text, not a fixed-
/// length numeric array).
std::string dump_tree(const std::string& structure, double paired, double unpaired) {
  validate_nonempty(structure);
  validate_well_nested(structure);
  validate_no_empty_loop(structure);
  char* json = rnadraw_oracle_dump_tree(structure.c_str(), paired, unpaired);
  if (json == nullptr) {
    throw std::runtime_error("rnadraw_oracle_dump_tree failed on structure of length " +
                             std::to_string(structure.size()));
  }
  std::string result(json);
  std::free(json);
  return result;
}

extern "C" char* rnadraw_oracle_dump_detections(const char* structure, double paired,
                                                double unpaired);
extern "C" char* rnadraw_oracle_dump_change_trace(const char* structure, double paired,
                                                  double unpaired, int max_config_changes);

/// The ordered SIBLING-resolver config-change trace (Milestone A step 7):
/// calls `rnadraw_oracle_dump_change_trace` (`vendor_instrument.c`), which
/// returns a malloc'd JSON string; copy it into a `std::string` and free
/// the buffer (same pattern as `dump_tree`/`dump_detections`).
///
/// CALLER WARNING (see `vendor_instrument.c`'s doc comment on the C
/// function): the vendored resolver does not terminate on every structure
/// under checkSiblingIntersections=1/checkAncestorIntersections=0/
/// optimize=0 -- do not call this on a structure not already known to
/// terminate (`benchmarks/hard_set_sibling_only_oracle_hangs.json` lists
/// the known exceptions on `hard_set.json`) without an external timeout.
std::string dump_change_trace(const std::string& structure, double paired, double unpaired,
                              int max_config_changes) {
  validate_nonempty(structure);
  validate_well_nested(structure);
  validate_no_empty_loop(structure);
  char* json =
      rnadraw_oracle_dump_change_trace(structure.c_str(), paired, unpaired, max_config_changes);
  if (json == nullptr) {
    throw std::runtime_error("rnadraw_oracle_dump_change_trace failed on structure of length " +
                             std::to_string(structure.size()));
  }
  std::string result(json);
  std::free(json);
  return result;
}

/// The intersection detection set (Milestone A step 5): calls
/// `rnadraw_oracle_dump_detections` (`vendor_instrument.c`), which returns
/// a malloc'd JSON string of `[[node1_id, node2_id, "type"], ...]`; copy it
/// into a `std::string` and free the buffer (same pattern as `dump_tree`).
std::string dump_detections(const std::string& structure, double paired, double unpaired) {
  validate_nonempty(structure);
  validate_well_nested(structure);
  validate_no_empty_loop(structure);
  char* json = rnadraw_oracle_dump_detections(structure.c_str(), paired, unpaired);
  if (json == nullptr) {
    throw std::runtime_error("rnadraw_oracle_dump_detections failed on structure of length " +
                             std::to_string(structure.size()));
  }
  std::string result(json);
  std::free(json);
  return result;
}
#endif

/// The compiled-against ViennaRNA version, as a string (e.g. `"2.7.0"`).
std::string version() { return VRNA_VERSION; }

/// The compiled-against ViennaRNA version, as `(major, minor, patch)`.
std::tuple<int, int, int> abi_version() {
  return {VRNA_VERSION_MAJOR, VRNA_VERSION_MINOR, VRNA_VERSION_PATCH};
}

/// `sizeof(vrna_plot_options_puzzler_t)` -- an ABI drift guard. A header
/// bump that reorders/extends the options struct must be caught here, not
/// silently mis-laid-out (see plan Step 2).
size_t sizeof_puzzler_options() { return sizeof(vrna_plot_options_puzzler_t); }

}  // namespace

PYBIND11_MODULE(_vienna_layout, m) {
  m.doc() =
      "In-process bindings to ViennaRNA's puzzler/turtle layout algorithms, "
      "compiled from a vendored, standalone copy of the layout core (no "
      "libRNA.a link). Both engines are reentrant.";

  m.def("version", &version, "Compiled-against ViennaRNA version string.");
  m.def("abi_version", &abi_version,
        "Compiled-against ViennaRNA version as (major, minor, patch).");
  m.def("sizeof_puzzler_options", &sizeof_puzzler_options,
        "sizeof(vrna_plot_options_puzzler_t), for ABI-drift detection.");

  m.def("plot_coords_puzzler", &plot_coords_puzzler, py::arg("structure"),
        "Lay out a dot-bracket structure with RNApuzzler; returns (x, y).");
  m.def("plot_coords_puzzler_opts", &plot_coords_puzzler_opts, py::arg("structure"),
        py::arg("allow_flipping") = false, py::arg("max_config_changes") = 0,
        py::arg("clearance") = 0.0,
        "Lay out a dot-bracket structure with RNApuzzler, exposing resolver "
        "levers: allow_flipping (exterior-branch flip heuristic), "
        "max_config_changes (config-change search budget; <= 0 uses the "
        "engine's default), and clearance (>1 scales the intersection "
        "clearance so puzzler resolves near-touches rna_draw's checker "
        "flags; <= 0 or 1.0 == stock). Defaults match plot_coords_puzzler "
        "exactly. Returns (x, y).");
  m.def("plot_coords_puzzler_resolver_off", &plot_coords_puzzler_resolver_off, py::arg("structure"),
        "Lay out a dot-bracket structure with RNApuzzler, with "
        "checkSiblingIntersections/checkAncestorIntersections/optimize all "
        "false (checkExteriorIntersections and allowFlipping left at their "
        "stock defaults) -- the vendored-side counterpart of the native "
        "rna_layout::layout_puzzler resolver-off finalization path "
        "(Milestone A step 6). Returns (x, y).");
  m.def("plot_coords_puzzler_sibling_only", &plot_coords_puzzler_sibling_only, py::arg("structure"),
        "Lay out a dot-bracket structure with RNApuzzler, with "
        "checkSiblingIntersections true and checkAncestorIntersections/"
        "optimize both false (checkExteriorIntersections and allowFlipping "
        "left at their stock defaults) -- the vendored-side counterpart of "
        "the native rna_layout::layout_puzzler SIBLING-only resolver path "
        "(Milestone A step 7). Returns (x, y).");
  m.def("plot_coords_turtle", &plot_coords_turtle, py::arg("structure"),
        "Lay out a dot-bracket structure with RNAturtle; returns (x, y).");
  m.def("dump_turtle", &plot_coords_turtle, py::arg("structure"),
        "T0 turtle-coordinate dump for parity testing against the native "
        "rna_layout core's dump_turtle; identical to plot_coords_turtle "
        "(turtle has no separate pre-tree dump point -- it terminates "
        "before any tree is built).");

#ifdef RNA_DRAW_BUILD_ORACLE
  m.def("oracle_instrumentation_version", &rnadraw_oracle_instrumentation_version,
        "Marker string proving the RNA_DRAW_BUILD_ORACLE instrumentation TU "
        "(vendor_instrument.c) is compiled in; see that file for the "
        "mechanism dump_detections/dump_change_trace entry points will use.");
  m.def("dump_tree", &dump_tree, py::arg("structure"), py::arg("paired") = 35.0,
        py::arg("unpaired") = 25.0,
        "T1 config-tree/bounding-box dump (post updateBoundingBoxes, "
        "pre-resolver), as a JSON string -- json.loads() it and compare "
        "against the native rna_layout core's dump_tree (a list[dict] of "
        "the same shape); only built when RNA_DRAW_BUILD_ORACLE is on.");
  m.def("dump_detections", &dump_detections, py::arg("structure"), py::arg("paired") = 35.0,
        py::arg("unpaired") = 25.0,
        "The full intersection detection set over the T1 tree (Milestone A "
        "step 5), as a JSON string of [[node1_id, node2_id, \"type\"], "
        "...] -- json.loads() it and compare against the native "
        "rna_layout core's dump_detections; only built when "
        "RNA_DRAW_BUILD_ORACLE is on.");
  m.def("dump_change_trace", &dump_change_trace, py::arg("structure"), py::arg("paired") = 35.0,
        py::arg("unpaired") = 25.0, py::arg("max_config_changes") = 25000,
        "The ordered SIBLING-resolver config-change trace (Milestone A step "
        "7, checkSiblingIntersections=1/checkAncestorIntersections=0/"
        "optimize=0), as a JSON string -- json.loads() it and compare "
        "against the native rna_layout core's dump_change_trace. WARNING: "
        "does not terminate on every structure (see "
        "benchmarks/hard_set_sibling_only_oracle_hangs.json); only built "
        "when RNA_DRAW_BUILD_ORACLE is on.");
#endif
}
