// pybind11 bindings for ViennaRNA's in-process secondary-structure layout
// algorithms (puzzler, naview, turtle), exposed to Python as
// `rna_draw._vienna_layout`.
//
// REENTRANCY: RNApuzzler and RNAturtle are reentrant (no mutable file-scope
// state). naview is NOT -- `naview.o` in libRNA.a reads/writes file-scope
// mutable BSS globals (`_bases`, `_nbase`, `_loops`, `_loop_count`,
// `_regions`, `_root`, `_lencut`). Two concurrent calls to
// `vrna_plot_coords_naview` in one address space corrupt each other and
// silently return garbage coordinates. Callers MUST NOT invoke the naview
// binding from more than one thread at a time in this process; a benchmark
// harness that parallelizes naview must use process-based parallelism
// (`ProcessPoolExecutor`/`multiprocessing`), never a thread pool.

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
// not match the plain C symbols in libRNA.a -- wrap the include ourselves.
//
// vrna_config.h must be included before plotting/layouts.h: layouts.h
// gates its naview include behind `#ifdef VRNA_WITH_NAVIEW_LAYOUT`, which
// vrna_config.h is what defines.
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
/// Verified (M5.1 spike): `vrna_plot_coords_naview("")` ABORTS (SIGABRT,
/// not a clean 0 return) -- puzzler/turtle already return 0 cleanly for
/// this case (caught by the `n == 0` check after the call), but naview
/// does not, so every binding rejects it uniformly up front instead of
/// relying on that post-call check.
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
/// input (naview alone tolerates it). This is the same degenerate shape
/// that already crashes `render_rna.py`'s tree recursion (see
/// `rna_draw.layout.base.has_empty_loop`); guard uniformly across all
/// three engines rather than special-case per algorithm.
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

/// Call `vrna_plot_coords_naview(structure, &x, &y)`.
///
/// NOT REENTRANT (see file header) -- caller must not invoke this from more
/// than one thread at a time in this process. Unlike puzzler/turtle,
/// naview tolerates a bare empty loop `"()"` (verified: returns sane
/// coordinates), so `validate_no_empty_loop` is deliberately NOT applied
/// here -- it would needlessly reject input naview can actually handle.
CoordVectors plot_coords_naview(const std::string& structure) {
  validate_nonempty(structure);
  validate_well_nested(structure);
  MallocBuffer<float> x;
  MallocBuffer<float> y;
  int n = vrna_plot_coords_naview(structure.c_str(), &x.ptr, &y.ptr);
  if (n == 0 || static_cast<size_t>(n) != structure.size()) {
    throw std::runtime_error("vrna_plot_coords_naview failed on structure of length " +
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

/// Call `vrna_plot_coords_puzzler(structure, &x, &y, NULL, options)` with an
/// explicit options struct exposing the two resolver levers:
/// `allow_flipping` (RNApuzzler's exterior-branch flip heuristic) and
/// `max_config_changes` (the config-change search budget; <= 0 falls back
/// to the engine's 25000 default -- see the vendored `RNApuzzler.c` edit
/// that makes this budget caller-respectable instead of hardcoded).
CoordVectors plot_coords_puzzler_opts(const std::string& structure, bool allow_flipping,
                                      int max_config_changes) {
  validate_nonempty(structure);
  validate_well_nested(structure);
  validate_no_empty_loop(structure);
  PuzzlerOptions options;
  options.ptr->allowFlipping = allow_flipping ? 1 : 0;
  options.ptr->maximumNumberOfConfigChangesAllowed = max_config_changes;
  MallocBuffer<float> x;
  MallocBuffer<float> y;
  int n = vrna_plot_coords_puzzler(structure.c_str(), &x.ptr, &y.ptr, nullptr, options.ptr);
  if (n == 0 || static_cast<size_t>(n) != structure.size()) {
    throw std::runtime_error("vrna_plot_coords_puzzler_opts failed on structure of length " +
                             std::to_string(structure.size()));
  }
  return to_vectors(x, y, n);
}

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
      "In-process bindings to ViennaRNA's puzzler/naview/turtle layout "
      "algorithms. naview is NOT reentrant (file-scope globals in "
      "naview.o) -- callers must serialize naview calls within a process "
      "and use process-based (not thread-based) parallelism across "
      "structures.";

  m.def("version", &version, "Compiled-against ViennaRNA version string.");
  m.def("abi_version", &abi_version,
        "Compiled-against ViennaRNA version as (major, minor, patch).");
  m.def("sizeof_puzzler_options", &sizeof_puzzler_options,
        "sizeof(vrna_plot_options_puzzler_t), for ABI-drift detection.");

  m.def("plot_coords_puzzler", &plot_coords_puzzler, py::arg("structure"),
        "Lay out a dot-bracket structure with RNApuzzler; returns (x, y).");
  m.def("plot_coords_puzzler_opts", &plot_coords_puzzler_opts, py::arg("structure"),
        py::arg("allow_flipping") = false, py::arg("max_config_changes") = 0,
        "Lay out a dot-bracket structure with RNApuzzler, exposing two "
        "resolver levers: allow_flipping (exterior-branch flip heuristic) "
        "and max_config_changes (config-change search budget; <= 0 uses "
        "the engine's default). Defaults match plot_coords_puzzler exactly. "
        "Returns (x, y).");
  m.def("plot_coords_naview", &plot_coords_naview, py::arg("structure"),
        "Lay out a dot-bracket structure with naview; returns (x, y). "
        "NOT reentrant -- see module docstring.");
  m.def("plot_coords_turtle", &plot_coords_turtle, py::arg("structure"),
        "Lay out a dot-bracket structure with RNAturtle; returns (x, y).");
}
