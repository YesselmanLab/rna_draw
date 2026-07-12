// ctest: `rna_layout::layout_puzzler` (`puzzler.hpp`), Milestone A step 6
// ("finalization with resolver OFF"). Parity against the vendored oracle is
// a Python-level concern (`tests/test_native_parity.py`); this proves the
// native resolver-off pipeline runs, guards its resolver-on precondition,
// and produces sane, deterministic output.

#include "rna_layout/puzzler.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

using rna_layout::Coords;
using rna_layout::layout_puzzler;
using rna_layout::PuzzlerOptions;

namespace {

int g_failures = 0;

void expect(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

PuzzlerOptions resolver_off_options() {
  PuzzlerOptions opts;
  opts.check_sibling = false;
  opts.check_ancestor = false;
  opts.optimize = false;
  // `check_exterior` deliberately left at its default (`true`): see
  // `puzzler.cpp`'s `run_config_tree_pipeline` doc comment -- this keeps
  // `update_bounding_boxes` in the pipeline (matching the vendored default
  // options), while `checkAndFixIntersections` is still provably a no-op.
  return opts;
}

void test_resolver_on_options_throw() {
  bool threw = false;
  try {
    const Coords unused =
        layout_puzzler(std::string("((((....))))"), PuzzlerOptions{});  // defaults are resolver-ON
    (void)unused;
  } catch (const std::logic_error&) {
    threw = true;
  }
  expect(threw, "layout_puzzler: resolver-on options (the default PuzzlerOptions) throw");
}

void test_each_resolver_flag_individually_throws() {
  for (bool sibling : {true, false}) {
    for (bool ancestor : {true, false}) {
      for (bool optimize : {true, false}) {
        if (!sibling && !ancestor && !optimize) {
          continue;  // the one allowed combination
        }
        PuzzlerOptions opts = resolver_off_options();
        opts.check_sibling = sibling;
        opts.check_ancestor = ancestor;
        opts.optimize = optimize;

        bool threw = false;
        try {
          const Coords unused = layout_puzzler(std::string("((((....))))"), opts);
          (void)unused;
        } catch (const std::logic_error&) {
          threw = true;
        }
        expect(threw, "layout_puzzler: any single resolver flag set throws");
      }
    }
  }
}

void test_resolver_off_produces_sane_coordinates() {
  const std::string structure = "((((...)))(((...)))(((...))))";
  const Coords coords = layout_puzzler(structure, resolver_off_options());

  expect(coords.x.size() == structure.size(), "layout_puzzler: output length matches input length");
  expect(coords.y.size() == structure.size(),
         "layout_puzzler: output length matches input length (y)");

  for (std::size_t i = 0; i < coords.x.size(); ++i) {
    expect(std::isfinite(coords.x[i]), "layout_puzzler: every x coordinate is finite");
    expect(std::isfinite(coords.y[i]), "layout_puzzler: every y coordinate is finite");
  }
}

void test_resolver_off_is_deterministic() {
  const std::string structure = "((((..((((....))))..))))";
  const Coords first = layout_puzzler(structure, resolver_off_options());
  const Coords second = layout_puzzler(structure, resolver_off_options());
  expect(first.x == second.x, "layout_puzzler: repeated calls produce identical x coordinates");
  expect(first.y == second.y, "layout_puzzler: repeated calls produce identical y coordinates");
}

void test_malformed_input_throws() {
  bool threw = false;
  try {
    const Coords unused = layout_puzzler(std::string(""), resolver_off_options());
    (void)unused;
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  expect(threw, "layout_puzzler: empty structure throws invalid_argument");

  threw = false;
  try {
    const Coords unused = layout_puzzler(std::string("()"), resolver_off_options());
    (void)unused;
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  expect(threw, "layout_puzzler: empty-loop structure throws invalid_argument");
}

}  // namespace

int main() {
  test_resolver_on_options_throw();
  test_each_resolver_flag_individually_throws();
  test_resolver_off_produces_sane_coordinates();
  test_resolver_off_is_deterministic();
  test_malformed_input_throws();

  if (g_failures > 0) {
    std::cerr << "puzzler_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
