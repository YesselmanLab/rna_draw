// ctest: `rna_layout::layout_turtle` sanity checks (`turtle.hpp`) --
// output-shape/finiteness/determinism properties (parity criterion 5),
// plus the malformed-input guards. Coordinate parity against the vendored
// oracle (criterion 1) is checked at the Python level
// (`tests/test_native_parity.py`), where both compiled extensions are
// importable together.

#include "rna_layout/turtle.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

using rna_layout::Coords;
using rna_layout::layout_turtle;

namespace {

int g_failures = 0;

void expect(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

bool all_finite(const Coords& coords) {
  for (std::size_t i = 0; i < coords.x.size(); ++i) {
    if (!std::isfinite(coords.x[i]) || !std::isfinite(coords.y[i])) {
      return false;
    }
  }
  return true;
}

void test_output_length_matches_structure() {
  const std::string structure = "((((....))))";
  const Coords coords = layout_turtle(structure);
  expect(coords.x.size() == structure.size(), "turtle output length == structure length (x)");
  expect(coords.y.size() == structure.size(), "turtle output length == structure length (y)");
}

void test_all_finite() {
  const std::string structure =
      "..(((((((((((......)).)))))))))...((....))................(((......((((....))"
      "))....)))..";
  expect(all_finite(layout_turtle(structure)), "turtle coordinates are finite (no NaN/Inf)");
}

void test_deterministic() {
  const std::string structure = "(((...(((...)))...(((...)))...)))";
  const Coords first = layout_turtle(structure);
  const Coords second = layout_turtle(structure);
  expect(first.x == second.x && first.y == second.y,
         "turtle is deterministic: two calls on the same structure agree exactly");
}

void test_various_motifs_do_not_throw() {
  // One hand structure per motif named in the parity harness corpus
  // (`.claude/plans/current-plan.md`): hairpin, bulge (single unpaired,
  // both strands), internal loop, multiloop, stacked helices, exterior
  // dangles, multi-branch exterior.
  const char* structures[] = {
      "((((....))))",                  // hairpin
      "(.((....)))",                   // bulge, unpaired on the near strand
      "(((....)).)",                   // bulge, unpaired on the far strand
      "((((..((((....))))..))))",      // internal loop
      "((((...)))(((...))))",          // multiloop (two-way)
      "((((((((....))))))))",          // stacked helices
      "..((((....))))..",              // exterior dangles both sides
      "((....))..((....))..((....))",  // multi-branch exterior
  };
  for (const char* structure : structures) {
    try {
      const Coords coords = layout_turtle(std::string(structure));
      expect(coords.x.size() == std::string(structure).size(),
             "motif structure output length matches input length");
      expect(all_finite(coords), "motif structure coordinates are finite");
    } catch (const std::exception& e) {
      std::cerr << "FAIL structure " << structure << " threw: " << e.what() << "\n";
      ++g_failures;
    }
  }
}

void test_malformed_input_throws() {
  bool threw = false;
  try {
    (void)layout_turtle(std::string("(("));  // unbalanced
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  expect(threw, "unbalanced structure throws std::invalid_argument");

  threw = false;
  try {
    (void)layout_turtle(std::string(""));  // empty
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  expect(threw, "empty structure throws std::invalid_argument");

  threw = false;
  try {
    (void)layout_turtle(std::string("().()"));  // bare empty loop
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  expect(threw, "empty-loop structure throws std::invalid_argument");
}

}  // namespace

int main() {
  test_output_length_matches_structure();
  test_all_finite();
  test_deterministic();
  test_various_motifs_do_not_throw();
  test_malformed_input_throws();

  if (g_failures > 0) {
    std::cerr << "turtle_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
