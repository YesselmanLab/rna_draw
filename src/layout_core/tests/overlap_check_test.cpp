// ctest: `rna_layout::overlap` primitives + checker (`overlap_check.hpp`),
// port of `rna_draw/overlap.py:36-468`. Covers the PLAN-CRITIC R1 PrimId
// string-order lock, exclusion rules, hand-computed overlaps, and the
// C++-only self-oracle `check_overlaps == check_overlaps_bruteforce`
// (mirrors `tests/test_overlap.py`'s `TestHashEqualsBruteforce`). Plain
// assert-based `main()`, matching `geometry_test.cpp`'s style.

#include "rna_layout/overlap_check.hpp"

#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

using namespace rna_layout::overlap;  // NOLINT(google-build-using-namespace)

namespace {

int g_failures = 0;

void expect_true(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

void expect_throws_invalid_argument(const std::function<void()>& fn, const char* what) {
  try {
    fn();
  } catch (const std::invalid_argument&) {
    return;
  } catch (...) {
    std::cerr << "FAIL " << what << ": threw the wrong exception type\n";
    ++g_failures;
    return;
  }
  std::cerr << "FAIL " << what << ": did not throw\n";
  ++g_failures;
}

// -----------------------------------------------------------------------
// PLAN-CRITIC R1: PrimKind's enum-declaration order must reproduce
// Python's kind-STRING sort order ("bb" < "nt" < "pair"), not build order
// or alphabetical member-name order. This locks the invariant so a future
// accidental reordering of the `PrimKind` enum fails loudly here instead
// of silently flipping id_a/id_b in the differential harness.
// -----------------------------------------------------------------------
void test_prim_kind_matches_python_string_order() {
  expect_true(kind_str(PrimKind::bb) < kind_str(PrimKind::nt), R"(kind_str: "bb" < "nt")");
  expect_true(kind_str(PrimKind::nt) < kind_str(PrimKind::pair), R"(kind_str: "nt" < "pair")");
  expect_true(PrimKind::bb < PrimKind::nt, "PrimKind::bb < PrimKind::nt (matches string order)");
  expect_true(PrimKind::nt < PrimKind::pair,
              "PrimKind::nt < PrimKind::pair (matches string order)");

  const PrimId bb{PrimKind::bb, 5};
  const PrimId nt{PrimKind::nt, 5};
  const PrimId pair{PrimKind::pair, 5};
  expect_true(bb < nt, "PrimId: bb-kind sorts before nt-kind at equal index");
  expect_true(nt < pair, "PrimId: nt-kind sorts before pair-kind at equal index");
}

// -----------------------------------------------------------------------
// is_excluded (overlap.py:219-241)
// -----------------------------------------------------------------------
void test_is_excluded_adjacent_disks() {
  const Disk d0{PrimId{PrimKind::nt, 0}, 0, 0, 10, Ends{0, -1}};
  const Disk d1{PrimId{PrimKind::nt, 1}, 5, 0, 10, Ends{1, -1}};
  const std::vector<int> pair_map = {-1, -1};
  expect_true(is_excluded(d0, d1, pair_map), "adjacent disks (|i-j|==1) are excluded");
}

void test_is_excluded_paired_disks() {
  const Disk d0{PrimId{PrimKind::nt, 0}, 0, 0, 10, Ends{0, -1}};
  const Disk d5{PrimId{PrimKind::nt, 5}, 100, 0, 10, Ends{5, -1}};
  const std::vector<int> pair_map = {5, -1, -1, -1, -1, 0};
  expect_true(is_excluded(d0, d5, pair_map), "base-paired disks are excluded");
}

void test_is_excluded_disk_on_its_own_capsule() {
  const Disk disk{PrimId{PrimKind::nt, 0}, 0, 0, 10, Ends{0, -1}};
  const Capsule capsule{PrimId{PrimKind::bb, 0}, 0, 0, 20, 0, 5, Ends{0, 1}};
  const std::vector<int> pair_map = {-1, -1};
  expect_true(is_excluded(disk, capsule, pair_map), "disk excluded from its own capsule");
}

void test_is_excluded_capsules_sharing_endpoint() {
  const Capsule c1{PrimId{PrimKind::bb, 0}, 0, 0, 20, 0, 5, Ends{0, 1}};
  const Capsule c2{PrimId{PrimKind::bb, 1}, 20, 0, 40, 0, 5, Ends{1, 2}};
  const std::vector<int> pair_map = {-1, -1, -1};
  expect_true(is_excluded(c1, c2, pair_map), "capsules sharing an endpoint are excluded");
}

void test_is_excluded_distant_disks_not_excluded() {
  const Disk d0{PrimId{PrimKind::nt, 0}, 0, 0, 10, Ends{0, -1}};
  const Disk d5{PrimId{PrimKind::nt, 5}, 500, 0, 10, Ends{5, -1}};
  const std::vector<int> pair_map = {-1, -1, -1, -1, -1, -1};
  expect_true(!is_excluded(d0, d5, pair_map), "distant, non-adjacent, unpaired disks not excluded");
}

// -----------------------------------------------------------------------
// Hand-computed overlaps (mirrors tests/test_overlap.py's TestKnownBad)
// -----------------------------------------------------------------------

const Witness* find_witness(const OverlapResult& result, PrimId id_a, PrimId id_b) {
  for (const Witness& witness : result.witnesses) {
    const bool forward = witness.id_a == id_a && witness.id_b == id_b;
    const bool backward = witness.id_a == id_b && witness.id_b == id_a;
    if (forward || backward) {
      return &witness;
    }
  }
  return nullptr;
}

void test_identical_coords_far_in_index_flagged() {
  const std::vector<double> x = {0.0, 30.0, 60.0, 0.0};
  const std::vector<double> y = {0.0, 0.0, 0.0, 0.0};
  const std::vector<int> pair_map = {-1, -1, -1, -1};
  const OverlapResult result = check_overlaps(x, y, pair_map, OverlapParams{});
  const Witness* witness = find_witness(result, PrimId{PrimKind::nt, 0}, PrimId{PrimKind::nt, 3});
  expect_true(witness != nullptr, "identical far-in-index disks flagged");
  if (witness != nullptr) {
    expect_true(witness->kind == OverlapKind::disk_disk, "disk-disk kind");
  }
}

void test_disk_on_foreign_pair_capsule_flagged() {
  const std::vector<double> x = {0.0, 50.0, 75.0, 100.0};
  const std::vector<double> y = {0.0, 0.0, 50.0, 0.0};
  const std::vector<int> pair_map = {3, -1, -1, 0};
  const OverlapResult result = check_overlaps(x, y, pair_map, OverlapParams{});
  const Witness* witness = find_witness(result, PrimId{PrimKind::nt, 1}, PrimId{PrimKind::pair, 0});
  expect_true(witness != nullptr, "disk on a foreign pair capsule flagged");
  if (witness != nullptr) {
    expect_true(witness->kind == OverlapKind::disk_capsule, "disk-capsule kind");
  }
}

// PLAN-CRITIC R1's adversarial case: a disk overlapping a NON-ADJACENT
// backbone capsule -- the "bb" vs "nt" ordering flip bug this test would
// have caught. Nucleotide 5 sits squarely on the (0,1) backbone segment,
// far from it in index (not 0, 1, or adjacent), so this is a genuine,
// non-excluded disk-vs-backbone overlap.
void test_disk_vs_nonadjacent_backbone_flagged_with_bb_ordered_first() {
  std::vector<double> x(6, 0.0);
  std::vector<double> y(6, 0.0);
  x[0] = 0.0;
  x[1] = 100.0;  // backbone capsule ("bb", 0) spans (0,0)-(100,0)
  x[2] = 200.0;
  x[3] = 300.0;
  x[4] = 400.0;
  x[5] = 50.0;  // nucleotide 5 sits on the (0,1) backbone segment's midpoint
  y[5] = 0.0;
  const std::vector<int> pair_map(6, -1);
  const OverlapParams params{10.0, 5.0, 5.0, 1e-6};
  const OverlapResult result = check_overlaps(x, y, pair_map, params);

  const PrimId bb0{PrimKind::bb, 0};
  const PrimId nt5{PrimKind::nt, 5};
  const Witness* witness = find_witness(result, nt5, bb0);
  expect_true(witness != nullptr, "disk vs non-adjacent backbone capsule flagged");
  if (witness != nullptr) {
    expect_true(witness->kind == OverlapKind::disk_capsule, "disk-vs-backbone kind");
    // R1: Python's sorted((a.pid,b.pid)) picks "bb" as id_a ("bb" < "nt"
    // lexicographically) -- a naive {nt,bb,pair} enum would flip this.
    expect_true(witness->id_a == bb0, R"(id_a is the "bb" primitive (string-order R1 fix))");
    expect_true(witness->id_b == nt5, R"(id_b is the "nt" primitive (string-order R1 fix))");
  }
}

void test_crossing_pair_capsules_flagged() {
  const std::vector<double> x = {0.0, 10.0, 10.0, 0.0};
  const std::vector<double> y = {0.0, 10.0, 0.0, 10.0};
  const std::vector<int> pair_map = {1, 0, 3, 2};
  const OverlapParams params{0.01, 0.01, 1.0, 1e-6};
  const OverlapResult result = check_overlaps(x, y, pair_map, params);
  const Witness* witness =
      find_witness(result, PrimId{PrimKind::pair, 0}, PrimId{PrimKind::pair, 2});
  expect_true(witness != nullptr, "crossing pair capsules flagged");
  if (witness != nullptr) {
    expect_true(witness->kind == OverlapKind::capsule_capsule, "capsule-capsule kind");
  }
}

void test_exactly_touching_disks_not_flagged() {
  const std::vector<double> x = {0.0, 100.0, 20.0};
  const std::vector<double> y = {0.0, 100.0, 0.0};
  const std::vector<int> pair_map = {-1, -1, -1};
  const OverlapParams params{10.0, 0.0, 0.0, 1e-6};
  const OverlapResult result = check_overlaps(x, y, pair_map, params);
  expect_true(result.passed(), "exactly-touching disks are not flagged");
}

// -----------------------------------------------------------------------
// counts_by_kind index pinning (R1: disk_disk->0, disk_capsule->1,
// capsule_capsule->2)
// -----------------------------------------------------------------------
void test_counts_by_kind_pinned_indices() {
  const std::vector<double> x = {0.0, 30.0, 60.0, 0.0};
  const std::vector<double> y = {0.0, 0.0, 0.0, 0.0};
  const std::vector<int> pair_map = {-1, -1, -1, -1};
  const OverlapResult result = check_overlaps(x, y, pair_map, OverlapParams{});
  expect_true(result.counts_by_kind[static_cast<std::size_t>(OverlapKind::disk_disk)] > 0,
              "counts_by_kind[disk_disk] holds the disk-disk tally");
  expect_true(static_cast<int>(result.counts_by_kind[0] + result.counts_by_kind[1] +
                               result.counts_by_kind[2]) == result.num_overlaps(),
              "counts_by_kind sums to num_overlaps");
}

// -----------------------------------------------------------------------
// Input validation (R2)
// -----------------------------------------------------------------------
void test_validate_inputs_raises() {
  expect_throws_invalid_argument(
      [] { (void)check_overlaps({0.0, 1.0}, {0.0}, {-1, -1}, OverlapParams{}); },
      "mismatched lengths raise std::invalid_argument");
  expect_throws_invalid_argument([] { (void)check_overlaps({}, {}, {}, OverlapParams{}); },
                                 "empty input raises std::invalid_argument");
  expect_throws_invalid_argument(
      [] { (void)check_overlaps({0.0, 1.0}, {0.0, 1.0}, {1, -1}, OverlapParams{}); },
      "asymmetric pair_map raises std::invalid_argument");
}

// -----------------------------------------------------------------------
// C++-only self-oracle: check_overlaps == check_overlaps_bruteforce
// (mirrors tests/test_overlap.py's TestHashEqualsBruteforce) -- isolates
// spatial-hash candidate-set bugs from predicate bugs.
// -----------------------------------------------------------------------

/// A tiny deterministic PRNG (xorshift32) so this ctest needs no external
/// RNG dependency and is reproducible across runs/platforms.
class Xorshift32 {
 public:
  explicit Xorshift32(std::uint32_t seed) : state_(seed == 0 ? 1 : seed) {}

  std::uint32_t next() {
    state_ ^= state_ << 13;
    state_ ^= state_ >> 17;
    state_ ^= state_ << 5;
    return state_;
  }

  double uniform(double lo, double hi) {
    const double frac = static_cast<double>(next()) / static_cast<double>(0xFFFFFFFFu);
    return lo + frac * (hi - lo);
  }

 private:
  std::uint32_t state_;
};

using WitnessKey = std::tuple<int, int, int, int, int>;

std::set<WitnessKey> witness_key_set(const OverlapResult& result) {
  std::set<WitnessKey> keys;
  for (const Witness& witness : result.witnesses) {
    keys.emplace(static_cast<int>(witness.kind), static_cast<int>(witness.id_a.kind),
                 witness.id_a.index, static_cast<int>(witness.id_b.kind), witness.id_b.index);
  }
  return keys;
}

void random_case(Xorshift32& rng, int n, std::vector<double>& x, std::vector<double>& y,
                 std::vector<int>& pair_map) {
  const double box = 40.0;  // packed tightly relative to the default node_r=10
  x.resize(static_cast<std::size_t>(n));
  y.resize(static_cast<std::size_t>(n));
  pair_map.assign(static_cast<std::size_t>(n), -1);
  for (int i = 0; i < n; ++i) {
    x[static_cast<std::size_t>(i)] = rng.uniform(0.0, box);
    y[static_cast<std::size_t>(i)] = rng.uniform(0.0, box);
  }
  std::vector<int> shuffled(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    shuffled[static_cast<std::size_t>(i)] = i;
  }
  for (int i = n - 1; i > 0; --i) {
    const int j = static_cast<int>(rng.next() % static_cast<std::uint32_t>(i + 1));
    std::swap(shuffled[static_cast<std::size_t>(i)], shuffled[static_cast<std::size_t>(j)]);
  }
  const int num_pairs =
      n / 2 > 0 ? static_cast<int>(rng.next() % static_cast<std::uint32_t>(n / 2 + 1)) : 0;
  for (int k = 0; k < num_pairs; ++k) {
    const int i = shuffled[static_cast<std::size_t>(2 * k)];
    const int j = shuffled[static_cast<std::size_t>(2 * k + 1)];
    pair_map[static_cast<std::size_t>(i)] = j;
    pair_map[static_cast<std::size_t>(j)] = i;
  }
}

void test_hash_equals_bruteforce_on_random_layouts() {
  Xorshift32 rng(12345);
  for (const int n : {5, 20, 50, 120}) {
    for (int seed = 0; seed < 5; ++seed) {
      std::vector<double> x;
      std::vector<double> y;
      std::vector<int> pair_map;
      random_case(rng, n, x, y, pair_map);
      const OverlapResult hashed = check_overlaps(x, y, pair_map, OverlapParams{});
      const OverlapResult brute = check_overlaps_bruteforce(x, y, pair_map, OverlapParams{});
      expect_true(witness_key_set(hashed) == witness_key_set(brute),
                  "check_overlaps == check_overlaps_bruteforce on a random layout");
      expect_true(hashed.num_overlaps() == brute.num_overlaps(),
                  "check_overlaps and check_overlaps_bruteforce agree on num_overlaps");
    }
  }
}

void test_hash_equals_bruteforce_on_long_diagonal_pair_capsule() {
  // Direct regression for the spatial-hash tiling path: nucleotides 0/1
  // are paired far apart on a near-diagonal (insert_segment tiling
  // territory); nucleotide 2 sits on that diagonal's midpoint.
  const std::vector<double> x = {0.0, 950.0, 475.0};
  const std::vector<double> y = {0.0, 940.0, 470.0};
  const std::vector<int> pair_map = {1, 0, -1};
  const OverlapParams params{10.0, 1.0, 5.0, 1e-6};
  const OverlapResult hashed = check_overlaps(x, y, pair_map, params);
  const OverlapResult brute = check_overlaps_bruteforce(x, y, pair_map, params);
  expect_true(witness_key_set(hashed) == witness_key_set(brute),
              "hash == bruteforce across a long diagonal pair capsule");
  expect_true(hashed.num_overlaps() > 0, "the long-diagonal case has a genuine overlap");
}

}  // namespace

int main() {
  test_prim_kind_matches_python_string_order();
  test_is_excluded_adjacent_disks();
  test_is_excluded_paired_disks();
  test_is_excluded_disk_on_its_own_capsule();
  test_is_excluded_capsules_sharing_endpoint();
  test_is_excluded_distant_disks_not_excluded();
  test_identical_coords_far_in_index_flagged();
  test_disk_on_foreign_pair_capsule_flagged();
  test_disk_vs_nonadjacent_backbone_flagged_with_bb_ordered_first();
  test_crossing_pair_capsules_flagged();
  test_exactly_touching_disks_not_flagged();
  test_counts_by_kind_pinned_indices();
  test_validate_inputs_raises();
  test_hash_equals_bruteforce_on_random_layouts();
  test_hash_equals_bruteforce_on_long_diagonal_pair_capsule();

  if (g_failures > 0) {
    std::cerr << g_failures << " failure(s)\n";
    return 1;
  }
  std::cout << "overlap_check_test: all checks passed\n";
  return 0;
}
