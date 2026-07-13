#pragma once

/**
 * @file overlap_check.hpp
 * @brief Value types + the overlap checker itself. Port of
 *        `rna_draw/overlap.py:36-468`.
 *
 * `rna_draw/overlap.py` (together with `geometry.py`/`spatial_hash.py`) is
 * the frozen, definitional arbiter of "does this layout overlap" -- this
 * file's ONLY job is to agree with it exactly (`.claude/plans/
 * current-plan-checker.md`). Preserves the primitive build order
 * (disks -> backbone -> pairs, for index stability), the exclusion rules,
 * and `separation = required - depth` (never a `dist` shortcut).
 *
 * PLAN-CRITIC R1 (BLOCKING) -- PrimId ordering: Python's `PrimitiveId` is
 * `@dataclass(order=True)` over `(kind: str, index: int)`, so
 * `sorted((a.pid, b.pid))` (`overlap.py:294`) orders by the KIND STRING:
 * `"bb" < "nt" < "pair"` (lexicographic: 'b' < 'n' < 'p'). A naive C++
 * `enum class PrimKind {nt, bb, pair}` with the default integer
 * `operator<` would give `nt(0) < bb(1)`, which FLIPS `id_a`/`id_b` for
 * every disk-vs-backbone overlap relative to Python. FIX (this file):
 * `PrimKind` is declared in the SAME ORDER as the kind strings sort --
 * `{bb, nt, pair}` -- so the compiler-generated integer `operator<`
 * reproduces the Python string order exactly, with no separate string
 * comparison needed. `kind_str(PrimKind)` below is the inverse mapping
 * back to those same three strings, used only for display/binding
 * (`bindings.cpp`, the differential harness) -- it is NOT consulted for
 * ordering (ordering is the enum value itself, by construction).
 * `overlap_check_test.cpp` asserts this equivalence directly (locks the
 * invariant against a future accidental enum reordering).
 *
 * `OverlapKind`'s integer values are PINNED to match `overlap.py`'s
 * `OverlapKind` member declaration order -- `disk_disk=0, disk_capsule=1,
 * capsule_capsule=2` -- so `OverlapResult::counts_by_kind`'s array index
 * is unambiguous against the differential harness's Python-side
 * `counts_by_kind` dict (keyed by the same three kinds).
 */

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rna_layout/spatial_hash.hpp"
#include "rna_layout/types.hpp"

namespace rna_layout::overlap {

/// Primitive family. Declared in STRING-SORT order (`"bb" < "nt" <
/// "pair"`), not alphabetical-by-full-name or build order -- see this
/// file's header (R1). `std::uint8_t`-backed (mirrors `types.hpp`'s
/// `BaseType`): only 3 values are ever needed.
enum class PrimKind : std::uint8_t { bb, nt, pair };

/// The three-letter/four-letter kind tag Python's `PrimitiveId.kind`
/// stores, e.g. for the differential harness's witness-key tuples.
[[nodiscard]] std::string kind_str(PrimKind kind);

/**
 * Stable identity for a drawn primitive. Mirrors `PrimitiveId`
 * (`geometry.py:16-31`); orderable by `(kind, index)`, where `kind`'s
 * enum-value order reproduces Python's kind-STRING order (see this
 * file's header).
 */
struct PrimId {
  PrimKind kind = PrimKind::nt;
  int index = 0;
};

[[nodiscard]] bool operator==(const PrimId& lhs, const PrimId& rhs);
[[nodiscard]] bool operator<(const PrimId& lhs, const PrimId& rhs);

/**
 * The two (or one, for a disk) nucleotide indices a primitive touches,
 * replacing Python's `frozenset[int]` (`Disk.ends`/`Capsule.ends`).
 * Allocation-free: a disk sets only @ref a (singleton `{i}`); @ref b
 * stays the `-1` sentinel (never a legitimate nucleotide index).
 */
struct Ends {
  int a = -1;
  int b = -1;

  /// Whether @p index is one of this primitive's end(s).
  [[nodiscard]] bool contains(int index) const;

  /// Whether this primitive shares at least one nucleotide index with
  /// @p other -- mirrors Python's `bool(a.ends & b.ends)`.
  [[nodiscard]] bool intersects(const Ends& other) const;
};

/// A nucleotide glyph: a circle. Mirrors `Disk` (`geometry.py:34-61`).
struct Disk {
  PrimId pid;
  double cx = 0.0;
  double cy = 0.0;
  double radius = 0.0;
  Ends ends;
};

/// A backbone or base-pair connector. Mirrors `Capsule` (`geometry.py:64-96`).
struct Capsule {
  PrimId pid;
  double x0 = 0.0;
  double y0 = 0.0;
  double x1 = 0.0;
  double y1 = 0.0;
  double half_width = 0.0;
  Ends ends;
};

/// Either drawn shape. Mirrors Python's `Primitive = Union[Disk, Capsule]`.
using Primitive = std::variant<Disk, Capsule>;

/// This primitive's identity, regardless of which shape it is.
[[nodiscard]] PrimId pid_of(const Primitive& primitive);

/// This primitive's nucleotide-index endpoint set, regardless of shape.
[[nodiscard]] Ends ends_of(const Primitive& primitive);

/// This primitive's axis-aligned bounding box. Mirrors `Disk.aabb`
/// (`geometry.py:54-61`) / `Capsule.aabb` (`geometry.py:88-96`).
[[nodiscard]] Aabb aabb_of(const Primitive& primitive);

/// The three kinds of primitive-pair overlap. Values PINNED to match
/// `overlap.py`'s `OverlapKind` declaration order (see this file's
/// header) -- `counts_by_kind` indexes on these integer values directly.
enum class OverlapKind : std::uint8_t { disk_disk = 0, disk_capsule = 1, capsule_capsule = 2 };

/// `OverlapKind`'s Python-side string value (`overlap.py:36-41`), e.g.
/// for the differential harness's witness-key tuples.
[[nodiscard]] std::string kind_str(OverlapKind kind);

/// Evidence of a single overlapping primitive pair. Mirrors `Witness`
/// (`overlap.py:44-62`).
struct Witness {
  OverlapKind kind = OverlapKind::disk_disk;
  PrimId id_a;
  PrimId id_b;
  double separation = 0.0;
  double overlap_depth = 0.0;
};

/// Geometry knobs the checker judges a layout against. Mirrors
/// `OverlapParams` (`overlap.py:89-117`), same defaults.
struct OverlapParams {
  double node_r = 10.0;
  double backbone_half_width = 7.5;
  double pair_half_width = 7.5;
  double tol = 1e-6;
};

/// Result of an overlap check. Mirrors `OverlapReport` (`overlap.py:65-86`).
struct OverlapResult {
  std::vector<Witness> witnesses;
  /// Indexed by `OverlapKind`'s pinned integer value (see this file's
  /// header): `[disk_disk, disk_capsule, capsule_capsule]`.
  std::array<int, 3> counts_by_kind{0, 0, 0};

  [[nodiscard]] int num_overlaps() const;
  [[nodiscard]] bool passed() const;
};

/// One disk per nucleotide. Mirrors `build_disks` (`overlap.py:120-134`).
[[nodiscard]] std::vector<Disk> build_disks(const std::vector<double>& x,
                                            const std::vector<double>& y, double node_r);

/// One capsule per consecutive backbone connection. Mirrors
/// `build_backbone_capsules` (`overlap.py:137-161`).
[[nodiscard]] std::vector<Capsule> build_backbone_capsules(const std::vector<double>& x,
                                                           const std::vector<double>& y,
                                                           double half_width);

/// One capsule per base pair `(i, j)`, `i < j`. Mirrors
/// `build_pair_capsules` (`overlap.py:164-191`).
[[nodiscard]] std::vector<Capsule> build_pair_capsules(const std::vector<double>& x,
                                                       const std::vector<double>& y,
                                                       const std::vector<int>& pair_map,
                                                       double half_width);

/// Every primitive a renderer draws: disks, then backbone capsules, then
/// pair capsules (index stability). Mirrors `build_primitives`
/// (`overlap.py:194-216`).
[[nodiscard]] std::vector<Primitive> build_primitives(const std::vector<double>& x,
                                                      const std::vector<double>& y,
                                                      const std::vector<int>& pair_map,
                                                      const OverlapParams& params);

/// Whether a primitive pair is *supposed* to touch (backbone-adjacent
/// disks, paired disks, or any pair sharing a nucleotide endpoint).
/// Mirrors `is_excluded` (`overlap.py:219-241`).
[[nodiscard]] bool is_excluded(const Primitive& a, const Primitive& b,
                               const std::vector<int>& pair_map);

/// Test one primitive pair for overlap. Mirrors `test_pair`
/// (`overlap.py:260-297`), split into three same-kind helpers
/// (`overlap_check.cpp`) per the project's <=30-line-function style note.
[[nodiscard]] std::optional<Witness> test_pair(const Primitive& a, const Primitive& b, double tol);

/**
 * Check a rendered layout for overlapping primitives via the spatial
 * hash. Mirrors `check_overlaps` (`overlap.py:409-437`).
 *
 * @throws std::invalid_argument If `x`/`y`/`pair_map` lengths mismatch,
 *     are empty, or `pair_map` is not symmetric (mirrors `_validate_inputs`,
 *     `overlap.py:300-319` -- PLAN-CRITIC R2: this entry point raises
 *     rather than leaving a degenerate call UB-adjacent; pybind11 maps
 *     this to a Python `ValueError`, matching Python's own contract).
 */
[[nodiscard]] OverlapResult check_overlaps(const std::vector<double>& x,
                                           const std::vector<double>& y,
                                           const std::vector<int>& pair_map,
                                           const OverlapParams& params);

/**
 * Reference/oracle overlap check: O(n^2), no spatial hash. Mirrors
 * `check_overlaps_bruteforce` (`overlap.py:440-468`) -- the C++ self-oracle
 * `overlap_check_test.cpp` uses to prove the spatial-hash path drops no
 * true overlap. Same input-validation contract as `check_overlaps`.
 */
[[nodiscard]] OverlapResult check_overlaps_bruteforce(const std::vector<double>& x,
                                                      const std::vector<double>& y,
                                                      const std::vector<int>& pair_map,
                                                      const OverlapParams& params);

}  // namespace rna_layout::overlap
