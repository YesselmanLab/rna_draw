"""Differential parity harness: the native C++ overlap-checker twin
(`rna_draw._layout_core.check_overlaps_count`/`check_overlaps_report`) vs
the frozen Python arbiter (`rna_draw.overlap.check_overlaps`).

`rna_draw/overlap.py` (+ `geometry.py` + `spatial_hash.py`) is the
DEFINITIONAL ARBITER and is NEVER modified by this module or by the C++
port it exercises. This is the load-bearing acceptance gate for that port
(`.claude/plans/current-plan-checker.md`): for every input, the Python and
C++ checkers must agree EXACTLY --

1. pass/fail EXACT: `py_report.passed == (cpp_count == 0)`.
2. count EXACT: total AND per-kind (`counts_by_kind`).
3. witness key-set EQUAL: `(kind, id_a.kind, id_a.index, id_b.kind,
   id_b.index)` tuples, as an unordered SET (no missing, no extra).
4. witness coords within `WITNESS_TOL = 1e-9` for every shared key
   (`separation`/`overlap_depth`) -- well above `std::hypot` vs
   `math.hypot` ULP noise (~1e-13 at these coordinate magnitudes), well
   below the `tol=1e-6` decision margin (see `overlap_geometry.hpp`'s file
   header).

Interpretation policy: a mismatch in any of the four checks is a PORT BUG
(a fidelity slip: reordered ops, `sqrt` vs `hypot`, a `separation=dist`
shortcut) -- fix the C++, not this test. The ONE exception is a genuine
cross-libm `hypot`/`sqrt` last-ULP tie surviving at the EXACT
`required - tol` boundary after verified expression-for-expression
fidelity; such a case is logged in `KNOWN_TIES` (case name -> ULP
argument) and SKIPPED here rather than loosened -- `rna_draw.overlap`
stays the arbiter for that class (the STOP criterion). `KNOWN_TIES` is
empty: no case in this port's corpora required one.

PLAN-CRITIC R3: "exactly required" adversarial center distances use
axis-aligned offsets (`dx = required, dy = 0`) so `hypot(dx, 0) == |dx|`
is bit-identical across libms -- a diagonal offset would make the "exact
boundary" itself FP-fuzzy.
"""

from __future__ import annotations

import os
import random
from pathlib import Path

import pytest

from benchmarks.hard_gate import CORPUS, parse_dbn
from rna_draw.layout.base import is_pseudoknot_free
from rna_draw.layout.fallback import SafeFallbackEngine
from rna_draw.overlap import OverlapKind, OverlapParams, check_overlaps
from rna_draw.overlap_native import _native_available
from rna_draw.render_rna import RNARenderer, get_pairmap_from_secstruct

requires_native = pytest.mark.skipif(
    not _native_available(), reason="rna_draw._layout_core not importable"
)

HARD_SET_JSON = Path(__file__).parent.parent / "benchmarks" / "hard_set.json"
WORST_SET_JSON = Path(__file__).parent.parent / "benchmarks" / "worst_set.json"

WITNESS_TOL = 1e-9

# Empty by design (the goal, per the plan's STOP criterion) -- see this
# module's docstring. Keyed by case name -> a one-line ULP argument; a
# populated entry causes `assert_parity` to SKIP (not silently pass) that
# one case, so a future addition here is always visible in test output.
KNOWN_TIES: dict[str, str] = {}


# ---------------------------------------------------------------------
# The core differential assertion
# ---------------------------------------------------------------------


def _cpp_report(
    x: list[float], y: list[float], pair_map: list[int], params: OverlapParams
) -> list[tuple]:
    """Raw witness rows from the native C++ twin's full report entry point."""
    from rna_draw import _layout_core

    return _layout_core.check_overlaps_report(
        list(x),
        list(y),
        list(pair_map),
        params.node_r,
        params.backbone_half_width,
        params.pair_half_width,
        params.tol,
    )


def _py_key(kind: str, kind_a: str, idx_a: int, kind_b: str, idx_b: int) -> tuple:
    return (kind, kind_a, idx_a, kind_b, idx_b)


def assert_parity(
    x: list[float],
    y: list[float],
    pair_map: list[int],
    params: OverlapParams | None = None,
    case_name: str = "case",
) -> None:
    """Assert the native C++ twin agrees EXACTLY with the Python arbiter.

    See this module's docstring for the four checks and the `KNOWN_TIES`
    escape hatch.
    """
    if case_name in KNOWN_TIES:
        pytest.skip(f"{case_name}: documented float-tie ({KNOWN_TIES[case_name]})")

    params = params or OverlapParams()
    py_report = check_overlaps(x, y, pair_map, params)
    cpp_rows = _cpp_report(x, y, pair_map, params)
    cpp_count = len(cpp_rows)

    # 1. pass/fail EXACT.
    assert py_report.passed == (cpp_count == 0), (
        f"{case_name}: pass/fail mismatch (py passed={py_report.passed}, cpp_count={cpp_count})"
    )

    # 2. count EXACT (total + per-kind, pinned index order disk_disk,
    # disk_capsule, capsule_capsule -- matches overlap.py's OverlapKind
    # declaration order, which is also this dict's iteration order).
    assert py_report.num_overlaps == cpp_count, (
        f"{case_name}: total count mismatch (py={py_report.num_overlaps}, cpp={cpp_count})"
    )
    py_counts = [
        py_report.counts_by_kind[OverlapKind.DISK_DISK],
        py_report.counts_by_kind[OverlapKind.DISK_CAPSULE],
        py_report.counts_by_kind[OverlapKind.CAPSULE_CAPSULE],
    ]
    kind_index = {"disk_disk": 0, "disk_capsule": 1, "capsule_capsule": 2}
    cpp_counts = [0, 0, 0]
    for row in cpp_rows:
        cpp_counts[kind_index[row[0]]] += 1
    assert py_counts == cpp_counts, (
        f"{case_name}: per-kind counts mismatch: py={py_counts} cpp={cpp_counts}"
    )

    # 3. witness key-set EQUAL.
    py_keys = {
        _py_key(w.kind.value, w.id_a.kind, w.id_a.index, w.id_b.kind, w.id_b.index)
        for w in py_report.witnesses
    }
    cpp_keys = {_py_key(row[0], row[1], row[2], row[3], row[4]) for row in cpp_rows}
    assert py_keys == cpp_keys, (
        f"{case_name}: witness key-set mismatch: py-only={py_keys - cpp_keys} "
        f"cpp-only={cpp_keys - py_keys}"
    )

    # 4. witness coords within WITNESS_TOL for every shared key.
    py_by_key = {
        _py_key(w.kind.value, w.id_a.kind, w.id_a.index, w.id_b.kind, w.id_b.index): w
        for w in py_report.witnesses
    }
    cpp_by_key = {_py_key(row[0], row[1], row[2], row[3], row[4]): row for row in cpp_rows}
    for key in py_keys:
        w = py_by_key[key]
        row = cpp_by_key[key]
        sep_diff = abs(w.separation - row[5])
        depth_diff = abs(w.overlap_depth - row[6])
        assert sep_diff <= WITNESS_TOL, f"{case_name}: separation diff {sep_diff} at {key}"
        assert depth_diff <= WITNESS_TOL, f"{case_name}: overlap_depth diff {depth_diff} at {key}"


# ---------------------------------------------------------------------
# (a) Production corpus / hard_set + a size-stratified corpus sample
# ---------------------------------------------------------------------

_BUCKETS = (("<300", 0, 300), ("300-1200", 300, 1200), (">1200", 1200, 10**9))
# Deliberately small default (this module runs inside the standard pytest
# suite, which must stay fast): a "full" pre-merge soak sets
# RNA_DRAW_PARITY_SAMPLE_N to a few hundred, per the plan's own knob.
_SAMPLE_N_PER_BUCKET = int(os.environ.get("RNA_DRAW_PARITY_SAMPLE_N", "3"))


def _corpus_sample(n_per_bucket: int) -> list[str]:
    """Deterministic size-stratified sample of pseudoknot-free structures."""
    if not CORPUS.exists():
        return []
    filled: dict[str, list[str]] = {b[0]: [] for b in _BUCKETS}
    for path in sorted(CORPUS.glob("*.dbn"), key=lambda p: p.name):
        parsed = parse_dbn(path)
        if parsed is None:
            continue
        _seq, struct = parsed
        if not is_pseudoknot_free(struct):
            continue
        for name, lo, hi in _BUCKETS:
            if lo <= len(struct) < hi and len(filled[name]) < n_per_bucket:
                filled[name].append(struct)
                break
        if all(len(v) >= n_per_bucket for v in filled.values()):
            break
    return [s for bucket in filled.values() for s in bucket]


@requires_native
class TestProductionCorpusParity:
    """Section (a): the real `layout_guaranteed` pipeline at its own node_r."""

    def test_layout_guaranteed_on_corpus_sample(self) -> None:
        from rna_draw.layout.pipeline import layout_guaranteed

        structures = _corpus_sample(_SAMPLE_N_PER_BUCKET)
        if not structures:
            pytest.skip("corpus directory not available in this environment")
        for i, secstruct in enumerate(structures):
            result = layout_guaranteed(secstruct)
            pair_map = result.pair_map or get_pairmap_from_secstruct(secstruct)
            params = OverlapParams(
                node_r=result.node_r,
                backbone_half_width=0.75 * result.node_r,
                pair_half_width=0.75 * result.node_r,
            )
            assert_parity(result.x, result.y, pair_map, params, case_name=f"corpus[{i}]")


# ---------------------------------------------------------------------
# (b) The checker's own existing fixtures
# ---------------------------------------------------------------------


def _named_structure(set_path: Path, name: str) -> str:
    import json

    with set_path.open() as handle:
        structures = json.load(handle)
    return next(entry["structure"] for entry in structures if entry["name"] == name)


def _hairpin_coords(secstruct: str) -> tuple[list[float], list[float], list[int]]:
    renderer = RNARenderer()
    renderer.setup_tree(secstruct, NODE_R=10, PRIMARY_SPACE=20, PAIR_SPACE=23)
    x, y = list(renderer.xarray_), list(renderer.yarray_)
    return x, y, get_pairmap_from_secstruct(secstruct)


def _random_case(n: int, seed: int) -> tuple[list[float], list[float], list[int]]:
    rng = random.Random(seed)
    box = 40.0
    x = [rng.uniform(0, box) for _ in range(n)]
    y = [rng.uniform(0, box) for _ in range(n)]
    pair_map = [-1] * n
    shuffled = list(range(n))
    rng.shuffle(shuffled)
    num_pairs = rng.randint(0, n // 2)
    for k in range(num_pairs):
        i, j = shuffled[2 * k], shuffled[2 * k + 1]
        pair_map[i] = j
        pair_map[j] = i
    return x, y, pair_map


@requires_native
class TestExistingFixturesParity:
    """Section (b): the checker's own trusted test inputs."""

    def test_clean_hairpin(self) -> None:
        x, y, pair_map = _hairpin_coords("((((....))))")
        assert_parity(x, y, pair_map, case_name="clean_hairpin")

    @pytest.mark.parametrize("n", [5, 20, 50, 120])
    @pytest.mark.parametrize("seed", range(3))
    def test_random_layouts(self, n: int, seed: int) -> None:
        x, y, pair_map = _random_case(n, seed)
        assert_parity(x, y, pair_map, case_name=f"random[{n},{seed}]")

    def test_crossing_pair_capsules(self) -> None:
        x, y = [0.0, 10.0, 10.0, 0.0], [0.0, 10.0, 0.0, 10.0]
        pair_map = [1, 0, 3, 2]
        params = OverlapParams(node_r=0.01, backbone_half_width=0.01, pair_half_width=1.0)
        assert_parity(x, y, pair_map, params, case_name="crossing_pair_capsules")

    @pytest.mark.parametrize(
        "name", ["bpRNA_CRW_5083.dbn", "bpRNA_SRP_476.dbn", "bpRNA_RFAM_3401.dbn"]
    )
    def test_safe_fallback_circle_named_structures(self, name: str) -> None:
        if not HARD_SET_JSON.exists():
            pytest.skip("benchmarks/hard_set.json not available")
        struct = _named_structure(HARD_SET_JSON, name)
        x, y = SafeFallbackEngine().layout(struct)
        pair_map = get_pairmap_from_secstruct(struct)
        assert_parity(x, y, pair_map, case_name=f"safe_fallback[{name}]")

    def test_safe_fallback_circle_worst_set_684nt(self) -> None:
        if not WORST_SET_JSON.exists():
            pytest.skip("benchmarks/worst_set.json not available")
        struct = _named_structure(WORST_SET_JSON, "bpRNA_RFAM_35409.dbn")
        x, y = SafeFallbackEngine().layout(struct)
        pair_map = get_pairmap_from_secstruct(struct)
        assert_parity(x, y, pair_map, case_name="safe_fallback_worst_684nt")


# ---------------------------------------------------------------------
# (c) Adversarial near-boundary generator (axis-aligned per R3)
# ---------------------------------------------------------------------

_EPS_VALUES = (1e-4, 1e-7, 1e-10, 1e-13)


def _disk_disk_distances(required: float, tol: float) -> list[float]:
    """Center distances straddling the `required - tol` knife-edge, axis-aligned (R3)."""
    boundary = required - tol
    distances = [required, boundary]
    for eps in _EPS_VALUES:
        distances += [boundary - eps, boundary + eps]
    return distances


@requires_native
class TestDiskDiskBoundary:
    @pytest.mark.parametrize("distance", _disk_disk_distances(required=20.0, tol=1e-6))
    def test_axis_aligned_disk_pair(self, distance: float) -> None:
        # dx = distance, dy = 0 (R3): hypot(dx, 0) == |dx| exactly, so the
        # "exact boundary" case really does sit at the literal FP boundary.
        x, y = [0.0, distance], [0.0, 0.0]
        pair_map = [-1, -1]
        params = OverlapParams(node_r=10.0, backbone_half_width=0.0, pair_half_width=0.0)
        assert_parity(x, y, pair_map, params, case_name=f"disk_disk_boundary[{distance!r}]")


@requires_native
class TestDiskCapsuleBoundary:
    @pytest.mark.parametrize(
        "distance", _disk_disk_distances(required=10.0 + 5.0, tol=1e-6)
    )
    def test_axis_aligned_disk_vs_vertical_backbone(self, distance: float) -> None:
        # A vertical backbone capsule (100,0)-(100,100); a disk offset
        # purely horizontally by `distance` at the capsule's interior
        # midpoint, so point_segment_distance's clamped-t branch reduces
        # to the exact horizontal offset (R3).
        x = [100.0, 100.0, 100.0 + distance, 250.0]
        y = [0.0, 100.0, 50.0, 50.0]
        pair_map = [-1, -1, -1, -1]
        params = OverlapParams(node_r=10.0, backbone_half_width=5.0, pair_half_width=0.0)
        assert_parity(x, y, pair_map, params, case_name=f"disk_capsule_boundary[{distance!r}]")


@requires_native
class TestCapsuleCapsuleBoundary:
    @pytest.mark.parametrize("distance", _disk_disk_distances(required=10.0, tol=1e-6))
    def test_axis_aligned_parallel_capsules(self, distance: float) -> None:
        # Two horizontal, x-overlapping backbone capsules -- (0,1) and
        # (2,3), two separate 2-nucleotide chains -- separated purely in y
        # by `distance`: their closest points are axis-aligned (R3), so
        # segment_segment_distance's closest-points solve reduces to |dy|.
        x = [0.0, 100.0, 0.0, 100.0]
        y = [0.0, 0.0, distance, distance]
        pair_map = [-1, -1, -1, -1]
        params = OverlapParams(node_r=0.01, backbone_half_width=5.0, pair_half_width=0.0)
        assert_parity(x, y, pair_map, params, case_name=f"capsule_capsule_boundary[{distance!r}]")

    def test_near_tangent_capsules(self) -> None:
        # A slight tilt so the two segments are close to tangent but not
        # perfectly parallel -- exercises closest_params' general (non-
        # denom==0) branch near its boundary.
        x = [0.0, 100.0, 0.0, 100.0]
        y = [0.0, 0.0, 5.0, 4.999999]
        pair_map = [-1, -1, -1, -1]
        params = OverlapParams(node_r=0.01, backbone_half_width=2.5, pair_half_width=0.0)
        assert_parity(x, y, pair_map, params, case_name="near_tangent_capsules")

    def test_near_parallel_capsules(self) -> None:
        x = [0.0, 100.0, 1.0, 101.0]
        y = [0.0, 0.0, 20.0, 20.001]
        pair_map = [-1, -1, -1, -1]
        params = OverlapParams(node_r=0.01, backbone_half_width=1.0, pair_half_width=0.0)
        assert_parity(x, y, pair_map, params, case_name="near_parallel_capsules")

    def test_degenerate_zero_length_backbone_segment(self) -> None:
        # Nucleotides 1 and 2 coincide exactly: the (1,2) backbone capsule
        # degenerates to a point (drives point_segment_distance's / the
        # segment-segment `seg_degenerate` guards).
        x = [0.0, 50.0, 50.0, 100.0, 52.0]
        y = [0.0, 0.0, 0.0, 0.0, 0.0]
        pair_map = [-1, -1, -1, -1, -1]
        params = OverlapParams(node_r=5.0, backbone_half_width=3.0, pair_half_width=0.0)
        assert_parity(x, y, pair_map, params, case_name="degenerate_zero_length_backbone")


@requires_native
class TestExcludedPairEdgeCases:
    """Both checkers must exclude identically (no witness on either side)."""

    def test_backbone_adjacent_disks_excluded_both_sides(self) -> None:
        x, y = [0.0, 5.0], [0.0, 0.0]
        pair_map = [-1, -1]
        assert_parity(x, y, pair_map, case_name="excluded_backbone_adjacent")

    def test_base_paired_disks_excluded_both_sides(self) -> None:
        x, y = [0.0, 0.0], [0.0, 0.0]
        pair_map = [1, 0]
        assert_parity(x, y, pair_map, case_name="excluded_base_paired")

    def test_disk_on_its_own_capsule_endpoint_excluded_both_sides(self) -> None:
        x, y = [0.0, 20.0, 40.0], [0.0, 0.0, 0.0]
        pair_map = [-1, -1, -1]
        assert_parity(x, y, pair_map, case_name="excluded_disk_own_capsule")


@requires_native
class TestRandomJitterBand:
    """Perturb a clean layout into the boundary band; pass/fail must still agree."""

    @pytest.mark.parametrize("seed", range(8))
    @pytest.mark.parametrize("magnitude", [1e-6, 1e-9, 1e-13])
    def test_jittered_hairpin(self, seed: int, magnitude: float) -> None:
        rng = random.Random(seed * 7919 + int(magnitude * 1e15))
        x, y, pair_map = _hairpin_coords("((((....))))")
        x = [xi + rng.uniform(-magnitude, magnitude) for xi in x]
        y = [yi + rng.uniform(-magnitude, magnitude) for yi in y]
        assert_parity(x, y, pair_map, case_name=f"jitter[{seed},{magnitude}]")


# ---------------------------------------------------------------------
# Batch API: counts must match looping check_overlaps_count; a malformed
# element must yield the sentinel, never raise (plan step 5).
# ---------------------------------------------------------------------


@requires_native
class TestBatchMatchesLooping:
    def test_batch_counts_match_looping_single_calls(self) -> None:
        from rna_draw.overlap_native import check_overlaps_batch_native, check_overlaps_native

        cases = [
            (([0.0, 30.0, 60.0, 0.0], [0.0, 0.0, 0.0, 0.0]), [-1, -1, -1, -1], 10.0),
            (([0.0, 20.0], [0.0, 0.0]), [-1, -1], 10.0),
            (([0.0, 100.0, 20.0], [0.0, 100.0, 0.0]), [-1, -1, -1], 10.0),
        ]
        expected = [
            check_overlaps_native(xy[0], xy[1], pm, OverlapParams(node_r=node_r))
            for xy, pm, node_r in cases
        ]
        batched = check_overlaps_batch_native(
            [xy for xy, _pm, _r in cases],
            [pm for _xy, pm, _r in cases],
            [node_r for _xy, _pm, node_r in cases],
        )
        assert batched == expected

    def test_batch_isolates_a_malformed_element(self) -> None:
        from rna_draw.overlap_native import check_overlaps_batch_native

        good = (([0.0, 30.0], [0.0, 0.0]), [-1, -1], 10.0)
        malformed = (([0.0, 1.0], [0.0]), [-1, -1], 10.0)  # length mismatch
        cases = [good, malformed, good]
        batched = check_overlaps_batch_native(
            [xy for xy, _pm, _r in cases],
            [pm for _xy, pm, _r in cases],
            [node_r for _xy, _pm, node_r in cases],
        )
        assert batched[1] == -1
        assert batched[0] == batched[2] != -1

    def test_batch_mismatched_outer_lengths_raises_value_error(self) -> None:
        """CPP-REVIEWER WARNING fix: differing outer-list lengths (a
        differing number of structures across coords_list/pair_maps/
        node_rs) must raise `ValueError` up front, not read past the end of
        the shorter list inside a worker thread."""
        from rna_draw.overlap_native import check_overlaps_batch_native

        coords_list = [(([0.0, 30.0], [0.0, 0.0])), (([0.0, 60.0], [0.0, 0.0]))]
        pair_maps = [[-1, -1]]  # one fewer than coords_list/node_rs
        node_rs = [10.0, 10.0]
        with pytest.raises(ValueError):
            check_overlaps_batch_native(coords_list, pair_maps, node_rs)
