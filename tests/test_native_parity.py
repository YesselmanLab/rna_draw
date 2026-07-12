"""Differential parity harness: native `rna_layout` turtle vs the vendored
RNAturtle oracle (`rna_draw._vienna_layout.plot_coords_turtle`).

Implements parity criterion 1 from `.claude/plans/current-plan.md`'s
"Parity oracle harness" (turtle base, deterministic, tight): after optimal
rigid alignment (Kabsch, allow reflection), per-nucleotide max abs diff
`<= 1e-6`.

FLOAT-PRECISION NOTE (found running this harness, not anticipated by the
plan): the vendored `vrna_plot_coords_turtle_pt` computes in `double`
internally but its PUBLIC C API returns `float*` -- `_vienna_layout`'s
oracle output is therefore float32-truncated, independent of any
algorithmic divergence. Comparing raw doubles directly hits diffs up to
~1.5e-5 on some hand structures (float32 relative precision at the ~100-300
unit coordinate magnitude turtle produces), which would fail a naive 1e-6
gate for reasons that have nothing to do with correctness. Verified this is
*purely* the oracle's own precision ceiling: rounding the native double
output to float32 (`numpy.float32` round-trip) before comparing yields
EXACT (0.0) agreement on every structure tried -- i.e. the two
implementations compute bit-identical doubles; only the oracle's `float*`
return type differs. This harness therefore rounds native output to
float32 before the parity gate (the correct like-for-like comparison), and
separately reports the raw double-vs-oracle-float diff for transparency
(gated at a documented, precision-ceiling-justified 1e-4).
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

import rna_draw._layout_core as native_layout
import rna_draw._vienna_layout as vienna_layout
from rna_draw.layout.base import has_empty_loop, is_pseudoknot_free

TIGHT_TOL = 1e-6  # after float32-rounding native output; plan criterion 1
RAW_DIFF_SAFETY_FACTOR = 8.0  # headroom over 1 ULP of float32 for the raw (informational) check

# Plan criterion 2 ("Config tree + boxes (deterministic)"): config radius/
# angle rel-tol, box-field abs-tol. Both sides dump full doubles as text
# (`%.17g` on the vendored side, pybind11's native double marshalling on the
# native side) -- unlike criterion 1's turtle coords, there is no float32
# return-type ceiling here, so these tolerances are the actual gate, not a
# precision-floor workaround. Empirically (see the module docstring's sibling
# note in the handoff report) every hard_set/hand structure currently comes
# back EXACTLY equal (0.0 diff); the tolerances are kept at the plan's
# documented values rather than tightened further, to avoid a brittle gate
# against legitimate tiny transcendental (asin/sqrt) ULP differences.
TREE_CFG_REL_TOL = 1e-9
TREE_BOX_ABS_TOL = 1e-6

HAND_STRUCTURES: dict[str, str] = {
    "hairpin": "((((....))))",
    "bulge_near_strand": "(.((....)))",
    "bulge_far_strand": "(((....)).)",
    "internal_loop": "((((..((((....))))..))))",
    "multiloop_two_way": "((((...)))(((...))))",
    "multiloop_three_way": "(((((...)))(((...)))(((...)))))",
    "stacked_helices": "((((((((....))))))))",
    "exterior_dangles_both_sides": "..((((....))))..",
    "multi_branch_exterior": "((....))..((....))..((....))",
    "single_dangle_5prime": ".((((....))))",
    "single_dangle_3prime": "((((....)))).",
    "long_unpaired_run": "((((....))))" + "." * 20,
    "two_helix": "((((....))))((((....))))((((....))))",
    "deep_nesting": "(" * 15 + "." * 4 + ")" * 15,
}

HARD_SET_PATH = Path(__file__).parent.parent / "benchmarks" / "hard_set.json"


def _turtle_usable(structure: str) -> bool:
    """Whether both engines' shared malformed-input guards accept `structure`."""
    return is_pseudoknot_free(structure) and not has_empty_loop(structure)


def _hard_set_structures() -> list[str]:
    import json

    entries = json.loads(HARD_SET_PATH.read_text())
    return [e["structure"] for e in entries if _turtle_usable(e["structure"])]


def _kabsch_align(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rigidly align `source` onto `target` (rotation + reflection +
    translation) minimizing RMSD, via the Kabsch algorithm.

    Args:
        source: `(n, 2)` points to align.
        target: `(n, 2)` reference points.

    Returns:
        `source`, rigidly transformed to best fit `target`.
    """
    source_mean = source.mean(axis=0)
    target_mean = target.mean(axis=0)
    source_centered = source - source_mean
    target_centered = target - target_mean

    covariance = source_centered.T @ target_centered
    u, _, vt = np.linalg.svd(covariance)
    # No determinant-sign correction: the plan explicitly allows reflection
    # here (turtle's clockwise/counter-clockwise convention is otherwise
    # arbitrary), so the raw SVD rotation (possibly improper) is used as-is.
    rotation = vt.T @ u.T
    return (source_centered @ rotation.T) + target_mean


def _max_aligned_diff(native_xy: np.ndarray, oracle_xy: np.ndarray) -> float:
    """Max per-nucleotide abs coordinate diff after Kabsch alignment."""
    aligned = _kabsch_align(native_xy, oracle_xy)
    return float(np.max(np.abs(aligned - oracle_xy)))


def _coords(structure: str, engine: ModuleType) -> np.ndarray:
    x, y = engine.plot_coords_turtle(structure)
    return np.column_stack([x, y])


def _assert_turtle_parity(structure: str) -> float:
    """Compare native vs vendored turtle on `structure`; return the tight
    (float32-rounded) aligned max-abs-diff, after asserting it is within
    `TIGHT_TOL` -- the actual parity gate (plan criterion 1).

    Also computes the raw (un-rounded) double-vs-oracle-float diff and
    asserts it against a coordinate-magnitude-SCALED tolerance (float32's
    relative precision, `~1.2e-7`, applied to the largest coordinate
    magnitude in the structure): a FIXED absolute tolerance is unsound here
    because turtle's coordinate magnitudes grow with structure size (long
    backbones walk far from the origin), and float32 rounding error is
    relative, not absolute -- see module docstring.
    """
    native_xy = _coords(structure, native_layout)
    oracle_xy = _coords(structure, vienna_layout)

    raw_diff = _max_aligned_diff(native_xy, oracle_xy)
    float32_eps = float(np.finfo(np.float32).eps)
    raw_tol = float32_eps * float(np.max(np.abs(oracle_xy))) * RAW_DIFF_SAFETY_FACTOR
    assert raw_diff <= raw_tol, (
        f"{structure!r}: raw double-vs-oracle-float diff {raw_diff} exceeds "
        f"the magnitude-scaled tolerance {raw_tol} (the oracle's float32 "
        "return-type precision ceiling -- see module docstring)"
    )

    native_f32_xy: np.ndarray = native_xy.astype(np.float32).astype(np.float64)
    tight_diff = _max_aligned_diff(native_f32_xy, oracle_xy)
    assert tight_diff <= TIGHT_TOL, (
        f"{structure!r}: float32-rounded aligned diff {tight_diff} exceeds "
        f"the {TIGHT_TOL} parity gate"
    )
    return tight_diff


class TestHandCorpusParity:
    """Fast tier: hand structures covering every motif (plan corpus spec)."""

    @pytest.mark.parametrize("name", sorted(HAND_STRUCTURES))
    def test_turtle_parity(self, name: str) -> None:
        _assert_turtle_parity(HAND_STRUCTURES[name])

    @pytest.mark.parametrize("name", sorted(HAND_STRUCTURES))
    def test_native_turtle_is_deterministic(self, name: str) -> None:
        structure = HAND_STRUCTURES[name]
        first = native_layout.plot_coords_turtle(structure)
        second = native_layout.plot_coords_turtle(structure)
        assert first == second

    @pytest.mark.parametrize("name", sorted(HAND_STRUCTURES))
    def test_vendored_turtle_is_deterministic(self, name: str) -> None:
        structure = HAND_STRUCTURES[name]
        first = vienna_layout.plot_coords_turtle(structure)
        second = vienna_layout.plot_coords_turtle(structure)
        assert first == second


class TestHardSetParity:
    """Fast tier: the frozen 450-structure hard set (turtle has no resolver
    iteration, so this runs in well under a second -- see `benchmarks/
    hard_set.json`; ~14 of 450 are pseudoknots or contain an empty loop and
    are skipped, matching both engines' shared input guards).
    """

    def test_turtle_parity_over_hard_set(self) -> None:
        structures = _hard_set_structures()
        assert len(structures) > 400, "expected most of the hard set to be turtle-usable"

        diffs = [_assert_turtle_parity(s) for s in structures]

        assert max(diffs) <= TIGHT_TOL
        # Report-worthy summary for the handoff message; not itself a gate.
        print(
            f"\nhard_set turtle parity: n={len(structures)}, "
            f"max_diff={max(diffs):.3e}, mean_diff={sum(diffs) / len(diffs):.3e}"
        )


def _tree_dumps(structure: str, paired: float = 35.0, unpaired: float = 25.0) -> tuple[list, list]:
    """`(native, vendored)` T1 tree dumps for `structure`, both as
    `list[dict]` (the vendored side is JSON text over the wire; decode it).
    """
    native_tree = native_layout.dump_tree(structure, paired, unpaired)
    vendored_tree = json.loads(vienna_layout.dump_tree(structure, paired, unpaired))
    return native_tree, vendored_tree


def _assert_tree_topology_exact(native_tree: list, vendored_tree: list, structure: str) -> None:
    """Plan criterion 2's topology half: identical node count, and for each
    node (matched by DFS pre-order position -- both dumps use the same
    order, see `debug_dump.hpp`) identical id/parent_id/loop_start/
    stem_start.
    """
    assert len(native_tree) == len(vendored_tree), (
        f"{structure!r}: node count differs (native={len(native_tree)}, "
        f"vendored={len(vendored_tree)})"
    )
    for native_node, vendored_node in zip(native_tree, vendored_tree):
        for key in ("id", "parent_id", "loop_start", "stem_start"):
            assert native_node[key] == vendored_node[key], (
                f"{structure!r} node {native_node['id']}: {key} differs "
                f"(native={native_node[key]}, vendored={vendored_node[key]})"
            )
        assert (native_node["cfg"] is None) == (vendored_node["cfg"] is None), (
            f"{structure!r} node {native_node['id']}: cfg presence differs"
        )
        assert (native_node["lbox"] is None) == (vendored_node["lbox"] is None), (
            f"{structure!r} node {native_node['id']}: lbox presence differs"
        )
        assert (native_node["sbox"] is None) == (vendored_node["sbox"] is None), (
            f"{structure!r} node {native_node['id']}: sbox presence differs"
        )


def _assert_rel_close(actual: float, expected: float, tol: float, what: str) -> None:
    rel = abs(actual - expected) / max(abs(expected), 1e-12)
    assert rel <= tol, (
        f"{what}: rel diff {rel} exceeds {tol} (native={actual}, vendored={expected})"
    )


def _assert_abs_close(actual: float, expected: float, tol: float, what: str) -> None:
    diff = abs(actual - expected)
    assert diff <= tol, (
        f"{what}: abs diff {diff} exceeds {tol} (native={actual}, vendored={expected})"
    )


def _assert_tree_numeric_parity(native_tree: list, vendored_tree: list, structure: str) -> None:
    """Plan criterion 2's numeric half: config radius/angle rel-tol
    `<= TREE_CFG_REL_TOL`, box fields abs-tol `<= TREE_BOX_ABS_TOL`.
    """
    for native_node, vendored_node in zip(native_tree, vendored_tree):
        node_id = native_node["id"]
        if native_node["cfg"] is not None:
            for key in ("radius", "min_radius", "default_radius"):
                _assert_rel_close(
                    native_node["cfg"][key],
                    vendored_node["cfg"][key],
                    TREE_CFG_REL_TOL,
                    f"{structure!r} node {node_id} cfg.{key}",
                )
            native_arcs = native_node["cfg"]["arcs"]
            vendored_arcs = vendored_node["cfg"]["arcs"]
            assert len(native_arcs) == len(vendored_arcs), (
                f"{structure!r} node {node_id}: arc count differs "
                f"(native={len(native_arcs)}, vendored={len(vendored_arcs)})"
            )
            for native_arc, vendored_arc in zip(native_arcs, vendored_arcs):
                assert native_arc["segments"] == vendored_arc["segments"], (
                    f"{structure!r} node {node_id}: arc segment count differs"
                )
                _assert_rel_close(
                    native_arc["angle"],
                    vendored_arc["angle"],
                    TREE_CFG_REL_TOL,
                    f"{structure!r} node {node_id} arc.angle",
                )
        if native_node["lbox"] is not None:
            for key in ("cx", "cy", "r"):
                _assert_abs_close(
                    native_node["lbox"][key],
                    vendored_node["lbox"][key],
                    TREE_BOX_ABS_TOL,
                    f"{structure!r} node {node_id} lbox.{key}",
                )
        if native_node["sbox"] is not None:
            for key in ("ax", "ay", "bx", "by", "cx", "cy", "ex", "ey", "bulge_dist"):
                _assert_abs_close(
                    native_node["sbox"][key],
                    vendored_node["sbox"][key],
                    TREE_BOX_ABS_TOL,
                    f"{structure!r} node {node_id} sbox.{key}",
                )
            assert native_node["sbox"]["bulge_count"] == vendored_node["sbox"]["bulge_count"], (
                f"{structure!r} node {node_id}: bulge_count differs"
            )


def _assert_tree_parity(structure: str) -> None:
    native_tree, vendored_tree = _tree_dumps(structure)
    _assert_tree_topology_exact(native_tree, vendored_tree, structure)
    _assert_tree_numeric_parity(native_tree, vendored_tree, structure)


class TestTreeParity:
    """Plan criterion 2 ("Config tree + boxes"): `dump_tree`@T1 (post
    `update_bounding_boxes`, pre-resolver) topology EXACT + config/box
    numeric parity, over the hand/motif corpus + the hard set.
    """

    @pytest.mark.parametrize("name", sorted(HAND_STRUCTURES))
    def test_tree_parity_on_hand_corpus(self, name: str) -> None:
        _assert_tree_parity(HAND_STRUCTURES[name])

    @pytest.mark.parametrize("name", sorted(HAND_STRUCTURES))
    def test_native_tree_is_deterministic(self, name: str) -> None:
        structure = HAND_STRUCTURES[name]
        first = native_layout.dump_tree(structure)
        second = native_layout.dump_tree(structure)
        assert first == second

    def test_tree_parity_over_hard_set(self) -> None:
        structures = _hard_set_structures()
        assert len(structures) > 400, "expected most of the hard set to be turtle-usable"

        for structure in structures:
            _assert_tree_parity(structure)

        print(f"\nhard_set tree parity: n={len(structures)}, topology exact + numeric parity OK")


@pytest.mark.skipif(
    not os.environ.get("RNA_DRAW_PARITY_BROAD"),
    reason="broad ~2k dbnFiles sample is on-demand (set RNA_DRAW_PARITY_BROAD=1); "
    "the hand corpus + hard set already cover this in every default run",
)
class TestBroadSampleParity:
    """On-demand tier: ~2000 pseudoknot-free structures from the user's
    bpRNA corpus (`~/Downloads/dbnFiles/`), per the plan's corpus spec.
    """

    def test_turtle_parity_over_broad_sample(self) -> None:
        from benchmarks.hard_gate import CORPUS, parse_dbn

        files = sorted(CORPUS.glob("*.dbn"), key=lambda p: p.name)
        sample_stride = max(1, len(files) // 2000)
        structures = []
        for path in files[::sample_stride]:
            parsed = parse_dbn(path)
            if parsed is None:
                continue
            _, structure = parsed
            if _turtle_usable(structure):
                structures.append(structure)

        assert len(structures) > 500, "expected a substantial broad sample"
        diffs = [_assert_turtle_parity(s) for s in structures]
        print(
            f"\nbroad-sample turtle parity: n={len(structures)}, "
            f"max_diff={max(diffs):.3e}, mean_diff={sum(diffs) / len(diffs):.3e}"
        )
