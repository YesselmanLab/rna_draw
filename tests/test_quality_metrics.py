"""Tests for `benchmarks.quality_metrics.compute_readability`.

The harness itself is a process-pool driver script (like `pipeline_qc.py`,
which has no unit test); the one piece worth pinning is the pure metric
math, since the readability index becomes the yardstick for judging future
layout engines. These lock its definition against a hand-computed geometry.
"""

from __future__ import annotations

from math import isclose

from benchmarks.quality_metrics import compute_readability
from rna_draw.layout.pipeline import layout_guaranteed
from rna_draw.render_rna import get_pairmap_from_secstruct

NODE_R = 10.0


def test_unit_square_metrics():
    # Four nucleotides on a 20x20 square; nt0 pairs nt3 (a single rung of len 20).
    x = [0.0, 20.0, 20.0, 0.0]
    y = [0.0, 0.0, 20.0, 20.0]
    pair_map = [3, -1, -1, 0]
    m = compute_readability(x, y, pair_map, NODE_R)

    assert m["n"] == 4
    # bbox 20*20=400; per nt 100; / node_r**2 (100) == 1.0
    assert isclose(m["area_per_nt_norm"], 1.0)
    # closest pair is an edge of the square, length 20; / node_r == 2.0
    assert isclose(m["min_nn_norm"], 2.0)
    # the single rung (nt0-nt3) has length 20; / node_r == 2.0
    assert isclose(m["min_rung_norm"], 2.0)
    # min feature 20 / sqrt(area)=20 -> 1.0
    assert isclose(m["readability_index"], 1.0)


def test_no_pairs_has_no_rung():
    x = [0.0, 30.0, 15.0]
    y = [0.0, 0.0, 40.0]
    m = compute_readability(x, y, [-1, -1, -1], NODE_R)
    assert m["min_rung_norm"] is None
    assert m["min_nn_norm"] is not None


def test_collinear_layout_has_no_readability():
    # Zero-area bounding box -> readability index is undefined (None), not a crash.
    x = [0.0, 10.0, 20.0]
    y = [0.0, 0.0, 0.0]
    m = compute_readability(x, y, [-1, -1, -1], NODE_R)
    assert m["area_per_nt_norm"] == 0.0
    assert m["readability_index"] is None
    assert isclose(m["min_nn_norm"], 1.0)  # closest gap 10 / node_r


def test_on_real_pipeline_output():
    # A small hairpin lays out cleanly; its readability index is a real ratio.
    secstruct = "((((....))))"
    result = layout_guaranteed(secstruct)
    pair_map = get_pairmap_from_secstruct(secstruct)
    m = compute_readability(result.x, result.y, pair_map, result.node_r)
    assert m["n"] == len(secstruct)
    assert m["min_nn_norm"] > 0
    assert m["readability_index"] is None or 0.0 < m["readability_index"] <= 2.0
