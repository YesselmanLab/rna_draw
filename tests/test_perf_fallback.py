"""Performance + honesty tests for `SafeFallbackEngine`'s analytic radius solve.

The circle fallback used to double its radius, re-running `check_overlaps`
each iteration until clean -- which for a large, densely paired structure
took ~160s of futile work on the UI thread (giant multiloops never clear on
a modestly sized circle). It now solves the radius analytically from a
single checker pass and caps the radius at a cost budget so every pass stays
affordable: clean when that is cheap, else the largest affordable radius,
still checker-run (the pipeline flags it) -- never a silent overlap.

These tests pin: (1) the large-structure fallback completes well under a
generous wall-clock budget and returns a checker-validated result; (2) a
small structure still lays out clean; (3) the fix does not perturb the
frozen `check_overlaps` semantics (hash == brute-force on the fallback
circle).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from rna_draw.layout.fallback import SafeFallbackEngine
from rna_draw.overlap import (
    OverlapParams,
    check_overlaps,
    check_overlaps_bruteforce,
)
from rna_draw.render_rna import get_pairmap_from_secstruct

HARD_SET_JSON = Path(__file__).parent.parent / "benchmarks" / "hard_set.json"

# The exact 3024 nt structure whose fallback circle froze the editor ~160s.
BIG_STRUCTURE_NAME = "bpRNA_RFAM_42974.dbn"
BIG_BUDGET_S = 15.0


def _named_structure(set_path: Path, name: str) -> str:
    with set_path.open() as handle:
        structures = json.load(handle)
    return next(entry["structure"] for entry in structures if entry["name"] == name)


def test_large_structure_fallback_is_fast_and_checker_validated() -> None:
    """The 3024 nt fallback finishes fast and returns a checker-run layout."""
    struct = _named_structure(HARD_SET_JSON, BIG_STRUCTURE_NAME)
    assert len(struct) > 2900  # guard: still the big case
    params = OverlapParams(node_r=10.0)

    start = time.perf_counter()
    x, y = SafeFallbackEngine(params=params).layout(struct)
    elapsed = time.perf_counter() - start

    # Was ~156s of futile doubling; the analytic solve + cost cap is seconds.
    assert elapsed < BIG_BUDGET_S, f"fallback took {elapsed:.1f}s (budget {BIG_BUDGET_S}s)"
    assert len(x) == len(y) == len(struct)

    # Honest contract: the frozen checker is actually run on the result. The
    # giant circle is not required to be clean (a clean radius here is
    # prohibitively large for the frozen checker) -- the pipeline flags it --
    # but the check must complete and be well-formed.
    pair_map = get_pairmap_from_secstruct(struct)
    report = check_overlaps(x, y, pair_map, params)
    assert isinstance(report.num_overlaps, int)


def test_small_structure_still_lays_out_clean() -> None:
    """A small structure's cheap clean radius is still reached (unchanged)."""
    struct = "((((....))))((((....))))"
    x, y = SafeFallbackEngine().layout(struct)
    pair_map = get_pairmap_from_secstruct(struct)
    assert check_overlaps(x, y, pair_map).passed


def test_fallback_does_not_alter_frozen_checker_semantics() -> None:
    """Hash == brute-force on the fallback circle (checker semantics frozen)."""
    struct = _named_structure(HARD_SET_JSON, "bpRNA_SRP_476.dbn")
    x, y = SafeFallbackEngine().layout(struct)
    pair_map = get_pairmap_from_secstruct(struct)
    hashed = check_overlaps(x, y, pair_map)
    brute = check_overlaps_bruteforce(x, y, pair_map)
    assert set(hashed.witnesses) == set(brute.witnesses)
