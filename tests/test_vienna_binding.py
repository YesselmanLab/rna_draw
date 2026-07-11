"""Tests for the in-process ViennaRNA bindings (`rna_draw._vienna_layout`)
and the three `LayoutEngine`s built on top of them
(`rna_draw.layout.vienna`).

Covers: coordinate parity against the subprocess `PuzzlerEngine` oracle
(including >=2 structures puzzler lays out DIRTY, so the parity gate is
non-vacuous), malformed/empty/pseudoknot input handling (must raise, never
segfault), the ABI version/sizeof drift guard, memory-safety stress loops
(proving the RAII `free()` path runs), and that all three engines slot
into `layout_guaranteed` honestly (clean or flagged, never a silent
overlap).
"""

from __future__ import annotations

import resource
import shutil
import subprocess
import sys

import pytest
from conftest import random_structure

import rna_draw._vienna_layout as vienna_layout
from rna_draw.layout.base import EngineError, EngineUnavailableError
from rna_draw.layout.pipeline import (
    EXPECTED_PUZZLER_OPTIONS_SIZEOF,
    EXPECTED_VIENNA_ABI_VERSION,
    layout_guaranteed,
    resolve_engine,
)
from rna_draw.layout.puzzler import PuzzlerEngine
from rna_draw.layout.vienna import ViennaNaviewEngine, ViennaPuzzlerEngine, ViennaTurtleEngine
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct

ATOL_COORD = 1e-4  # in-process vs subprocess call the SAME algorithm -- near-exact
ATOL_BBOX = 1e-4

# Curated bpRNA_CRW/SRP structures (300-800 nt) that ViennaRNA's puzzler
# itself lays out DIRTY at the default checker geometry (non-zero M2
# witnesses at node_r=10) -- found by sampling the M5.1 dbnFiles corpus.
# Non-vacuous parity requires this: a parity check on only clean hairpins
# would pass 0==0 regardless of whether the in-process binding is correct.
DIRTY_STRUCTURES = {
    "bpRNA_CRW_15054": (
        "....................(((((((.....((((((((..(((((((....))))))).((((....)))))))"
        ")).))).((.((((..(((((((((...(((((((((....)))..((((......))))..)))))).....((("
        "(.(((((((...((..((......))))....)))))))..((((((((.....)))))))).....))))....)"
        "))).)))...))))))))....)))))))...................(((((((.....(((..((...(((..."
        ".)))...))....))).....)))))))......(...((((((((........))))))))...).........."
        ".....((((((((.......))))))))................................................"
        "...................."
    ),
    "bpRNA_CRW_9308": (
        "((((((....)))..((((......))))..)))........((((.(((((((...((..((......))))..."
        ".)))))))..((((((((.....)))))))).....))))...................................."
        ".................(((((((.....(((..((...(((....)))...))....))).....)))))))..."
        "...(...((((((((........))))))))...)...............((((((((.......))))))))..."
        "..................................(((((((..((((((((((((....((((((.((((..((.."
        "...)).))))))))))...))))))))))))..)))))))"
    ),
    "bpRNA_SRP_623": (
        "(.(((((((((((((((..))))))))))))).....(.((((((....(((((....)))))....)))))).))"
        ").)(((((.(((((.((((((..(((.((((((.((.((((((((((.(((((..((((((((((((.(((((((("
        "((.((....))))))))))))..)))).(...).(((((.....((((....(((....)))....))))..))))"
        ").))))))))....)))))))...)))).))))..)).))))))...)))))))))...)))))))).))......"
        "..............."
    ),
}

CLEAN_STRUCTURES = {
    "hairpin": "((((....))))",
    "two_helix": "((((....))))((((....))))((((....))))",
    "curated_puzzler_clean": (
        ".(((((((....(((..........))))))))))(((((((.((((.....))))..(((((......)))))"
        "(((....)))......)))))))..\n"
    ).strip(),
}

HAVE_VIENNARNA = shutil.which("RNAfold") is not None and shutil.which("RNAplot") is not None
requires_viennarna = pytest.mark.skipif(
    not HAVE_VIENNARNA, reason="RNAfold/RNAplot (ViennaRNA) not found on PATH"
)


def _rnaplot_version() -> str:
    """The `RNAplot --version` string, or `""` if the binary is absent."""
    if shutil.which("RNAplot") is None:
        return ""
    result = subprocess.run(
        ["RNAplot", "--version"], capture_output=True, text=True, check=True
    )
    return result.stdout.strip()


def _bbox_aspect_ratio(x: list[float], y: list[float]) -> float:
    """`min(width, height) / max(width, height)` -- a rigid-transform (up
    to axis swap) shape signature independent of absolute position/scale.
    """
    width = max(x) - min(x)
    height = max(y) - min(y)
    lo, hi = sorted((width, height))
    return lo / hi if hi else 0.0


class TestModuleAbiFacts:
    def test_version_matches_expected(self) -> None:
        assert vienna_layout.version() == "2.7.0"

    def test_abi_version_matches_expected(self) -> None:
        assert vienna_layout.abi_version() == EXPECTED_VIENNA_ABI_VERSION

    def test_sizeof_puzzler_options_matches_recorded_constant(self) -> None:
        assert vienna_layout.sizeof_puzzler_options() == EXPECTED_PUZZLER_OPTIONS_SIZEOF


class TestMalformedInputBindingLevel:
    """The raw `_vienna_layout` functions must raise, never segfault."""

    BINDING_FNS = [
        vienna_layout.plot_coords_puzzler,
        vienna_layout.plot_coords_naview,
        vienna_layout.plot_coords_turtle,
    ]

    @pytest.mark.parametrize("secstruct", ["(", ")", "((", "))", ")(", "[.]"])
    @pytest.mark.parametrize("fn_index", range(3))
    def test_malformed_input_raises_value_error(self, fn_index: int, secstruct: str) -> None:
        with pytest.raises(ValueError):
            self.BINDING_FNS[fn_index](secstruct)

    @pytest.mark.parametrize("fn_index", range(3))
    def test_empty_structure_raises_runtime_error(self, fn_index: int) -> None:
        with pytest.raises(RuntimeError):
            self.BINDING_FNS[fn_index]("")

    @pytest.mark.parametrize("fn_index", range(3))
    def test_single_nucleotide_succeeds(self, fn_index: int) -> None:
        x, y = self.BINDING_FNS[fn_index](".")
        assert len(x) == 1
        assert len(y) == 1


class TestEmptyLoopGuard:
    """A bare empty loop (`"()"` with nothing between the pair) is a
    verified hang/crash trigger, not merely "malformed": M5.1 spike found
    `vrna_plot_coords_puzzler("().()")` LOOPS (does not return; its
    iterative intersection-resolution never converges on a
    zero-nucleotide loop) and `vrna_plot_coords_turtle` SEGFAULTS on the
    same input -- naview alone tolerates it. `@pytest.mark.timeout` caps
    each case so a regression that reintroduces the hang fails fast
    instead of hanging the whole suite.
    """

    @pytest.mark.timeout(5)
    def test_puzzler_raises_value_error_not_hangs(self) -> None:
        with pytest.raises(ValueError, match="empty loop"):
            vienna_layout.plot_coords_puzzler("().()")

    @pytest.mark.timeout(5)
    def test_turtle_raises_value_error_not_segfaults(self) -> None:
        with pytest.raises(ValueError, match="empty loop"):
            vienna_layout.plot_coords_turtle("().()")

    @pytest.mark.timeout(5)
    def test_naview_tolerates_empty_loop(self) -> None:
        x, y = vienna_layout.plot_coords_naview("().()")
        assert len(x) == 5
        assert len(y) == 5

    @pytest.mark.timeout(5)
    @pytest.mark.parametrize(
        "engine_factory", [ViennaPuzzlerEngine, ViennaTurtleEngine]
    )
    def test_engine_layout_raises_engine_error(self, engine_factory: type) -> None:
        # `is_pseudoknot_free` (well-nested check) does not reject "().()"
        # -- it IS well-nested -- so this exercises the binding's
        # ValueError being caught and re-raised as EngineError, not the
        # earlier is_pseudoknot_free short-circuit (EngineUnavailableError).
        with pytest.raises(EngineError):
            engine_factory().layout("().()")

    @pytest.mark.timeout(5)
    @pytest.mark.parametrize(
        "engine_factory", [ViennaPuzzlerEngine, ViennaTurtleEngine]
    )
    def test_pipeline_falls_back_instead_of_hanging(self, engine_factory: type) -> None:
        result = layout_guaranteed("().()", engine=engine_factory())
        assert result.flagged is True
        assert result.engine_name == "fallback"

    @pytest.mark.timeout(5)
    def test_naview_pipeline_succeeds_directly(self) -> None:
        result = layout_guaranteed("().()", engine=ViennaNaviewEngine())
        assert result.flagged is False
        assert result.engine_name == "naview"


class TestPuzzlerOptsResolverLevers:
    """`plot_coords_puzzler_opts` exposes two resolver levers on top of the
    vendored (editable) RNApuzzler: `allow_flipping` and
    `max_config_changes`. Its defaults must reproduce `plot_coords_puzzler`
    exactly -- see the vendored `RNApuzzler.c` edit that makes the
    25000-change budget caller-respectable (only a fallback for <= 0)
    instead of hardcoded, which this binding's default (0) still triggers.
    """

    STRUCTURES = [
        CLEAN_STRUCTURES["hairpin"],
        CLEAN_STRUCTURES["two_helix"],
        CLEAN_STRUCTURES["curated_puzzler_clean"],
        DIRTY_STRUCTURES["bpRNA_CRW_9308"],
    ]

    @pytest.mark.parametrize("secstruct", STRUCTURES)
    def test_defaults_match_plot_coords_puzzler(self, secstruct: str) -> None:
        want_x, want_y = vienna_layout.plot_coords_puzzler(secstruct)
        got_x, got_y = vienna_layout.plot_coords_puzzler_opts(secstruct)
        assert got_x == pytest.approx(want_x, abs=ATOL_COORD)
        assert got_y == pytest.approx(want_y, abs=ATOL_COORD)

    def test_large_budget_does_not_change_small_clean_structure(self) -> None:
        secstruct = CLEAN_STRUCTURES["hairpin"]
        default_x, default_y = vienna_layout.plot_coords_puzzler_opts(secstruct)
        big_budget_x, big_budget_y = vienna_layout.plot_coords_puzzler_opts(
            secstruct, max_config_changes=1_000_000
        )
        assert big_budget_x == pytest.approx(default_x, abs=ATOL_COORD)
        assert big_budget_y == pytest.approx(default_y, abs=ATOL_COORD)

    def test_allow_flipping_returns_valid_finite_coords(self) -> None:
        secstruct = CLEAN_STRUCTURES["curated_puzzler_clean"]
        x, y = vienna_layout.plot_coords_puzzler_opts(secstruct, allow_flipping=True)
        assert len(x) == len(secstruct)
        assert len(y) == len(secstruct)
        assert all(v == v and abs(v) < 1e8 for v in x)  # v == v excludes NaN
        assert all(v == v and abs(v) < 1e8 for v in y)

    @pytest.mark.parametrize("secstruct", ["(", ")", "((", "))", ")(", "[.]"])
    def test_malformed_input_raises_value_error(self, secstruct: str) -> None:
        with pytest.raises(ValueError):
            vienna_layout.plot_coords_puzzler_opts(secstruct)

    def test_empty_structure_raises_runtime_error(self) -> None:
        with pytest.raises(RuntimeError):
            vienna_layout.plot_coords_puzzler_opts("")

    def test_empty_loop_raises_value_error_not_hangs(self) -> None:
        with pytest.raises(ValueError, match="empty loop"):
            vienna_layout.plot_coords_puzzler_opts("().()")


class TestMalformedInputEngineLevel:
    """The `LayoutEngine` wrappers guard with `is_pseudoknot_free` first."""

    ENGINES = [ViennaPuzzlerEngine(), ViennaNaviewEngine(), ViennaTurtleEngine()]

    @pytest.mark.parametrize("secstruct", ["([)]", "((", "(]"])
    @pytest.mark.parametrize("engine_index", range(3))
    def test_pseudoknot_or_unbalanced_raises_unavailable(
        self, engine_index: int, secstruct: str
    ) -> None:
        with pytest.raises(EngineUnavailableError):
            self.ENGINES[engine_index].layout(secstruct)

    @pytest.mark.parametrize("engine_index", range(3))
    def test_empty_structure_returns_empty_lists(self, engine_index: int) -> None:
        assert self.ENGINES[engine_index].layout("") == ([], [])

    @pytest.mark.parametrize("engine_index", range(3))
    def test_single_nucleotide_returns_origin(self, engine_index: int) -> None:
        assert self.ENGINES[engine_index].layout(".") == ([0.0], [0.0])


@requires_viennarna
class TestParityVsSubprocessOracle:
    """Non-vacuous parity: >=2 DIRTY structures + clean ones, coordinate
    match + witness-count match (non-zero on the dirty cases) + a
    rigid-transform-invariant bbox-aspect-ratio shape signature.
    """

    def test_oracle_is_same_viennarna_version_as_the_linked_archive(self) -> None:
        # If RNAplot on PATH were a different ViennaRNA build, the two
        # puzzler paths could legitimately diverge and parity would be
        # meaningless -- fail loudly here rather than silently comparing
        # apples to oranges.
        assert _rnaplot_version() == "RNAplot 2.7.0"

    @pytest.mark.parametrize("name", list(DIRTY_STRUCTURES))
    def test_dirty_structures_match_with_nonzero_witness_count(self, name: str) -> None:
        secstruct = DIRTY_STRUCTURES[name]
        self._assert_parity(secstruct, expect_nonzero_witnesses=True)

    @pytest.mark.parametrize("name", list(CLEAN_STRUCTURES))
    def test_clean_structures_match(self, name: str) -> None:
        secstruct = CLEAN_STRUCTURES[name]
        self._assert_parity(secstruct, expect_nonzero_witnesses=False)

    def _assert_parity(self, secstruct: str, expect_nonzero_witnesses: bool) -> None:
        oracle_x, oracle_y = PuzzlerEngine().layout(secstruct)
        inproc_x, inproc_y = ViennaPuzzlerEngine().layout(secstruct)

        for got, want in zip(inproc_x + inproc_y, oracle_x + oracle_y):
            assert got == pytest.approx(want, abs=ATOL_COORD)

        pair_map = get_pairmap_from_secstruct(secstruct)
        params = OverlapParams(node_r=DrawParameters().NODE_R)
        oracle_report = check_overlaps(oracle_x, oracle_y, pair_map, params)
        inproc_report = check_overlaps(inproc_x, inproc_y, pair_map, params)
        assert inproc_report.num_overlaps == oracle_report.num_overlaps
        if expect_nonzero_witnesses:
            assert oracle_report.num_overlaps > 0

        oracle_bbox = _bbox_aspect_ratio(oracle_x, oracle_y)
        inproc_bbox = _bbox_aspect_ratio(inproc_x, inproc_y)
        assert inproc_bbox == pytest.approx(oracle_bbox, abs=ATOL_BBOX)


def _rss_bytes() -> int:
    """Peak RSS so far, normalized to bytes (Linux reports KB, macOS bytes)."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak if sys.platform == "darwin" else peak * 1024


def _no_empty_loop_structure(n: int) -> str:
    """A deterministic length-`n` structure with no bare `"()"` loop.

    `random_structure` reliably contains a `"()"` empty loop for n>=~50
    (statistically near-certain -- every seed tried for n=500 had one), and
    that shape hangs puzzler / segfaults turtle (see `TestEmptyLoopGuard`).
    The memory-safety stress loops below need a large structure that is
    NOT that degenerate case, so build one from repeated 4-nt-loop
    hairpins instead.
    """
    unit = "((((....))))"
    repeats, remainder = divmod(n, len(unit))
    return unit * repeats + "." * remainder


class TestMemorySafety:
    """Stress loops proving the RAII `free()` path runs (no leak growth
    beyond a small bound) for the reentrant engines, and that turtle's
    verified upstream leak (see below) is at least BOUNDED, not unbounded.
    """

    STRUCTURE = _no_empty_loop_structure(500)

    def test_puzzler_no_leak_over_2000_calls(self) -> None:
        self._assert_bounded_growth(
            vienna_layout.plot_coords_puzzler, iterations=2000, bound_bytes=5_000_000
        )

    def test_naview_no_leak_over_2000_calls(self) -> None:
        self._assert_bounded_growth(
            vienna_layout.plot_coords_naview, iterations=2000, bound_bytes=5_000_000
        )

    def test_turtle_leak_is_bounded_known_upstream_issue(self) -> None:
        # KNOWN UPSTREAM BUG (verified against libRNA.a with a standalone C
        # harness that frees x/y/arc_coords correctly, outside Python/pybind
        # entirely): `vrna_plot_coords_turtle` leaks internally,
        # proportional to structure length, independent of this binding's
        # free() calls (which do run -- puzzler/naview with the identical
        # call pattern do not leak). Not fixable from this binding; tracked
        # here with a generous bound so a REGRESSION (leak rate getting much
        # worse) is still caught, without blocking M5.1 on an upstream fix.
        self._assert_bounded_growth(
            vienna_layout.plot_coords_turtle, iterations=200, bound_bytes=10_000_000
        )

    def _assert_bounded_growth(self, fn, iterations: int, bound_bytes: int) -> None:
        fn(self.STRUCTURE)  # warm up any one-time allocation before measuring
        before = _rss_bytes()
        for _ in range(iterations):
            fn(self.STRUCTURE)
        after = _rss_bytes()
        assert after - before < bound_bytes, (
            f"RSS grew {after - before} bytes over {iterations} calls "
            f"(bound {bound_bytes}) -- possible leak"
        )


class TestAbiDriftGuard:
    """The startup ABI assert must actually FIRE on a recorded mismatch,
    not just pass on the current (matching) build.
    """

    def test_wrong_version_raises_unavailable(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr(vienna_layout, "abi_version", lambda: (9, 9, 9))
        with pytest.raises(EngineUnavailableError):
            resolve_engine("vienna_puzzler")
        # The guard is in the engine's __init__, so direct construction is
        # covered too (M5.2 constructs engines directly for option sweeps).
        with pytest.raises(EngineUnavailableError):
            ViennaPuzzlerEngine()

    def test_wrong_sizeof_raises_unavailable(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr(vienna_layout, "sizeof_puzzler_options", lambda: -1)
        with pytest.raises(EngineUnavailableError):
            resolve_engine("naview")

    def test_matching_abi_does_not_raise(self) -> None:
        engine = resolve_engine("turtle")
        assert engine is not None
        assert engine.name == "turtle"


class TestNaviewSingleThreaded:
    """naview is NOT reentrant (see `bindings.cpp`/`vienna.py` docstrings);
    this test only exercises it single-threaded, and checks the coordinates
    it returns are sane (finite, non-degenerate).
    """

    def test_produces_sane_finite_coords(self) -> None:
        secstruct = random_structure(seed=7, n=200)
        x, y = ViennaNaviewEngine().layout(secstruct)
        assert len(x) == len(secstruct)
        assert all(v == v and abs(v) < 1e8 for v in x)  # v == v excludes NaN
        assert all(v == v and abs(v) < 1e8 for v in y)


class TestResolveEngineNewNames:
    @pytest.mark.parametrize(
        "name, expected_class",
        [
            ("vienna_puzzler", ViennaPuzzlerEngine),
            ("naview", ViennaNaviewEngine),
            ("turtle", ViennaTurtleEngine),
        ],
    )
    def test_resolves_to_expected_engine(self, name: str, expected_class: type) -> None:
        assert isinstance(resolve_engine(name), expected_class)


class TestEnginesNeverSilentlyOverlapInPipeline:
    """Every in-process engine must slot into `layout_guaranteed` honestly:
    checker-clean, or explicitly flagged -- never a silent overlap.

    `random_structure` can independently generate other bare-empty-loop
    (or otherwise degenerate) shapes besides the curated `"().()"` case
    above; `@pytest.mark.timeout` is a safety net so any such case fails
    fast rather than hanging this sweep.
    """

    @pytest.mark.timeout(10)
    @pytest.mark.parametrize(
        "engine_factory", [ViennaPuzzlerEngine, ViennaNaviewEngine, ViennaTurtleEngine]
    )
    @pytest.mark.parametrize("seed", range(5))
    @pytest.mark.parametrize("n", [5, 20, 60])
    def test_never_silent_overlap(self, engine_factory: type, seed: int, n: int) -> None:
        secstruct = random_structure(seed, n)
        result = layout_guaranteed(secstruct, engine=engine_factory())
        assert result.flagged or result.report.passed

    def test_curated_puzzler_clean_case_is_clean_and_unflagged_for_vienna_puzzler(self) -> None:
        # This structure was hand-picked (see `test_puzzler_engine.py`) as
        # one the SUBPROCESS puzzler lays out checker-clean; only assert
        # the equivalent in-process engine matches that -- naview/turtle
        # use a different algorithm and are not expected to agree (the
        # general never-silent-overlap property above already covers them).
        secstruct = CLEAN_STRUCTURES["curated_puzzler_clean"]
        result = layout_guaranteed(secstruct, engine=ViennaPuzzlerEngine())
        assert result.flagged is False
        assert result.report.passed
        assert result.engine_name == "vienna_puzzler"
