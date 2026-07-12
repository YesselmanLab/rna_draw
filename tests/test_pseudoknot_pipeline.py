"""Tests for the M3 pipeline/renderer wiring: `pipeline._try_pseudoknot` /
`layout_guaranteed`'s new tier, and `draw.py`'s pseudoknot-aware pair/bounds
helpers (MUST-FIX #1/#2).
"""

from __future__ import annotations

from rna_draw import settings
from rna_draw.colorer import COLORS
from rna_draw.draw import RNADrawer, _drawn_pairs, _figure_bounds, _shift_routed_lines
from rna_draw.layout.base import LayoutResult, RoutedLine, empty_report
from rna_draw.layout.pipeline import _try_pseudoknot, layout_guaranteed
from rna_draw.overlap import OverlapParams
from rna_draw.render_rna import RNARenderer, get_pairmap_from_secstruct

H_TYPE = "((((....[[[[....))))....]]]]"


class TestTryPseudoknot:
    def test_returns_pseudoknot_engine_result(self) -> None:
        result = _try_pseudoknot(H_TYPE, OverlapParams())
        assert result is not None
        assert result.engine_name == "pseudoknot"
        assert result.pair_map is not None

    def test_never_a_silent_overlap(self) -> None:
        result = _try_pseudoknot(H_TYPE, OverlapParams())
        assert result is not None
        assert result.flagged or result.report.passed


class TestLayoutGuaranteedPseudoknotTier:
    def test_pseudoknot_input_reaches_pseudoknot_tier(self) -> None:
        result = layout_guaranteed(H_TYPE)
        assert result.engine_name == "pseudoknot"
        assert result.flagged or result.report.passed

    def test_pk_free_input_never_reaches_pseudoknot_tier(self) -> None:
        result = layout_guaranteed("((((....))))")
        assert result.engine_name != "pseudoknot"

    def test_all_crossing_pairs_present_in_full_pair_map(self) -> None:
        result = layout_guaranteed(H_TYPE)
        assert result.pair_map is not None
        for i, j in result.crossing_pairs:
            assert result.pair_map[i] == j
            assert result.pair_map[j] == i


class TestDrawnPairsHelper:
    """`draw._drawn_pairs`: MUST-FIX #1 -- prefer the pk path's FULL
    `pair_map` over the `()`-only view, and color crossing pairs distinctly.
    """

    def test_non_pk_result_uses_get_pairmap_from_secstruct(self) -> None:
        ss = "((((....))))"
        result = LayoutResult(
            x=[0.0] * len(ss),
            y=[0.0] * len(ss),
            engine_name="legacy",
            report=empty_report(),
            flagged=False,
            node_r=10.0,
        )
        pairs = _drawn_pairs(result, ss)
        expected_pairmap = get_pairmap_from_secstruct(ss)
        expected_pairs = {(i, partner) for i, partner in enumerate(expected_pairmap) if partner > i}
        assert {(p["from"], p["to"]) for p in pairs} == expected_pairs
        assert all(p["color"] == COLORS["e"] for p in pairs)

    def test_pk_result_uses_full_pair_map_not_dot_bracket_view(self) -> None:
        # The dot-bracket string's OWN ()-only view differs from the pk
        # path's `pair_map` on purpose here -- MUST-FIX #1's regression
        # target: a naive `get_pairmap_from_secstruct(ss)` read would miss
        # this pair entirely (index 0/3 are `.` in `ss`).
        ss = "...."
        pair_map = [3, -1, -1, 0]
        result = LayoutResult(
            x=[0.0, 1.0, 2.0, 3.0],
            y=[0.0, 0.0, 0.0, 0.0],
            engine_name="pseudoknot",
            report=empty_report(),
            flagged=True,
            node_r=10.0,
            pair_map=pair_map,
            crossing_pairs=[(0, 3)],
        )
        pairs = _drawn_pairs(result, ss)
        assert pairs == [{"from": 0, "to": 3, "p": 1.0, "color": COLORS["o"]}]

    def test_non_crossing_pk_pair_uses_default_color(self) -> None:
        ss = "...."
        result = LayoutResult(
            x=[0.0, 1.0, 2.0, 3.0],
            y=[0.0, 0.0, 0.0, 0.0],
            engine_name="pseudoknot",
            report=empty_report(),
            flagged=False,
            node_r=10.0,
            pair_map=[3, -1, -1, 0],
            crossing_pairs=[],
        )
        pairs = _drawn_pairs(result, ss)
        assert pairs[0]["color"] == COLORS["e"]


class TestFigureBoundsIncludesRoutedLines:
    """MUST-FIX #2: an outward-routed PK-A line's points must extend the
    figure bounds, or its apex is invisibly clipped.
    """

    def test_line_point_outside_nucleotide_hull_extends_bounds(self) -> None:
        r = RNARenderer()
        r.set_coords([0.0, 10.0], [0.0, 0.0], 10.0)
        line = RoutedLine(i=0, j=1, points=[(0.0, 0.0), (500.0, 500.0), (10.0, 0.0)])
        xs, ys = _figure_bounds(r, [line])
        assert max(xs) >= 500.0
        assert max(ys) >= 500.0

    def test_no_lines_reduces_to_nucleotide_only_bounds(self) -> None:
        r = RNARenderer()
        r.set_coords([0.0, 10.0], [0.0, 0.0], 10.0)
        xs, ys = _figure_bounds(r, [])
        assert xs == list(r.xarray_)
        assert ys == list(r.yarray_)


class TestShiftRoutedLines:
    def test_shifts_every_point(self) -> None:
        line = RoutedLine(i=0, j=1, points=[(10.0, 20.0), (30.0, 40.0)])
        shifted = _shift_routed_lines([line], shift_x=5.0, shift_y=2.0)
        assert shifted[0].points == [(5.0, 18.0), (25.0, 38.0)]
        assert shifted[0].i == 0
        assert shifted[0].j == 1


class TestEndToEndDrawPseudoknot:
    """Regression: a pseudoknot no longer routes to the bare all-pairs-
    deleted circle fallback -- it produces a PNG with real crossing
    connectors, via the public `RNADrawer.draw` API.
    """

    def test_pseudoknot_renders_a_png(self) -> None:
        out = settings.Paths.UNITTEST_PATH + "test_renders/test_pseudoknot_h_type"
        drawer = RNADrawer()
        seq = "GCGCAUGC" * 4
        seq = seq[: len(H_TYPE)]
        fig = drawer.draw(ss=H_TYPE, seq=seq, filename=out)
        assert fig is not None
