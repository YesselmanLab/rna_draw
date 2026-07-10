"""Integration tests for the `-engine` CLI flag and `RNARenderer.set_coords`.

Covers the render-side wiring introduced in M3: coord injection into
`RNARenderer`, the `-engine` argparse flag, and end-to-end PNG output for
both an explicit engine and the default `auto` selection.
"""

from __future__ import annotations

from pathlib import Path

import rna_draw as rd
from rna_draw.draw import get_parser
from rna_draw.render_rna import RNARenderer


class TestSetCoords:
    """`RNARenderer.set_coords` injects externally computed coordinates."""

    def test_sets_shifted_arrays_size_and_node_r(self) -> None:
        renderer = RNARenderer()
        renderer.set_coords([0.0, 20.0, 40.0], [0.0, 0.0, 0.0], 10.0)

        assert renderer.NODE_R == 10.0
        assert min(renderer.xarray_) >= 0.0
        assert min(renderer.yarray_) >= 0.0
        assert renderer.size_ == [60.0, 20.0]

    def test_shift_preserves_relative_spacing(self) -> None:
        renderer = RNARenderer()
        renderer.set_coords([5.0, 25.0], [-3.0, -3.0], 10.0)

        assert renderer.xarray_[1] - renderer.xarray_[0] == 20.0
        assert renderer.yarray_[1] - renderer.yarray_[0] == 0.0


class TestEngineFlagParsing:
    """The `-engine` CLI flag defaults to `auto` and accepts overrides."""

    def test_default_engine_is_auto(self) -> None:
        args = get_parser().parse_args(["-ss", "...."])
        assert args.engine == "auto"

    def test_explicit_legacy_engine_flag(self) -> None:
        args = get_parser().parse_args(["-ss", "....", "-engine", "legacy"])
        assert args.engine == "legacy"


class TestRnaDrawEngineSelection:
    """`rna_draw(...)` writes a PNG for both explicit and default engines."""

    def test_legacy_engine_writes_png(self, tmp_path: Path) -> None:
        out = tmp_path / "legacy_out"
        rd.rna_draw(
            ss="((((....))))",
            seq="GGGGAAAACCCC",
            out=str(out),
            engine="legacy",
        )
        assert (tmp_path / "legacy_out.png").exists()

    def test_default_auto_engine_writes_png(self, tmp_path: Path) -> None:
        out = tmp_path / "auto_out"
        rd.rna_draw(ss="((((....))))", seq="GGGGAAAACCCC", out=str(out))
        assert (tmp_path / "auto_out.png").exists()
