"""Tests for the P1 persistence layer: `Document`/`StylePreset` save/load,
the layered never-silent load gate (RC1 scaled params, RC2 reused routed-line
validator), and `draw.py`'s B2 flagged-warning surfacing.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from rna_draw import parameters
from rna_draw.document import Document
from rna_draw.document_render import (
    SourceInputs,
    document_from_layout,
    draw_document,
    relayout,
    resolve_colors,
    resolve_preset,
    resolve_style,
)
from rna_draw.document_schema import ColoringIntent, DataIntent
from rna_draw.draw import RNADrawer
from rna_draw.layout import RoutedLine, layout_guaranteed
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.style import StylePreset, default_preset
from rna_draw.validate import LayoutOverlapError, check_routed_lines, never_silent_gate

NESTED_SS = "((((....))))"
NESTED_SEQ = "GGGGAAAACCCC"
# Same H-type pseudoknot used throughout `test_pseudoknot_*.py`; produces
# real PK-A `crossing_lines` at the default target radius (verified below).
PK_SS = "((((....[[[[....))))....]]]]"
# A structure whose adaptive search shrinks `node_r` below the 10.0 target
# AND produces a routed PK-A line -- exercises RC1 (scaled half-widths) and
# RC2 (reused `polyline_is_clean`) together on reload.
PK_SHRINKING_SS = "..<<<...[..(..>>>.]...{...).........}"

# The PK dot-brackets embedded in `test_pseudoknot_engine.py`/
# `test_pseudoknot_pipeline.py`/`test_pseudoknot_extraction.py`/
# `test_pseudoknot_proximity.py` -- the plan's B1 corpus.
PK_CORPUS = [
    PK_SS,
    "((([[[)))]]]",
    ".<<<<..[[(>>>>]]....{{{{).....}}}}",
    "((([[[[.....<<{]]]]))).>>}",
    "(((<<<)))>>>",
    "((((..[[..))))((((..]]..))))",
    PK_SHRINKING_SS,
    "((((...((((....))))...[[[[....))))....]]]]",
]


class FakeShrinkingEngine:
    """A `LayoutEngine` whose fixed coords force the adaptive search to
    shrink `node_r` below the target -- mirrors
    `test_layout_pipeline.TestAdaptiveRenderRadius`'s fixture.
    """

    name = "fake_shrink"

    # nt1 sits perpendicular (never near-touches); nt0/nt2 are 18 units
    # apart -- overlapping at the default node_r=10 (required 20 > 18) but
    # clean at node_r=9 (required 18 == 18).
    X = [0.0, 0.0, 18.0]
    Y = [0.0, 200.0, 0.0]

    def layout(self, secstruct: str) -> tuple[list[float], list[float]]:
        return list(self.X), list(self.Y)


def _build_document(ss: str, seq: str | None = None) -> Document:
    """Build a `Document` from a real `layout_guaranteed` result."""
    seq = seq or " " * len(ss)
    result = layout_guaranteed(ss, params=OverlapParams())
    source = SourceInputs(ss=ss, seq=seq)
    return document_from_layout(source, result)


class TestRoundTripFidelity:
    def test_nested_round_trip_is_idempotent(self, tmp_path: Path) -> None:
        doc = _build_document(NESTED_SS, NESTED_SEQ)
        path1 = tmp_path / "a.rnadoc.json"
        path2 = tmp_path / "b.rnadoc.json"
        doc.save(path1)
        loaded = Document.load(path1)
        loaded.save(path2)
        assert path1.read_text() == path2.read_text()

    def test_nested_coords_survive_exactly(self, tmp_path: Path) -> None:
        doc = _build_document(NESTED_SS, NESTED_SEQ)
        assert doc.derived.layout is not None
        path = tmp_path / "a.rnadoc.json"
        doc.save(path)
        loaded = Document.load(path)
        assert loaded.derived.layout is not None
        assert loaded.derived.layout.coords == doc.derived.layout.coords

    def test_pseudoknot_round_trip_is_idempotent(self, tmp_path: Path) -> None:
        doc = _build_document(PK_SS)
        path1 = tmp_path / "a.rnadoc.json"
        path2 = tmp_path / "b.rnadoc.json"
        doc.save(path1)
        Document.load(path1).save(path2)
        assert path1.read_text() == path2.read_text()

    def test_pseudoknot_crossing_data_survives_exactly(self, tmp_path: Path) -> None:
        doc = _build_document(PK_SS)
        assert doc.derived.layout is not None
        assert doc.derived.layout.crossing_pairs or doc.derived.layout.crossing_lines
        path = tmp_path / "a.rnadoc.json"
        doc.save(path)
        loaded = Document.load(path)
        assert loaded.derived.layout is not None
        assert loaded.derived.layout.crossing_pairs == doc.derived.layout.crossing_pairs
        assert loaded.derived.layout.crossing_lines == doc.derived.layout.crossing_lines


class TestUnknownFieldPassthrough:
    def test_extra_layout_field_survives_load_and_save(self, tmp_path: Path) -> None:
        doc = _build_document(NESTED_SS, NESTED_SEQ)
        path = tmp_path / "a.rnadoc.json"
        doc.save(path)

        raw = json.loads(path.read_text())
        raw["derived"]["layout"]["edited_indices"] = [1, 2]
        path.write_text(json.dumps(raw, indent=2))

        loaded = Document.load(path)
        path2 = tmp_path / "b.rnadoc.json"
        loaded.save(path2)

        resaved = json.loads(path2.read_text())
        assert resaved["derived"]["layout"]["edited_indices"] == [1, 2]


class TestNeverSilentOnLoad:
    def test_unmodified_pk_doc_loads_without_raising(self, tmp_path: Path) -> None:
        """Concern (a): an UNMODIFIED pipeline-produced PK doc must load
        clean -- guards the RC1 gate-consistency regression directly.
        """
        doc = _build_document(PK_SS)
        path = tmp_path / "a.rnadoc.json"
        doc.save(path)
        fig = draw_document(Document.load(path), filename=str(tmp_path / "out"))
        assert fig is not None

    def test_hand_edited_overlapping_coords_raises(self, tmp_path: Path) -> None:
        doc = _build_document(NESTED_SS, NESTED_SEQ)
        assert doc.derived.layout is not None
        coords = list(doc.derived.layout.coords)
        # Stack nucleotide 0 directly on top of nucleotide 6 (both disks at
        # the same point is an unmistakable disk-disk overlap).
        coords[0] = coords[6]
        doc.derived.layout.coords = coords

        with pytest.raises(LayoutOverlapError):
            draw_document(doc, filename=str(tmp_path / "out"))

    def test_routed_line_overlap_raises_but_frozen_checker_alone_passes(
        self, tmp_path: Path
    ) -> None:
        doc = _build_document(PK_SS)
        layout = doc.derived.layout
        assert layout is not None
        line = next(line for line in layout.crossing_lines if len(line.points) >= 2)
        midpoint = (
            (line.points[0][0] + line.points[1][0]) / 2,
            (line.points[0][1] + line.points[1][1]) / 2,
        )
        # Move some OTHER nucleotide (not the line's own endpoints) exactly
        # onto the routed line's segment.
        victim = next(i for i in range(len(layout.coords)) if i not in (line.i, line.j))
        coords = list(layout.coords)
        coords[victim] = midpoint
        layout.coords = coords

        scaled_report = check_overlaps(
            [pt[0] for pt in coords],
            [pt[1] for pt in coords],
            layout.pair_map,
            OverlapParams(node_r=layout.node_r),
        )
        assert scaled_report.passed, "the frozen disk/capsule checker alone must stay clean"

        with pytest.raises(LayoutOverlapError):
            draw_document(doc, filename=str(tmp_path / "out"))

    def test_check_routed_lines_direct_unit(self) -> None:
        line = RoutedLine(i=0, j=1, points=[(0.0, 0.0), (100.0, 0.0)])
        params = OverlapParams()
        # A disk sitting on the line's own endpoint is excluded.
        clean = check_routed_lines([0.0, 100.0], [0.0, 0.0], [-1, -1], [line], params)
        assert clean == []
        # A disk placed ON the line's interior segment is a real overlap.
        dirty = check_routed_lines(
            [0.0, 100.0, 50.0], [0.0, 0.0, 0.0], [-1, -1, -1], [line], params
        )
        assert dirty == [line]

    def test_advisory_cache_override_does_not_win(self, tmp_path: Path) -> None:
        doc = _build_document(NESTED_SS, NESTED_SEQ)
        assert doc.derived.layout is not None
        coords = list(doc.derived.layout.coords)
        coords[0] = coords[6]
        doc.derived.layout.coords = coords
        doc.derived.layout.checker = {"verdict": "passed", "node_r": doc.derived.layout.node_r}

        with pytest.raises(LayoutOverlapError):
            draw_document(doc, filename=str(tmp_path / "out"))


class TestRC1ScaledParamsRoundTrip:
    """RC1: the gate MUST use the pipeline's own scaled half-widths, not
    the library default, or a shrunk-`node_r` layout falsely raises.
    """

    def test_routed_line_pk_doc_round_trip_never_raises(self, tmp_path: Path) -> None:
        doc = _build_document(PK_SS)
        assert doc.derived.layout is not None
        assert doc.derived.layout.crossing_lines
        path = tmp_path / "a.rnadoc.json"
        doc.save(path)
        fig = draw_document(Document.load(path), filename=str(tmp_path / "out"))
        assert fig is not None

    def test_node_r_shrinking_structure_round_trip_never_raises(self, tmp_path: Path) -> None:
        result = layout_guaranteed("...", engine=FakeShrinkingEngine(), params=OverlapParams())
        assert result.node_r < 10.0  # confirms the adaptive search actually shrank it

        source = SourceInputs(ss="...", seq="AAA")
        doc = document_from_layout(source, result)
        path = tmp_path / "a.rnadoc.json"
        doc.save(path)
        fig = draw_document(Document.load(path), filename=str(tmp_path / "out"))
        assert fig is not None

    def test_pk_shrinking_structure_round_trip_never_raises(self, tmp_path: Path) -> None:
        """A real pseudoknot-tier structure that both shrinks `node_r` AND
        carries a routed PK-A line -- RC1 and RC2 exercised together.
        """
        doc = _build_document(PK_SHRINKING_SS)
        assert doc.derived.layout is not None
        assert doc.derived.layout.node_r < 10.0
        assert doc.derived.layout.crossing_lines
        path = tmp_path / "a.rnadoc.json"
        doc.save(path)
        fig = draw_document(Document.load(path), filename=str(tmp_path / "out"))
        assert fig is not None

    def test_gate_uses_scaled_not_default_params(self) -> None:
        """Direct unit proof of RC1's arithmetic: a disk sitting between the
        SCALED pair-capsule clearance (clean) and the library-DEFAULT
        clearance (a false overlap) at a shrunk `node_r=8.0` (target 10.0,
        ratio 0.75 -> half-width 6.0 scaled vs. the unscaled default 7.5).
        """
        x = [-50.0, 50.0, 0.0]
        y = [0.0, 0.0, 14.5]  # required: scaled 8+6.0=14.0 (clean); default 8+7.5=15.5 (overlap)
        pair_map = [1, 0, -1]
        node_r, target = 8.0, 10.0

        clean, report, dirty = never_silent_gate(x, y, pair_map, [], node_r, target)
        assert clean and report.passed and not dirty

        default_report = check_overlaps(x, y, pair_map, OverlapParams(node_r=node_r))
        assert not default_report.passed, "unscaled default params falsely flag the shrunk layout"


class TestB2FlaggedWarning:
    def test_flagged_layout_warns_on_normal_draw(self, tmp_path: Path) -> None:
        rd = RNADrawer()
        with pytest.warns(UserWarning):
            rd.draw(".((A)).", filename=str(tmp_path / "flagged_out"))


class TestDefaultPresetParity:
    def test_matches_draw_parameters_field_for_field(self) -> None:
        dp = default_preset().to_draw_parameters()
        ref = parameters.DrawParameters()
        assert dp.NODE_R == ref.NODE_R
        assert dp.PRIMARY_SPACE == ref.PRIMARY_SPACE
        assert dp.PAIR_SPACE == ref.PAIR_SPACE
        assert dp.CELL_PADDING == ref.CELL_PADDING
        assert dp.TEXT_SIZE == ref.TEXT_SIZE
        assert dp.RENDER_IN_LETTERS == ref.RENDER_IN_LETTERS
        assert dp.output_width == ref.output_width
        assert dp.output_height == ref.output_height

    def test_resolve_preset_default_name(self) -> None:
        preset = resolve_preset("default")
        assert preset.name == "default"
        assert preset is not resolve_preset("default")  # fresh instance each call


class TestBackCompat:
    def test_rna_draw_kwargs_shim_still_works(self, tmp_path: Path) -> None:
        import rna_draw as rd

        fig = rd.rna_draw(ss="(....)", seq="CUUCGG", out=str(tmp_path / "test_1"))
        assert fig is not None

        fig = rd.rna_draw(
            ss="(....)", seq="CUUCGG", out=str(tmp_path / "test_2"), color_str="1-6:r"
        )
        assert fig is not None

        fig = rd.rna_draw(
            ss="(....)",
            seq="CUUCGG",
            out=str(tmp_path / "test_3"),
            data_str="0;1;2;3;4;5",
        )
        assert fig is not None


class TestKwargsShimObjectPath:
    def test_rna_draw_doc_kwarg_accepts_a_path_string(self, tmp_path: Path) -> None:
        import rna_draw as rd

        doc = _build_document(NESTED_SS, NESTED_SEQ)
        doc_path = tmp_path / "a.rnadoc.json"
        doc.save(doc_path)

        fig = rd.rna_draw(doc=str(doc_path), out=str(tmp_path / "out1"))
        assert fig is not None

    def test_rna_draw_doc_kwarg_accepts_an_inline_dict(self, tmp_path: Path) -> None:
        import rna_draw as rd

        doc = _build_document(NESTED_SS, NESTED_SEQ)
        fig = rd.rna_draw(doc=doc.to_dict(), style="default", out=str(tmp_path / "out2"))
        assert fig is not None

    def test_rna_draw_doc_kwarg_accepts_a_document_instance(self, tmp_path: Path) -> None:
        import rna_draw as rd

        doc = _build_document(NESTED_SS, NESTED_SEQ)
        fig = rd.rna_draw(doc=doc, out=str(tmp_path / "out3"))
        assert fig is not None


class TestCliDocAndSaveDocFlags:
    """Exercises the true argparse `-doc`/`-style`/`-save_doc` path
    (`__rna_draw_from_args`/`_draw_from_document_args`/`_write_save_doc`),
    distinct from `rna_draw(doc=..., style=...)`'s object-kwargs shortcut.
    """

    def test_save_doc_flag_writes_a_document(self, tmp_path: Path, monkeypatch) -> None:
        from rna_draw import draw

        out_stem = str(tmp_path / "cli_out")
        doc_path = str(tmp_path / "cli.rnadoc.json")
        argv = [
            "rna_draw",
            "-ss",
            "(....)",
            "-seq",
            "CUUCGG",
            "-out",
            out_stem,
            "-save_doc",
            doc_path,
        ]
        monkeypatch.setattr("sys.argv", argv)
        draw.main()

        saved = Document.load(doc_path)
        assert saved.derived.layout is not None

    def test_doc_flag_draws_from_a_saved_document(self, tmp_path: Path, monkeypatch) -> None:
        from rna_draw import draw

        doc = _build_document(NESTED_SS, NESTED_SEQ)
        doc_path = tmp_path / "a.rnadoc.json"
        doc.save(doc_path)
        out_stem = str(tmp_path / "cli_doc_out")
        argv = ["rna_draw", "-doc", str(doc_path), "-style", "default", "-out", out_stem]
        monkeypatch.setattr("sys.argv", argv)
        draw.main()

        assert Path(out_stem + ".png").exists()

    def test_missing_ss_and_doc_raises(self) -> None:
        from rna_draw import draw

        # `__rna_draw_from_args` is module-private (leading double-underscore
        # at module scope, NOT name-mangled there); `getattr` by string
        # avoids this test class's OWN mangling of a literal `__`-prefixed
        # identifier.
        rna_draw_from_args = getattr(draw, "__rna_draw_from_args")
        args = draw.get_parser().parse_args([])
        with pytest.raises(ValueError):
            rna_draw_from_args(args)


class TestB1CorpusMeasurement:
    def test_corpus_measurement_is_well_defined(self, capsys: pytest.CaptureFixture) -> None:
        """Measurement, not a gate: count how many corpus structures have a
        `report.passed` frozen checker yet the RC2 routed-line test still
        finds a witness at the SAME scaled params the pipeline used.
        """
        target = default_preset().layout_defaults.node_r
        count = 0
        for ss in PK_CORPUS:
            result = layout_guaranteed(ss, params=OverlapParams(node_r=target))
            if not result.report.passed:
                continue
            _, _, dirty = never_silent_gate(
                result.x,
                result.y,
                result.pair_map or [-1] * len(result.x),
                result.crossing_lines,
                result.node_r,
                target,
            )
            if dirty:
                count += 1
        print(
            f"B1 corpus measurement: {count}/{len(PK_CORPUS)} clean-checker structures "
            f"still have a routed-line blind-spot witness"
        )
        assert count >= 0


class TestStylePresetRoundTrip:
    def test_save_load_round_trip(self, tmp_path: Path) -> None:
        preset = default_preset()
        preset.name = "custom"
        path = tmp_path / "a.rnastyle.json"
        preset.save(path)
        loaded = StylePreset.load(path)
        assert loaded.name == "custom"
        assert loaded.to_draw_parameters().NODE_R == preset.to_draw_parameters().NODE_R
        assert loaded.palette == preset.palette

    def test_load_wrong_schema_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.rnastyle.json"
        path.write_text(json.dumps({"schema": "not/a/style", "version": 1}))
        with pytest.raises(ValueError):
            StylePreset.load(path)


class TestResolvePresetBranches:
    def test_passes_through_an_existing_preset_instance(self) -> None:
        preset = default_preset()
        assert resolve_preset(preset) is preset

    def test_builds_from_an_inline_dict(self) -> None:
        preset = resolve_preset(default_preset().to_dict())
        assert preset.name == "default"

    def test_loads_from_a_path(self, tmp_path: Path) -> None:
        path = tmp_path / "a.rnastyle.json"
        default_preset().save(path)
        preset = resolve_preset(str(path))
        assert preset.name == "default"


class TestResolveStyleOverrides:
    def test_style_overrides_render_in_letters(self) -> None:
        doc = _build_document(NESTED_SS, NESTED_SEQ)
        doc.source.style_overrides = {"render_in_letters": True}
        draw_params, _ = resolve_style(default_preset(), doc)
        assert draw_params.RENDER_IN_LETTERS is True


class TestResolveColorsWithDataIntent:
    def test_data_intent_colors_are_resolved(self) -> None:
        doc = _build_document(NESTED_SS, NESTED_SEQ)
        doc.source.coloring = ColoringIntent(
            data=DataIntent(values=[float(i) for i in range(len(NESTED_SS))], palette="Reds")
        )
        colors = resolve_colors(doc, default_preset())
        assert len(colors) == len(NESTED_SS)

    def test_data_intent_round_trips(self, tmp_path: Path) -> None:
        doc = _build_document(NESTED_SS, NESTED_SEQ)
        doc.source.coloring = ColoringIntent(
            data=DataIntent(values=[1.0, 2.0, 3.0], vmin=0.0, vmax=5.0, ignore_restype="GU")
        )
        path = tmp_path / "a.rnadoc.json"
        doc.save(path)
        loaded = Document.load(path)
        assert loaded.source.coloring.data == doc.source.coloring.data


class TestRelayoutAndSourceOnlyDocument:
    def test_relayout_replaces_the_layout_band(self) -> None:
        doc = _build_document(NESTED_SS, NESTED_SEQ)
        original_node_r = doc.derived.layout.node_r  # type: ignore[union-attr]
        relaid = relayout(doc, default_preset())
        assert relaid.derived.layout is not None
        assert relaid.derived.layout.node_r == original_node_r

    def test_draw_document_relayouts_a_source_only_document(self, tmp_path: Path) -> None:
        from rna_draw.document import DerivedBand, SourceBand, StructureIntent

        source = SourceBand(structure=StructureIntent(ss=NESTED_SS, seq=NESTED_SEQ))
        doc = Document(source=source, derived=DerivedBand(layout=None))
        fig = draw_document(doc, filename=str(tmp_path / "out"))
        assert fig is not None
