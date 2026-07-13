import argparse
import warnings
from collections.abc import Callable

from rna_draw import render_rna, parameters, colorer
from rna_draw.colorer import *
from rna_draw.data import Data
from rna_draw.layout import LayoutResult, RoutedLine, layout_guaranteed, resolve_engine
from rna_draw.layout.base import params_at_node_r
from rna_draw.overlap import OverlapParams
from rna_draw.validate import check_routed_lines

# Figure-size model fit once from four reference (area, figsize) points -- see
# git history for the original per-render scipy.optimize.curve_fit call this
# replaces: areas = [3781, 126207, 1150472, 4286761], figsize = [25, 30, 35, 40],
# model a * (x - b) ** c. Precomputing makes figure sizing deterministic across
# scipy versions (the live fit was environment-sensitive, b ~= -15161 is an
# ill-conditioned local minimum) and drops the pandas/scipy runtime deps.
_FIGSIZE_A = 10.834742003709371
_FIGSIZE_B = -15161.400300804866
_FIGSIZE_C = 0.08504100128169866

# Pseudoknot crossing-connector colors (rna_draw.layout.pseudoknot, M3):
# PK-B in-plane connectors and PK-A routed lines are drawn in distinct,
# already-defined swatches (reuse before add) so they read apart from the
# ordinary gray ("e") nested pairs at a glance.
PK_CONNECTOR_COLOR = COLORS["o"]  # PK-B: in-plane straight crossing helix
PK_LINE_COLOR = COLORS["r"]  # PK-A: routed non-overlapping line


def get_parser():
    parser = argparse.ArgumentParser()
    # `required=False`: `-doc` bypasses -ss/-seq/etc. entirely (see
    # `__rna_draw_from_args`'s manual guard for the non-doc path).
    parser.add_argument("-ss", help="secondary structure in dot bracket notation", required=False)
    parser.add_argument("-seq", help="rna sequence", required=False)
    parser.add_argument("-out", help="output png file", required=False, default="secstruct")
    parser.add_argument(
        "-color_str",
        help="description of coloring, see docs for options",
        required=False,
    )
    parser.add_argument(
        "-render_type",
        help="scheme to color by: res_type,paired,motif,none",
        required=False,
    )
    parser.add_argument("-default_color", help="the color used when no other color is supplied")

    parser.add_argument("-data_str", help="data values by res seperated by ;", required=False)
    parser.add_argument("-data_file", help="path to data to color by", required=False)
    parser.add_argument("-data_palette", help="matplotlib color palette", required=False)
    parser.add_argument(
        "-data_vmin",
        help="min data value everything lower than this will be set to this value",
        required=False,
    )
    parser.add_argument(
        "-data_vmax",
        help="max data value everything greater than this will be set to this value",
        required=False,
    )
    parser.add_argument(
        "-data_ignore_restype",
        help="data values will be ignored for these restypes, e.g. G and U for DMS",
        required=False,
    )
    parser.add_argument(
        "-engine",
        default="auto",
        help=(
            "layout engine: auto|legacy|puzzler|production|constructive|"
            "native_puzzler|native_turtle|vienna_puzzler|turtle "
            "(vienna_puzzler/turtle need a -DRNA_DRAW_BUILD_ORACLE=ON build)"
        ),
        required=False,
    )
    parser.add_argument(
        "-doc",
        help="path to a .rnadoc.json document to load and draw (bypasses -ss/-seq/etc.)",
        required=False,
    )
    parser.add_argument(
        "-style",
        help="style preset for -doc: a registry name (e.g. 'default') or a path to .rnastyle.json",
        required=False,
    )
    parser.add_argument(
        "-save_doc",
        help="path to write a .rnadoc.json snapshot of what was just drawn",
        required=False,
    )
    return parser


def parse_args():
    parser = get_parser()
    args = parser.parse_args()
    return args


class RNADrawer(object):
    def __init__(self):
        self.__colorer = colorer.Colorer()
        self.__draw_params = parameters.DrawParameters()

    def draw(
        self,
        ss,
        seq=None,
        filename="secstruct",
        color_str=None,
        render_type=None,
        default_color=None,
        data=None,
        draw_params=None,
        engine="auto",
        on_layout=None,
    ):
        """Draw `ss`/`seq` to `filename`.png.

        Args:
            on_layout: Optional callback invoked with the raw
                `LayoutResult` right after layout, before any coordinate
                mutation (`RNARenderer.set_coords` shifts coords in
                place) -- lets a caller snapshot the pre-shift layout for
                persistence (`rna_draw.document_render.document_from_layout`)
                without changing this method's return type.
        """
        self.__setup(draw_params)

        if seq is None:
            seq = " " * len(ss)

        final_color_rbgs = self.__colorer.get_rgb_colors(
            seq, ss, color_str, data, render_type, default_color
        )

        return self.__render(
            seq, ss, final_color_rbgs, filename, self.__draw_params, engine, on_layout
        )

    def __setup(self, draw_params):
        if draw_params is not None:
            self.__draw_params = draw_params

    def __render(self, seq, ss, colors, filename, params, engine="auto", on_layout=None):
        r = render_rna.RNARenderer()

        gate = OverlapParams(node_r=params.NODE_R)  # target radius; adaptive search may shrink it
        result = layout_guaranteed(ss, engine=resolve_engine(engine), params=gate)
        _warn_if_flagged(result, params.NODE_R)
        if on_layout is not None:
            on_layout(result)
        pairs = _drawn_pairs(result, ss)
        # `set_coords` mutates result.x/y IN PLACE (shifts to a non-negative
        # bounding box); compute the same shift here first so any PK-A
        # `crossing_lines` points (computed in the ORIGINAL, unshifted
        # coordinate space) can be shifted to match (MUST-FIX #2).
        min_x = min(result.x) - result.node_r
        min_y = min(result.y) - result.node_r
        r.set_coords(result.x, result.y, result.node_r)
        lines = _shift_routed_lines(result.crossing_lines, min_x, min_y)

        r.ax.axis("off")
        _set_bounds_and_size(r, params, lines)
        r.draw(
            params.CELL_PADDING,
            params.CELL_PADDING,
            colors,
            pairs,
            seq,
            params.RENDER_IN_LETTERS,
        )
        r.draw_routed_lines(lines, params.CELL_PADDING, params.CELL_PADDING, PK_LINE_COLOR)
        r.fig.savefig(filename + ".png")
        return r.fig


def _warn_if_flagged(result: LayoutResult, target_node_r: float) -> None:
    """Surface a flagged/blind-spot layout as a warning (B2) -- never raise.

    `flagged` is a legitimate contract (PK/fallback tiers MAY be non-clean,
    see `LayoutResult`'s docstring); only the document-LOAD path
    (`document_render.draw_document`) raises on an unclean layout. This
    also runs the RC2 layered routed-line test (`validate.check_routed_lines`,
    reusing `pseudoknot.validate.polyline_is_clean`) at the RC1 SCALED
    params (`params_at_node_r(OverlapParams(node_r=target_node_r),
    result.node_r)`) so the pre-existing PK-A blind spot (a routed line
    crossing a disk/backbone/pair capsule) is visible in the ordinary
    batch-draw path too, without changing routing behavior.

    Args:
        result: The `layout_guaranteed` result, BEFORE `set_coords`
            mutates its coordinates.
        target_node_r: The `DrawParameters.NODE_R` this layout's adaptive
            search was gated against (the search ceiling).
    """
    if result.flagged:
        warnings.warn(
            f"layout for this structure is flagged (engine={result.engine_name}); "
            "it may contain overlaps",
            stacklevel=3,
        )
    pair_map = result.pair_map if result.pair_map is not None else [-1] * len(result.x)
    scaled = params_at_node_r(OverlapParams(node_r=target_node_r), result.node_r)
    dirty = check_routed_lines(result.x, result.y, pair_map, result.crossing_lines, scaled)
    if dirty:
        warnings.warn(
            "layout has PK-A routed line(s) crossing a nucleotide/backbone/pair "
            "capsule (pre-existing blind spot); rendered anyway",
            stacklevel=3,
        )


def _drawn_pairs(result, ss):
    """The base-pair connectors `__render` should draw.

    MUST-FIX #1: on the pseudoknot path (`result.pair_map` set), the
    max-nested extraction can retain some `[]{}<>` pairs and drop some
    `()` pairs, so the drawn pairs can differ from the plain `()`-only
    view -- `result.pair_map` (the FULL nested + PK-B augmented map) is
    used instead of `render_rna.get_pairmap_from_secstruct(ss)` whenever
    it is present. PK-B crossing connectors are colored distinctly.

    Args:
        result: The `layout_guaranteed` result.
        ss: The original dot-bracket string (used only on the non-pk path).

    Returns:
        `render_rna.RNARenderer.draw`'s `pairs` argument: one dict per
        pair with `from`/`to`/`p`/`color`.
    """
    pair_map = result.pair_map
    if pair_map is None:
        pair_map = render_rna.get_pairmap_from_secstruct(ss)
    crossing = {frozenset(pair) for pair in result.crossing_pairs}
    pairs = []
    for i, partner in enumerate(pair_map):
        if partner <= i:
            continue
        color = PK_CONNECTOR_COLOR if frozenset((i, partner)) in crossing else COLORS["e"]
        pairs.append({"from": i, "to": partner, "p": 1.0, "color": color})
    return pairs


def _shift_routed_lines(lines, shift_x, shift_y):
    """Translate every `RoutedLine`'s points by `-shift_x, -shift_y`.

    Mirrors the shift `RNARenderer.set_coords` applies to nucleotide
    coordinates in place, so a PK-A line stays aligned with the (now
    shifted) disks/pairs it connects.

    Args:
        lines: `RoutedLine`s in the ORIGINAL (unshifted) coordinate space.
        shift_x: Amount subtracted from every point's x.
        shift_y: Amount subtracted from every point's y.

    Returns:
        New `RoutedLine`s with every point shifted; `i`/`j` unchanged.
    """
    return [
        RoutedLine(
            i=line.i,
            j=line.j,
            points=[(px - shift_x, py - shift_y) for px, py in line.points],
        )
        for line in lines
    ]


def _figure_bounds(r, lines):
    """All x/y extents that must stay inside the drawn figure.

    Nucleotide positions plus any PK-A routed-line points -- MUST-FIX #2:
    an outward-routed line's ring/apex points can sit OUTSIDE the
    nucleotide hull, and must not be clipped by the figure bounds.

    Args:
        r: The `RNARenderer`, already `set_coords`-shifted.
        lines: `RoutedLine`s already shifted to match `r.xarray_`/`yarray_`.

    Returns:
        `(xs, ys)`: every x and y value the figure must cover.
    """
    xs = list(r.xarray_) + [px for line in lines for px, _ in line.points]
    ys = list(r.yarray_) + [py for line in lines for _, py in line.points]
    return xs, ys


def _set_bounds_and_size(r, params, lines):
    """Set the axes limits and figure size to cover every drawn point.

    A direct generalization of the original nucleotide-only bounds: with
    no `lines`, `_figure_bounds` reduces to exactly `r.xarray_`/`yarray_`,
    so the non-pseudoknot path is byte-for-byte unchanged.

    Args:
        r: The `RNARenderer`, already `set_coords`-shifted.
        params: `DrawParameters` (`CELL_PADDING` is the render offset).
        lines: `RoutedLine`s already shifted to match `r.xarray_`/`yarray_`.
    """
    xs, ys = _figure_bounds(r, lines)
    pad = params.CELL_PADDING
    r.ax.set_xlim([min(xs) + pad - 15, max(xs) + pad + 15])
    r.ax.set_ylim([min(ys) + pad - 15, max(ys) + pad + 15])

    if min(r.xarray_) == max(r.yarray_):
        return
    area = (max(xs) - min(xs)) * (max(ys) - min(ys))
    denom = _FIGSIZE_A * (area - _FIGSIZE_B) ** _FIGSIZE_C
    r.fig.set_size_inches((max(xs) - min(xs)) / denom, (max(ys) - min(ys)) / denom)


def _get_data_from_args(args):
    """Build a `Data` object from CLI-parsed `-data_*` flags, or `None`."""
    if args.data_str is not None or args.data_file is not None:
        return Data(
            args.data_str,
            args.data_file,
            args.data_palette,
            args.data_vmin,
            args.data_vmax,
            args.data_ignore_restype,
        )
    else:
        return None


def _get_render_type(render_type_name):
    """Map a `-render_type`/`Document.source.coloring.render_type` string
    to `colorer.RenderType` -- the single shared implementation both the
    CLI flag path (`__rna_draw_from_args`) and the document coloring
    bridge (`document_render.resolve_colors`) call, per the plan's concern
    (b) (one copy, no duplication).
    """
    if render_type_name is None:
        return None
    render_type_name = render_type_name.lower()
    if render_type_name == "res_type":
        return colorer.RenderType.RES_TYPE
    elif render_type_name == "paired":
        return colorer.RenderType.PAIRED
    elif render_type_name == "none":
        return None
    else:
        raise ValueError("unknown render type: " + render_type_name)


def _capture_layout(sink: dict) -> Callable[[LayoutResult], None]:
    """Build an `on_layout` callback snapshotting a layout before `set_coords`
    mutates its coordinates in place.

    Args:
        sink: A dict the callback stores the snapshot into (under
            `"result"`) -- the side channel `-save_doc` uses to get the
            pre-shift `LayoutResult` back out without changing
            `RNADrawer.draw`'s `Figure` return type.

    Returns:
        A callable suitable for `RNADrawer.draw(..., on_layout=...)`.
    """

    def _capture(result: LayoutResult) -> None:
        sink["result"] = LayoutResult(
            x=list(result.x),
            y=list(result.y),
            engine_name=result.engine_name,
            report=result.report,
            flagged=result.flagged,
            node_r=result.node_r,
            pair_map=list(result.pair_map) if result.pair_map is not None else None,
            crossing_pairs=list(result.crossing_pairs),
            crossing_lines=list(result.crossing_lines),
        )

    return _capture


def _write_save_doc(args, result: LayoutResult) -> None:
    """Write a `.rnadoc.json` snapshot of a just-drawn `-save_doc` layout."""
    from rna_draw.document_render import SourceInputs, document_from_layout

    source = SourceInputs(
        ss=args.ss,
        seq=args.seq or " " * len(args.ss),
        render_type=args.render_type,
        color_str=args.color_str,
        default_color=args.default_color,
        engine=args.engine,
    )
    document_from_layout(source, result).save(args.save_doc)


def _draw_from_document_args(args):
    """Load `-doc`/`-style` and draw via `document_render.draw_document`."""
    from rna_draw.document import Document
    from rna_draw.document_render import draw_document, resolve_preset

    doc = Document.load(args.doc)
    preset = resolve_preset(args.style)
    return draw_document(doc, preset=preset, filename=args.out)


def __rna_draw_from_args(args):
    if args.doc:
        return _draw_from_document_args(args)
    if args.ss is None:
        raise ValueError("either -ss or -doc is required")

    data = _get_data_from_args(args)
    render_type = _get_render_type(args.render_type)
    default_color = None
    if args.default_color is not None:
        default_color = colorer.parse_color_code(args.default_color)

    rd = RNADrawer()
    captured: dict = {}
    on_layout = _capture_layout(captured) if args.save_doc else None
    fig = rd.draw(
        args.ss,
        args.seq,
        args.out,
        args.color_str,
        render_type,
        default_color,
        data,
        engine=args.engine,
        on_layout=on_layout,
    )
    if args.save_doc:
        _write_save_doc(args, captured["result"])
    return fig


def _rna_draw_from_document(doc, style, kwargs):
    """Object-typed `rna_draw(doc=..., style=...)` path -- bypasses argparse
    entirely so a `Document`/`StylePreset` object can be passed directly.

    Args:
        doc: A `Document`, a path/registry-name string, or a raw dict.
        style: A `StylePreset`, registry name, path, dict, or `None`.
        kwargs: Remaining `rna_draw(**kwargs)` arguments; only `out` (the
            output filename stem) is currently consumed.

    Returns:
        The matplotlib `Figure` `document_render.draw_document` returns.
    """
    from rna_draw.document import Document
    from rna_draw.document_render import draw_document, resolve_preset

    if isinstance(doc, str):
        doc = Document.load(doc)
    elif isinstance(doc, dict):
        doc = Document.from_dict(doc)
    preset = resolve_preset(style)
    filename = kwargs.pop("out", "secstruct")
    return draw_document(doc, preset=preset, filename=filename)


def rna_draw(**kwargs):
    doc = kwargs.pop("doc", None)
    style = kwargs.pop("style", None)
    if doc is not None:
        return _rna_draw_from_document(doc, style, kwargs)

    parser = get_parser()
    args = []
    for k, v in kwargs.items():
        args.append("-" + k)
        args.append(v)
    args = parser.parse_args(args)
    return __rna_draw_from_args(args)


def main():
    args = parse_args()
    # Draw for its side effect (writes <out>.png). Deliberately return None:
    # the `rna_draw` console script runs `sys.exit(main())`, so returning the
    # matplotlib Figure would make Python print its repr and exit with code 1
    # on an otherwise-successful run. The library entry point `rna_draw(**kwargs)`
    # still returns the Figure for programmatic use.
    __rna_draw_from_args(args)


if __name__ == "__main__":
    main()
