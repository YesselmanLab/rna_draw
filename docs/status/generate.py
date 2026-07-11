"""Regenerate the layout-rebuild status report (HTML for the Claude Artifact).

Renders the before/after showcase structures with the current engines, then
fills `template.html` with the build-trend chart, milestone bars, algorithm
stages, showcase, and attempts table.

Usage (from repo root, in the ViennaRNA-enabled env):
    micromamba run -n py3 python docs/status/generate.py [OUTDIR]

Writes OUTDIR/status.html (default: docs/status/). To update the published
Artifact, add a row to BUILDS below when a new build moves the hard-set
numbers, rerun this, and republish the SAME status.html file path (keeps the
existing Artifact URL). Keep docs/layout_algorithm.md's progression in sync.

Artifact: https://claude.ai/code/artifact/eb56f105-7ace-427b-817b-ef9011f02fc7
"""

from __future__ import annotations

import base64
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

# --- build series: (label, short-key, clean %, clean count / 450, overlaps) ---
# APPEND a row here when a new build changes the hard-set clean rate/overlaps.
BUILDS = [
    ("Stock puzzler", "baseline", 30.4, 137, 2943),
    ("+ intersection clearance 1.5×", "clearance", 42.9, 193, 2008),
    ("+ per-structure escalating ladder", "escalation", 52.2, 235, 1896),
    ("+ rigid loop-inflation post-pass", "post-pass", 72.9, 328, 1739),
    ("+ cap / move-budget tuning", "tuning", 87.3, 393, 945),
]

# before/after showcase structures (must be in benchmarks/hard_set.json)
SHOWCASE = ["bpRNA_CRW_5583.dbn", "bpRNA_RFAM_35409.dbn", "bpRNA_CRW_10035.dbn"]

ATTEMPTS = [
    ("Portfolio: best of puzzler / naview / turtle", "The three engines fail on different structures, so per-structure best-of-three should beat any one.", "Clean count identical (137). naview and turtle never rescue a structure puzzler leaves dirty.", "dead-end"),
    ("Enable puzzler's allowFlipping resolver knob", "Flipping loops lets the resolver escape intersections rotation alone can't.", "Worse — 329 non-clean vs 313.", "dead-end"),
    ("Raise the config-change budget 25k → 1M", "The resolver is giving up because it runs out of iterations on big structures.", "Byte-identical output. The resolver converges; it is not budget-limited.", "dead-end"),
    ("Global loop-radius inflation", "Bigger loops give crowded branches more room.", "Produced NaN coordinates past 1.5× and slowed the resolver catastrophically.", "dead-end"),
    ("Scale puzzler's intersection clearance", "Diagnosis: puzzler resolves every intersection IT detects, but its clearance is more lenient than our checker, so it calls near-touches clean.", "2943 → 2008 overlaps; clean 30% → 43%. First lever that moved the number.", "win"),
    ("Per-structure escalating clearance ladder", "Different structures need different clearance; escalate only where still dirty.", "Clean 43% → 52% (finer 1.0–1.5× ladder).", "win"),
    ("Rigid loop-inflation post-pass (monotone)", "Diagnosis: the residual overlaps are small-loop crowding of unpaired bases, not branch collisions. Inflate the offending loop; accept only if the checker reports strictly fewer overlaps.", "Clean 52% → 73%. The single biggest jump.", "big-win"),
    ("Raise post-pass cap 6 → 20, move budget 12 → 40", "Most 7–20-overlap structures are the same crowding the inflation move clears; the cap was skipping them.", "Clean 73% → 87%.", "win"),
    ("Productionize as the default engine", "Make the whole stack rna_draw's real output, behind the honest clean-or-flagged contract.", "default_engine() renders a stock-dirty 684nt structure clean in 0.39s; empty-loop inputs flag safely.", "shipped"),
    ("Spatial-hash O(L²) → O(L) for long chords", "Fallback renders took minutes because diagonal base-pair chords blew up the hash's bounding-box insertion.", "Checker ~14× faster on the fallback circle; byte-identical output (proven vs brute force).", "shipped"),
]

STAGES = [
    ("1", "Clearance escalation", "Lay out with RNApuzzler at an increasing intersection-clearance factor (1.0 → 1.5×), scaling puzzler's own overlap-detection margin so it resolves the near-touches our stricter checker flags. Keep the first clean result.", "in-process · hang-safe ≤1.5×"),
    ("2", "Rigid loop-inflation post-pass", "For whatever crowding remains, inflate the offending loop — spread its unpaired bases, translate child branches rigidly so helices stay straight. Apply a move only if the frozen checker reports strictly fewer overlaps. Never worse; wall-clock bounded.", "monotone · engine-agnostic"),
    ("3", "Checker-gated radius + guaranteed fallback", "Render at the largest disk radius that passes the checker. If no conventional layout clears, fall back to a guaranteed-clean circle, flagged. Every output is checker-clean or explicitly flagged — never a silent overlap.", "honest contract"),
]

VMAP = {"dead-end": ("dead end", "v-bad"), "win": ("win", "v-win"), "big-win": ("big win", "v-big"), "shipped": ("shipped", "v-ship")}


def render_showcase() -> dict:
    """Render each showcase structure stock-dirty vs production-clean; return base64 data URIs."""
    sys.path.insert(0, str(REPO))
    from PIL import Image  # noqa: PLC0415
    from rna_draw import render_rna  # noqa: PLC0415
    from rna_draw.colorer import COLORS, Colorer  # noqa: PLC0415
    from rna_draw.layout.pipeline import layout_guaranteed, resolve_engine  # noqa: PLC0415
    from rna_draw.layout.vienna import ViennaPuzzlerEngine  # noqa: PLC0415
    from rna_draw.overlap import check_overlaps  # noqa: PLC0415
    from rna_draw.parameters import DrawParameters  # noqa: PLC0415

    hard = {s["name"]: s for s in json.loads((REPO / "benchmarks/hard_set.json").read_text())}
    params, colorer = DrawParameters(), Colorer()
    tmp = HERE / "_renders"
    tmp.mkdir(exist_ok=True)

    def draw(ss, x, y, node_r, path):
        r = render_rna.RNARenderer()
        pm = render_rna.get_pairmap_from_secstruct(ss)
        pairs = [{"from": i, "to": pm[i], "p": 1.0, "color": COLORS["e"]} for i in range(len(pm)) if pm[i] > i]
        seq = " " * len(ss)
        colors = colorer.get_rgb_colors(seq, ss, None, None, None, None)
        r.set_coords(x, y, node_r)
        r.ax.axis("off")
        r.ax.set_xlim([min(r.xarray_) + 25, max(r.xarray_) + 55])
        r.ax.set_ylim([min(r.yarray_) + 25, max(r.yarray_) + 55])
        area = (max(r.xarray_) - min(r.xarray_)) * (max(r.yarray_) - min(r.yarray_))
        denom = 10.834742003709371 * (area + 15161.400300804866) ** 0.08504100128169866
        r.fig.set_size_inches(max(2, (max(r.xarray_) - min(r.xarray_)) / denom), max(2, (max(r.yarray_) - min(r.yarray_)) / denom))
        r.draw(params.CELL_PADDING, params.CELL_PADDING, colors, pairs, seq, False)
        r.fig.savefig(str(path), dpi=85, bbox_inches="tight")

    def encode(path):
        im = Image.open(path).convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        bg.alpha_composite(im)
        im = bg.convert("RGB")
        w, h = im.size
        s = min(1.0, 640 / w)
        im = im.resize((int(w * s), int(h * s)), Image.LANCZOS)
        im.save(path, optimize=True)
        return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()

    out = {}
    for name in SHOWCASE:
        ss = hard[name]["structure"]
        pm = render_rna.get_pairmap_from_secstruct(ss)
        sx, sy = ViennaPuzzlerEngine().layout(ss)
        s_ov = check_overlaps(sx, sy, pm).num_overlaps
        draw(ss, sx, sy, params.NODE_R, tmp / f"{name}_s.png")
        res = layout_guaranteed(ss, engine=resolve_engine("production"))
        draw(ss, res.x, res.y, res.node_r, tmp / f"{name}_p.png")
        out[name] = {"len": len(ss), "stock_ov": s_ov, "stock": encode(tmp / f"{name}_s.png"), "prod": encode(tmp / f"{name}_p.png")}
    return out


def chart_svg() -> str:
    """A build-over-time chart: overlaps (area, left axis) + clean rate (line, right axis)."""
    w, h, pl, pr, pt, pb = 720, 344, 60, 56, 30, 60
    pw, ph = w - pl - pr, h - pt - pb
    n = len(BUILDS)
    ov = [b[4] for b in BUILDS]
    cl = [b[2] for b in BUILDS]
    keys = [b[1] for b in BUILDS]
    ovmax = max(3000, max(ov))
    xs = [pl + (i / (n - 1)) * pw for i in range(n)]
    yov = lambda v: pt + (1 - v / ovmax) * ph  # noqa: E731
    ycl = lambda v: pt + (1 - v / 100) * ph  # noqa: E731
    bottom = pt + ph
    p = []
    for gv in range(0, ovmax + 1, 1000):
        y = yov(gv)
        p.append(f'<line class="grid" x1="{pl}" y1="{y:.1f}" x2="{pl+pw}" y2="{y:.1f}"/>')
        p.append(f'<text class="ytick ytick-ov" x="{pl-10}" y="{y+4:.1f}" text-anchor="end">{gv}</text>')
    for cv in (0, 25, 50, 75, 100):
        p.append(f'<text class="ytick ytick-cl" x="{pl+pw+10}" y="{ycl(cv)+4:.1f}" text-anchor="start">{cv}%</text>')
    p.append(f'<text class="axlabel" x="{pl-38}" y="{pt-12:.1f}" text-anchor="start">OVERLAPS</text>')
    p.append(f'<text class="axlabel" x="{pl+pw+2:.1f}" y="{pt-12:.1f}" text-anchor="start">CLEAN</text>')
    area = f"M {xs[0]:.1f},{yov(ov[0]):.1f} " + " ".join(f"L {xs[i]:.1f},{yov(ov[i]):.1f}" for i in range(1, n))
    area += f" L {xs[-1]:.1f},{bottom:.1f} L {xs[0]:.1f},{bottom:.1f} Z"
    p.append(f'<path class="area" d="{area}"/>')
    p.append('<polyline class="line-cl" points="' + " ".join(f"{xs[i]:.1f},{ycl(cl[i]):.1f}" for i in range(n)) + '"/>')
    p.append('<polyline class="line-ov" points="' + " ".join(f"{xs[i]:.1f},{yov(ov[i]):.1f}" for i in range(n)) + '"/>')
    for i in range(n):
        p.append(f'<circle class="dot-cl" cx="{xs[i]:.1f}" cy="{ycl(cl[i]):.1f}" r="3.5"/>')
        p.append(f'<circle class="dot-ov" cx="{xs[i]:.1f}" cy="{yov(ov[i]):.1f}" r="4.5"/>')
        p.append(f'<text class="val" x="{xs[i]:.1f}" y="{yov(ov[i])-12:.1f}" text-anchor="middle">{ov[i]}</text>')
        p.append(f'<text class="xtick" x="{xs[i]:.1f}" y="{bottom+22:.1f}" text-anchor="middle">{keys[i]}</text>')
        p.append(f'<text class="axlabel" x="{xs[i]:.1f}" y="{bottom+37:.1f}" text-anchor="middle">build {i+1}</text>')
    aria = f"Total overlaps falling from {ov[0]} to {ov[-1]} across {n} builds; clean rate rising from {cl[0]:.0f} to {cl[-1]:.0f} percent"
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{aria}">' + "".join(p) + "</svg>"


def build(outdir: Path) -> Path:
    imgs = render_showcase()
    ms = ""
    for i, (label, _key, pct, clean, ov) in enumerate(BUILDS):
        delta = "" if i == 0 else f'<span class="ms-delta">+{round(pct-BUILDS[i-1][2],1)}</span>'
        ms += (f'<div class="ms-row"><div class="ms-label"><span class="ms-name">{label}</span></div>'
               f'<div class="ms-track"><div class="ms-fill" style="width:{round(pct,1)}%"></div>'
               f'<span class="ms-pct">{pct}%{delta}</span></div>'
               f'<div class="ms-meta"><span>{clean}/450 clean</span><span class="ms-ov">{ov} overlaps</span></div></div>')
    show = ""
    for name in SHOWCASE:
        d = imgs[name]
        disp = name.replace(".dbn", "")
        show += (f'<figure class="shot"><figcaption class="shot-cap"><span class="shot-name">{disp}</span>'
                 f'<span class="shot-len">{d["len"]} nt</span></figcaption><div class="shot-pair">'
                 f'<div class="shot-cell"><div class="shot-tag tag-bad">stock · {d["stock_ov"]} overlaps</div>'
                 f'<img loading="lazy" alt="{disp} stock puzzler, {d["stock_ov"]} overlaps" src="{d["stock"]}"></div>'
                 f'<div class="shot-arrow" aria-hidden="true">→</div>'
                 f'<div class="shot-cell"><div class="shot-tag tag-good">production · 0 overlaps</div>'
                 f'<img loading="lazy" alt="{disp} production engine, overlap-free" src="{d["prod"]}"></div></div></figure>')
    att = ""
    for i, (lever, hyp, res, verd) in enumerate(ATTEMPTS, 1):
        vt, vc = VMAP[verd]
        att += (f'<tr><td class="a-num">{i:02d}</td><td class="a-lever"><div class="a-lever-t">{lever}</div>'
                f'<div class="a-hyp">{hyp}</div></td><td class="a-res">{res}</td>'
                f'<td class="a-verd"><span class="verd {vc}">{vt}</span></td></tr>')
    stg = ""
    for num, title, body, tag in STAGES:
        stg += (f'<div class="stage"><div class="stage-n">{num}</div><div class="stage-b"><h3>{title}</h3>'
                f'<p>{body}</p><span class="stage-tag">{tag}</span></div></div>')
    html = (HERE / "template.html").read_text()
    html = (html.replace("<!--MILES-->", ms).replace("<!--SHOW-->", show)
            .replace("<!--ATT-->", att).replace("<!--STAGES-->", stg).replace("<!--CHART-->", chart_svg()))
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "status.html"
    out.write_text(html)
    return out


if __name__ == "__main__":
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
    written = build(dest)
    print(f"wrote {written} ({written.stat().st_size // 1024} KB)")
