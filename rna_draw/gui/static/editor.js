// rna_draw interactive editor view (anywidget, dependency-free ESM).
//
// The kernel is authoritative. This view only: renders a Scene, sends
// {type:"select"} on click and {type:"rotate", angle} on drag-release. The
// drag preview is cosmetic; the scene pushed back after a move is the REAL
// checker's result (overlaps tinted red).

const SVG_NS = "http://www.w3.org/2000/svg";
const ACCENT = "#14b8a6"; // teal: selection
const OVERLAP = "#ef4444"; // red: flagged overlap
const NESTED = "#9aa0a6";
const CROSSING = "#f59e0b";
const BACKBONE = "#c4c9d0";
const MAX_PX = 640;

function el(name, attrs) {
  const node = document.createElementNS(SVG_NS, name);
  for (const k in attrs) node.setAttribute(k, attrs[k]);
  return node;
}

export function render({ model, el: root }) {
  root.innerHTML = "";
  root.style.fontFamily = "system-ui, sans-serif";

  const wrap = document.createElement("div");
  wrap.style.cssText = "display:inline-block;border:1px solid #e0e0e0;border-radius:8px;padding:6px;background:#fff;";
  const svg = document.createElementNS(SVG_NS, "svg");
  const status = document.createElement("div");
  status.style.cssText = "margin-top:6px;font-size:13px;color:#374151;min-height:18px;";
  wrap.appendChild(svg);
  wrap.appendChild(status);
  root.appendChild(wrap);

  // group holding all flipped-y geometry; recreated each draw.
  let gGeom = null;
  let vp = { min_x: 0, min_y: 0, w: 1, h: 1 };

  // drag state
  let drag = null; // {pivot:[lx,ly], startAngle, elems:[...]}

  function localOf(nt) {
    return [nt.x - vp.min_x, nt.y - vp.min_y];
  }

  function pointerLocal(evt) {
    const pt = svg.createSVGPoint();
    pt.x = evt.clientX;
    pt.y = evt.clientY;
    const ctm = gGeom.getScreenCTM();
    if (!ctm) return null;
    const p = pt.matrixTransform(ctm.inverse());
    return [p.x, p.y];
  }

  function draw() {
    const scene = model.get("scene") || {};
    vp = scene.viewport || { min_x: 0, min_y: 0, w: 1, h: 1 };
    const nts = scene.nucleotides || [];
    const pairs = scene.pairs || [];
    const routed = scene.routed_lines || [];

    const scale = Math.min(MAX_PX / vp.w, MAX_PX / vp.h, 4);
    svg.setAttribute("width", (vp.w * scale).toFixed(1));
    svg.setAttribute("height", (vp.h * scale).toFixed(1));
    svg.setAttribute("viewBox", `0 0 ${vp.w} ${vp.h}`);
    svg.style.touchAction = "none";
    svg.innerHTML = "";

    // flip y: engine y-up -> screen y-down
    gGeom = el("g", { transform: `translate(0,${vp.h}) scale(1,-1)` });
    svg.appendChild(gGeom);

    // backbone polyline through nucleotides in index order
    if (nts.length > 1) {
      const pts = nts.map((n) => localOf(n).map((v) => v.toFixed(2)).join(",")).join(" ");
      gGeom.appendChild(el("polyline", { points: pts, fill: "none", stroke: BACKBONE, "stroke-width": 2 }));
    }

    // base-pair lines
    for (const p of pairs) {
      const a = nts[p.i], b = nts[p.j];
      if (!a || !b) continue;
      const [ax, ay] = localOf(a), [bx, by] = localOf(b);
      gGeom.appendChild(el("line", {
        x1: ax, y1: ay, x2: bx, y2: by,
        stroke: p.kind === "crossing" ? CROSSING : NESTED,
        "stroke-width": p.kind === "crossing" ? 2 : 3,
        "stroke-dasharray": p.kind === "crossing" ? "4,3" : "none",
      }));
    }

    // routed (PK-A) polylines
    for (const r of routed) {
      const pts = (r.points || []).map((pt) => [pt[0] - vp.min_x, pt[1] - vp.min_y].map((v) => v.toFixed(2)).join(",")).join(" ");
      gGeom.appendChild(el("polyline", { points: pts, fill: "none", stroke: CROSSING, "stroke-width": 2, "stroke-dasharray": "2,3" }));
    }

    // nucleotide disks
    for (const n of nts) {
      const [cx, cy] = localOf(n);
      const c = el("circle", {
        id: `nt-${n.id}`, cx, cy, r: n.r, fill: n.fill,
        stroke: "#5f6368", "stroke-width": 0.75, "data-id": n.id, cursor: "pointer",
      });
      gGeom.appendChild(c);
      c.addEventListener("pointerdown", (e) => onPointerDown(e, n.id));
    }

    updateHighlight();
  }

  function updateHighlight() {
    const scene = model.get("scene") || {};
    const sel = new Set(model.get("selection") || []);
    const overlaps = new Set(scene.overlaps || []);
    const nts = scene.nucleotides || [];
    for (const n of nts) {
      const c = svg.querySelector(`#nt-${n.id}`);
      if (!c) continue;
      if (scene.flagged && overlaps.has(n.id)) {
        c.setAttribute("fill", OVERLAP);
      } else {
        c.setAttribute("fill", n.fill);
      }
      if (sel.has(n.id)) {
        c.setAttribute("stroke", ACCENT);
        c.setAttribute("stroke-width", 2.5);
      } else {
        c.setAttribute("stroke", "#5f6368");
        c.setAttribute("stroke-width", 0.75);
      }
    }
  }

  function updateStatus() {
    status.textContent = model.get("status") || "";
  }

  // -- interaction --------------------------------------------------------

  function onPointerDown(e, id) {
    const sel = model.get("selection") || [];
    const scene = model.get("scene") || {};
    if (sel.includes(id) && scene.pivot) {
      // begin drag-rotate of the already-selected helix
      const pivotLocal = [scene.pivot[0] - vp.min_x, scene.pivot[1] - vp.min_y];
      const p = pointerLocal(e);
      if (!p) return;
      const elems = sel.map((k) => svg.querySelector(`#nt-${k}`)).filter(Boolean);
      drag = {
        pivotLocal,
        startAngle: Math.atan2(p[1] - pivotLocal[1], p[0] - pivotLocal[0]),
        delta: 0,
        elems,
      };
      e.target.setPointerCapture(e.pointerId);
      e.preventDefault();
    } else {
      // plain click -> select
      model.send({ type: "select", index: id });
    }
  }

  function onPointerMove(e) {
    if (!drag) return;
    const p = pointerLocal(e);
    if (!p) return;
    const ang = Math.atan2(p[1] - drag.pivotLocal[1], p[0] - drag.pivotLocal[0]);
    drag.delta = ang - drag.startAngle;
    const deg = (drag.delta * 180) / Math.PI;
    const [px, py] = drag.pivotLocal;
    for (const c of drag.elems) c.setAttribute("transform", `rotate(${deg}, ${px}, ${py})`);
    status.textContent = `Rotating ${deg >= 0 ? "+" : ""}${deg.toFixed(1)} deg (release to commit)`;
  }

  function onPointerUp() {
    if (!drag) return;
    const delta = drag.delta;
    for (const c of drag.elems) c.removeAttribute("transform");
    drag = null;
    // authoritative move: kernel rotates + re-checks + pushes back scene
    model.send({ type: "rotate", angle: delta });
  }

  svg.addEventListener("pointermove", onPointerMove);
  window.addEventListener("pointerup", onPointerUp);

  model.on("change:scene", draw);
  model.on("change:selection", updateHighlight);
  model.on("change:status", updateStatus);

  draw();
  updateStatus();

  return () => {
    window.removeEventListener("pointerup", onPointerUp);
  };
}

export default { render };
