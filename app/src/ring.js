// The annotations: what is exchanging signal with this body, each label standing
// at the projected position of the 3D point it is about.
//
// This used to be a ring. Labels were dealt onto a screen-space ellipse in the
// order they happened to rank, and every leader ran to ONE anchor -- the
// highlighted system's -- so a label's position said nothing about where its
// structure was, and only one system's 26 strongest edges could be shown at all.
//
// Now every annotation carries its own point in canonical metres, the same point
// the geometry is drawn from. Each frame that point is projected and the label is
// placed where it lands, so the labels orbit with the body because they ARE on
// the body. That is what makes the count affordable: the boxes are measured once
// when they are built and moved by transform afterwards, never re-measured, so
// 150 annotations cost 150 projections and no layout at all.
//
// **The constraint that decides every other decision here: arrow width IS the
// measured magnitude.** Most of these paths are near zero. Only 0.03-0.07% of a
// sensory drive reaches the disjoint motor region, and severing the cortical
// kernel changes the motor command by nothing. An arrow glowing from sensory to
// motor cortex would assert something this project has measured to be false, so
// there is no flow animation anywhere in this file, no log scale that would
// flatter a dead path into legibility, and no default width for an edge whose
// magnitude nobody measured. What there is instead:
//
//   measured    a solid arrow, width linear in the measurement
//   live        the same, read from the running body this tick
//   declared    a gain or latency, drawn as a bracket and never as width
//   unmeasured  a dashed hairline that says so
//
// Two honesty rules the 3D placement adds:
//
//   * An annotation whose edge carries no point of its own falls back to its
//     system's anchor and says so (`data-anchor="system"`). Its label is near the
//     right structure and not ON its own, and that distinction is visible rather
//     than implied.
//   * A point behind the camera or off the viewport is not drawn. Pinning such a
//     label to the nearest edge of the screen would put a label for the back of
//     the body over the front of it.

const SVG = "http://www.w3.org/2000/svg";

// A hairline floor so a dead edge is still hittable, and a ceiling that is the
// widest a measurement of 1.0 draws. Between them the map is LINEAR: a log
// scale would draw 0.0005 at a third of full width, which is the exact false
// impression this interface exists not to give. The floor is a rendering
// minimum and never a claim -- the number is always printed beside it.
export const HAIRLINE_PX = 0.6;
export const FULL_PX = 14;

// How far a label sits from its own point. Crowded labels take another direction
// around the same point rather than a longer leader: past about four times this
// the line is long enough that which point it came from stops being obvious, and
// a label that cannot be placed within it is hidden and counted instead.
const LEADER_PX = 34;

export function arrowWidth(magnitude, kind) {
  if (kind === "declared" || kind === "unmeasured" || magnitude == null) return HAIRLINE_PX;
  const value = Math.max(0, Math.min(1, Number(magnitude)));
  if (!Number.isFinite(value)) return HAIRLINE_PX;
  return HAIRLINE_PX + value * (FULL_PX - HAIRLINE_PX);
}

/** How a magnitude reads in the panel, in the units it was measured in. */
export function magnitudeText(edge) {
  if (edge.magnitude_kind === "unmeasured") return "not measured";
  if (edge.magnitude_kind === "declared")
    return Number.isFinite(edge.declared_gain)
      ? `gain ${edge.declared_gain >= 0 ? "+" : ""}${edge.declared_gain.toFixed(2)} · declared`
      : "declared parameter";
  if (edge.magnitude == null) return "not measured";
  const percent = edge.magnitude * 100;
  const value = percent >= 1 ? `${percent.toFixed(1)}%`
    : percent >= 0.001 ? `${percent.toFixed(3)}%`
    : `${percent.toExponential(1)}%`;
  return edge.magnitude_kind === "live" ? `${value} · live this tick` : `${value} · measured`;
}

/** Merge the per-side records of one nerve into one annotation.
 *
 * Two labels reading "left median nerve" and "right median nerve" are two labels
 * nobody needs; the sides are kept on the item so the panel can still name both,
 * and the delays are identical by construction only when the measured lengths
 * are, which they are not, so both are carried.
 */
function mergeKey(edge) {
  if (edge.kind !== "route") return edge.id;
  return edge.id.replace(/^peripheral-nerve-(left|right)-/, "route:");
}

/** Everything touching `systemId`, ranked by what it is currently carrying.
 *
 * `state.edges` overrides the graph's static magnitude where the running body
 * has one. With no body running nothing is overridden and the measured values
 * stand -- labelled measured, not live.
 */
export function contributors(graph, state, systemId) {
  const overrides = state?.edges || {};
  const items = new Map();
  for (const raw of graph?.edges || []) {
    const touches = raw.system === systemId || raw.target === systemId;
    if (!touches) continue;
    const override = overrides[raw.id] || {};
    const edge = { ...raw, ...override };
    edge.direction = raw.target === systemId ? "into" : "out of";
    const key = mergeKey(raw);
    const held = items.get(key);
    if (!held) { items.set(key, { ...edge, key, members: [edge] }); continue; }
    held.members.push(edge);
    // Two sides of one nerve: keep the larger measured magnitude and say the
    // label covers both. The anchor stays the first side's own point, so the
    // label stands on one of the two structures it names rather than between
    // them, which is a place neither of them is.
    if ((edge.magnitude ?? -1) > (held.magnitude ?? -1)) {
      held.magnitude = edge.magnitude;
      held.magnitude_kind = edge.magnitude_kind;
      held.measurement = edge.measurement;
    }
    held.label = held.label.replace(/^(left|right)\s+/i, "");
    held.side = "both";
  }
  const order = { live: 0, measured: 1, declared: 2, unmeasured: 3 };
  return [...items.values()].sort((a, b) => {
    const byMagnitude = (b.magnitude ?? -1) - (a.magnitude ?? -1);
    if (byMagnitude) return byMagnitude;
    const byKind = (order[a.magnitude_kind] ?? 9) - (order[b.magnitude_kind] ?? 9);
    return byKind || a.label.localeCompare(b.label);
  });
}

/** Delay per fibre class, longest last: one route, N delays, which is the join. */
export function fibreRows(item) {
  const rows = [];
  for (const member of item.members || [item])
    for (const fibre of member.fibres || [])
      rows.push({ side: member.side, ...fibre });
  return rows.sort((a, b) => a.delay_s - b.delay_s);
}

/** Screen-space placement for already-projected annotations.
 *
 * Exported because it is the part worth testing without a browser: given points
 * and box sizes it returns where each label goes, and it is a pure function of
 * its inputs. Labels are placed at their own point, offset to whichever side has
 * room, and nudged vertically only when they would overlap one already placed.
 * Ordered by magnitude, so the largest measurement keeps the position it earned
 * and the hairlines move around it.
 */
export function placeLabels(entries, width, height, bounds) {
  const area = {
    left: bounds?.left ?? 0, right: bounds?.right ?? width,
    top: bounds?.top ?? 0, bottom: bounds?.bottom ?? height,
  };
  const placed = [];
  const out = [];
  // Candidate directions, tried in order: straight out from the middle of the
  // view first, then progressively steeper, then back across the body as a last
  // resort. Going OUTWARD matters -- a label placed toward the centre covers the
  // structure it is naming, and with a hundred of them the body disappears under
  // its own annotations.
  const ANGLES = [0, 22, -22, 45, -45, 68, -68, 90, -90, 112, -112, 135, -135, 158, -158, 180];
  const RADII = [LEADER_PX, LEADER_PX * 1.8, LEADER_PX * 2.8, LEADER_PX * 4];
  const overlaps = (cx, cy, w, h) => placed.some((p) =>
    Math.abs(p.cx - cx) * 2 < p.w + w + 6 && Math.abs(p.cy - cy) * 2 < p.h + h + 4);
  for (const entry of entries) {
    const { x, y, w, h } = entry;
    const outward = x < width / 2 ? -1 : 1;
    let spot = null;
    for (const radius of RADII) {
      for (const angle of ANGLES) {
        const t = (angle * Math.PI) / 180;
        const dx = outward * Math.cos(t) * radius;
        const dy = Math.sin(t) * radius;
        // The box hangs off the anchor on the side the direction points, so its
        // near edge is what the leader meets and the far edge is what grows.
        const cx = x + dx + (dx >= 0 ? w / 2 : -w / 2);
        const cy = y + dy;
        // Reject rather than clamp: a clamped position collapses every label
        // pushed past an edge onto one pixel, where the overlap test no longer
        // sees them and they stack silently. The bounds are the clear middle of
        // the view, not the viewport -- the columns of panes overlay it, and a
        // label under one of them is a label nobody can read.
        if (cx - w / 2 < area.left || cx + w / 2 > area.right) continue;
        if (cy - h / 2 < area.top || cy + h / 2 > area.bottom) continue;
        if (overlaps(cx, cy, w, h)) continue;
        spot = { cx, cy, side: dx >= 0 ? "r" : "l" };
        break;
      }
      if (spot) break;
    }
    if (!spot) { out.push({ ...entry, hidden: true }); continue; }
    placed.push({ cx: spot.cx, cy: spot.cy, w, h });
    out.push({ ...entry, ...spot, hidden: false });
  }
  return out;
}

const el = (tag, cls, text) => {
  const node = document.createElement(tag);
  if (cls) node.className = cls;
  if (text !== undefined) node.textContent = text;
  return node;
};
const svgEl = (tag, attributes = {}) => {
  const node = document.createElementNS(SVG, tag);
  for (const key in attributes) node.setAttribute(key, attributes[key]);
  return node;
};
const escaped = (s) => String(s ?? "");

/**
 * @param host      the element the annotations and their controls live in
 * @param project   (point_m) => [x, y] or [x, y, ndcZ] in host pixels, or null.
 *                  An ndcZ above 1 is behind the camera and is not drawn.
 * @param onSystem  called when the focused system changes
 */
export function mountRing(host, { project, onSystem, panelHost } = {}) {
  const layer = el("div", "ring-layer");
  layer.id = "ring";
  const lines = svgEl("svg", { class: "ring-leaders", "aria-hidden": "true" });
  const panel = el("aside", "ring-panel");
  panel.id = "ring-panel";
  panel.hidden = true;
  const chooser = el("div", "ring-systems");
  chooser.id = "ring-systems";
  chooser.setAttribute("role", "group");
  chooser.setAttribute("aria-label", "Annotation groups");
  const note = el("p", "note ring-basis");
  note.id = "ring-basis";
  note.setAttribute("role", "status");
  const more = el("p", "note ring-more");
  more.id = "ring-more";
  more.hidden = true;
  const controls = el("div", "ring-controls");
  controls.id = "ring-controls";
  controls.append(chooser, note, more);
  host.append(lines, layer, controls);
  // The panel is a reading surface, not an overlay: it goes in the column of
  // panes with every other readout.
  (panelHost || host).append(panel);

  let graph = null, state = null, system = null, selected = null, items = [];
  let width = 0, height = 0, visible = true;
  const groups = new Set();

  function setVisible(on) {
    visible = !!on;
    layer.hidden = !visible;
    lines.style.display = visible ? "" : "none";
    note.hidden = !visible;
    drawChooser();
    if (visible) { build(); } else { panel.hidden = true; }
  }

  function setBasis() {
    if (!graph) { note.textContent = "Signal graph unavailable."; return; }
    const live = state?.live;
    note.textContent = live
      ? `Live · t ${(state.time_s ?? 0).toFixed(2)} s${state.severed ? " · cortical kernel severed" : ""}. Arrow width is the magnitude this tick.`
      : "Measured, not live. No body is running, so every width below is a stored measurement with its receipt — nothing is animated to stand in for a signal that is not flowing.";
    note.dataset.live = String(!!live);
  }

  function drawChooser() {
    chooser.replaceChildren();
    const fold = el("button", "ring-fold", visible ? "▾ Annotations" : "▸ Annotations");
    fold.type = "button";
    fold.id = "ring-fold";
    fold.setAttribute("aria-expanded", String(visible));
    fold.onclick = () => { setVisible(!visible); };
    chooser.append(fold);
    if (!visible) return;
    // Groups are a multi-select: any number of them can be on at once, which is
    // the whole reason the labels had to leave the ellipse.
    for (const entry of graph?.systems || []) {
      const on = groups.has(entry.id);
      const count = countFor(entry.id);
      const button = el("button", "ring-system", `${entry.label} · ${count}`);
      button.type = "button";
      button.dataset.system = entry.id;
      button.dataset.group = entry.id;
      button.setAttribute("aria-pressed", String(on));
      button.title = entry.note || "";
      // A plain click is the toggle: off turns the group on and focuses it, on
      // turns it off. Clicking a group that is already shown but not focused
      // focuses it instead of hiding it, so reading the panel for one group
      // never costs you the annotations of another. No modifier key decides
      // anything here -- a control nobody can find is a control nobody has.
      button.onclick = () => {
        if (!groups.has(entry.id)) { groups.add(entry.id); system = entry.id; selected = null; onSystem?.(entry.id); }
        else if (system !== entry.id) { system = entry.id; selected = null; onSystem?.(entry.id); }
        else {
          groups.delete(entry.id);
          system = [...groups][0] || null;
        }
        drawChooser();
        build();
      };
      chooser.append(button);
    }
    const all = el("button", "ring-system ring-all", groups.size === (graph?.systems || []).length ? "none" : "all");
    all.type = "button";
    all.id = "ring-all";
    all.onclick = () => {
      const every = (graph?.systems || []).map((s) => s.id);
      if (groups.size === every.length) groups.clear();
      else for (const id of every) groups.add(id);
      drawChooser();
      build();
    };
    chooser.append(all);
  }

  function countFor(id) {
    return contributors(graph, state, id).length;
  }

  function toggleGroup(id) {
    if (groups.has(id)) groups.delete(id); else groups.add(id);
    if (!groups.has(system)) system = [...groups][0] || null;
    drawChooser();
    build();
  }

  function focusGroup(id) {
    groups.clear();
    groups.add(id);
    system = id;
    selected = null;
    drawChooser();
    build();
    onSystem?.(id);
  }

  /** One annotation per contributor of every visible group, built once. */
  function build() {
    layer.replaceChildren();
    lines.replaceChildren();
    items = [];
    if (!graph || !groups.size) { panel.hidden = true; return; }
    const systemAnchor = new Map((graph.systems || []).map((s) => [s.id, s.anchor_m]));
    const seen = new Set();
    for (const id of groups) {
      for (const item of contributors(graph, state, id)) {
        // One structure can touch two visible groups; it gets one label.
        if (seen.has(item.key)) continue;
        seen.add(item.key);
        const own = item.anchor_m || (item.members || []).map((m) => m.anchor_m).find(Boolean);
        const point = own || systemAnchor.get(id);
        if (!point) continue;
        const button = el("button", "ring-item");
        button.type = "button";
        button.dataset.edge = item.key;
        button.dataset.kind = item.magnitude_kind;
        button.dataset.group = id;
        button.dataset.anchor = own ? "own" : "system";
        if (!own) button.title = "No point of its own: placed at this system's anchor.";
        button.setAttribute("aria-pressed", String(selected === item.key));
        const swatch = svgEl("svg", { class: "ring-gauge", viewBox: "0 0 20 16", "aria-hidden": "true" });
        const stroke = arrowWidth(item.magnitude, item.magnitude_kind);
        swatch.append(svgEl("line", {
          x1: "1", y1: "8", x2: "19", y2: "8",
          "stroke-width": stroke.toFixed(2),
          "stroke-dasharray": item.magnitude_kind === "unmeasured" ? "2 2" : "",
        }));
        button.append(swatch, el("span", "ring-label", item.label),
          el("span", "ring-magnitude", magnitudeText(item)));
        button.onclick = () => choose(item.key);
        layer.append(button);
        const leader = svgEl("path", { class: "ring-leader", "stroke-width": stroke.toFixed(2) });
        if (item.magnitude_kind === "unmeasured") leader.setAttribute("stroke-dasharray", "3 4");
        lines.append(leader);
        // Measured once, here. Everything after this is a transform.
        items.push({ item, button, leader, point, group: id,
                     w: button.offsetWidth || 180, h: button.offsetHeight || 20 });
      }
    }
    follow();
    drawPanel();
  }

  /** Project every anchor, place the labels, run the leaders. Per frame. */
  function follow() {
    width = host.clientWidth || 0;
    height = host.clientHeight || 0;
    if (!width || !height) return;
    lines.setAttribute("viewBox", `0 0 ${width} ${height}`);
    const candidates = [];
    for (const entry of items) {
      const p = project ? project(entry.point) : null;
      const behind = p && p.length > 2 && !(p[2] <= 1);
      const onScreen = p && Number.isFinite(p[0]) && Number.isFinite(p[1])
        && p[0] >= -80 && p[0] <= width + 80 && p[1] >= -80 && p[1] <= height + 80;
      if (!p || behind || !onScreen) { hide(entry); continue; }
      if (!entry.w) { entry.w = entry.button.offsetWidth || 180; entry.h = entry.button.offsetHeight || 20; }
      candidates.push({ entry, x: p[0], y: p[1], w: entry.w, h: entry.h });
    }
    // The band the labels may use: clear of the columns of panes on either side,
    // of the annotation controls above, and of the transport and prompt bar
    // below. All measured rather than assumed, so a wrapped row of group
    // buttons never lands a label under the text that explains what it means.
    const style = getComputedStyle(host);
    const hostTop = host.getBoundingClientRect?.().top ?? 0;
    const bounds = {
      left: (parseFloat(style.getPropertyValue("--left-inset")) || 12) + 8,
      right: width - (parseFloat(style.getPropertyValue("--right-inset")) || 12) - 8,
      top: Math.max(8, (controls.getBoundingClientRect?.().bottom ?? 0) - hostTop + 10),
      bottom: height - Math.min((parseFloat(style.getPropertyValue("--prompt-inset")) || 48) + 70,
                                Math.max(90, height * 0.22)),
    };
    let dropped = 0;
    for (const placedEntry of placeLabels(candidates, width, height, bounds)) {
      const { entry, x, y, cx, cy, side, hidden } = placedEntry;
      if (hidden) { hide(entry); dropped += 1; continue; }
      entry.button.hidden = false;
      entry.button.dataset.side = side;
      entry.button.style.transform = `translate3d(${(cx - entry.w / 2).toFixed(1)}px, ${(cy - entry.h / 2).toFixed(1)}px, 0)`;
      const edgeX = side === "r" ? cx - entry.w / 2 - 3 : cx + entry.w / 2 + 3;
      entry.leader.style.display = "";
      entry.leader.setAttribute("d",
        `M${edgeX.toFixed(1)},${cy.toFixed(1)} Q${((edgeX + x) / 2).toFixed(1)},${cy.toFixed(1)} ${x.toFixed(1)},${y.toFixed(1)}`);
    }
    const off = items.length - candidates.length;
    more.hidden = !(off || dropped);
    if (!more.hidden)
      more.textContent = [
        off ? `${off} out of view` : "",
        dropped ? `${dropped} too crowded to place` : "",
      ].filter(Boolean).join(" · ") + `, of ${items.length} annotations.`;
    paintSelection();
  }

  function hide(entry) {
    entry.button.hidden = true;
    entry.leader.style.display = "none";
  }

  // Dimming is a selection change, not a camera change, so it is its own pass
  // and does not ride along on the per-frame one.
  function paintSelection() {
    for (const entry of items) {
      const dim = selected && selected !== entry.item.key;
      entry.leader.style.opacity = dim ? "0.12" : "";
      entry.button.classList.toggle("is-dim", !!dim);
      entry.button.classList.toggle("is-active", selected === entry.item.key);
    }
  }

  function drawPanel() {
    if (!graph || !system) { panel.hidden = true; return; }
    panel.hidden = false;
    panel.replaceChildren();
    const entry = (graph.systems || []).find((s) => s.id === system);
    const all = contributors(graph, state, system);
    const item = selected
      ? all.find((x) => x.key === selected)
        || [...groups].flatMap((g) => contributors(graph, state, g)).find((x) => x.key === selected)
      : null;

    const head = el("div", "ring-head");
    head.append(el("div", "eyebrow", item ? "Contribution" : "Everything feeding this"));
    head.append(el("h3", null, item ? item.label : entry?.label || system));
    if (item) {
      const back = el("button", "ring-back", "← all contributors");
      back.type = "button";
      back.onclick = () => choose(null);
      head.append(back);
    }
    panel.append(head);

    if (!item) {
      panel.append(el("p", "note", entry?.note || ""));
      const list = el("ol", "ring-rank");
      for (const row of all) {
        const line = el("li");
        line.dataset.kind = row.magnitude_kind;
        const gauge = svgEl("svg", { class: "ring-gauge wide", viewBox: "0 0 120 16", "aria-hidden": "true" });
        const stroke = arrowWidth(row.magnitude, row.magnitude_kind);
        gauge.setAttribute("preserveAspectRatio", "none");
        gauge.append(svgEl("line", {
          x1: "1", y1: "8", x2: "119", y2: "8", "stroke-width": stroke.toFixed(2),
          "stroke-dasharray": row.magnitude_kind === "unmeasured" ? "2 2" : "",
        }));
        const name = el("button", "ring-rank-name", row.label);
        name.type = "button";
        name.onclick = () => choose(row.key);
        line.append(name, el("span", "ring-magnitude", magnitudeText(row)), gauge);
        list.append(line);
      }
      panel.append(list);
      panel.append(el("p", "note",
        "Ranked by what each is carrying. Width is the measurement, linearly — a path measured near zero is drawn near zero."));
      return;
    }

    const gauge = svgEl("svg", { class: "ring-gauge hero", viewBox: "0 0 260 20", "aria-hidden": "true" });
    const stroke = arrowWidth(item.magnitude, item.magnitude_kind);
    gauge.append(svgEl("line", {
      x1: "2", y1: "10", x2: "246", y2: "10", "stroke-width": stroke.toFixed(2),
      "stroke-dasharray": item.magnitude_kind === "unmeasured" ? "3 3" : "",
    }));
    gauge.append(svgEl("path", { d: "M246,3 L258,10 L246,17 Z", class: "ring-arrowhead" }));
    panel.append(gauge);

    const figures = el("dl", "ring-figures");
    const figure = (term, value) => { figures.append(el("dt", null, term), el("dd", null, value)); };
    figure("Magnitude", magnitudeText(item));
    figure("Direction", `${item.direction} ${entry?.label || system}`);
    if (item.muscle_ports) figure("Muscle ports written", String(item.muscle_ports));
    if (item.receptor_ports) figure("Receptor ports read", String(item.receptor_ports));
    if (Number.isFinite(item.path_length_m))
      figure("Measured route", `${(item.path_length_m * 1000).toFixed(0)} mm · ${escaped(item.path_length_scope || "")}`);
    if (Number.isFinite(item.loop_delay_s))
      figure("Loop delay", `${(item.loop_delay_s * 1000).toFixed(0)} ms`);
    if (item.evidence_kind) figure("Evidence", escaped(item.evidence_kind));
    if (item.runtime_support) figure("Runtime", escaped(item.runtime_support));
    panel.append(figures);

    const rows = fibreRows(item);
    if (rows.length) {
      panel.append(el("h4", null, "Fibre classes and the delay on each"));
      const table = el("table", "ring-fibres");
      const body = el("tbody");
      for (const row of rows) {
        const line = el("tr");
        line.append(el("td", null, row.fibre_class), el("td", null, `${row.velocity_m_s} m/s`),
          el("td", null, `${(row.delay_s * 1000).toFixed(1)} ms`),
          el("td", null, row.side || ""));
        body.append(line);
      }
      table.append(body);
      panel.append(table);
      panel.append(el("p", "note",
        `One route, ${rows.length} delays. Composition is IBM-1's declaration; the length is this body's measured route.`));
    }

    if (item.measurement) panel.append(el("p", "ring-measurement", item.measurement));
    if (item.source) panel.append(el("p", "note", `Source: ${item.source}`));
    for (const limitation of item.limitations || []) panel.append(el("p", "note", limitation));
    if (item.magnitude_kind === "unmeasured")
      panel.append(el("p", "ring-warn",
        "No measurement exists for this edge. That is not the same as a measured zero, and this line makes no claim about what it is carrying."));
  }

  function choose(key) {
    selected = key;
    for (const entry of items)
      entry.button.setAttribute("aria-pressed", String(selected === entry.item.key));
    paintSelection();
    drawPanel();
  }

  function highlight(id) {
    if (!graph) return;
    focusGroup(id);
  }

  return {
    element: layer,
    panel,
    setGraph(next) {
      graph = next;
      system = system || next?.systems?.[0]?.id || null;
      if (!groups.size && system) groups.add(system);
      drawChooser();
      setBasis();
      build();
    },
    setState(next) { state = next; setBasis(); build(); },
    setVisible,
    highlight,
    toggleGroup,
    resize: follow,
    follow,
    get system() { return system; },
    get groups() { return [...groups]; },
    get selected() { return selected; },
    get count() { return items.length; },
  };
}
