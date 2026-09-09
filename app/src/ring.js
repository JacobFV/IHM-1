// The ring: what is exchanging signal with the thing you are looking at, and
// how much of it, drawn at the size it was measured.
//
// The IBM-1 site rings a brain with its materializations and runs leader lines
// to the structures each one touches. The interaction is reused here; the ring
// is not. On the site the ring is a menu of models. Here the signals are live
// and there is a body on the other end, so the ring is a ring of the things
// currently exchanging signal, and clicking one answers a question you cannot
// otherwise ask: what is this contributing, right now, to the thing I am
// looking at?
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
// The one genuinely load-bearing result -- severing the cortical dynamics drops
// the body at 2.90 s while the intact one holds at 0.17 mm -- is magnitude 1.0
// and comes out as the thick arrow it actually is. Everything else is thin
// because everything else is thin.

const SVG = "http://www.w3.org/2000/svg";

// A hairline floor so a dead edge is still hittable, and a ceiling that is the
// widest a measurement of 1.0 draws. Between them the map is LINEAR: a log
// scale would draw 0.0005 at a third of full width, which is the exact false
// impression this interface exists not to give. The floor is a rendering
// minimum and never a claim -- the number is always printed beside it.
export const HAIRLINE_PX = 0.6;
export const FULL_PX = 14;

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

/** Merge the per-side records of one nerve into one ring item.
 *
 * A ring with "left median nerve" and "right median nerve" as separate labels
 * is a ring nobody can read; the sides are kept on the item so the panel can
 * still name both, and the delays are identical by construction only when the
 * measured lengths are, which they are not, so both are carried.
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
    // label covers both.
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
 * @param host      the element the ring's buttons and panel live in
 * @param project   (point_m) => [x, y] in host pixels, or null when unprojectable
 * @param onSystem  called when the highlighted system changes
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
  chooser.setAttribute("aria-label", "Highlighted system");
  const note = el("p", "note ring-basis");
  note.id = "ring-basis";
  note.setAttribute("role", "status");
  const more = el("p", "note ring-more");
  more.id = "ring-more";
  more.hidden = true;
  // The chooser and the basis line are one stacked control: the basis is a
  // statement about the widths below it and has to sit under them, not float
  // over whichever row the buttons happen to wrap onto.
  const controls = el("div", "ring-controls");
  controls.id = "ring-controls";
  controls.append(chooser, note, more);
  host.append(lines, layer, controls);
  // The panel is a reading surface, not an overlay: it goes in the column of
  // panes with every other readout, which leaves the whole middle of the view
  // to the ring itself.
  (panelHost || host).append(panel);

  let graph = null, state = null, system = null, selected = null, items = [];
  let width = 0, height = 0, visible = true;

  function setVisible(on) {
    visible = !!on;
    layer.hidden = !visible;
    lines.style.display = visible ? "" : "none";
    note.hidden = !visible;
    drawChooser();
    if (visible) { drawRing(); } else { panel.hidden = true; }
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
    // The ring is the centrepiece, so it is on; it also covers the body, so it
    // folds away to its own title the way every pane here does.
    const fold = el("button", "ring-fold", visible ? "▾ Ring" : "▸ Ring");
    fold.type = "button";
    fold.id = "ring-fold";
    fold.setAttribute("aria-expanded", String(visible));
    fold.onclick = () => { setVisible(!visible); };
    chooser.append(fold);
    if (!visible) return;
    for (const entry of graph?.systems || []) {
      const button = el("button", "ring-system", entry.label);
      button.type = "button";
      button.dataset.system = entry.id;
      button.setAttribute("aria-pressed", String(entry.id === system));
      button.title = entry.note || "";
      button.onclick = () => highlight(entry.id);
      chooser.append(button);
    }
  }

  function drawRing() {
    layer.replaceChildren();
    lines.replaceChildren();
    items = [];
    if (!graph || !system) return;
    const all = contributors(graph, state, system);
    // A ring nobody can read is not a ring. The tail is kept, counted, and
    // reachable from the panel rather than silently dropped.
    const shown = all.slice(0, 26);
    const hidden = all.length - shown.length;
    shown.forEach((item, i) => {
      const button = el("button", "ring-item");
      button.type = "button";
      button.dataset.edge = item.key;
      button.dataset.kind = item.magnitude_kind;
      button.setAttribute("aria-pressed", String(selected === item.key));
      const swatch = svgEl("svg", { class: "ring-gauge", viewBox: "0 0 20 16", "aria-hidden": "true" });
      const stroke = arrowWidth(item.magnitude, item.magnitude_kind);
      const rule = svgEl("line", {
        x1: "1", y1: "8", x2: "19", y2: "8",
        "stroke-width": stroke.toFixed(2),
        "stroke-dasharray": item.magnitude_kind === "unmeasured" ? "2 2" : "",
      });
      swatch.append(rule);
      button.append(swatch, el("span", "ring-label", item.label),
        el("span", "ring-magnitude", magnitudeText(item)));
      button.onclick = () => choose(item.key);
      layer.append(button);
      const leader = svgEl("path", { class: "ring-leader", "stroke-width": stroke.toFixed(2) });
      if (item.magnitude_kind === "unmeasured") leader.setAttribute("stroke-dasharray", "3 4");
      lines.append(leader);
      items.push({ item, button, leader, index: i, count: shown.length });
    });
    more.hidden = true;
    if (hidden > 0) {
      more.textContent =
        `${hidden} further route${hidden === 1 ? "" : "s"} touch this system and are listed in the panel.`;
      more.hidden = false;
    }
    layout();
    drawPanel();
  }

  // The ellipse, and one leader per item from its inner edge to the system's
  // own anchor projected into the scene. Same construction as the site's ring:
  // the labels own the ellipse, the leaders own the middle.
  function layout() {
    width = host.clientWidth || 0;
    height = host.clientHeight || 0;
    if (!width || !height || !items.length) return;
    const style = getComputedStyle(host);
    const left = parseFloat(style.getPropertyValue("--left-inset")) || 12;
    const right = parseFloat(style.getPropertyValue("--right-inset")) || 12;
    // The band the ring may use: clear of the system chooser above and of the
    // transport and prompt bar below, both of which are measured rather than
    // assumed so a five-entry retrieval never pushes labels under the bar.
    // The band starts below the chooser and its basis line, measured, so an
    // extra row of system buttons never lands a ring label under the text that
    // explains what its width means.
    const hostTop = host.getBoundingClientRect?.().top ?? 0;
    const controlsBottom = controls.getBoundingClientRect?.().bottom ?? 0;
    const above = Math.max(120, controlsBottom - hostTop + 16);
    // The prompt panel overlays rather than reflows: letting a five-entry
    // retrieval squeeze the ring would pile 26 labels on top of each other.
    const below = Math.min((parseFloat(style.getPropertyValue("--prompt-inset")) || 48) + 116,
                           Math.max(160, height * 0.30));
    const cx = (left + (width - right)) / 2;
    const cy = (above + (height - below)) / 2;
    const rx = Math.max(140, (width - left - right) / 2 - 150);
    const ry = Math.max(90, (height - above - below) / 2);
    items.forEach(({ button, index, count }, i) => {
      // Start at the top and walk the ellipse; the sides carry the most items,
      // which is where the labels have room.
      const t = -Math.PI / 2 + (index / count) * Math.PI * 2;
      let x = cx + rx * Math.cos(t);
      const y = cy + ry * Math.sin(t);
      const side = Math.abs(Math.cos(t)) < 0.12 ? "c" : Math.cos(t) < 0 ? "l" : "r";
      // A left label hangs to the LEFT of its anchor, so the ellipse's own edge
      // is not where the text ends. Clamp by the widest label the CSS allows,
      // or the ring runs under the columns on both sides.
      const reach = 250;
      if (side === "l") x = Math.max(x, left + reach);
      if (side === "r") x = Math.min(x, width - right - reach);
      button.style.left = `${x}px`;
      button.style.top = `${y}px`;
      button.dataset.side = side;
      items[i].x = x;
      items[i].y = y;
      items[i].side = side;
    });
    drawLeaders();
  }

  function anchor() {
    const entry = (graph?.systems || []).find((s) => s.id === system);
    const point = entry?.anchor_m;
    const projected = point && project ? project(point) : null;
    if (projected && Number.isFinite(projected[0]) && Number.isFinite(projected[1])) return projected;
    const left = parseFloat(getComputedStyle(host).getPropertyValue("--left-inset")) || 12;
    const right = parseFloat(getComputedStyle(host).getPropertyValue("--right-inset")) || 12;
    return [(left + (width - right)) / 2, height / 2];
  }

  function drawLeaders() {
    if (!width || !height) return;
    lines.setAttribute("viewBox", `0 0 ${width} ${height}`);
    const [ax, ay] = anchor();
    for (const entry of items) {
      const { button, leader, item } = entry;
      const dim = selected && selected !== item.key;
      leader.style.opacity = dim ? "0.12" : "";
      button.classList.toggle("is-dim", !!dim);
      button.classList.toggle("is-active", selected === item.key);
      const box = button.getBoundingClientRect?.();
      const hostBox = host.getBoundingClientRect?.();
      let x0 = entry.x, y0 = entry.y;
      if (box && hostBox && box.width) {
        y0 = box.top - hostBox.top + box.height / 2;
        x0 = entry.side === "l" ? box.right - hostBox.left + 4 : box.left - hostBox.left - 4;
        if (entry.side === "c") x0 = box.left - hostBox.left + box.width / 2;
      }
      const mx = (x0 + ax) / 2, my = (y0 + ay) / 2;
      leader.setAttribute("d", `M${x0.toFixed(1)},${y0.toFixed(1)} Q${mx.toFixed(1)},${y0.toFixed(1)} ${ax.toFixed(1)},${ay.toFixed(1)}`);
    }
  }

  function drawPanel() {
    if (!graph || !system) { panel.hidden = true; return; }
    panel.hidden = false;
    panel.replaceChildren();
    const entry = (graph.systems || []).find((s) => s.id === system);
    const all = contributors(graph, state, system);
    const item = selected ? all.find((x) => x.key === selected) : null;

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
      // Selecting the system rather than a contributor: the converse view.
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
    drawLeaders();
    drawPanel();
  }

  function highlight(id) {
    if (!graph) return;
    system = id;
    selected = null;
    for (const button of chooser.querySelectorAll("button"))
      button.setAttribute("aria-pressed", String(button.dataset.system === id));
    drawRing();
    onSystem?.(id);
  }

  return {
    element: layer,
    panel,
    setGraph(next) {
      graph = next;
      system = system || next?.systems?.[0]?.id || null;
      drawChooser();
      setBasis();
      drawRing();
    },
    setState(next) { state = next; setBasis(); drawRing(); },
    setVisible,
    highlight,
    resize: layout,
    follow: drawLeaders,
    get system() { return system; },
    get selected() { return selected; },
  };
}
