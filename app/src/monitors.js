// Every card owns one existing, live widget. Removing a card parks its DOM;
// simulation subscriptions and user inputs survive without duplicated controls.
export const MONITORS = [
  {id:'live',title:'Live body',description:'Current unified native physiology and neural/mechanical clock.',tags:['live','physiology','signals']},
  ...[1,2,3].map(i=>({id:`live-signal-${i}`,title:`Live signal ${i}`,description:'An independently selected signal from the active body, using actual recorded live samples.',tags:['live','physiology','signals']})),
  {id:'microvascular',title:'Local microvessels',description:'Source-conditioned vastus lateralis vessel graph with physical diameters, pressures and flows at local zoom.',tags:['anatomy','vascular','evidence','inspect']},
  {id:'temporal-spectrum',title:'Live Laplace spectrum',description:'Finite-window Laplace magnitude from actual live signal samples, with elapsed seconds, native units and selectable damping.',tags:['live','physiology','signals','spectra','analysis']},
  {id:'intake',title:'Food & drink',description:'Schedule native nutrient and water intake on the active body clock; inspect queued, issued and accepted events.',tags:['live','controls','physiology','food','drink']},
  {id:'motor',title:'Motor & skin inputs',description:'Named muscle descending drive, selective sensory/motor blocks and whole-Skin pressure.',tags:['live','controls','neural','muscle','skin']},
  {
    id: "selection",
    title: "Selection",
    description: "Anatomy, source and uncertainty for the selected structure.",
    tags: ["anatomy", "evidence", "inspect"],
  },
  {
    id: "signals",
    title: "Signal graph",
    description:
      "Recorded trajectories, Fourier power and finite Laplace spectra.",
    tags: ["signals", "physiology", "spectra", "analysis"],
  },
  {
    id: "experiments",
    title: "Experiments",
    description: "Choose computed body studies or run a native protocol.",
    tags: ["physiology", "protocol", "mechanics", "controls"],
  },
  {
    id: "playback",
    title: "Playback",
    description: "Inspect actual recorded spatial frames and playback time.",
    tags: ["time", "controls", "signals", "mechanics"],
  },
  {
    id: "systemic",
    title: "Systemic state",
    description:
      "Shared state and executing pathways of the selected whole-body experiment.",
    tags: ["physiology", "signals", "evidence"],
  },
  {
    id: "contact",
    title: "Local contact",
    description:
      "Computed shorts-panel forces, tissue state and material evidence.",
    tags: ["mechanics", "clothing", "signals", "evidence"],
  },
  {
    id: "evidence",
    title: "Evidence",
    description:
      "Measurement fits, system coverage, exchange and vascular audits.",
    tags: ["evidence", "calibration", "vascular", "anatomy"],
  },
  {
    id: "scene",
    title: "Body interaction",
    description:
      "Applied forces, current body time and environment limits.",
    tags: ["mechanics", "environment", "controls"],
  },
];
export function filterMonitors(query = "", tag = "all") {
  const words = query.toLowerCase().trim().split(/\s+/).filter(Boolean);
  return MONITORS.filter(
    (m) =>
      (tag === "all" || m.tags.includes(tag)) &&
      words.every((word) =>
        `${m.title} ${m.description} ${m.tags.join(" ")}`
          .toLowerCase()
          .includes(word),
      ),
  );
}
const $ = (id) => document.getElementById(id);
const node = (tag, className, text) => {
  const n = document.createElement(tag);
  if (className) n.className = className;
  if (text) n.textContent = text;
  return n;
};
function section(title, nodes, { open = false } = {}) {
  const fold = node("details", "view-section");
  fold.open = open;
  fold.append(node("summary", "", title));
  const body = node("div", "view-section-body");
  body.append(...nodes.filter(Boolean));
  fold.append(body);
  return fold;
}
export function mountMonitors() {
  const library = $("library-panel"),
    right = $("inspector-panel"),
    viewport = $("viewport");
  const toolbar = document.querySelector(".workspace-tools");
  viewport.append(toolbar);
  $("toggle-library").innerHTML = "◧ <span>View</span>";
  $("toggle-inspector").innerHTML = "◨ <span>Monitors</span>";
  $("toggle-signals").hidden = true;
  document.querySelector("#app > header").remove();
  library.setAttribute("aria-label", "View controls");
  const camera = document.querySelector(".view-tools"),
    display = document.querySelector(".display-controls");
  const renderSettings = node("div", "render-settings");
  for (const label of [...display.querySelectorAll("label")])
    renderSettings.append(label);
  renderSettings.append($("flow-field"));
  const playback = node("div", "playback-settings");
  playback.append($("play"), $("time"), $("time-value"));
  display.remove();
  const experiments = node("div", "experiment-widgets");
  experiments.append(
    $("regional-controls"),
    $("systemic-controls"),
    document.querySelector(".scenario"),
  );
  const heading = node("div", "sidebar-heading");
  heading.append(
    node("h1", "", "View"),
    node("span", "sidebar-identity", "IHM–1"),
    $("total"),
  );
  const anatomyTitle = document.querySelector(".results-heading"),
    layersTitle = document.querySelector(".layer-title");
  const anatomy = section(
    "Anatomy",
    [
      $("anatomy-view"),
      $("layer-note"),
      document.querySelector(".search"),
      layersTitle,
      $("systems"),
      anatomyTitle,
      $("structures"),
    ],
    { open: true },
  );
  const source = section("Sources & body status", [
    $("canonical-body"),
    $("body-status"),
    $("model-note"),
    $("source-inspection"),
  ]);
  const sceneControls = node("div");
  sceneControls.id = "scene-controls";
  const clothing = $("clothing-components");
  library.replaceChildren(
    heading,
    section("Camera & rendering", [camera, renderSettings], { open: true }),
    sceneControls,
    section("Clothing", [clothing]),
    anatomy,
    source,
  );
  const signals = $("signals-panel");
  signals.querySelector(".signal-header h2").textContent =
    "Recorded physiology";
  const sceneMonitor = node("div");
  sceneMonitor.id = "scene-monitor";
  sceneMonitor.append(
    node(
      "p",
      "muted",
      "Select a mechanics interaction to inspect its computed response.",
    ),
  );
  const systemic = node("div");
  systemic.append($("systemic-monitor"));
  systemic.append(
    node(
      "p",
      "inactive-monitor-note",
      "Select a whole-body experiment in Experiments to inspect its shared state.",
    ),
  );
  const contact = node("div");
  contact.append($("garment-contact-monitor"));
  contact.append(
    node(
      "p",
      "inactive-monitor-note",
      "Select the shorts-panel study in Experiments to inspect computed local mechanics.",
    ),
  );
  const liveContents={};
  for(const id of ['live','motor','intake','temporal-spectrum','microvascular','live-signal-1','live-signal-2','live-signal-3']){
    const mount=node('div');mount.id=id==='live'?'live-body-monitor':id==='motor'?'live-motor-monitor':`${id}-monitor`;
    mount.append(node('p','muted','Start Body to inspect its current computed state.'));
    liveContents[id]=mount;
  }
  const contents = {
    ...liveContents,
    selection: $("details"),
    signals,
    experiments,
    playback,
    systemic,
    contact,
    evidence: document.querySelector(".evidence"),
    scene: sceneMonitor,
  };
  right.setAttribute("aria-label", "Monitor cards");
  const rightHeading = node("div", "monitor-sidebar-heading");
  rightHeading.append(node("h2", "", "Monitors"));
  const add = node("button", "new-pane", "+ New pane");
  add.setAttribute("aria-label", "New pane");
  rightHeading.append(add);
  const stack = node("div", "monitor-stack");
  stack.id = "monitor-stack";
  const storage = node("div");
  storage.id = "monitor-storage";
  storage.hidden = true;
  right.replaceChildren(rightHeading, stack, storage);
  let saved;
  try {
    saved = JSON.parse(localStorage.getItem("ihm.monitors.v1"));
  } catch {}
  const defaults = ["live", "live-signal-1", "selection", "motor"];
  let active = Array.isArray(saved?.active)
    ? [...new Set(saved.active)].filter((id) => contents[id])
    : defaults;
  const collapsed = new Set(
      Array.isArray(saved?.collapsed) ? saved.collapsed : [],
    ),
    compact = new Set(Array.isArray(saved?.compact) ? saved.compact : []);
  const cards = new Map();
  const save = () => {
    try {
      localStorage.setItem(
        "ihm.monitors.v1",
        JSON.stringify({
          active,
          collapsed: [...collapsed],
          compact: [...compact],
        }),
      );
    } catch {}
  };
  function render() {
    const focused = document.activeElement;
    for (const m of MONITORS) {
      const card = cards.get(m.id);
      (active.includes(m.id) ? stack : storage).append(card);
      card.dataset.density = compact.has(m.id) ? "compact" : "comfortable";
      card.querySelector(".monitor-body").hidden = collapsed.has(m.id);
      const button = card.querySelector(".monitor-collapse");
      button.textContent = collapsed.has(m.id) ? "▸" : "▾";
      button.setAttribute(
        "aria-label",
        `${collapsed.has(m.id) ? "Expand" : "Collapse"} ${m.title}`,
      );
      button.setAttribute("aria-expanded", String(!collapsed.has(m.id)));
    }
    for (const id of active) stack.append(cards.get(id));
    stack.classList.toggle("empty-stack", !active.length);
    save();
    if (focused?.isConnected && !focused.closest("[hidden]"))
      focused.focus({ preventScroll: true });
  }
  for (const m of MONITORS) {
    const card = node("section", "monitor-card");
    card.dataset.monitor = m.id;
    card.setAttribute("aria-label", m.title);
    const head = node("div", "monitor-heading"),
      collapse = node("button", "monitor-collapse", "▾"),
      title = node("h2", "", m.title),
      configure = node("button", "monitor-configure", "⋯"),
      remove = node("button", "monitor-remove", "×");
    configure.setAttribute("aria-label", `Configure ${m.title}`);
    remove.setAttribute("aria-label", `Remove ${m.title}`);
    head.append(collapse, title, configure, remove);
    const config = node("div", "monitor-options");
    config.hidden = true;
    const dense = node("label");
    const check = node("input");
    check.type = "checkbox";
    check.checked = compact.has(m.id);
    check.setAttribute("aria-label", `Compact ${m.title}`);
    dense.append(check, document.createTextNode(" Compact spacing"));
    const up = node("button", "", "Move up"),
      down = node("button", "", "Move down");
    config.append(dense, up, down);
    up.onclick = () => move(-1);
    down.onclick = () => move(1);
    function move(direction) {
      const index = active.indexOf(m.id),
        next = index + direction;
      if (next >= 0 && next < active.length) {
        [active[index], active[next]] = [active[next], active[index]];
        render();
      }
    }
    check.onchange = () => {
      check.checked ? compact.add(m.id) : compact.delete(m.id);
      render();
    };
    collapse.onclick = () => {
      collapsed.has(m.id) ? collapsed.delete(m.id) : collapsed.add(m.id);
      render();
    };
    configure.onclick = () => {
      config.hidden = !config.hidden;
      configure.setAttribute("aria-expanded", String(!config.hidden));
    };
    remove.onclick = () => {
      active = active.filter((id) => id !== m.id);
      render();
      add.focus();
    };
    const body = node("div", "monitor-body");
    body.append(contents[m.id]);
    card.append(head, config, body);
    cards.set(m.id, card);
  }
  const dialog = node("dialog", "monitor-catalog");
  dialog.setAttribute("aria-label", "Monitor catalog");
  const dialogHead = node("div", "catalog-heading");
  dialogHead.append(node("h2", "", "Add a monitor"));
  const close = node("button", "", "×");
  close.setAttribute("aria-label", "Close monitor catalog");
  dialogHead.append(close);
  const search = node("input");
  search.type = "search";
  search.placeholder = "Search signals, anatomy, mechanics…";
  search.setAttribute("aria-label", "Search monitor catalog");
  const tags = node("div", "catalog-tags");
  const results = node("div", "catalog-results");
  let selectedTag = "all";
  function list() {
    results.replaceChildren();
    for (const m of filterMonitors(search.value, selectedTag)) {
      const entry = node("div", "catalog-entry");
      entry.append(
        node("h3", "", m.title),
        node("p", "", m.description),
        node("small", "", m.tags.join(" · ")),
      );
      const button = node(
        "button",
        "",
        active.includes(m.id) ? "Added" : "Add pane",
      );
      button.disabled = active.includes(m.id);
      button.setAttribute("aria-label", `Add ${m.title}`);
      button.onclick = () => {
        active.push(m.id);
        render();
        dialog.close();
        cards.get(m.id).scrollIntoView({ block: "nearest" });
      };
      entry.append(button);
      results.append(entry);
    }
    if (!results.children.length)
      results.append(node("p", "muted", "No monitors match this search."));
  }
  for (const tag of ["all", ...new Set(MONITORS.flatMap((m) => m.tags))]) {
    const b = node("button", "", tag);
    b.setAttribute("aria-pressed", String(tag === "all"));
    b.onclick = () => {
      selectedTag = tag;
      for (const c of tags.children)
        c.setAttribute("aria-pressed", String(c === b));
      list();
    };
    tags.append(b);
  }
  dialog.append(dialogHead, search, tags, results);
  document.body.append(dialog);
  search.oninput = list;
  add.onclick = () => {
    list();
    dialog.showModal();
    search.focus();
  };
  close.onclick = () => dialog.close();
  dialog.addEventListener("close", () => add.focus());
  dialog.addEventListener("keydown", (e) => {
    if (e.key === "Escape") {
      e.preventDefault();
      dialog.close();
    }
  });
  render();
  return {
    show(id) {
      if (!contents[id] || (active.includes(id) && !collapsed.has(id))) return;
      if (!active.includes(id)) active.push(id);
      collapsed.delete(id);
      render();
    },
    cards,
  };
}
