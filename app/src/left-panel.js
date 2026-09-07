// The left column: Materialization, Layers, Clothing, Environment, Simulation.
// Nothing here asks which source, candidate or execution owner to use.
import { mountTiles } from "./tiles.js";

const COLORS = {
  skeletal: "#d4c7ac", cardiac: "#ba6870", muscular: "#b66f67", nervous: "#d2ae57",
  arterial: "#bf716c", venous: "#8390b8", respiratory: "#7aa9a8", digestive: "#b59680",
  urinary: "#c19678", lymphatic: "#8cab78", integumentary: "#c4a18d", hair: "#8a6540",
  microvascular: "#cf6679", connective: "#9db3a4", sensory: "#c7b26a",
  endocrine: "#b58ab5", reproductive: "#c98fa2",
};

// Checked-and-locked dynamics are integral to the one native engine: it has no
// switch for them. Only the three that the runtime can actually gate are live.
export const DYNAMICS = [
  { id: "cardiovascular", label: "Cardiovascular", fixed: true },
  { id: "respiratory", label: "Respiratory", fixed: true },
  { id: "digestive", label: "Digestive", option: "intake_mass" },
  { id: "urinary", label: "Renal", fixed: true },
  { id: "endocrine", label: "Endocrine & metabolic", fixed: true },
  { id: "nervous", label: "Nervous & motor", fixed: true },
  { id: "integumentary", label: "Skin · regional", option: "regional_skin" },
  { id: "thermal", label: "Thermoregulation", fixed: true },
  { id: "hair", label: "Hair · elastic strands", option: "hair" },
];

const el = (tag, className, text) => {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
};
function section(title, ...children) {
  const node = el("section", "column-section");
  node.append(el("h2", null, title), ...children);
  return node;
}

export function mountLeftColumn(host, hooks) {
  const systems = new Set();
  const hidden = new Set();
  const expanded = new Set();
  const options = { regional_skin: false, intake_mass: false, hair: false };
  let rows = [], query = "";

  // Materialization -------------------------------------------------------
  const materialization = el("select");
  materialization.id = "materialization";
  materialization.setAttribute("aria-label", "Materialization");
  materialization.onchange = () => hooks.onMaterialization(materialization.value);
  const materializationNote = el("p", "note", "Loading materializations…");

  // Layers ----------------------------------------------------------------
  const search = el("input");
  search.type = "search";
  search.id = "layer-search";
  search.placeholder = "Find a member…";
  search.setAttribute("aria-label", "Find a member");
  const layers = el("div", "layers");
  layers.id = "layers";
  search.oninput = () => {
    query = search.value.trim().toLowerCase();
    renderLayers();
  };

  // Clothing --------------------------------------------------------------
  const clothingHost = el("div");
  clothingHost.id = "clothing-tiles";
  const clothing = mountTiles(clothingHost, { onChange: (ids) => hooks.onClothing(ids) });
  const clothingSection = section("Clothing", clothingHost);
  clothingSection.id = "clothing-section";

  // Environment -----------------------------------------------------------
  const environmentHost = el("div");
  environmentHost.id = "environment-tiles";
  const environment = mountTiles(environmentHost, { onChange: (ids) => hooks.onEnvironment(ids[0] || null) });
  const environmentSection = section("Environment", environmentHost);
  environmentSection.id = "environment-section";

  // Simulation ------------------------------------------------------------
  const dynamics = el("div", "dynamics");
  dynamics.id = "dynamics";
  const boxes = new Map();
  for (const item of DYNAMICS) {
    const label = el("label", "check-row");
    const box = el("input");
    box.type = "checkbox";
    box.checked = true;
    box.value = item.id;
    if (item.fixed) {
      box.disabled = true;
      label.title = "Integral to the one native engine; it has no switch for this.";
    } else {
      box.checked = options[item.option];
      box.onchange = () => {
        options[item.option] = box.checked;
        hooks.onSimulation(item.option, box.checked);
      };
    }
    label.append(box, el("span", null, item.label));
    boxes.set(item.id, box);
    dynamics.append(label);
  }
  const run = el("button", "run", "Start body");
  run.type = "button";
  run.id = "run-body";
  run.onclick = () => hooks.onRun();
  const runNote = el("p", "note", "");
  runNote.id = "run-note";
  runNote.setAttribute("role", "status");

  host.replaceChildren(
    section("Materialization", materialization, materializationNote),
    section("Layers", search, layers),
    clothingSection,
    environmentSection,
    section("Simulation", dynamics, run, runNote),
  );

  // -----------------------------------------------------------------------
  function renderLayers() {
    const names = [...new Set(rows.map((r) => r.system))].sort();
    layers.replaceChildren();
    for (const name of names) {
      const members = rows.filter((r) => r.system === name);
      const matching = query
        ? members.filter((m) => m.name.toLowerCase().includes(query))
        : members;
      if (query && !matching.length && !name.includes(query)) continue;
      const item = el("div", "layer");
      const row = el("div", "layer-row");
      const open = expanded.has(name) || (query && matching.length);
      const disclose = el("button", "disclose", open ? "▾" : "▸");
      disclose.type = "button";
      disclose.setAttribute("aria-expanded", String(!!open));
      disclose.setAttribute("aria-label", `${open ? "Collapse" : "Expand"} ${name}`);
      disclose.onclick = () => {
        expanded.has(name) ? expanded.delete(name) : expanded.add(name);
        renderLayers();
      };
      const label = el("label", "check-row");
      const box = el("input");
      box.type = "checkbox";
      box.value = name;
      box.checked = systems.has(name);
      box.onchange = () => {
        box.checked ? systems.add(name) : systems.delete(name);
        hooks.onLayers();
        renderLayers();
      };
      const swatch = el("i");
      swatch.style.background = COLORS[name] || "#91abb0";
      label.append(box, swatch, el("span", null, name.replaceAll("_", " ")));
      row.append(disclose, label, el("small", null, String(members.length)));
      item.append(row);
      if (open) {
        const list = el("div", "members");
        for (const member of matching.slice(0, 250)) {
          // The checkbox owns visibility; the name selects the structure.
          const row = el("div", "check-row member");
          const memberBox = el("input");
          memberBox.type = "checkbox";
          memberBox.value = member.id;
          memberBox.checked = !hidden.has(member.id);
          memberBox.disabled = !systems.has(name);
          memberBox.setAttribute("aria-label", `Show ${member.name}`);
          memberBox.onchange = () => {
            memberBox.checked ? hidden.delete(member.id) : hidden.add(member.id);
            hooks.onLayers();
          };
          const label = el("button", "member-name", member.name);
          label.type = "button";
          label.onclick = () => hooks.onMember?.(member.id);
          row.append(memberBox, label);
          list.append(row);
        }
        if (matching.length > 250)
          list.append(el("p", "note", "Refine the search to see more members."));
        if (!matching.length) list.append(el("p", "note", "No members match."));
        item.append(list);
      }
      layers.append(item);
    }
    if (!layers.children.length) layers.append(el("p", "note", "No members match."));
  }

  return {
    get systems() { return systems; },
    get hidden() { return hidden; },
    get options() { return { ...options }; },
    get environment() { return environment.active()[0] || null; },
    get garments() { return clothing.active(); },
    setGarments(items, initial) { clothing.setItems(items, initial); },
    setEnvironments(items, initial) { environment.setItems(items, initial); },
    setMaterializations(items, current) {
      materialization.replaceChildren();
      for (const item of items) {
        const option = document.createElement("option");
        option.value = item.value;
        option.textContent = item.label;
        materialization.append(option);
      }
      if (current) materialization.value = current;
    },
    setMaterializationNote(text) { materializationNote.textContent = text; },
    setStructures(next, initialSystems) {
      rows = next;
      if (initialSystems) {
        systems.clear();
        for (const name of initialSystems) systems.add(name);
      }
      renderLayers();
    },
    setWholeBody(on) {
      clothingSection.hidden = !on;
      environmentSection.hidden = !on;
    },
    setRun(label, disabled = false) {
      run.textContent = label;
      run.disabled = disabled;
    },
    setRunNote(text) { runNote.textContent = text; },
    lockDynamics(locked) {
      for (const item of DYNAMICS)
        if (!item.fixed) boxes.get(item.id).disabled = locked;
    },
    refresh: renderLayers,
  };
}
