// The left column: Materialization, Layers, Clothing, Environment, Simulation.
// Nothing here asks which source, candidate or execution owner to use.
import { mountTiles } from "./tiles.js";
import { DEFAULT_PALETTE, DEFAULT_SKIN_TONE, readPaletteChoice, writePaletteChoice } from "./palette.js";

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
// The same disclosure as the panes and the Layers list: one triangle, one
// meaning, on both sides of the screen, and it survives a reload.
const STORE = "ihm.sections.v1";
const shutSections = (() => {
  try {
    const saved = JSON.parse(localStorage.getItem(STORE));
    return new Set(Array.isArray(saved) ? saved : []);
  } catch {
    return new Set();
  }
})();
function saveSections() {
  try { localStorage.setItem(STORE, JSON.stringify([...shutSections])); } catch {}
}
// Per-system opacity. It used to be an inline slider on every row and was cut as
// clutter; it comes back behind one control per row, so a row at rest still
// carries only its triangle, checkbox, swatch, name and count. The value is what
// the renderer paints the system at, which is why skin starts translucent.
export const DEFAULT_OPACITY = { integumentary: 0.22 };
const OPACITY_STORE = "ihm.layer-opacity.v1";
const savedOpacity = new Map(Object.entries((() => {
  try {
    const saved = JSON.parse(localStorage.getItem(OPACITY_STORE));
    return saved && typeof saved === "object" ? saved : {};
  } catch {
    return {};
  }
})()).filter(([, v]) => Number.isFinite(v) && v >= 0 && v <= 1));
function saveOpacity() {
  try { localStorage.setItem(OPACITY_STORE, JSON.stringify(Object.fromEntries(savedOpacity))); } catch {}
}
export const layerOpacity = (name) =>
  savedOpacity.has(name) ? savedOpacity.get(name) : (DEFAULT_OPACITY[name] ?? 1);
function section(id, title, ...children) {
  const node = el("section", "column-section");
  node.dataset.section = id;
  const head = el("div", "section-head");
  const disclose = el("button", "disclose");
  disclose.type = "button";
  const body = el("div", "section-body");
  body.append(...children);
  head.append(disclose, el("h2", null, title));
  node.append(head, body);
  const paint = () => {
    const shut = shutSections.has(id);
    body.hidden = shut;
    node.dataset.collapsed = String(shut);
    disclose.textContent = shut ? "▸" : "▾";
    disclose.setAttribute("aria-expanded", String(!shut));
    disclose.setAttribute("aria-label", `${shut ? "Expand" : "Collapse"} ${title}`);
  };
  disclose.onclick = () => {
    shutSections.has(id) ? shutSections.delete(id) : shutSections.add(id);
    paint();
    saveSections();
  };
  paint();
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

  // One popover, reused by every row, anchored to whichever settings control
  // opened it. It is fixed to the window rather than parented to the row, so a
  // long Layers list scrolling under it cannot clip it, and it flips above the
  // row when there is no room below.
  const opacityButtons = new Map();
  let openOpacity = null;
  const popover = el("div", "float opacity-popover");
  popover.id = "layer-opacity";
  popover.hidden = true;
  popover.setAttribute("role", "dialog");
  const opacityTitle = el("p", "opacity-title");
  const opacityRange = el("input");
  opacityRange.type = "range";
  opacityRange.id = "layer-opacity-range";
  opacityRange.min = "0";
  opacityRange.max = "100";
  opacityRange.step = "1";
  const opacityValue = el("output", "opacity-value");
  opacityValue.id = "layer-opacity-value";
  popover.append(opacityTitle, opacityRange, opacityValue);
  document.body.append(popover);

  function placeFloat(node, button) {
    const box = button.getBoundingClientRect();
    node.style.visibility = "hidden";
    node.hidden = false;
    const width = node.offsetWidth, height = node.offsetHeight;
    const below = box.bottom + 6;
    const top = below + height > window.innerHeight - 8
      ? Math.max(8, box.top - height - 6) : below;
    const left = Math.max(8, Math.min(box.right - width, window.innerWidth - width - 8));
    node.style.top = `${top}px`;
    node.style.left = `${left}px`;
    node.style.visibility = "";
  }
  const placeOpacity = (button) => placeFloat(popover, button);
  function closeOpacity() {
    if (!openOpacity) return;
    openOpacity = null;
    popover.hidden = true;
    popover.removeAttribute("data-system");
    for (const button of opacityButtons.values()) button.setAttribute("aria-expanded", "false");
  }
  function showOpacity(name) {
    const button = opacityButtons.get(name);
    if (!button) return closeOpacity();
    // Only one at a time: opening a second row's control closes the first.
    for (const [id, other] of opacityButtons) other.setAttribute("aria-expanded", String(id === name));
    openOpacity = name;
    popover.dataset.system = name;
    opacityTitle.textContent = `${name.replaceAll("_", " ")} opacity`;
    popover.setAttribute("aria-label", `${name.replaceAll("_", " ")} opacity`);
    opacityRange.setAttribute("aria-label", `${name.replaceAll("_", " ")} opacity`);
    const percent = Math.round(layerOpacity(name) * 100);
    opacityRange.value = String(percent);
    opacityValue.textContent = `${percent}%`;
    placeOpacity(button);
    opacityRange.focus();
  }
  // Opacity is not visibility: this never touches the checkbox, and a system at
  // zero opacity is still a shown system.
  opacityRange.oninput = () => {
    if (!openOpacity) return;
    const percent = Number(opacityRange.value);
    opacityValue.textContent = `${percent}%`;
    savedOpacity.set(openOpacity, percent / 100);
    saveOpacity();
    hooks.onOpacity?.(openOpacity, percent / 100);
  };
  document.addEventListener("pointerdown", (event) => {
    if (!openOpacity) return;
    if (popover.contains(event.target) || opacityButtons.get(openOpacity)?.contains(event.target)) return;
    closeOpacity();
  });
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape" && openOpacity) {
      const button = opacityButtons.get(openOpacity);
      closeOpacity();
      button?.focus();
    }
  });
  for (const [target, type] of [[window, "resize"], [host, "scroll"]])
    target.addEventListener(type, () => {
      if (openOpacity) placeOpacity(opacityButtons.get(openOpacity));
      if (paletteOpen) placeFloat(palettePopover, paletteSettings);
    }, { passive: true });

  // Tissue colour. One control for the whole Layers section rather than one per
  // row, but the same idiom as the opacity gear beside it: a gear, an anchored
  // popover fixed to the window, and a choice that survives a reload. Palette
  // and skin tone are two axes and are two controls; the server publishes them
  // separately and nothing here folds one into the other.
  const paletteRow = el("div", "layer-row palette-row");
  paletteRow.id = "palette-row";
  const paletteCurrent = el("small", "palette-current", "Loading palettes…");
  paletteCurrent.id = "palette-current";
  const paletteSettings = el("button", "layer-settings", "⚙");
  paletteSettings.type = "button";
  paletteSettings.id = "palette-settings";
  paletteSettings.disabled = true;
  paletteSettings.setAttribute("aria-haspopup", "dialog");
  paletteSettings.setAttribute("aria-expanded", "false");
  paletteSettings.setAttribute("aria-label", "Tissue colour palette");
  paletteSettings.title = "Tissue colour palette";
  paletteRow.append(el("span", "palette-label", "Tissue colour"), paletteCurrent, paletteSettings);

  const palettePopover = el("div", "float palette-popover");
  palettePopover.id = "tissue-palette";
  palettePopover.hidden = true;
  palettePopover.setAttribute("role", "dialog");
  palettePopover.setAttribute("aria-label", "Tissue colour palette");
  const paletteSelect = el("select");
  paletteSelect.id = "palette-choice";
  paletteSelect.setAttribute("aria-label", "Palette");
  const paletteNote = el("p", "note");
  paletteNote.id = "palette-note";
  const skinToneSelect = el("select");
  skinToneSelect.id = "skin-tone-choice";
  skinToneSelect.setAttribute("aria-label", "Skin tone");
  const skinToneNote = el("p", "note");
  skinToneNote.id = "skin-tone-note";
  const palettePolicy = el("p", "note");
  palettePolicy.id = "palette-policy";
  palettePopover.append(
    el("p", "opacity-title", "Palette"), paletteSelect, paletteNote,
    el("p", "opacity-title", "Skin tone"), skinToneSelect, skinToneNote, palettePolicy);
  document.body.append(palettePopover);

  let paletteIndex = null, paletteOpen = false, skinToneApplies = false, paletteFallbacks = 0;
  let paletteError = "";
  const choice = readPaletteChoice();
  const option = (value, label) => {
    const node = document.createElement("option");
    node.value = value;
    node.textContent = label;
    return node;
  };
  function paintPalette() {
    const entry = (paletteIndex?.palettes || []).find((p) => p.id === choice.palette);
    const tone = paletteIndex?.skin_tone_options?.[choice.skin_tone];
    paletteCurrent.textContent = paletteError || entry?.label || choice.palette;
    palettePopover.dataset.palette = choice.palette;
    palettePopover.dataset.skinTone = choice.skin_tone;
    // The palette's own words. Each record states what it represents and what it
    // does not, and restating it here in shorter words is how a claim gets made
    // that the data did not make.
    const about = [entry?.description];
    if (paletteFallbacks)
      about.push(`${paletteFallbacks} structure${paletteFallbacks === 1 ? "" : "s"} in this scene are not in this palette and keep the colour they already had.`);
    if (paletteError) about.push(paletteError);
    paletteNote.textContent = about.filter(Boolean).join(" ");
    skinToneSelect.disabled = !skinToneApplies;
    const says = [tone?.description];
    if (tone?.cohort_as_reported) says.push(`Cohort as reported: ${tone.cohort_as_reported}.`);
    if (paletteIndex && !skinToneApplies)
      says.push(`${entry?.label || choice.palette} declares no skin tone, so this choice is not applied to it.`);
    skinToneNote.textContent = says.filter(Boolean).join(" ");
  }
  function closePalette() {
    if (!paletteOpen) return;
    paletteOpen = false;
    palettePopover.hidden = true;
    paletteSettings.setAttribute("aria-expanded", "false");
  }
  function openPalette() {
    if (!paletteIndex) return;
    closeOpacity();
    paletteOpen = true;
    paletteSettings.setAttribute("aria-expanded", "true");
    placeFloat(palettePopover, paletteSettings);
    paletteSelect.focus();
  }
  paletteSettings.onclick = () => (paletteOpen ? closePalette() : openPalette());
  const takeChoice = () => {
    choice.palette = paletteSelect.value;
    choice.skin_tone = skinToneSelect.value;
    writePaletteChoice(choice);
    paintPalette();
    hooks.onPalette?.(choice.palette, choice.skin_tone);
  };
  paletteSelect.onchange = takeChoice;
  skinToneSelect.onchange = takeChoice;
  document.addEventListener("pointerdown", (event) => {
    if (!paletteOpen) return;
    if (palettePopover.contains(event.target) || paletteSettings.contains(event.target)) return;
    closePalette();
  });
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape" && paletteOpen) { closePalette(); paletteSettings.focus(); }
  });

  // Clothing --------------------------------------------------------------
  const clothingHost = el("div");
  clothingHost.id = "clothing-tiles";
  const clothing = mountTiles(clothingHost, { onChange: (ids) => hooks.onClothing(ids) });
  const clothingNote = el("p", "note", "Loading the wardrobe…");
  clothingNote.id = "clothing-note";
  const clothingSection = section("clothing", "Clothing", clothingHost, clothingNote);
  clothingSection.id = "clothing-section";

  // Environment -----------------------------------------------------------
  const environmentHost = el("div");
  environmentHost.id = "environment-tiles";
  const environment = mountTiles(environmentHost, {
    onChange: (ids, selection, objects) => hooks.onEnvironment(selection, objects),
  });
  const environmentSection = section("environment", "Environment", environmentHost);
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
    section("materialization", "Materialization", materialization, materializationNote),
    section("layers", "Layers", paletteRow, search, layers),
    clothingSection,
    environmentSection,
    section("simulation", "Simulation", dynamics, run, runNote),
  );

  // -----------------------------------------------------------------------
  function renderLayers() {
    const names = [...new Set(rows.map((r) => r.system))].sort();
    layers.replaceChildren();
    opacityButtons.clear();
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
      // One small control at the right edge, and nothing else added to the row.
      const settings = el("button", "layer-settings", "\u2699");
      settings.type = "button";
      settings.dataset.system = name;
      settings.setAttribute("aria-haspopup", "dialog");
      settings.setAttribute("aria-expanded", String(openOpacity === name));
      settings.setAttribute("aria-label", `Opacity for ${name.replaceAll("_", " ")}`);
      settings.title = "Opacity";
      settings.onclick = () => (openOpacity === name ? closeOpacity() : showOpacity(name));
      opacityButtons.set(name, settings);
      row.append(disclose, label, el("small", null, String(members.length)), settings);
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
    // The rows were rebuilt: re-anchor an open popover to its new control, or
    // close it if its row is no longer listed.
    if (openOpacity) {
      const button = opacityButtons.get(openOpacity);
      if (!button) closeOpacity();
      else { button.setAttribute("aria-expanded", "true"); placeOpacity(button); }
    }
  }

  return {
    get systems() { return systems; },
    get hidden() { return hidden; },
    opacity(system) { return layerOpacity(system); },
    get options() { return { ...options }; },
    get environment() { return environment.selection.environment || null; },
    get environmentSelection() { return environment.selection; },
    get sceneObjects() { return environment.instances; },
    get garments() { return clothing.active(); },
    setGarments(items, options) { clothing.setItems(items, options); },
    setClothingNote(text) { clothingNote.textContent = text; clothingNote.hidden = !text; },
    setEnvironments(items, options) { environment.setItems(items, { title: "Environment", ...options }); },
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
    // The published index fills both lists and settles the default: the saved
    // choice is honoured only if the index still names it, and the index's own
    // default_palette wins otherwise. Returns the choice the caller should paint.
    setPalettes(index) {
      paletteIndex = index;
      paletteError = "";
      const palettes = index?.palettes || [];
      const tones = index?.skin_tone_options || {};
      paletteSelect.replaceChildren(...palettes.map((p) => option(p.id, p.label || p.id)));
      skinToneSelect.replaceChildren(...Object.entries(tones).map(([id, t]) => option(id, t.label || id)));
      if (!palettes.some((p) => p.id === choice.palette))
        choice.palette = index?.default_palette || DEFAULT_PALETTE;
      if (!(choice.skin_tone in tones)) choice.skin_tone = DEFAULT_SKIN_TONE;
      paletteSelect.value = choice.palette;
      skinToneSelect.value = choice.skin_tone;
      palettePolicy.textContent = index?.default_policy || "";
      paletteSettings.disabled = !palettes.length;
      paintPalette();
      return { ...choice };
    },
    // What actually happened when the choice was painted, reported back so the
    // popover can say it rather than assert it in advance.
    setPaletteApplied({ skinToneApplies: applies = false, fallbacks = 0, error = "" } = {}) {
      skinToneApplies = applies;
      paletteFallbacks = fallbacks;
      paletteError = error;
      paintPalette();
    },
    get paletteChoice() { return { ...choice }; },
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
