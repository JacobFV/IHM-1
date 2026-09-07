// Panes float directly on the renderer. There is no panel, no sidebar surface
// and no backdrop: each pane carries its own minimal ground so it stays legible
// over arbitrary 3D content. Which panes are on screen, and which are collapsed
// to their title, is the reader's choice and survives a reload.
const STORE = "ihm.panes.v1";

function load() {
  try {
    const saved = JSON.parse(localStorage.getItem(STORE));
    return {
      shown: Array.isArray(saved?.shown) ? saved.shown : null,
      collapsed: Array.isArray(saved?.collapsed) ? saved.collapsed : [],
    };
  } catch {
    return { shown: null, collapsed: [] };
  }
}

export function mountPanes(host, definitions) {
  const panes = new Map();
  const saved = load();
  const collapsed = new Set(saved.collapsed);
  const save = () => {
    try {
      localStorage.setItem(STORE, JSON.stringify({
        shown: definitions.filter((d) => !panes.get(d.id).pane.hidden).map((d) => d.id),
        collapsed: [...collapsed],
      }));
    } catch {}
  };

  const add = document.createElement("button");
  add.type = "button";
  add.className = "pane-add";
  add.id = "add-pane";
  add.textContent = "＋";
  add.setAttribute("aria-label", "Choose monitors");
  add.setAttribute("aria-expanded", "false");
  const picker = document.createElement("div");
  picker.className = "pane-picker";
  picker.id = "pane-picker";
  picker.hidden = true;
  const bar = document.createElement("div");
  bar.className = "pane-bar";
  bar.append(add, picker);
  host.append(bar);

  function renderPicker() {
    picker.replaceChildren();
    const missing = definitions.filter((d) => panes.get(d.id).pane.hidden);
    if (!missing.length) {
      const empty = document.createElement("p");
      empty.className = "note";
      empty.textContent = "Every monitor is already shown.";
      picker.append(empty);
      return;
    }
    for (const definition of missing) {
      const button = document.createElement("button");
      button.type = "button";
      button.textContent = definition.title;
      button.onclick = () => {
        show(definition.id);
        panes.get(definition.id).pane.scrollIntoView({ block: "nearest" });
        renderPicker();
      };
      picker.append(button);
    }
  }
  add.onclick = () => {
    const open = picker.hidden;
    if (open) renderPicker();
    picker.hidden = !open;
    add.setAttribute("aria-expanded", String(open));
  };
  document.addEventListener("pointerdown", (event) => {
    if (!picker.hidden && !bar.contains(event.target)) {
      picker.hidden = true;
      add.setAttribute("aria-expanded", "false");
    }
  });

  for (const definition of definitions) {
    const pane = document.createElement("section");
    pane.className = "pane";
    pane.dataset.pane = definition.id;
    pane.setAttribute("aria-label", definition.title);
    pane.hidden = saved.shown ? !saved.shown.includes(definition.id) : !definition.open;

    const head = document.createElement("div");
    head.className = "pane-head";
    // The same triangle as the Layers list, meaning the same thing.
    const disclose = document.createElement("button");
    disclose.className = "disclose";
    disclose.type = "button";
    const title = document.createElement("h2");
    title.textContent = definition.title;
    const dismiss = document.createElement("button");
    dismiss.className = "pane-dismiss";
    dismiss.type = "button";
    dismiss.textContent = "×";
    dismiss.setAttribute("aria-label", `Remove ${definition.title}`);
    dismiss.onclick = () => {
      pane.hidden = true;
      save();
      if (!picker.hidden) renderPicker();
    };
    head.append(disclose, title, dismiss);

    const body = document.createElement("div");
    body.className = "pane-body";
    body.append(definition.content);
    pane.append(head, body);
    host.append(pane);
    panes.set(definition.id, { pane, body, disclose });

    const paint = () => {
      const shut = collapsed.has(definition.id);
      body.hidden = shut;
      pane.dataset.collapsed = String(shut);
      disclose.textContent = shut ? "▸" : "▾";
      disclose.setAttribute("aria-expanded", String(!shut));
      disclose.setAttribute("aria-label", `${shut ? "Expand" : "Collapse"} ${definition.title}`);
    };
    disclose.onclick = () => {
      collapsed.has(definition.id) ? collapsed.delete(definition.id) : collapsed.add(definition.id);
      paint();
      save();
    };
    paint();
  }

  function show(id) {
    const entry = panes.get(id);
    if (!entry) return;
    entry.pane.hidden = false;
    save();
  }
  return {
    show(id) {
      show(id);
      if (!picker.hidden) renderPicker();
    },
    hide(id) {
      const entry = panes.get(id);
      if (entry) { entry.pane.hidden = true; save(); }
    },
    visible: (id) => !panes.get(id)?.pane.hidden,
    body: (id) => panes.get(id)?.body,
  };
}
