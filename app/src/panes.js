// Panes float directly on the renderer. There is no panel, no sidebar surface
// and no backdrop: each pane carries its own minimal ground so it stays legible
// over arbitrary 3D content. Which panes are on screen is the reader's choice.
export function mountPanes(host, definitions) {
  const panes = new Map();
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
        panes.get(definition.id).pane.hidden = false;
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
    pane.hidden = !definition.open;
    const head = document.createElement("div");
    head.className = "pane-head";
    const title = document.createElement("h2");
    title.textContent = definition.title;
    const dismiss = document.createElement("button");
    dismiss.className = "pane-dismiss";
    dismiss.type = "button";
    dismiss.textContent = "×";
    dismiss.setAttribute("aria-label", `Remove ${definition.title}`);
    dismiss.onclick = () => {
      pane.hidden = true;
      if (!picker.hidden) renderPicker();
    };
    head.append(title, dismiss);
    const body = document.createElement("div");
    body.className = "pane-body";
    body.append(definition.content);
    pane.append(head, body);
    host.append(pane);
    panes.set(definition.id, { pane, body });
  }
  return {
    show(id) {
      const entry = panes.get(id);
      if (entry) entry.pane.hidden = false;
      if (!picker.hidden) renderPicker();
    },
    hide(id) {
      const entry = panes.get(id);
      if (entry) entry.pane.hidden = true;
    },
    visible: (id) => !panes.get(id)?.pane.hidden,
    body: (id) => panes.get(id)?.body,
  };
}
