// Picture-and-label tiles, driven entirely by a catalog: id, label, the slots an
// entry occupies, thumbnail and requirements. Exclusivity and dependency are
// data, never front-end rules.
//
// Exclusivity is slot-set intersection. An entry occupies one or more slots and
// displaces whatever already holds any of them, so a dress that occupies
// torso_base and legs turns off both the shirt and the trousers, while a hat and
// a pair of gloves never touch each other. A single-slot catalog is the same
// rule with sets of one, so environments behave exactly as before.
const svg = (body) =>
  `<svg viewBox="0 0 64 48" aria-hidden="true" focusable="false">${body}</svg>`;

const GENERIC = svg(
  `<rect x="16" y="10" width="32" height="30" rx="4" fill="#48605f" stroke="#9fbdb8" stroke-width="1.4"/>
   <path d="M16 20 H48" stroke="#9fbdb8" stroke-width="1.2"/>`,
);

const escape = (value) =>
  String(value ?? "").replace(/[&<>"']/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);

function picture(item) {
  if (item.thumbnail_url)
    return `<img src="${escape(item.thumbnail_url)}" alt="" loading="lazy">`;
  return GENERIC;
}

const el = (tag, className, text) => {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
};

export function mountTiles(host, { onChange }) {
  host.classList.add("tile-grid");
  let items = [], slots = [], title = "", instances = [], serial = 0;
  const chosen = new Map(); // slot -> id
  // A record may declare several slots; one is the degenerate case.
  const slotsOf = (item) =>
    (item.slots && item.slots.length ? item.slots : [item.slot || `slot:${item.id}`]);
  const homeOf = (item) => slotsOf(item)[0];
  const slotSpec = (id) => slots.find((s) => s.id === id) || { id, exclusive: true };
  const isOn = (item) => slotsOf(item).some((slot) => chosen.get(slot) === item.id);

  const met = (requires) =>
    (requires || []).every((rule) => (rule.any_of || []).includes(chosen.get(rule.slot)));
  const available = (item) =>
    met(item.requires) && slotsOf(item).every((slot) => met(slotSpec(slot).requires));

  const clear = (id) => {
    for (const [slot, held] of [...chosen]) if (held === id) chosen.delete(slot);
  };
  function wear(item) {
    // Everything the new entry's slots already hold comes off whole: a garment
    // is never left occupying half the slots it declares.
    for (const slot of slotsOf(item)) {
      const held = chosen.get(slot);
      if (held && held !== item.id) clear(held);
    }
    for (const slot of slotsOf(item)) chosen.set(slot, item.id);
  }

  // A choice can invalidate another slot's choice; drop what is no longer
  // available until nothing else falls, rather than leaving a stale selection.
  function settle() {
    for (let pass = 0; pass < slots.length + items.length; pass++) {
      const stale = items.find((item) => isOn(item) && !available(item));
      if (!stale) break;
      clear(stale.id);
    }
    instances = instances.filter((instance) => {
      const item = items.find((i) => i.id === instance.id);
      return item && available(item);
    });
    // A slot that has just become reachable takes the default the catalog
    // declares for it, rather than opening empty.
    for (const spec of slots) {
      if (chosen.has(spec.id) || !spec.default || !met(spec.requires)) continue;
      const fallback = items.find((i) => i.id === spec.default);
      if (fallback && available(fallback)) wear(fallback);
    }
  }

  // Three interaction kinds, told apart by the catalog rather than by anything
  // hardcoded: exclusive selections are tiles, per-slot configuration is a
  // dropdown that appears only when its requirements are met, and objects are
  // additive inserts reached through one menu.
  let insertOpen = false;
  function render() {
    host.replaceChildren();
    const named = slots.filter((s) => s.label).length > 1;
    const pick = items.filter((i) => i.kind !== "component" && i.kind !== "object");
    const config = items.filter((i) => i.kind === "component");
    const objects = items.filter((i) => i.kind === "object");

    const order = named ? slots.map((s) => s.id) : [null];
    for (const slot of order) {
      const group = slot === null ? pick : pick.filter((item) => homeOf(item) === slot);
      if (!group.length) continue;
      const spec = slot === null ? {} : slotSpec(slot);
      if (spec.label && spec.label.toLowerCase() !== (title || "").toLowerCase())
        host.append(el("p", "tile-slot", spec.label));
      const grid = el("div", "tile-row");
      for (const item of group) grid.append(tile(item));
      host.append(grid);
    }

    for (const spec of slots) {
      const group = config.filter((item) => homeOf(item) === spec.id);
      if (!group.length || !met(spec.requires)) continue;
      const usable = group.filter(available);
      if (!usable.length) continue;
      const row = el("label", "config-row");
      row.append(el("span", null, spec.label || spec.id));
      const select = document.createElement("select");
      select.dataset.slot = spec.id;
      if (!spec.required) select.append(new Option("—", ""));
      for (const item of usable) select.append(new Option(item.label, item.id));
      select.value = chosen.get(spec.id) || "";
      select.onchange = () => {
        const next = items.find((i) => i.id === select.value);
        chosen.delete(spec.id);
        if (next) wear(next);
        settle();
        render();
        emit();
      };
      row.append(select);
      host.append(row);
    }

    if (objects.length) renderInserts(objects);
  }

  // One plus button, not a row of them: it opens a menu of what can be put in
  // the scene, and choosing an entry inserts it. What is already there is a
  // list underneath, each line removable on its own.
  function renderInserts(objects) {
    const usable = objects.filter(available);
    if (usable.length) {
      const head = el("div", "insert-head");
      head.append(el("p", "tile-slot", slotSpec("objects").label || "Objects"));
      const add = el("button", "insert-add", "＋");
      add.type = "button";
      add.id = "insert-object";
      add.setAttribute("aria-label", "Insert an object into the scene");
      add.setAttribute("aria-expanded", String(insertOpen));
      const menu = el("div", "insert-menu");
      menu.id = "insert-menu";
      menu.hidden = !insertOpen;
      menu.setAttribute("role", "menu");
      // The catalog says which objects a session can actually instantiate and
      // which it can only draw. Both go in the menu, under headings, so the
      // difference is stated once rather than repeated on every line.
      const groups = usable.some((i) => i.insertable === false)
        ? [["Simulated", usable.filter((i) => i.insertable !== false)],
           ["Drawn only", usable.filter((i) => i.insertable === false)]]
        : [[null, usable]];
      for (const [heading, group] of groups) {
        if (!group.length) continue;
        if (heading) menu.append(el("p", "insert-group", heading));
        for (const item of group) {
          const entry = el("button", null, item.label);
          entry.type = "button";
          entry.dataset.insert = item.id;
          entry.setAttribute("role", "menuitem");
          entry.onclick = () => {
            instances.push({ uid: `${item.id}:${++serial}`, id: item.id });
            insertOpen = false;
            render();
            emit();
          };
          menu.append(entry);
        }
      }
      add.onclick = () => {
        insertOpen = !insertOpen;
        menu.hidden = !insertOpen;
        add.setAttribute("aria-expanded", String(insertOpen));
      };
      const wrap = el("div", "insert-wrap");
      wrap.append(add, menu);
      head.append(wrap);
      host.append(head);
    } else insertOpen = false;
    if (instances.length) {
      const list = el("ul", "instance-list");
      for (const instance of instances) {
        const entry = el("li");
        entry.dataset.instance = instance.uid;
        entry.append(el("span", null, items.find((i) => i.id === instance.id)?.label || instance.id));
        const remove = el("button", "instance-remove", "×");
        remove.type = "button";
        remove.setAttribute("aria-label", `Remove ${instance.id}`);
        remove.onclick = () => {
          instances = instances.filter((x) => x !== instance);
          render();
          emit();
        };
        entry.append(remove);
        list.append(entry);
      }
      host.append(list);
    }
  }
  function tile(item) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "tile";
    button.dataset.tile = item.id;
    button.dataset.slot = slotsOf(item).join(" ");
    button.innerHTML = `${picture(item)}<span>${escape(item.label)}</span>`;
    const rule = slotSpec(homeOf(item));
    const usable = available(item);
    button.disabled = !usable;
    if (!usable && rule.note) button.title = rule.note;
    const on = isOn(item);
    button.setAttribute("aria-pressed", String(on));
    button.classList.toggle("on", on);
    button.onclick = () => {
      if (!available(item)) return;
      if (on && !rule.required) clear(item.id);
      else wear(item);
      settle();
      render();
      emit();
    };
    return button;
  }
  const emit = () => onChange(active(), Object.fromEntries(chosen), instances.map((i) => i.id));

  const active = () => [...new Set(chosen.values())];

  return {
    active,
    get selection() { return Object.fromEntries(chosen); },
    get instances() { return instances.map((i) => ({ ...i })); },
    setItems(next, options = {}) {
      items = next;
      slots = options.slots || [];
      title = options.title || "";
      instances = [];
      insertOpen = false;
      chosen.clear();
      // Slot defaults are applied in catalog order, and only where the slot's
      // own requirements are already satisfied by the slots before it.
      for (const id of options.initial || []) {
        const item = items.find((i) => i.id === id);
        if (item && available(item)) wear(item);
      }
      settle();
      render();
    },
    select(ids) {
      chosen.clear();
      for (const id of ids) {
        const item = items.find((i) => i.id === id);
        if (item && available(item)) wear(item);
      }
      settle();
      render();
    },
  };
}
