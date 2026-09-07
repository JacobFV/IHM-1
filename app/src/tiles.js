// Picture-and-label tiles, driven entirely by a catalog: id, label, slot,
// thumbnail and requirements. Exclusivity and dependency are data, never
// front-end rules — tiles in one slot replace each other, different slots
// combine, and a tile whose requirements are unmet is shown unavailable rather
// than selectable-and-then-refused.
const svg = (body) =>
  `<svg viewBox="0 0 64 48" aria-hidden="true" focusable="false">${body}</svg>`;

// Representative drawings, used only while a catalog entry has no image of its own.
export const THUMBNAILS = {
  shirt: svg(
    `<path d="M22 9 L32 14 L42 9 L52 15 L47 22 L44 20 L44 41 L20 41 L20 20 L17 22 L12 15 Z" fill="#5f8f8c" stroke="#9fd3cd" stroke-width="1.4" stroke-linejoin="round"/>
     <path d="M25 10 q7 8 14 0" fill="none" stroke="#cfeae5" stroke-width="1.2"/>`,
  ),
  shorts: svg(
    `<path d="M18 10 H46 L44 40 H35 L32 24 L29 40 H20 Z" fill="#33506b" stroke="#9dc3e0" stroke-width="1.4" stroke-linejoin="round"/>
     <path d="M18 15 H46" fill="none" stroke="#cfe3f2" stroke-width="1.2"/>`,
  ),
};
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
  return THUMBNAILS[item.thumbnail || item.id] || GENERIC;
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
  const slotOf = (item) => item.slot || `slot:${item.id}`;
  const slotSpec = (id) => slots.find((s) => s.id === id) || { id, exclusive: true };

  const met = (requires) =>
    (requires || []).every((rule) => (rule.any_of || []).includes(chosen.get(rule.slot)));
  const available = (item) => met(item.requires) && met(slotSpec(slotOf(item)).requires);

  // A choice can invalidate another slot's choice; drop what is no longer
  // available until nothing else falls, rather than leaving a stale selection.
  function settle() {
    for (let pass = 0; pass < slots.length + items.length; pass++) {
      const stale = items.find((item) => chosen.get(slotOf(item)) === item.id && !available(item));
      if (!stale) break;
      chosen.delete(slotOf(stale));
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
      if (fallback && available(fallback)) chosen.set(spec.id, spec.default);
    }
  }

  // Three interaction kinds, told apart by the catalog rather than by anything
  // hardcoded: exclusive selections are tiles, per-slot configuration is a
  // dropdown that appears only when its requirements are met, and objects are
  // additive inserts in the one non-exclusive slot.
  function render() {
    host.replaceChildren();
    const named = slots.filter((s) => s.label).length > 1;
    const pick = items.filter((i) => i.kind !== "component" && i.kind !== "object");
    const config = items.filter((i) => i.kind === "component");
    const objects = items.filter((i) => i.kind === "object");

    const order = named ? slots.map((s) => s.id) : [null];
    for (const slot of order) {
      const group = slot === null ? pick : pick.filter((item) => slotOf(item) === slot);
      if (!group.length) continue;
      const spec = slot === null ? {} : slotSpec(slot);
      if (spec.label && spec.label.toLowerCase() !== (title || "").toLowerCase())
        host.append(el("p", "tile-slot", spec.label));
      const grid = el("div", "tile-row");
      for (const item of group) grid.append(tile(item));
      host.append(grid);
    }

    for (const spec of slots) {
      const group = config.filter((item) => slotOf(item) === spec.id);
      if (!group.length || !met(spec.requires)) continue;
      const available_ = group.filter(available);
      if (!available_.length) continue;
      const row = el("label", "config-row");
      row.append(el("span", null, spec.label || spec.id));
      const select = document.createElement("select");
      select.dataset.slot = spec.id;
      if (!spec.required) select.append(new Option("—", ""));
      for (const item of available_) select.append(new Option(item.label, item.id));
      select.value = chosen.get(spec.id) || "";
      select.onchange = () => {
        select.value ? chosen.set(spec.id, select.value) : chosen.delete(spec.id);
        settle();
        render();
        emit();
      };
      row.append(select);
      host.append(row);
    }

    if (objects.length) {
      const usable = objects.filter(available);
      if (usable.length) {
        host.append(el("p", "tile-slot", slotSpec("objects").label || "Objects"));
        const actions = el("div", "insert-actions");
        for (const item of usable) {
          const button = el("button", "insert", `+ Insert ${item.label}`);
          button.type = "button";
          button.dataset.insert = item.id;
          button.onclick = () => {
            instances.push({ uid: `${item.id}:${++serial}`, id: item.id });
            render();
            emit();
          };
          actions.append(button);
        }
        host.append(actions);
      }
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
  }
  function tile(item) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "tile";
    button.dataset.tile = item.id;
    button.dataset.slot = slotOf(item);
    button.innerHTML = `${picture(item)}<span>${escape(item.label)}</span>`;
    const home = slotOf(item), rule = slotSpec(home);
    const usable = available(item);
    button.disabled = !usable;
    if (!usable && rule.note) button.title = rule.note;
    const on = chosen.get(home) === item.id;
    button.setAttribute("aria-pressed", String(on));
    button.classList.toggle("on", on);
    button.onclick = () => {
      if (!available(item)) return;
      if (chosen.get(home) === item.id && !rule.required) chosen.delete(home);
      else chosen.set(home, item.id);
      settle();
      render();
      emit();
    };
    return button;
  }
  const emit = () => onChange(active(), Object.fromEntries(chosen), instances.map((i) => i.id));

  const active = () => [...chosen.values()];

  return {
    active,
    get selection() { return Object.fromEntries(chosen); },
    get instances() { return instances.map((i) => ({ ...i })); },
    setItems(next, options = {}) {
      items = next;
      slots = options.slots || [];
      title = options.title || "";
      instances = [];
      chosen.clear();
      // Slot defaults are applied in catalog order, and only where the slot's
      // own requirements are already satisfied by the slots before it.
      for (const id of options.initial || []) {
        const item = items.find((i) => i.id === id);
        if (item && available(item)) chosen.set(slotOf(item), id);
      }
      settle();
      render();
    },
    select(ids) {
      chosen.clear();
      for (const id of ids) {
        const item = items.find((i) => i.id === id);
        if (item && available(item)) chosen.set(slotOf(item), id);
      }
      settle();
      render();
    },
  };
}
