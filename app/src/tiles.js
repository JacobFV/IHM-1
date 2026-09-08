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
  `<rect x="16" y="10" width="32" height="30" rx="4" fill="#1c272b" stroke="#8ea3a9" stroke-width="1.4"/>
   <path d="M16 20 H48" stroke="#8ea3a9" stroke-width="1.2"/>`,
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

  // `state` is a slot -> id map, defaulting to the live selection. Everything
  // that asks whether an entry is wearable takes one, so the same rules can be
  // run against a hypothetical selection without disturbing the real one.
  const met = (requires, state = chosen) =>
    (requires || []).every((rule) => (rule.any_of || []).includes(state.get(rule.slot)));
  const available = (item, state = chosen) =>
    met(item.requires, state) && slotsOf(item).every((slot) => met(slotSpec(slot).requires, state));

  const clearFrom = (state, id) => {
    for (const [slot, held] of [...state]) if (held === id) state.delete(slot);
  };
  const clear = (id) => clearFrom(chosen, id);
  function wearInto(state, item) {
    // Everything the new entry's slots already hold comes off whole: a garment
    // is never left occupying half the slots it declares.
    for (const slot of slotsOf(item)) {
      const held = state.get(slot);
      if (held && held !== item.id) clearFrom(state, held);
    }
    for (const slot of slotsOf(item)) state.set(slot, item.id);
  }
  const wear = (item) => wearInto(chosen, item);

  // A requirement is a route, not a wall.
  //
  // Four of the six scenes name `floor` as the environment they sit on, and the
  // app opens on `bed`, so four of six tiles were drawn disabled with no
  // indication that one press elsewhere would light them. They are not
  // unavailable; they are one declared step away. This walks that step: the
  // entries the tile needs are collected, deepest first, on a copy of the
  // selection, and the press applies them in order before the tile itself.
  //
  // A rule that offers more than one way to be satisfied is left alone. Two
  // options is a choice the catalogue has declined to make, and making it here
  // would be the front end inventing a rule again.
  function route(item, state = new Map(chosen), seen = new Set()) {
    if (available(item, state)) return [];
    if (seen.has(item.id)) return null;
    seen.add(item.id);
    const rules = [...(item.requires || []),
                   ...slotsOf(item).flatMap((slot) => slotSpec(slot).requires || [])];
    const steps = [];
    for (const rule of rules) {
      if (met([rule], state)) continue;
      const options = (rule.any_of || []).map((id) => items.find((i) => i.id === id)).filter(Boolean);
      if (options.length !== 1) return null;
      const inner = route(options[0], state, seen);
      if (inner === null) return null;
      steps.push(...inner, options[0]);
      wearInto(state, options[0]);
    }
    return available(item, state) ? steps : null;
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
      const usable = group.filter((item) => available(item));
      if (!usable.length) continue;
      const row = el("label", "config-row");
      row.append(el("span", null, spec.label || spec.id));
      const select = document.createElement("select");
      select.dataset.slot = spec.id;
      if (!spec.required) select.append(new Option("—", ""));
      for (const item of usable) {
        const option = new Option(item.label, item.id);
        option.disabled = item.live_supported === false;
        if (option.disabled) option.text += " · unavailable for live body";
        select.append(option);
      }
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
    const usable = objects.filter((item) => available(item));
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
    const steps = route(item);
    const usable = steps !== null;
    button.disabled = !usable;
    if (!usable && rule.note) button.title = rule.note;
    // What the press will do before it does it, named, so the environment tile
    // changing under the reader is something they were told about.
    else if (steps.length) button.title = `Also selects ${steps.map((s) => s.label || s.id).join(", ")}`;
    const on = isOn(item);
    button.setAttribute("aria-pressed", String(on));
    button.classList.toggle("on", on);
    button.onclick = () => {
      const path = route(item);
      if (path === null) return;
      if (on && !rule.required) clear(item.id);
      else { for (const step of path) wear(step); wear(item); }
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
