// Picture-and-label tiles, driven entirely by a catalog: id, label, slot and an
// optional thumbnail URL. Exclusivity is data, not a front-end rule — tiles that
// share a slot are mutually exclusive, different slots combine freely.
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
  bed: svg(
    `<rect x="8" y="24" width="48" height="12" rx="3" fill="#5c7570" stroke="#a9cec5" stroke-width="1.4"/>
     <rect x="12" y="18" width="16" height="8" rx="3" fill="#cfe3dc" stroke="#a9cec5" stroke-width="1.2"/>
     <path d="M10 36 v6 M54 36 v6" stroke="#a9cec5" stroke-width="1.6" stroke-linecap="round"/>`,
  ),
  floor: svg(
    `<path d="M4 34 H60" stroke="#a9cec5" stroke-width="1.8" stroke-linecap="round"/>
     <path d="M10 34 L20 42 M24 34 L34 42 M38 34 L48 42" stroke="#5c7570" stroke-width="1.3"/>
     <circle cx="32" cy="18" r="6" fill="#5c7570" stroke="#a9cec5" stroke-width="1.3"/>
     <path d="M32 24 v8" stroke="#a9cec5" stroke-width="1.6" stroke-linecap="round"/>`,
  ),
  studio: svg(
    `<circle cx="32" cy="24" r="8" fill="#4b6b78" stroke="#a9cfdd" stroke-width="1.4"/>
     <circle cx="12" cy="12" r="1.6" fill="#cfe6ee"/><circle cx="52" cy="14" r="1.4" fill="#cfe6ee"/>
     <circle cx="48" cy="38" r="1.6" fill="#cfe6ee"/><circle cx="16" cy="36" r="1.3" fill="#cfe6ee"/>
     <circle cx="32" cy="8" r="1.2" fill="#cfe6ee"/>`,
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

export function mountTiles(host, { onChange }) {
  host.classList.add("tile-grid");
  let items = [];
  const selected = new Set();
  const buttons = new Map();
  const slotOf = (item) => item.slot || `slot:${item.id}`;

  function render() {
    host.replaceChildren();
    buttons.clear();
    for (const item of items) {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "tile";
      button.dataset.tile = item.id;
      button.innerHTML = `${picture(item)}<span>${escape(item.label)}</span>`;
      button.onclick = () => {
        const on = !selected.has(item.id);
        if (on)
          for (const other of items)
            if (other !== item && slotOf(other) === slotOf(item)) selected.delete(other.id);
        on ? selected.add(item.id) : selected.delete(item.id);
        paint();
        onChange(active());
      };
      buttons.set(item.id, button);
      host.append(button);
    }
    paint();
  }
  function paint() {
    for (const [id, button] of buttons) {
      const on = selected.has(id);
      button.setAttribute("aria-pressed", String(on));
      button.classList.toggle("on", on);
    }
  }
  const active = () => items.filter((i) => selected.has(i.id)).map((i) => i.id);
  return {
    active,
    setItems(next, initial) {
      items = next;
      const known = new Set(items.map((i) => i.id));
      for (const id of [...selected]) if (!known.has(id)) selected.delete(id);
      if (initial) {
        selected.clear();
        for (const id of initial) if (known.has(id)) selected.add(id);
      }
      render();
    },
    select(ids) {
      selected.clear();
      for (const id of ids) if (buttons.has(id)) selected.add(id);
      paint();
    },
  };
}
