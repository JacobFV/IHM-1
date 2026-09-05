// Drawer visibility changes layout, never physiology or camera state.
export function mountWorkspace() {
  const $ = (id) => document.getElementById(id),
    workspace = $("workspace"),
    media = matchMedia("(max-width: 1000px)");
  const states = new Map(),
    scroll = new Map();
  function state() {
    const key = media.matches ? "compact" : "desktop";
    if (!states.has(key)) {
      let saved;
      try {
        saved = JSON.parse(localStorage.getItem(`ihm.workspace.cards.${key}`));
      } catch {}
      const value = {
        library: !media.matches,
        inspector: !media.matches,
        focus: false,
      };
      for (const k of Object.keys(value))
        if (typeof saved?.[k] === "boolean") value[k] = saved[k];
      if (media.matches && value.library && value.inspector)
        value.inspector = false;
      states.set(key, value);
    }
    return states.get(key);
  }
  function render() {
    const s = state();
    for (const name of ["library", "inspector"]) {
      const panel = $(`${name}-panel`),
        show = s[name] && !s.focus;
      if (!panel.hidden && !show) scroll.set(name, panel.scrollTop);
      const before = panel.hidden;
      panel.hidden = !show;
      workspace.dataset[name] = String(show);
      $(`toggle-${name}`).setAttribute("aria-expanded", String(show));
      if (show && before) panel.scrollTop = scroll.get(name) || 0;
    }
    $("toggle-focus").setAttribute("aria-pressed", String(s.focus));
    try {
      localStorage.setItem(
        `ihm.workspace.cards.${media.matches ? "compact" : "desktop"}`,
        JSON.stringify(s),
      );
    } catch {}
  }
  for (const name of ["library", "inspector"])
    $(`toggle-${name}`).onclick = () => {
      const s = state();
      s[name] = s.focus || !s[name];
      s.focus = false;
      if (media.matches && s[name])
        s[name === "library" ? "inspector" : "library"] = false;
      render();
    };
  $("toggle-focus").onclick = () => {
    state().focus = !state().focus;
    render();
  };
  document.addEventListener("keydown", (e) => {
    if (
      e.key !== "Escape" ||
      !media.matches ||
      document.querySelector("dialog[open]")
    )
      return;
    const name = ["library", "inspector"].find((n) => !$(`${n}-panel`).hidden);
    if (name) {
      state()[name] = false;
      render();
      $(`toggle-${name}`).focus();
    }
  });
  media.addEventListener("change", render);
  render();
  for (const name of ["library", "inspector"]) {
    const panel = $(`${name}-panel`),
      handle = document.createElement("div");
    handle.className = `panel-resizer resize-${name}`;
    handle.tabIndex = 0;
    handle.setAttribute("role", "separator");
    handle.setAttribute(
      "aria-label",
      `Resize ${name === "library" ? "view" : "monitors"} panel`,
    );
    handle.setAttribute("aria-orientation", "vertical");
    workspace.append(handle);
    function size(width) {
      const next = Math.max(220, Math.min(workspace.clientWidth * 0.4, width));
      workspace.style.setProperty(`--${name}-preferred-width`, `${next}px`);
      handle.setAttribute("aria-valuenow", String(Math.round(next)));
    }
    handle.onkeydown = (e) => {
      if (!["ArrowLeft", "ArrowRight"].includes(e.key)) return;
      e.preventDefault();
      size(
        panel.clientWidth +
          (e.key === "ArrowRight" ? 16 : -16) * (name === "library" ? 1 : -1),
      );
    };
    handle.onpointerdown = (e) => {
      e.preventDefault();
      handle.setPointerCapture(e.pointerId);
      const x = e.clientX,
        width = panel.clientWidth;
      const move = (e) =>
        size(width + (e.clientX - x) * (name === "library" ? 1 : -1));
      handle.addEventListener("pointermove", move);
      handle.addEventListener(
        "lostpointercapture",
        () => handle.removeEventListener("pointermove", move),
        { once: true },
      );
      handle.addEventListener(
        "pointerup",
        (e) => handle.releasePointerCapture(e.pointerId),
        { once: true },
      );
    };
  }
}
