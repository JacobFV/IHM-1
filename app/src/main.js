import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import {
  filterStructures,
  geometryArrays,
  chartPath,
  scenarioInput,
  spectralSeries,
  trajectoryFromChannels,
  scalarCoordinates,
} from "./state.js";
import "./style.css";
const $ = (id) => document.getElementById(id),
  esc = (s) =>
    String(s ?? "").replace(
      /[&<>"']/g,
      (c) =>
        ({
          "&": "&amp;",
          "<": "&lt;",
          ">": "&gt;",
          '"': "&quot;",
          "'": "&#39;",
        })[c],
    );
document.querySelector("#app").innerHTML =
  `<header><div class="brand"><span class="brandmark">ih<span>•</span></span><div>INTEGRATED HUMAN<small>Computational physiology workbench</small></div></div><div class="header-right"><span class="live-dot"></span> LOCAL RESEARCH ENVIRONMENT <span class="version">v0.3</span></div></header><main><aside class="library"><div class="section-heading"><span>ANATOMY LIBRARY</span><span id="total">—</span></div><label class="field-label" for="model">Source model</label><select id="model"><option>Loading source families…</option></select><p class="muted" id="model-note">Reading the source manifest.</p><div class="search"><span>⌕</span><input id="search" placeholder="Find a structure…" aria-label="Search anatomy"></div><div class="section-heading layer-title">SYSTEM LAYERS <button id="all" class="text-button">Show all</button></div><div id="systems"></div><div class="section-heading results-heading">STRUCTURES <span id="count">0</span></div><div id="structures" class="structures"></div><div class="library-footer"><span class="live-dot"></span> Source geometry · independent frames</div></aside><section class="center"><div class="viewport" id="viewport"><div class="view-heading"><div class="eyebrow">SPATIAL EXPLORER</div><h1 id="view-title">Human anatomy</h1><p id="frame-label">Source-defined coordinates</p></div><div class="view-tools"><button id="reset" title="Reset camera">↺ <span>Reset view</span></button><button id="posture" title="Rigid display rotation only; native posture is not established">Supine view</button><button id="front">Anterior</button><button id="side">Lateral</button></div><div id="scene-status" role="status">Loading anatomical assets…</div><div id="flow-legend" class="flow-legend" hidden></div><div class="orientation">Y <span>↑</span><br><small>X →</small></div><div class="view-footer"><span id="render-count">0 structures visible</span><span>Drag to orbit · Scroll to zoom · Click to inspect</span></div></div><div class="display-controls"><label>Opacity <input id="opacity" type="range" min=".05" max="1" step=".05" value="1"><output id="opacity-value">100%</output></label><label>Section <input id="clip" type="range" min="0" max="100" value="100"><output id="clip-value">Off</output></label><select id="flow-field" aria-label="Vascular field"><option value="velocity">Velocity</option><option value="pressure">Pressure</option></select><button id="play" disabled aria-label="Play spatial time series">▶</button><input id="time" type="range" min="0" max="0" value="0" aria-label="Spatial time frame" disabled><output id="time-value">No flow frames</output></div><section class="signals"><div class="signal-header"><div><div class="eyebrow">TEMPORAL OBSERVATORY</div><h2>Physiology & dynamics</h2></div><div class="tabs"><button id="tab-phys" class="active">Trajectory</button><button id="tab-spectral">Spectrum</button></div></div><div id="spectral-controls" class="spectral-controls" hidden><select id="spectral-run" aria-label="Spectral evidence source"></select><select id="spectral-mode" aria-label="Spectral representation"><option value="psd">Fourier power density</option><option value="laplace">Finite Laplace magnitude</option></select><select id="sigma" aria-label="Laplace damping" hidden></select></div><div id="trajectory-controls" class="spectral-controls"><select id="trajectory-run" aria-label="Recorded physiology run"><option value="baseline">Recorded resting baseline</option></select></div><div class="signal-select"><select id="variable" aria-label="Physiology variable"></select><span id="signal-source">Awaiting simulation data</span></div><div id="chart" class="chart"><p class="empty">Loading physiological trajectories…</p></div><p id="chart-note" class="muted">Native simulation evidence is distinct from human measurements.</p></section></section><aside class="inspector"><div class="section-heading">STRUCTURE INSPECTOR <span>↗</span></div><div id="details"><div class="inspector-symbol">◎</div><h2>Explore the body</h2><p class="muted">Select a structure in the viewport or library to inspect its source, coordinate frame, and evidence status.</p></div><div class="scenario"><div class="eyebrow">NATIVE PHYSIOLOGY</div><h2>Run a scenario</h2><p class="muted">Execute the upstream BioGears engine and inspect its recorded response.</p><form id="scenario-form"><label class="field-label" for="engine-variant">Native implementation</label><select id="engine-variant"><option value="">Loading implementations…</option></select><p id="variant-note" class="muted">Source corrections do not establish clinical calibration.</p><label class="field-label" for="patient">Upstream patient</label><select id="patient"><option value="StandardMale">StandardMale</option></select><label class="field-label" for="scenario">Protocol</label><select id="scenario"><option value="baseline">Resting baseline</option><option value="exercise">Exercise & recovery</option><option value="hemorrhage_saline">Hemorrhage, saline & recovery</option></select><p id="protocol-note" class="muted">Record the upstream resting state.</p><label class="field-label" for="duration">Duration · seconds</label><input id="duration" type="number" min="1" max="600" value="60" required><label class="field-label" for="ambient">Ambient temperature · °C · optional</label><input id="ambient" type="number" min="10" max="35" step="any" placeholder="Upstream conditions"><label class="field-label" for="clothing">Clothing · clo · optional</label><input id="clothing" type="number" min="0" max="3" step="any" placeholder="Upstream conditions"><p class="muted">Blank values preserve upstream environment and clothing.</p><button class="primary" id="run" type="submit">Run native simulation <span>↗</span></button></form><div id="run-status" class="run-status" role="status">Checking native backend…</div></div><div class="evidence"><div class="section-heading">EVIDENCE BOUNDARY</div><p>Geometry describes anatomy. Simulation describes model behavior. Neither alone establishes patient-specific calibration.</p><div id="evidence-status" class="tag">Uncertainty remains explicit</div><details class="evidence-fold"><summary>Human measurement fit</summary><div id="calibration-details"><p>Loading calibration evidence…</p></div></details><details class="evidence-fold"><summary>System coverage</summary><div id="coverage-details"><p>Coverage data unavailable.</p></div></details><details class="evidence-fold"><summary>Conservative exchange</summary><div id="coupling-details"><p>Exchange data unavailable.</p></div></details><details class="evidence-fold"><summary>Vascular CFD audit</summary><div id="vascular-audit"><p>Audit unavailable.</p></div></details></div></aside></main>`;
let manifest,
  modelId,
  systems = new Set(),
  selected = null,
  objects = new Map(),
  generation = 0,
  loadController,
  modelBounds = new THREE.Box3(),
  activeRun = "baseline",
  phys = {},
  temporal = {},
  reproductive = null,
  csfAvailable = false,
  thermalIndex = null,
  spectral = false,
  flowFrames = 0,
  playing = false,
  lastTick = 0;
const colors = {
  skeletal: "#d4c7ac",
  cardiac: "#ba6870",
  muscular: "#b66f67",
  nervous: "#d2ae57",
  circulatory: "#bf716c",
  cardiovascular: "#bf716c",
  respiratory: "#7aa9a8",
  digestive: "#b59680",
  urinary: "#c19678",
  lymphatic: "#8cab78",
  integumentary: "#c4a18d",
};
const viewport = $("viewport");
let renderer, scene, camera, controls, group, webglError;
try {
  renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.localClippingEnabled = true;
  renderer.setClearColor(0x10191c, 0);
  viewport.prepend(renderer.domElement);
  scene = new THREE.Scene();
  camera = new THREE.PerspectiveCamera(35, 1, 0.01, 100000);
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  scene.add(new THREE.HemisphereLight(0xd7f0ee, 0x26333d, 2));
  const light = new THREE.DirectionalLight(0xfff1da, 3);
  light.position.set(3, 4, 5);
  scene.add(light);
  const rim = new THREE.DirectionalLight(0x87c5ce, 2);
  rim.position.set(-3, 2, -4);
  scene.add(rim);
  group = new THREE.Group();
  scene.add(group);
  new ResizeObserver(() => {
    const w = viewport.clientWidth,
      h = viewport.clientHeight;
    renderer.setSize(w, h);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }).observe(viewport);
  renderer.setAnimationLoop((now) => {
    controls.update();
    if (playing && now - lastTick > 160) {
      $("time").value = (Number($("time").value) + 1) % flowFrames;
      updateFrame();
      lastTick = now;
    }
    renderer.render(scene, camera);
  });
  const ray = new THREE.Raycaster();
  let pointer;
  renderer.domElement.addEventListener("pointerdown", (e) => {
    pointer = [e.clientX, e.clientY];
  });
  renderer.domElement.addEventListener("pointerup", (e) => {
    if (
      !pointer ||
      Math.hypot(e.clientX - pointer[0], e.clientY - pointer[1]) > 5
    )
      return;
    const r = renderer.domElement.getBoundingClientRect();
    ray.setFromCamera(
      new THREE.Vector2(
        ((e.clientX - r.left) / r.width) * 2 - 1,
        (-(e.clientY - r.top) / r.height) * 2 + 1,
      ),
      camera,
    );
    ray.params.Line.threshold =
      modelBounds.getSize(new THREE.Vector3()).length() * 0.003;
    const hit = ray.intersectObjects(
      [...objects.values()].filter((x) => x.visible),
      true,
    )[0];
    if (hit) {
      let obj = hit.object;
      while (!obj.userData.structure && obj.parent) obj = obj.parent;
      selectStructure(obj.userData.structure);
      if (obj.userData.cellIds && hit.face) {
        obj.userData.selectedCell = obj.userData.cellIds[hit.face.a];
        updateFrame();
      }
    }
  });
} catch (e) {
  webglError = e.message;
  $("scene-status").textContent =
    `WebGL unavailable: ${e.message}. The anatomy library and data panels remain available.`;
}
async function api(path, options) {
  const r = await fetch(path, options);
  if (!r.ok) {
    let msg;
    try {
      msg = await r.json();
    } catch {}
    throw Error(msg?.error || `Request failed (${r.status})`);
  }
  return r.json();
}
function model() {
  return manifest?.models.find((x) => x.id === modelId);
}
function resetCamera(direction = "front") {
  if (!camera || modelBounds.isEmpty()) return;
  const center = modelBounds.getCenter(new THREE.Vector3()),
    size = modelBounds.getSize(new THREE.Vector3()),
    distance = Math.max(size.y, size.x, size.z) * 1.85;
  controls.target.copy(center);
  camera.position
    .copy(center)
    .add(
      direction === "side"
        ? new THREE.Vector3(distance, 0, 0)
        : new THREE.Vector3(0, 0, distance),
    );
  camera.near = distance / 1000;
  camera.far = distance * 20;
  camera.updateProjectionMatrix();
  controls.update();
}
function rebuildSystems() {
  const rows = manifest.structures.filter((x) => x.model_id === modelId);
  const names = [...new Set(rows.map((x) => x.system))];
  systems = new Set(
    names.filter((x) => /skelet|respirat|cardi|circul|digest|urinar/.test(x)),
  );
  if (!systems.size || modelId.startsWith("opensim")) systems = new Set(names);
  $("systems").innerHTML = names
    .map(
      (name) =>
        `<label class="system-row"><input type="checkbox" value="${esc(name)}" ${systems.has(name) ? "checked" : ""}><i style="background:${colors[name] || "#91abb0"}"></i><span>${esc(name.replaceAll("_", " "))}</span><small>${rows.filter((x) => x.system === name).length}</small></label>`,
    )
    .join("");
  $("systems")
    .querySelectorAll("input")
    .forEach(
      (el) =>
        (el.onchange = () => {
          el.checked ? systems.add(el.value) : systems.delete(el.value);
          refresh();
        }),
    );
}
function chooseModel() {
  modelId = $("model").value;
  selected = null;
  loadController?.abort();
  objects.forEach(dispose);
  objects.clear();
  group?.clear();
  if (group) group.rotation.x = 0;
  $("posture").textContent = "Supine view";
  $("details").innerHTML =
    '<h2>Explore this model</h2><p class="muted">Select a structure to inspect source evidence.</p>';
  const m = model();
  $("opacity").value = modelId.startsWith("vascular") ? ".25" : "1";
  $("clip").value = "100";
  $("view-title").textContent = m.name;
  $("frame-label").textContent =
    `${m.frame || "Source frame"} · ${m.display_units || m.source_units || "source units"} display`;
  $("model-note").textContent =
    m.description ||
    "Independent source family. Display placement is not physiological registration.";
  modelBounds.makeEmpty();
  if (m.bounds) {
    const b = m.bounds;
    modelBounds.set(
      new THREE.Vector3(...(b.min || b[0])),
      new THREE.Vector3(...(b.max || b[1])),
    );
  }
  rebuildSystems();
  resetCamera();
  refresh();
}
function dispose(o) {
  o.traverse((x) => {
    x.geometry?.dispose();
    if (x.material) {
      (Array.isArray(x.material) ? x.material : [x.material]).forEach((m) =>
        m.dispose(),
      );
    }
  });
  o.removeFromParent();
}
function refresh() {
  const rows = filterStructures(
    manifest.structures,
    modelId,
    systems,
    $("search").value,
  );
  $("count").textContent = rows.length;
  $("structures").innerHTML =
    rows
      .slice(0, 350)
      .map(
        (s) =>
          `<button class="structure ${selected?.id === s.id ? "selected" : ""}" data-id="${esc(s.id)}"><span>${esc(s.name)}</span><small>↗</small></button>`,
      )
      .join("") +
    (rows.length > 350
      ? '<p class="muted">Refine search to see more structures.</p>'
      : "");
  $("structures")
    .querySelectorAll("button")
    .forEach(
      (b) =>
        (b.onclick = () =>
          selectStructure(
            manifest.structures.find((x) => x.id === b.dataset.id),
          )),
    );
  loadVisible(rows);
}
async function loadVisible(rows) {
  if (!group) {
    $("render-count").textContent = "3D rendering unavailable";
    $("scene-status").textContent =
      `WebGL unavailable: ${webglError || "No graphics context"}. Use the library to inspect anatomy.`;
    return;
  }
  const current = ++generation;
  loadController?.abort();
  loadController = new AbortController();
  const signal = loadController.signal,
    wanted = new Set(rows.map((x) => x.id));
  objects.forEach((o, id) => {
    o.visible = wanted.has(id);
  });
  let queue = rows
      .filter((x) => !objects.has(x.id))
      .sort(
        (a, b) =>
          Number(
            /femur|tibia|humerus|hip bone|sternum|skull|lung|liver|stomach|kidney/i.test(
              b.name,
            ),
          ) -
          Number(
            /femur|tibia|humerus|hip bone|sternum|skull|lung|liver|stomach|kidney/i.test(
              a.name,
            ),
          ),
      ),
    done = rows.length - queue.length,
    failures = 0;
  function status() {
    if (current !== generation) return;
    $("render-count").textContent =
      `${done} / ${rows.length} structures visible`;
    $("scene-status").textContent =
      done === rows.length
        ? ""
        : `Loading anatomy · ${done} / ${rows.length}${failures ? ` · ${failures} unavailable` : ""}`;
  }
  status();
  await Promise.all(
    Array.from({ length: 6 }, async () => {
      while (queue.length && current === generation) {
        const s = queue.shift();
        try {
          const g = await api(
            s.geometry_url || `/api/geometry/${encodeURIComponent(s.id)}`,
            { signal },
          );
          if (current !== generation) return;
          const object = createGeometry(s, g);
          if (object && group) {
            objects.set(s.id, object);
            group.add(object);
            if (!model()?.bounds) {
              modelBounds.expandByObject(object);
              if (done === 0) resetCamera();
            }
          }
          done++;
          status();
        } catch (e) {
          if (e.name === "AbortError") return;
          failures++;
          status();
        }
      }
    }),
  );
  if (current === generation) {
    if (!model()?.bounds) resetCamera();
    $("scene-status").textContent = failures
      ? `${failures} structures unavailable. Other source geometry remains interactive.`
      : rows.length
        ? ""
        : "No structures match these filters.";
    updateDisplay();
    setupFrames();
  }
}
function createGeometry(s, g) {
  if (!group) return null;
  const { positions, indices } = geometryArrays(g),
    geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  if (indices) geometry.setIndex(new THREE.BufferAttribute(indices, 1));
  const color = s.color || colors[s.system] || "#9eb7b7";
  let object;
  if (s.kind === "vectors") {
    object = new THREE.Group();
    const vectors = (g.velocities || g.vectors || []).flat();
    const scale = g.vector_scale || 1;
    for (let i = 0; i < positions.length; i += 3) {
      const v = new THREE.Vector3(...vectors.slice(i, i + 3));
      if (!Number.isFinite(v.length()) || !v.length()) continue;
      const arrow = new THREE.ArrowHelper(
        v.clone().normalize(),
        new THREE.Vector3(...positions.slice(i, i + 3)),
        v.length() * scale,
        color,
      );
      arrow.userData.nodeIndex = i / 3;
      object.add(arrow);
    }
    geometry.dispose();
    object.userData.vectorData = { positions, vectors, scale };
  } else if (s.kind === "lines") {
    object = new (g.topology === "polyline" ? THREE.Line : THREE.LineSegments)(
      geometry,
      new THREE.LineBasicMaterial({ color, transparent: true }),
    );
  } else if (s.kind === "scalar_mesh") {
    geometry.setAttribute(
      "color",
      new THREE.Float32BufferAttribute(new Float32Array(positions.length), 3),
    );
    object = new THREE.Mesh(
      geometry,
      new THREE.MeshBasicMaterial({
        vertexColors: true,
        side: THREE.DoubleSide,
      }),
    );
    object.userData.cellIds = g.cell_ids;
    object.userData.scalarFields = g.scalar_fields;
  } else {
    geometry.computeVertexNormals();
    object = new THREE.Mesh(
      geometry,
      new THREE.MeshStandardMaterial({
        color,
        roughness: 0.65,
        metalness: 0.05,
        side: THREE.DoubleSide,
        transparent: true,
      }),
    );
  }
  object.userData.structure = s;
  object.userData.frames = g.frames || [];
  object.userData.times = g.times || [];
  object.userData.pressure = g.pressure || [];
  let pmin = Infinity,
    pmax = -Infinity;
  for (const frame of [g, ...(g.frames || [])])
    for (const p of frame.pressure || []) {
      if (Number.isFinite(p)) {
        pmin = Math.min(pmin, p);
        pmax = Math.max(pmax, p);
      }
    }
  object.userData.pressureRange = [pmin, pmax];
  object.userData.units = g.units || "source units";
  return object;
}
function selectStructure(s) {
  if (!s) return;
  selected = s;
  objects.forEach((obj, id) =>
    obj.traverse((x) => {
      if (x.material?.emissive)
        x.material.emissive.set(id === s.id ? "#1d514f" : "#000000");
    }),
  );
  const source = s.source || {};
  $("details").innerHTML =
    `<div class="eyebrow">${esc(s.system)} / ${esc(s.kind || "mesh")}</div><h2>${esc(s.name)}</h2><span class="tag">${esc(s.calibration_status || "Not independently calibrated")}</span><dl><dt>Source</dt><dd>${esc(source.label || model()?.name)}</dd><dt>Coordinate frame</dt><dd>${esc(source.frame || model()?.frame || "Source-defined")}</dd><dt>Specimen</dt><dd>${esc(source.specimen || "See source metadata")}</dd><dt>Units</dt><dd>${esc(source.units || model()?.source_units || "Source units")}</dd><dt>Evidence status</dt><dd>${esc(source.status || s.status || "Source geometry; physiological registration unverified")}</dd>${s.path_interpretation ? `<dt>Muscle path interpretation</dt><dd>${esc(s.path_interpretation)}</dd><dt>Wrapping solution</dt><dd>${s.wrap_solved ? "Solved by source" : "Unsolved; attachment chords only"}</dd>` : ""}${
      s.native_mechanics
        ? `<dt>Native mechanics · source simulation</dt><dd>Activation ${esc(s.native_mechanics.activation)} · path length ${Number(s.native_mechanics.length_m).toFixed(4)} m<br>Tendon force ${Number(s.native_mechanics.tendon_force_N).toFixed(3)} N<br>Active / passive fiber force ${Number(s.native_mechanics.active_fiber_force_N).toFixed(3)} / ${Number(s.native_mechanics.passive_fiber_force_N).toFixed(3)} N</dd><dt>Joint moment arms</dt><dd>${Object.entries(
            s.native_mechanics.moment_arms_m || {},
          )
            .map(([k, v]) => `${esc(k)}: ${Number(v).toFixed(5)} m`)
            .join(
              "<br>",
            )}</dd><dt>Mechanical validity</dt><dd>Source default pose; external force balance ${s.native_mechanics.external_force_balance_solved ? "solved" : "not solved"}. No subject calibration.</dd>`
        : ""
    }${source.sha256 ? `<dt>SHA-256</dt><dd class="hash">${esc(source.sha256)}</dd>` : ""}</dl>${/^https?:\/\//.test(source.url || "") ? `<a class="source-link" href="${esc(source.url)}" target="_blank" rel="noopener noreferrer">Open original source ↗</a>` : ""}`;
  $("structures")
    .querySelectorAll("button")
    .forEach((b) => b.classList.toggle("selected", b.dataset.id === s.id));
  if (s.kind === "scalar_mesh") updateFrame();
}
function updateDisplay() {
  const opacity = Number($("opacity").value),
    clip = Number($("clip").value),
    size = modelBounds.getSize(new THREE.Vector3());
  const plane = new THREE.Plane(
    new THREE.Vector3(-1, 0, 0),
    modelBounds.min.x + (size.x * clip) / 100,
  );
  objects.forEach((obj) =>
    obj.traverse((o) => {
      if (o.material) {
        const alpha = obj.userData.structure?.kind === "vectors" ? 1 : opacity;
        o.material.opacity = alpha;
        o.material.transparent = alpha < 1 || !!o.material.isLineBasicMaterial;
        o.material.depthWrite = alpha > 0.5;
        o.material.clippingPlanes = clip < 100 ? [plane] : [];
      }
    }),
  );
  $("opacity-value").textContent = `${Math.round(opacity * 100)}%`;
  $("clip-value").textContent = clip === 100 ? "Off" : `${clip}%`;
}
function setupFrames() {
  const scalar = [...objects.values()].find(
    (o) => o.visible && o.userData.scalarFields,
  );
  $("flow-field").innerHTML = scalar
    ? scalar.userData.scalarFields
        .map(
          (f) =>
            `<option value="${esc(f.id)}">${esc(f.label)} · ${esc(f.unit)}</option>`,
        )
        .join("")
    : '<option value="velocity">Velocity</option><option value="pressure">Pressure</option>';
  $("flow-field").setAttribute("aria-label", "Spatial field");
  flowFrames = Math.max(
    0,
    ...[...objects.values()]
      .filter((o) => o.visible)
      .map((o) => o.userData.frames.length),
  );
  $("play").disabled = !flowFrames;
  $("time").disabled = !flowFrames;
  $("time").max = Math.max(0, flowFrames - 1);
  $("time").value = 0;
  playing = false;
  $("play").textContent = "▶";
  updateFrame();
}
function updateFrame() {
  $("flow-legend").hidden = !flowFrames;
  const index = Number($("time").value);
  let label = "";
  objects.forEach((o) => {
    const f = o.userData.frames[index];
    if (!f) return;
    label =
      o.userData.times[index] !== undefined
        ? `${Number(o.userData.times[index]).toFixed(3)} s`
        : `Frame ${index + 1} / ${flowFrames}`;
    if (f.positions && o.geometry) {
      o.geometry.setAttribute(
        "position",
        new THREE.Float32BufferAttribute(f.positions.flat(), 3),
      );
      o.geometry.computeBoundingSphere();
    }
    if (o.userData.scalarFields) {
      const field = o.userData.scalarFields.find(
        (x) => x.id === $("flow-field").value,
      );
      if (!field) return;
      const values = f.values[field.id];
      const normalized = scalarCoordinates(
        o.userData.cellIds,
        values,
        field.range,
      );
      const colors = o.geometry.getAttribute("color"),
        color = new THREE.Color();
      normalized.forEach((v, i) => {
        color.setHSL(0.66 * (1 - v), 0.75, 0.55);
        colors.setXYZ(i, color.r, color.g, color.b);
      });
      colors.needsUpdate = true;
      $("flow-legend").textContent =
        `${field.label} · blue ${field.range[0].toPrecision(4)} → red ${field.range[1].toPrecision(4)} ${field.unit} · fixed source range`;
      if (selected?.id === o.userData.structure.id) {
        let panel = $("cell-values");
        if (!panel) {
          panel = document.createElement("div");
          panel.id = "cell-values";
          $("details").append(panel);
        }
        const cell = o.userData.selectedCell;
        panel.innerHTML =
          `<h3>Native tissue fields</h3><p>212 planar solver cells · ${esc(label)}. Transmembrane voltage (V) differs from extracellular wound field (V/m).</p>` +
          (cell === undefined
            ? "<p>Click a cell to inspect its native values.</p>"
            : `<strong>Cell ${cell}</strong><dl>${o.userData.scalarFields.map((x) => `<dt>${esc(x.label)}</dt><dd>${Number(f.values[x.id][cell]).toPrecision(6)} ${esc(x.unit)}</dd>`).join("")}</dl>`) +
          '<p class="muted">Generic computational tissue; no human skin registration or wound-healing validation.</p>';
      }
    }
    const v = f.velocities || f.vectors;
    if (v && o.userData.vectorData) {
      const data = o.userData.vectorData,
        values = v.flat();
      const pressure = f.pressure || o.userData.pressure,
        [low, high] = o.userData.pressureRange;
      $("flow-legend").textContent =
        $("flow-field").value === "pressure"
          ? `Pressure · blue ${low.toPrecision(3)} → red ${high.toPrecision(3)} · source units unconfirmed`
          : "Velocity vectors · source units unconfirmed · archived CFD";
      o.children.forEach((arrow) => {
        const i = arrow.userData.nodeIndex,
          vec = new THREE.Vector3(...values.slice(i * 3, i * 3 + 3));
        if ($("flow-field").value === "pressure" && pressure.length) {
          arrow.setColor(
            new THREE.Color().setHSL(
              0.66 * (1 - (pressure[i] - low) / (high - low || 1)),
              0.7,
              0.6,
            ),
          );
        } else arrow.setColor(o.userData.structure.color || "#73dfd4");
        if (vec.length() > 0) {
          arrow.setDirection(vec.clone().normalize());
          arrow.setLength(vec.length() * data.scale);
        }
      });
    }
  });
  $("time-value").textContent = label || "No flow frames";
  $("time-value").title =
    modelId === "betse-tissue"
      ? "Native solver time; fixed field ranges across recorded frames"
      : $("flow-field").value === "pressure"
        ? "Pressure: blue = archive minimum, red = archive maximum; fixed scale across frames; source units unconfirmed"
        : "Velocity direction and magnitude from archived solver states";
}
function spectralRun() {
  return (
    temporal.runs?.find((r) => r.id === $("spectral-run").value) ||
    temporal.runs?.[0]
  );
}
function updateVariables() {
  $("spectral-controls").hidden = !spectral;
  $("trajectory-controls").hidden = spectral;
  const vars = spectral
    ? (spectralRun()?.variables || []).map((v) => ({
        id: v.id,
        label: v.label || v.id,
      }))
    : Object.keys(phys.values || {})
        .sort((a, b) =>
          a.startsWith("ArterialPressure")
            ? -1
            : b.startsWith("ArterialPressure")
              ? 1
              : 0,
        )
        .map((id) => ({ id, label: id }));
  const old = $("variable").value;
  $("variable").innerHTML = vars
    .map((v) => `<option value="${esc(v.id)}">${esc(v.label)}</option>`)
    .join("");
  if (vars.some((v) => v.id === old)) $("variable").value = old;
  drawChart();
}
function drawChart() {
  const id = $("variable").value,
    v = spectralRun()?.variables.find((x) => x.id === id),
    series = spectralSeries(
      spectralRun(),
      id,
      $("spectral-mode").value,
      Number($("sigma").value) || 0,
    );
  const x = spectral ? series.x : phys.time_axis || phys.time_s,
    y = spectral ? series.y : phys.values?.[id];
  $("signal-source").textContent = spectral
    ? spectralRun()?.source_kind || "No spectral evidence"
    : phys.metadata?.source_kind || "Native simulation";
  $("chart-note").textContent = spectral
    ? $("spectral-mode").value === "laplace"
      ? "Finite-horizon transform magnitude |L(σ + iω)| · mean removed · not an infinite-time transfer function."
      : spectralRun()?.sample_rate_hz <= 1
        ? "1 Hz source output · Nyquist 0.5 Hz; instantaneous arterial pressure and lung volume excluded as aliased. 900 s windows resolve 1/900 Hz. Drift does not establish a physiological cycle."
        : spectralRun()?.limitations?.[0] ||
          temporal.limitations?.[0] ||
          "Finite observation horizon. Power spectra do not establish causality."
    : activeRun === "reproductive"
      ? `Schlosser–Selgrade source model · ${phys.channel_types?.[id]?.replaceAll("_", " ") || "source channel"} · 0–10 day example; E2/P4/inhibin are prescribed, not a generated full cycle.`
      : activeRun.startsWith("csf:")
        ? "Ursino–Lodi source simulation; rounded published initial state. " +
          (activeRun === "csf:native_map_driven"
            ? "One-way native MAP forcing; unmatched subjects and initial state; no ICP feedback."
            : "Separate literature model, not measured or patient-calibrated ICP.")
        : activeRun.startsWith("thermal:")
          ? "JOS-3 source supine model · 85 thermal nodes / 17 regions · one hour after neutral standing initialization. " +
            (activeRun === "thermal:lying_default"
              ? "Source default environment."
              : "Published whole-body bedding resistance applied uniformly; regional contact measurements unavailable.") +
            " Separate source subject; no clinical bed-rest calibration."
          : "Upstream native simulation · model output, not a human recording.";
  if (!x?.length || !y?.length) {
    $("chart").innerHTML =
      '<p class="empty">No recorded data available for this view.</p>';
    return;
  }
  const valid = y.filter(Number.isFinite),
    min = Math.min(...valid),
    max = Math.max(...valid);
  $("chart").innerHTML =
    `<div class="chart-values"><strong>${Number(y.at(-1)).toPrecision(5)}</strong><span>${esc(spectral ? series.unit || "Source units" : phys.units?.[id] || id.match(/\(([^)]+)\)/)?.[1] || "Source units")}</span></div><svg viewBox="0 0 600 120" preserveAspectRatio="none" aria-label="${esc(id)} ${spectral ? "power spectrum" : "trajectory"}" role="img"><defs><linearGradient id="fade" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#8ec9bd" stop-opacity=".2"/><stop offset="1" stop-color="#8ec9bd" stop-opacity="0"/></linearGradient></defs><path d="M0 0H600 M0 40H600 M0 80H600 M0 120H600" stroke="#ffffff0b" fill="none"/><path d="${chartPath(x, y)}" fill="none" stroke="#8ec9bd" stroke-width="2" vector-effect="non-scaling-stroke"/></svg><div class="chart-axis"><span>${Number(x[0]).toPrecision(3)} ${spectral ? "Hz" : phys.time_unit || "s"}</span><span>range ${min.toPrecision(4)} — ${max.toPrecision(4)}</span><span>${Number(x.at(-1)).toPrecision(3)} ${spectral ? "Hz" : phys.time_unit || "s"}</span></div>`;
}
async function pollRuns() {
  try {
    const data = await api("/api/scenarios");
    const runs = data.runs || [];
    const selectedVariant = $("engine-variant").value;
    $("engine-variant").innerHTML = (
      data.engine_variants || [
        { id: "upstream", label: "Original upstream", available: true },
      ]
    )
      .map(
        (v) =>
          `<option value="${esc(v.id)}" ${v.available ? "" : "disabled"}>${esc(v.label)}</option>`,
      )
      .join("");
    $("engine-variant").value =
      selectedVariant || data.default_engine_variant || "upstream";
    const selectedPatient = $("patient").value;
    $("patient").innerHTML = (data.patients || ["StandardMale"])
      .map(
        (p) =>
          `<option value="${esc(p)}">${esc(p)}${p === "StandardFemale" && $("engine-variant").value === "upstream" ? " · upstream initialization fails" : ""}</option>`,
      )
      .join("");
    if ((data.patients || []).includes(selectedPatient))
      $("patient").value = selectedPatient;
    else if ((data.patients || []).includes("StandardMale"))
      $("patient").value = "StandardMale";
    const latest = runs.at(-1);
    $("trajectory-run").innerHTML =
      '<option value="baseline">Recorded resting baseline</option>' +
      (reproductive
        ? '<option value="reproductive">Gonadotropin regulation · source model · 0–10 days</option>'
        : "") +
      (csfAvailable
        ? '<option value="csf:baseline">CSF · source baseline</option><option value="csf:native_map_driven">CSF · one-way native MAP · unmatched initial state</option><option value="csf:hypotension">CSF · source hypotension ramp</option>'
        : "") +
      (thermalIndex?.runs || [])
        .map(
          (r) =>
            `<option value="thermal:${esc(r.id)}">${esc(r.label)}</option>`,
        )
        .join("") +
      runs
        .map(
          (r) =>
            `<option value="${esc(r.id)}" ${["completed", "complete", "succeeded"].includes(r.status) ? "" : "disabled"}>${esc(r.id)} · ${esc(r.status)}</option>`,
        )
        .join("");
    $("trajectory-run").value = activeRun;
    const chosen = runs.find((r) => r.id === activeRun);
    $("run-status").textContent = latest
      ? `${latest.id} · ${latest.status}${latest.error ? ` · ${latest.error}` : ""}`
      : data.available === false
        ? "Native backend unavailable"
        : "Native engine ready · bounded local execution";
    if (
      chosen &&
      ["completed", "complete", "succeeded"].includes(chosen.status) &&
      phys._run !== chosen.id
    ) {
      phys = await api(`/api/physiology?run=${encodeURIComponent(chosen.id)}`);
      phys._run = chosen.id;
      updateVariables();
    }
  } catch (e) {
    $("run-status").textContent = e.message;
  }
}
$("trajectory-run").onchange = async () => {
  activeRun = $("trajectory-run").value;
  try {
    phys =
      activeRun === "reproductive"
        ? trajectoryFromChannels(reproductive)
        : activeRun.startsWith("thermal:")
          ? trajectoryFromChannels(
              await api(
                `/api/thermal?run=${encodeURIComponent(activeRun.slice(8))}`,
              ),
            )
          : activeRun.startsWith("csf:")
            ? trajectoryFromChannels(
                await api(
                  `/api/csf?run=${encodeURIComponent(activeRun.slice(4))}`,
                ),
              )
            : await api(
                activeRun === "baseline"
                  ? "/api/physiology"
                  : `/api/physiology?run=${encodeURIComponent(activeRun)}`,
              );
    phys._run = activeRun;
    updateVariables();
  } catch (e) {
    $("chart").innerHTML = `<p class="empty">${esc(e.message)}</p>`;
  }
};
$("engine-variant").onchange = () => {
  for (const option of $("patient").options)
    option.textContent =
      option.value +
      (option.value === "StandardFemale" &&
      $("engine-variant").value === "upstream"
        ? " · upstream initialization fails"
        : "");
  $("variant-note").textContent =
    $("engine-variant").value === "upstream"
      ? "Original upstream build; female initialization failure is preserved in run history."
      : "Explicit source bounds correction; heatflux variant also corrects diagnostic output. These corrections do not establish clinical calibration.";
};
$("scenario").onchange = () => {
  $("protocol-note").textContent =
    $("scenario").value === "exercise"
      ? "Native exercise intensity 0.1, then stop and observe recovery."
      : $("scenario").value === "hemorrhage_saline"
        ? "Right-leg hemorrhage 10 mL/min, then stop; saline 20 mL/min, then stop and observe. This does not guarantee return to baseline."
        : "Record the upstream resting state.";
};
$("scenario-form").onsubmit = async (e) => {
  e.preventDefault();
  $("run").disabled = true;
  try {
    const config = scenarioInput(
      $("scenario").value,
      $("duration").value,
      $("patient").value,
      {
        ambient_temperature_c: $("ambient").value,
        clothing_clo: $("clothing").value,
        engine_variant: $("engine-variant").value,
      },
    );
    const run = await api("/api/scenarios", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(config),
    });
    activeRun = run.id;
    $("run-status").textContent = `${run.id} · ${run.status}`;
    await pollRuns();
  } catch (e) {
    $("run-status").textContent = e.message;
  } finally {
    $("run").disabled = false;
  }
};
$("posture").onclick = () => {
  if (!group) return;
  group.rotation.x = group.rotation.x === 0 ? -Math.PI / 2 : 0;
  $("posture").textContent = group.rotation.x ? "Upright view" : "Supine view";
  modelBounds.setFromObject(group);
  resetCamera();
  $("frame-label").textContent =
    `${group.rotation.x ? "Supine" : "Upright"} display · rigid rotation only · native posture not established`;
  updateDisplay();
};
$("model").onchange = chooseModel;
$("search").oninput = () => manifest && refresh();
$("all").onclick = () => {
  if (!manifest) return;
  const boxes = [...$("systems").querySelectorAll("input")],
    show = boxes.some((b) => !b.checked);
  boxes.forEach((b) => {
    b.checked = show;
    show ? systems.add(b.value) : systems.delete(b.value);
  });
  $("all").textContent = show ? "Hide all" : "Show all";
  refresh();
};
$("opacity").oninput = updateDisplay;
$("clip").oninput = updateDisplay;
$("reset").onclick = () => resetCamera();
$("front").onclick = () => resetCamera();
$("side").onclick = () => resetCamera("side");
$("time").oninput = updateFrame;
$("flow-field").onchange = updateFrame;
$("play").onclick = () => {
  playing = !playing;
  $("play").textContent = playing ? "Ⅱ" : "▶";
};
$("variable").onchange = drawChart;
for (const [id, value] of [
  ["tab-phys", false],
  ["tab-spectral", true],
])
  $(id).onclick = () => {
    spectral = value;
    $("tab-phys").classList.toggle("active", !value);
    $("tab-spectral").classList.toggle("active", value);
    updateVariables();
  };
function updateSigma() {
  $("sigma").innerHTML = (spectralRun()?.laplace?.sigma_per_s || [])
    .map((x, i) => `<option value="${i}">σ = ${x} s⁻¹</option>`)
    .join("");
}
$("spectral-run").onchange = () => {
  updateSigma();
  updateVariables();
};
$("spectral-mode").onchange = () => {
  $("sigma").hidden = $("spectral-mode").value !== "laplace";
  drawChart();
};
$("sigma").onchange = drawChart;
async function loadEvidencePanels() {
  try {
    const data = await api("/api/calibration");
    $("calibration-details").innerHTML =
      (data.models || [])
        .map(
          (m) =>
            `<p><strong>Lateral skin electric field</strong><br>Human observations · ${esc(m.unit)}</p><dl><dt>Training / held-out observations</dt><dd>${esc(m.metrics?.train?.n)} / ${esc(m.metrics?.holdout?.n)}</dd><dt>Held-out prediction RMSE</dt><dd>${Number(m.metrics?.holdout?.rmse_V_per_m).toFixed(2)} V/m</dd><dt>Identifiable parameter rank</dt><dd>${esc(m.rank)} / ${esc(m.parameters_count)}</dd></dl><p>${esc(m.limitations?.[0])}</p><p>No external cohort validation or microscopic conductance calibration.</p>`,
        )
        .join("") || "<p>No measurement fits available.</p>";
  } catch {
    $("calibration-details").innerHTML =
      "<p>Calibration evidence unavailable.</p>";
  }

  const panels = await Promise.allSettled([
    api("/api/coverage"),
    api("/api/coupling"),
    api("/api/coupling/skin-lymph"),
    api("/api/vascular/audit"),
  ]);
  if (panels[0].status === "fulfilled") {
    const data = panels[0].value,
      summary = data.summary || {};
    $("coverage-details").innerHTML =
      `<p>${esc(summary.system_domains)} system domains · ${esc(summary.declared_components)} declared components · ${esc(summary.population_measured_components)} population measurement targets.</p><ul class="coverage-list">${(data.systems || []).map((s) => `<li><span>${esc(s.name)}</span><small>${s.measured_components?.length || 0} measured / ${s.components_count} components</small></li>`).join("")}</ul><p>Coverage and anatomical naming do not establish calibrated interactions or cross-specimen registration.</p>`;
  }
  if (panels[1].status === "fulfilled") {
    const data = panels[1].value,
      summary = data.summary || {},
      model = data.models?.find((m) => m.id === "native_skin_lymph_circuit");
    const duration = (t) =>
      t >= 3600 ? `${(t / 3600).toFixed(2)} h` : `${t.toFixed(2)} s`;
    $("coupling-details").innerHTML =
      `<p>${esc(summary.nodes)} native circuit nodes · ${esc(summary.paths)} paths · ${esc(summary.compartments)} compartments.</p>${
        model
          ? `<dl><dt>Frozen skin–lymph response times</dt><dd>${[
              ...(model.time_constants_s || []),
            ]
              .sort((a, b) => a - b)
              .map(duration)
              .join(
                " · ",
              )}</dd><dt>Zero modes</dt><dd>${esc(model.zero_modes)} · long-horizon stability not established</dd></dl><p>${esc(model.limitations?.[0])}</p>`
          : ""
      }<p>${esc(summary.cross_source_coupling)}</p>`;
    if (panels[2].status === "fulfilled") {
      const step = panels[2].value.baseline_step;
      const residual = step?.balance?.max_abs_free_node_residual_m3_s;
      if (Number.isFinite(residual))
        $("coupling-details").innerHTML +=
          `<dl><dt>Maximum free-node flow residual</dt><dd>${residual.toExponential(2)} m³/s</dd><dt>Frozen-step checks</dt><dd>${step.gate_violations?.length || 0} gate violations · ${step.negative_volume_nodes?.length || 0} negative volume nodes</dd></dl><p>Numerical balance at the saved operating point is not independent physiological calibration.</p>`;
    }
  }
  if (panels[3].status === "fulfilled") {
    const data = panels[3].value,
      conservation = data.conservation || {},
      periodic = data.periodic_endpoint_check;
    $("vascular-audit").innerHTML =
      `<p>${data.clock?.audited_frames} archived states · ${esc(data.archived_job_status)}</p><dl><dt>Maximum global mass imbalance</dt><dd>${(100 * conservation.maximum_relative_global_mass_imbalance).toPrecision(3)}%</dd><dt>Same-phase velocity difference</dt><dd>${(100 * (periodic?.velocity_nodal_difference?.relative_rmse_to_expected_rms || 0)).toFixed(1)}%</dd><dt>Wall velocity</dt><dd>${esc(conservation.maximum_wall_nodal_speed)} source units</dd></dl><p>Good global balance does not establish periodic convergence. Units and archived pressure gauge remain unconfirmed.</p><a class="source-link" href="/api/vascular/audit" target="_blank" rel="noopener">Open complete numerical audit ↗</a>`;
  }
}
async function start() {
  try {
    manifest = await api("/api/manifest");
    $("total").textContent = manifest.structures.length.toLocaleString();
    $("model").innerHTML = manifest.models
      .map((m) => `<option value="${esc(m.id)}">${esc(m.name)}</option>`)
      .join("");
    chooseModel();
  } catch (e) {
    $("scene-status").textContent = `Anatomy unavailable · ${e.message}`;
    $("model").innerHTML = "<option>Source manifest unavailable</option>";
    $("model-note").textContent =
      "Start the local API and regenerate visualization assets. Reload to retry.";
  }
  const results = await Promise.allSettled([
    api("/api/physiology"),
    api("/api/temporal"),
    api("/api/evidence"),
    api("/api/reproductive"),
    api("/api/csf/index"),
    api("/api/thermal/index"),
  ]);
  if (results[0].status === "fulfilled") phys = results[0].value;
  if (results[1].status === "fulfilled") {
    temporal = results[1].value;
    $("spectral-run").innerHTML = (temporal.runs || [])
      .map(
        (r) =>
          `<option value="${esc(r.id)}">${esc(r.id)} · ${esc(r.source_kind)}</option>`,
      )
      .join("");
    updateSigma();
  }
  if (results[2].status === "fulfilled")
    $("evidence-status").textContent = results[2].value.whole_body_calibrated
      ? "See domain calibration evidence"
      : "Whole-body calibration incomplete";
  if (results[3].status === "fulfilled") reproductive = results[3].value;
  csfAvailable = results[4].status === "fulfilled";
  if (results[5].status === "fulfilled") thermalIndex = results[5].value;
  updateVariables();
  await pollRuns();
  setInterval(pollRuns, 5000);
}
start();
loadEvidencePanels();
