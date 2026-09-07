import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import {
  filterStructures,
  geometryArrays,
  bodyTransform,
  validBodyTrajectory,
  deformSkinVertices,
} from "./state.js";
import "./style.css";
import { attachedHairPositions } from "./hair_motion.js";
import { attachElasticHair, updateElasticHair } from "./hair-view.js";
import { WardrobeView } from "./clothing.js";
import { DomainView, ROLE_LABELS } from "./domains.js";
import { mountSceneInteraction } from "./scene-interaction.js";
import { mountLeftColumn } from "./left-panel.js";
import { mountPanes } from "./panes.js";
import { mountProvenance } from "./provenance.js";
import { mountGimbal } from "./gimbal.js";
import { mountSurround } from "./surround.js";
import { mountMicrovascularDetail } from "./microvascular-detail.js";

const MODEL_ID = "ihm-body";
const SPEEDS = [0.25, 0.5, 1, 2, 5, 10];
const esc = (s) =>
  String(s ?? "").replace(/[&<>"']/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);
const el = (tag, className, text) => {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
};
const $ = (id) => document.getElementById(id);

document.querySelector("#app").innerHTML = `
<div id="viewport">
  <button id="toggle-left" class="edge-toggle left" type="button" aria-label="Collapse controls" aria-expanded="true">◧</button>
  <button id="toggle-right" class="edge-toggle right" type="button" aria-label="Collapse panes" aria-expanded="true">◨</button>
  <aside id="left-column" aria-label="Body controls"></aside>
  <canvas id="gimbal" width="108" height="108" aria-label="Body plane views · drag the scene to turn it"></canvas>
  <aside id="pane-column" aria-label="Panes"></aside>
  <div id="transport">
    <button id="play" type="button" aria-label="Play recorded body" disabled>▶</button>
    <button id="speed" type="button" aria-expanded="false" aria-label="Playback speed">1x</button>
    <div id="speed-menu" role="group" aria-label="Playback rate" hidden></div>
  </div>
  <p id="scene-status" role="status">Loading anatomy…</p>
</div>`;

// ---------------------------------------------------------------- scene ----
const viewport = $("viewport");
let renderer, scene, camera, controls, group, webglError;
const modelBounds = new THREE.Box3();
try {
  renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(1);
  renderer.localClippingEnabled = true;
  renderer.setClearColor(0x0d1416, 1);
  renderer.domElement.id = "scene";
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
    const w = viewport.clientWidth, h = viewport.clientHeight;
    if (!w || !h) return;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }).observe(viewport);
} catch (e) {
  webglError = e.message;
  $("scene-status").textContent = `WebGL unavailable: ${e.message}`;
}

async function api(path, options) {
  const r = await fetch(path, options);
  if (!r.ok) {
    let message;
    try { message = await r.json(); } catch {}
    throw Error(message?.error || `Request failed (${r.status})`);
  }
  return r.json();
}

// ---------------------------------------------------------------- panes ----
const paneHost = $("pane-column");
const selectionBody = el("div");
selectionBody.id = "details";
selectionBody.append(el("p", "note", "Click a structure to identify it."));
const microvascularHost = el("div");
microvascularHost.id = "microvascular-monitor";
const sceneMonitor = el("div");
sceneMonitor.id = "scene-monitor";
const sceneControls = el("div");
sceneControls.id = "scene-controls";
const bodyInteraction = el("div");
bodyInteraction.append(sceneControls, sceneMonitor);
const liveHosts = {};
for (const id of ["live-body", "live-motor", "intake", "intake-mass",
  "skin-voltage", "temporal-spectrum", "live-signal-1", "live-signal-2", "live-signal-3"]) {
  const host = el("div");
  host.id = `${id}-monitor`;
  host.append(el("p", "note", "Start the body to inspect its computed state."));
  liveHosts[id] = host;
}
const playbackBody = el("div");
playbackBody.innerHTML =
  `<input id="time" type="range" min="0" max="0" value="0" aria-label="Recorded body frame" disabled><p id="time-value" class="note">No recorded frames</p>`;

const panes = mountPanes(paneHost, [
  { id: "selection", title: "Selection", content: selectionBody, open: true },
  { id: "live", title: "Live body", content: liveHosts["live-body"] },
  { id: "live-signal-1", title: "Live signal 1", content: liveHosts["live-signal-1"] },
  { id: "live-signal-2", title: "Live signal 2", content: liveHosts["live-signal-2"] },
  { id: "live-signal-3", title: "Live signal 3", content: liveHosts["live-signal-3"] },
  { id: "motor", title: "Motor & skin inputs", content: liveHosts["live-motor"] },
  { id: "intake", title: "Food & drink", content: liveHosts["intake"] },
  { id: "intake-mass", title: "Mechanical mass", content: liveHosts["intake-mass"] },
  { id: "skin-voltage", title: "Skin membrane voltage", content: liveHosts["skin-voltage"] },
  { id: "temporal-spectrum", title: "Live Laplace spectrum", content: liveHosts["temporal-spectrum"] },
  { id: "microvessels", title: "Local microvessels", content: microvascularHost },
  { id: "scene", title: "Body interaction", content: bodyInteraction },
  { id: "playback", title: "Playback", content: playbackBody, open: true },
  { id: "domain", title: "Material owners", content: (() => {
      const node = el("div");
      node.innerHTML = `<div id="domain-roles"></div><p id="domain-selected" class="note">Click a surface to identify the owner that carries it.</p>`;
      return node;
    })() },
]);
const provenanceView = mountProvenance(() => $("details"), api);
const microvascularDetail = mountMicrovascularDetail(microvascularHost);

// ------------------------------------------------------------- left column -
let manifest = null, structures = [];
let objects = new Map(), generation = 0, loadController;
let selected = null, clothingView = null, clothingRequest = 0, clothingCatalog = null;
let domainView = null, domainRequest = 0, materialization = "body";
let bodyTrajectory = null, bodyError = "Body trajectory unavailable";
let liveFrame = null, sceneInteraction = null;
let playing = false, speed = 1, lastTick = 0, lastRender = 0, lastLiveVisual = 0;
let hairDynamics = false;

const left = mountLeftColumn($("left-column"), {
  onMaterialization: chooseMaterialization,
  onLayers: () => { refresh(); },
  onMember: (id) => selectStructure(structures.find((s) => s.id === id)),
  onClothing: (ids) => applyGarments(ids),
  onEnvironment: (selection, objects) => {
    applyEnvironment(selection.environment);
    applySurround(selection, objects);
  },
  onSimulation: (option, on) => {
    if (option === "hair") { hairDynamics = on; updateFrame(); return; }
    sceneInteraction?.setOptions(left.options).then(syncRun).catch((e) => left.setRunNote(e.message));
  },
  onRun: () => {
    if (!sceneInteraction) return;
    if (!sceneInteraction.started) sceneInteraction.start().then(syncRun).catch((e) => left.setRunNote(e.message));
    else if (sceneInteraction.running) { sceneInteraction.pause(); syncRun(); }
    else sceneInteraction.start().then(syncRun).catch((e) => left.setRunNote(e.message));
    syncRun();
  },
});
// The declared surround: every environment and scene carries its sky, its
// ground treatment and its enclosing geometry under `world`, and names the
// parts to cut away so an interior can be seen into. All of it is scenery — the
// catalogue's own world.physics says no collider exists for any of it — so it
// is drawn into a group of its own that anatomy picking never reaches.
let sceneCatalog = null, surround = null, surroundRequest = 0, framedScene = null;
function applySurround(selection = {}, objects = []) {
  if (!surround || !sceneCatalog) return;
  const request = ++surroundRequest;
  const environment = sceneCatalog.environments?.find((e) => e.id === selection.environment);
  const scene = sceneCatalog.scenes?.find((s) => s.id === selection.scene);
  const entry = scene || environment;
  surround.apply(entry, environment, objects, sceneCatalog.objects || [])
    .then((result) => {
      if (request !== surroundRequest) return;
      if (result.note) left.setRunNote(result.note);
      // Choosing a world is a request to see it. The camera goes to the view the
      // catalogue itself renders its scene tiles from, so what appears is what
      // the tile promised; leaving the scene leaves the camera where it is.
      if (scene && scene.id !== framedScene) { framedScene = scene.id; frameScene(); }
      if (!scene) framedScene = null;
    })
    .catch((error) => { if (request === surroundRequest) left.setRunNote("Surround unavailable: " + error.message); });
}
// camera.view_direction_canonical points from the scene centre toward the
// catalogue's camera; image_right_canonical is its right vector, which is how
// this reading was checked rather than guessed.
function frameScene() {
  const view = sceneCatalog?.camera;
  if (!camera || !view?.view_direction_canonical) return;
  const half = view.half_extent_m?.scene_tiles;
  if (!Number.isFinite(half)) return;
  const centre = new THREE.Vector3(...(view.centre_m || [0, 0, 0]));
  const direction = new THREE.Vector3(...view.view_direction_canonical).normalize();
  const distance = half / Math.tan((camera.fov * Math.PI) / 360) * 1.15;
  controls.target.copy(centre);
  camera.position.copy(centre).addScaledVector(direction, distance);
  camera.up.set(...(view.image_up_canonical || [0, 1, 0]));
  camera.near = distance / 1000;
  camera.far = Math.max(camera.far, distance * 40);
  camera.updateProjectionMatrix();
  controls.update();
  gimbal?.select(null);
}

// POST /api/embodied/sessions accepts one environment and nothing else, so the
// component slots and scene objects are held as declared state until the
// runtime accepts them.
function applyEnvironment(id) {
  if (!id) return;
  sceneInteraction?.setEnvironment(id).then(syncRun).catch((e) => left.setRunNote(e.message));
}
function syncRun() {
  if (!sceneInteraction) { left.setRun("Start body", true); return; }
  left.setRun(!sceneInteraction.started ? "Start body" : sceneInteraction.running ? "Pause body" : "Resume body");
  left.lockDynamics(sceneInteraction.started);
}

// -------------------------------------------------------------- rendering --
function visibleRows() {
  return filterStructures(structures, MODEL_ID, left.systems, "")
    .filter((s) => !left.hidden.has(s.id));
}
function refresh() {
  if (!manifest || materialization !== "body") return;
  loadVisible(visibleRows());
}
async function loadVisible(rows) {
  if (!group) { $("scene-status").textContent = `WebGL unavailable: ${webglError || "no graphics context"}`; return; }
  const current = ++generation;
  loadController?.abort();
  loadController = new AbortController();
  const signal = loadController.signal, wanted = new Set(rows.map((x) => x.id));
  objects.forEach((object, id) => { object.visible = wanted.has(id); });
  let queue = rows.filter((x) => !objects.has(x.id)), done = rows.length - queue.length, failures = 0;
  const status = () => {
    if (current !== generation || domainView) return;
    $("scene-status").textContent = done === rows.length
      ? "" : `Loading anatomy · ${done} / ${rows.length}${failures ? ` · ${failures} unavailable` : ""}`;
  };
  status();
  await Promise.all(Array.from({ length: 2 }, async () => {
    while (queue.length && current === generation) {
      const s = queue.shift();
      try {
        const g = await api(s.geometry_url || `/api/geometry/${encodeURIComponent(s.id)}`, { signal });
        if (current !== generation) return;
        const object = createGeometry(s, g);
        if (object && group) {
          objects.set(s.id, object);
          group.add(object);
          if (!modelBounds.isEmpty()) modelBounds.expandByObject(object);
        }
        done++; status();
      } catch (e) {
        if (e.name === "AbortError") return;
        failures++; status();
      }
    }
  }));
  if (current === generation) { status(); updateFrame(); }
}
function createGeometry(s, g) {
  const { positions, indices } = geometryArrays(g), geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  if (indices) geometry.setIndex(new THREE.BufferAttribute(indices, 1));
  const color = s.color || "#9eb7b7";
  let object;
  if (s.kind === "lines") {
    object = new (g.topology === "polyline" ? THREE.Line : THREE.LineSegments)(
      geometry, new THREE.LineBasicMaterial({ color, transparent: true }));
  } else {
    geometry.computeVertexNormals();
    object = new THREE.Mesh(geometry, new THREE.MeshStandardMaterial({
      color, roughness: 0.65, metalness: 0.05, side: THREE.DoubleSide,
      transparent: true, opacity: s.system === "integumentary" ? 0.22 : 1 }));
  }
  object.userData.structure = s;
  attachElasticHair(object, g);
  if (g.attachment?.kind === "MaterialPoint" && g.attachment.reference_triangles_m) {
    object.userData.hairAttachment = g.attachment;
    object.userData.hairReference = positions.slice();
  }
  return object;
}
function selectStructure(s) {
  if (!s) return;
  selected = s;
  microvascularDetail.updateSelection(s);
  objects.forEach((object, id) => object.traverse((x) => {
    if (x.material?.emissive) x.material.emissive.set(id === s.id ? "#1d514f" : "#000000");
  }));
  const source = s.source || {};
  $("details").innerHTML =
    `<div class="eyebrow">${esc(s.system)}</div><h3>${esc(s.name)}</h3>` +
    `<dl><dt>Source</dt><dd>${esc(source.label || "IHM · one generic human")}</dd>` +
    `<dt>Evidence</dt><dd>${esc(source.status || s.calibration_status || "Source geometry; physiological registration unverified")}</dd></dl>`;
  // Seam for the uniform per-structure provenance record.
  provenanceView.show(s);
  panes.show("selection");
}

// --------------------------------------------------------------- camera ----
function resetCamera(plane = "coronal") {
  gimbal?.select(plane);
  if (domainView) { domainView.resetCamera(plane === "sagittal" ? "side" : "front"); return; }
  if (!camera || modelBounds.isEmpty()) return;
  const center = modelBounds.getCenter(new THREE.Vector3()),
    size = modelBounds.getSize(new THREE.Vector3()),
    distance = Math.max(size.y, size.x, size.z) * 1.85;
  controls.target.copy(center);
  const offset = plane === "sagittal" ? new THREE.Vector3(distance, 0, 0)
    : plane === "transverse" ? new THREE.Vector3(0, distance, 0.001)
      : new THREE.Vector3(0, 0, distance);
  camera.position.copy(center).add(offset);
  camera.up.set(0, plane === "transverse" ? 0 : 1, plane === "transverse" ? -1 : 0);
  camera.near = distance / 1000;
  camera.far = distance * 20;
  camera.updateProjectionMatrix();
  controls.update();
}
if (scene) surround = mountSurround(scene, { api });
let gimbal = null;
if (camera) {
  try {
    gimbal = mountGimbal($("gimbal"), { camera, onSelect: (plane) => resetCamera(plane) });
    gimbal.select("coronal");
  } catch (error) {
    $("gimbal").hidden = true;
    gimbal = null;
  }
}

// ------------------------------------------------------------- transport ---
$("speed-menu").replaceChildren(...SPEEDS.map((value) => {
  const button = el("button", null, `${value}x`);
  button.type = "button";
  button.onclick = () => {
    speed = value;
    $("speed").textContent = `${value}x`;
    $("speed-menu").hidden = true;
    $("speed").setAttribute("aria-expanded", "false");
  };
  return button;
}));
$("speed").onclick = () => {
  const open = $("speed-menu").hidden;
  $("speed-menu").hidden = !open;
  $("speed").setAttribute("aria-expanded", String(open));
};
$("play").onclick = () => {
  if (liveFrame) return;
  playing = !playing;
  lastTick = performance.now();
  $("play").textContent = playing ? "❚❚" : "▶";
};
$("time").oninput = updateFrame;

// -------------------------------------------------------------- toggles ----
function syncInsets() {
  for (const [side, name] of [["left-column", "--left-inset"], ["pane-column", "--right-inset"]]) {
    const column = $(side);
    viewport.style.setProperty(name,
      column.classList.contains("collapsed") ? "12px" : `${column.offsetWidth + 22}px`);
  }
}
for (const [id, side] of [["toggle-left", "left-column"], ["toggle-right", "pane-column"]])
  $(id).onclick = () => {
    const target = $(side), collapsed = !target.classList.contains("collapsed");
    target.classList.toggle("collapsed", collapsed);
    $(id).setAttribute("aria-expanded", String(!collapsed));
    syncInsets();
  };
syncInsets();
new ResizeObserver(syncInsets).observe(document.documentElement);

// ------------------------------------------------------------- clothing ----
// The wardrobe is 33 registered, cloth-simulated garments served whole by the
// API. They are already in this body's own frame, so there is no fitting step
// here: the view binds them to the skin entity and the body's transform carries
// them. What is worn is the catalog's exclusivity model, resolved by the tiles.
function loadClothing() {
  const request = ++clothingRequest;
  clothingView?.dispose(); clothingView = null;
  if (materialization !== "body" || !group) return;
  const skin = structures.find((s) => s.id === "body-bp3d-FJ2810");
  if (!skin || request !== clothingRequest) return;
  clothingView = new WardrobeView(group, { fetchGeometry: (url) => api(url) });
  clothingView.bind(skin);
  if (clothingCatalog) clothingView.setCatalog(clothingCatalog);
  applyGarments(left.garments);
}

function applyGarments(ids) {
  if (!clothingView) return;
  clothingView.setActive(ids).then((failures) => {
    if (failures.length) left.setClothingNote("Garment geometry unavailable · " + failures.join("; "));
    updateFrame();
  });
}

// The wardrobe declares no outfit, so the workbench opens with one garment in
// each everyday slot — the first the catalog lists for it — and leaves hats,
// gloves, outerwear and the rest for the reader to add. Exclusivity is still
// the catalog's: these are ordinary selections a click can undo.
const OPENING_SLOTS = ["underwear_bottom", "underwear_top", "torso_base", "legs", "feet_outer"];
function openingOutfit(catalog) {
  const worn = [];
  for (const slot of OPENING_SLOTS) {
    const garment = catalog.garments.find((g) => (g.slots || []).includes(slot));
    if (garment && !worn.some((id) => {
      const held = catalog.garments.find((g) => g.id === id);
      return (held?.slots || []).some((s) => (garment.slots || []).includes(s));
    })) worn.push(garment.id);
  }
  return worn;
}

// -------------------------------------------------------- materialization --
async function chooseMaterialization(value) {
  const request = ++domainRequest;
  materialization = value;
  playing = false; $("play").textContent = "▶";
  try { if (sceneInteraction) await sceneInteraction.reset(); }
  catch (error) { left.setMaterializationNote("Body cleanup must finish first: " + error.message); return; }
  if (request !== domainRequest) return;
  syncRun();
  domainView?.close(); domainView = null;
  panes.hide("domain");
  left.setWholeBody(value === "body");
  if (value === "body") {
    $("scene-status").textContent = "";
    left.setMaterializationNote("The whole assembled body. Conforming tetrahedral domains are listed above.");
    refresh(); loadClothing(); setupFrames();
    return;
  }
  left.setMaterializationNote("Loading the conforming tetrahedral domain…");
  $("scene-status").textContent = "Loading the conforming tetrahedral domain…";
  try {
    const data = await api("/api/body/experiments/" + encodeURIComponent(value));
    if (request !== domainRequest) return;
    if (!renderer) throw Error("A conforming domain needs WebGL");
    domainView = new DomainView({ scene, camera, controls, bodyGroup: group });
    domainView.open(data);
    renderDomainRoles();
    panes.show("domain");
    $("scene-status").textContent = "";
    left.setMaterializationNote(data.summary);
    setupFrames();
  } catch (error) {
    if (request !== domainRequest) return;
    domainView?.close(); domainView = null;
    left.setMaterializationNote(error.message);
    $("scene-status").textContent = "";
  }
}
function renderDomainRoles() {
  if (!domainView) return;
  $("domain-roles").innerHTML = domainView.roles().map((r) =>
    `<label class="check-row"><input type="checkbox" data-role="${esc(r.role)}"${r.visible ? " checked" : ""}><i style="background:${esc(domainView.meshes.find((m) => m.name === r.role).material.color.getStyle())}"></i><span>${esc(ROLE_LABELS[r.role] || r.role)}</span><small>${r.structures.toLocaleString()}</small></label>`).join("");
  for (const box of $("domain-roles").querySelectorAll("input"))
    box.onchange = () => { domainView.setRoleVisible(box.dataset.role, box.checked); renderDomainRoles(); };
}

// ---------------------------------------------------------------- frames ---
function setupFrames() {
  const frames = domainView ? 0 : bodyTrajectory?.frames.length || 0;
  const live = !!liveFrame;
  $("play").disabled = live || frames < 2;
  $("time").disabled = live || !frames;
  $("time").max = Math.max(0, frames - 1);
  if (!live) $("time").value = 0;
  playing = false;
  $("play").textContent = "▶";
  updateFrame();
}
function updateFrame() {
  if (domainView) {
    $("time-value").textContent = "Static conforming volume · nothing is integrated in time here.";
    return;
  }
  if (liveFrame) {
    const reference = Object.fromEntries(Object.entries(liveFrame.entities).map(([id, state]) =>
      [id, state.centroid_m.map((x, i) => x - (state.translation_m?.[i] || 0))]));
    applyBodyFrame(liveFrame, { centroids_m: reference });
    $("time-value").textContent = `Live body · ${liveFrame.time_s.toFixed(3)} s`;
    return;
  }
  const frame = bodyTrajectory?.frames[Number($("time").value)];
  applyBodyFrame(frame, bodyTrajectory);
  $("time-value").textContent = frame ? `${Number(frame.time_s).toFixed(3)} s` : bodyError;
}
function applyBodyFrame(frame, trajectory) {
  const skinField = frame?.respiration?.skin_field;
  const skinIds = new Set(skinField?.entity_ids || []);
  objects.forEach((object, id) => {
    if (updateElasticHair(object, {
      frame, referenceCentroids: trajectory?.centroids_m,
      recordKey: liveFrame ? "live" : "recorded", enabled: hairDynamics,
      visible: object.visible && !document.hidden,
      gravity_m_s2: frame?.environment?.gravity || [0, -9.81, 0],
    })) return;
    const positions = object.geometry?.getAttribute("position");
    const hair = object.userData.hairAttachment;
    if (positions && hair) {
      attachedHairPositions(object.userData.hairReference, hair,
        skinIds.has(hair.skin_entity_id) ? skinField : null, positions.array);
      positions.needsUpdate = true;
      object.geometry.computeVertexNormals();
      object.geometry.computeBoundingSphere();
    }
    if (positions && (skinIds.has(id) || object.userData.skinDeformationActive)) {
      object.userData.skinReference ||= positions.array.slice();
      deformSkinVertices(object.userData.skinReference, skinIds.has(id) ? skinField : null, positions.array);
      positions.needsUpdate = true;
      object.geometry.computeVertexNormals();
      object.geometry.computeBoundingSphere();
      object.userData.skinDeformationActive = skinIds.has(id);
    }
    const motionId = hair?.skin_entity_id || id;
    object.matrixAutoUpdate = false;
    object.matrix.set(...bodyTransform(frame?.entities?.[motionId], trajectory?.centroids_m?.[motionId]));
    object.matrixWorldNeedsUpdate = true;
  });
  clothingView?.update(frame, trajectory?.centroids_m);
}

// ------------------------------------------------------------------ loop ---
if (renderer) {
  renderer.setAnimationLoop((now) => {
    if (document.hidden || now - lastRender < 1000 / 30) return;
    lastRender = now;
    controls.update();
    gimbal?.update();
    surround?.follow(camera);
    sceneInteraction?.update(now);
    const frames = bodyTrajectory?.frames;
    if (playing && frames?.length && !liveFrame && !domainView) {
      const index = Number($("time").value);
      const interval = 1000 * ((frames[index + 1]?.time_s ?? frames[index].time_s + 0.1) - frames[index].time_s) / speed;
      if (now - lastTick >= interval) {
        $("time").value = (index + 1) % frames.length;
        updateFrame();
        lastTick = now;
      }
    }
    renderer.render(scene, camera);
  });
  const ray = new THREE.Raycaster();
  let pointer;
  renderer.domElement.addEventListener("pointerdown", (e) => { pointer = [e.clientX, e.clientY]; });
  renderer.domElement.addEventListener("pointerup", (e) => {
    if (!pointer || Math.hypot(e.clientX - pointer[0], e.clientY - pointer[1]) > 5) return;
    const r = renderer.domElement.getBoundingClientRect();
    ray.setFromCamera(new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1,
      (-(e.clientY - r.top) / r.height) * 2 + 1), camera);
    if (domainView) {
      const found = domainView.pick(ray);
      $("domain-selected").textContent = found
        ? `${found.structure.name} · ${ROLE_LABELS[found.role] || found.role} · ${found.structure.tets.toLocaleString()} owned tetrahedra`
        : "No visible owner under the cursor.";
      panes.show("domain");
      return;
    }
    const hit = ray.intersectObjects([...objects.values()].filter((x) => x.visible), true)[0];
    if (!hit) return;
    let object = hit.object;
    while (!object.userData.structure && object.parent) object = object.parent;
    selectStructure(object.userData.structure);
  });
}
document.addEventListener("visibilitychange", () => { lastRender = 0; lastTick = performance.now(); });

// ------------------------------------------------------------------ boot ---
function mountBody() {
  if (!renderer) return;
  sceneInteraction = mountSceneInteraction({
    scene, camera, renderer, controls, group,
    getObjects: () => objects,
    onSelect: selectStructure,
    mount: sceneControls, monitor: sceneMonitor,
    onStatus: (text) => left.setRunNote(text),
    onPauseReplay: () => { playing = false; $("play").textContent = "▶"; $("play").disabled = true; $("time").disabled = true; },
    onFrame: (frame) => {
      const first = frame && !liveFrame;
      liveFrame = frame;
      if (frame) {
        if (first) for (const id of ["live", "live-signal-1", "motor", "scene"]) panes.show(id);
        const now = performance.now();
        if (!document.hidden && (first || now - lastLiveVisual >= 200)) { lastLiveVisual = now; updateFrame(); }
      } else setupFrames();
      syncRun();
    },
  });
  syncRun();
}

async function start() {
  try {
    manifest = await api("/api/manifest");
  } catch (e) {
    $("scene-status").textContent = `Anatomy unavailable · ${e.message}`;
    left.setMaterializationNote("Start the local API and rebuild visualization assets.");
    return;
  }
  structures = manifest.structures.filter((s) => s.model_id === MODEL_ID);
  const model = manifest.models.find((m) => m.id === MODEL_ID);
  if (model?.bounds) {
    const b = model.bounds;
    modelBounds.set(new THREE.Vector3(...(b.min || b[0])), new THREE.Vector3(...(b.max || b[1])));
  }
  const defaults = new Set(structures.filter((s) => s.default_visible).map((s) => s.system));
  defaults.add("integumentary");
  defaults.add("hair");
  left.setStructures(structures, defaults);
  left.setWholeBody(true);
  left.setMaterializations([{ value: "body", label: "Whole body" }], "body");
  left.setMaterializationNote("The whole assembled body.");
  resetCamera();
  refresh();
  await loadClothing();
  mountBody();
  api("/api/clothing")
    .then((catalog) => {
      clothingCatalog = catalog;
      clothingView?.setCatalog(catalog);
      left.setGarments(catalog.garments, { initial: openingOutfit(catalog) });
      left.setClothingNote("");
      applyGarments(left.garments);
    })
    .catch((error) => {
      left.setGarments([], {});
      left.setClothingNote("Wardrobe unavailable: " + error.message);
    });
  // The scene catalog owns identity, label, slot, thumbnail, requirements and
  // per-slot defaults, for environments and their components alike. Nothing
  // about exclusivity or dependency is decided here.
  api("/api/scene/catalog")
    .then((catalog) => {
      sceneCatalog = catalog;
      const tiles = catalog.tiles || (catalog.environments || []).map((e) => ({
        id: e.id, label: e.label, slot: e.slot || "environment", kind: e.kind || "environment",
        thumbnail_url: e.thumbnail_url ?? null, requires: e.requires || [],
      }));
      const slots = catalog.slots || [{ id: "environment", exclusive: true, required: true, default: "bed" }];
      left.setEnvironments(tiles, {
        slots,
        initial: slots.map((s) => s.default).filter(Boolean),
      });
      applyEnvironment(left.environment);
      applySurround(left.environmentSelection, left.sceneObjects.map((i) => i.id));
    })
    .catch((error) => {
      left.setEnvironments([], {});
      left.setRunNote("Scene catalog unavailable: " + error.message);
    });
  api("/api/body/experiments/conforming-domains").then((index) => {
    left.setMaterializations(
      [{ value: "body", label: "Whole body" },
        ...index.domains.map((d) => ({ value: "conforming-domain-" + d.id, label: d.title }))],
      materialization);
  }).catch((error) => left.setMaterializationNote("Conforming domains unavailable: " + error.message));
  try {
    const data = await api("/api/body/trajectory?view=display");
    if (!validBodyTrajectory(data)) throw Error("No valid computed body frames available");
    bodyTrajectory = data;
    bodyError = "";
  } catch (error) { bodyError = error.message; }
  setupFrames();
}
start();
