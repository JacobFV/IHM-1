import {configureSurfaceAssets} from './surface-assets.js';
import {updateSegmentSurface} from "./surface-binding.js";
import * as THREE from "three";
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
import { mountLeftColumn, readMaterialization, writeMaterialization } from "./left-panel.js";
import { labToHex, resolvePalette } from "./palette.js";
import { mountPanes } from "./panes.js";
import { mountProvenance } from "./provenance.js";
import { mountGimbal } from "./gimbal.js";
import { mountCameraOrbit } from "./camera-orbit.js";
import { mountSurround } from "./surround.js";
import { mountMicrovascularDetail } from "./microvascular-detail.js";
import { transportState, recordingCatalog, recordingHeader, RECORDING } from "./transport.js";
import { mountPromptBar } from "./prompt-bar.js";
import { mountRing } from "./ring.js";
import { tissueMaterial, tissueUVs, applyOpacity } from "./tissue-materials.js";

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
  <div id="recording-header" role="status" hidden>
    <span class="eyebrow">Recording</span>
    <b id="recording-name"></b>
    <span id="recording-detail"></span>
    <button id="recording-close" type="button">Close recording</button>
  </div>
  <div id="ring-host"></div>
  <div id="transport" role="group" aria-label="Transport">
    <span id="transport-label">Simulation</span>
    <button id="play" type="button" aria-label="Run or hold the simulation">▶</button>
    <input id="time" type="range" min="0" max="0" value="0" aria-label="Scrub the open recording" hidden disabled>
    <button id="speed" type="button" aria-expanded="false" aria-label="Playback speed" hidden>1x</button>
    <div id="speed-menu" role="group" aria-label="Playback rate" hidden></div>
    <button id="open-recording" type="button" aria-expanded="false" aria-label="Open a recording">Open recording…</button>
    <div id="recording-menu" role="group" aria-label="Recordings" hidden></div>
  </div>
  <div id="live-readout" role="status" hidden><i></i><span id="live-readout-text"></span></div>
  <div id="prompt-host"></div>
  <div id="orientation" aria-hidden="true"><span data-edge="top"></span><span data-edge="bottom"></span><span data-edge="left"></span><span data-edge="right"></span></div>
  <div id="scale-bar" aria-hidden="true" hidden><i></i><span></span></div>
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
  renderer.setClearColor(0x0a0e11, 1);
  renderer.domElement.id = "scene";
  viewport.prepend(renderer.domElement);
  scene = new THREE.Scene();
  camera = new THREE.PerspectiveCamera(35, 1, 0.02, 200);
  // A pole-free orbit rather than a turntable. The canonical frame is x left, y
  // superior, z anterior; a turntable orbits about one fixed axis, which was y,
  // so the arc over the head -- around the coronal plane -- ran into a stop and
  // the transverse view sat on a pole where no drag moved anything. See
  // camera-orbit.js.
  controls = mountCameraOrbit(camera, renderer.domElement);
  controls.enableDamping = true;
  scene.add(new THREE.HemisphereLight(0xd7f0ee, 0x26333d, 2));
  const light = new THREE.DirectionalLight(0xfff1da, 3);
  light.position.set(3, 4, 5);
  light.castShadow = true;
  light.shadow.mapSize.set(2048, 2048);
  Object.assign(light.shadow.camera, {left:-4, right:4, top:4, bottom:-4, near:.1, far:20});
  light.shadow.normalBias = .012;
  // Front faces write the shadow map now that solids are single sided, so the
  // depth it stores is the lit surface itself and a constant bias is enough to
  // keep a surface from shadowing itself.
  light.shadow.bias = -0.0004;
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = .85;
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
// A recording is a document you open, not a mode you switch into. This pane
// says which one is open and what it is; the scrubbing itself belongs to the
// transport, because while a recording is open the transport IS the
// recording's — the way a video player's scrubber is obviously not a camera.
const playbackBody = el("div");
playbackBody.innerHTML =
  `<p id="time-value" class="note">No recording is open. The simulation is live.</p>` +
  `<p id="recording-source" class="note"></p>` +
  `<p class="note">Opening a recording holds the simulation and hands the transport to the recording. Close it to hand the transport back.</p>`;

const ringPanelHost = el("div");
ringPanelHost.id = "ring-pane";
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
  { id: "signal-graph", title: "Signal graph", content: ringPanelHost, open: true },
  { id: "playback", title: "Recording", content: playbackBody },
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
// The implicit model keeps every structure; a materialization draws a subset
// of it. Which subset, what it is called and what it costs are all the
// server's, derived from the manifest and the colocation audits. Nothing
// here knows a structure count or an identifier.
let fidelityIndex = null, fidelityOmits = new Set(), conformingDomains = [];
let bodyTrajectory = null, bodyError = "Body trajectory unavailable";
let liveFrame = null, sceneInteraction = null;
// A recording is open or it is not; that one fact decides what the transport is.
let recording = null, recordingData = null, recordingError = "";
let playing = false, speed = 1, lastTick = 0, lastRender = 0, lastLiveVisual = 0;
let ring = null, brainGraph = null, brainStateTimer = 0, lastLeaders = 0;
let hairDynamics = false;

const left = mountLeftColumn($("left-column"), {
  onMaterialization: chooseMaterialization,
  onLayers: () => { refresh(); },
  onOpacity: (system, value) => applySystemOpacity(system, value),
  onPalette: (id, skinTone) => applyPalette(id, skinTone),
  onMember: (id) => selectStructure(structures.find((s) => s.id === id)),
  onClothing: (ids) => applyGarments(ids),
  onEnvironment: (selection, objects) => {
    applyEnvironment(selection.environment, selection, objects);
    applySurround(selection, objects);
  },
  onSimulation: (option, on) => {
    if (option === "hair") { hairDynamics = on; updateFrame(); return; }
    sceneInteraction?.setOptions(left.options).then(syncRun).catch((e) => left.setRunNote(e.message));
  },
});
// Catalogue geometry and server-owned environment state share canonical coordinates.
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
// Both clipping planes, set together from how far the camera is standing off,
// and refreshed every frame because the wheel is not bounded. A near plane a
// thousandth of the viewing distance spends almost the whole depth buffer on
// the first few centimetres, and what is left cannot separate a mattress from
// the sheet on it; a far plane that only ever grew kept the ratio at whatever
// the widest scene had ever needed. A fiftieth and fortyfold hold a whole room
// at a ratio of two thousand, which coincident surfaces survive.
function clipTo(distance) {
  if (!Number.isFinite(distance) || distance <= 0) return;
  camera.near = distance / 50;
  camera.far = distance * 40;
}
// The orbit dollies without limit, so the standoff that set the planes at the
// last snap is not the standoff now. Recomputing costs one subtraction and only
// touches the projection when the reader has actually moved.
function clipToView() {
  if (!camera || !controls) return;
  const distance = camera.position.distanceTo(controls.target);
  const near = distance / 50;
  if (Math.abs(near - camera.near) < camera.near * 0.05) return;
  clipTo(distance);
  camera.updateProjectionMatrix();
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
  const up = new THREE.Vector3(...(view.image_up_canonical || [0, 1, 0]));
  const gravityFrame = view.gravity_frames?.[left.environment];
  if (gravityFrame?.rotation) {
    const r = gravityFrame.rotation;
    const transform = new THREE.Matrix3().set(...r[0], ...r[1], ...r[2]).transpose();
    direction.applyMatrix3(transform); up.applyMatrix3(transform);
    if (gravityFrame.centre) centre.fromArray(gravityFrame.centre).applyMatrix3(transform);
  }
  const distance = half / Math.tan((camera.fov * Math.PI) / 360) * 1.15;
  controls.target.copy(centre);
  camera.position.copy(centre).addScaledVector(direction, distance);
  camera.up.copy(up);
  // A catalogue scene carries its own up, and in a gravity frame that is not the
  // superior axis. It is what a settled view should keep vertical from here on.
  controls.upReference.copy(up).normalize();
  clipTo(distance);
  camera.updateProjectionMatrix();
  controls.update();
  gimbal?.select(null);
}

function applyEnvironment(id, selection = left.environmentSelection, objects = left.sceneObjects.map(i => i.id)) {
  if (!id) return;
  const configuration = {objects};
  for (const slot of ["scene", "bed_support_model", "mattress_material", "ambient_thermal"])
    if (selection[slot]) configuration[slot] = selection[slot];
  sceneInteraction?.setEnvironment(id, configuration).then(syncRun).catch((e) => left.setRunNote(e.message));
}
function syncRun() {
  if (!sceneInteraction) { left.setRun("Body view not mounted"); syncTransport(); return; }
  left.setRun(!sceneInteraction.started ? "Body stopped · press play in the transport"
    : sceneInteraction.running ? "Body running" : "Body held");
  left.lockDynamics(sceneInteraction.started);
  syncTransport();
}

// -------------------------------------------------------------- rendering --
function visibleRows() {
  return filterStructures(structures, MODEL_ID, left.systems, "")
    .filter((s) => !left.hidden.has(s.id) && !fidelityOmits.has(s.id));
}
function refresh() {
  if (!manifest || !wholeBody()) return;
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
  // Whatever the reader last set this system to, which starts at the value the
  // system was always drawn at -- skin translucent, everything else solid.
  const opacity = left.opacity(s.system);
  let object;
  if (s.kind === "lines") {
    object = new (g.topology === "polyline" ? THREE.Line : THREE.LineSegments)(
      geometry, applyOpacity(new THREE.LineBasicMaterial({ color }), opacity));
  } else {
    geometry.computeVertexNormals();
    // Tissue is drawn as tissue: the palette still decides the colour, and the
    // material adds the grain and relief that colour alone cannot carry. See
    // tissue-materials.js for why the maps are grey.
    const material = tissueMaterial(s, { color, opacity });
    tissueUVs(geometry, material.userData.tissueFamily);
    object = new THREE.Mesh(geometry, material);
  }
  object.castShadow = object.isMesh && s.system === "integumentary";
  object.receiveShadow = object.isMesh;
  object.userData.structure = s;
  attachElasticHair(object, g);
  if (g.attachment?.kind === "MaterialPoint" && g.attachment.reference_triangles_m) {
    object.userData.hairAttachment = g.attachment;
    object.userData.hairReference = positions.slice();
  }
  return object;
}
// Opacity is not visibility: this changes how a system is painted and leaves
// every checkbox, and every object's `visible`, exactly where it was.
function applySystemOpacity(system, value) {
  objects.forEach((object) => {
    if (object.userData.structure?.system !== system) return;
    object.traverse((x) => {
      if (!x.material) return;
      for (const material of Array.isArray(x.material) ? x.material : [x.material]) {
        const wasTransparent = material.transparent;
        applyOpacity(material, value);
        // Crossing the boundary between opaque and translucent changes which
        // pass the material is compiled for, and three only notices if told.
        if (material.transparent !== wasTransparent) material.needsUpdate = true;
      }
    });
  });
}

// Choosing a palette repaints what is already in the scene and rewrites the
// colour on the structure record, so geometry that loads later comes up in the
// same palette. No geometry is refetched: the only thing that changes is the
// material colour, and the record the renderer reads it from.
const paletteCache = new Map();
let paletteIndex = null, paletteBase = null, paletteRequest = 0;
async function applyPalette(id, skinTone) {
  if (!structures.length) return;
  // The colour every structure started at, so a palette that does not name a
  // structure leaves it exactly where it was rather than blanking it.
  paletteBase ||= new Map(structures.map((s) => [s.id, s.color]));
  const request = ++paletteRequest;
  let palette;
  try {
    if (!paletteCache.has(id))
      paletteCache.set(id, await api(`/api/body/palettes/${encodeURIComponent(id)}`));
    palette = paletteCache.get(id);
  } catch (error) {
    paletteCache.delete(id);
    if (request === paletteRequest) left.setPaletteApplied({ error: `Palette unavailable: ${error.message}` });
    return;
  }
  if (request !== paletteRequest) return;
  const tone = paletteIndex?.skin_tone_options?.[skinTone]?.colour;
  const toneHex = labToHex(tone?.lab, tone?.illuminant_observer);
  const { colours, missing, skinToneApplies } = resolvePalette({ palette, structures, base: paletteBase, toneHex });
  for (const s of structures) {
    const hex = colours.get(s.id);
    if (!hex) continue;
    s.color = hex;
    const object = objects.get(s.id);
    if (!object?.material) continue;
    for (const material of Array.isArray(object.material) ? object.material : [object.material])
      material.color?.set(hex);
  }
  left.setPaletteApplied({ skinToneApplies, fallbacks: missing.length });
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
  // The transverse view used to be nudged a millimetre off the axis because a
  // turntable's azimuth degenerates on its pole; the orbit no longer has one, so
  // the snap is exactly on axis. A snap also restores the superior axis as the
  // upright reference, which is how a view carried over the head gets back up.
  const offset = plane === "sagittal" ? new THREE.Vector3(distance, 0, 0)
    : plane === "transverse" ? new THREE.Vector3(0, distance, 0)
      : new THREE.Vector3(0, 0, distance);
  camera.position.copy(center).add(offset);
  camera.up.set(0, plane === "transverse" ? 0 : 1, plane === "transverse" ? -1 : 0);
  controls.upReference.set(0, 1, 0);
  clipTo(distance);
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

// A pose seam. How far a drag turns the camera is a number in degrees, and no
// screenshot reports it, so the browser specs and app/test/orbit-isotropy.mjs
// read the pose here and put the camera at a known start before dragging.
// Nothing in the app calls this; it only reads and writes the camera.
globalThis.__ihmCamera = {
  pose() {
    if (!camera || !controls) return null;
    const direction = camera.position.clone().sub(controls.target).normalize();
    return {
      position: camera.position.toArray(),
      target: controls.target.toArray(),
      up: camera.up.toArray(),
      quaternion: camera.quaternion.toArray(),
      direction: direction.toArray(),
    };
  },
  setPose({ position, target, up }) {
    if (!camera || !controls) return null;
    if (target) controls.target.fromArray(target);
    if (position) camera.position.fromArray(position);
    if (up) camera.up.fromArray(up).normalize();
    camera.lookAt(controls.target);
    controls.update();
    return this.pose();
  },
  snap(plane) { resetCamera(plane); return this.pose(); },
  // Damping spreads one drag over many frames. This runs the same update the
  // render loop runs, so a measurement reads a finished camera rather than one
  // still coasting.
  settle(frames = 600) {
    for (let i = 0; i < frames; i++) controls?.update();
    return this.pose();
  },
};

// What the renderer is actually painting a system at, read back from the
// materials rather than from the panel that asked for it. The browser spec uses
// it to tell a slider that moved from a body that changed.
globalThis.__ihmLayers = {
  opacity(system) {
    const values = [];
    objects.forEach((object) => {
      if (object.userData.structure?.system !== system) return;
      object.traverse((x) => { if (x.material?.opacity !== undefined) values.push(x.material.opacity); });
    });
    return values.length ? values.reduce((a, b) => a + b, 0) / values.length : null;
  },
  // The same readback for colour: what one structure's material is actually
  // painted, not what the palette record said it should be.
  colour(id) {
    const material = objects.get(id)?.material;
    const first = Array.isArray(material) ? material[0] : material;
    return first?.color ? `#${first.color.getHexString()}` : null;
  },
};

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
// ONE control. What it runs or holds is whatever the transport is attached to:
// the simulation, always, unless a recording is open. There is no shared button
// and no mode toggle -- opening a recording is opening a document.
$("play").onclick = () => {
  const model = transportModel();
  if (model.disabled) return;
  if (model.owner === RECORDING) {
    playing = !playing;
    lastTick = performance.now();
    syncTransport();
    return;
  }
  if (!sceneInteraction) { left.setRunNote("The body view is not mounted."); return; }
  if (sceneInteraction.running) { sceneInteraction.pause(); syncRun(); return; }
  sceneInteraction.start().then(syncRun).catch((e) => left.setRunNote(e.message));
  syncRun();
};
$("time").oninput = () => { updateFrame(); syncTransport(); };

// -- opening a recording -----------------------------------------------------
$("open-recording").onclick = async () => {
  const menu = $("recording-menu");
  const open = menu.hidden;
  menu.hidden = !open;
  $("open-recording").setAttribute("aria-expanded", String(open));
  if (!open) return;
  menu.replaceChildren(el("p", "note", "Reading what has been computed…"));
  const scenarios = await api("/api/scenarios").catch(() => null);
  const entries = recordingCatalog({
    canonical: { available: !bodyError, detail: bodyTrajectory
      ? `${bodyTrajectory.frames.length} computed frames` : "The materialized native run for this body." },
    scenarios,
  });
  menu.replaceChildren();
  if (!entries.length) {
    menu.append(el("p", "note", bodyError
      ? `No recording is available: ${bodyError}`
      : "No computed recording is available yet."));
    return;
  }
  for (const entry of entries) {
    const button = el("button", null, entry.label);
    button.type = "button";
    button.append(el("small", null, entry.detail));
    button.onclick = () => { menu.hidden = true; $("open-recording").setAttribute("aria-expanded", "false"); openRecording(entry); };
    menu.append(button);
  }
};
document.addEventListener("pointerdown", (event) => {
  const menu = $("recording-menu");
  if (!menu.hidden && !$("transport").contains(event.target)) {
    menu.hidden = true;
    $("open-recording").setAttribute("aria-expanded", "false");
  }
});
$("recording-close").onclick = () => closeRecording();

async function openRecording(entry) {
  // Opening a recording holds the simulation. It does not close the session:
  // the body is still there, held at its last accepted state, and closing the
  // recording hands the transport straight back to it.
  if (sceneInteraction?.running) sceneInteraction.pause();
  syncRun();
  playing = false;
  recording = entry;
  recordingError = "";
  $("recording-header").hidden = false;
  $("recording-name").textContent = entry.label;
  $("recording-detail").textContent = "loading…";
  try {
    if (entry.id === "canonical" && bodyTrajectory) recordingData = bodyTrajectory;
    else {
      const data = await api(entry.url);
      if (!validBodyTrajectory(data)) throw Error("No valid computed body frames in this recording");
      recordingData = data;
    }
  } catch (error) {
    recordingData = null;
    recordingError = error.message;
  }
  $("recording-source").textContent = recordingError
    ? recordingError : `${entry.label} · ${entry.detail}`;
  setupFrames();
}
function closeRecording() {
  recording = null;
  recordingData = null;
  recordingError = "";
  playing = false;
  $("recording-header").hidden = true;
  $("recording-source").textContent = "";
  setupFrames();
}

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
// The prompt panel grows upward from the bottom edge, and the transport sits
// above it. What "above" is has to be measured, not guessed at, or a five-entry
// retrieval covers the play button.
new ResizeObserver(() => {
  viewport.style.setProperty("--prompt-inset", `${$("prompt-host").offsetHeight || 48}px`);
  ring?.resize();
}).observe($("prompt-host"));

// ------------------------------------------------------------- clothing ----
// The wardrobe is 33 registered, cloth-simulated garments served whole by the
// API. They are already in this body's own frame, so there is no fitting step
// here: the view binds them to the same native segment supports as the skin
// when supplied by the live owner. What is worn is the catalog's exclusivity model, resolved by the tiles.
function loadClothing() {
  const request = ++clothingRequest;
  clothingView?.dispose(); clothingView = null;
  if (!wholeBody() || !group) return;
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
// each everyday slot — the first the catalog lists for it — and leaves bra tops,
// hats, gloves, outerwear and the rest for the reader to add. Exclusivity is still
// the catalog's: these are ordinary selections a click can undo.
const OPENING_SLOTS = ["underwear_bottom", "torso_base", "legs", "feet_outer"];
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
// A fidelity tier is a whole-body materialization; a conforming domain is not.
const fidelityTier = (value = materialization) =>
  fidelityIndex?.tiers?.find((t) => t.value === value) || null;
const wholeBody = (value = materialization) => value === "body" || !!fidelityTier(value);
// What a tier costs, in the numbers it was measured at through this app. A
// materialization the reader is warned about vaguely is one they cannot judge.
function fidelityPrice(tier) {
  const cost = tier?.cost || {};
  return Number.isFinite(cost.time_to_ready_s) && Number.isFinite(cost.js_heap_mb)
    ? `~${cost.time_to_ready_s} s to load, ~${cost.js_heap_mb} MB` : "";
}
// Every tier is offered, the complete one included, each labelled with the count
// it draws and with what that cost when it was measured.
function materializationOptions() {
  const tiers = fidelityIndex?.tiers?.map((t) => {
    const price = fidelityPrice(t);
    return {
      value: t.value,
      label: `${t.label} · ${t.structures.toLocaleString()} structures${price ? " · " + price : ""}`,
    };
  }) || [{ value: "body", label: "Whole body" }];
  return [...tiers, ...conformingDomains];
}
// Switching tiers is a change of membership, not a reload: what the scene
// already holds stays loaded and only what the tier newly asks for is fetched.
function applyFidelity(value) {
  const tier = fidelityTier(value);
  fidelityOmits = new Set(tier?.omits || []);
  // The closed menu truncates, so the tier that is actually in force says here
  // what it draws and what it cost.
  const price = fidelityPrice(tier);
  left.setMaterializationNote(tier
    ? `${tier.note}${price ? ` Measured here at ${price}.` : ""}`
    : "The whole assembled body. Conforming tetrahedral domains are listed above.");
}
async function chooseMaterialization(value) {
  const request = ++domainRequest;
  materialization = value;
  // A conforming tetrahedral volume is static and a recording is a trajectory;
  // they cannot both own the view, and the volume is what was just asked for.
  if (recording && !wholeBody(value)) closeRecording();
  playing = false;
  try { if (sceneInteraction) await sceneInteraction.reset(); }
  catch (error) { left.setMaterializationNote("Body cleanup must finish first: " + error.message); return; }
  if (request !== domainRequest) return;
  syncRun();
  domainView?.close(); domainView = null;
  panes.hide("domain");
  left.setWholeBody(wholeBody(value));
  if (wholeBody(value)) {
    $("scene-status").textContent = "";
    applyFidelity(value);
    if (fidelityTier(value)) writeMaterialization(value);
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

// The bottom centre is ONE control with one meaning: run or hold whatever the
// transport is attached to. It is attached to the simulation, because the
// simulation is the app; it is attached to a recording only while one is open,
// and then the header above says which one. Nothing about this is a tooltip.
function transportModel() {
  return transportState({
    open: !!recording,
    recording,
    frames: recordingData?.frames.length || 0,
    index: Number($("time").value) || 0,
    playing,
    recordingError,
    started: !!sceneInteraction?.started,
    running: !!sceneInteraction?.running,
    busy: !!sceneInteraction && sceneInteraction.started && !sceneInteraction.state,
    staticVolume: !!domainView,
  });
}
function syncTransport() {
  const model = transportModel();
  const transport = $("transport");
  transport.dataset.owner = model.owner;
  transport.dataset.empty = String(model.disabled);
  $("transport-label").textContent = model.label;
  $("play").textContent = model.symbol;
  $("play").disabled = model.disabled;
  $("play").title = model.title;
  $("play").setAttribute("aria-label", model.title);
  $("time").hidden = !model.scrub;
  $("time").disabled = !model.scrub || model.frames < 2;
  $("time").max = Math.max(0, model.frames - 1);
  $("speed").hidden = !model.speed;
  if (!model.speed) { $("speed-menu").hidden = true; $("speed").setAttribute("aria-expanded", "false"); }
  $("open-recording").textContent = recording ? "Open another…" : "Open recording…";

  const header = recordingHeader(recording, {
    frames: model.frames, index: model.index,
    time_s: recordingData?.frames[model.index]?.time_s,
    error: recordingError,
  });
  $("recording-header").hidden = !header;
  // The header owns the top strip while it is up, so what sits under it has to
  // know how tall it actually is rather than assume.
  viewport.style.setProperty("--header-inset",
    header ? `${($("recording-header").offsetHeight || 32) + 8}px` : "0px");
  if (header) {
    $("recording-name").textContent = header.name;
    $("recording-detail").textContent = header.detail;
    $("recording-header").dataset.failed = String(header.failed);
  }
  // The live badge is the body's own clock and stays readable while a recording
  // is open, because the body is still there, held.
  const live = !!liveFrame;
  $("live-readout").hidden = !live;
  if (live) {
    const running = !!sceneInteraction?.running;
    $("live-readout").dataset.running = String(running);
    $("live-readout-text").textContent =
      `${running ? "Live" : "Held"} · ${liveFrame.time_s.toFixed(2)} s computed`;
  }
}

// ---------------------------------------------------------------- frames ---
function setupFrames() {
  const frames = recording ? recordingData?.frames.length || 0 : 0;
  $("time").max = Math.max(0, frames - 1);
  if (!recording) $("time").value = 0;
  playing = false;
  syncTransport();
  updateFrame();
}
function updateFrame() {
  if (domainView) {
    $("time-value").textContent = "Static conforming volume · nothing is integrated in time here.";
    return;
  }
  // An open recording is what you are looking at. Nothing else can claim the
  // view while it is open, which is the whole point of opening it.
  if (recording) {
    const frame = recordingData?.frames[Number($("time").value)];
    applyBodyFrame(frame, recordingData);
    $("time-value").textContent = frame
      ? `${recording.label} · recorded · ${Number(frame.time_s).toFixed(3)} s`
      : recordingError || "This recording holds no readable frames.";
    return;
  }
  if (liveFrame) {
    const reference = Object.fromEntries(Object.entries(liveFrame.entities).map(([id, state]) =>
      [id, state.centroid_m.map((x, i) => x - (state.translation_m?.[i] || 0))]));
    applyBodyFrame(liveFrame, { centroids_m: reference });
    $("time-value").textContent = `No recording is open. Live body · ${liveFrame.time_s.toFixed(3)} s.`;
    return;
  }
  applyBodyFrame(null, null);
  $("time-value").textContent = "No recording is open. The simulation is live.";
}
function applyBodyFrame(frame, trajectory) {
  const skinField = frame?.respiration?.skin_field;
  const skinIds = new Set(skinField?.entity_ids || []);
  objects.forEach((object, id) => {
    if (updateSegmentSurface(object,frame,id)) return;
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

// ------------------------------------------------------------ annotation ---
// What an imaging console prints on every frame: the anatomical direction at
// each edge of the view, and a scale bar true at the orbit target's distance.
// The canonical frame is x left, y superior, z anterior, so the letter at an
// edge names the direction of the body that lies that way on the screen.
const ORIENTATION = [["L", "R"], ["S", "I"], ["A", "P"]];
const SCALE_STEPS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5];
const orientationEdges = Object.fromEntries(
  [...document.querySelectorAll("#orientation span")].map((node) => [node.dataset.edge, node]));
const orientationScratch = { axis: new THREE.Vector3(), turn: new THREE.Quaternion(), text: "" };
function orientationLetter(axis) {
  const parts = [axis.x, axis.y, axis.z];
  const index = parts.map(Math.abs).indexOf(Math.max(...parts.map(Math.abs)));
  return Math.abs(parts[index]) < 0.35 ? "" : ORIENTATION[index][parts[index] >= 0 ? 0 : 1];
}
function annotate() {
  const { axis, turn } = orientationScratch;
  group.getWorldQuaternion(turn).invert().multiply(camera.quaternion);
  const letters = {
    right: orientationLetter(axis.set(1, 0, 0).applyQuaternion(turn)),
    left: orientationLetter(axis.set(-1, 0, 0).applyQuaternion(turn)),
    top: orientationLetter(axis.set(0, 1, 0).applyQuaternion(turn)),
    bottom: orientationLetter(axis.set(0, -1, 0).applyQuaternion(turn)),
  };
  const text = letters.left + letters.right + letters.top + letters.bottom;
  if (text !== orientationScratch.text) {
    orientationScratch.text = text;
    for (const edge in letters) orientationEdges[edge].textContent = letters[edge];
  }
  const bar = $("scale-bar");
  const height = viewport.clientHeight;
  const distance = camera.position.distanceTo(controls.target);
  if (!height || !(distance > 0) || modelBounds.isEmpty()) { bar.hidden = true; return; }
  const metresPerPixel = 2 * distance * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2) / height;
  const metres = [...SCALE_STEPS].reverse().find((step) => step / metresPerPixel <= 150) || SCALE_STEPS[0];
  const pixels = Math.round(metres / metresPerPixel);
  if (bar.dataset.pixels !== String(pixels)) {
    bar.dataset.pixels = String(pixels);
    bar.firstElementChild.style.width = `${pixels}px`;
    bar.lastElementChild.textContent = metres >= 1 ? `${metres} m` : metres >= 0.01 ? `${Math.round(metres * 100)} cm` : `${Math.round(metres * 1000)} mm`;
  }
  bar.hidden = pixels < 12;
}

// ------------------------------------------------------------------ loop ---
if (renderer) {
  renderer.setAnimationLoop((now) => {
    if (document.hidden || now - lastRender < 1000 / 30) return;
    lastRender = now;
    controls.update();
    clipToView();
    gimbal?.update();
    annotate();
    surround?.follow(camera);
    sceneInteraction?.update(now);
    // Only an OPEN recording advances here. The simulation advances in
    // sceneInteraction.update, because it is a simulation and not a film.
    const frames = recording ? recordingData?.frames : null;
    if (playing && frames?.length && !domainView) {
      const index = Number($("time").value);
      const interval = 1000 * ((frames[index + 1]?.time_s ?? frames[index].time_s + 0.1) - frames[index].time_s) / speed;
      if (now - lastTick >= interval) {
        $("time").value = (index + 1) % frames.length;
        updateFrame();
        syncTransport();
        lastTick = now;
      }
    }
    // The leaders track the orbit, but re-measuring 26 label boxes every frame
    // is layout work for a line that moves a pixel; ten times a second is
    // indistinguishable and costs a fortieth of it.
    if (now - lastLeaders >= 100) { lastLeaders = now; ring?.follow(); }
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

// ------------------------------------------------------- ring & prompt bar -
// The ring is the ablation made continuous: everything exchanging signal with
// the highlighted system, each arrow drawn at the width its contribution was
// measured at. See ring.js for why there is no flow animation in it.
const projectScratch = new THREE.Vector3();
function projectCanonical(point) {
  if (!camera || !group || !Array.isArray(point)) return null;
  projectScratch.set(point[0], point[1], point[2]);
  group.localToWorld(projectScratch);
  projectScratch.project(camera);
  const w = viewport.clientWidth, h = viewport.clientHeight;
  if (!w || !h) return null;
  return [(projectScratch.x * 0.5 + 0.5) * w, (-projectScratch.y * 0.5 + 0.5) * h];
}
function mountRingAndPrompt() {
  ring = mountRing($("ring-host"), { project: projectCanonical, panelHost: ringPanelHost });
  mountPromptBar($("prompt-host"), {
    post: async (body) => {
      const response = await fetch("/api/brain/prompt", {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      });
      let payload = null;
      try { payload = await response.json(); } catch {}
      if (!response.ok) {
        const error = Error(payload?.error || `Prompt request failed (${response.status})`);
        error.httpStatus = response.status;
        throw error;
      }
      return payload;
    },
  });
  new ResizeObserver(() => ring?.resize()).observe(viewport);
  api("/api/brain/graph")
    .then((graph) => { brainGraph = graph; ring.setGraph(graph); refreshBrainState(); })
    .catch((error) => {
      $("ring-host").append(el("p", "note ring-basis",
        `Signal graph unavailable: ${error.message}. Nothing is drawn in its place.`));
    });
}
// Per-edge magnitudes. With no body running this comes back saying so, and the
// ring shows the measured values labelled measured — it never fills the gap.
async function refreshBrainState() {
  if (!brainGraph || !ring) return;
  const session = sceneInteraction?.session;
  try {
    ring.setState(await api(`/api/brain/state${session ? `?session=${encodeURIComponent(session)}` : ""}`));
  } catch { /* the graph's own measured values stand */ }
}

// ------------------------------------------------------------------ boot ---
function mountBody() {
  configureSurfaceAssets({onReady:()=>updateFrame(),onError:error=>left.setRunNote(error.message)});
  if (!renderer) return;
  mountRingAndPrompt();
  sceneInteraction = mountSceneInteraction({
    scene, camera, renderer, controls, group,
    getObjects: () => objects, getEnvironmentObjects: () => surround?.interactiveObjects||[],
    onSelect: selectStructure,
    mount: sceneControls, monitor: sceneMonitor,
    onStatus: (text) => left.setRunNote(text),
    // A body coming up no longer has to disable anything: the transport it
    // shares with nothing is already the simulation's.
    onPauseReplay: () => { syncTransport(); },
    onFrame: (frame) => {
      const first = frame && !liveFrame;
      liveFrame = frame;
      viewport.dataset.live = String(!!frame);
      surround?.update(frame?.environment_state);
      if (frame) {
        if (first) for (const id of ["live", "live-signal-1", "motor", "scene"]) panes.show(id);
        const now = performance.now();
        if (!document.hidden && (first || now - lastLiveVisual >= 200)) { lastLiveVisual = now; updateFrame(); syncTransport(); }
        if (now - brainStateTimer >= 1000) { brainStateTimer = now; refreshBrainState(); }
      } else { setupFrames(); refreshBrainState(); }
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
  // The palette index names the default and the app keeps it: didactic is what
  // the manifest already paints, so restoring it is a no-op rather than a
  // recolour. Another palette is an explicit choice, and it survives a reload.
  api("/api/body/palettes")
    .then((index) => {
      paletteIndex = index;
      const choice = left.setPalettes(index);
      return applyPalette(choice.palette, choice.skin_tone);
    })
    .catch((error) => left.setPaletteApplied({ error: `Palettes unavailable: ${error.message}` }));
  left.setWholeBody(true);
  resetCamera();
  // Which structures this materialization draws is settled before the first one
  // is fetched, so a saved choice opens at its own fidelity rather than loading
  // one body and then another.
  try {
    fidelityIndex = await api("/api/body/materializations");
    const saved = readMaterialization();
    materialization = fidelityTier(saved) ? saved : fidelityIndex.default;
    applyFidelity(materialization);
    left.setMaterializations(materializationOptions(), materialization);
  } catch (error) {
    left.setMaterializations(materializationOptions(), materialization);
    left.setMaterializationNote("Materialization fidelity unavailable · " + error.message);
  }
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
    conformingDomains = index.domains.map((d) => ({ value: "conforming-domain-" + d.id, label: d.title }));
    left.setMaterializations(materializationOptions(), materialization);
  }).catch((error) => left.setMaterializationNote("Conforming domains unavailable: " + error.message));
  // The recorded trajectory is tens of megabytes and nothing has asked to play
  // it yet, so the scene is ready without it. setupFrames already reports an
  // absent trajectory rather than substituting anything for it.
  api("/api/body/trajectory?view=display")
    .then((data) => {
      if (!validBodyTrajectory(data)) throw Error("No valid computed body frames available");
      bodyTrajectory = data;
      bodyError = "";
    })
    .catch((error) => { bodyError = error.message; })
    .finally(setupFrames);
}
start();
