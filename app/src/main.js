import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import {
  filterStructures,
  anatomyViews,
  anatomyView,
  geometryArrays,
  chartPath,
  scenarioInput,
  spectralSeries,
  trajectoryFromChannels,
  scalarCoordinates,
  defaultModelId,
  bodyTransform,
  validBodyTrajectory,
  deformSkinVertices,
} from "./state.js";
import "./style.css";
import { attachedHairPositions } from "./hair_motion.js";
import { RegionalView, ElectricRegionalView } from "./regional.js";
import { ClothingView } from "./clothing.js";
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
  `<header><div class="brand"><span class="brandmark">ih<span>•</span></span><div>INTEGRATED HUMAN<small>Computational physiology workbench</small></div></div><nav class="workspace-tools" aria-label="Workspace panels"><button id="toggle-library" aria-label="Anatomy panel" aria-controls="library-panel" aria-expanded="true"><span aria-hidden="true">◧</span> Anatomy</button><button id="toggle-inspector" aria-label="Inspector panel" aria-controls="inspector-panel" aria-expanded="true"><span aria-hidden="true">◨</span> Inspector</button><button id="toggle-signals" aria-label="Signals panel" aria-controls="signals-panel" aria-expanded="true"><span aria-hidden="true">▤</span> Signals</button><button id="toggle-focus" aria-label="Focus mode" aria-pressed="false" title="Hide panels; click again to restore your workspace"><span aria-hidden="true">⛶</span> Focus</button></nav></header><main id="workspace"><aside class="library" id="library-panel" aria-label="Anatomy library" tabindex="0"><div class="section-heading"><span>GENERIC HUMAN</span><span id="total">—</span></div><button id="canonical-body" class="text-button">Explore the assembled body</button><p id="body-status" class="muted">Loading the canonical body state.</p><details id="source-inspection" class="evidence-fold"><summary>Advanced · source evidence inspection</summary><label class="field-label" for="model">Body or source model</label><select id="model"><option>Loading source families…</option></select></details><p class="muted" id="model-note">Reading the anatomical assembly.</p><label class="field-label" for="anatomy-view">Explore anatomy</label><select id="anatomy-view" aria-label="Anatomy view"></select><p id="layer-note" class="muted"></p><div class="search"><span>⌕</span><input id="search" placeholder="Find a structure…" aria-label="Search anatomy"></div><div class="section-heading layer-title">SYSTEM LAYERS <button id="all" class="text-button">Show all</button></div><div id="systems"></div><div class="section-heading results-heading">STRUCTURES <span id="count">0</span></div><div id="structures" class="structures"></div><div class="library-footer"><span class="live-dot"></span> One body · traceable sources and assumptions</div></aside><section class="center"><div class="viewport" id="viewport"><div class="view-heading"><div class="eyebrow">SPATIAL EXPLORER</div><h1 id="view-title">Human anatomy</h1><p id="frame-label">Source-defined coordinates</p></div><div class="view-tools"><button id="reset" title="Reset camera">↺ <span>Reset view</span></button><button id="posture" title="Rigid display rotation only; native posture is not established">Supine view</button><button id="front">Anterior</button><button id="side">Lateral</button></div><div id="scene-status" role="status">Loading anatomical assets…</div><div id="flow-legend" class="flow-legend" hidden></div><div class="orientation">Y <span>↑</span><br><small>X →</small></div><div class="view-footer"><span id="render-count">0 structures visible</span><span>Drag to orbit · Scroll to zoom · Click to inspect</span></div></div><div class="display-controls"><label>Opacity <input id="opacity" type="range" min=".05" max="1" step=".05" value="1"><output id="opacity-value">100%</output></label><label>Section <input id="clip" type="range" min="0" max="100" value="100"><output id="clip-value">Off</output></label><select id="flow-field" aria-label="Vascular field"><option value="velocity">Velocity</option><option value="pressure">Pressure</option></select><button id="play" disabled aria-label="Play spatial time series">▶</button><input id="time" type="range" min="0" max="0" value="0" aria-label="Spatial time frame" disabled><output id="time-value">No flow frames</output></div><section class="signals" id="signals-panel" aria-label="Physiology signals" tabindex="0"><div class="signal-header"><div><div class="eyebrow">TEMPORAL OBSERVATORY</div><h2>Physiology & dynamics</h2></div><div class="tabs"><button id="tab-phys" class="active">Trajectory</button><button id="tab-spectral">Spectrum</button></div></div><div id="spectral-controls" class="spectral-controls" hidden><select id="spectral-run" aria-label="Spectral evidence source"></select><select id="spectral-mode" aria-label="Spectral representation"><option value="psd">Fourier power density</option><option value="laplace">Finite Laplace magnitude</option></select><select id="sigma" aria-label="Laplace damping" hidden></select></div><div id="trajectory-controls" class="spectral-controls"><select id="trajectory-run" aria-label="Recorded physiology run"><option value="baseline">Recorded resting baseline</option></select></div><div class="signal-select"><select id="variable" aria-label="Physiology variable"></select><span id="signal-source">Awaiting simulation data</span></div><div id="chart" class="chart"><p class="empty">Loading physiological trajectories…</p></div><p id="chart-note" class="muted">Native simulation evidence is distinct from human measurements.</p></section></section><aside class="inspector" id="inspector-panel" aria-label="Structure inspector and scenarios" tabindex="0"><div class="section-heading">STRUCTURE INSPECTOR <span>↗</span></div><div id="details"><div class="inspector-symbol">◎</div><h2>Explore the body</h2><p class="muted">Select a structure in the viewport or library to inspect its source, coordinate frame, and evidence status.</p></div><div class="scenario"><div class="eyebrow">NATIVE PHYSIOLOGY</div><h2>Run a scenario</h2><p id="scenario-description" class="muted">Execute the body model and inspect its recorded response.</p><form id="scenario-form"><label class="field-label" for="engine-variant">Native implementation</label><select id="engine-variant"><option value="">Loading implementations…</option></select><p id="variant-note" class="muted">Source corrections do not establish clinical calibration.</p><label class="field-label" for="patient">Model profile</label><select id="patient"><option value="StandardMale">StandardMale</option></select><label class="field-label" for="scenario">Protocol</label><select id="scenario"><option value="baseline">Resting baseline</option><option value="exercise">Exercise & recovery</option><option value="hemorrhage_saline">Hemorrhage, saline & recovery</option></select><p id="protocol-note" class="muted">Record the upstream resting state.</p><label class="field-label" for="duration">Duration · seconds</label><input id="duration" type="number" min="1" max="600" value="60" required><label class="field-label" for="ambient">Ambient temperature · °C · optional</label><input id="ambient" type="number" min="10" max="35" step="any" placeholder="Upstream conditions"><label class="field-label" for="clothing">Clothing · clo · optional</label><input id="clothing" type="number" min="0" max="3" step="any" placeholder="Upstream conditions"><p class="muted">Blank values preserve upstream environment and clothing.</p><button class="primary" id="run" type="submit">Run native simulation <span>↗</span></button></form><div id="run-status" class="run-status" role="status">Checking native backend…</div></div><div class="evidence"><div class="section-heading">EVIDENCE BOUNDARY</div><p>Geometry describes anatomy. Simulation describes model behavior. Neither alone establishes patient-specific calibration.</p><div id="evidence-status" class="tag">Uncertainty remains explicit</div><details class="evidence-fold"><summary>Human measurement fit</summary><div id="calibration-details"><p>Loading calibration evidence…</p></div></details><details class="evidence-fold"><summary>System coverage</summary><div id="coverage-details"><p>Coverage data unavailable.</p></div></details><details class="evidence-fold"><summary>Conservative exchange</summary><div id="coupling-details"><p>Exchange data unavailable.</p></div></details><details class="evidence-fold"><summary>Vascular CFD audit</summary><div id="vascular-audit"><p>Audit unavailable.</p></div></details></div></aside></main>`;
// Keep panels mounted so filters, forms and disclosure state survive.
const compactWorkspace = matchMedia("(max-width: 1000px)");
const workspaceStates = new Map();
const panelScroll = new Map();
function workspaceState() {
  const mode = compactWorkspace.matches ? "compact" : "desktop";
  if (!workspaceStates.has(mode)) {
    let saved;
    try { saved = JSON.parse(localStorage.getItem(`ihm.workspace.${mode}`)); } catch {}
    const state = { library: !compactWorkspace.matches, inspector: !compactWorkspace.matches, signals: !compactWorkspace.matches, focus: false };
    for (const key of Object.keys(state))
      if (typeof saved?.[key] === "boolean") state[key] = saved[key];
    if (compactWorkspace.matches && state.library && state.inspector) state.inspector = false;
    workspaceStates.set(mode, state);
  }
  return workspaceStates.get(mode);
}
function renderWorkspace() {
  const state = workspaceState();
  for (const name of ["library", "inspector", "signals"]) {
    const panel = $(`${name}-panel`);
    const visible = state[name] && !state.focus;
    if (!panel.hidden && !visible) panelScroll.set(name, panel.scrollTop);
    const wasHidden = panel.hidden;
    $("workspace").dataset[name] = String(visible);
    panel.hidden = !visible;
    if (visible && wasHidden) panel.scrollTop = panelScroll.get(name) || 0;
    $(`toggle-${name}`).setAttribute("aria-expanded", String(visible));
  }
  $("toggle-focus").setAttribute("aria-pressed", String(state.focus));
}
function saveWorkspace() {
  try {
    localStorage.setItem(`ihm.workspace.${compactWorkspace.matches ? "compact" : "desktop"}`, JSON.stringify(workspaceState()));
  } catch { /* Storage may be unavailable; current-session controls still work. */ }
  renderWorkspace();
}
for (const name of ["library", "inspector", "signals"]) {
  $(`toggle-${name}`).onclick = () => {
    const state = workspaceState();
    state[name] = state.focus || !state[name];
    state.focus = false;
    if (compactWorkspace.matches && state[name] && name !== "signals")
      state[name === "library" ? "inspector" : "library"] = false;
    saveWorkspace();
  };
}
$("toggle-focus").onclick = () => {
  workspaceState().focus = !workspaceState().focus;
  saveWorkspace();
};
document.addEventListener("keydown", event => {
  if (event.key !== "Escape" || !compactWorkspace.matches) return;
  const name = ["library", "inspector"].find(name => !$(`${name}-panel`).hidden);
  if (!name) return;
  workspaceState()[name] = false;
  saveWorkspace();
  $(`toggle-${name}`).focus();
});
compactWorkspace.addEventListener("change", renderWorkspace);
renderWorkspace();
import { mountPanelResizers } from './workspace.js';
mountPanelResizers();
const bodyControls=document.createElement('div');
bodyControls.id='body-controls';
bodyControls.innerHTML=`<label class="field-label" for="chest-compliance">Chest compliance · L/cmH₂O · optional</label><input id="chest-compliance" type="number" min=".05" max="1" step=".01" placeholder="Native constitutive value"><p class="muted">Changes the physiological pressure–flow solve. Tissue parameters remain uncertain.</p><label class="field-label" for="body-protocol">Body perturbation</label><select id="body-protocol"><option value="none">No somatic stimulus</option value="touch">Left palm pressure pulse</option value="motor">Left tibialis anterior command</option value="block">Same command with nerve block</option></select><p class="muted">Explicit test inputs; no inferred voluntary motor policy.</p>`;
$('run').before(bodyControls);
const regionalControls=document.createElement('div');
regionalControls.id='regional-controls';
regionalControls.innerHTML='<label class="field-label" for="regional-study">Active materialization</label><select id="regional-study"><option value="body">Whole-body replay</option><option value="forearm-touch">Forearm · contact & IBM receptors</option><option value="skin-electric">Skin · non-neural electricity</option></select><select id="regional-condition" aria-label="Regional electrical intervention" hidden><option value="wound_shunt">Barrier shunt</option><option value="baseline">Intact baseline</option><option value="electrode">Electrode pair</option><option value="membrane_perturbation">Membrane perturbation</option></select><p id="regional-note" class="muted">Detailed studies share this body’s material coordinates.</p>';
$('anatomy-view').before(regionalControls);
const clothingControls=document.createElement('div');
clothingControls.id='clothing-components';
clothingControls.innerHTML='<div class="section-heading">CLOTHING</div><label class="system-row"><input id="garment-shirt" type="checkbox" checked><span>Sleeveless shirt</span></label><label class="system-row"><input id="garment-shorts" type="checkbox" checked><span>Shorts</span></label><p id="clothing-status" class="muted">Fitting garments to the assembled body…</p><details class="evidence-fold"><summary>Garment construction & mechanics</summary><p>Engineered patterns fitted to canonical skin sections with assumed ease. Shirt motion follows the computed thoracic field. This display has no solved fabric friction or deformable genital contact. Garment visibility is independent of the native thermal clothing input.</p></details>';
$('layer-note').after(clothingControls);
let clothingView=null,clothingRequest=0;
for(const id of ['shirt','shorts'])$('garment-'+id).onchange=()=>{
  clothingView?.setEnabled(id,$('garment-'+id).checked);updateClothingStatus();
};
function updateClothingStatus(){
  const count=[...clothingView?.meshes.values()||[]].filter(m=>m.visible).length;
  if(clothingView?.meshes.size)$('clothing-status').textContent=`Fitted to canonical skin · ${count} garment${count===1?'':'s'} shown · geometric pattern prior`;
}
async function loadClothing(){
  const request=++clothingRequest;
  clothingView?.dispose();clothingView=null;
  $('clothing-components').hidden=modelId!=='ihm-body';
  if(modelId!=='ihm-body'||!group)return;
  const skin=manifest.structures.find(s=>s.model_id===modelId&&s.id==='body-bp3d-FJ2810');
  if(!skin){$('clothing-status').textContent='Garment fit unavailable: canonical skin reference is missing.';return;}
  $('clothing-status').textContent='Fitting garments to the assembled body…';
  try {
    const geometry=await api(skin.geometry_url||`/api/geometry/${encodeURIComponent(skin.id)}`);
    if(request!==clothingRequest||modelId!=='ihm-body')return;
    clothingView=new ClothingView(group);clothingView.fit(geometry,skin);
    for(const id of ['shirt','shorts'])clothingView.setEnabled(id,$('garment-'+id).checked);
    updateClothingStatus();updateDisplay();updateFrame();
  }catch(error){
    if(request!==clothingRequest)return;
    clothingView?.dispose();clothingView=null;
    $('clothing-status').textContent=`Garment fit unavailable: ${error.message}`;
  }
}

let manifest,
  modelId,
  bodyTrajectory = null,
  bodySummary = null,
  bodyError = "Body trajectory unavailable",
  bodyRequest = 0,
  canonicalRuns = new Set(),
  systems = new Set(),
  layerOpacity = new Map(),
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
  hair: "#5a4030",
  microvascular: "#cf6679",
};
const viewport = $("viewport");
let renderer, scene, camera, controls, group, webglError;
let regionalView, regionalRequest=0, regionalSpectra=[], regionalReturnPhys=null,regionalReturnLabels=null;
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
    if (!w || !h) return;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }).observe(viewport);
  renderer.setAnimationLoop((now) => {
    controls.update();
    const frameIndex = Number($("time").value);
    const bodyInterval = regionalView?.active ? (regionalView instanceof ElectricRegionalView ? 1000*regionalView.data.clock.dt_s : 33) : modelId === "ihm-body" && bodyTrajectory
      ? 1000 * ((bodyTrajectory.frames[frameIndex+1]?.time_s ?? bodyTrajectory.frames[frameIndex].time_s + .1) - bodyTrajectory.frames[frameIndex].time_s)
      : 160;
    if (playing && flowFrames && now - lastTick >= bodyInterval) {
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
    if (regionalView?.active) return;
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
  if(regionalView?.active){regionalView.resetCamera(direction);return;}
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
function syncLayers() {
  $("systems").querySelectorAll('input[type="checkbox"]').forEach(el=>{
    el.checked=systems.has(el.value);
    const slider=$("systems").querySelector(`[data-opacity="${el.value}"]`);
    slider.disabled=!el.checked;
    slider.value=layerOpacity.get(el.value) ?? 1;
  });
  const total=manifest.structures.filter(s=>s.model_id===modelId).length;
  const selectedCount=manifest.structures.filter(s=>s.model_id===modelId && systems.has(s.system)).length;
  $("layer-note").textContent=`${selectedCount} of ${total} anatomical structures enabled. Each layer has its own opacity.`;
  $("all").textContent=systems.size===$("systems").querySelectorAll('input[type="checkbox"]').length?'Hide all':'Show all';
}
function applyAnatomyView() {
  if ($("anatomy-view").value==='custom') return;
  const available=[...new Set(manifest.structures.filter(s=>s.model_id===modelId).map(s=>s.system))];
  if($('anatomy-view').value==='clothed'){
    ({systems, opacity:layerOpacity}=anatomyView('core',available));
    if(available.includes('integumentary'))systems.add('integumentary');
  }else ({systems, opacity:layerOpacity}=anatomyView($("anatomy-view").value,available));
  $("search").value='';
  $("opacity").value='1';
  syncLayers();
  refresh();
  updateDisplay();
}
function rebuildSystems() {
  const rows = manifest.structures.filter((x) => x.model_id === modelId);
  const names = [...new Set(rows.map((x) => x.system))];
  systems = new Set(rows.filter(s=>s.default_visible).map(s=>s.system));
  if (!systems.size || modelId.startsWith("opensim")) systems = new Set(names);
  layerOpacity=new Map(names.map(s=>[s,1]));
  $("anatomy-view").innerHTML=anatomyViews.map(v=>{
    const matching=anatomyView(v.id,names).systems;
    const count=rows.filter(s=>matching.has(s.system)).length;
    return `<option value="${v.id}" ${count?'':'disabled'}>${esc(v.label)} · ${count}</option>`;
  }).join('')+(modelId==='ihm-body'?'<option value="clothed">Dressed human · skin & internal anatomy</option>':'')+'<option value="custom">Custom layers</option>';
  $("anatomy-view").value='custom';
  $("systems").innerHTML = names.map(name=>
    `<div class="system-item"><label class="system-row"><input type="checkbox" value="${esc(name)}" ${systems.has(name) ? "checked" : ""}><i style="background:${colors[name] || "#91abb0"}"></i><span>${esc(name.replaceAll("_", " "))}</span><small>${rows.filter(x=>x.system===name).length}</small></label><input class="layer-opacity" type="range" min=".05" max="1" step=".01" value="1" data-opacity="${esc(name)}" aria-label="${esc(name)} layer opacity"></div>`).join('');
  $("systems").querySelectorAll('input[type="checkbox"]').forEach(el=>el.onchange=()=>{
    el.checked ? systems.add(el.value) : systems.delete(el.value);
    $("anatomy-view").value='custom';syncLayers();refresh();
  });
  $("systems").querySelectorAll('[data-opacity]').forEach(el=>el.oninput=()=>{
    layerOpacity.set(el.dataset.opacity,Number(el.value));updateDisplay();
  });
  syncLayers();
}
function chooseModel() {
  regionalRequest++;
  closeRegional();
  $('regional-study').value='body';
  $('regional-condition').hidden=true;
  modelId = $("model").value;
  $("total").textContent = manifest.structures.filter(s => s.model_id === modelId).length.toLocaleString();
  playing = false;
  $("patient").disabled = modelId === "ihm-body";
  $("body-controls").hidden = modelId !== "ihm-body";
  $('regional-controls').hidden = modelId !== 'ihm-body';
  $("scenario-description").textContent = modelId === "ihm-body"
    ? "Run this generic body's physiology, mechanics and brain model. Playback uses computed states."
    : "Run a native source profile and inspect its recorded response.";
  syncCanonicalProfile();
  if (modelId === "ihm-body" && activeRun === "baseline") activeRun = "body";
  if (modelId === "ihm-body" && activeRun === "body" && bodyTrajectory) useBodyPhysiology();
  selected = null;
  loadController?.abort();
  objects.forEach(dispose);
  objects.clear();
  clothingRequest++;
  clothingView?.dispose();clothingView=null;
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
  if (modelId === "ihm-body") {
    const available = [...new Set(manifest.structures.filter(x => x.model_id === modelId).map(x => x.system))];
    ({systems, opacity: layerOpacity} = anatomyView("core", available));
    if(available.includes('integumentary'))systems.add('integumentary');
    $("anatomy-view").value = "clothed";
    syncLayers();
  }
  setupFrames();
  resetCamera();
  refresh();
  loadClothing();
}
function syncCanonicalProfile() {
  if (modelId !== "ihm-body") return;
  if (![...$("patient").options].some(o => o.value === "IHMGenericMale"))
    $("patient").add(new Option("IHMGenericMale · this body", "IHMGenericMale"));
  $("patient").value = "IHMGenericMale";
}
function useBodyPhysiology() {
  if (!bodyTrajectory) return;
  const frames = bodyTrajectory.frames;
  const channels = [...new Set(frames.flatMap(f => Object.keys(f.physiology || {})))].filter(key => frames.some(f => Number.isFinite(f.physiology?.[key])));
  receivePhysiology({time_s: frames.map(f => f.time_s), values: Object.fromEntries(channels.map(key => [key, frames.map(f => Number.isFinite(f.physiology?.[key]) ? f.physiology[key] : null)])), metadata: {source_kind: "Computed generic body physiology"}});
  updateVariables();
}
function receivePhysiology(value){
  if(regionalReturnPhys!==null)regionalReturnPhys=value;
  else phys=value;
}
function closeRegional(){
  regionalView?.close();
  if(regionalReturnPhys!==null){phys=regionalReturnPhys;regionalReturnPhys=null;}
  if(regionalReturnLabels){
    $('view-title').textContent=regionalReturnLabels.title;$('frame-label').textContent=regionalReturnLabels.frame;
    $('render-count').textContent=`${[...objects.values()].filter(o=>o.visible).length} structures visible`;
    regionalReturnLabels=null;
  }
  $('scene-status').hidden=false;
  $('regional-condition').hidden=true;
  updateSigma();updateVariables();
}
async function loadBodyTrajectory(run) {
  const request = ++bodyRequest;
  try {
    const data = await api("/api/body/trajectory?view=display" + (run ? `&run=${encodeURIComponent(run)}` : ""));
    if (request !== bodyRequest) return;
    if (!validBodyTrajectory(data)) throw Error("No valid computed body frames available");
    bodyTrajectory = data;
    bodyError = "";
    if (activeRun === "body") useBodyPhysiology();
  } catch (error) {
    if (request !== bodyRequest) return;
    bodyTrajectory = null;
    bodyError = error.message;
  }
  if (modelId === "ihm-body") setupFrames();
  $("body-status").textContent = bodySummary
    ? `${bodySummary.entity_count?.toLocaleString() || "Canonical"} anatomical entities · ${bodyTrajectory ? "computed body playback available" : "reference anatomy; body trajectory unavailable"}`
    : bodyTrajectory ? "Computed generic body state available" : "Reference anatomy available; body trajectory unavailable";
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
    if (current !== generation || regionalView?.active) return;
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
  if (g.attachment?.kind === "MaterialPoint" && g.attachment.reference_triangles_m) {
    object.userData.hairAttachment = g.attachment;
    object.userData.hairReference = positions.slice();
  }
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
function canonicalEvidence(s) {
  if (s.model_id !== "ihm-body") return "";
  const u = s.uncertainty || {};
  const registration = model()?.registrations?.[u.registration_or_synthesis?.registration_id];
  const assumptions = (s.assumptions || []).map(id => {
    const item = model()?.assumption_ledger?.find(a => a.id === id);
    return `<li><strong>${esc(id)}</strong><br>${esc(item?.statement || "See canonical assumption ledger")}</li>`;
  }).join("");
  const numerical = u.display_numerics || {};
  const triangleCount = s.display_geometry?.source_faces;
  return `<section class="canonical-evidence"><h3>What this representation supports</h3><dl><dt>Biological uncertainty</dt><dd>${esc(u.biological?.status || "Biological variation is not quantified for this generic body.")} No calibrated confidence percentage.</dd><dt>Registration & synthesis</dt><dd>${esc((s.evidence_kind || u.registration_or_synthesis?.kind || "Source-backed representation").replaceAll("_", " "))}.${registration?.held_out_rms_m !== undefined ? ` Held-out geometric fit: ${(registration.held_out_rms_m*1000).toFixed(2)} mm RMS. This is a registration diagnostic, not biological accuracy.` : " See the explicit assumptions below."}</dd><dt>Display & numerical precision</dt><dd>${triangleCount === undefined ? "Structural geometry" : `${Number(triangleCount).toLocaleString()} acquired triangles`}; ${numerical.display_reduction_applied === false ? "no additional display reduction" : "see representation metadata"}. Coordinates displayed in meters; source surface resolution and computed transform precision do not establish biological certainty.</dd><dt>Canonical frame</dt><dd>${esc(model()?.frame)} · meters</dd></dl>${assumptions ? `<details class="evidence-fold"><summary>Source dependencies & modeling assumptions</summary><ul>${assumptions}</ul></details>` : ""}</section>`;
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
    `<div class="eyebrow">${esc(s.system)} / ${esc(s.kind || "mesh")}</div><h2>${esc(s.name)}</h2><span class="tag">${esc(s.calibration_status || "Not independently calibrated")}</span><dl><dt>Source</dt><dd>${esc(source.label || model()?.name)}</dd><dt>Original source frame</dt><dd>${esc(source.frame || model()?.frame || "Source-defined")}</dd><dt>Specimen</dt><dd>${esc(source.specimen || "See source metadata")}</dd><dt>Original source units</dt><dd>${esc(source.units || model()?.source_units || "Source units")}</dd><dt>Evidence status</dt><dd>${esc(source.status || s.status || "Source geometry; physiological registration unverified")}</dd>${s.path_interpretation ? `<dt>Muscle path interpretation</dt><dd>${esc(s.path_interpretation)}</dd><dt>Wrapping solution</dt><dd>${s.wrap_solved ? "Solved by source" : "Unsolved; attachment chords only"}</dd>` : ""}${
      s.native_mechanics
        ? `<dt>Native mechanics · source simulation</dt><dd>Activation ${esc(s.native_mechanics.activation)} · path length ${Number(s.native_mechanics.length_m).toFixed(4)} m<br>Tendon force ${Number(s.native_mechanics.tendon_force_N).toFixed(3)} N<br>Active / passive fiber force ${Number(s.native_mechanics.active_fiber_force_N).toFixed(3)} / ${Number(s.native_mechanics.passive_fiber_force_N).toFixed(3)} N</dd><dt>Joint moment arms</dt><dd>${Object.entries(
            s.native_mechanics.moment_arms_m || {},
          )
            .map(([k, v]) => `${esc(k)}: ${Number(v).toFixed(5)} m`)
            .join(
              "<br>",
            )}</dd><dt>Mechanical validity</dt><dd>Source default pose; external force balance ${s.native_mechanics.external_force_balance_solved ? "solved" : "not solved"}. No subject calibration.</dd>`
        : ""
    }${s.classification ? `<dt>Anatomical grouping</dt><dd>${esc(s.classification.classification_basis || s.classification.status)} · ${esc((s.classification.systems || []).join(", "))}${s.classification.ambiguous_primary ? " · multiple supported primary groups" : ""}</dd>` : ""}${s.source_collections ? `<dt>Source collections</dt><dd>${esc(s.source_collections.join(", "))}</dd>` : ""}${s.evaluation_warnings?.length ? `<dt>Source evaluation caveat</dt><dd>${esc(s.evaluation_warnings.join("; "))}</dd>` : ""}${source.geometry_stage ? `<dt>Geometry stage</dt><dd>${esc(source.geometry_stage)}</dd>` : ""}${s.display_geometry ? `<dt>Surface geometry</dt><dd>${esc(s.display_geometry.source_faces)} original triangles · display budget ${esc(s.display_geometry.target_faces ?? s.display_geometry.source_faces)}. Original mesh retained.</dd>` : ""}${source.license ? `<dt>Asset license</dt><dd>${esc(source.license)}</dd>` : ""}${source.sha256 ? `<dt>SHA-256</dt><dd class="hash">${esc(source.sha256)}</dd>` : ""}</dl>${/^https?:\/\//.test(source.url || "") ? `<a class="source-link" href="${esc(source.url)}" target="_blank" rel="noopener noreferrer">Open original source ↗</a>` : ""}`;
  $("details").querySelector(":scope > dl")?.insertAdjacentHTML("beforebegin", canonicalEvidence(s));
  if(s.id==='body-detail-hair'||s.id==='body-detail-microvascular') {
    api('/api/body/experiments/details').then(data=>{
      if(selected?.id!==s.id)return;
      const content=s.id==='body-detail-hair'
        ? `<h3>Follicle population and display</h3><p>${data.hair.full_count.toLocaleString()} material-attached follicles in the retained regional sample; ${data.hair.display_count.toLocaleString()} displayed. Physical radii are not enlarged.</p><p>${esc(data.hair.morphology_prior)}</p><p>Regional densities derive from a 15-person human subset; their SD describes that sample, not this generic body's posterior certainty. Glabrous and unresolved masks remain excluded.</p>`
        : `<h3>Perfused reference graphs</h3><p>${data.microvascular.unit_count} bilateral forearm units · ${data.microvascular.capillary_count} capillary links · ${data.microvascular.edge_count} total edges.</p><p>Paired supply/return graphs solve pressure and flow with prescribed boundaries. Native blood storage is not connected; radii, topology and apparent viscosity remain priors.</p><p>Maximum internal flow residual: ${data.microvascular.maximum_internal_residual_m3_s.toExponential(2)} m³/s.</p>`;
      $('details').querySelector('.canonical-evidence').innerHTML=content;
    }).catch(error=>{if(selected?.id===s.id)$('details').insertAdjacentHTML('beforeend',`<p class="muted">${esc(error.message)}</p>`);});
  }
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
        const alpha = obj.userData.structure?.kind === "vectors" ? 1 : opacity * (layerOpacity.get(obj.userData.structure?.system) ?? 1);
        o.material.opacity = alpha;
        o.material.transparent = alpha < 1 || !!o.material.isLineBasicMaterial;
        o.material.depthWrite = alpha > 0.5;
        o.material.clippingPlanes = clip < 100 ? [plane] : [];
      }
    }),
  );
  clothingView?.setClipping(clip < 100 ? [plane] : []);
  $("opacity-value").textContent = `${Math.round(opacity * 100)}%`;
  $("clip-value").textContent = clip === 100 ? "Off" : `${clip}%`;
}
function setupFrames() {
  $('posture').disabled=!!regionalView?.active;
  if (modelId === "ihm-body") {
    flowFrames = regionalView?.active ? regionalView.data.frames.length : bodyTrajectory?.frames.length || 0;
    $("flow-field").hidden = true;
    $("play").disabled = flowFrames < 2;
    $("time").disabled = !flowFrames;
    $("time").max = Math.max(0, flowFrames - 1);
    $("time").value = 0;
    playing = false;
    $("play").textContent = "▶";
    updateFrame();
    return;
  }
  $("flow-field").hidden = false;
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
  if (regionalView?.active) {
    const frame=regionalView.draw(Number($('time').value));
    $('flow-legend').hidden=false;
    $('flow-legend').textContent=regionalView instanceof ElectricRegionalView
      ? `Non-neural electrical state · Vm ${(frame.membrane_voltage_V[4]*1000).toFixed(2)} mV · apical ${(frame.apical_voltage_V[4]*1000).toFixed(2)} mV · peak |E| ${Math.max(...frame.edge_field_V_m.map(Math.abs)).toFixed(2)} V/m`
      : `Computed contact · ${(frame.indentation_m*1e6).toFixed(1)} µm · ${(frame.reaction_n*1e3).toFixed(3)} mN · signed IBM response; slow playback`;
    $('time-value').textContent=`${frame.time_s.toFixed(3)} s`;
    $('time-value').title=regionalView instanceof ElectricRegionalView ? 'Recorded electrical clock; fixed ion reservoirs. Vm colors −80 to20 mV; apical colors −40 to10 mV.' : 'Recorded quasistatic mechanical endpoints and causal receptor state; playback slowed to inspect adaptation.';
    return;
  }
  if (modelId === "ihm-body") {
    const frame = bodyTrajectory?.frames[Number($("time").value)];
    const skinField = frame?.respiration?.skin_field;
    const skinIds = new Set(skinField?.entity_ids || []);
    let deformedSkins = 0;
    objects.forEach((object, id) => {
      const positions = object.geometry?.getAttribute("position");
      const hair = object.userData.hairAttachment;
      if (positions && hair) {
        attachedHairPositions(object.userData.hairReference, hair,
          skinIds.has(hair.skin_entity_id) ? skinField : null, positions.array);
        positions.needsUpdate = true;
        object.geometry.computeVertexNormals();
        object.geometry.computeBoundingSphere();
      }
      if (positions && (skinIds.has(id) || object.userData.skinReference)) {
        object.userData.skinReference ||= positions.array.slice();
        deformSkinVertices(object.userData.skinReference, skinIds.has(id) ? skinField : null, positions.array);
        positions.needsUpdate = true;
        object.geometry.computeVertexNormals();
        object.geometry.computeBoundingSphere();
        if (skinIds.has(id)) deformedSkins++;
      }
      const motionId = hair?.skin_entity_id || id;
      const transform = bodyTransform(frame?.entities?.[motionId], bodyTrajectory?.centroids_m[motionId]);
      object.matrixAutoUpdate = false;
      object.matrix.set(...transform);
      object.matrixWorldNeedsUpdate = true;
    });
    clothingView?.update(frame,bodyTrajectory?.centroids_m);
    $("flow-legend").hidden = !frame;
    const p = frame?.physiology || {};
    const values = [["HeartRate(1/min)", "HR", "/min"], ["MeanArterialPressure(mmHg)", "MAP", "mmHg"]]
      .filter(([key]) => Number.isFinite(p[key])).map(([key,label,unit]) => `${label} ${Number(p[key]).toFixed(1)} ${unit}`);
    $("flow-legend").textContent = frame ? `Computed body state · ${Object.keys(frame.entities).length} tissue transforms${deformedSkins ? ` · ${deformedSkins} thoracic skin field${deformedSkins === 1 ? "" : "s"}` : ""}${values.length ? " · " + values.join(" · ") : ""}` : "";
    $("time-value").textContent = frame ? `${Number(frame.time_s).toFixed(3)} s` : "Body trajectory unavailable";
    $("time-value").title = frame ? "Computed snapshots; omitted tissues retain reference geometry. No interpolation or extrapolation." : bodyError;
    return;
  }
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
  if(regionalView?.active){
    const id=regionalView instanceof ElectricRegionalView ? 'skin-electric:'+regionalView.condition : 'forearm-touch';
    return regionalSpectra.find(r=>r.id===id);
  }
  return (
    temporal.runs?.find((r) => r.id === $("spectral-run").value) ||
    temporal.runs?.[0]
  );
}
function updateVariables() {
  $('spectral-run').disabled=!!regionalView?.active;
  $('trajectory-run').disabled=!!regionalView?.active;
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
    : regionalView?.active
      ? 'Computed regional materialization at a pinned body site. Explicit parameter priors; inspect the regional evidence and charge/force audit.'
    : activeRun === "body"
      ? "Computed physiology for IHMGenericMale on the body playback clock. Recorded model states; empirical calibration remains incomplete."
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
    syncCanonicalProfile();
    canonicalRuns = new Set(runs.filter(r => r.canonical_body).map(r => r.id));
    const latest = runs.at(-1);
    $("trajectory-run").innerHTML =
      '<option value="body">Computed generic body state</option><option value="baseline">Recorded source resting baseline</option>' +
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
      (regionalReturnPhys || phys)._run !== chosen.id
    ) {
      const received = await api(`/api/physiology?run=${encodeURIComponent(chosen.id)}`);
      received._run = chosen.id;receivePhysiology(received);
      if (chosen.canonical_body) await loadBodyTrajectory(chosen.id);
      updateVariables();
    }
  } catch (e) {
    $("run-status").textContent = e.message;
  }
}
$("trajectory-run").onchange = async () => {
  activeRun = $("trajectory-run").value;
  const requestedRun=activeRun;
  try {
    if (activeRun === "body") { useBodyPhysiology(); return; }
    const received =
      requestedRun === "reproductive"
        ? trajectoryFromChannels(reproductive)
        : requestedRun.startsWith("thermal:")
          ? trajectoryFromChannels(
              await api(
                `/api/thermal?run=${encodeURIComponent(requestedRun.slice(8))}`,
              ),
            )
          : requestedRun.startsWith("csf:")
            ? trajectoryFromChannels(
                await api(
                  `/api/csf?run=${encodeURIComponent(requestedRun.slice(4))}`,
                ),
              )
            : await api(
                requestedRun === "baseline"
                  ? "/api/physiology"
                  : `/api/physiology?run=${encodeURIComponent(requestedRun)}`,
              );
    if(requestedRun!==activeRun)return;
    received._run = requestedRun;receivePhysiology(received);
    if (canonicalRuns.has(requestedRun)) await loadBodyTrajectory(requestedRun);
    updateVariables();
  } catch (e) {
    if(requestedRun===activeRun && !regionalView?.active)$("chart").innerHTML = `<p class="empty">${esc(e.message)}</p>`;
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
      modelId === "ihm-body" ? "IHMGenericMale" : $("patient").value,
      {
        ambient_temperature_c: $("ambient").value,
        clothing_clo: $("clothing").value,
        engine_variant: $("engine-variant").value,
      },
    );
    if(modelId==='ihm-body') {
      if($('chest-compliance').value!=='')config.chest_compliance_l_cmH2O=Number($('chest-compliance').value);
      const protocol=$('body-protocol').value;
      if(protocol!=='none') {
        const event={start_s:config.seconds*.2,end_s:config.seconds*.5};
        if(protocol==='touch')event.stimuli={'peripheral-skin-left-palm':{pressure_pa:20000}};
        else {
          event.motor_commands={'body-muscle-opensim-tibant_l':.5};
          if(protocol==='block') {
            const data=await api('/api/body/peripheral');
            const binding=data.muscle_bindings.find(b=>b.muscle_id==='body-muscle-opensim-tibant_l');
            if(!binding)throw Error('Requested motor pathway is unavailable');
            event.blocked_nerves=[binding.nerve_id];
          }
        }
        config.body_interventions=[event];
      }
    }
    const run = await api(modelId === "ihm-body" ? "/api/body/scenarios" : "/api/scenarios", {
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
  if (!group || regionalView?.active) return;
  group.rotation.x = group.rotation.x === 0 ? -Math.PI / 2 : 0;
  $("posture").textContent = group.rotation.x ? "Upright view" : "Supine view";
  modelBounds.setFromObject(group);
  resetCamera();
  $("frame-label").textContent =
    `${group.rotation.x ? "Supine" : "Upright"} display · rigid rotation only · native posture not established`;
  updateDisplay();
};
$("model").onchange = chooseModel;
$('regional-study').onchange=async()=>{
  const request=++regionalRequest,kind=$('regional-study').value;
  playing=false;
  closeRegional();
  if(kind==='body'){
    updateSigma();
    $('regional-note').textContent='Detailed studies share this body’s material coordinates.';
    setupFrames();return;
  }
  regionalReturnPhys=phys;
  regionalReturnLabels={title:$('view-title').textContent,frame:$('frame-label').textContent};
  $('regional-note').textContent='Loading source-pinned regional experiment…';
  try {
    const data=await api('/api/body/experiments/'+encodeURIComponent(kind));
    if(request!==regionalRequest)return;
    if(!renderer)throw Error('Regional mesh needs WebGL');
    const View=kind==='skin-electric'?ElectricRegionalView:RegionalView;
    regionalView = new View({scene,camera,controls,bodyGroup:group});
    regionalView.open(data);
    $('scene-status').hidden=true;
    $('view-title').textContent=kind==='skin-electric'?'Skin · non-neural electricity':'Forearm · contact & sensation';
    $('frame-label').textContent='Pinned left forearm coordinates · meters · regional reference model';
    $('render-count').textContent=kind==='skin-electric'?'9 membrane sites · 9 apical sites':`${data.geometry.tetrahedra.length.toLocaleString()} volume tetrahedra · boundary surface shown`;
    try {regionalSpectra=(await api('/api/body/experiments/spectra')).runs;updateSigma();}
    catch {regionalSpectra=[];}
    if(request!==regionalRequest)return;
    if(kind==='skin-electric') {
      $('regional-condition').hidden=false;$('regional-condition').value='wound_shunt';
      $('regional-note').textContent='Separate membrane and apical networks · fixed basal bath. Physical site spacing; marker radii illustrate sites.';
      spectral=false;$('tab-phys').click();useElectricalSignals();setupFrames();
      $('details').innerHTML=`<h2>Non-neural skin electricity</h2><p>${esc(data.voltage_reference)}</p><dl><dt>Body attachment</dt><dd>${esc(data.anchor.entity_id)} · face ${data.anchor.face_index}</dd><dt>Parameter evidence</dt><dd>${esc(data.parameter_status)}</dd><dt>Charge balance</dt><dd>${data.experiments.wound_shunt.charge_audit.max_node_residual_c.toExponential(2)} C maximum node residual</dd><dt>Color scales</dt><dd>Membrane −80 to20 mV; apical −40 to10 mV. Blue → red within each separate scale.</dd></dl><p class="muted">${esc(data.limitations.join(' '))}</p>`;
      return;
    }
    $('regional-note').textContent='Left forearm · synthesized reference volume. Measured stiffness transfer prior; no whole-body force feedback.';
    spectral=false;$('tab-phys').click();
    phys={time_s:data.frames.map(f=>f.time_s),values:Object.fromEntries(['indentation_m','reaction_n','elastic_energy_j','rapid_response','slow_response'].map(k=>[k,data.frames.map(f=>f[k])])),units:{indentation_m:'m',reaction_n:'N',elastic_energy_j:'J',rapid_response:'donor response',slow_response:'donor response'},metadata:{source_kind:'Computed contact → pinned IBM causal transfer'}};
    updateVariables();setupFrames();
    $('details').innerHTML=`<h2>Forearm contact & sensation</h2><p>Body entity ${esc(data.anchor.entity_id)} · material face ${esc(data.anchor.face_index)}</p><dl><dt>Volume mesh</dt><dd>${data.geometry.tetrahedra.length} tetrahedra · physical meters</dd><dt>Stiffness</dt><dd>Apparent human shear modulus 2.8 ± 0.8 kPa; transferred from 20 adults aged55–70. Bulk response and local thickness remain priors.</dd><dt>Neural dynamics</dt><dd>Exact causal IBM rapid/slow transfer; signed response, no calibrated firing-rate conversion.</dd><dt>Numerical force residual</dt><dd>${data.states.loaded.force_balance_residual_n.toExponential(2)} N</dd></dl><p class="muted">${esc(data.limitations.join(' '))}</p>`;
  }catch(error){if(request===regionalRequest){closeRegional();$('regional-note').textContent=error.message;$('regional-study').value='body';setupFrames();}}
};
function useElectricalSignals(){
  if(!(regionalView instanceof ElectricRegionalView)||!regionalView.active)return;
  const frames=regionalView.data.frames;
  phys={time_s:frames.map(f=>f.time_s),values:{'Central membrane voltage':frames.map(f=>f.membrane_voltage_V[4]),'Central apical voltage':frames.map(f=>f.apical_voltage_V[4]),'Maximum absolute lateral field':frames.map(f=>Math.max(...f.edge_field_V_m.map(Math.abs)))},units:{'Central membrane voltage':'V','Central apical voltage':'V','Maximum absolute lateral field':'V/m'},metadata:{source_kind:'Computed non-neural RC network · '+regionalView.condition}};
  updateVariables();
}
$('regional-condition').onchange=()=>{
  if(!(regionalView instanceof ElectricRegionalView)||!regionalView.active)return;
  regionalView.selectCondition($('regional-condition').value);updateSigma();useElectricalSignals();setupFrames();
};
$("canonical-body").onclick = () => {
  if (!manifest?.models.some(m => m.id === "ihm-body")) return;
  $("model").value = "ihm-body";
  $("source-inspection").open = false;
  chooseModel();
};
$("anatomy-view").onchange = applyAnatomyView;
$("search").oninput = () => manifest && refresh();
$("all").onclick = () => {
  if (!manifest) return;
  const boxes = [...$("systems").querySelectorAll('input[type="checkbox"]')],
    show = boxes.some((b) => !b.checked);
  boxes.forEach((b) => {
    b.checked = show;
    show ? systems.add(b.value) : systems.delete(b.value);
  });
  $("anatomy-view").value="custom";
  if(show) layerOpacity.set("integumentary",.18);
  syncLayers();
  refresh();
  updateDisplay();
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
  lastTick = performance.now();
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
    $("model").value = defaultModelId(manifest.models);
    chooseModel();
    api("/api/body").then(data => { bodySummary = data; }).catch(() => {}).finally(() => loadBodyTrajectory());
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
  if (results[0].status === "fulfilled") receivePhysiology(results[0].value);
  if (activeRun === "body" && bodyTrajectory) useBodyPhysiology();
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
