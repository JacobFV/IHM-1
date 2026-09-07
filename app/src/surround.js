import * as THREE from "three";

// The world a scene declares: its enclosing geometry, the sky above it, the
// ground under it and the objects standing in it. Every field read here is
// served by /api/scene/catalog under `world`; nothing about a room's shape,
// colour or extent is decided in this file.
//
// All of it is scenery. The catalogue says so in world.physics — no collider
// exists for a wall, a floor sheet, a window or a piece of furniture — so every
// mesh built here refuses to be raycast and can never be selected as anatomy.

const UP = new THREE.Vector3(0, 1, 0);

function scenery(mesh) {
  mesh.raycast = () => {};
  mesh.userData.scenery = true;
  return mesh;
}

// The gravity frame the environment declares: the sky gradient runs along it and
// the ground plane lies across it.
function upAxis(environment) {
  const gravity = environment?.gravity;
  if (Array.isArray(gravity) && gravity.length === 3 && gravity.some((v) => v)) {
    return new THREE.Vector3(...gravity).normalize().negate();
  }
  const axis = environment?.axis;
  return Number.isInteger(axis)
    ? new THREE.Vector3(axis === 0 ? 1 : 0, axis === 1 ? 1 : 0, axis === 2 ? 1 : 0)
    : UP.clone();
}

function sampleStops(stops, t) {
  if (!stops?.length) return new THREE.Color(1, 1, 1);
  const sorted = [...stops].sort((a, b) => a.t - b.t);
  if (t <= sorted[0].t) return new THREE.Color(...sorted[0].rgb);
  const last = sorted[sorted.length - 1];
  if (t >= last.t) return new THREE.Color(...last.rgb);
  for (let i = 1; i < sorted.length; i++) {
    const a = sorted[i - 1], b = sorted[i];
    if (t > b.t) continue;
    const f = b.t === a.t ? 0 : (t - a.t) / (b.t - a.t);
    return new THREE.Color(...a.rgb).lerp(new THREE.Color(...b.rgb), f);
  }
  return new THREE.Color(...last.rgb);
}

// The sky is a declared gradient between stops, painted on the inside of a
// sphere large enough to sit behind everything, so the stops run along the
// environment's own gravity axis rather than up the screen.
function skyMesh(sky, up, radius) {
  if (!sky || sky.type === "none" || !sky.stops?.length) return null;
  const geometry = new THREE.SphereGeometry(radius, 32, 24);
  const position = geometry.getAttribute("position");
  const colors = new Float32Array(position.count * 3);
  const vertex = new THREE.Vector3();
  for (let i = 0; i < position.count; i++) {
    vertex.fromBufferAttribute(position, i);
    const t = (vertex.dot(up) / radius + 1) / 2;
    const color = sampleStops(sky.stops, t);
    colors[3 * i] = color.r; colors[3 * i + 1] = color.g; colors[3 * i + 2] = color.b;
  }
  geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  const mesh = new THREE.Mesh(geometry, new THREE.MeshBasicMaterial({
    vertexColors: true, side: THREE.BackSide,
    depthWrite: false, depthTest: false, fog: false,
  }));
  // The dome rides with the camera and is painted before anything else, so it
  // reads as sky at any distance instead of a ball the reader can leave behind.
  mesh.renderOrder = -1000;
  mesh.userData.sky = true;
  mesh.frustumCulled = false;
  return scenery(mesh);
}

// The engine plane, drawn as a grid. It is the only surface the engine solves,
// so it is drawn where the environment says it is and nowhere else.
function groundGrid(ground, up, extent) {
  if (!ground || ground.type !== "ideal_half_space" || !Number.isFinite(ground.level_m)) return null;
  const grid = new THREE.GridHelper(extent, Math.round(extent / 0.25), 0x4c6e70, 0x22383a);
  grid.material.transparent = true;
  grid.material.opacity = 0.55;
  grid.material.depthWrite = false;
  grid.quaternion.setFromUnitVectors(UP, up);
  grid.position.copy(up).multiplyScalar(ground.level_m);
  grid.renderOrder = -900;
  return scenery(grid);
}

// A surround is one indexed mesh with a colour per face and a flat list of
// parts. Faces are emitted part by part, which is checked here rather than
// assumed: if the blocks do not line up the whole enclosure is drawn and
// nothing is cut away, so a reader never silently loses a wall.
function surroundMesh(data, cutaway) {
  const parts = data.parts || [];
  const faces = data.indices.length / 3;
  const blocked = parts.length > 0 && faces % parts.length === 0 &&
    data.positions.length / 3 % parts.length === 0 &&
    parts.every((part, index) => {
      const vertices = data.positions.length / 3 / parts.length;
      const perFace = faces / parts.length;
      for (let f = index * perFace; f < (index + 1) * perFace; f++)
        for (let k = 0; k < 3; k++) {
          const vertex = data.indices[3 * f + k];
          if (vertex < index * vertices || vertex >= (index + 1) * vertices) return false;
        }
      return true;
    });
  const cut = new Set(blocked ? cutaway || [] : []);
  const perFace = parts.length ? faces / parts.length : faces;
  const keep = [];
  for (let f = 0; f < faces; f++) {
    const part = parts[Math.floor(f / perFace)];
    if (part && cut.has(part.name)) continue;
    keep.push(f);
  }
  const positions = new Float32Array(keep.length * 9);
  const colors = new Float32Array(keep.length * 9);
  const hasColors = data.face_colours_rgb?.length === faces * 3;
  keep.forEach((f, n) => {
    for (let k = 0; k < 3; k++) {
      const vertex = data.indices[3 * f + k];
      for (let c = 0; c < 3; c++) {
        positions[9 * n + 3 * k + c] = data.positions[3 * vertex + c];
        colors[9 * n + 3 * k + c] = hasColors ? data.face_colours_rgb[3 * f + c] : 0.7;
      }
    }
  });
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  geometry.computeVertexNormals();
  const mesh = new THREE.Mesh(geometry, new THREE.MeshStandardMaterial({
    vertexColors: true, roughness: 0.94, metalness: 0, side: THREE.DoubleSide,
  }));
  mesh.name = data.id;
  mesh.userData.cutaway = blocked ? [...cut] : [];
  mesh.userData.cutawayApplied = blocked;
  return scenery(mesh);
}

function objectMesh(data, colour, offset) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(new Float32Array(data.positions), 3));
  geometry.setIndex(new THREE.BufferAttribute(new Uint32Array(data.indices), 1));
  geometry.computeVertexNormals();
  const mesh = new THREE.Mesh(geometry, new THREE.MeshStandardMaterial({
    color: Array.isArray(colour) ? new THREE.Color(...colour) : 0x9aa7ad,
    roughness: 0.86, metalness: 0.02, side: THREE.DoubleSide,
  }));
  if (Array.isArray(offset) && offset.length === 3) mesh.position.fromArray(offset);
  mesh.name = data.id;
  return scenery(mesh);
}

export function mountSurround(scene, { api }) {
  const root = new THREE.Group();
  root.name = "surround";
  scene.add(root);
  const geometryCache = new Map();
  let generation = 0;

  const load = (url) => {
    if (!geometryCache.has(url)) geometryCache.set(url, api(url));
    return geometryCache.get(url);
  };
  function empty() {
    root.traverse((node) => {
      if (node.geometry) node.geometry.dispose();
      if (node.material) node.material.dispose();
    });
    root.clear();
  }

  return {
    group: root,
    // `entry` is the selected scene, or the environment when no scene is chosen;
    // `environment` is always the base environment record, because the sky axis
    // and the engine plane belong to it and a scene never changes them.
    async apply(entry, environment, instances = [], objectRecords = []) {
      const current = ++generation;
      empty();
      const world = entry?.world || environment?.world;
      if (!world) return { drawn: 0, note: "This environment declares no world." };
      const up = upAxis(environment);
      const extent = world.enclosure_dimensions_m
        ? Math.max(world.enclosure_dimensions_m.width, world.enclosure_dimensions_m.depth) * 1.2
        : (world.ground?.extent_m || 12);
      const sky = skyMesh(world.sky, up, 5);
      if (sky) root.add(sky);
      const grid = groundGrid(world.ground, up, Math.max(6, extent));
      if (grid) root.add(grid);
      let drawn = sky ? 1 : 0;
      const notes = [];
      if (world.surround_url) {
        try {
          const data = await load(world.surround_url);
          if (current !== generation) return { drawn, note: "" };
          const mesh = surroundMesh(data, world.cutaway_parts);
          if (!mesh.userData.cutawayApplied && (world.cutaway_parts || []).length)
            notes.push("The surround's parts do not line up with its faces, so the whole enclosure is drawn.");
          root.add(mesh);
          drawn++;
        } catch (error) { notes.push("Surround geometry unavailable: " + error.message); }
      }
      // A scene's own placements, plus whatever the reader has inserted. Both
      // are drawn the same way, because the engine instantiates neither.
      const records = new Map(objectRecords.map((o) => [o.id, o]));
      const placements = [
        ...(entry?.placements || []).map((p) => ({ id: p.object, offset: p.offset_m })),
        ...instances.map((id) => ({ id, offset: [0, 0, 0] })),
      ];
      await Promise.all(placements.map(async (placement) => {
        const record = records.get(placement.id);
        if (!record?.geometry_url) { notes.push("No geometry for " + placement.id); return; }
        try {
          const data = await load(record.geometry_url);
          if (current !== generation) return;
          root.add(objectMesh(data, record.colour_rgb, placement.offset));
          drawn++;
        } catch (error) { notes.push(placement.id + ": " + error.message); }
      }));
      if (current !== generation) return { drawn, note: "" };
      return { drawn, note: notes.join(" ") };
    },
    // Called every frame: the sky dome keeps the camera at its centre.
    follow(camera) {
      for (const child of root.children)
        if (child.userData.sky) child.position.copy(camera.position);
    },
    dispose() { generation++; empty(); root.removeFromParent(); },
  };
}
