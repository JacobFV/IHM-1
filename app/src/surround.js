import * as THREE from "three";
import { RoundedBoxGeometry } from "three/addons/geometries/RoundedBoxGeometry.js";
import { environmentMaterial, materialUVs } from "./environment-materials.js";
import { canopyLeaves } from "./environment-foliage.js";

// Rendering follows catalogue geometry and accepted server environment frames.
// Cutaways affect visibility only: hidden walls remain physical on the server.

const UP = new THREE.Vector3(0, 1, 0);

function scenery(mesh) {
  mesh.raycast = () => {};
  mesh.castShadow = true; mesh.receiveShadow = true;
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
  const geometry = new THREE.PlaneGeometry(extent, extent);
  geometry.rotateX(-Math.PI / 2);
  const mesh = new THREE.Mesh(materialUVs(geometry, "vinyl"), environmentMaterial("vinyl", 0xb9b7ac));
  mesh.quaternion.setFromUnitVectors(UP, up);
  mesh.position.copy(up).multiplyScalar(ground.level_m);
  return scenery(mesh);
}

// A surround is one indexed mesh with a colour per face and a flat list of
// parts. Faces are emitted part by part, which is checked here rather than
// assumed: if the blocks do not line up the whole enclosure is drawn and
// nothing is cut away, so a reader never silently loses a wall.
function surroundMesh(data, cutaway) {
  const parts = data.parts || [];
  const faces = data.indices.length / 3;
  const blocked = parts.length > 0 && parts.every(part => part.primitive !== "tiled_plane") && faces % parts.length === 0 &&
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
  const materials = parts.map(part => environmentMaterial(part.colour || part.name, 0xffffff, {vertexColors:true, doubleSided:true}));
  if (!materials.length) materials.push(environmentMaterial("paint", 0xffffff, {vertexColors:true, doubleSided:true}));
  const faceParts = new Uint16Array(faces);
  if (blocked) {
    for (let f = 0; f < faces; f++) faceParts[f] = Math.floor(f / perFace);
  } else if (parts.filter(p => p.primitive === "tiled_plane").length === 1 && parts.every(p => ["box", "tiled_plane"].includes(p.primitive))) {
    let offset = 0;
    for (let i = 0; i < parts.length; i++) {
      const count = parts[i].primitive === "box" ? 12 : faces - 12 * (parts.length - 1);
      faceParts.fill(i, offset, offset + count); offset += count;
    }
  }
  // One draw per material run, not one draw per triangle (outdoor ground has thousands).
  let start = 0, materialIndex = keep.length ? faceParts[keep[0]] : 0;
  for (let n = 1; n <= keep.length; n++) {
    const next = n === keep.length ? -1 : faceParts[keep[n]];
    if (next !== materialIndex) {
      geometry.addGroup(start * 3, (n - start) * 3, materialIndex);
      start = n; materialIndex = next;
    }
  }
  materialUVs(geometry, parts.some(p => /wood/.test(p.colour)) ? "wood" : "paint");
  const mesh = new THREE.Mesh(geometry, materials);
  mesh.name = data.id;
  mesh.userData.cutaway = blocked ? [...cut] : [];
  mesh.userData.cutawayApplied = blocked;
  return scenery(mesh);
}

function objectMesh(data, colour, offset) {
  const root = new THREE.Group(); root.name = data.id;
  // Beveled edges preserve the authored extents and produce real highlight rolloff.
  for (const part of data.parts || []) {
    if (part.name.includes("canopy")) { root.add(canopyLeaves(part)); continue; }
    let geometry, center;
    if (part.primitive === "box") {
      const low = new THREE.Vector3(...part.min_m), high = new THREE.Vector3(...part.max_m);
      const size = high.clone().sub(low); center = high.add(low).multiplyScalar(.5);
      geometry = new RoundedBoxGeometry(size.x, size.y, size.z, 3, Math.min(.012, size.x / 5, size.y / 5, size.z / 5));
    } else if (part.primitive === "cylinder") {
      geometry = new THREE.CylinderGeometry(part.radius_m, part.radius_m, part.height_m, 32);
      if (part.axis === 2) geometry.rotateX(Math.PI / 2);
      if (part.axis === 0) geometry.rotateZ(Math.PI / 2);
      center = new THREE.Vector3(...part.centre_m);
    } else {
      geometry = new THREE.SphereGeometry(part.radius_m, 40, 28);
      center = new THREE.Vector3(...part.centre_m);
    }
    geometry.translate(...center.toArray());
    const name = part.colour || (part.name.includes("canopy") ? "canopy" : data.id);
    const mesh = new THREE.Mesh(materialUVs(geometry, name), environmentMaterial(name, part.name === "trunk" ? 0x67513b : colour));
    root.add(scenery(mesh));
  }
  if (!root.children.length) {
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.Float32BufferAttribute(data.positions, 3));
    geometry.setIndex(data.indices); geometry.computeVertexNormals();
    root.add(scenery(new THREE.Mesh(materialUVs(geometry, data.id), environmentMaterial(data.id, colour))));
  }
  if (offset?.length === 3) root.position.fromArray(offset);
  return root;
}

export function mountSurround(scene, { api }) {
  const root = new THREE.Group();
  root.name = "surround";
  scene.add(root);
  const geometryCache = new Map();
  let generation = 0, latestState = null;
  const instancesById = new Map();

  const load = (url) => {
    if (!geometryCache.has(url)) geometryCache.set(url, api(url));
    return geometryCache.get(url);
  };
  function empty() {
    root.traverse((node) => {
      if (node.geometry) node.geometry.dispose();
      if (node.material) for (const material of Array.isArray(node.material) ? node.material : [node.material]) material.dispose();
    });
    root.clear(); instancesById.clear();
  }

  function update(state) {
    latestState = state;
    if (!state) return;
    for (const item of state.objects || []) {
      const instance = instancesById.get(item.id);
      if (!instance) continue;
      if (item.kind === "rigid") {
        const matrix = new THREE.Matrix4().set(
          ...item.rotation_matrix[0], 0, ...item.rotation_matrix[1], 0,
          ...item.rotation_matrix[2], 0, 0, 0, 0, 1);
        instance.quaternion.setFromRotationMatrix(matrix);
        const origin = new THREE.Vector3(...item.origin_m).applyQuaternion(instance.quaternion);
        instance.position.fromArray(item.position_m).sub(origin);
      } else {
        let mesh = instance.userData.deformed;
        if (!mesh) {
          instance.traverse(node => {
            if (node.geometry) node.geometry.dispose();
            if (node.material) node.material.dispose();
          });
          instance.clear(); instance.position.set(0, 0, 0);
          const geometry = new THREE.BufferGeometry();
          geometry.setAttribute("position", new THREE.Float32BufferAttribute(item.positions, 3));
          geometry.setIndex(item.indices); geometry.computeVertexNormals();
          const name = item.kind === "cloth" ? "blanket" : "pillow";
          mesh = scenery(new THREE.Mesh(materialUVs(geometry, name), environmentMaterial(name, item.kind === "cloth" ? 0x627e99 : 0xe9e3d7, {doubleSided:true})));
          instance.add(mesh); instance.userData.deformed = mesh;
        }
        mesh.geometry.getAttribute("position").array.set(item.positions);
        mesh.geometry.getAttribute("position").needsUpdate = true;
        mesh.geometry.computeVertexNormals(); mesh.geometry.computeBoundingSphere();
      }
    }
  }
  return {
    update,
    group: root,
    // `entry` is the selected scene, or the environment when no scene is chosen;
    // `environment` is always the base environment record, because the sky axis
    // and the engine plane belong to it and a scene never changes them.
    async apply(entry, environment, instances = [], objectRecords = []) {
      const current = ++generation;
      empty(); latestState = null;
      const world = entry?.world || environment?.world;
      if (!world) return { drawn: 0, note: "This environment declares no world." };
      const up = upAxis(environment);
      scene.traverse(node => {
        if (node.isDirectionalLight && node.castShadow) {
          node.position.copy(new THREE.Vector3(3, 4, 5).applyQuaternion(new THREE.Quaternion().setFromUnitVectors(UP, up)));
        }
      });
      const extent = world.enclosure_dimensions_m
        ? Math.max(world.enclosure_dimensions_m.width, world.enclosure_dimensions_m.depth) * 1.2
        : (world.ground?.extent_m || 12);
      const sky = skyMesh(world.sky, up, 5);
      if (sky) root.add(sky);
      const grid = world.surround_url ? null : groundGrid(world.ground, up, Math.max(6, extent));
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
      // Stable instance identities match the server, including repeated inserts.
      const records = new Map(objectRecords.map((o) => [o.id, o]));
      const placements = [
        ...(entry?.placements || (environment?.id === "bed" ? ["bed-frame", "bed-mattress", "pillow", "blanket"].map(object => ({object})) : [])).map((p) => ({ id: p.object, offset: p.offset_m })),
        ...instances.map((id) => ({ id, offset: [0, 0, 0] })),
      ];
      const counts = new Map();
      for (const p of placements) { counts.set(p.id, (counts.get(p.id) || 0) + 1); p.instance = p.id + "-" + counts.get(p.id); }
      await Promise.all(placements.map(async (placement) => {
        const record = records.get(placement.id);
        if (!record?.geometry_url) { notes.push("No geometry for " + placement.id); return; }
        try {
          const data = await load(record.geometry_url);
          if (current !== generation) return;
          const mesh = objectMesh(data, record.colour_rgb, placement.offset);
          if (record.initial_mesh) {
            mesh.traverse(node => { node.geometry?.dispose(); node.material?.dispose(); }); mesh.clear();
            const geometry = new THREE.BufferGeometry();
            geometry.setAttribute("position", new THREE.Float32BufferAttribute(record.initial_mesh.positions, 3));
            geometry.setIndex(record.initial_mesh.indices); geometry.computeVertexNormals();
            mesh.add(scenery(new THREE.Mesh(materialUVs(geometry, placement.id), environmentMaterial(placement.id, record.colour_rgb, {doubleSided:true}))));
          }
          mesh.userData.environmentInstance = placement.instance;
          instancesById.set(placement.instance, mesh); root.add(mesh);
          drawn++;
        } catch (error) { notes.push(placement.id + ": " + error.message); }
      }));
      if (current !== generation) return { drawn, note: "" };
      if (latestState) update(latestState);
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
