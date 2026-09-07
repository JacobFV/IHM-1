import * as THREE from "three";
import model from "./gimbal-body.json";

// A view gizmo in the CAD sense, with this body instead of a cube. It lives in
// its own scene with its own fixed camera, so the main camera's position, zoom
// and clipping never reach it; only its orientation is copied, every frame, so
// the little body turns to show which way the real body is facing.
//
// The mesh is this body's own outer envelope, decimated offline by
// scripts/bake_gimbal_body.py. It carries orientation, never measurement.
export const PLANES = [
  // Plane geometry lies in XY with +Z normal; each rotation puts it in place.
  { id: "sagittal", rotation: [0, Math.PI / 2, 0], direction: [1, 0, 0], up: [0, 1, 0] },
  { id: "coronal", rotation: [0, 0, 0], direction: [0, 0, 1], up: [0, 1, 0] },
  { id: "transverse", rotation: [-Math.PI / 2, 0, 0], direction: [0, 1, 0], up: [0, 0, -1] },
];

export function mountGimbal(canvas, { camera, onSelect, size = 108 }) {
  const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(2, globalThis.devicePixelRatio || 1));
  renderer.setSize(size, size, false);
  renderer.setClearColor(0x000000, 0);

  const scene = new THREE.Scene();
  const view = new THREE.PerspectiveCamera(30, 1, 0.1, 20);
  view.position.set(0, 0, 3.1);
  scene.add(new THREE.HemisphereLight(0xdff2ee, 0x2b3a3d, 2.4));
  const key = new THREE.DirectionalLight(0xfff3e2, 2.2);
  key.position.set(1.4, 1.8, 2.4);
  scene.add(key);
  const root = new THREE.Group();
  scene.add(root);

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.Float32BufferAttribute(model.positions, 3));
  geometry.setIndex(model.indices);
  geometry.computeVertexNormals();
  const body = new THREE.Mesh(
    geometry,
    new THREE.MeshStandardMaterial({ color: 0xd8e8e2, roughness: 0.62, metalness: 0.04 }),
  );
  body.renderOrder = 0;
  root.add(body);

  // Planes are drawn after the body and never write depth, so the parts behind
  // the body are occluded and the silhouette stays readable through the rest.
  const quads = new Map();
  for (const plane of PLANES) {
    const group = new THREE.Group();
    group.rotation.set(...plane.rotation);
    const fill = new THREE.Mesh(
      new THREE.PlaneGeometry(0.78, 0.78),
      new THREE.MeshBasicMaterial({
        color: 0x86ccbb, transparent: true, opacity: 0.1,
        side: THREE.DoubleSide, depthWrite: false,
      }),
    );
    fill.renderOrder = 1;
    fill.userData.plane = plane.id;
    const edge = new THREE.LineSegments(
      new THREE.EdgesGeometry(fill.geometry),
      new THREE.LineBasicMaterial({ color: 0xbedeD6, transparent: true, opacity: 0.5, depthWrite: false }),
    );
    edge.renderOrder = 2;
    edge.userData.plane = plane.id;
    group.add(fill, edge);
    root.add(group);
    quads.set(plane.id, { fill, edge });
  }

  let hovered = null, selected = null, disposed = false;
  function paint() {
    for (const [id, quad] of quads) {
      const on = id === selected, hot = id === hovered;
      quad.fill.material.opacity = on ? 0.3 : hot ? 0.22 : 0.1;
      quad.fill.material.color.set(on || hot ? 0x9fe0cd : 0x86ccbb);
      quad.edge.material.opacity = on ? 1 : hot ? 0.85 : 0.5;
      quad.edge.material.color.set(on || hot ? 0x9fe0cd : 0xbeded6);
    }
    canvas.dataset.hover = hovered || "";
    canvas.dataset.selected = selected || "";
    canvas.style.cursor = hovered ? "pointer" : "";
  }
  paint();

  const ray = new THREE.Raycaster();
  function planeAt(event) {
    const rect = canvas.getBoundingClientRect();
    if (!rect.width || !rect.height) return null;
    ray.setFromCamera(
      new THREE.Vector2(
        ((event.clientX - rect.left) / rect.width) * 2 - 1,
        (-(event.clientY - rect.top) / rect.height) * 2 + 1,
      ),
      view,
    );
    // Faces first, then outlines: a plane seen edge-on has no face area left,
    // and its border is the only thing the reader can still aim at.
    ray.params.Line.threshold = 0.03;
    const targets = [...quads.values()].flatMap((q) => [q.fill, q.edge]);
    const hit = ray.intersectObjects(targets, false)[0];
    if (!hit) return null;
    // A plane sliver hidden behind the body is not a target the reader can see.
    const occluder = ray.intersectObject(body, false)[0];
    return occluder && occluder.distance < hit.distance - 1e-3 ? null : hit.object.userData.plane;
  }
  canvas.addEventListener("pointermove", (event) => {
    const found = planeAt(event);
    if (found !== hovered) { hovered = found; paint(); render(); }
  });
  canvas.addEventListener("pointerleave", () => {
    if (hovered) { hovered = null; paint(); render(); }
  });
  canvas.addEventListener("click", (event) => {
    const found = planeAt(event);
    if (found) { selected = found; paint(); render(); onSelect(found); }
  });

  const orientation = new THREE.Quaternion();
  function render() { renderer.render(scene, view); }
  return {
    // The gizmo shows the scene as the main camera sees it: the inverse of the
    // camera's own rotation, at a fixed distance of its own.
    update() {
      if (disposed) return;
      const next = camera.quaternion.clone().invert();
      if (next.angleTo(orientation) < 1e-4) return;
      orientation.copy(next);
      root.quaternion.copy(next);
      render();
    },
    select(id) { selected = id; paint(); render(); },
    get hovered() { return hovered; },
    dispose() {
      disposed = true;
      geometry.dispose();
      body.material.dispose();
      for (const quad of quads.values()) {
        quad.fill.geometry.dispose(); quad.fill.material.dispose();
        quad.edge.geometry.dispose(); quad.edge.material.dispose();
      }
      renderer.dispose();
    },
  };
}
