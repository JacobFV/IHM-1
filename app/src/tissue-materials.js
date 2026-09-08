import * as THREE from "three";
import { fbm } from "./procedural-noise.js";

// Relief and grain for tissue, at the scale the tissue actually has it: dermal
// furrows around a millimetre, muscle fascicles around two, trabecular mottle
// around a centimetre. Everything here is deterministic and generated from a
// fixed seed, so two runs paint the same body.
//
// The one rule the rest of this file obeys: a texture modulates luminance and
// never hue. The palette lane measured, transferred or synthesized a colour per
// structure and says which of the three each one is; a texture that tinted it
// would quietly invalidate that audit. So every map here is grey, centred near
// white, and multiplies the structure's own colour.

// A tissue family is what the surface is made of, which is not the same as the
// system that owns it: an artery and a vein are both vessel wall, an oesophagus
// and a bronchus are both mucosa. The family is read from the system first,
// because that is the field every structure record carries, and refined by name
// only where one system holds two materials.
export const TISSUE_FAMILIES = ["skin", "muscle", "bone", "cartilage", "vessel", "mucosa", "nerve", "fibrous", "hair"];
const BY_SYSTEM = {
  integumentary: "skin", muscular: "muscle", skeletal: "bone", hair: "hair",
  arterial: "vessel", venous: "vessel", microvascular: "vessel", cardiac: "muscle",
  lymphatic: "vessel", nervous: "nerve", sensory: "nerve",
  digestive: "mucosa", respiratory: "mucosa", urinary: "mucosa",
  reproductive: "mucosa", endocrine: "mucosa", connective: "fibrous",
};
export function tissueFamily(structure = {}) {
  const name = String(structure.name || structure.id || "").toLowerCase();
  const system = String(structure.system || "").toLowerCase();
  // Cartilage and the fibrous sheets are named, not systematised: they sit in
  // the skeletal and connective lists beside bone and beside fat.
  if (/cartilage|meniscus|disc of|intervertebral/.test(name)) return "cartilage";
  if (system === "skeletal" && /ligament|tendon|membrane|capsule|fascia/.test(name)) return "fibrous";
  if (system === "muscular" && /tendon|aponeurosis|fascia/.test(name)) return "fibrous";
  return BY_SYSTEM[system] || "fibrous";
}

// Metres per texture tile. Set from the feature, not from taste: dermal furrows
// are sub-millimetre so 256 texels cover 3 cm, trabecular mottle is coarse so
// the same 256 cover 12 cm.
const TILE_M = { skin: 0.03, muscle: 0.05, bone: 0.12, cartilage: 0.09, vessel: 0.04, mucosa: 0.06, nerve: 0.03, fibrous: 0.05, hair: 0.02 };

// Surface response per family, in the terms MeshPhysicalMaterial takes. The
// wet surfaces -- serosa, mucosa, a fresh muscle belly, sebum on skin -- are
// the ones that get a clearcoat, because a thin film over a rough substrate is
// two lobes and one roughness cannot be both. Dry periosteum and dry collagen
// get none, which is not only right but is what keeps the extra lobe out of the
// shader for most of the seven thousand structures in this body.
const FINISH = {
  skin: { roughness: 0.58, clearcoat: 0.18, clearcoatRoughness: 0.6, bumpScale: 0.0007 },
  muscle: { roughness: 0.42, clearcoat: 0.3, clearcoatRoughness: 0.35, bumpScale: 0.0006 },
  bone: { roughness: 0.66, clearcoat: 0, clearcoatRoughness: 0.5, bumpScale: 0.0005 },
  cartilage: { roughness: 0.34, clearcoat: 0.35, clearcoatRoughness: 0.25, bumpScale: 0.00025 },
  vessel: { roughness: 0.4, clearcoat: 0.28, clearcoatRoughness: 0.3, bumpScale: 0.0003 },
  mucosa: { roughness: 0.3, clearcoat: 0.42, clearcoatRoughness: 0.22, bumpScale: 0.0004 },
  nerve: { roughness: 0.46, clearcoat: 0, clearcoatRoughness: 0.4, bumpScale: 0.0003 },
  fibrous: { roughness: 0.62, clearcoat: 0, clearcoatRoughness: 0.55, bumpScale: 0.00045 },
  hair: { roughness: 0.4, clearcoat: 0.3, clearcoatRoughness: 0.3, bumpScale: 0.0002 },
};

const SIZE = 256;
const cache = new Map();
// `albedo` returns luminance about 1 and `relief` returns height in 0..1; the
// pair is what distinguishes one tissue from another, so each family is written
// as those two functions of the tile coordinate and nothing else.
const PATTERN = {
  skin(u, v, n) {
    // Dermal furrows: two crossed families of shallow creases at about 60
    // degrees, over the coarse mottle of vascularity beneath. The phase of each
    // family is pushed around by more than a whole cycle of noise, because a
    // furrow field that stays in step reads as a printed grid rather than as
    // skin.
    const a = Math.sin((u * 34 + v * 19 + 5 * n(u, v)) * Math.PI * 2);
    const b = Math.sin((u * 21 - v * 31 + 5 * n(v, u)) * Math.PI * 2);
    const furrow = Math.min(1, 0.5 * (Math.abs(a) + Math.abs(b)));
    const pore = n(u * 4, v * 4);
    return { albedo: 0.88 + 0.09 * furrow + 0.06 * n(u, v) - 0.05 * (pore > 0.84 ? 1 : 0),
             relief: 0.3 + 0.6 * furrow - 0.3 * (pore > 0.84 ? 1 : 0) };
  },
  muscle(u, v, n) {
    // Fascicles run along the tile's u axis: bundles a few millimetres across,
    // with the fibre striation inside them an order finer.
    const bundle = Math.abs(Math.sin((v * 7 + 0.35 * n(u, v)) * Math.PI));
    const fibre = 0.5 + 0.5 * Math.sin((v * 46 + 1.6 * n(u, v)) * Math.PI * 2);
    return { albedo: 0.86 + 0.12 * bundle + 0.05 * fibre * bundle,
             relief: 0.2 + 0.6 * bundle + 0.2 * fibre };
  },
  bone(u, v, n) {
    // Periosteal surface: broad mottle with the vascular foramina picked out as
    // small dark pits rather than drawn as a pattern.
    const mottle = n(u, v), fine = n(u * 4, v * 4);
    const pit = fine > 0.9 ? 1 : 0;
    return { albedo: 0.9 + 0.09 * mottle - 0.14 * pit, relief: 0.5 + 0.35 * mottle - 0.5 * pit };
  },
  cartilage(u, v, n) {
    return { albedo: 0.95 + 0.05 * n(u, v), relief: 0.45 + 0.2 * n(u * 2, v * 2) };
  },
  vessel(u, v, n) {
    // Adventitia: fine longitudinal fibre, and no cross structure, because that
    // is what reads as a vessel wall rather than as a tube.
    // Constant along u and varying across it, because u is the vessel's own
    // long axis: a pattern that varied along u would ring the tube instead.
    const fibre = 0.5 + 0.5 * Math.sin((v * 58 + 1.2 * n(v, u)) * Math.PI * 2);
    return { albedo: 0.92 + 0.06 * fibre + 0.03 * n(u, v), relief: 0.4 + 0.35 * fibre };
  },
  mucosa(u, v, n) {
    // Rugae: a coarse folded relief, wet, with no fine detail to speak of.
    const fold = Math.abs(Math.sin((u * 9 + 2.8 * n(u, v)) * Math.PI));
    return { albedo: 0.9 + 0.09 * fold + 0.04 * n(u * 3, v * 3), relief: 0.25 + 0.7 * fold };
  },
  nerve(u, v, n) {
    const fascicle = 0.5 + 0.5 * Math.sin((v * 26 + 0.9 * n(u, v)) * Math.PI * 2);
    return { albedo: 0.93 + 0.06 * fascicle, relief: 0.4 + 0.4 * fascicle };
  },
  fibrous(u, v, n) {
    // Crossed collagen, the same construction as a weave but far less regular.
    const warp = 0.5 + 0.5 * Math.sin((u * 30 + 2 * n(u, v)) * Math.PI * 2);
    const weft = 0.5 + 0.5 * Math.sin((v * 24 + 2 * n(v, u)) * Math.PI * 2);
    return { albedo: 0.88 + 0.1 * Math.max(warp, weft), relief: 0.25 + 0.6 * Math.max(warp, weft) };
  },
  hair(u, v, n) {
    const strand = 0.5 + 0.5 * Math.sin((v * 90 + 1.4 * n(v, u)) * Math.PI * 2);
    return { albedo: 0.85 + 0.15 * strand, relief: 0.3 + 0.6 * strand };
  },
};

function textureSet(family) {
  if (cache.has(family)) return cache.get(family);
  const pattern = PATTERN[family] || PATTERN.fibrous;
  const n = fbm([8, 16, 32], 1729 + family.length * 131);
  const rgb = new Uint8Array(SIZE * SIZE * 4), relief = new Uint8Array(SIZE * SIZE * 4);
  for (let y = 0; y < SIZE; y++) for (let x = 0; x < SIZE; x++) {
    const { albedo, relief: height } = pattern(x / SIZE, y / SIZE, n);
    const i = 4 * (y * SIZE + x);
    const a = Math.max(0, Math.min(255, Math.round(255 * albedo)));
    const h = Math.max(0, Math.min(255, Math.round(255 * height)));
    rgb[i] = rgb[i + 1] = rgb[i + 2] = a;
    relief[i] = relief[i + 1] = relief[i + 2] = h;
    rgb[i + 3] = relief[i + 3] = 255;
  }
  const make = (data, colour) => {
    const texture = new THREE.DataTexture(data, SIZE, SIZE, THREE.RGBAFormat);
    texture.wrapS = texture.wrapT = THREE.RepeatWrapping;
    texture.magFilter = THREE.LinearFilter;
    texture.minFilter = THREE.LinearMipmapLinearFilter;
    texture.generateMipmaps = true;
    texture.anisotropy = 8;
    if (colour) texture.colorSpace = THREE.SRGBColorSpace;
    texture.needsUpdate = true;
    return texture;
  };
  // One set per family, shared by every structure in it: the body is thousands
  // of meshes and each one owning a 256-square pair would cost hundreds of
  // megabytes for nine distinct images.
  const set = { map: make(rgb, true), bumpMap: make(relief, false) };
  cache.set(family, set);
  return set;
}

// Families whose pattern has a direction. Fascicles run along a muscle belly,
// collagen along a ligament, adventitia along a vessel; skin, bone and mucosa
// have no such axis and are projected on the world frame.
const DIRECTIONAL = new Set(["muscle", "nerve", "vessel", "fibrous", "hair"]);

// The structure's own long axis, by power iteration on the vertex covariance.
// A muscle belly's fibres run along it to first order, which is the whole point:
// projecting every structure on the world frame banded the entire body in one
// direction at once, so the abdomen and the arm were striped the same way and
// neither was striped the way its fibres run.
function principalAxis(position) {
  const count = position.count;
  let mx = 0, my = 0, mz = 0;
  for (let i = 0; i < count; i++) { mx += position.getX(i); my += position.getY(i); mz += position.getZ(i); }
  mx /= count; my /= count; mz /= count;
  const c = [0, 0, 0, 0, 0, 0]; // xx, yy, zz, xy, xz, yz
  for (let i = 0; i < count; i++) {
    const x = position.getX(i) - mx, y = position.getY(i) - my, z = position.getZ(i) - mz;
    c[0] += x * x; c[1] += y * y; c[2] += z * z; c[3] += x * y; c[4] += x * z; c[5] += y * z;
  }
  let v = new THREE.Vector3(1, 0.7, 0.3);
  for (let step = 0; step < 12; step++) {
    const x = c[0] * v.x + c[3] * v.y + c[4] * v.z;
    const y = c[3] * v.x + c[1] * v.y + c[5] * v.z;
    const z = c[4] * v.x + c[5] * v.y + c[2] * v.z;
    const length = Math.hypot(x, y, z);
    if (!(length > 0)) return new THREE.Vector3(0, 1, 0);
    v = new THREE.Vector3(x / length, y / length, z / length);
  }
  return v;
}

// Body geometry arrives as positions and indices in metres, with no texture
// coordinates, so they are projected here. The tile lies in a frame the family
// chooses -- the structure's own long axis where the pattern has a direction,
// the world frame where it does not -- and which of the frame's two cross
// planes a vertex lands on is decided by its dominant normal, which keeps the
// tile square wherever the surface faces one way and stretches it only across
// the seams where the normal turns over.
export function tissueUVs(geometry, family) {
  const position = geometry.getAttribute("position");
  const normal = geometry.getAttribute("normal");
  if (!position || !normal || !position.count) return geometry;
  const scale = 1 / (TILE_M[family] || 0.05);
  const along = DIRECTIONAL.has(family) ? principalAxis(position) : new THREE.Vector3(1, 0, 0);
  // Two directions across the fibre, orthogonal to it and to each other.
  const across = new THREE.Vector3(0, 0, 0);
  const seed = Math.abs(along.x) < 0.9 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0);
  across.crossVectors(along, seed).normalize();
  const third = new THREE.Vector3().crossVectors(across, along);
  const uv = new Float32Array(position.count * 2);
  const p = new THREE.Vector3(), n = new THREE.Vector3();
  for (let i = 0; i < position.count; i++) {
    p.fromBufferAttribute(position, i);
    n.fromBufferAttribute(normal, i);
    // The fibre coordinate is the same everywhere on the structure; only the
    // coordinate across it switches planes, so a band never breaks mid-belly.
    uv[2 * i] = p.dot(along) * scale;
    uv[2 * i + 1] = (Math.abs(n.dot(across)) >= Math.abs(n.dot(third)) ? p.dot(third) : p.dot(across)) * scale;
  }
  geometry.setAttribute("uv", new THREE.BufferAttribute(uv, 2));
  return geometry;
}

// Opacity is a property of how the reader has set the system, not of the
// tissue, so it is passed in. A surface at full opacity is opaque and writes
// depth; below that it stops writing depth, because a translucent skin that
// writes depth hides whatever is drawn after it -- which, in a scene sorted
// back to front every frame, is a different set of organs each time the camera
// moves.
export function applyOpacity(material, opacity) {
  material.opacity = opacity;
  material.transparent = opacity < 1;
  material.depthWrite = opacity >= 1;
  return material;
}

export function tissueMaterial(structure, { color, opacity = 1 } = {}) {
  const family = tissueFamily(structure);
  const finish = FINISH[family] || FINISH.fibrous;
  const material = new THREE.MeshPhysicalMaterial({
    color: color ?? 0x9eb7b7,
    // Anatomical surfaces are closed shells whose winding the source atlas gets
    // right; drawing their back faces as well doubles the cost, and under a
    // translucent skin it is the back faces that read as a second body inside
    // the first.
    side: THREE.FrontSide,
    metalness: 0,
    roughness: finish.roughness,
    clearcoat: finish.clearcoat,
    clearcoatRoughness: finish.clearcoatRoughness,
    bumpScale: finish.bumpScale,
    ...textureSet(family),
  });
  material.userData.tissueFamily = family;
  return applyOpacity(material, opacity);
}
