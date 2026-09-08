import * as THREE from 'three';
import { fbm } from './procedural-noise.js';

// Deterministic, metre-scaled albedo and relief. No screen-space decoration:
// grain, weave and paint stay attached to each object's material coordinates.
const textures = new Map();
export function materialFamily(name = '') {
  if (/window-pane/.test(name)) return 'glass';
  if (/rail|iv-stand|metal/.test(name)) return 'metal';
  if (/blanket|pillow|mattress|carpet|fabric|cloth|garment|shirt/.test(name)) return 'fabric';
  if (/wood|frame|table|chair|nightstand|trunk|tree/.test(name)) return 'wood';
  if (/grass|canopy/.test(name)) return 'grass';
  if (/patio|vinyl/.test(name)) return 'stone';
  if (/ball/.test(name)) return 'rubber';
  return 'paint';
}
function textureSet(family) {
  if (textures.has(family)) return textures.get(family);
  const size = 256, rgb = new Uint8Array(size * size * 4), relief = new Uint8Array(size * size * 4);
  // Lattice noise, not one random number per texel. A floor runs twenty metres
  // to the far wall, and white noise minified that far has no mipmap to fall
  // back on: it crawls, and a camera that moves a few millimetres repaints the
  // whole floor. See procedural-noise.js.
  const field = fbm([8, 16, 32], 1729 + family.length * 131);
  for (let y = 0; y < size; y++) for (let x = 0; x < size; x++) {
    const u = x / size, v = y / size, noise = field(u, v);
    let value = .94 + .06 * noise, height = noise * .25;
    if (family === 'wood') {
      const grain = Math.sin(2 * Math.PI * (u * 38 + .45 * Math.sin(v * Math.PI * 2) + .18 * Math.sin(v * 12 * Math.PI)));
      const fine = Math.sin(u * 420 * Math.PI + Math.sin(v * Math.PI * 4));
      const seam = x < 2 || (y < 2 && x > size / 2);
      value = seam ? .70 : .88 + .055 * grain + .018 * fine + .03 * noise;
      height = seam ? 0 : .55 + .1 * grain;
    } else if (family === 'fabric') {
      const warp = Math.pow(Math.sin(x * Math.PI / 4), 2), weft = Math.pow(Math.sin(y * Math.PI / 4), 2);
      const over = (Math.floor(x / 4) + Math.floor(y / 4)) % 2;
      height = .3 + .5 * (over ? warp : weft);
      value = .82 + .14 * height + .05 * noise;
    } else if (family === 'stone') {
      const joint = x < 3 || y < 3;
      value = joint ? .56 : .89 + .09 * noise; height = joint ? 0 : .6 + .05 * noise;
    } else if (family === 'grass') {
      // Blades in clumps rather than static: a fine directional streak over the
      // patchiness of the lawn beneath it.
      const blade = .5 + .5 * Math.sin((v * 64 + 3 * noise) * Math.PI * 2);
      value = .6 + .28 * noise + .12 * blade; height = .35 * noise + .6 * blade;
    } else if (family === 'paint') {
      value = .95 + .04 * noise + .01 * Math.sin(x * .08) * Math.sin(y * .11); height = noise * .4;
    }
    const i = 4 * (y * size + x);
    for (let c = 0; c < 3; c++) { rgb[i + c] = Math.round(255 * value); relief[i + c] = Math.round(255 * height); }
    rgb[i + 3] = relief[i + 3] = 255;
  }
  const make = (data, color) => {
    const texture = new THREE.DataTexture(data, size, size, THREE.RGBAFormat);
    texture.wrapS = texture.wrapT = THREE.RepeatWrapping;
    texture.magFilter = THREE.LinearFilter; texture.minFilter = THREE.LinearMipmapLinearFilter;
    texture.generateMipmaps = true; texture.anisotropy = 8;
    if (color) texture.colorSpace = THREE.SRGBColorSpace;
    texture.needsUpdate = true; return texture;
  };
  const result = {map: make(rgb, true), bumpMap: make(relief, false)};
  textures.set(family, result); return result;
}
// `doubleSided` is asked for, never assumed. A wall, a cut-open enclosure and a
// blanket are open sheets and have to be drawn from both sides; a bed frame, a
// ball and a table are closed solids, and drawing their interiors costs twice
// the fill, lights the back faces as though they were front ones, and writes
// the far side of the box into the shadow map, which is where the crawling
// speckle on the floor came from.
export function environmentMaterial(name, color, {vertexColors = false, doubleSided = false} = {}) {
  const family = materialFamily(name);
  const material = new THREE.MeshPhysicalMaterial({
    color: Array.isArray(color) ? new THREE.Color(...color) : color ?? 0xffffff,
    vertexColors, side: doubleSided ? THREE.DoubleSide : THREE.FrontSide,
    roughness: family === 'metal' ? .28 : family === 'wood' ? .48 : .88,
    metalness: family === 'metal' ? .85 : 0,
    clearcoat: family === 'wood' ? .22 : 0, clearcoatRoughness: .45,
    sheen: family === 'fabric' ? .65 : 0, sheenRoughness: .85,
    sheenColor: new THREE.Color(0xc9c1ad),
    ...textureSet(family), bumpScale: family === 'fabric' ? .0012 : family === 'paint' ? .0004 : .0007,
  });
  if (family === 'glass') {
    material.transparent = true; material.opacity = .28; material.roughness = .12;
    material.depthWrite = false; material.bumpScale = 0;
  }
  material.userData.environmentFamily = family;
  return material;
}

export function materialUVs(geometry, name) {
  const p = geometry.getAttribute('position'), n = geometry.getAttribute('normal');
  const uv = new Float32Array(p.count * 2);
  const family = materialFamily(name), scale = family === 'fabric' ? 6 : family === 'wood' ? 2 : 1.8;
  for (let i = 0; i < p.count; i++) {
    const ax = Math.abs(n.getX(i)), ay = Math.abs(n.getY(i)), az = Math.abs(n.getZ(i));
    const a = ax > ay && ax > az ? p.getZ(i) : p.getX(i);
    const b = az >= ax && az >= ay ? p.getY(i) : ay >= ax ? p.getZ(i) : p.getY(i);
    uv[2 * i] = a * scale; uv[2 * i + 1] = b * scale;
  }
  geometry.setAttribute('uv', new THREE.BufferAttribute(uv, 2));
  return geometry;
}
