import test from 'node:test';
import assert from 'node:assert/strict';
import * as THREE from 'three';
import {tissueFamily,tissueMaterial,tissueUVs,applyOpacity,TISSUE_FAMILIES} from '../src/tissue-materials.js';
import {fbm,valueNoise} from '../src/procedural-noise.js';

test('tissue is read from the system, and from the name only where one system holds two materials',()=>{
 assert.equal(tissueFamily({system:'integumentary',name:'Skin'}),'skin');
 assert.equal(tissueFamily({system:'arterial',name:'Aorta'}),'vessel');
 assert.equal(tissueFamily({system:'venous',name:'Superior vena cava'}),'vessel');
 assert.equal(tissueFamily({system:'digestive',name:'Stomach'}),'mucosa');
 assert.equal(tissueFamily({system:'skeletal',name:'Femur'}),'bone');
 assert.equal(tissueFamily({system:'skeletal',name:'Costal cartilage of rib 3'}),'cartilage');
 assert.equal(tissueFamily({system:'skeletal',name:'Anterior longitudinal ligament'}),'fibrous');
 assert.equal(tissueFamily({system:'muscular',name:'Biceps brachii'}),'muscle');
 assert.equal(tissueFamily({system:'muscular',name:'Calcaneal tendon'}),'fibrous');
 // A system nobody anticipated still gets a surface rather than a crash.
 assert.ok(TISSUE_FAMILIES.includes(tissueFamily({system:'not-a-system',name:'x'})));
});

test('every family carries grain and relief, and no family tints the palette colour',()=>{
 for(const system of ['integumentary','muscular','skeletal','arterial','digestive','nervous','connective']){
  const material=tissueMaterial({system,name:system},{color:'#8899aa'});
  assert.ok(material.map,`${system} albedo`);
  assert.ok(material.bumpMap,`${system} relief`);
  assert.ok(material.bumpScale>0);
  assert.equal(`#${material.color.getHexString()}`,'#8899aa','the palette still decides the colour');
  // Grey maps only: a coloured albedo would silently restate the measured,
  // transferred and synthesized colours the palette lane audits.
  const data=material.map.image.data;
  for(let i=0;i<data.length;i+=4)
   assert.ok(data[i]===data[i+1]&&data[i+1]===data[i+2],`${system} albedo texel ${i/4} is not grey`);
 }
});

test('one texture set per family, shared by every structure in it',()=>{
 const a=tissueMaterial({system:'muscular',name:'Deltoid'});
 const b=tissueMaterial({system:'muscular',name:'Soleus'});
 assert.equal(a.map,b.map);
 assert.notEqual(a.map,tissueMaterial({system:'skeletal',name:'Tibia'}).map);
});

test('full opacity is opaque and writes depth; anything less does neither',()=>{
 const solid=tissueMaterial({system:'skeletal',name:'Femur'},{opacity:1});
 assert.equal(solid.transparent,false);
 assert.equal(solid.depthWrite,true);
 const skin=tissueMaterial({system:'integumentary',name:'Skin'},{opacity:.22});
 assert.equal(skin.transparent,true);
 assert.equal(skin.depthWrite,false,'a translucent skin that wrote depth would hide the organs behind it');
 applyOpacity(skin,1);
 assert.equal(skin.transparent,false);
 assert.equal(skin.depthWrite,true);
});

test('closed anatomical shells are drawn from the front only',()=>{
 assert.equal(tissueMaterial({system:'skeletal',name:'Femur'}).side,THREE.FrontSide);
});

test('projected coordinates put the tile at the family scale on any facing',()=>{
 const geometry=new THREE.BufferGeometry();
 geometry.setAttribute('position',new THREE.Float32BufferAttribute([0,0,0, .03,0,0, 0,.03,0],3));
 geometry.setAttribute('normal',new THREE.Float32BufferAttribute([0,0,1, 0,0,1, 0,0,1],3));
 tissueUVs(geometry,'skin');
 const uv=geometry.getAttribute('uv');
 // Skin tiles every 3 cm, so a 3 cm edge is exactly one tile.
 assert.ok(Math.abs(uv.getX(1)-uv.getX(0)-1)<1e-6);
 assert.ok(Math.abs(uv.getY(2)-uv.getY(0)-1)<1e-6);
});

test('noise is lattice noise: tileable, smooth and the same on every run',()=>{
 const field=valueNoise(8,1729);
 assert.ok(Math.abs(field(0,0)-field(1,1))<1e-12,'wraps, so a repeating texture has no seam');
 // Neighbouring texels of a 256-square differ by far less than white noise
 // would: this is the property that lets the mipmap chain average correctly
 // instead of crawling when the surface is minified.
 let worst=0;
 for(let i=0;i<256;i++)worst=Math.max(worst,Math.abs(field(i/256,.5)-field((i+1)/256,.5)));
 assert.ok(worst<.1,`neighbouring texels differ by ${worst}`);
 assert.equal(fbm([8,16,32],11)(.25,.75),fbm([8,16,32],11)(.25,.75));
});

test('a directional tissue is striped along its own long axis, not along the world',()=>{
 // A belly twice as long as it is wide, lying on the world diagonal.
 const axis=new THREE.Vector3(1,1,0).normalize();
 const positions=[],normals=[];
 for(const t of [-.1,-.05,0,.05,.1])for(const off of [-.01,.01]){
  positions.push(axis.x*t,axis.y*t+off,off);normals.push(0,0,1);
 }
 const geometry=new THREE.BufferGeometry();
 geometry.setAttribute('position',new THREE.Float32BufferAttribute(positions,3));
 geometry.setAttribute('normal',new THREE.Float32BufferAttribute(normals,3));
 tissueUVs(geometry,'muscle');
 const uv=geometry.getAttribute('uv');
 // The fibre coordinate follows the belly: the two ends are a full 0.2 m apart
 // along the axis, which at the muscle tile of 5 cm is four tiles.
 assert.ok(Math.abs(Math.abs(uv.getX(8)-uv.getX(0))-4)<.2,`fibre coordinate spans ${uv.getX(8)-uv.getX(0)}`);
 // And the coordinate across it barely moves between the two ends, because
 // across is across however the belly is turned in the world.
 assert.ok(Math.abs(uv.getY(8)-uv.getY(0))<.5);
});
