import {test} from 'node:test';
import assert from 'node:assert/strict';
import * as THREE from 'three';
import fs from 'node:fs';
import zlib from 'node:zlib';
import {prepareGarmentContact,GarmentContactView} from '../src/garment-contact.js';

function fixture(){
  const p=[[0,0,0],[1,0,0],[0,1,0]],q=p.map(v=>[v[0],v[1],v[2]+.01]),zero=p.map(()=>[0,0,0]);
  return {schema_version:1,frame:{units:'m',axes:{x:'left',y:'superior',z:'anterior'}},time_s:[0,.005],
    tissue:{positions_m:[p,q],triangles:[[0,1,2]],triangle_material_index:[0],material_regions:[{source_id:'owner',name:'Tissue'}],source_volume_node_indices:[0,1,2],contact_force_n:[zero,[[0,0,-1],[0,0,0],[0,0,0]]]},
    panel:{positions_m:[p,q],triangles:[[0,2,1]],source_node_indices:[0,1,2],source_triangle_indices:[0],source_positions_m:p,source_triangles:[[0,1,2]],source_garment_id:'shorts',contact_force_n:[zero,[[0,0,1],[0,0,0],[0,0,0]]]},
    frames:[{time_s:.005,minimum_jacobian:.99,contact_dissipation_j:.001}],configuration:{},report:{limitations:['Local fixture']}};
}
const garment=()=>({id:'shorts',positions:[0,0,0,1,0,0,0,1,0,1,1,0],indices:[0,1,2,1,3,2]});
test('exact recorded positions and actual clock survive materialization; unrecorded initial metrics remain missing',()=>{
  const data=fixture(),prepared=prepareGarmentContact(data,garment(),new Set(['owner']));
  assert.deepEqual(prepared.frames.map(f=>f.time_s),[0,.005]);
  assert.deepEqual(prepared.remainder_indices,[1,3,2]);
  assert.equal(prepared.physiology.values['Panel contact resultant'][1],1);
  assert.equal(prepared.physiology.values['Minimum tissue Jacobian'][0],null);
  assert.equal(prepared.physiology.values['Maximum tissue displacement'][1],.01);
  assert.equal(data.tissue.positions_m[1][0][2],.01);
});
test('rejects mismatched garment topology, unrelated owners, nonfinite positions and inconsistent clocks',()=>{
  for(const mutate of [d=>d.panel.source_triangle_indices=[1],d=>d.tissue.material_regions[0].source_id='unrelated',d=>d.tissue.positions_m[1][0][0]=NaN,d=>d.time_s[1]=0,d=>d.frames[0].time_s=.004,d=>d.panel.positions_m.pop(),d=>d.frame.units='mm']){
    const data=fixture();mutate(data);assert.throws(()=>prepareGarmentContact(data,garment(),new Set(['owner'])));
  }
});
test('replaces only participating owners and exact panel faces, restores indices on close, follows clothes visibility',()=>{
  const bodyGroup=new THREE.Group(),objects=new Map(),reference=garment();
  const owner=new THREE.Mesh(new THREE.BufferGeometry(),new THREE.MeshBasicMaterial());objects.set('owner',owner);bodyGroup.add(owner);
  const unrelated=new THREE.Mesh(new THREE.BufferGeometry(),new THREE.MeshBasicMaterial());objects.set('bone',unrelated);bodyGroup.add(unrelated);
  const shorts=new THREE.Mesh(new THREE.BufferGeometry(),new THREE.MeshBasicMaterial());shorts.geometry.setIndex(reference.indices);bodyGroup.add(shorts);
  const clothingView={materializations:[reference],meshes:new Map([['shorts',shorts]])};
  const view=new GarmentContactView({bodyGroup,objects,clothingView,ownerVisible:()=>true});
  view.open(fixture());view.draw(1);
  assert.equal(bodyGroup.visible,true);assert.equal(owner.visible,false);assert.equal(unrelated.visible,true);
  assert.deepEqual(Array.from(shorts.geometry.index.array),[1,3,2]);
  assert.ok(Math.abs(view.tissueMeshes[0].geometry.attributes.position.array[2]-.01)<1e-9);
  shorts.visible=false;view.draw(1);assert.equal(view.panelMesh.visible,false);
  view.close();assert.equal(owner.visible,true);assert.deepEqual(Array.from(shorts.geometry.index.array),reference.indices);
  assert.equal(bodyGroup.children.length,3);
});

const held=new URL('../../data/derived/garment-tissue-display-v2/display.json.gz',import.meta.url);
test('held computed boundary export maps every panel face to its retained canonical garment',{skip:!fs.existsSync(held)},()=>{
  const data=JSON.parse(zlib.gunzipSync(fs.readFileSync(held)));
  const garment=JSON.parse(fs.readFileSync(new URL('../../data/derived/clothing/garments.json',import.meta.url))).garments.find(g=>g.id==='shorts');
  const prepared=prepareGarmentContact(data,garment,new Set(data.replaced_body_entity_ids));
  assert.equal(prepared.frames.length,49);assert.equal(data.tissue.positions_m[0].length,1327);
  assert.equal(prepared.remainder_indices.length,garment.indices.length-225*3);
  assert.equal(data.report.whole_garment_containment_validated,false);
  assert.ok(prepared.physiology.values['Panel contact resultant'][7]>.035);
});
