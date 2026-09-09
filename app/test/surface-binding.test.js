import {test} from 'node:test';
import assert from 'node:assert/strict';
import * as THREE from 'three';
import {bindSurfaceVertices,poseSurfaceVertices,updateSegmentSurface,surfacePickBone} from '../src/surface-binding.js';
import {WardrobeView} from '../src/clothing.js';
const I=[[1,0,0],[0,1,0],[0,0,1]];
const binding={schema:'ihm.segment-surface-binding.v1',coordinate_frame:'canonical_current_world',rule:'nearest_named_bone_envelope',reference_coordinate_frame:'canonical_rest',tie_break:'lexicographic_native_segment_id',binding_identity:'fixture-v1',surface_entity_ids:['skin'],segments:[{id:'right',bone_id:'bone-r',bounds_min_m:[0,-1,-1],bounds_max_m:[2,1,1],reference_centroid_m:[1,0,0]},{id:'left',bone_id:'bone-l',bounds_min_m:[-2,-1,-1],bounds_max_m:[0,1,1],reference_centroid_m:[-1,0,0]}]};
const transforms={left:{centroid_m:[-1,2,0],rotation_matrix:I},right:{centroid_m:[1,0,3],rotation_matrix:[[0,-1,0],[1,0,0],[0,0,1]]}};
const frame={mechanics:{surface_binding:binding,surface_transforms:transforms},entities:{skin:{centroid_m:[90,0,0],translation_m:[90,0,0]}}};
const reference=new Float32Array([-1,0,0,1,1,0,0,0,0]);
function mesh(){const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.BufferAttribute(reference.slice(),3));g.setIndex([0,1,2]);return new THREE.Mesh(g,new THREE.MeshBasicMaterial());}
test('hard envelope support matches analytic segment motion, ties by native name and never accumulates',()=>{
 const bound=bindSurfaceVertices(reference,binding),out=new Float32Array(9);
 assert.deepEqual([...bound.owners],[0,1,0]);poseSurfaceVertices(reference,bound,transforms,out);
 assert.deepEqual([...out],[-1,2,0,0,0,3,0,2,0]);
 poseSurfaceVertices(reference,bound,transforms,out);assert.deepEqual([...out],[-1,2,0,0,0,3,0,2,0]);
 assert.throws(()=>poseSurfaceVertices(reference,bound,{left:transforms.left},out),/Missing/);
});
test('skin uses canonical per-vertex transforms with identity mesh matrix and restores rest on reset',()=>{
 const m=mesh();m.matrix.makeTranslation(90,0,0);assert.equal(updateSegmentSurface(m,frame,'skin'),true);
 assert.deepEqual([...m.geometry.attributes.position.array],[-1,2,0,0,0,3,0,2,0]);assert.equal(m.matrix.equals(new THREE.Matrix4()),true);
 assert.equal(surfacePickBone(m,{a:0,b:1,c:2},[0,0,3]),'bone-r');
 assert.equal(updateSegmentSurface(m,null,'skin'),false);assert.deepEqual([...m.geometry.attributes.position.array],[...reference]);
});
test('wardrobe follows the same posed support, independent owned cloth wins without a second transform',async()=>{
 const view=new WardrobeView(new THREE.Group(),{fetchGeometry:async()=>({positions:[...reference],indices:[0,1,2]})});
 view.bind({id:'skin'});view.setCatalog({garments:[{id:'shirt',geometry_url:'/shirt'}]});await view.setActive(['shirt']);
 view.update(frame,{skin:[0,0,0]});const m=view.meshes.get('shirt');assert.deepEqual([...m.geometry.attributes.position.array],[-1,2,0,0,0,3,0,2,0]);assert.equal(m.matrix.equals(new THREE.Matrix4()),true);
 view.update({...frame,mechanics:{...frame.mechanics,garment_mechanics:{garments:{shirt:{positions_m:[[3,0,0],[3,1,0],[3,0,1]]}}}}},{skin:[0,0,0]});
 assert.deepEqual([...m.geometry.attributes.position.array],[3,0,0,3,1,0,3,0,1]);assert.equal(m.matrix.equals(new THREE.Matrix4()),true);view.dispose();
});

test('continuous owned skin waits for verified sidecar, then renders the shared blended material support',async()=>{
 const {surfaceAssets}=await import('../src/surface-assets.js');const {createHash}=await import('node:crypto');const {surfacePickAnchor}=await import('../src/surface-binding.js');const {continuousMaterialPorts}=await import('../src/continuous-surface.js');
 const payload={schema:'ihm.continuous-surface-binding.v1',binding_identity:'continuous-fixture',segments:binding.segments,reference_positions_m:[...reference].reduce((a,v,i)=>{if(i%3===0)a.push([]);a.at(-1).push(v);return a;},[]),weights:[[0,1],[1,0],[.5,.5]]};
 const raw=JSON.stringify(payload),digest=createHash('sha256').update(raw).digest('hex'),meta={...binding,schema:payload.schema,binding_identity:payload.binding_identity,rule:'graph_regularized_linear_blend',weights_sha256:digest,weights_url:'/api/surface-binding/'+digest};
 const previous=surfaceAssets.fetcher;surfaceAssets.fetcher=async()=>new Response(raw);
 try{
  const f={mechanics:{surface_binding:meta,surface_transforms:transforms}},m=mesh();
  assert.equal(updateSegmentSurface(m,f,'skin'),true);assert.equal(m.visible,false);
  await surfaceAssets.ensure(meta);assert.equal(updateSegmentSurface(m,f,'skin'),true);assert.equal(m.visible,true);
  assert.deepEqual([...m.geometry.attributes.position.array],[-1,2,0,0,0,3,.5,.5,1.5]);assert.equal(m.matrix.equals(new THREE.Matrix4()),true);
  const anchor=surfacePickAnchor(m,{a:0,b:1,c:2},[.5,.5,1.5]),ports=continuousMaterialPorts(anchor,transforms,[0,0,2]);
  assert.deepEqual(ports.point_m,[.5,.5,1.5]);assert.equal(ports.ports.length,2);assert.deepEqual(ports.ports.map(p=>p.force_n),[[0,0,1],[0,0,1]]);
 }finally{surfaceAssets.fetcher=previous;surfaceAssets.records.delete(digest);}
});
