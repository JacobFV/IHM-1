import {test} from 'node:test';
import assert from 'node:assert/strict';
import * as THREE from 'three';
import {mountSurround} from '../src/surround.js';
import {canonicalEnvironmentObject,environmentWorldMatrix} from '../src/world-frame.js';
const state={object_coordinate_frame:'gravity_aligned_world',world_frame:{world_to_canonical:[[0,-1,0,2],[1,0,0,3],[0,0,1,4],[0,0,0,1]]},objects:[]};
test('world frame renders static and dynamic geometry once and resets preview',async()=>{
 const scene=new THREE.Scene(),surround=mountSurround(scene,{api:async()=>({id:'ball',parts:[{name:'ball',primitive:'sphere',radius_m:.1,centre_m:[1,0,0]}]})});
 await surround.apply({world:{},placements:[{object:'ball'}]}, {},[],[{id:'ball',geometry_url:'/ball'}]);
 surround.update({...state,objects:[{id:'ball-1',kind:'rigid',origin_m:[1,0,0],position_m:[5,0,0],rotation_matrix:[[1,0,0],[0,1,0],[0,0,1]]}]});
 const object=surround.group.children[0];surround.group.updateMatrixWorld(true);
 assert.deepEqual(object.localToWorld(new THREE.Vector3(1,0,0)).toArray().map(v=>Math.round(v*1e10)/1e10),[2,8,4]);
 assert.deepEqual(surround.group.localToWorld(new THREE.Vector3(1,0,0)).toArray().map(v=>Math.round(v*1e10)/1e10),[2,4,4]);
 surround.update(null);assert.deepEqual(surround.group.position.toArray(),[0,0,0]);surround.dispose();
});
test('canonical picking transforms rigid orientation and every cloth node without modifying material origin',()=>{
 const raw={position_m:[1,2,3],origin_m:[7,8,9],rotation_matrix:[[1,0,0],[0,1,0],[0,0,1]],positions:[1,2,3,0,0,0],fixed_nodes:[1]};
 const result=canonicalEnvironmentObject(state,raw);
 assert.deepEqual(result.position_m,[0,4,7]);assert.deepEqual(result.positions,[0,4,7,2,3,4]);assert.deepEqual(result.rotation_matrix,[[0,-1,0],[1,0,0],[0,0,1]]);assert.deepEqual(result.origin_m,[7,8,9]);assert.deepEqual(raw.position_m,[1,2,3]);
 assert.throws(()=>environmentWorldMatrix({object_coordinate_frame:'gravity_aligned_world'}),/transform/);
 assert.throws(()=>environmentWorldMatrix({object_coordinate_frame:'other'}),/Unknown/);
});
