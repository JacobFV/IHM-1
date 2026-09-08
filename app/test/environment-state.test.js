import test from 'node:test';
import assert from 'node:assert/strict';
import * as THREE from 'three';
import {mountSurround} from '../src/surround.js';
import {environmentMaterial,materialFamily} from '../src/environment-materials.js';

const geometry={id:'blanket',positions:[0,0,0,1,0,0,0,1,0],indices:[0,1,2]};
test('accepted server cloth positions replace the preview; rendering never integrates physics',async()=>{
 const scene=new THREE.Scene(),view=mountSurround(scene,{api:async()=>geometry});
 const entry={world:{},placements:[{object:'blanket'}]};
 await view.apply(entry,{id:'studio'},[],[{id:'blanket',geometry_url:'/cloth',colour_rgb:[.3,.4,.5]}]);
 view.update({objects:[{id:'blanket-1',kind:'cloth',positions:[0,0,.1,1,0,.2,0,1,.3],indices:[0,1,2]}]});
 const instance=view.group.children[0],mesh=instance.userData.deformed;
 assert.ok(mesh);assert.ok(Math.abs(mesh.geometry.attributes.position.getZ(2)-.3)<1e-6);
 const before=Array.from(mesh.geometry.attributes.position.array);
 view.follow(new THREE.PerspectiveCamera());assert.deepEqual(Array.from(mesh.geometry.attributes.position.array),before);
 assert.equal(mesh.material.userData.environmentFamily,'fabric');assert.ok(mesh.material.map);assert.ok(mesh.material.bumpMap);
 view.dispose();assert.equal(scene.children.length,0);
});
test('duplicate rigid instances preserve server identities and rotation about their true centre',async()=>{
 const view=mountSurround(new THREE.Scene(),{api:async()=>geometry});
 await view.apply({world:{},placements:[]},{id:'studio'},['ball-small','ball-small'],[{id:'ball-small',geometry_url:'/ball'}]);
 view.update({objects:[{id:'ball-small-2',kind:'rigid',position_m:[2,3,4],origin_m:[1,0,0],rotation_matrix:[[0,-1,0],[1,0,0],[0,0,1]]}]});
 const second=view.group.children.find(o=>o.userData.environmentInstance==='ball-small-2');
 assert.ok(second);assert.ok(second.position.distanceTo(new THREE.Vector3(2,2,4))<1e-10);
 view.dispose();
});
test('material families provide fabric relief, painted surfaces and metallic rails',()=>{
 assert.equal(materialFamily('bed-rails'),'metal');assert.equal(materialFamily('wall-clinical'),'paint');
 for(const name of ['blanket','floor-wood','wall-clinical','bed-rails']){
  const material=environmentMaterial(name,0xffffff);assert.ok(material.map);assert.ok(material.bumpMap);material.dispose();
 }
});
test('unequal ground parts keep paving material confined to the patio and batch draw calls',async()=>{
 const positions=[],indices=[];
 for(let f=0;f<36;f++){const x=f*.1;positions.push(x,0,0,x+.05,0,0,x,0,.05);indices.push(f*3,f*3+1,f*3+2);}
 const data={id:'ground',positions,indices,parts:[{name:'sheet',primitive:'tiled_plane',colour:'grass'},{name:'patio',primitive:'box',colour:'patio'}]};
 const view=mountSurround(new THREE.Scene(),{api:async()=>data});
 await view.apply({world:{surround_url:'/ground'}},{id:'floor',gravity:[0,-9.81,0]});
 const mesh=view.group.children[0];
 assert.deepEqual(mesh.geometry.groups,[{start:0,count:72,materialIndex:0},{start:72,count:36,materialIndex:1}]);
 assert.equal(mesh.material[0].userData.environmentFamily,'grass');assert.equal(mesh.material[1].userData.environmentFamily,'stone');
 view.dispose();
});
