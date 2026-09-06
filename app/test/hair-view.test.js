import test from 'node:test';
import assert from 'node:assert/strict';
import * as THREE from 'three';
import {attachElasticHair,updateElasticHair as updateView} from '../src/hair-view.js';
import {createHairController} from '../src/hair_dynamics.js';
const updateElasticHair=(object,options)=>updateView(object,{...options,controllerFactory:createHairController});
function fixture(){
 const g={strands:{centerlines_m:[0,0,0,.005,0,0,.01,0,0],strand_offsets:[0,3],radius_m:[40e-6],tensile_modulus_pa:7.11e9,bending_modulus_pa:5.7e9,density_kg_m3:1312},
 attachment:{kind:'ElasticHairMaterialPoint',skin_entity_id:'skin',barycentric:[[1,0,0]],reference_triangles_m:[0,0,0,0,1,0,0,0,1]}};
 const geometry=new THREE.BufferGeometry();geometry.setAttribute('position',new THREE.BufferAttribute(new Float32Array(36),3));
 const object=new THREE.Mesh(geometry);return {g,object};
}
test('elastic hair receives skin motion once and retains an identity object matrix',()=>{
 const {g,object}=fixture();assert.equal(attachElasticHair(object,g),true);
 const frame={time_s:0,entities:{skin:{translation_m:[.1,.2,.3]}}};
 assert.equal(updateElasticHair(object,{frame,referenceCentroids:{skin:[0,0,0]},recordKey:'fixture',visible:true,enabled:true}),true);
 const positions=object.geometry.attributes.position.array;
 const root=[0,1,2].map(k=>(positions[k]+positions[3+k]+positions[6+k]+positions[9+k])/4);
 assert.ok(root.every((v,i)=>Math.abs(v-[.1,.2,.3][i])<1e-7));
 assert.deepEqual(object.matrix.elements,new THREE.Matrix4().elements);
 assert.equal(object.matrixAutoUpdate,false);
 assert.equal(object.userData.hairDiagnostics.physical_radius_multiplier,1);
});
test('hidden hair does no GPU update or time advance; rewind remains explicit',()=>{
 const {g,object}=fixture();attachElasticHair(object,g);
 updateElasticHair(object,{frame:{time_s:0,entities:{}},recordKey:'a',visible:true,enabled:true});
 const version=object.geometry.attributes.position.version;
 updateElasticHair(object,{frame:{time_s:.05,entities:{}},recordKey:'a',visible:false,enabled:true});
 assert.equal(object.geometry.attributes.position.version,version);
 assert.equal(object.userData.elasticHair.system.time,0);
 updateElasticHair(object,{frame:{time_s:.1,entities:{}},recordKey:'a',visible:true,enabled:true});
 assert.equal(object.userData.hairDiagnostics.reset_reason,'resume_without_catchup');
 updateElasticHair(object,{frame:{time_s:0,entities:{}},recordKey:'a',visible:true,enabled:true});
 assert.equal(object.userData.hairDiagnostics.reset_reason,'rewind');
});
test('default static strands never construct or advance an elastic controller',()=>{
 const {g,object}=fixture();const reference=object.geometry.attributes.position.array.slice();
 attachElasticHair(object,g);
 for(const time_s of [0,.1,1,1000])updateElasticHair(object,{frame:{time_s,entities:{}},visible:true});
 assert.equal(object.userData.elasticHair,undefined);
 assert.deepEqual(object.geometry.attributes.position.array,reference);
 assert.equal(object.geometry.attributes.position.version,0);
 assert.equal(object.userData.hairDiagnostics.mode,'static');
 updateElasticHair(object,{frame:{time_s:1001,entities:{skin:{translation_m:[.1,.2,.3]}}},referenceCentroids:{skin:[0,0,0]}});
 assert.deepEqual(object.matrix.elements.slice(12,15),[.1,.2,.3]);
 assert.equal(object.visible,true);
 assert.deepEqual(object.geometry.attributes.position.array,reference);
 updateElasticHair(object,{frame:{time_s:0,entities:{}},visible:true,enabled:true});
 assert.ok(object.userData.elasticHair);
 updateElasticHair(object,{frame:{time_s:.1,entities:{}},visible:true,enabled:false});
 assert.equal(object.userData.elasticHair,undefined);
 assert.deepEqual(object.geometry.attributes.position.array,reference);
});

test('worker results update the view asynchronously and cannot revive disabled dynamics',()=>{
 const {g,object}=fixture();attachElasticHair(object,g);
 let deliver,disposed=0,notified=0,submitted;
 const controllerFactory=(_,callback)=>{deliver=callback;return {update:request=>{submitted=request;},dispose:()=>disposed++};};
 const frame={time_s:0,entities:{skin:{translation_m:[.1,.2,.3]},unrelated:{large:'omitted'}}};
 const referenceCentroids={skin:[0,0,0]};
 updateView(object,{frame,referenceCentroids,enabled:true,controllerFactory,onDiagnostics:()=>notified++});
 assert.equal(object.geometry.attributes.position.version,0);
 assert.equal(submitted.frame.entities.unrelated,undefined);
 assert.equal(object.userData.hairDiagnostics.worker_pending,true);
 const result=createHairController(g).update({frame,referenceCentroids,simulationTime:0});
 deliver(result);assert.equal(notified,1);assert.deepEqual(object.geometry.attributes.position.array,result.positions);
 updateView(object,{frame,referenceCentroids,enabled:false});const version=object.geometry.attributes.position.version;
 deliver(result);assert.equal(disposed,1);assert.equal(object.geometry.attributes.position.version,version);
 assert.equal(object.userData.hairDiagnostics.mode,'static');
});
