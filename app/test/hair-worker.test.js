import test from 'node:test';
import assert from 'node:assert/strict';
import {createHairWorkerPool} from '../src/hair-worker-client.js';
import {createHairWorkerHandler} from '../src/hair-worker.js';
import {createHairController} from '../src/hair_dynamics.js';
import {Worker} from 'node:worker_threads';

const fixture=()=>({strands:{centerlines_m:[0,0,0,.005,0,0,.01,0,0],strand_offsets:[0,3],radius_m:[40e-6],tensile_modulus_pa:7.11e9,bending_modulus_pa:5.7e9,density_kg_m3:1312},attachment:{skin_entity_id:'skin',barycentric:[[1,0,0]],reference_triangles_m:[0,0,0,0,1,0,0,0,1]}});
const request=time=>({simulationTime:time,frame:{time_s:time,entities:{}},recordKey:'a',active:true,visible:true});
function harness(){
 const messages=[],results=[];let workers=0,terminated=0;
 const worker={postMessage:m=>messages.push(m),terminate:()=>terminated++};
 const pool=createHairWorkerPool(()=>{workers++;return worker;});
 const handler=createHairWorkerHandler(data=>worker.onmessage({data}));
 return {pool,results,messages,worker,flush(){const message=messages.shift();handler({data:message});},get workers(){return workers;},get terminated(){return terminated;}};
}
test('one worker bounds outstanding work and coalesces rendering requests without changing submitted time',()=>{
 const h=harness(),a=h.pool.createController(fixture(),r=>h.results.push(r)),b=h.pool.createController(fixture(),r=>h.results.push(r));
 a.update(request(0));a.update(request(.05));a.update(request(.1));b.update(request(.1));
 assert.equal(h.workers,1);assert.equal(h.messages.length,1);
 h.flush();assert.equal(h.messages.length,1);assert.equal(h.messages[0].request.simulationTime,.1);
 h.flush();h.flush();assert.equal(h.messages.length,0);
 assert.equal(h.results[1].diagnostics.coalesced_updates,1);
 assert.equal(h.results[1].diagnostics.integrated_interval_s,.1);
 a.dispose();b.dispose();assert.equal(h.terminated,1);
});
test('worker coordinates equal synchronous controller and disclose unintegrated gaps',()=>{
 const geometry=fixture(),controller=createHairController(geometry),outputs=[];
 const handler=createHairWorkerHandler(result=>outputs.push(result));
 for(const [sequence,time] of [0,.1,.2,1].entries()){
  const input={...request(time),frame:{time_s:time,entities:{skin:{translation_m:[time*.001,0,0]}}},referenceCentroids:{skin:[0,0,0]}};
  handler({data:{id:1,sequence,geometry:sequence===0?geometry:undefined,request:input}});
  const expected=controller.update(input),actual=outputs.at(-1);
  assert.deepEqual(actual.positions,expected.positions);
 }
 assert.equal(outputs.at(-1).diagnostics.reset_reason,'time_gap_not_simulated');
 assert.equal(outputs.at(-1).diagnostics.integrated_interval_s,0);
 assert.equal(outputs.at(-1).diagnostics.unsimulated_interval_s,.8);
 assert.equal(outputs.at(-1).diagnostics.physical_interval_status,'not_integrated');
});
test('hidden and disposed handles cannot apply stale worker results',()=>{
 const h=harness(),controller=h.pool.createController(fixture(),r=>h.results.push(r));
 controller.update(request(0));controller.update({...request(.1),visible:false});h.flush();
 assert.equal(h.results.length,0);
 controller.update(request(.2));assert.equal(h.messages[0].pausedSinceLast,true);h.flush();
 assert.equal(h.results[0].diagnostics.reset_reason,'resume_without_catchup');
 controller.update(request(.3));controller.dispose();h.flush();assert.equal(h.results.length,1);
});
test('worker failure is explicit and never executes a synchronous fallback',()=>{
 const h=harness(),controller=h.pool.createController(fixture(),r=>h.results.push(r));
 controller.update(request(0));h.worker.onerror({message:'fixture failure'});
 assert.equal(h.results[0].diagnostics.worker_error,'fixture failure');
 controller.update(request(.1));assert.equal(h.results.at(-1).diagnostics.worker_pending,false);
 assert.equal(h.terminated,1);
});

test('a solver error stops that region without dropping the other region queue',()=>{
 const h=harness(),bad=h.pool.createController(fixture(),r=>h.results.push(r)),good=h.pool.createController(fixture(),r=>h.results.push(r));
 bad.update(request(NaN));bad.update(request(.1));good.update(request(0));
 h.flush();assert.match(h.results[0].diagnostics.worker_error,/finite mechanical clock/);
 h.flush();assert.equal(h.results[1].diagnostics.updated,true);
 bad.update(request(.2));assert.equal(h.messages.length,0);assert.ok(h.results.at(-1).diagnostics.worker_error);
 bad.dispose();good.dispose();
});

test('actual isolated worker transfers the same fixture coordinates',{timeout:10000},async()=>{
 const moduleUrl=new URL('../src/hair-worker.js',import.meta.url).href;
 const worker=new Worker(`const {parentPort}=require('node:worker_threads');import(${JSON.stringify(moduleUrl)}).then(({createHairWorkerHandler})=>{const handler=createHairWorkerHandler((data,transfer)=>parentPort.postMessage(data,transfer));parentPort.on('message',data=>handler({data}));});`,{eval:true,resourceLimits:{maxOldGenerationSizeMb:64}});
 try {
  const geometry=fixture(),input=request(0);
  const result=await new Promise((resolve,reject)=>{worker.once('message',resolve);worker.once('error',reject);worker.postMessage({id:1,sequence:1,geometry,request:input});});
  assert.deepEqual(result.positions,createHairController(geometry).update(input).positions);
  assert.equal(result.diagnostics.physical_interval_status,'not_integrated');
 }finally{await worker.terminate();}
});
