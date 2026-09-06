// Tiny fixture only: node --max-old-space-size=256 scripts/benchmark_hair_worker.mjs
import {Worker} from 'node:worker_threads';
import {performance} from 'node:perf_hooks';
import assert from 'node:assert/strict';
import {createHairController} from '../app/src/hair_dynamics.js';
import {createHairWorkerPool} from '../app/src/hair-worker-client.js';

const strands={centerlines_m:[],strand_offsets:[0],radius_m:[],tensile_modulus_pa:7.11e9,bending_modulus_pa:5.7e9,density_kg_m3:1312};
const attachment={skin_entity_id:'skin',barycentric:[],reference_triangles_m:[]};
const guides=128,nodes=7,updates=40;
for(let s=0;s<guides;s++){
 for(let j=0;j<nodes;j++)strands.centerlines_m.push(j*.005,0,s*.001);
 strands.strand_offsets.push((s+1)*nodes);strands.radius_m.push(40e-6);
 attachment.barycentric.push([1,0,0]);attachment.reference_triangles_m.push(0,0,s*.001,0,1,s*.001,0,0,s*.001+1);
}
const geometry={strands,attachment},request=i=>({simulationTime:i*.1,frame:{time_s:i*.1,entities:{}},recordKey:'fixture',visible:true,active:true});
const direct=createHairController(geometry);direct.update(request(0));
let expected;
const syncStart=performance.now();
for(let i=1;i<=updates;i++)expected=direct.update(request(i));
const syncMs=performance.now()-syncStart;
const moduleUrl=new URL('../app/src/hair-worker.js',import.meta.url).href;
let rawWorker,deliver,dispatchMs=0,heartbeats=0;
const pool=createHairWorkerPool(()=>{
 rawWorker=new Worker(`const {parentPort}=require('node:worker_threads');import(${JSON.stringify(moduleUrl)}).then(({createHairWorkerHandler})=>{const handler=createHairWorkerHandler((data,transfer)=>parentPort.postMessage(data,transfer));parentPort.on('message',data=>handler({data}));});`,{eval:true,resourceLimits:{maxOldGenerationSizeMb:64}});
 const adapter={postMessage:message=>rawWorker.postMessage(message),terminate:()=>rawWorker.terminate()};
 rawWorker.on('message',data=>adapter.onmessage({data}));rawWorker.on('error',error=>adapter.onerror(error));
 return adapter;
});
const controller=pool.createController(geometry,result=>deliver(result));
const advance=i=>new Promise((resolve,reject)=>{
 deliver=result=>result.diagnostics.worker_error?reject(Error(result.diagnostics.worker_error)):resolve(result);
 const start=performance.now();controller.update(request(i));dispatchMs+=performance.now()-start;
});
const deadline=setTimeout(()=>{controller.dispose();throw Error('Hair fixture exceeded 10 s');},10000);
try {
 await advance(0);dispatchMs=0;
 const heartbeat=setInterval(()=>heartbeats++,1),start=performance.now();let actual;
 try{for(let i=1;i<=updates;i++)actual=await advance(i);}finally{clearInterval(heartbeat);}
 const workerMs=performance.now()-start;
 assert.deepEqual(actual.positions,expected.positions);
 assert.equal(actual.diagnostics.unsimulated_total_s,0);
 console.log(JSON.stringify({guides,nodes,updates,synchronous_main_thread_ms:syncMs,
  worker_main_thread_dispatch_ms:dispatchMs,worker_roundtrip_ms:workerMs,main_thread_heartbeat_ticks:heartbeats,
  rss_bytes:process.memoryUsage().rss,exact_tube_coordinates_match:true,
  limitation:'Fixture only; excludes GPU normals/bounds/upload and real-viewer latency. One-way beam model, no contact or reaction impulses. Startup excluded from worker timing.'},null,2));
}finally{clearTimeout(deadline);controller.dispose();await rawWorker.terminate();}
