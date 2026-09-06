// Bounded synthetic guide benchmark; no assets, UI, native code or workers.
// node --max-old-space-size=256 --expose-gc scripts/benchmark_hair_dynamics.mjs [source-file]
// A saved pre-change hair_dynamics.js can be supplied for before/after runs.
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {resolve} from 'node:path';
import {performance} from 'node:perf_hooks';

const current=new URL('../app/src/hair_dynamics.js',import.meta.url);
const sourcePath=process.argv[2]?pathToFileURL(resolve(process.argv[2])):current;
const source=readFileSync(sourcePath,'utf8').replace("'./state.js'",JSON.stringify(new URL('../app/src/state.js',import.meta.url).href));
const {ElasticHairSystem}=await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
const guides=128,nodes=7,updates=60,dt=.1;
const data={centerlines_m:[],strand_offsets:[0],radius_m:[],tensile_modulus_pa:7.11e9,bending_modulus_pa:5.7e9,density_kg_m3:1312};
for(let s=0;s<guides;s++){
 for(let j=0;j<nodes;j++)data.centerlines_m.push(j*.005,0,s*.001);
 data.strand_offsets.push((s+1)*nodes);data.radius_m.push(40e-6);
}
const system=new ElasticHairSystem(data);
let coordinateCopyAllocations=0;
system.positions.slice=function(...args){coordinateCopyAllocations++;return Float64Array.prototype.slice.apply(this,args);};
const roots=system.root.slice(),directions=system.tangents.slice(),forces=new Float64Array(system.positions.length);
for(let s=0;s<guides;s++)forces[3*((s+1)*nodes-1)+2]=1e-8;
const advance=(i)=>{
 for(let s=0;s<guides;s++){roots[3*s]=1e-4*Math.sin(i*.03);directions[3*s+1]=.001*Math.sin(i*.02);}
 system.step(dt,{roots_m:roots,root_tangents:directions,external_forces_n:forces});
};
for(let i=0;i<8;i++)advance(i);
coordinateCopyAllocations=0;
global.gc?.();
const before=process.memoryUsage(),start=performance.now();
for(let i=0;i<updates;i++)advance(i+8);
const elapsed=performance.now()-start,after=process.memoryUsage();
const hash=createHash('sha256');
for(const array of [system.positions,system.velocities,system.root,system.tangents])hash.update(Buffer.from(array.buffer,array.byteOffset,array.byteLength));
console.log(JSON.stringify({source:sourcePath.pathname,guides,nodes,updates,substeps_per_update:24,
 elapsed_ms:elapsed,ms_per_update:elapsed/updates,
 peak_rss_kib:process.resourceUsage().maxRSS,rss_before_bytes:before.rss,rss_after_bytes:after.rss,
 heap_used_delta_bytes:after.heapUsed-before.heapUsed,
 full_coordinate_copy_allocations:coordinateCopyAllocations,
 full_coordinate_copy_allocated_bytes:coordinateCopyAllocations*system.positions.byteLength,
 retained_substep_scratch_bytes:system.beforeSubstep?.byteLength||0,
 final_state_sha256:hash.digest('hex'),
 limitations:'Synthetic guide solver only; wall time and RSS are noisy. Excludes root attachment, interpolation, tube uploads, UI and collision (none). Dynamics remain opt-in.',
},null,2));
