// Lightweight retained-asset verification; no browser or full population dynamics.
import fs from 'node:fs';
import zlib from 'node:zlib';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {createHairController,materialHairRoots} from '../app/src/hair_dynamics.js';
const directory=process.argv[2]||'data/derived/hair/elastic_v3';
const manifest=JSON.parse(fs.readFileSync(`${directory}/manifest_fragment.json`));
const hash=b=>crypto.createHash('sha256').update(b).digest('hex');
const results=[];
for(const structure of manifest.structures){
 const bytes=fs.readFileSync(structure.geometry_path);assert.equal(hash(bytes),structure.geometry_sha256);
 const g=JSON.parse(zlib.gunzipSync(bytes));const controller=createHairController(g);
 for(const [strands,attachment] of [[g.strands,g.attachment],[g.render_strands,g.render_attachment]]){
  const roots=materialHairRoots(attachment);for(let i=0;i<strands.radius_m.length;i++)for(let k=0;k<3;k++)assert.ok(Math.abs(roots.roots_m[3*i+k]-strands.centerlines_m[3*strands.strand_offsets[i]+k])<1e-9);
 }
 const start=controller.update({simulationTime:0});assert.equal(start.positions.length,g.positions.length);assert.ok(start.positions.every(Number.isFinite));
 const before=controller.system.positions.slice();const t=performance.now();const moved=controller.update({simulationTime:.1});const elapsed=performance.now()-t;
 assert.ok(moved.positions.every(Number.isFinite));assert.ok(moved.diagnostics.within_small_deflection);assert.ok(controller.system.maxStrain()<.01);
 let displacement=0;for(let i=0;i<before.length;i++)displacement=Math.max(displacement,Math.abs(before[i]-controller.system.positions[i]));assert.ok(displacement>0);
 const cached=moved.positions;assert.equal(controller.update({simulationTime:.11}).positions,cached);
 assert.equal(controller.update({simulationTime:1,visible:false}).diagnostics.paused,true);
 results.push({id:structure.id,maximum_guide_displacement_m:displacement,one_cold_update_ms:elapsed,...moved.diagnostics});
}
const report={schema:'hair_strand_verification_v1',passed:true,manifest_sha256:hash(fs.readFileSync(`${directory}/manifest_fragment.json`)),module_sha256:hash(fs.readFileSync('app/src/hair_dynamics.js')),results,limits:'One cold CPU update is not a frame-rate benchmark; no browser or collision acceptance.'};
fs.writeFileSync(`${directory}/verification.json`,JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report,null,2));
