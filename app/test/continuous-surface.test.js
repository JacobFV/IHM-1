import test from 'node:test';
import assert from 'node:assert/strict';
import {prepareContinuousAsset,bindContinuousVertices,poseContinuousVertices,continuousMaterialAnchor,continuousMaterialPorts} from '../src/continuous-surface.js';
const I=[[1,0,0],[0,1,0],[0,0,1]],R=[[0,-1,0],[1,0,0],[0,0,1]];
function payload(){return {schema:'ihm.continuous-surface-binding.v1',binding_identity:'fixture',segments:[{id:'a',bone_id:'bone-a',reference_centroid_m:[0,0,0]},{id:'b',bone_id:'bone-b',reference_centroid_m:[1,0,0]}],reference_positions_m:[[1,2,3],[-1,0,0],[1,0,0]],weights:[[.25,.75],[1,0],[0,1]]};}
const transforms={a:{centroid_m:[2,0,0],rotation_matrix:R},b:{centroid_m:[1,3,0],rotation_matrix:I}};
const near=(a,b)=>{assert.equal(a.length,b.length);a.forEach((v,i)=>assert.ok(Math.abs(v-b[i])<1e-12,`${i}: ${v} != ${b[i]}`));};
const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
test('continuous attachment blends native stations with unrounded weights',()=>{
 const asset=prepareContinuousAsset(payload()),reference=asset.reference.slice(),bound=bindContinuousVertices(reference,asset,{skin:true});
 const output=poseContinuousVertices(reference,bound,transforms,new Float64Array(reference.length));
 near([...output.slice(0,3)],[.75,4,3]);
 const common={a:{centroid_m:[5,6,7],rotation_matrix:R},b:{centroid_m:[5,7,7],rotation_matrix:R}};
 near([...poseContinuousVertices(reference,bound,common,new Float64Array(reference.length)).slice(0,3)],[3,7,10]);
 const p=payload();p.weights[0]=[.123456789123456,.876543210876544];
 assert.equal(prepareContinuousAsset(p).weights[0],p.weights[0][0]);
});
test('continuous material ports conserve force, moment and virtual work',()=>{
 const asset=prepareContinuousAsset(payload()),reference=asset.reference.slice(),bound=bindContinuousVertices(reference,asset,{skin:true}),anchor=continuousMaterialAnchor(bound,reference,0),force=[2,-3,4];
 const {point_m,ports}=continuousMaterialPorts(anchor,transforms,force);near(point_m,[.75,4,3]);assert.equal(ports.length,2);
 near([0,1,2].map(k=>ports.reduce((s,p)=>s+p.force_n[k],0)),force);
 near([0,1,2].map(k=>ports.reduce((s,p)=>s+cross(p.point_m,p.force_n)[k],0)),cross(point_m,force));
 const velocities=[[1,2,3],[-2,3,1]],omegas=[[.1,.2,.3],[.3,-.2,.1]];
 const stationVelocity=ports.map((p,j)=>{const v=cross(omegas[j],p.point_m.map((x,k)=>x-transforms[j?'b':'a'].centroid_m[k]));return v.map((x,k)=>x+velocities[j][k]);});
 const blended=[0,1,2].map(k=>stationVelocity.reduce((s,v,j)=>s+anchor.weights[j]*v[k],0));
 assert.ok(Math.abs(dot(force,blended)-ports.reduce((s,p,j)=>s+dot(p.force_n,stationVelocity[j]),0))<1e-12);
 assert.equal(continuousMaterialPorts(continuousMaterialAnchor(bound,reference,1),transforms,force).ports.length,1);
 assert.equal(continuousMaterialPorts(anchor,transforms,[0,0,0]).ports.length,0);
});
test('garments transfer nearest rest vertex and exact ties use lowest source index',()=>{
 const asset=prepareContinuousAsset(payload()),reference=new Float64Array([0,0,0,1,0,0]),bound=bindContinuousVertices(reference,asset,{skin:false});
 assert.deepEqual([...bound.sourceIndices],[1,2]);
 poseContinuousVertices(reference,bound,transforms,new Float64Array(reference.length));
 assert.deepEqual([...bound.sourceIndices],[1,2]);
 const moved=reference.slice();moved[0]=9;
 assert.throws(()=>poseContinuousVertices(moved,bound,transforms,new Float64Array(reference.length)),/reference changed/);
});
test('full skin uses exact source index even at duplicated coordinates; rejects different topology',()=>{
 const p=payload();p.reference_positions_m[2]=[...p.reference_positions_m[1]];
 const asset=prepareContinuousAsset(p),reference=Float32Array.from(asset.reference),bound=bindContinuousVertices(reference,asset,{skin:true});
 assert.deepEqual([...bound.sourceIndices],[0,1,2]);assert.equal(bound.weights[5],1);
 const changed=reference.slice();changed[0]+=.001;assert.throws(()=>bindContinuousVertices(changed,asset,{skin:true}),/do not match/);
 assert.throws(()=>bindContinuousVertices(reference.slice(0,3),asset,{skin:true}),/do not match/);
 assert.throws(()=>poseContinuousVertices(reference,bound,transforms,reference),/output buffer/);
});
test('rejects invalid sidecars and missing transforms',()=>{
 for(const mutate of [p=>p.schema='other',p=>p.segments[1].id='a',p=>p.segments[0].reference_centroid_m[0]=NaN,p=>p.reference_positions_m[0][0]=Infinity,p=>p.weights[0][0]=-.1,p=>p.weights[0][0]=.5,p=>p.weights[0]=[1]]){const p=payload();mutate(p);assert.throws(()=>prepareContinuousAsset(p));}
 const asset=prepareContinuousAsset(payload()),bound=bindContinuousVertices(asset.reference,asset,{skin:true});
 assert.throws(()=>poseContinuousVertices(asset.reference,bound,{},new Float64Array(asset.reference.length)),/transform/);
});
test('KD transfer agrees with exhaustive search including duplicate ties',()=>{
 const p=payload();p.reference_positions_m=Array.from({length:500},(_,i)=>[Math.sin(i*3),Math.cos(i*7),Math.sin(i*11)]);p.reference_positions_m[123]=[...p.reference_positions_m[12]];p.weights=p.reference_positions_m.map(()=>[.2,.8]);
 const asset=prepareContinuousAsset(p),reference=Float64Array.from(Array.from({length:200},(_,i)=>[Math.cos(i),Math.sin(i*2),Math.cos(i*3)]).flat());
 const bound=bindContinuousVertices(reference,asset);
 for(let i=0;i<200;i++){let best=Infinity,index=-1;for(let j=0;j<500;j++){const d=Math.sqrt([0,1,2].reduce((s,k)=>s+(reference[3*i+k]-asset.reference[3*j+k])**2,0));if(d<best){best=d;index=j;}}assert.equal(bound.sourceIndices[i],index);}
 const duplicate=bindContinuousVertices(Float64Array.from(p.reference_positions_m[123]),asset);assert.equal(duplicate.sourceIndices[0],12);
});
test('Float32 rendering input retains authoritative Float64 skin material coordinates',()=>{
 const p=payload();p.reference_positions_m[0]=[.1,.2,.3];
 const asset=prepareContinuousAsset(p),reference=Float32Array.from(asset.reference),bound=bindContinuousVertices(reference,asset,{skin:true});
 const rest={a:{centroid_m:[0,0,0],rotation_matrix:I},b:{centroid_m:[1,0,0],rotation_matrix:I}};
 const output=poseContinuousVertices(reference,bound,rest,new Float64Array(reference.length));near([...output.slice(0,3)],[.1,.2,.3]);
 assert.notEqual(reference[0],bound.reference[0]);
 near(continuousMaterialPorts(continuousMaterialAnchor(bound,reference,0),rest,[1,0,0]).point_m,[.1,.2,.3]);
});
