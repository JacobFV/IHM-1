import test from 'node:test';
import assert from 'node:assert/strict';
import {LiveBodyHistory,bodyEndpoint,bodyEnvironment,bodyCommand,signalInfo,frameScope,materialOffset,materialPoint} from '../src/embodied-live.js';
const frame=(time,sequence=0)=>({schema:'ihm.embodied-frame.v1',time_s:time,sequence,entities:{},mechanics:{muscles:{bra_r:{}}},physiology:{values:{heart_rate_per_min:72,unknown:null},units:{heart_rate_per_min:'1/min'}}});
test('default live owner never silently selects reduced scene',()=>{
 assert.equal(bodyEndpoint('embodied'),'/api/embodied/sessions');
 assert.equal(bodyEnvironment('bed','embodied'),'supine');
 assert.equal(bodyEnvironment('floor','embodied'),'upright');
 assert.throws(()=>bodyEndpoint('fake'));assert.throws(()=>bodyEnvironment('fake','embodied'));
 assert.deepEqual(bodyCommand('embodied',3,[],{descending:{bra_r:.3},sensory_blocks:['bra_r'],motor_blocks:[],skin_compression_pa:100}),{seconds:.02,sequence:3,forces:[],descending:{bra_r:.3},sensory_blocks:['bra_r'],motor_blocks:[],skin_compression_pa:100});
 assert.throws(()=>bodyCommand('embodied',1,[],{skin_compression_pa:NaN}));
});
test('live histories are bounded, preserve real time and never blend session owners',()=>{
 const h=new LiveBodyHistory(3);h.push('a',frame(0));h.push('a',frame(.02,1));h.push('a',frame(.04,2));h.push('a',frame(.06,3));
 assert.deepEqual(h.series('heart_rate_per_min').time_s,[.02,.04,.06]);
 assert.equal(h.push('a',frame(.06,3)),false);
 h.push('b',frame(0));assert.deepEqual(h.series('heart_rate_per_min').time_s,[0]);
 assert.equal(h.series('unknown').values[0],null);
 assert.throws(()=>h.push('b',frame(-1,1)));
 assert.equal(signalInfo(frame(0),'heart_rate_per_min').unit,'1/min');
 assert.match(frameScope(frame(0)),/calibration|uncalibrated/i);
});
test('grab point follows current affine and rigid material transform',()=>{
 const entity={centroid_m:[1,2,3],rotation_matrix:[[0,-1,0],[1,0,0],[0,0,1]],deformation_gradient:[[2,0,0],[0,1,0],[0,0,1]]};
 const offset=materialOffset([1,2.2,3],entity);assert.ok(Math.abs(offset[0]-.1)<1e-12);
 const next={...entity,centroid_m:[2,2,3]};assert.ok(Math.abs(materialPoint(offset,next)[1]-2.2)<1e-12);
});
import {MotorInputs} from '../src/embodied-panels.js';
test('muscle inputs retain selective targets and never create unknown effectors',()=>{
 const m=new MotorInputs();m.bind(['bra_r','bra_l']);m.set('bra_r',.4,true,false);
 assert.deepEqual(m.snapshot(),{descending:{bra_r:.4},sensory_blocks:['bra_r'],motor_blocks:[],skin_compression_pa:0});
 assert.throws(()=>m.set('fake',1,false,false));m.set('bra_l',0,false,true);
 assert.equal(m.snapshot().motor_blocks[0],'bra_l');m.bind(['bra_l']);assert.deepEqual(m.snapshot().descending,{});assert.deepEqual(m.snapshot().sensory_blocks,[]);
 m.clear();assert.deepEqual(m.snapshot().motor_blocks,[]);
});
import {createBodyOwner,closeBodyOwner} from '../src/embodied-live.js';
test('lost creation recovers only a new actor without a repeated create',async()=>{
 const calls=[];let lists=0;
 const result=await createBodyOwner(async(path,data)=>{calls.push([path,data]);if(data)throw Error('Lost response');if(path==='/body')return {sessions:lists++?[{id:'new',closed:false}]:[]};return {id:'new',status:'initializing',closed:false};},'/body',{environment:'supine'});
 assert.equal(result.id,'new');assert.equal(calls.filter(c=>c[1]).length,1);
 await assert.rejects(createBodyOwner(async(path,data)=>{if(data)throw Error('Resource occupied');return {sessions:[{id:'old',closed:false}]};},'/body',{}),/Resource occupied/);
});
test('source handoff waits for confirmed pending body cleanup',async()=>{
 let polls=0,waits=0;
 const result=await closeBodyOwner(async(path,data)=>data?{status:'closing',closed:false}:++polls<2?{status:'closing',closed:false}:{status:'closed',closed:true},'/body/id',async()=>waits++);
 assert.equal(result.closed,true);assert.equal(waits,2);
 await assert.rejects(closeBodyOwner(async()=>{throw Error('Offline');},'/body/id',async()=>{}),/Offline/);
});
