import test from 'node:test';
import assert from 'node:assert/strict';
import {LiveBodyHistory,bodyEndpoint,bodyEnvironment,bodyCommand,signalInfo,frameScope,materialOffset,materialPoint} from '../src/embodied-live.js';
const frame=(time,sequence=0)=>({schema:'ihm.embodied-frame.v1',time_s:time,sequence,entities:{},mechanics:{muscles:{bra_r:{}}},physiology:{values:{heart_rate_per_min:72,unknown:null},units:{heart_rate_per_min:'1/min'}}});
test('default live owner never silently selects reduced scene',()=>{
 assert.equal(bodyEndpoint('embodied'),'/api/embodied/sessions');
 assert.equal(bodyEnvironment('bed','embodied'),'supine');
 assert.equal(bodyEnvironment('floor','embodied'),'upright');
 assert.throws(()=>bodyEndpoint('fake'));assert.throws(()=>bodyEnvironment('fake','embodied'));
 assert.deepEqual(bodyCommand('embodied',3,[],{descending:{bra_r:.3},sensory_blocks:['bra_r'],motor_blocks:[],skin_compression_pa:100},{input_capabilities:{skin_pressure:{whole_skin:true,regional_ids:[]}}}),{seconds:.02,sequence:3,forces:[],descending:{bra_r:.3},sensory_blocks:['bra_r'],motor_blocks:[],skin_compression_pa:100});
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
test('schedule sequence updates do not fabricate physiology time samples',()=>{
 const h=new LiveBodyHistory(3);h.push('body',frame(0));h.push('body',frame(.02,1));
 const scheduled={...frame(.02,2),intake_schedule:{events:[{event_id:'drink'}]}};
 assert.equal(h.push('body',scheduled),false);
 assert.equal(h.frame,scheduled);assert.equal(h.sequence,2);
 assert.deepEqual(h.series('heart_rate_per_min').time_s,[0,.02]);
 h.push('body',frame(.04,3));assert.deepEqual(h.series('heart_rate_per_min').time_s,[0,.02,.04]);
 assert.throws(()=>h.push('body',frame(.02,4)),/clock reversed/);
});
import {scheduleBodyIntakes} from '../src/embodied-live.js';
test('intake scheduling uses the current body sequence and recovers without replay',async()=>{
 const events=[{event_id:'water',time_s:.02,meal:{name:'water',water_ml:250}}];
 const accepted={...frame(.02,5),intake_schedule:{events}};const calls=[];
 const result=await scheduleBodyIntakes(async(path,data)=>{
  calls.push({path,data});if(data)throw Error('Lost response');return accepted;
 },'/body/id',4,events);
 assert.deepEqual(calls,[{path:'/body/id/intakes',data:{sequence:4,events}},{path:'/body/id',data:undefined}]);
 assert.equal(result.frame,accepted);assert.equal(result.recovered,true);
 assert.equal(result.frame.time_s,.02);
});
test('intake double network failure never retries POST and invalid sequences never send',async()=>{
 let posts=0;const request=async(path,data)=>{if(data)posts++;throw Error('Offline');};
 await assert.rejects(scheduleBodyIntakes(request,'/body/id',0,[{event_id:'water'}]),/Offline/);assert.equal(posts,1);
 await assert.rejects(scheduleBodyIntakes(request,'/body/id',-1,[{}]),/Missing current/);assert.equal(posts,1);
});
test('only pre-mutation 400 intake rejection is marked correctable while refreshing the owner',async()=>{
 for(const httpStatus of [400,409,422,503,undefined]){
  const error=Object.assign(Error('Request failed'),{httpStatus});const calls=[];
  const result=await scheduleBodyIntakes(async(path,data)=>{calls.push(path);if(data)throw error;return frame(.02,2);},'/body/id',1,[{event_id:'water'}]);
  assert.equal(result.error.definitelyRejected===true,httpStatus===400);
  assert.deepEqual(calls,['/body/id/intakes','/body/id']);
 }
});
test('known rejection remains correctable if read recovery is offline',async()=>{
 const error=Object.assign(Error('Past intake time'),{httpStatus:400});let posts=0;
 await assert.rejects(scheduleBodyIntakes(async(path,data)=>{if(data){posts++;throw error;}throw Error('Offline');},'/body/id',2,[{event_id:'water'}]),e=>e===error&&e.definitelyRejected===true);
 assert.equal(posts,1);
});

test('controller telemetry distinguishes ablation from evidence of cortical control',async()=>{
 const {controllerReadouts}=await import('../src/embodied-live.js');
 const rows=Object.fromEntries(controllerReadouts({controller:{kind:'implicit',sever:true,checkpoint_sha256:'abc'},neural:{motor_excitations:{a:.1,b:.4},arc_max:{stretch:0,renshaw:.2}}}));
 assert.equal(rows.Controller,'IBM fused implicit kernel');
 assert.equal(rows.Ablation,'Cortical kernel severed');
 assert.equal(rows['Motor excitation peak'],.4);
 assert.equal(rows['Stretch arc peak'],0);
 assert.equal(rows['Renshaw arc peak'],.2);
 assert.equal(rows['Autogenic arc peak'],null);
 assert.match(rows['Cortical contribution'],/matched full \/ severed/);
 assert.equal(Object.fromEntries(controllerReadouts({}))['Motor excitation peak'],null);
 assert.equal(Object.fromEntries(controllerReadouts({})).Controller,'Unavailable');
});

test('deformable object grab follows one free material node through motion',async()=>{
 const {softObjectAnchor}=await import('../src/embodied-live.js');
 const object={positions:[0,0,0,.1,0,0,.2,0,0],fixed_nodes:[0]};
 const first=softObjectAnchor(object,[0,0,0]);
 assert.equal(first.nodeIndex,1);
 assert.deepEqual(first.point,[.1,0,0]);
 const moved={...object,positions:[0,0,0,.3,.4,0,.1,0,0]};
 assert.deepEqual(softObjectAnchor(moved,[.1,0,0],first.nodeIndex).point,[.3,.4,0]);
 assert.throws(()=>softObjectAnchor({positions:[0,0,0],fixed_nodes:[0]},[0,0,0]),/No movable/);
});

test('learned primitive configuration restricts targets and incompatible ablations',async()=>{
 const {controllerConfiguration,controllerReadouts}=await import('../src/embodied-live.js');
 assert.deepEqual(controllerConfiguration('implicit_cortical_ankle','no-cord',.1),{kind:'implicit_cortical_ankle',sever:false,no_cord:true,target_rad:.1});
 assert.deepEqual(controllerConfiguration('implicit_ankle_primitive','sever',-.1),{kind:'implicit_ankle_primitive',sever:true,target_rad:-.1});
 assert.throws(()=>controllerConfiguration('implicit_ankle_primitive','no-cord'),/Unsupported/);
 assert.throws(()=>controllerConfiguration('implicit_ankle_primitive','full',NaN),/target/);
 assert.throws(()=>controllerConfiguration('implicit_ankle_primitive','full',.3),/target/);
 const rows=Object.fromEntries(controllerReadouts({controller:{kind:'implicit_ankle_primitive',target_rad:.12},neural:{motor_primitive:{angle_rad:.01}}}));
 assert.equal(rows['Ankle angle · rad'],.01);assert.match(rows['Ankle feedback'],/privileged/);
});

test('engineered balance has no ablations and never claims cortical control',async()=>{
 const {controllerConfiguration,controllerReadouts}=await import('../src/embodied-live.js');
 assert.deepEqual(controllerConfiguration('engineering_stance'),{kind:'engineering_stance'});
 assert.throws(()=>controllerConfiguration('engineering_stance','sever'),/Unsupported/);
 const rows=Object.fromEntries(controllerReadouts({controller:{kind:'engineering_stance',motor_owner:'engineered-postural-feedback'},neural:{motor_excitations:{soleus_r:.3}}}));
 assert.equal(rows.Controller,'Engineered balance');
 assert.equal(rows['Motor owner'],'engineered-postural-feedback');
 assert.equal(rows.Ablation,'Not applicable');
 assert.equal(rows['Cortical contribution'],'Not used by this controller');
 assert.equal(rows['Stretch arc peak'],'Not provided by this controller');
 assert.match(rows.Task,/walking not established/);
});

test('learned cortical stance preserves explicit severing and fixed no-cord topology',async()=>{
 const {controllerConfiguration,controllerReadouts}=await import('../src/embodied-live.js');
 assert.deepEqual(controllerConfiguration('implicit_cortical_stance','full'),{kind:'implicit_cortical_stance',sever:false,no_cord:false});
 assert.deepEqual(controllerConfiguration('implicit_cortical_stance','sever'),{kind:'implicit_cortical_stance',sever:true,no_cord:false});
 assert.throws(()=>controllerConfiguration('implicit_cortical_stance','no-cord'),/Unsupported/);
 const rows=Object.fromEntries(controllerReadouts({neural:{controller:{kind:'implicit_cortical_stance',sites:1024,muscle_count:98,motor_owner:'trained1024site-persistent-IBM-EI-cortex',sampling_interval_s:.01},cortical_stance:{max_cortical_correction:.02}}}));
 assert.match(rows.Controller,/Experimental.*1024/);assert.equal(rows['Controlled muscles'],98);assert.equal(rows['Cortical correction peak'],.02);
 assert.equal(rows['Stretch arc peak'],'Not provided by this controller');assert.match(rows['Balance feedback'],/Privileged/);assert.match(rows['Cortical contribution'],/matched/);assert.match(rows.Task,/walking not established/);
});

test('curriculum16 kernel selections keep their own ablations and disclose the kernel lineage',async()=>{
 const {controllerConfiguration,controllerReadouts}=await import('../src/embodied-live.js');
 // The retained-kernel raw path keeps the cord ablation the trained stance path cannot have.
 assert.deepEqual(controllerConfiguration('implicit_curriculum16','no-cord'),{kind:'implicit_curriculum16',sever:false,no_cord:true});
 assert.deepEqual(controllerConfiguration('implicit_curriculum16_stance','sever'),{kind:'implicit_curriculum16_stance',sever:true,no_cord:false});
 assert.throws(()=>controllerConfiguration('implicit_curriculum16_stance','no-cord'),/Unsupported/);
 const lineage='IBM-1 ckpt/ibm1_curriculum16.pt, the 16-objective consolidated kernel, retained locally';
 const rows=Object.fromEntries(controllerReadouts({neural:{controller:{kind:'implicit_curriculum16_stance',sites:1024,muscle_count:98,sampling_interval_s:.01,kernel_identity:lineage,checkpoint_sha256:'a'.repeat(64)},cortical_stance:{max_cortical_correction:.01}}}));
 assert.match(rows.Controller,/16-objective/);assert.equal(rows['Association kernel'],lineage);
 assert.equal(rows['Base kernel SHA-256'],'a'.repeat(64));assert.equal(rows['Controlled muscles'],98);
 assert.match(rows['Cortical contribution'],/matched/);assert.match(rows.Walking,/Not established/);
 const raw=Object.fromEntries(controllerReadouts({neural:{controller:{kind:'implicit_curriculum16',sites:1024,no_cord:true,kernel_name:'kernel.pt'},arc_max:{stretch:0}}}));
 assert.match(raw.Controller,/16-objective/);assert.equal(raw['Stretch arc peak'],'Inactive for this controller');
});
