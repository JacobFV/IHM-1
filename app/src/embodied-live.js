import {rotateOffset} from './scene-forces.js';
import {systemicLabel} from './systemic.js';
export function bodyEndpoint(kind='embodied') {
 if(kind==='embodied')return '/api/embodied/sessions';
 if(kind==='reduced')return '/api/scene/sessions';
 throw Error('Unknown body execution owner');
}
export function bodyEnvironment(view,kind='embodied') {
 const mapped={studio:'free',floor:'upright',bed:'supine'};
 if(!mapped[view])throw Error('Unknown body environment');
 return kind==='embodied'?mapped[view]:kind==='reduced'?view:bodyEndpoint(kind);
}
export function skinPressureCapability(frame) {
 const value=frame?.input_capabilities?.skin_pressure;
 const unavailable={topology:'unavailable',regional_ids:[],message:'Skin pressure capability is unavailable. Reconnect to an owner that reports supported pressure inputs.'};
 if(!value||typeof value.whole_skin!=='boolean'||!Array.isArray(value.regional_ids)||value.regional_ids.some(id=>typeof id!=='string'||!id)||new Set(value.regional_ids).size!==value.regional_ids.length||value.whole_skin&&value.regional_ids.length)return unavailable;
 if(value.whole_skin)return {topology:'whole',regional_ids:[],message:'Whole-skin pressure boundary is available.'};
 if(value.regional_ids.length)return {topology:'regional',regional_ids:[...value.regional_ids],message:'Whole-skin pressure is unavailable for this regional owner. Reset whole-skin pressure to 0; use named regional boundaries through the regional input interface.'};
 return unavailable;
}
export function bodyCommand(kind,sequence,forces,inputs={},frame=null) {
 bodyEndpoint(kind);
 if(!Number.isInteger(sequence)||sequence<0)throw Error('Missing current body sequence');
 const command={seconds:.02,sequence,forces};
 if(kind==='embodied') {
  const descending={...inputs.descending},pressure=inputs.skin_compression_pa??0;
  if(!Number.isFinite(pressure)||pressure<0||pressure>5000||Object.values(descending).some(v=>!Number.isFinite(v)||v<0||v>1))throw Error('Invalid body input');
  Object.assign(command,{descending,sensory_blocks:[...(inputs.sensory_blocks||[])],motor_blocks:[...(inputs.motor_blocks||[])]});
  const capability=skinPressureCapability(frame),regional=inputs.regional_skin_pressures;
  if(pressure!==0&&capability.topology!=='whole')throw Error('Whole-skin pressure is unavailable. '+capability.message);
  if(regional!==undefined) {
   if(capability.topology!=='regional')throw Error('Regional skin pressure is unavailable for this owner.');
   if(!regional||typeof regional!=='object'||Array.isArray(regional)||Object.entries(regional).some(([id,v])=>!capability.regional_ids.includes(id)||!Number.isFinite(v)||v<0||v>5000))throw Error('Invalid regional skin pressure boundary');
   command.regional_skin_pressures={...regional};
  } else if(capability.topology==='whole')command.skin_compression_pa=pressure;
 }
 return command;
}
export function signalInfo(frame,key) {
 const meta=frame?.physiology?.signal_metadata?.[key]||frame?.signal_metadata?.[key]||{};
 return {key,label:meta.label||systemicLabel(key),unit:meta.unit??frame?.physiology?.units?.[key]??'unit unavailable',
  owner:meta.owner||'native physiology',value:Number.isFinite(frame?.physiology?.values?.[key])?frame.physiology.values[key]:null};
}
export class LiveBodyHistory {
 constructor(limit=600){if(!Number.isInteger(limit)||limit<2||limit>10000)throw Error('Invalid live sample budget');this.limit=limit;this.clear();}
 clear(){this.owner=null;this.samples=[];this.sequence=-1;this.frame=null;}
 push(owner,frame) {
  if(!owner||frame?.schema!=='ihm.embodied-frame.v1'||!Number.isFinite(frame.time_s)||frame.time_s<0||!Number.isInteger(frame.sequence)||frame.sequence<0||!frame.physiology?.values)throw Error('Invalid live body frame');
  if(this.owner!==owner){this.clear();this.owner=owner;}
  if(frame.sequence===this.sequence){if(frame.time_s!==this.frame.time_s)throw Error('Inconsistent duplicate body sequence');return false;}
  if(frame.sequence<this.sequence||this.samples.length&&frame.time_s<this.samples.at(-1).time_s)throw Error('Live body clock reversed');
  if(this.frame&&frame.time_s===this.frame.time_s){this.sequence=frame.sequence;this.frame=frame;return false;}
  const values=Object.fromEntries(Object.entries(frame.physiology.values).map(([k,v])=>[k,Number.isFinite(v)?v:null]));
  this.samples.push({time_s:frame.time_s,values});if(this.samples.length>this.limit)this.samples.shift();
  this.sequence=frame.sequence;this.frame=frame;return true;
 }
 series(key){return {time_s:this.samples.map(s=>s.time_s),values:this.samples.map(s=>s.values[key]??null)};}
}
export function frameScope(frame) {
 if(frame?.schema==='ihm.embodied-frame.v1')return [
  'Live articulated muscles ↔ pinned IBM/controller ↔ native physiology. Generic source registration and decoder calibration remain incomplete.',
  frame.mechanics?.body_environment?.scope,
  frame.environment_state ? `Environment: ${frame.environment_state.contact_count} active body contacts. ${frame.environment_state.scope}` : 'Native supports are computational contact proxies; no scene contact owner is attached.',...(frame.mechanics?.limitations||[]),frame.coupling?.metabolic_law,
 ].filter(Boolean).join(' ');
 const s=frame?.scope||{};return [s.body_mechanics,s.body_rotations,s.body_gravity,s.body_environment,
  s.body_object_contact===false?'Body–object contact is not coupled.':null,
  s.clothing_contact===false?'Clothing contact is not coupled to this experiment.':null,
  s.physiology_feedback===false?'Reduced experiment: forces do not feed physiology.':null].filter(Boolean).join(' ');
}
const identity=[[1,0,0],[0,1,0],[0,0,1]];
function inverse(m){
 if(m?.length!==3||m.some(r=>r.length!==3||r.some(v=>!Number.isFinite(v))))throw Error('Invalid material transform');
 const [[a,b,c],[d,e,f],[g,h,i]]=m,det=a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
 if(!Number.isFinite(det)||det<=1e-12)throw Error('Singular material transform');
 return [[e*i-f*h,c*h-b*i,b*f-c*e],[f*g-d*i,a*i-c*g,c*d-a*f],[d*h-e*g,b*g-a*h,a*e-b*d]].map(row=>row.map(x=>x/det));
}
export function materialOffset(point,entity) {
 const delta=point.map((v,i)=>v-entity.centroid_m[i]);
 return rotateOffset(rotateOffset(delta,entity.rotation_matrix||identity,true),inverse(entity.deformation_gradient||identity));
}
export function materialPoint(offset,entity) {
 return rotateOffset(rotateOffset(offset,entity.deformation_gradient||identity),entity.rotation_matrix||identity).map((v,i)=>v+entity.centroid_m[i]);
}
export async function createBodyOwner(request,path,configuration) {
 const before=await request(path),previous=new Set((before.sessions||[]).map(s=>s.id));
 try{return await request(path,configuration);}
 catch(error){
  if(error.httpStatus>=400&&error.httpStatus<500)throw error;
  const after=await request(path),fresh=(after.sessions||[]).filter(s=>!previous.has(s.id)&&!s.closed);
  if(fresh.length!==1)throw error;
  return await request(path+'/'+fresh[0].id);
 }
}
export async function closeBodyOwner(request,path,wait) {
 let result=await request(path+'/close',{});
 while(result.closed!==true){
  if(result.error)throw Error(result.error);
  if(!['closing','initializing','ready'].includes(result.status))throw Error('Body cleanup outcome is unknown');
  await wait();result=await request(path);
 }
 return result;
}

// Recover an uncertain intake POST by reading the owner, never replaying events.
export async function scheduleBodyIntakes(request,path,sequence,events) {
 if(!Number.isInteger(sequence)||sequence<0||!Array.isArray(events)||!events.length)throw Error('Missing current body sequence or intake events');
 try{return {frame:await request(path+'/intakes',{sequence,events}),recovered:false};}
 catch(error){
  // The intake API guarantees 400 validation failures occur before mutation.
  // Refresh the owning sequence, but leave this proven rejection correctable.
  if(error.httpStatus===400)error.definitelyRejected=true;
  try{return {frame:await request(path),recovered:true,error};}catch{throw error;}
 }
}
