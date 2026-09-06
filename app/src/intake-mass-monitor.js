const finite=value=>typeof value==='number'&&Number.isFinite(value)?value:null;
const mass=value=>value===null?'Unavailable':`${Number(value).toPrecision(6)} kg`;
const boundary=value=>value===null?'None reported':value===undefined?'Unavailable':mass(finite(value.mass_kg));
export function intakeMassSnapshot(frame){
 if(frame?.schema!=='ihm.embodied-frame.v1')return null;
 const audit=frame.intake_mass,enabled=audit?.enabled===true,bridge=enabled?audit.bridge:null;
 return {mass_kg:finite(frame.mechanics?.effective_native_body_mass_kg),enabled:audit?.enabled===false?false:enabled?true:null,
  binding:enabled?audit.binding:null,failed:enabled&&bridge?.failed===true,applied_mass_kg:finite(bridge?.applied_mass_kg),capacity_kg:finite(bridge?.capacity_kg),
  pending_boundary:enabled?bridge?.pending_boundary:undefined,last_boundary:enabled?bridge?.last_boundary:undefined,
  mechanical_sequence:bridge?.mechanical_sequence,mechanical_reference_id:bridge?.mechanical_reference_id,checkpoint_scope:bridge?.checkpoint_scope};
}
export function mountIntakeMassMonitor(host,{document:doc=globalThis.document}={}){
 const el=(tag,text)=>{const node=doc.createElement(tag);if(text)node.textContent=text;return node;};
 const status=el('p'),clock=el('p'),readouts=el('dl'),details=el('details'),source=el('div');status.setAttribute('role','status');
 readouts.className='live-vitals';clock.className='muted';details.append(el('summary','Mass owner & assumptions'),source);host.replaceChildren(status,clock,readouts,details);
 let owner=null,signature='',disposed=false;
 function clear(){owner=null;signature='';status.textContent='No live body mass frame.';clock.textContent='';readouts.replaceChildren();source.replaceChildren();}
 function row(label,value){readouts.append(el('dt',label),el('dd',value));}
 clear();
 return {clear,dispose(){clear();disposed=true;},update(nextOwner,frame){
  if(disposed)return;const snapshot=intakeMassSnapshot(frame);if(!snapshot||!nextOwner){clear();return;}
  if(owner!==nextOwner){clear();owner=nextOwner;}
  clock.textContent=`Body time ${finite(frame.time_s)===null?'unavailable':frame.time_s.toFixed(2)+' s'} · sequence ${Number.isInteger(frame.sequence)?frame.sequence:'unavailable'}`;
  const next=JSON.stringify(snapshot);if(signature===next)return;signature=next;
  status.textContent=snapshot.failed?'Intake mass transfer failed; pending mass is unconfirmed.':snapshot.enabled===true?'Intake mass coupling enabled.':snapshot.enabled===false?'Intake mass coupling disabled for this body.':'Intake mass coupling unavailable in this frame.';
  readouts.replaceChildren();row('Effective native body mass',mass(snapshot.mass_kg));
  row('Applied intake mass',mass(snapshot.applied_mass_kg));row('Validated intake capacity',mass(snapshot.capacity_kg));
  row('Pending consumed boundary',boundary(snapshot.pending_boundary));row('Last consumed boundary',boundary(snapshot.last_boundary));
  source.replaceChildren(el('p','Values come from the accepted native mechanical frame and intake bridge receipts. Scheduled meals do not establish applied mass.'));
  if(snapshot.enabled){
   const binding=snapshot.binding||{};source.append(el('p',`Native body: ${binding.body||'unavailable'} · canonical binding: ${binding.canonical_entity_id||'unavailable'}.`));
   const station=binding.station_source_m;source.append(el('p',`Binding station · source coordinates: ${Array.isArray(station)&&station.length===3&&station.every(v=>finite(v)!==null)?station.map(v=>Number(v).toPrecision(5)).join(', ')+' m':'unavailable'}.`));
   if(binding.scope)source.append(el('p',String(binding.scope)));
   if(binding.incoming_velocity_basis)source.append(el('p',`Incoming velocity basis: ${String(binding.incoming_velocity_basis).replaceAll('_',' ')}.`));
   for(const [label,value] of [['Last',snapshot.last_boundary],['Pending',snapshot.pending_boundary]])if(value){
    source.append(el('p',`${label} boundary · ticks ${Number.isInteger(value.interval_start_tick)?value.interval_start_tick:'unavailable'}–${Number.isInteger(value.interval_end_tick)?value.interval_end_tick:'unavailable'} · mechanical transfer ${value.mechanical_transfer_applied===true?'confirmed':value.mechanical_transfer_applied===false?'not applied':'unconfirmed'} · boundary accounted ${value.boundary_accounted===true?'yes':value.boundary_accounted===false?'no':'unconfirmed'}.`));
   }
   source.append(el('p',`Mechanical transfer sequence: ${Number.isInteger(snapshot.mechanical_sequence)?snapshot.mechanical_sequence:'unavailable'}.`));
   if(snapshot.checkpoint_scope)source.append(el('p',String(snapshot.checkpoint_scope)));
   const provenance=el('p',`Canonical reference: ${binding.canonical_sha256||'unavailable'} · registration: ${binding.registration_sha256||'unavailable'} · native mass reference: ${snapshot.mechanical_reference_id||'unavailable'}.`);provenance.style.overflowWrap='anywhere';source.append(provenance);
  }
 }};
}
