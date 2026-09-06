// Queue configuration and authoritative receipts only; no delivery simulation.
export const MEAL_FIELDS = Object.freeze([
 ['carbohydrate_g','Carbohydrate · g'],['protein_g','Protein · g'],['fat_g','Fat · g'],
 ['sodium_g','Sodium · g'],['calcium_mg','Calcium · mg'],['water_ml','Water · mL'],
]);
const bounded=(n,max)=>typeof n==='number'&&Number.isFinite(n)&&n>=0&&n<=max;
export function validateIntakeEvent(event){
 if(!event||typeof event.event_id!=='string'||!(/^[A-Za-z0-9_-]{1,64}$/).test(event.event_id))throw Error('Event ID must use 1–64 letters, digits, underscores or hyphens.');
 if(!bounded(event.time_s,86400)||Math.abs(event.time_s*50-Math.round(event.time_s*50))>1e-7)throw Error('Time must be 0–86400 seconds in 0.02 second steps.');
 const meal=event.meal;
 if(!meal||typeof meal.name!=='string'||!(/^[A-Za-z0-9_]{1,64}$/).test(meal.name))throw Error('Meal name must use 1–64 letters, digits or underscores.');
 const keys=['name',...MEAL_FIELDS.map(([key])=>key)];
 if(Object.keys(meal).some(key=>!keys.includes(key)))throw Error('Unsupported native Meal field.');
 const result={name:meal.name};
 for(const [key,label] of MEAL_FIELDS){const amount=meal[key]===undefined?0:meal[key];if(!bounded(amount,10000))throw Error(`${label} must be a finite amount from 0 to 10000.`);result[key]=amount;}
 if(!MEAL_FIELDS.some(([key])=>result[key]>0))throw Error('Meal must contain a positive amount.');
 return {event_id:event.event_id,time_s:event.time_s,meal:result};
}
export function mountIntakeMonitor(host,{submit,document:doc=globalThis.document}={}){
 if(!host||typeof submit!=='function')throw Error('Intake monitor requires a host and submit callback.');
 const el=(tag,text)=>{const n=doc.createElement(tag);if(text)n.textContent=text;return n;};
 let connected=false,busy=false,locked=false,disposed=false,records=new Map();
 const form=el('form'),list=el('ol'),status=el('p','No live intake connection.');status.setAttribute('role','status');
 const fields={};
 function input(name,label,value,type='number'){
  const row=el('label',label),control=el('input');control.name=name;control.type=type;control.value=value;control.required=true;control.setAttribute('aria-label',label);
  if(type==='number'){control.min='0';control.max=name==='time_s'?'86400':'10000';control.step=name==='time_s'?'0.02':'any';}else control.maxLength=64;
  row.append(control);form.append(row);fields[name]=control;
 }
 input('event_id','Event ID','','text');input('time_s','Elapsed body time · s','0');input('name','Meal name','mixed_meal','text');
 for(const [key,label] of MEAL_FIELDS)input(key,label,'0');
 const button=el('button','Queue intake');button.type='submit';form.append(button);
 host.replaceChildren(el('p','Schedule food or water using native Meal quantities. Acceptance acknowledges the command; absorption remains native physiology.'),form,status,list);
 function controls(){for(const control of [...Object.values(fields),button])control.disabled=disposed||!connected||busy||locked;}
 function update(snapshot){
  if(disposed)return;
  if(!snapshot||!Array.isArray(snapshot.events)||snapshot.events.length>10000)throw Error('Invalid intake snapshot.');
  const next=new Map();
  for(const row of snapshot.events){
   const event=validateIntakeEvent(row);if(next.has(event.event_id))throw Error('Duplicate snapshot event ID.');
   if(!['queued','issued','accepted','uncertain'].includes(row.state))throw Error('Invalid intake state.');
   const prior=records.get(event.event_id);
   const permitted={queued:['queued','issued','accepted','uncertain'],issued:['issued','accepted','uncertain'],accepted:['accepted'],uncertain:['uncertain']};
   if(prior&&(!permitted[prior.state].includes(row.state)||JSON.stringify(validateIntakeEvent(prior))!==JSON.stringify(event)))throw Error('Intake identity or state cannot regress.');
   for(const key of ['issued_tick','native_sequence'])if(row[key]!=null&&(!Number.isSafeInteger(row[key])||row[key]<0))throw Error('Invalid receipt integer.');
   if(row.reason!=null&&(typeof row.reason!=='string'||row.reason.length>1024))throw Error('Invalid receipt reason.');
   next.set(event.event_id,{...event,state:row.state,issued_tick:row.issued_tick,native_sequence:row.native_sequence,reason:row.reason});
  }
  for(const id of records.keys())if(!next.has(id))throw Error('Intake snapshot cannot discard existing identities.');
  records=next;locked=locked||[...records.values()].some(row=>row.state==='uncertain');
  list.replaceChildren();
  for(const row of records.values()){
   const amounts=MEAL_FIELDS.filter(([key])=>row.meal[key]>0).map(([key,label])=>`${label}: ${row.meal[key]}`).join(', ');
   const receipt=row.issued_tick==null?'':` · issued at ${(row.issued_tick*.02).toFixed(2)} s`;
   list.append(el('li',`${row.event_id} · ${row.state} · scheduled ${row.time_s.toFixed(2)} s · ${row.meal.name} · ${amounts}${receipt}${row.native_sequence==null?'':` · native sequence ${row.native_sequence}`}${row.reason?` · ${row.reason}`:''}`));
  }
  if(locked)status.textContent='Intake outcome uncertain. Further requests are disabled; no retry is available.';
  else if(!busy)status.textContent=connected?'Intake schedule synchronized.':'No live intake connection.';
  controls();
 }
 form.onsubmit=async e=>{
  e.preventDefault();if(disposed||!connected||busy||locked)return;
  let event;
  try{
   const numeric=key=>fields[key].value.trim()===''?NaN:Number(fields[key].value);
   event=validateIntakeEvent({event_id:fields.event_id.value,time_s:numeric('time_s'),meal:{name:fields.name.value,...Object.fromEntries(MEAL_FIELDS.map(([key])=>[key,numeric(key)]))}});
   if(records.has(event.event_id))throw Error('Event ID already exists.');
  }catch(error){status.textContent=error.message;return;}
  busy=true;controls();status.textContent='Submitting schedule request; native delivery is not yet established.';
  try{const snapshot=await submit([event]);if(disposed)return;
   if(!snapshot?.events?.some(row=>JSON.stringify(validateIntakeEvent(row))===JSON.stringify(event)))throw Error('Response does not contain the submitted event.');
   update(snapshot);status.textContent=locked?'Intake outcome uncertain. Further requests are disabled; no retry is available.':'Schedule request acknowledged. Native command status is listed below.';}
  catch(error){if(disposed)return;if(error?.definitelyRejected===true){status.textContent=`Schedule request rejected: ${error.message}`;}else{locked=true;status.textContent=`Schedule request outcome uncertain. No retry is available. ${error?.message||'No valid response.'}`;}}
  finally{busy=false;controls();}
 };
 controls();
 return {update,dispose(){disposed=true;connected=false;controls();},setConnected(value){if(disposed)return;connected=value===true;if(!busy&&!locked)status.textContent=connected?'Ready to queue intake.':'No live intake connection.';controls();}};
}
