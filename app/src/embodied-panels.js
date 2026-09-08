import {LiveBodyHistory,signalInfo,skinPressureCapability} from './embodied-live.js';
const $=id=>document.getElementById(id);
const el=(tag,text,cls)=>{const n=document.createElement(tag);if(text)n.textContent=text;if(cls)n.className=cls;return n;};
const format=v=>Number.isFinite(v)?Number(v).toPrecision(5):'Unavailable';
export class MotorInputs {
 constructor(){this.clear();}
 clear(){this.keys=new Set();this.drives={};this.sensory=new Set();this.motor=new Set();this.pressure=0;}
 bind(keys){this.keys=new Set(keys);for(const key of Object.keys(this.drives))if(!this.keys.has(key))delete this.drives[key];for(const set of [this.sensory,this.motor])for(const key of set)if(!this.keys.has(key))set.delete(key);}
 set(key,drive,sensory,motor){if(!this.keys.has(key)||!Number.isFinite(drive)||drive<0||drive>1||typeof sensory!=='boolean'||typeof motor!=='boolean')throw Error('Invalid selected muscle input');if(drive)this.drives[key]=drive;else delete this.drives[key];sensory?this.sensory.add(key):this.sensory.delete(key);motor?this.motor.add(key):this.motor.delete(key);}
 snapshot(){return {descending:{...this.drives},sensory_blocks:[...this.sensory],motor_blocks:[...this.motor],skin_compression_pa:this.pressure};}
}
export function mountEmbodiedPanels() {
 const history=new LiveBodyHistory(600),inputs=new MotorInputs();let frame=null,lastDraw=-Infinity,signature='',muscleSignature='';
 const live=$('live-body-monitor'),motor=$('live-motor-monitor');
 live.replaceChildren(el('p','Start Body to initialize a unified live state.','live-body-status'),el('dl',null,'live-vitals'),el('p','Source calibration and integrated homeostasis remain under verification.','muted'));
 const readouts=live.querySelector('dl');
 const select=el('select');select.setAttribute('aria-label','Live muscle effector');
 const search=el('input');search.type='search';search.placeholder='Find muscle…';search.setAttribute('aria-label','Search live muscles');
 const drive=el('input');drive.type='number';drive.min='0';drive.max='1';drive.step='.05';drive.value='0';drive.setAttribute('aria-label','Requested descending muscle drive');
 const sensory=el('input');sensory.type='checkbox';const block=el('input');block.type='checkbox';
 const pressure=el('input');pressure.type='number';pressure.min='0';pressure.max='5000';pressure.step='10';pressure.value='0';pressure.setAttribute('aria-label','Whole skin compression in pascals');
 const row=(text,input)=>{const label=el('label',text,'live-input-label');label.append(input);return label;};
 const muscleState=el('p','Actual muscle activation and force appear after initialization.','muted');
 const note=el('p','Inputs are held for the next live body tick. Paused bodies apply changes on resume.','muted');
 const release=el('button','Release all drives, blocks & pressure');release.type='button';
 const pressureStatus=el('p',null,'muted');pressureStatus.setAttribute('role','status');
 const pressureReset=el('button','Clear pending whole-skin pressure');pressureReset.type='button';
 function syncPressure(){
  const capability=skinPressureCapability(frame);pressure.disabled=capability.topology!=='whole';
  pressureReset.disabled=inputs.pressure===0;
  pressureStatus.textContent=capability.topology==='whole'?'Uniform native Skin compartment pressure; not local pressure beneath the cursor.':capability.message;
  if(inputs.pressure!==0&&capability.topology!=='whole')pressureStatus.textContent+=` Pending ${inputs.pressure} Pa is preserved and blocks advance. Clear pending whole-skin pressure explicitly to continue.`;
 }
 pressureReset.onclick=()=>{inputs.pressure=0;pressure.value=0;syncPressure();note.textContent='Pending whole-skin pressure cleared. Other requested inputs are unchanged.';};
 motor.replaceChildren(search,select,row('Requested descending drive · 0–1',drive),row('Block selected sensory pathway',sensory),row('Block selected motor pathway',block),muscleState,
  row('Whole-Skin compression · Pa',pressure),pressureStatus,pressureReset,release,note);
 const syncControls=()=>{const key=select.value;drive.value=inputs.drives[key]??0;sensory.checked=inputs.sensory.has(key);block.checked=inputs.motor.has(key);};
 const save=()=>{try{inputs.set(select.value,drive.value===''?NaN:Number(drive.value),sensory.checked,block.checked);note.textContent='Requested inputs updated · applied on the next live body tick.';}catch(error){note.textContent=error.message;syncControls();}};
 drive.onchange=save;sensory.onchange=save;block.onchange=save;select.onchange=()=>{syncControls();drawMuscle();};
 pressure.onchange=()=>{if(skinPressureCapability(frame).topology!=='whole'){pressure.value=inputs.pressure;syncPressure();return;}const p=pressure.value===''?NaN:Number(pressure.value);if(!Number.isFinite(p)||p<0||p>5000){pressure.value=inputs.pressure;note.textContent='Skin pressure must be between 0 and 5000 Pa.';}else{inputs.pressure=p;note.textContent='Whole-Skin pressure requested for the next live body tick.';syncPressure();}};
 release.onclick=()=>{const keys=[...inputs.keys];inputs.clear();inputs.bind(keys);syncControls();pressure.value=0;syncPressure();note.textContent='Release requested · muscle activation still decays through its native law.';};
 function options(){const current=select.value;select.replaceChildren();for(const key of [...inputs.keys].sort())if(key.toLowerCase().includes(search.value.trim().toLowerCase())){const o=el('option',key);o.value=key;select.append(o);}if([...select.options].some(o=>o.value===current))select.value=current;syncControls();drawMuscle();}
 search.oninput=options;
 function drawMuscle(){const m=frame?.mechanics?.muscles?.[select.value];muscleState.textContent=m?`${select.value} · activation ${format(m.activation)} · excitation ${format(m.excitation)} · tendon ${format(m.tendon_force_n)} N · fiber ${format(m.fiber_length_m)} m. ${m.sensor_basis||'Native muscle state'}`:'No matching live muscle selected.';}
 const graphDefaults=['heart_rate_per_min','arterial_co2_mmhg','core_temperature_c'];
 const graphs=graphDefaults.map((preferred,i)=>{
  const host=$(`live-signal-${i+1}-monitor`),chooser=el('select');chooser.setAttribute('aria-label',`Live signal ${i+1}`);
  const value=el('p','Unavailable','live-signal-value'),plot=el('div',null,'live-signal-plot'),axis=el('p','No live samples.','muted');
  host.replaceChildren(chooser,value,plot,axis);chooser.onchange=()=>drawGraph({chooser,value,plot,axis});
  return {chooser,value,plot,axis,preferred};
 });
 function drawGraph(g){
  const info=signalInfo(frame,g.chooser.value),s=history.series(g.chooser.value);g.value.textContent=`${format(info.value)} ${info.unit}`;
  g.value.title=`${info.label} · ${info.owner}`;const good=s.values.filter(Number.isFinite);
  if(s.time_s.length<2||!good.length){g.plot.replaceChildren();g.axis.textContent='Waiting for two actual live samples.';return;}
  let low=Math.min(...good),high=Math.max(...good);const span=Math.max(high-low,Math.abs(high)*1e-6,1e-9);const dt=s.time_s.at(-1)-s.time_s[0];
  let d='',pen=false;for(let i=0;i<s.values.length;i++){if(!Number.isFinite(s.values[i])){pen=false;continue;}const x=400*(s.time_s[i]-s.time_s[0])/dt,y=82-74*(s.values[i]-low)/span;d+=`${pen?'L':'M'}${x.toFixed(2)},${y.toFixed(2)} `;pen=true;}
  const svg=document.createElementNS('http://www.w3.org/2000/svg','svg'),path=document.createElementNS('http://www.w3.org/2000/svg','path');
  svg.setAttribute('viewBox','0 0 400 90');svg.setAttribute('role','img');svg.setAttribute('aria-label',`${info.label}, live body samples`);path.setAttribute('d',d);path.setAttribute('fill','none');path.setAttribute('stroke','currentColor');path.setAttribute('stroke-width','1.5');svg.append(path);g.plot.replaceChildren(svg);
  g.axis.textContent=`Live t ${s.time_s[0].toFixed(2)}–${s.time_s.at(-1).toFixed(2)} s · ${s.values.length} stored samples · range ${format(low)}–${format(high)} ${info.unit}`;
 }
 function disabled(value){for(const control of motor.querySelectorAll('input,select,button'))control.disabled=value;}
 syncPressure();disabled(true);
 return {
  get inputs(){return inputs.snapshot();},
  update(owner,next){
   if(next?.schema!=='ihm.embodied-frame.v1'){this.clear();return;}
   if(history.owner&&history.owner!==owner)this.clear();
   history.push(owner,next);frame=next;disabled(false);syncPressure();
   const keys=Object.keys(next.physiology.values).sort(),sig=keys.join('|');
   if(signature!==sig){signature=sig;for(const g of graphs){const old=g.chooser.value;g.chooser.replaceChildren();for(const key of keys){const o=el('option',signalInfo(next,key).label);o.value=key;g.chooser.append(o);}g.chooser.value=keys.includes(old)?old:keys.includes(g.preferred)?g.preferred:keys[0]||'';}}
   const muscles=Object.keys(next.mechanics?.muscles||{}).sort(),ms=muscles.join('|');if(ms!==muscleSignature){muscleSignature=ms;inputs.bind(muscles);options();}
   if(next.time_s-lastDraw<.2&&lastDraw>=0)return;lastDraw=next.time_s;
   live.querySelector('.live-body-status').textContent=`Live body · t ${next.time_s.toFixed(2)} s · sequence ${next.sequence} · ${muscles.length} muscle effectors`;
   readouts.replaceChildren();for(const key of ['heart_rate_per_min','mean_arterial_pressure_mmhg','oxygen_saturation','core_temperature_c','arterial_co2_mmhg','lung_volume_ml']){
    const info=signalInfo(next,key);readouts.append(el('dt',info.label),el('dd',`${format(info.value)} ${info.unit}`));
   }
   drawMuscle();for(const g of graphs)drawGraph(g);
  },
  status(text){live.querySelector('.live-body-status').textContent=text;},
  clear(){frame=null;history.clear();inputs.clear();signature='';muscleSignature='';lastDraw=-Infinity;select.replaceChildren();search.value='';drive.value='0';sensory.checked=false;block.checked=false;pressure.value='0';syncPressure();disabled(true);readouts.replaceChildren();live.querySelector('.live-body-status').textContent='No live body owner. Recorded experiments are separate.';for(const g of graphs){g.chooser.replaceChildren();g.value.textContent='Unavailable';g.plot.replaceChildren();g.axis.textContent='No live samples.';}},
 };
}
