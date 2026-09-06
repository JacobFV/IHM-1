import {LiveBodyHistory} from './embodied-live.js';
export const SKIN_VOLTAGE_KEY='tissue.Skin.bulk_membrane_potential_v';
export function skinVoltageMv(frame){
 const value=frame?.schema==='ihm.embodied-frame.v1'?frame.physiology?.values?.[SKIN_VOLTAGE_KEY]:null;
 return typeof value==='number'&&Number.isFinite(value)&&Number.isFinite(value*1000)?value*1000:null;
}
export function mountSkinVoltageMonitor(host,{document:doc=globalThis.document}={}){
 const el=(tag,text)=>{const n=doc.createElement(tag);if(text)n.textContent=text;return n;};
 const value=el('p','Unavailable'),status=el('p','No native Skin tissue voltage frame.'),plot=el('div'),axis=el('p'),details=el('details');
 value.className='live-signal-value';plot.className='live-signal-plot';axis.className='muted';status.setAttribute('role','status');
 details.append(el('summary','Native quantity & evidence'),el('p','BioGears SkinTissue intracellular potential relative to extracellular potential. This is bulk tissue Vm, not epidermal TEP or a measured keratinocyte voltage.'),el('p',`Read-only native field: ${SKIN_VOLTAGE_KEY}. Source volts are converted to millivolts (×1000). Missing observations remain unavailable; source-model output carries no implied empirical calibration.`));
 host.replaceChildren(value,status,plot,axis,details);
 const history=new LiveBodyHistory(300);let disposed=false,lastDraw=-Infinity;
 function clear(){history.clear();lastDraw=-Infinity;value.textContent='Unavailable';status.textContent='No native Skin tissue voltage frame.';plot.replaceChildren();axis.textContent='No live voltage samples.';}
 clear();
 return {clear,dispose(){clear();disposed=true;},update(owner,frame){
  if(disposed)return;if(!owner||frame?.schema!=='ihm.embodied-frame.v1'){clear();return;}
  if(history.owner&&history.owner!==owner)clear();
  const mv=skinVoltageMv(frame);history.push(owner,{...frame,physiology:{values:{[SKIN_VOLTAGE_KEY]:mv}}});
  value.textContent=mv===null?'Unavailable':`${mv.toPrecision(6)} mV`;
  status.textContent=mv===null?'Native bulk Skin tissue Vm is unavailable in this frame.':`Bulk Skin tissue Vm · native model · t ${frame.time_s.toFixed(2)} s`;
  if(mv===null){plot.replaceChildren();axis.textContent='No current native voltage observation.';lastDraw=-Infinity;return;}
  if(frame.time_s-lastDraw<.2)return;lastDraw=frame.time_s;
  const samples=history.series(SKIN_VOLTAGE_KEY),good=samples.values.filter(Number.isFinite),duration=samples.time_s.at(-1)-samples.time_s[0];
  if(good.length<2||!(duration>0)){plot.replaceChildren();axis.textContent='Waiting for two actual native voltage samples.';return;}
  const low=Math.min(...good),high=Math.max(...good),span=Math.max(high-low,1e-6);let path='',pen=false;
  for(let i=0;i<samples.values.length;i++){
   if(!Number.isFinite(samples.values[i])){pen=false;continue;}
   const x=400*(samples.time_s[i]-samples.time_s[0])/duration,y=80-70*(samples.values[i]-low)/span;
   path+=`${pen?'L':'M'}${x.toFixed(3)},${y.toFixed(3)} `;pen=true;
  }
  const svg=doc.createElementNS('http://www.w3.org/2000/svg','svg'),line=doc.createElementNS('http://www.w3.org/2000/svg','path');
  svg.setAttribute('viewBox','0 0 400 90');svg.setAttribute('role','img');svg.setAttribute('aria-label','Actual native bulk Skin tissue Vm samples');
  line.setAttribute('d',path);line.setAttribute('stroke','#8ec9bd');line.setAttribute('stroke-width','1.5');line.setAttribute('fill','none');svg.append(line);plot.replaceChildren(svg);
  axis.textContent=`${samples.time_s[0].toFixed(2)}–${samples.time_s.at(-1).toFixed(2)} s · ${good.length} finite samples · ${low.toPrecision(5)}–${high.toPrecision(5)} mV`;
 }};
}
