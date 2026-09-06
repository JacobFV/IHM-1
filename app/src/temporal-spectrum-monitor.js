import {LiveBodyHistory,signalInfo} from './embodied-live.js';

const finite=Number.isFinite;
const unavailable=reason=>({available:false,reason});
// Bounded quadrature, no FFT, resampling, tail extrapolation or pole fitting.
export function finiteWindowSpectrum(series,{unit,window_s=10,sigma=0,bins=65,maxFrequency=5}={}) {
 if(!finite(window_s)||window_s<=0||window_s>60||!finite(sigma)||sigma<0||sigma>10||!Number.isInteger(bins)||bins<2||bins>129||!finite(maxFrequency)||maxFrequency<=0||maxFrequency>100)throw Error('Invalid finite-spectrum settings');
 const times=series?.time_s,values=series?.values;
 if(!Array.isArray(times)||!Array.isArray(values)||times.length!==values.length||times.length>512)throw Error('Invalid or oversized spectrum samples');
 if(typeof unit!=='string'||!unit.trim()||unit==='unit unavailable')return unavailable('Signal unit unavailable.');
 if(times.length<4)return unavailable('Need at least four live samples.');
 if(!times.every(finite)||times.some((t,i)=>i&&t<=times[i-1]))return unavailable('Recorded clock must strictly increase.');
 const end=times.at(-1),first=times.findIndex(t=>t>=end-window_s),t=times.slice(first),x=values.slice(first);
 if(t.length<4)return unavailable('Need at least four samples in the selected window.');
 if(!x.every(finite))return unavailable('Missing signal samples in the selected window.');
 const dt=(t.at(-1)-t[0])/(t.length-1),steps=t.slice(1).map((v,i)=>v-t[i]);
 if(steps.some(v=>Math.abs(v-dt)>Math.max(1e-9,dt*1e-5)))return unavailable('Irregular sample clock; spectrum unavailable without resampling.');
 const top=Math.min(maxFrequency,.5/Math.max(...steps));
 const mean=x.reduce((sum,v)=>sum+v/x.length,0),elapsed=t.map(v=>v-t[0]);
 const frequency_hz=[],real=[],imag=[],magnitude=[];
 for(let k=0;k<bins;k++){
  const f=top*k/(bins-1);let re=0,im=0,priorRe=0,priorIm=0;
  for(let i=0;i<x.length;i++){
   const weight=(x[i]-mean)*Math.exp(-sigma*elapsed[i]),angle=2*Math.PI*f*elapsed[i],r=weight*Math.cos(angle),q=-weight*Math.sin(angle);
   if(i){re+=(priorRe+r)*steps[i-1]/2;im+=(priorIm+q)*steps[i-1]/2;}priorRe=r;priorIm=q;
  }
  if(!finite(re)||!finite(im)||!finite(Math.hypot(re,im)))return unavailable('Transform exceeds finite numerical range.');
  frequency_hz.push(f);real.push(re);imag.push(im);magnitude.push(Math.hypot(re,im));
 }
 return {available:true,frequency_hz,real,imag,magnitude,unit:`${unit}·s`,time_interval_s:[t[0],end],samples:t.length,
  sample_interval_s:dt,sigma_per_s:sigma,removed_sample_mean:mean,detrend:'subtract retained-window sample mean',
  interpretation:'Finite-window transform, not physiological poles or a causal transfer model.'};
}

export function mountTemporalSpectrumMonitor(host,{document:doc=globalThis.document}={}){
 if(!host||!doc)throw Error('Spectrum monitor requires a host and document');
 const history=new LiveBodyHistory(512),el=(tag,text)=>{const n=doc.createElement(tag);if(text)n.textContent=text;return n;};
 let disposed=false,frame=null,signature='',lastRender=-Infinity;
 const signal=el('select'),window=el('select'),damping=el('select');
 signal.setAttribute('aria-label','Spectrum signal');window.setAttribute('aria-label','Spectrum window seconds');damping.setAttribute('aria-label','Spectrum damping per second');
 for(const value of [2,5,10]){const o=el('option',`${value} s window`);o.value=String(value);window.append(o);}window.value='10';
 for(const value of [0,.1,1]){const o=el('option',`${value} s⁻¹ damping`);o.value=String(value);damping.append(o);}damping.value='0';
 const status=el('p','No live samples.'),plot=el('div'),note=el('p','Finite-window transform, not physiological poles or a causal transfer model. Retained-window mean removed; elapsed seconds from its first sample. 512 sample budget; no source anti-alias filter established. Float64 quadrature; tiny products may underflow.');
 status.setAttribute('role','status');host.replaceChildren(signal,window,damping,status,plot,note);
 function render(){
  plot.replaceChildren();if(!frame){status.textContent='No live samples.';return;}
  const info=signalInfo(frame,signal.value),r=finiteWindowSpectrum(history.series(signal.value),{unit:info.unit,window_s:Number(window.value),sigma:Number(damping.value)});
  if(!r.available){status.textContent=r.reason;return;}
  status.textContent=`${info.label} · ${r.unit} · ${r.samples} samples · ${r.time_interval_s[0].toFixed(3)}–${r.time_interval_s[1].toFixed(3)} s · σ ${r.sigma_per_s} s⁻¹ · 0–${r.frequency_hz.at(-1).toPrecision(3)} Hz`;
  const svg=doc.createElementNS('http://www.w3.org/2000/svg','svg'),line=doc.createElementNS('http://www.w3.org/2000/svg','polyline'),peak=Math.max(...r.magnitude);
  svg.setAttribute('viewBox','0 0 320 120');svg.setAttribute('role','img');svg.setAttribute('aria-label',`Finite Laplace magnitude; linear scale 0 to ${peak.toPrecision(3)} ${r.unit}; frequency 0 to ${r.frequency_hz.at(-1).toPrecision(3)} Hz`);
  line.setAttribute('points',r.magnitude.map((v,i)=>`${10+300*i/(r.magnitude.length-1)},${110-(peak?100*v/peak:0)}`).join(' '));line.setAttribute('fill','none');line.setAttribute('stroke','currentColor');line.setAttribute('stroke-width','2');svg.append(line);plot.append(svg,el('p',`Linear magnitude: 0–${peak.toPrecision(3)} ${r.unit}. Frequency: 0–${r.frequency_hz.at(-1).toPrecision(3)} Hz.`));
 }
 for(const control of [signal,window,damping])control.onchange=()=>{if(!disposed)render();};
 function clear(){history.clear();frame=null;signature='';lastRender=-Infinity;signal.replaceChildren();render();}
 return {clear,dispose(){clear();disposed=true;for(const c of [signal,window,damping])c.disabled=true;},update(owner,next){
  if(disposed)return;
  const keys=Object.keys(next?.physiology?.values||{}),sig=JSON.stringify(keys.map(key=>[key,signalInfo(next,key).unit]));
  if(history.owner!==owner||signature!==sig){clear();signature=sig;for(const key of keys){const o=el('option',signalInfo(next,key).label);o.value=key;signal.append(o);}signal.value=keys[0]||'';}
  const added=history.push(owner,next);frame=next;
  if(added&&(next.time_s-lastRender>=.5||history.samples.length<=4)){lastRender=next.time_s;render();}
 }};
}
