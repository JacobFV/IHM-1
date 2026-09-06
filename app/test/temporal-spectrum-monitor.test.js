import test from 'node:test';
import assert from 'node:assert/strict';
import {finiteWindowSpectrum,mountTemporalSpectrumMonitor} from '../src/temporal-spectrum-monitor.js';
const sine=(start=0)=>({time_s:Array.from({length:401},(_,i)=>start+i/20),values:Array.from({length:401},(_,i)=>Math.sin(2*Math.PI*.5*i/20))});
test('finite observed transform preserves clock origin, units, damping and analytic sine amplitude',()=>{
 const options={unit:'mV',window_s:20,sigma:0,bins:21,maxFrequency:1};
 const a=finiteWindowSpectrum(sine(),options),b=finiteWindowSpectrum(sine(100),options);
 assert.equal(a.available,true);assert.equal(a.unit,'mV·s');assert.deepEqual(a.time_interval_s,[0,20]);
 assert.ok(Math.abs(a.imag[10]+10)<1e-10);assert.ok(Math.abs(a.real[10])<1e-10);
 a.magnitude.forEach((v,i)=>assert.ok(Math.abs(v-b.magnitude[i])<1e-10));
 const damped=finiteWindowSpectrum(sine(),{...options,sigma:1});assert.ok(damped.magnitude[10]<a.magnitude[10]);
 assert.equal(damped.sigma_per_s,1);assert.match(a.interpretation,/not.*poles/);
});
test('window truncation uses recorded endpoints and caps frequencies at sampled Nyquist',()=>{
 const r=finiteWindowSpectrum(sine(),{unit:'V',window_s:2,sigma:0,maxFrequency:100});
 assert.deepEqual(r.time_interval_s,[18,20]);assert.equal(r.samples,41);assert.ok(r.frequency_hz.at(-1)<=10);
 assert.equal(r.detrend,'subtract retained-window sample mean');
});
test('missing values, irregular clocks and unknown units never fabricate a spectrum',()=>{
 const s=sine();s.values[399]=null;assert.equal(finiteWindowSpectrum(s,{unit:'V'}).available,false);
 const t=sine();t.time_s[399]+=.001;assert.equal(finiteWindowSpectrum(t,{unit:'V'}).available,false);
 assert.equal(finiteWindowSpectrum(sine(),{unit:'unit unavailable'}).available,false);
 assert.throws(()=>finiteWindowSpectrum(sine(),{unit:'V',sigma:-1}));
 assert.throws(()=>finiteWindowSpectrum(sine(),{unit:'V',bins:100000}));
 const oversized={time_s:Array.from({length:513},(_,i)=>i),values:Array(513).fill(0)};
 assert.throws(()=>finiteWindowSpectrum(oversized,{unit:'V'}));
});
class Element{
 constructor(tag){this.tagName=tag;this.children=[];this.value='';this._text='';this.attributes={};}
 set textContent(t){this._text=String(t);this.children=[];}get textContent(){return this._text+this.children.map(c=>c.textContent).join('');}
 append(...n){this.children.push(...n);}replaceChildren(...n){this._text='';this.children=n;}
 setAttribute(k,v){this.attributes[k]=v;}
}
const frame=(sequence,unit='mV')=>({schema:'ihm.embodied-frame.v1',sequence,time_s:sequence*.02,physiology:{values:{ECG:Math.sin(sequence*.1)},units:{ECG:unit}}});
test('stream monitor bounds history and clears owner and unit changes without duplicate-time samples',()=>{
 const host=new Element('div'),doc={createElement:t=>new Element(t),createElementNS:(_,t)=>new Element(t)};
 const m=mountTemporalSpectrumMonitor(host,{document:doc});
 for(let i=0;i<550;i++)m.update('one',frame(i));
 assert.match(host.textContent,/512 sample budget/);assert.match(host.textContent,/mV·s/);
 const prior=host.textContent;m.update('one',frame(549));assert.equal(host.textContent,prior);
 assert.throws(()=>m.update('one',frame(548)),/reversed/);
 m.update('two',frame(0));assert.match(host.textContent,/at least four/i);assert.doesNotMatch(host.textContent,/mV·s/);
 for(let i=1;i<30;i++)m.update('two',frame(i));
 m.update('two',frame(30,'V'));assert.match(host.textContent,/at least four/i);
 m.clear();assert.match(host.textContent,/No live samples/);m.dispose();m.update('one',frame(500));assert.match(host.textContent,/No live samples/);
});
