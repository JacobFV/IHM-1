import test from 'node:test';
import assert from 'node:assert/strict';
import {mountEmbodiedPanels} from '../src/embodied-panels.js';
// Minimal DOM surface used by this widget; no browser, geometry or network.
class Element {
 constructor(tag){this.tagName=tag;this.children=[];this.attributes={};this._value=undefined;this._text='';}
 set textContent(text){this._text=String(text);this.children=[];}
 get textContent(){return this._text+this.children.map(n=>n.textContent).join('');}
 append(...nodes){this.children.push(...nodes);}
 replaceChildren(...nodes){this._text='';this.children=[...nodes];if(this.tagName==='select')this._value=undefined;}
 setAttribute(key,value){this.attributes[key]=value;}
 get value(){return this._value??(this.tagName==='select'?this.children[0]?.value:'')??'';}
 set value(value){this._value=String(value);}
 get options(){return this.children;}
 querySelectorAll(selector){const choices=selector.split(',');const matches=n=>choices.some(s=>s.startsWith('.')?n.className===s.slice(1):n.tagName===s);return this.children.flatMap(n=>[...(matches(n)?[n]:[]),...n.querySelectorAll(selector)]);}
 querySelector(selector){return this.querySelectorAll(selector)[0]||null;}
}
test('real monitor event handlers preserve independent signals and selective body inputs',()=>{
 const previous=globalThis.document;const ids=['live-body-monitor','live-motor-monitor',...['1','2','3'].map(i=>`live-signal-${i}-monitor`)];const hosts=Object.fromEntries(ids.map(id=>[id,new Element('div')]));
 globalThis.document={getElementById:id=>hosts[id],createElement:t=>new Element(t),createElementNS:(_,t)=>new Element(t)};
 try{
  const panels=mountEmbodiedPanels();const native=(t,s)=>({schema:'ihm.embodied-frame.v1',time_s:t,sequence:s,input_capabilities:{skin_pressure:{whole_skin:true,regional_ids:[]}},entities:{},mechanics:{muscles:{bra_r:{activation:.1,tendon_force_n:3,fiber_length_m:.1}}},physiology:{values:{heart_rate_per_min:72+s,core_temperature_c:37+s},signal_metadata:{heart_rate_per_min:{unit:'1/min'},core_temperature_c:{unit:'C'}}}});
  panels.update('one',native(0,0));panels.update('one',native(.2,1));
  const first=hosts['live-signal-1-monitor'].querySelector('select'),second=hosts['live-signal-2-monitor'].querySelector('select');
  first.value='core_temperature_c';first.onchange();second.value='heart_rate_per_min';second.onchange();
  assert.match(hosts['live-signal-1-monitor'].querySelector('.live-signal-value').textContent,/38.* C/);
  assert.match(hosts['live-signal-2-monitor'].querySelector('.live-signal-value').textContent,/73.*1\/min/);
  const motor=hosts['live-motor-monitor'];const drive=motor.querySelectorAll('input').find(n=>n.attributes['aria-label']==='Requested descending muscle drive');
  drive.value='.4';drive.onchange();assert.deepEqual(panels.inputs.descending,{bra_r:.4});
  const check=motor.querySelectorAll('input').find(n=>n.type==='checkbox');check.checked=true;check.onchange();assert.deepEqual(panels.inputs.sensory_blocks,['bra_r']);
  const pressure=motor.querySelectorAll('input').find(n=>n.attributes['aria-label']==='Whole skin compression in pascals');pressure.value='800';pressure.onchange();assert.equal(panels.inputs.skin_compression_pa,800);
  const regional=native(.22,2);regional.input_capabilities.skin_pressure={whole_skin:false,regional_ids:['region_a','residual']};
  panels.update('one',regional);assert.equal(pressure.disabled,true);assert.equal(panels.inputs.skin_compression_pa,800);
  assert.match(motor.textContent,/Clear pending whole-skin pressure/);
  pressure.value='300';pressure.onchange();assert.equal(panels.inputs.skin_compression_pa,800);
  const reset=motor.querySelectorAll('button').find(n=>n.textContent==='Clear pending whole-skin pressure');reset.onclick();
  assert.equal(panels.inputs.skin_compression_pa,0);assert.deepEqual(panels.inputs.descending,{bra_r:.4});assert.equal(reset.disabled,true);
  motor.querySelectorAll('button').find(n=>n.textContent==='Release all drives, blocks & pressure').onclick();assert.deepEqual(panels.inputs.descending,{});assert.equal(panels.inputs.skin_compression_pa,0);
  panels.clear();assert.ok(motor.querySelectorAll('input').every(n=>n.disabled));assert.match(hosts['live-body-monitor'].textContent,/No live body owner/);
 }finally{globalThis.document=previous;}
});
