import test from 'node:test';
import assert from 'node:assert/strict';
import {SKIN_VOLTAGE_KEY,skinVoltageMv,mountSkinVoltageMonitor} from '../src/skin-voltage-monitor.js';
class Element{
 constructor(tag){this.tagName=tag;this.children=[];this.attributes={};this._text='';}
 set textContent(text){this._text=String(text);this.children=[];}get textContent(){return this._text+this.children.map(n=>n.textContent).join('');}
 append(...nodes){this.children.push(...nodes);}replaceChildren(...nodes){this._text='';this.children=nodes;}setAttribute(k,v){this.attributes[k]=v;}
 querySelectorAll(tag){return this.children.flatMap(n=>[...(n.tagName===tag?[n]:[]),...n.querySelectorAll(tag)]);}
}
const frame=(time,sequence,volts)=>({schema:'ihm.embodied-frame.v1',time_s:time,sequence,physiology:{values:{[SKIN_VOLTAGE_KEY]:volts}}});
const setup=()=>{const host=new Element('div');return {host,view:mountSkinVoltageMonitor(host,{document:{createElement:t=>new Element(t),createElementNS:(_,t)=>new Element(t)}})};};
test('bulk Skin voltage uses only the explicit native volts field and converts to mV',()=>{
 assert.equal(skinVoltageMv(frame(0,0,-.083)),-83);assert.equal(skinVoltageMv(frame(0,0,0)),0);
 for(const value of [null,undefined,NaN,Infinity,'-.083',false])assert.equal(skinVoltageMv(frame(0,0,value)),null);
 assert.equal(skinVoltageMv({schema:'ihm.embodied-frame.v1',physiology:{values:{epidermal_tep_mv:-83}}}),null);
});
test('missing native values show unavailable with no fabricated graph',()=>{
 const {host,view}=setup();view.update('body',frame(0,0,null));view.update('body',frame(.2,1,null));
 assert.match(host.textContent,/Unavailable/);assert.equal(host.querySelectorAll('svg').length,0);
 assert.match(host.textContent,/not epidermal TEP or a measured keratinocyte voltage/);
});
test('plot uses actual timestamps, suppresses duplicate-time schedule samples and preserves gaps',()=>{
 const {host,view}=setup();view.update('body',frame(0,0,-.083));view.update('body',frame(.2,1,-.082));view.update('body',frame(.2,2,-.082));view.update('body',frame(.4,3,-.081));
 assert.match(host.textContent,/3 finite samples/);assert.match(host.textContent,/0.00–0.40 s/);assert.equal(host.querySelectorAll('path')[0].attributes.d.split('L').length-1,2);
 view.update('body',frame(.6,4,null));assert.equal(host.querySelectorAll('svg').length,0);
 view.update('body',frame(.8,5,-.08));assert.equal(host.querySelectorAll('path')[0].attributes.d.split('M').length-1,2);
});
test('owner switch, reset and disposal clear old voltage samples',()=>{
 const {host,view}=setup();view.update('a',frame(0,0,-.083));view.update('a',frame(.2,1,-.082));
 view.update('b',frame(0,0,-.081));assert.equal(host.querySelectorAll('svg').length,0);assert.match(host.textContent,/-81.0000 mV/);
 view.clear();assert.match(host.textContent,/No native Skin tissue voltage frame/);view.dispose();view.update('b',frame(.2,1,-.08));assert.doesNotMatch(host.textContent,/-80.0000/);
});
