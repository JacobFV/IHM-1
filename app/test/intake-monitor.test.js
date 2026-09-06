import test from 'node:test';
import assert from 'node:assert/strict';
import {mountIntakeMonitor,validateIntakeEvent} from '../src/intake-monitor.js';
class Element {
 constructor(tag){this.tagName=tag;this.children=[];this.attributes={};this.value='';this._text='';}
 set textContent(t){this._text=String(t);this.children=[];} get textContent(){return this._text+this.children.map(c=>c.textContent).join('');}
 append(...nodes){this.children.push(...nodes);} replaceChildren(...nodes){this._text='';this.children=[...nodes];}
 setAttribute(k,v){this.attributes[k]=v;}
 querySelectorAll(tag){return this.children.flatMap(c=>[...(c.tagName===tag?[c]:[]),...c.querySelectorAll(tag)]);}
}
const event={event_id:'water',time_s:.04,meal:{name:'drink',carbohydrate_g:0,protein_g:0,fat_g:0,sodium_g:0,calcium_mg:0,water_ml:250}};
function setup(submit){const host=new Element('div');const monitor=mountIntakeMonitor(host,{submit,document:{createElement:t=>new Element(t)}});return {host,monitor};}
function fill(host){for(const n of host.querySelectorAll('input'))n.value=String(n.name==='event_id'?event.event_id:n.name==='time_s'?event.time_s:event.meal[n.name]);}
test('native field validation rejects nonfinite, fractional ticks, empty food and coercion',()=>{
 assert.deepEqual(validateIntakeEvent(event),event);
 for(const value of [NaN,Infinity,-1,86400.02,.021,true,'0.04'])assert.throws(()=>validateIntakeEvent({...event,time_s:value}));
 for(const value of [NaN,Infinity,-1,10001,true,'250',null])assert.throws(()=>validateIntakeEvent({...event,meal:{...event.meal,water_ml:value}}));
 assert.throws(()=>validateIntakeEvent({...event,meal:{...event.meal,water_ml:0}}));
 assert.throws(()=>validateIntakeEvent({...event,meal:{...event.meal,energy_kcal:1}}));
});
test('DOM submit calls injected queue callback once and only reports authoritative state',async()=>{
 let resolve;const calls=[];const {host,monitor}=setup(events=>{calls.push(events);return new Promise(r=>resolve=r);});
 monitor.setConnected(true);fill(host);const form=host.querySelectorAll('form')[0];const pending=form.onsubmit({preventDefault(){}});
 await form.onsubmit({preventDefault(){}});assert.equal(calls.length,1);assert.doesNotMatch(host.textContent,/water · issued/);
 resolve({events:[{...event,state:'queued'}]});await pending;assert.match(host.textContent,/water · queued/);
 await form.onsubmit({preventDefault(){}});assert.equal(calls.length,1);assert.match(host.textContent,/already exists/);
 monitor.update({events:[{...event,state:'issued',issued_tick:2}]});assert.match(host.textContent,/water · issued/);
 monitor.update({events:[{...event,state:'accepted',issued_tick:2,native_sequence:4}]});assert.match(host.textContent,/water · accepted/);
 assert.throws(()=>monitor.update({events:[{...event,state:'queued'}]}));
});
test('uncertain callback locks further submission without inventing native acceptance',async()=>{
 let calls=0;const {host,monitor}=setup(async()=>{calls++;throw Error('response lost');});monitor.setConnected(true);fill(host);
 const form=host.querySelectorAll('form')[0];await form.onsubmit({preventDefault(){}});await form.onsubmit({preventDefault(){}});
 assert.equal(calls,1);assert.match(host.textContent,/request outcome uncertain/i);assert.doesNotMatch(host.textContent,/water · accepted/);
 monitor.setConnected(false);monitor.setConnected(true);await form.onsubmit({preventDefault(){}});assert.equal(calls,1);
});
test('authoritative uncertainty remains terminal and disables UI',()=>{
 const {host,monitor}=setup(async()=>({events:[]}));monitor.setConnected(true);monitor.update({events:[{...event,state:'uncertain',reason:'timeout'}]});
 assert.match(host.textContent,/water · uncertain/);assert.ok(host.querySelectorAll('button')[0].disabled);
 assert.throws(()=>monitor.update({events:[{...event,state:'accepted'}]}));
});

test('definite queue rejection permits correction while disposed owner ignores late response',async()=>{
 let calls=0;const {host,monitor}=setup(async()=>{calls++;throw Object.assign(Error('past time'),{definitelyRejected:true});});
 monitor.setConnected(true);fill(host);const form=host.querySelectorAll('form')[0];await form.onsubmit({preventDefault(){}});await form.onsubmit({preventDefault(){}});assert.equal(calls,2);assert.match(host.textContent,/request rejected/);
 let resolve;const old=setup(()=>new Promise(r=>resolve=r));old.monitor.setConnected(true);fill(old.host);const pending=old.host.querySelectorAll('form')[0].onsubmit({preventDefault(){}});old.monitor.dispose();resolve({events:[{...event,state:'accepted'}]});await pending;assert.doesNotMatch(old.host.textContent,/water · accepted/);assert.ok(old.host.querySelectorAll('button')[0].disabled);
});
