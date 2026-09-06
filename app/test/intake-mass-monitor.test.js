import test from 'node:test';
import assert from 'node:assert/strict';
import {intakeMassSnapshot,mountIntakeMassMonitor} from '../src/intake-mass-monitor.js';
class Element{
 constructor(tag){this.tagName=tag;this.children=[];this.attributes={};this.style={};this._text='';}
 set textContent(text){this._text=String(text);this.children=[];}get textContent(){return this._text+this.children.map(n=>n.textContent).join('');}
 append(...nodes){this.children.push(...nodes);}replaceChildren(...nodes){this._text='';this.children=nodes;}setAttribute(k,v){this.attributes[k]=v;}
}
const frame=()=>({schema:'ihm.embodied-frame.v1',time_s:.2,sequence:12,mechanics:{effective_native_body_mass_kg:82.1},intake_mass:{enabled:true,binding:{body:'torso',canonical_entity_id:'stomach',station_source_m:[.1,.2,.3],incoming_velocity_basis:'co_moving_at_ingestion_assumption',scope:'Inferred rigid torso; no internal deformation',registration_sha256:'retained'},bridge:{applied_mass_kg:.1,capacity_kg:.5,failed:false,pending_boundary:null,last_boundary:{mass_kg:.1,interval_start_tick:9,interval_end_tick:10,mechanical_transfer_applied:true,boundary_accounted:true},mechanical_sequence:1,mechanical_reference_id:'native'}},intake_schedule:{events:[{meal:{water_ml:100000}}]}});
const setup=()=>{const host=new Element('div');return {host,view:mountIntakeMassMonitor(host,{document:{createElement:t=>new Element(t)}})};};
test('mass snapshot reports accepted mechanical values and never calculates scheduled mass',()=>{
 const data=frame(),snapshot=intakeMassSnapshot(data);assert.equal(snapshot.mass_kg,82.1);assert.equal(snapshot.applied_mass_kg,.1);assert.equal(snapshot.capacity_kg,.5);
 data.intake_mass={enabled:false,bridge:{applied_mass_kg:999}};assert.equal(intakeMassSnapshot(data).applied_mass_kg,null);
 assert.equal(intakeMassSnapshot({schema:'replay'}),null);
});
test('disabled or missing coupling fields are unavailable while actual body mass remains shown',()=>{
 const {host,view}=setup(),data=frame();data.intake_mass={enabled:false};view.update('body',data);
 assert.match(host.textContent,/coupling disabled/);assert.match(host.textContent,/82.1000 kg/);assert.match(host.textContent,/Applied intake massUnavailable/);assert.doesNotMatch(host.textContent,/100000|0.500000 kg/);
 delete data.intake_mass;view.update('body',data);assert.match(host.textContent,/coupling unavailable/);
});
test('pending boundaries remain unconfirmed and explicit zero values remain actual zeros',()=>{
 const {host,view}=setup(),data=frame();data.intake_mass.bridge={applied_mass_kg:0,capacity_kg:.5,failed:true,pending_boundary:{mass_kg:.1,interval_start_tick:9,interval_end_tick:10},last_boundary:null};view.update('body',data);
 assert.match(host.textContent,/transfer failed/);assert.match(host.textContent,/Applied intake mass0.00000 kg/);assert.match(host.textContent,/Pending consumed boundary0.100000 kg/);assert.match(host.textContent,/mechanical transfer unconfirmed/);assert.doesNotMatch(host.textContent,/transfer confirmed/);
});
test('owner/reset/disposal clears receipts and compact details preserve actual source assumptions',()=>{
 const {host,view}=setup();view.update('body',frame());assert.match(host.textContent,/Incoming velocity basis: co moving at ingestion assumption/);assert.match(host.textContent,/0.10000, 0.20000, 0.30000 m/);assert.match(host.textContent,/mechanical transfer confirmed/);
 const next=frame();next.intake_mass={enabled:false};view.update('new',next);assert.doesNotMatch(host.textContent,/registration: retained/);
 view.clear();assert.match(host.textContent,/No live body mass frame/);view.dispose();view.update('body',frame());assert.doesNotMatch(host.textContent,/82.1000/);
});
