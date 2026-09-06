import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {mountMicrovascularDetail,projectMicrovascularPatch} from '../src/microvascular-detail.js';
class Element{
 constructor(tag){this.tagName=tag;this.children=[];this.attributes={};this.style={};this.value='';this._text='';}
 set textContent(text){this._text=String(text);this.children=[];}get textContent(){return this._text+this.children.map(n=>n.textContent).join('');}
 append(...nodes){this.children.push(...nodes);}replaceChildren(...nodes){this._text='';this.children=nodes;}
 setAttribute(k,v){this.attributes[k]=v;}querySelectorAll(tag){return this.children.flatMap(n=>[...(n.tagName===tag?[n]:[]),...n.querySelectorAll(tag)]);}
}
const fixture=()=>({scenario:{id:'declared',label:'Declared source scenario'},selection:{name:'Right vastus lateralis'},patch:{resolution_m:50e-6,registration:{translation_m:[1,2,3],local_to_body_rotation:[[1,0,0],[0,1,0],[0,0,1]],evidence:'Synthetic geometry; inferred fiber direction'},ownership:'No independent native blood store',boundary_scope:'Macro-vessel correspondence unresolved',graph:{conditioning:{source:'source A'}},zoom_edges:[{edge_id:'one',edge_index:0,endpoints_m:[[1,2,3],[1.0005,2.0002,3]],centerline_samples_m:[[1,2,3],[1.0005,2.0002,3]],radius_m:3e-6,endpoint_pressure_pa:[4000,2000],flow_m3_per_s:1e-12}]}});
const setup=request=>{const host=new Element('div');return {host,view:mountMicrovascularDetail(host,{request,document:{createElement:t=>new Element(t),createElementNS:(_,t)=>new Element(t)}})};};
test('projection uses canonical rotation and real vessel diameter, zoom changes only view scale',()=>{
 const data=fixture(),snapshot=JSON.stringify(data),fit=projectMicrovascularPatch(data),zoom=projectMicrovascularPatch(data,{zoom:2});
 assert.ok(Math.abs(fit.edges[0].diameterUnits*fit.umPerUnit-6)<1e-12);
 assert.equal(zoom.edges[0].diameterUnits,fit.edges[0].diameterUnits*2);assert.equal(JSON.stringify(data),snapshot);
 data.patch.registration.local_to_body_rotation=[[0,-1,0],[1,0,0],[0,0,1]];
 const rotated=projectMicrovascularPatch(data);assert.ok(rotated.edges[0].points[1][0]>rotated.edges[0].points[0][0]);assert.ok(rotated.edges[0].points[1][1]>rotated.edges[0].points[0][1]);
});
test('invalid radii and nonfinite physical quantities cannot render decorative substitutes',()=>{
 for(const radius of [null,0,-1,NaN]){const data=fixture();data.patch.zoom_edges[0].radius_m=radius;assert.throws(()=>projectMicrovascularPatch(data));}
 const data=fixture();data.patch.zoom_edges[0].endpoint_pressure_pa=[Infinity,1];assert.throws(()=>projectMicrovascularPatch(data));
});
test('monitor requests explicit region/scenario only on load, exposes authoritative units and sources',async()=>{
 const calls=[];const {host,view}=setup(async payload=>{calls.push(payload);return fixture();});assert.equal(calls.length,0);
 view.updateSelection({id:'body-bp3d-FJ1442M',name:'Left vastus lateralis'});
 await host.querySelectorAll('button')[0].onclick();assert.deepEqual(calls,[{entity_id:'body-bp3d-FJ1442M',resolution_m:50e-6,scenario_id:'human_quadriceps_1988_corrected_v1'}]);
 const line=host.querySelectorAll('polyline')[0];line.onclick();assert.match(host.textContent,/diameter 6\.0000 µm/);assert.match(host.textContent,/4000\.0 → 2000\.0 Pa/);assert.match(host.textContent,/1\.0000e-12 m³\/s/);assert.match(host.textContent,/source A/);assert.match(host.textContent,/Declared source scenario/);
 host.querySelectorAll('select')[3].value='2';host.querySelectorAll('select')[3].onchange();assert.equal(calls.length,1);
});
test('selection changes and disposal discard old requests without repainting another region',async()=>{
 let resolve;const {host,view}=setup(()=>new Promise(r=>resolve=r));const pending=host.querySelectorAll('button')[0].onclick();view.updateSelection({id:'body-bp3d-FJ1442M'});resolve(fixture());await pending;assert.equal(host.querySelectorAll('polyline').length,0);
 const again=host.querySelectorAll('button')[0].onclick();view.dispose();resolve(fixture());await again;assert.equal(host.querySelectorAll('polyline').length,0);assert.equal(host.querySelectorAll('button')[0].disabled,true);
});

test('tortuous vessel geometry and slider pressure follow retained arc length, not endpoint chord',async()=>{
 const data=fixture(),edge=data.patch.zoom_edges[0];edge.endpoints_m=[[1,2,3],[1.0001,2,3]];
 edge.centerline_samples_m=[[1,2,3],[1,2.0001,3],[1.0001,2.0001,3],[1.0001,2,3]];
 const projected=projectMicrovascularPatch(data).edges[0];assert.equal(projected.points.length,4);assert.ok(Math.abs(projected.length_m-300e-6)<1e-12);
 assert.ok(Math.abs(projected.arcFractions[1]-1/3)<1e-10);
 const {host}=setup(async()=>data);await host.querySelectorAll('button')[0].onclick();
 const line=host.querySelectorAll('polyline')[0];assert.equal(line.attributes.points.split(' ').length,4);line.onclick();
 assert.match(host.textContent,/At 150.00 µm along this fragment: 3000.0 Pa/);
 const position=host.querySelectorAll('input')[0];position.value='25';position.oninput();assert.match(host.textContent,/At 75.000 µm along this fragment: 3500.0 Pa/);
 const marker=host.querySelectorAll('circle')[0];assert.ok(Math.abs(Number(marker.attributes.cx)-projected.points[0][0])<1e-9);
});
test('kidney selection requests explicit control or injury scenario without muscle parameters',async()=>{
 const data=fixture();data.patch.containment={whole_patch_inside_authored_surface:true,cortical_location_verified:false};data.patch.boundary_scope='No native terminal registration';
 const calls=[],{host,view}=setup(async payload=>{calls.push(payload);return data;});view.updateSelection({id:'body-bp3d-FJ3147',name:'Right kidney'});
 await host.querySelectorAll('button')[0].onclick();assert.deepEqual(calls[0],{entity_id:'body-bp3d-FJ3147',resolution_m:50e-6,scenario_id:'human_kidney_control_v1'});
 assert.match(host.textContent,/containment: certified; cortical location remains unverified/);assert.match(host.textContent,/No native terminal registration/);
 const scenario=host.querySelectorAll('select')[1];scenario.value='human_kidney_injury_v1';scenario.onchange();await host.querySelectorAll('button')[0].onclick();assert.equal(calls[1].scenario_id,'human_kidney_injury_v1');
 assert.equal(host.querySelectorAll('polyline').length,1);
});
test('missing polyline geometry and samples inconsistent with endpoints are rejected',()=>{
 const data=fixture();delete data.patch.zoom_edges[0].centerline_samples_m;assert.throws(()=>projectMicrovascularPatch(data),/polyline samples/);
 data.patch.zoom_edges[0].centerline_samples_m=[[1,2,3],[1.2,2,3]];assert.throws(()=>projectMicrovascularPatch(data),/endpoints/);
});
test('skin selection requires named territory and exposes transferred priors and uncertainty',async()=>{
 const data=fixture();data.selection.name='Skin · left lower-leg anteromedial';
 data.patch.registration.territory_id='left_lower_leg_anteromedial';data.patch.registration.source_triangle_id=126303;
 data.patch.uncertainty={posterior_inference_performed:false,summary:'Forearm-to-lower-leg transfer; posterior uncertainty has not been inferred.'};
 data.patch.graph.source_conditioning={human_confocal:{site:'middle volar forearm'}};
 const calls=[],{host,view}=setup(async payload=>{calls.push(payload);return data;});view.updateSelection({id:'body-bp3d-FJ2810',name:'Skin'});
 await host.querySelectorAll('button')[0].onclick();
 assert.deepEqual(calls[0],{entity_id:'body-bp3d-FJ2810',territory_id:'left_lower_leg_anteromedial',patch_area_mm2:1,resolution_m:50e-6,scenario_id:'forearm_baseline_transfer_to_lower_leg_v1'});
 assert.match(host.textContent,/posterior uncertainty has not been inferred/);assert.match(host.textContent,/126303/);assert.match(host.textContent,/middle volar forearm/);
 const territory=host.querySelectorAll('select').find(e=>e.attributes['aria-label']==='Exact engineered skin territory');territory.value='right_lower_leg_posterolateral';territory.onchange();
 await host.querySelectorAll('button')[0].onclick();assert.equal(calls[1].territory_id,'right_lower_leg_posterolateral');
 assert.equal(host.querySelectorAll('polyline').length,1);
});

test('actual retained skin API fixture renders all physical loop vertices in tangent-depth view',async()=>{
 const data=JSON.parse(readFileSync(new URL('../../data/research/skin_microvascular/api_patch.json',import.meta.url),'utf8'));
 const projection=projectMicrovascularPatch(data,{plane:'02'}),caps=projection.edges.filter(x=>x.edge.vessel_class==='papillary_capillary_loop');
 assert.equal(caps.length,7);assert.equal(projection.edges.length,25);
 for(const cap of caps){assert.ok(cap.points.length>=33);assert.ok(Math.abs(cap.diameterUnits*projection.umPerUnit-9.59)<1e-10);assert.ok(Math.max(...cap.points.map(p=>p[1]))-Math.min(...cap.points.map(p=>p[1]))>1);}
 const {host,view}=setup(async()=>data);view.updateSelection({id:'body-bp3d-FJ2810'});
 await host.querySelectorAll('button')[0].onclick();assert.equal(host.querySelectorAll('polyline').length,25);
 assert.match(host.textContent,/zero additional native blood volume/);assert.match(host.textContent,/posterior uncertainty has not been inferred/);
 host.querySelectorAll('polyline')[5].onclick();assert.match(host.textContent,/papillary capillary loop/);
});
