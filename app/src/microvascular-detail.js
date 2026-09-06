const REGIONS=[['body-bp3d-FJ1442','Right vastus lateralis'],['body-bp3d-FJ1442M','Left vastus lateralis']];
const SCENARIO='human_quadriceps_1988_corrected_v1';
const finite3=v=>Array.isArray(v)&&v.length===3&&v.every(Number.isFinite);
const number=v=>Number.isFinite(v)?Number(v).toPrecision(5):'Unavailable';
export function projectMicrovascularPatch(response,{plane='01',zoom=1}={}){
 const patch=response?.patch,edges=patch?.zoom_edges,rotation=patch?.registration?.local_to_body_rotation,origin=patch?.registration?.translation_m;
 if(!Array.isArray(edges)||edges.length>256||!finite3(origin)||!Array.isArray(rotation)||rotation.length!==3||!rotation.every(finite3)||!['01','02','12'].includes(plane)||![1,2,4].includes(zoom))throw Error('Invalid bounded microvascular geometry.');
 if(!Number.isFinite(patch.resolution_m)||patch.resolution_m<=0)throw Error('Missing physical sampling resolution.');
 const axes=[...plane].map(Number);
 const projected=edges.map(edge=>{
  if(!Array.isArray(edge.endpoints_m)||edge.endpoints_m.length!==2||!edge.endpoints_m.every(finite3)||!Number.isFinite(edge.radius_m)||edge.radius_m<=0)throw Error('Vessel endpoints or physical radius unavailable.');
  if(edge.endpoint_pressure_pa!==null&&(!Array.isArray(edge.endpoint_pressure_pa)||edge.endpoint_pressure_pa.length!==2||!edge.endpoint_pressure_pa.every(Number.isFinite)))throw Error('Invalid pressure values.');
  if(edge.flow_m3_per_s!==null&&!Number.isFinite(edge.flow_m3_per_s))throw Error('Invalid flow value.');
  return {edge,points:edge.endpoints_m.map(point=>axes.map(axis=>point.reduce((sum,value,k)=>sum+rotation[k][axis]*(value-origin[k]),0)*1e6))};
 });
 if(!projected.length)return {edges:[],umPerUnit:null,axes};
 const coords=projected.flatMap(edge=>edge.points),bounds=[0,1].map(i=>[Math.min(...coords.map(p=>p[i])),Math.max(...coords.map(p=>p[i]))]);
 const radius=Math.max(...projected.map(p=>p.edge.radius_m))*1e6;
 const center=bounds.map(b=>(b[0]+b[1])/2),umPerUnit=Math.max((bounds[0][1]-bounds[0][0]+2*radius)/360,(bounds[1][1]-bounds[1][0]+2*radius)/260,1e-9)/zoom;
 return {axes,umPerUnit,edges:projected.map(({edge,points})=>({edge,points:points.map(p=>[200+(p[0]-center[0])/umPerUnit,150-(p[1]-center[1])/umPerUnit]),diameterUnits:2*edge.radius_m*1e6/umPerUnit}))};
}
async function fetchPatch(payload){
 const response=await fetch('/api/body/microvascular-patch',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)});
 const data=await response.json();if(!response.ok)throw Error(data.error||'Microvascular patch unavailable.');return data;
}
export function mountMicrovascularDetail(host,{request=fetchPatch,document:doc=globalThis.document}={}){
 const el=(tag,text)=>{const n=doc.createElement(tag);if(text)n.textContent=text;return n;};
 const svgEl=tag=>doc.createElementNS('http://www.w3.org/2000/svg',tag);
 let response=null,generation=0,disposed=false,busy=false;
 const selection=el('p','Default local region: right vastus lateralis.'),region=el('select');region.setAttribute('aria-label','Microvascular region');
 for(const [id,name] of REGIONS){const option=el('option',name);option.value=id;region.append(option);}region.value=REGIONS[0][0];
 const load=el('button','Load local vessel graph');load.type='button';
 const plane=el('select');plane.setAttribute('aria-label','Local vessel projection');
 for(const [value,label] of [['01','Longitudinal · axes 1–2'],['02','Longitudinal · axes 1–3'],['12','Cross section · axes 2–3']]){const o=el('option',label);o.value=value;plane.append(o);}plane.value='01';
 const zoom=el('select');zoom.setAttribute('aria-label','Local vessel zoom');for(const value of [1,2,4]){const o=el('option',`${value}× local view`);o.value=String(value);zoom.append(o);}zoom.value='1';
 const status=el('p','Load the region to inspect its returned hydraulic graph.'),plot=el('div'),scale=el('p'),readout=el('p','Select a vessel for physical quantities.'),scope=el('p'),evidence=el('details'),source=el('pre');status.setAttribute('role','status');
 source.style.whiteSpace='pre-wrap';source.style.overflowWrap='anywhere';evidence.append(el('summary','Scenario, source conditioning & registration'),source);
 host.replaceChildren(selection,region,load,plane,zoom,status,plot,scale,readout,scope,evidence);
 function invalidate(){generation++;response=null;plot.replaceChildren();scale.textContent='';scope.textContent='';source.textContent='';readout.textContent='Select a vessel for physical quantities.';busy=false;load.disabled=disposed;}
 region.onchange=()=>{invalidate();selection.textContent=`Local region: ${REGIONS.find(r=>r[0]===region.value)?.[1]||'unavailable'}.`;status.textContent='Load this region to inspect its graph.';};
 function draw(){
  if(disposed||!response)return;
  const projection=projectMicrovascularPatch(response,{plane:plane.value,zoom:Number(zoom.value)});plot.replaceChildren();
  if(!projection.edges.length){status.textContent='No vessels intersect this local view.';return;}
  const svg=svgEl('svg');svg.setAttribute('viewBox','0 0 400 300');svg.setAttribute('role','img');svg.setAttribute('aria-label','Local microvascular graph with physical vessel diameters');svg.style.width='100%';svg.style.height='300px';svg.style.background='#142225';
  const pressures=projection.edges.flatMap(({edge})=>edge.endpoint_pressure_pa||[]),low=Math.min(...pressures),high=Math.max(...pressures);
  for(const item of projection.edges){
   const {edge,points}=item,line=svgEl('line');for(const [key,value] of Object.entries({x1:points[0][0],y1:points[0][1],x2:points[1][0],y2:points[1][1],'stroke-width':item.diameterUnits}))line.setAttribute(key,String(value));
   const mean=edge.endpoint_pressure_pa?.reduce((a,b)=>a+b,0)/2,fraction=Number.isFinite(mean)?(mean-low)/(high-low||1):0;
   line.setAttribute('stroke',Number.isFinite(mean)?`hsl(${220-220*fraction} 65% 65%)`:'#9cb5b2');line.setAttribute('tabindex','0');line.setAttribute('role','button');line.setAttribute('aria-label',`Inspect vessel ${edge.edge_index}`);
   const description=`Vessel ${edge.edge_index} · diameter ${number(edge.radius_m*2e6)} µm · endpoint pressure ${edge.endpoint_pressure_pa?.map(number).join(' → ')||'Unavailable'} Pa · flow ${number(edge.flow_m3_per_s)} m³/s.`;
   const title=svgEl('title');title.textContent=description;line.append(title);line.onclick=()=>{readout.textContent=description;};line.onkeydown=e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();readout.textContent=description;}};svg.append(line);
  }
  const barUm=10**Math.floor(Math.log10(projection.umPerUnit*100)),bar=svgEl('path'),barLabel=svgEl('text');
  bar.setAttribute('d',`M20 278H${20+barUm/projection.umPerUnit}`);bar.setAttribute('stroke','#edf4ee');bar.setAttribute('stroke-width','2');
  barLabel.setAttribute('x','20');barLabel.setAttribute('y','270');barLabel.setAttribute('fill','#edf4ee');barLabel.setAttribute('font-size','11');barLabel.textContent=`${number(barUm)} µm`;svg.append(bar,barLabel);
  plot.append(svg);
  scale.textContent=`Local axes ${projection.axes.map(i=>i+1).join('–')} · ${number(projection.umPerUnit)} µm per SVG view unit · centerline sampling ${number(response.patch.resolution_m*1e6)} µm. Diameter uses physical scale; sampling is numerical and biological resolution remains source-limited. Color: mean endpoint pressure ${Number.isFinite(low)?number(low)+'–'+number(high)+' Pa':'unavailable'}.`;
  status.textContent=`${response.selection?.name||REGIONS.find(r=>r[0]===region.value)?.[1]} · ${projection.edges.length} returned vessels · ${response.scenario?.label||response.scenario?.id||SCENARIO}.`;
 }
 plane.onchange=draw;zoom.onchange=draw;
 load.onclick=async()=>{
  if(disposed||busy)return;invalidate();const token=generation;busy=true;load.disabled=true;status.textContent='Loading the bounded local hydraulic graph…';
  try{
   const next=await request({entity_id:region.value,resolution_m:50e-6,scenario_id:SCENARIO});if(disposed||token!==generation)return;
   projectMicrovascularPatch(next);response=next;
   scope.textContent=[next.patch.registration?.evidence,next.patch.ownership,next.patch.boundary_scope].filter(Boolean).join(' ');
   source.textContent=JSON.stringify({scenario:next.scenario,selection:next.selection,conditioning:next.patch.graph?.conditioning,radius_conditioning:next.patch.graph?.radius_conditioning,hydraulic_model:next.patch.graph?.hydraulic_model,registration:next.patch.registration,containment:next.patch.containment},null,2);
   draw();
  }catch(error){if(!disposed&&token===generation){response=null;plot.replaceChildren();status.textContent=`Local graph unavailable: ${error.message}`;}}
  finally{if(token===generation){busy=false;load.disabled=disposed;}}
 };
 return {updateSelection(structure){
  if(disposed)return;invalidate();const found=REGIONS.find(r=>r[0]===structure?.id);
  if(found){region.value=found[0];selection.textContent=`Selected local region: ${found[1]}.`;}
  else selection.textContent=`${structure?.name?structure.name+': ':''}This patch prior supports vastus lateralis. Region selector explicitly chooses the local muscle patch.`;
  status.textContent='Load the named region to inspect its graph.';
 },clear(){if(disposed)return;invalidate();selection.textContent='Choose a local vastus lateralis region.';status.textContent='No local graph loaded.';},dispose(){disposed=true;invalidate();region.disabled=true;plane.disabled=true;zoom.disabled=true;}};
}
