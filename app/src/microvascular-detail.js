const REGIONS=[['body-bp3d-FJ1442','Right vastus lateralis'],['body-bp3d-FJ1442M','Left vastus lateralis'],['body-bp3d-FJ3145','Left kidney'],['body-bp3d-FJ3147','Right kidney'],['body-bp3d-FJ2810','Skin · exact lower-leg territory']];
const skin=id=>id==='body-bp3d-FJ2810';
const SKIN_SCENARIO='forearm_baseline_transfer_to_lower_leg_v1';
const SKIN_TERRITORIES=['left','right'].flatMap(side=>['anteromedial','anterolateral','posteromedial','posterolateral'].map(part=>`${side}_lower_leg_${part}`));
const SCENARIO='human_quadriceps_1988_corrected_v1';
const kidney=id=>['body-bp3d-FJ3145','body-bp3d-FJ3147'].includes(id);
const finite3=v=>Array.isArray(v)&&v.length===3&&v.every(Number.isFinite);
const number=v=>Number.isFinite(v)?Number(v).toPrecision(5):'Unavailable';
export function projectMicrovascularPatch(response,{plane='01',zoom=1}={}){
 const patch=response?.patch,edges=patch?.zoom_edges,rotation=patch?.registration?.local_to_body_rotation,origin=patch?.registration?.translation_m;
 if(!Array.isArray(edges)||edges.length>256||!finite3(origin)||!Array.isArray(rotation)||rotation.length!==3||!rotation.every(finite3)||!['01','02','12'].includes(plane)||![1,2,4].includes(zoom))throw Error('Invalid bounded microvascular geometry.');
 if(!Number.isFinite(patch.resolution_m)||patch.resolution_m<=0)throw Error('Missing physical sampling resolution.');
 const axes=[...plane].map(Number);
 let sampleCount=0;
 const projected=edges.map(edge=>{
  if(!Array.isArray(edge.endpoints_m)||edge.endpoints_m.length!==2||!edge.endpoints_m.every(finite3)||!Number.isFinite(edge.radius_m)||edge.radius_m<=0)throw Error('Vessel endpoints or physical radius unavailable.');
  if(edge.endpoint_pressure_pa!==null&&(!Array.isArray(edge.endpoint_pressure_pa)||edge.endpoint_pressure_pa.length!==2||!edge.endpoint_pressure_pa.every(Number.isFinite)))throw Error('Invalid pressure values.');
  if(edge.flow_m3_per_s!==null&&!Number.isFinite(edge.flow_m3_per_s))throw Error('Invalid flow value.');
  const samples=edge.centerline_samples_m;
  if(!Array.isArray(samples)||samples.length<2||!samples.every(finite3)||(sampleCount+=samples.length)>20000)throw Error('Retained vessel polyline samples unavailable or over budget.');
  if(samples[0].some((v,i)=>Math.abs(v-edge.endpoints_m[0][i])>1e-9)||samples.at(-1).some((v,i)=>Math.abs(v-edge.endpoints_m[1][i])>1e-9))throw Error('Polyline endpoints do not match the returned fragment.');
  const arc=[0];for(let i=1;i<samples.length;i++)arc.push(arc.at(-1)+Math.hypot(...samples[i].map((v,k)=>v-samples[i-1][k])));
  const length_m=arc.at(-1);if(!(length_m>0))throw Error('Zero-length vessel fragment.');
  return {edge,length_m,arcFractions:arc.map(value=>value/length_m),points:samples.map(point=>axes.map(axis=>point.reduce((sum,value,k)=>sum+rotation[k][axis]*(value-origin[k]),0)*1e6))};
 });
 if(!projected.length)return {edges:[],umPerUnit:null,axes};
 const coords=projected.flatMap(edge=>edge.points),bounds=[0,1].map(i=>[Math.min(...coords.map(p=>p[i])),Math.max(...coords.map(p=>p[i]))]);
 const radius=Math.max(...projected.map(p=>p.edge.radius_m))*1e6;
 const center=bounds.map(b=>(b[0]+b[1])/2),umPerUnit=Math.max((bounds[0][1]-bounds[0][0]+2*radius)/360,(bounds[1][1]-bounds[1][0]+2*radius)/260,1e-9)/zoom;
 return {axes,umPerUnit,edges:projected.map(({edge,points,length_m,arcFractions})=>({edge,length_m,arcFractions,points:points.map(p=>[200+(p[0]-center[0])/umPerUnit,150-(p[1]-center[1])/umPerUnit]),diameterUnits:2*edge.radius_m*1e6/umPerUnit}))};
}
async function fetchPatch(payload){
 const response=await fetch('/api/body/microvascular-patch',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)});
 const data=await response.json();if(!response.ok)throw Error(data.error||'Microvascular patch unavailable.');return data;
}
export function mountMicrovascularDetail(host,{request=fetchPatch,document:doc=globalThis.document}={}){
 const el=(tag,text)=>{const n=doc.createElement(tag);if(text)n.textContent=text;return n;};
 const svgEl=tag=>doc.createElementNS('http://www.w3.org/2000/svg',tag);
 let response=null,generation=0,disposed=false,busy=false,selectedEdge=null,marker=null;
 const selection=el('p','Default local region: right vastus lateralis.'),region=el('select');region.setAttribute('aria-label','Microvascular region');
 for(const [id,name] of REGIONS){const option=el('option',name);option.value=id;region.append(option);}region.value=REGIONS[0][0];
 const scenario=el('select');scenario.setAttribute('aria-label','Microvascular scenario');
 function scenarios(){scenario.replaceChildren();for(const [id,name] of skin(region.value)?[[SKIN_SCENARIO,'Skin · forearm baseline transferred to lower leg']]:kidney(region.value)?[['human_kidney_control_v1','Kidney · control-conditioned scenario'],['human_kidney_injury_v1','Kidney · injury-conditioned scenario']]:[[SCENARIO,'Muscle · corrected quadriceps scenario']]){const option=el('option',name);option.value=id;scenario.append(option);}scenario.value=skin(region.value)?SKIN_SCENARIO:kidney(region.value)?'human_kidney_control_v1':SCENARIO;}
 scenarios();
 const territory=el('select');territory.setAttribute('aria-label','Exact engineered skin territory');
 for(const id of SKIN_TERRITORIES){const option=el('option',id.replaceAll('_',' '));option.value=id;territory.append(option);}territory.value=SKIN_TERRITORIES[0];territory.hidden=true;
 const area=el('input');area.type='number';area.min='0.15';area.max='8';area.step='0.05';area.value='1';area.setAttribute('aria-label','Skin patch area in square millimetres');area.hidden=true;
 const skinScope=el('p','Skin patch area (mm²). Named territories are engineered source-skin masks; they are not measured vascular watersheds.');skinScope.hidden=true;
 function skinControls(){territory.hidden=area.hidden=skinScope.hidden=!skin(region.value);}
 const load=el('button','Load local vessel graph');load.type='button';
 const plane=el('select');plane.setAttribute('aria-label','Local vessel projection');
 for(const [value,label] of [['01','Longitudinal · axes 1–2'],['02','Longitudinal · axes 1–3'],['12','Cross section · axes 2–3']]){const o=el('option',label);o.value=value;plane.append(o);}plane.value='01';
 const zoom=el('select');zoom.setAttribute('aria-label','Local vessel zoom');for(const value of [1,2,4]){const o=el('option',`${value}× local view`);o.value=String(value);zoom.append(o);}zoom.value='1';
 const position=el('input');position.type='range';position.min='0';position.max='100';position.step='1';position.value='50';position.disabled=true;position.setAttribute('aria-label','Position along selected vessel · percent of arc length');
 const status=el('p','Load the region to inspect its returned hydraulic graph.'),plot=el('div'),scale=el('p'),readout=el('p','Select a vessel for physical quantities.'),scope=el('p'),evidence=el('details'),source=el('pre');status.setAttribute('role','status');
 source.style.whiteSpace='pre-wrap';source.style.overflowWrap='anywhere';evidence.append(el('summary','Scenario, source conditioning & registration'),source);
 host.replaceChildren(selection,region,scenario,load,plane,zoom,territory,status,plot,scale,el('label','Position along selected vessel · percent of arc length'),position,skinScope,area,readout,scope,evidence);
 function invalidate(){generation++;response=null;selectedEdge=null;marker=null;position.disabled=true;plot.replaceChildren();scale.textContent='';scope.textContent='';source.textContent='';readout.textContent='Select a vessel for physical quantities.';busy=false;load.disabled=disposed;}
 territory.onchange=area.onchange=()=>{invalidate();status.textContent='Load this exact skin territory and bounded footprint.';};
 region.onchange=()=>{invalidate();scenarios();skinControls();plane.value=skin(region.value)?'02':'01';selection.textContent=`Local region: ${REGIONS.find(r=>r[0]===region.value)?.[1]||'unavailable'}.`;status.textContent='Load this region to inspect its graph.';};
 scenario.onchange=()=>{invalidate();status.textContent='Load the selected declared scenario.';};
 function inspect(item){
  selectedEdge=item;position.disabled=false;const t=Number(position.value)/100,{edge,points,arcFractions,length_m}=item;
  let i=1;while(i<arcFractions.length-1&&arcFractions[i]<t)i++;
  const local=(t-arcFractions[i-1])/(arcFractions[i]-arcFractions[i-1]||1),point=points[i].map((value,k)=>points[i-1][k]+local*(value-points[i-1][k]));
  marker?.setAttribute('cx',String(point[0]));marker?.setAttribute('cy',String(point[1]));marker?.setAttribute('visibility','visible');
  const pressure=edge.endpoint_pressure_pa?edge.endpoint_pressure_pa[0]+t*(edge.endpoint_pressure_pa[1]-edge.endpoint_pressure_pa[0]):null;
  readout.textContent=`Vessel ${edge.edge_index}${edge.vessel_class?' · '+edge.vessel_class.replaceAll('_',' '):''}${edge.fragment_index==null?'':` fragment ${edge.fragment_index}`} · diameter ${number(edge.radius_m*2e6)} µm · retained arc ${number(length_m*1e6)} µm · endpoint pressure ${edge.endpoint_pressure_pa?.map(number).join(' → ')||'Unavailable'} Pa · flow ${number(edge.flow_m3_per_s)} m³/s. At ${number(t*length_m*1e6)} µm along this fragment: ${number(pressure)} Pa (interpolated by arc length).`;
 }
 position.oninput=()=>{if(selectedEdge)inspect(selectedEdge);};
 function draw(){
  if(disposed||!response)return;
  const projection=projectMicrovascularPatch(response,{plane:plane.value,zoom:Number(zoom.value)});plot.replaceChildren();
  if(!projection.edges.length){status.textContent='No vessels intersect this local view.';return;}
  const svg=svgEl('svg');svg.setAttribute('viewBox','0 0 400 300');svg.setAttribute('role','img');svg.setAttribute('aria-label','Local microvascular graph with physical vessel diameters');svg.style.width='100%';svg.style.height='300px';svg.style.background='#0a0e11';
  const pressures=projection.edges.flatMap(({edge})=>edge.endpoint_pressure_pa||[]),low=Math.min(...pressures),high=Math.max(...pressures);
  for(const item of projection.edges){
   const {edge,points}=item,line=svgEl('polyline');line.setAttribute('points',points.map(point=>point.join(',')).join(' '));line.setAttribute('fill','none');line.setAttribute('stroke-width',String(item.diameterUnits));
   const mean=edge.endpoint_pressure_pa?.reduce((a,b)=>a+b,0)/2,fraction=Number.isFinite(mean)?(mean-low)/(high-low||1):0;
   line.setAttribute('stroke',Number.isFinite(mean)?`hsl(${220-220*fraction} 65% 65%)`:'#9cb5b2');line.setAttribute('tabindex','0');line.setAttribute('role','button');line.setAttribute('aria-label',`Inspect vessel ${edge.edge_index}`);
   const description=`Vessel ${edge.edge_index} · diameter ${number(edge.radius_m*2e6)} µm · endpoint pressure ${edge.endpoint_pressure_pa?.map(number).join(' → ')||'Unavailable'} Pa · flow ${number(edge.flow_m3_per_s)} m³/s.`;
   const title=svgEl('title');title.textContent=description;line.append(title);line.onclick=()=>inspect(item);line.onkeydown=e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();inspect(item);}};svg.append(line);
  }
  const barUm=10**Math.floor(Math.log10(projection.umPerUnit*100)),bar=svgEl('path'),barLabel=svgEl('text');
  bar.setAttribute('d',`M20 278H${20+barUm/projection.umPerUnit}`);bar.setAttribute('stroke','#edf4ee');bar.setAttribute('stroke-width','2');
  barLabel.setAttribute('x','20');barLabel.setAttribute('y','270');barLabel.setAttribute('fill','#edf4ee');barLabel.setAttribute('font-size','11');barLabel.textContent=`${number(barUm)} µm`;svg.append(bar,barLabel);
  marker=svgEl('circle');marker.setAttribute('r','4');marker.setAttribute('fill','#fff');marker.setAttribute('stroke','#111');marker.setAttribute('visibility','hidden');svg.append(marker);
  if(selectedEdge){const current=projection.edges.find(item=>item.edge.edge_id===selectedEdge.edge.edge_id);if(current)inspect(current);}
  plot.append(svg);
  scale.textContent=`Local axes ${projection.axes.map(i=>i+1).join('–')} · ${number(projection.umPerUnit)} µm per SVG view unit · centerline sampling ${number(response.patch.resolution_m*1e6)} µm. Diameter uses physical scale; sampling is numerical and biological resolution remains source-limited. Color: mean endpoint pressure ${Number.isFinite(low)?number(low)+'–'+number(high)+' Pa':'unavailable'}.`;
  status.textContent=`${response.selection?.name||REGIONS.find(r=>r[0]===region.value)?.[1]} · ${projection.edges.length} returned vessels · ${response.scenario?.label||response.scenario?.id||SCENARIO}.`;
 }
 plane.onchange=draw;zoom.onchange=draw;
 load.onclick=async()=>{
  if(disposed||busy)return;invalidate();const token=generation;busy=true;load.disabled=true;status.textContent='Loading the bounded local hydraulic graph…';
  try{
   const payload={entity_id:region.value,resolution_m:50e-6,scenario_id:scenario.value};
   if(skin(region.value)){payload.territory_id=territory.value;payload.patch_area_mm2=Number(area.value);if(!Number.isFinite(payload.patch_area_mm2)||payload.patch_area_mm2<.15||payload.patch_area_mm2>8)throw Error('Skin patch area must be 0.15–8 mm².');}
   const next=await request(payload);if(disposed||token!==generation)return;
   projectMicrovascularPatch(next);response=next;
   scope.textContent=[next.patch.registration?.evidence,next.patch.ownership,next.patch.boundary_scope,next.patch.uncertainty?.summary,next.patch.containment?.cortical_location_verified===false?`Whole-kidney authored-surface containment: ${next.patch.containment.whole_patch_inside_authored_surface===true?'certified':'unverified'}; cortical location remains unverified.`:null].filter(Boolean).join(' ');
   source.textContent=JSON.stringify({scenario:next.scenario,selection:next.selection,conditioning:next.patch.graph?.source_conditioning||next.patch.graph?.conditioning,uncertainty:next.patch.uncertainty,radius_conditioning:next.patch.graph?.radius_conditioning,hydraulic_model:next.patch.graph?.hydraulic_model,registration:next.patch.registration,containment:next.patch.containment,boundary_scope:next.patch.boundary_scope,ownership:next.patch.ownership},null,2);
   draw();
  }catch(error){if(!disposed&&token===generation){response=null;plot.replaceChildren();status.textContent=`Local graph unavailable: ${error.message}`;}}
  finally{if(token===generation){busy=false;load.disabled=disposed;}}
 };
 return {updateSelection(structure){
  if(disposed)return;invalidate();const found=REGIONS.find(r=>r[0]===structure?.id);
  if(found){if(region.value!==found[0]){region.value=found[0];scenarios();plane.value=skin(region.value)?'02':'01';}skinControls();selection.textContent=`Selected local region: ${found[1]}.`;}
  else selection.textContent=`${structure?.name?structure.name+': ':''}Available source-conditioned patches are vastus lateralis, kidney and exact skin territories. The region selector explicitly chooses the local patch.`;
  status.textContent='Load the named region to inspect its graph.';
 },clear(){if(disposed)return;invalidate();selection.textContent='Choose a local muscle, kidney or skin territory.';status.textContent='No local graph loaded.';},dispose(){disposed=true;invalidate();region.disabled=true;scenario.disabled=true;territory.disabled=true;area.disabled=true;plane.disabled=true;zoom.disabled=true;}};
}
