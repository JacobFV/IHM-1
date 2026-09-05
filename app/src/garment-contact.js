import * as THREE from 'three';

const demand=(condition,message)=>{if(!condition)throw Error('Garment contact: '+message);};
const vector=v=>Array.isArray(v)&&v.length===3&&v.every(Number.isFinite);
const indices=(v,n)=>Array.isArray(v)&&v.every(i=>Number.isInteger(i)&&i>=0&&i<n);
const resultant=forces=>Math.hypot(...forces.reduce((s,v)=>s.map((x,k)=>x+v[k]),[0,0,0]));
const displacement=(positions,reference)=>Math.max(...positions.map((v,i)=>Math.hypot(...v.map((x,k)=>x-reference[i][k]))));

export function prepareGarmentContact(data,garment,availableOwnerIds){
  demand(data?.schema_version===1&&data.frame?.units==='m'&&data.frame?.axes?.x==='left'&&data.frame?.axes?.y==='superior'&&data.frame?.axes?.z==='anterior','unsupported coordinate frame');
  const times=data.time_s;
  demand(Array.isArray(times)&&times.length>1&&times.every((t,i)=>Number.isFinite(t)&&(!i||t>times[i-1])),'invalid recorded clock');
  demand(data.frames?.length===times.length-1&&data.frames.every((f,i)=>f.time_s===times[i+1]),'metric clock differs from geometry clock');
  for(const name of ['tissue','panel']){
    const d=data[name],n=d?.positions_m?.[0]?.length;
    demand(n>2&&d.positions_m.length===times.length&&d.positions_m.every(p=>p.length===n&&p.every(vector)),name+' has invalid positions');
    demand(d.triangles?.length>0&&d.triangles.every(t=>t.length===3&&indices(t,n)&&new Set(t).size===3),name+' has invalid triangles');
    demand(d.contact_force_n?.length===times.length&&d.contact_force_n.every(f=>f.length===n&&f.every(vector)),name+' has invalid force observations');
  }
  const {tissue,panel}=data,owners=tissue.material_regions?.map(m=>m.source_id);
  demand(owners?.length>0&&new Set(owners).size===owners.length&&owners.every(id=>availableOwnerIds.has(id)),'source owners are absent or duplicated');
  demand(tissue.triangle_material_index?.length===tissue.triangles.length&&indices(tissue.triangle_material_index,owners.length)&&new Set(tissue.triangle_material_index).size===owners.length,'invalid surface ownership');
  demand(panel.source_garment_id===garment?.id,'source garment unavailable');
  demand(panel.source_node_indices?.length===panel.positions_m[0].length&&indices(panel.source_node_indices,garment.positions.length/3)&&new Set(panel.source_node_indices).size===panel.source_node_indices.length,'invalid garment node mapping');
  demand(panel.source_triangle_indices?.length===panel.triangles.length&&indices(panel.source_triangle_indices,garment.indices.length/3)&&new Set(panel.source_triangle_indices).size===panel.source_triangle_indices.length,'invalid garment face mapping');
  demand(panel.source_positions_m?.length===panel.source_node_indices.length&&panel.source_triangles?.length===panel.triangles.length,'missing retained garment correspondence');
  for(let i=0;i<panel.source_node_indices.length;i++){
    const source=panel.source_positions_m[i],id=panel.source_node_indices[i];
    demand(vector(source)&&source.every((v,k)=>Math.abs(v-garment.positions[3*id+k])<1e-12),'garment node coordinates changed');
  }
  for(let i=0;i<panel.triangles.length;i++){
    const original=garment.indices.slice(3*panel.source_triangle_indices[i],3*panel.source_triangle_indices[i]+3);
    demand(panel.source_triangles[i]?.length===3&&panel.source_triangles[i].every((v,k)=>v===original[k]),'garment face topology changed');
    demand([...panel.triangles[i].map(v=>panel.source_node_indices[v])].sort((a,b)=>a-b).join(',')===[...original].sort((a,b)=>a-b).join(','),'panel face does not correspond to its garment source');
  }
  const replaced=new Set(panel.source_triangle_indices),remainder_indices=garment.indices.filter((_,i)=>!replaced.has(Math.floor(i/3)));
  const frames=times.map((time_s,i)=>({time_s,...(i?data.frames[i-1]:{})}));
  const definitions=[
    ['Panel contact resultant','N',(_,i)=>resultant(panel.contact_force_n[i])],
    ['Tissue contact resultant','N',(_,i)=>resultant(tissue.contact_force_n[i])],
    ['Maximum tissue displacement','m',(_,i)=>displacement(tissue.positions_m[i],tissue.positions_m[0])],
    ['Maximum panel displacement','m',(_,i)=>displacement(panel.positions_m[i],panel.positions_m[0])],
    ['Minimum tissue Jacobian','1',f=>f.minimum_jacobian??null],
    ['Contact dissipation','J',f=>f.contact_dissipation_j??null],
    ['Numerical energy defect','J',f=>f.numerical_energy_defect_j??null],
    ['Tissue interface work','J',f=>f.tissue_interface_work_j??null],
    ['Panel interface work','J',f=>f.panel_interface_work_j??null],
    ['Unresolved edge samples','count',f=>f.unresolved_nearest_edge_samples??null],
  ];
  return {frames,owners,remainder_indices,physiology:{time_s:times,values:Object.fromEntries(definitions.map(([name,,read])=>[name,frames.map(read)])),units:Object.fromEntries(definitions.map(([name,unit])=>[name,unit])),metadata:{source_kind:'Computed local garment–tissue contact · explicit engineering priors'}}};
}

export class GarmentContactView {
  constructor({bodyGroup,objects,clothingView,ownerVisible,camera,controls,availableOwnerIds}){
    Object.assign(this,{bodyGroup,objects,clothingView,ownerVisible,camera,controls,availableOwnerIds});this.active=false;
  }
  open(source){
    const garment=this.clothingView?.materializations.find(g=>g.id===source.panel.source_garment_id);
    this.prepared=prepareGarmentContact(source,garment,this.availableOwnerIds??new Set(this.objects.keys()));
    this.source=source;this.data={...source,frames:this.prepared.frames};
    this.garmentMesh=this.clothingView.meshes.get(garment.id);
    demand(this.garmentMesh,'garment mesh unavailable');
    this.originalIndex=this.garmentMesh.geometry.index;
    this.garmentMesh.geometry.setIndex(this.prepared.remainder_indices);
    this.group=new THREE.Group();this.group.name='computed-garment-contact';this.bodyGroup.add(this.group);
    const make=(positions,triangles,color)=>{
      const geometry=new THREE.BufferGeometry();geometry.setAttribute('position',new THREE.Float32BufferAttribute(positions.flat(),3));geometry.setIndex(triangles.flat());geometry.computeVertexNormals();
      const mesh=new THREE.Mesh(geometry,new THREE.MeshStandardMaterial({color,roughness:.8,side:THREE.DoubleSide}));this.group.add(mesh);return mesh;
    };
    this.tissueMeshes=source.tissue.material_regions.map((region,i)=>{
      const mesh=make(source.tissue.positions_m[0],source.tissue.triangles.filter((_,f)=>source.tissue.triangle_material_index[f]===i),['#b47c80','#c49b83','#d0a49b'][i%3]);
      mesh.name=region.source_id;return mesh;
    });
    this.panelMesh=make(source.panel.positions_m[0],source.panel.triangles,'#486d83');this.panelMesh.name='computed-shorts-front-panel';
    this.savedVisibility=new Map(this.prepared.owners.map(id=>[id,this.objects.get(id)?.visible]));
    if(this.camera)this.savedCamera={position:this.camera.position.clone(),target:this.controls.target.clone(),near:this.camera.near,far:this.camera.far};
    this.active=true;this.draw(0);
  }
  syncOwners(){
    if(!this.active)return;
    this.prepared.owners.forEach((id,i)=>{
      const original=this.objects.get(id),visible=this.ownerVisible?.(id)??this.savedVisibility.get(id);
      if(original)original.visible=false;
      this.tissueMeshes[i].visible=visible;
    });
    this.panelMesh.visible=this.garmentMesh.visible;
  }
  draw(index){
    demand(Number.isInteger(index)&&index>=0&&index<this.data.frames.length,'unrecorded frame');
    for(const [meshes,positions] of [[this.tissueMeshes,this.source.tissue.positions_m[index]],[[this.panelMesh],this.source.panel.positions_m[index]]]){
      const flattened=positions.flat();for(const mesh of meshes){mesh.geometry.attributes.position.array.set(flattened);mesh.geometry.attributes.position.needsUpdate=true;mesh.geometry.computeVertexNormals();mesh.geometry.computeBoundingSphere();}
    }
    this.syncOwners();return this.data.frames[index];
  }
  resetCamera(direction='front'){
    if(!this.camera)return;
    const box=new THREE.Box3().setFromObject(this.group),center=box.getCenter(new THREE.Vector3()),distance=Math.max(.3,box.getSize(new THREE.Vector3()).length()*2);
    this.controls.target.copy(center);this.camera.position.copy(center).add(direction==='side'?new THREE.Vector3(distance,0,0):new THREE.Vector3(0,0,distance));
    this.camera.near=.0001;this.camera.far=10;this.camera.updateProjectionMatrix();this.controls.update();
  }
  setClipping(planes){if(this.active)this.group.children.forEach(m=>{m.material.clippingPlanes=planes;});}
  close(){
    if(!this.active)return;
    this.garmentMesh.geometry.setIndex(this.originalIndex);
    this.prepared.owners.forEach(id=>{const obj=this.objects.get(id);if(obj)obj.visible=this.ownerVisible?.(id)??this.savedVisibility.get(id);});
    this.group.children.forEach(m=>{m.geometry.dispose();m.material.dispose();});this.group.removeFromParent();
    if(this.savedCamera){this.camera.position.copy(this.savedCamera.position);this.controls.target.copy(this.savedCamera.target);this.camera.near=this.savedCamera.near;this.camera.far=this.savedCamera.far;this.camera.updateProjectionMatrix();this.controls.update();}
    this.active=false;
  }
}
