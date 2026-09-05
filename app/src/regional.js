import * as THREE from 'three';

export function tetraSurface(tetrahedra) {
  const faces=new Map();
  for(const tet of tetrahedra) for(const local of [[0,2,1],[0,1,3],[0,3,2],[1,2,3]]) {
    const face=local.map(i=>tet[i]),key=[...face].sort((a,b)=>a-b).join(',');
    const old=faces.get(key);faces.set(key,old?{...old,count:old.count+1}:{face,count:1});
  }
  if([...faces.values()].some(f=>f.count>2)) throw Error('Nonmanifold regional volume');
  return [...faces.values()].filter(f=>f.count===1).flatMap(f=>f.face);
}

// The regional mesh uses the body's pinned material anchor and its own computed
// mechanical state. Its mass is not added to the global replay a second time.
export class RegionalView {
  constructor({scene,camera,controls,bodyGroup}) {
    Object.assign(this,{scene,camera,controls,bodyGroup});this.active=false;
  }
  get viewDistance(){return .026;}
  resetCamera(direction='front'){
    const center=new THREE.Vector3(...this.data.anchor.origin_m),axes=this.data.anchor.local_axes;
    this.controls.target.copy(center);
    this.camera.position.copy(center).addScaledVector(new THREE.Vector3(...axes[direction==='side'?0:2]),this.viewDistance)
      .addScaledVector(new THREE.Vector3(...axes[1]),-.6*this.viewDistance);
    this.camera.near=this.viewDistance/2500;this.camera.far=10;this.camera.updateProjectionMatrix();this.controls.update();
  }
  open(data) {
    this.close();
    this.data=data;
    const positions=data.geometry.reference_positions_m.flat();
    const geometry=new THREE.BufferGeometry();
    geometry.setAttribute('position',new THREE.Float32BufferAttribute(positions,3));
    geometry.setIndex(tetraSurface(data.geometry.tetrahedra));geometry.computeVertexNormals();
    this.object=new THREE.Mesh(geometry,new THREE.MeshStandardMaterial({color:'#ce9980',roughness:.8,side:THREE.DoubleSide}));
    const edges=new THREE.LineSegments(new THREE.EdgesGeometry(geometry,1),new THREE.LineBasicMaterial({color:'#365b58'}));
    this.scene.add(this.object);this.edges=edges;this.scene.add(edges);
    this.saved={position:this.camera.position.clone(),target:this.controls.target.clone(),near:this.camera.near,far:this.camera.far};
    this.resetCamera();this.bodyGroup.visible=false;this.active=true;
    this.draw(0);
  }
  draw(index) {
    const frame=this.data.frames[index];if(!frame)return;
    if(this.state!==frame.mechanical_state) {
      const positions=this.data.states[frame.mechanical_state].positions_m.flat();
      this.object.geometry.attributes.position.array.set(positions);
      this.object.geometry.attributes.position.needsUpdate=true;
      this.object.geometry.computeVertexNormals();this.object.geometry.computeBoundingSphere();
      this.edges.geometry.dispose();this.edges.geometry=new THREE.EdgesGeometry(this.object.geometry,1);
      this.state=frame.mechanical_state;
    }
    return frame;
  }
  close() {
    for(const object of [this.object,this.edges]) if(object){this.scene.remove(object);object.geometry.dispose();object.material.dispose();}
    if(this.saved){this.camera.position.copy(this.saved.position);this.controls.target.copy(this.saved.target);this.camera.near=this.saved.near;this.camera.far=this.saved.far;this.camera.updateProjectionMatrix();this.controls.update();}
    this.bodyGroup.visible=true;this.object=this.edges=this.saved=this.state=null;this.active=false;
  }
}

export class ElectricRegionalView extends RegionalView {
  get viewDistance(){return .0025;}
  open(data) {
    this.close();this.source=data;
    this.data={...data,frames:data.experiments.wound_shunt.frames};
    this.saved={position:this.camera.position.clone(),target:this.controls.target.clone(),near:this.camera.near,far:this.camera.far};
    this.object=new THREE.Group();this.scene.add(this.object);
    this.nodes=[];
    for(const [kind,positions] of [['membrane',data.geometry.cell_positions_m],['apical',data.geometry.surface_positions_m]]) {
      const mesh=new THREE.InstancedMesh(new THREE.SphereGeometry(.000025,12,8),new THREE.MeshStandardMaterial({roughness:.6}),positions.length);
      for(let i=0;i<positions.length;i++)mesh.setMatrixAt(i,new THREE.Matrix4().makeTranslation(...positions[i]));
      this.object.add(mesh);this.nodes.push({kind,mesh});
    }
    const edgePositions=[];
    for(const [edges,positions] of [[data.geometry.gap_edges,data.geometry.cell_positions_m],[data.geometry.surface_edges,data.geometry.surface_positions_m]])
      for(const [a,b] of edges)edgePositions.push(...positions[a],...positions[b]);
    const edges=new THREE.LineSegments(new THREE.BufferGeometry().setAttribute('position',new THREE.Float32BufferAttribute(edgePositions,3)),new THREE.LineBasicMaterial({color:'#809b98'}));
    this.object.add(edges);
    this.resetCamera();
    this.bodyGroup.visible=false;this.active=true;this.condition='wound_shunt';this.draw(0);
  }
  selectCondition(condition){this.condition=condition;this.data.frames=this.source.experiments[condition].frames;}
  draw(index){
    const frame=this.data.frames[index];if(!frame)return;
    // Each electrical quantity has a fixed range across conditions and time.
    for(const {kind,mesh} of this.nodes){
      const values=frame[kind==='membrane'?'membrane_voltage_V':'apical_voltage_V'];
      const low=kind==='membrane'?-.08:-.04,high=kind==='membrane'?.02:.01;
      for(let i=0;i<values.length;i++)mesh.setColorAt(i,new THREE.Color().setHSL(.65*(1-Math.max(0,Math.min(1,(values[i]-low)/(high-low)))),.65,.6));
      mesh.instanceColor.needsUpdate=true;
    }
    return frame;
  }
  close(){
    if(this.object?.isGroup){this.object.traverse(o=>{o.geometry?.dispose();o.material?.dispose();});this.scene.remove(this.object);this.object=null;}
    super.close();
  }
}
