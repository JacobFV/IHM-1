import * as THREE from 'three';

// Material roles the domain builders record. Colour is a legend, not a measurement.
export const ROLE_COLORS={rigid_bone:'#e9e3d2',cartilage:'#bcd6e0',muscle:'#b0453f',ligament:'#d6c268',
  tendon:'#e2d9ab',connective_tissue:'#8fae95',body_envelope:'#d7a389',skin:'#d7a389',soft_organ:'#a86279',
  vascular:'#8f333f',nerve:'#e3d264',lymph_node_group:'#93a86b',fluid_cavity:'#6fa3b5',unclaimed_complement:'#55666f'};
export const ROLE_LABELS={rigid_bone:'Bone',cartilage:'Cartilage',muscle:'Muscle',ligament:'Ligament',
  tendon:'Tendon',connective_tissue:'Connective / tendon sheath',body_envelope:'Body envelope (interstitium)',
  skin:'Skin',soft_organ:'Soft organ',vascular:'Vascular',nerve:'Nerve',lymph_node_group:'Lymph node group',
  fluid_cavity:'Fluid cavity',unclaimed_complement:'Unclaimed complement'};
// The envelope and the complement enclose everything else, so they start hidden.
const ENCLOSING=new Set(['body_envelope','unclaimed_complement']);

// One conforming tetrahedral volume, shown as the boundary surface each material owner
// carries inside the shared vertex array. No interior tetrahedron reaches the browser and
// no vertex is moved, so what is drawn is exactly the recorded ownership partition.
export class DomainView {
  constructor({scene,camera,controls,bodyGroup}) {
    Object.assign(this,{scene,camera,controls,bodyGroup});this.active=false;this.meshes=[];
  }
  // Frame the whole domain: the bounding sphere has to fit the narrower of the two field angles.
  get viewDistance(){
    const radius=new THREE.Vector3(...this.data.anchor.extent_m).length()/2;
    const vertical=THREE.MathUtils.degToRad(this.camera.fov);
    return 1.08*radius/Math.sin(Math.min(vertical,2*Math.atan(Math.tan(vertical/2)*this.camera.aspect))/2);
  }
  resetCamera(direction='front'){
    const center=new THREE.Vector3(...this.data.anchor.origin_m),distance=this.viewDistance;
    this.controls.target.copy(center);
    this.camera.position.copy(center).add(direction==='side'?new THREE.Vector3(distance,.15*distance,0)
      :new THREE.Vector3(.15*distance,.15*distance,distance));
    this.camera.near=distance/5000;this.camera.far=100*distance;this.camera.updateProjectionMatrix();this.controls.update();
  }
  open(data) {
    this.close();
    this.data=data;
    // One position buffer for the whole volume: owners that meet in the mesh stay coincident.
    const position=new THREE.Float32BufferAttribute(data.geometry.positions_m,3);
    this.object=new THREE.Group();
    for(const group of data.groups) {
      const geometry=new THREE.BufferGeometry();
      geometry.setAttribute('position',position);
      geometry.setIndex(new THREE.BufferAttribute(new Uint32Array(group.triangles),1));
      geometry.computeVertexNormals();
      const mesh=new THREE.Mesh(geometry,new THREE.MeshStandardMaterial({
        color:ROLE_COLORS[group.role]||'#9aa5ac',roughness:.85,metalness:0,side:THREE.DoubleSide,
        transparent:true,opacity:1}));
      mesh.name=group.role;mesh.userData.group=group;
      mesh.visible=!ENCLOSING.has(group.role);
      this.object.add(mesh);this.meshes.push(mesh);
    }
    this.scene.add(this.object);
    this.saved={position:this.camera.position.clone(),target:this.controls.target.clone(),near:this.camera.near,far:this.camera.far};
    this.resetCamera();this.bodyGroup.visible=false;this.active=true;
  }
  roles(){return this.meshes.map(mesh=>({role:mesh.name,visible:mesh.visible,
    triangles:mesh.userData.group.triangle_count,structures:mesh.userData.group.structures.length}));}
  setRoleVisible(role,visible){const mesh=this.meshes.find(m=>m.name===role);if(mesh)mesh.visible=visible;}
  setOpacity(value){for(const mesh of this.meshes){mesh.material.opacity=value;mesh.material.transparent=value<1;mesh.material.depthWrite=value>.5;}}
  setClipping(planes){for(const mesh of this.meshes)mesh.material.clippingPlanes=planes;}
  // Which owner carries the triangle under the cursor.
  pick(ray){
    const hit=ray.intersectObjects(this.meshes.filter(m=>m.visible),false)[0];
    if(!hit)return null;
    const group=hit.object.userData.group;
    const found=group.structures.find(s=>hit.faceIndex>=s.triangle_start&&hit.faceIndex<s.triangle_start+s.triangle_count);
    return found?{role:group.role,structure:found,point:hit.point}:null;
  }
  draw(index){return this.data.frames[index];}
  close(){
    if(this.object){
      for(const mesh of this.meshes){mesh.geometry.dispose();mesh.material.dispose();}
      this.scene.remove(this.object);this.object=null;this.meshes=[];
    }
    if(this.saved){this.camera.position.copy(this.saved.position);this.controls.target.copy(this.saved.target);
      this.camera.near=this.saved.near;this.camera.far=this.saved.far;this.camera.updateProjectionMatrix();this.controls.update();}
    this.bodyGroup.visible=true;this.saved=null;this.data=null;this.active=false;
  }
}
