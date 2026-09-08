import * as THREE from 'three';

// Small geometric leaves fill the same declared canopy envelope used by contact.
// These are material detail; branches/canopy contact stays server-owned.
export function canopyLeaves(part) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position',new THREE.Float32BufferAttribute([
    0,0,0, -.018,.03,.005, .018,.03,.005, -.014,.07,0, .014,.07,0, 0,.105,-.009,
  ],3));
  geometry.setIndex([0,2,1,1,2,4,1,4,3,3,4,5]);geometry.computeVertexNormals();
  const material=new THREE.MeshPhysicalMaterial({color:0x65833a,roughness:.83,side:THREE.DoubleSide,sheen:.25,sheenColor:new THREE.Color(0xa1ad56)});
  const leaves=new THREE.InstancedMesh(geometry,material,1400);
  const object=new THREE.Object3D(),color=new THREE.Color();let seed=741;
  const random=()=>{seed=(1664525*seed+1013904223)>>>0;return seed/4294967296;};
  for(let i=0;i<leaves.count;i++) {
    const theta=random()*Math.PI*2,z=random()*2-1,r=part.radius_m*Math.cbrt(random())*.91;
    const xy=Math.sqrt(1-z*z);
    object.position.set(part.centre_m[0]+r*xy*Math.cos(theta),part.centre_m[1]+r*z,part.centre_m[2]+r*xy*Math.sin(theta));
    object.rotation.set(random()*Math.PI,random()*Math.PI*2,random()*Math.PI*2);
    object.scale.setScalar(.7+random()*.65);object.updateMatrix();leaves.setMatrixAt(i,object.matrix);
    color.setHSL(.22+random()*.06,.3+random()*.25,.22+random()*.15);leaves.setColorAt(i,color);
  }
  leaves.castShadow=true;leaves.receiveShadow=true;leaves.raycast=()=>{};
  return leaves;
}
