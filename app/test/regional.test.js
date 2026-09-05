import test from 'node:test';
import assert from 'node:assert/strict';
import {tetraSurface,RegionalView,ElectricRegionalView} from '../src/regional.js';
import * as THREE from 'three';
test('regional tetrahedral view removes shared interior faces without changing node IDs',()=>{
  const faces=tetraSurface([[0,1,2,3],[0,2,1,4]]);
  assert.equal(faces.length,18);
  for(let i=0;i<faces.length;i+=3)assert.notDeepEqual(faces.slice(i,i+3).sort((a,b)=>a-b),[0,1,2]);
  assert.throws(()=>tetraSurface([[0,1,2,3],[0,2,1,4],[0,1,2,5]]),/Nonmanifold/);
});
test('camera presets retain the active regional physical scale',()=>{
  for(const Type of [RegionalView,ElectricRegionalView]) {
    const camera=new THREE.PerspectiveCamera(),controls={target:new THREE.Vector3(),update(){}};
    const view=new Type({camera,controls});
    view.data={anchor:{origin_m:[.4,1.2,.1],local_axes:[[1,0,0],[0,1,0],[0,0,1]]}};
    for(const direction of ['front','side']){
      view.resetCamera(direction);
      assert.ok(camera.position.distanceTo(controls.target)<.04);
      assert.ok(camera.position.distanceTo(controls.target)>.001);
      assert.deepEqual(controls.target.toArray(),[.4,1.2,.1]);
    }
  }
});
