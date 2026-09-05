import test from 'node:test';
import assert from 'node:assert/strict';
import {attachedHairPositions} from '../src/hair_motion.js';
test('material hair transport preserves shaft length and restores exact reference',()=>{
  const root=[0,0,1],reference=new Float32Array([-0.01,0,1,0,.01,1,.01,0,1,0,-.01,1,-.01,0,1.1,0,.01,1.1,.01,0,1.1,0,-.01,1.1]);
  const attachment={vertices_per_shaft:8,barycentric:[[.25,.25,.5]],reference_triangles_m:[-1,-1,1,1,-1,1,0,1,1]};
  const field={center_m:[0,0,0],bounds_m:{min:[-2,-2,-2],max:[2,2,2]},reference_radii_m:[1,1],thorax_y_offsets_m:[-1,1],displacement_m:[.2,.3,.1]};
  const moved=attachedHairPositions(reference,attachment,field);
  const length=Math.hypot(...[0,1,2].map(k=>moved[12+k]-moved[k]));
  assert.ok(Math.abs(length-.1)<1e-6);
  assert.notDeepEqual(moved,reference);
  assert.deepEqual(attachedHairPositions(reference,attachment,null,moved),reference);
});
