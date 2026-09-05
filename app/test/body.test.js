import test from 'node:test';
import assert from 'node:assert/strict';
import { defaultModelId, bodyTransform, validBodyTrajectory } from '../src/state.js';

test('canonical body is the default while source families remain available', () => {
  assert.equal(defaultModelId([{id:'bodyparts3d'}, {id:'ihm-body'}]), 'ihm-body');
  assert.equal(defaultModelId([{id:'source-only'}]), 'source-only');
});
test('body deformation acts about the supplied tissue centroid in RF order', () => {
  const transform = bodyTransform({translation_m:[1,2,3], rotation_matrix:[[0,-1,0],[1,0,0],[0,0,1]], deformation_gradient:[[2,0,0],[0,1,0],[0,0,.5]]}, [4,5,6]);
  const apply = p => [0,1,2].map(i => transform[i*4]*p[0]+transform[i*4+1]*p[1]+transform[i*4+2]*p[2]+transform[i*4+3]);
  assert.deepEqual(apply([4,5,6]), [5,7,9]);
  assert.deepEqual(apply([5,5,6]), [5,9,9]);
});
test('sparse frames reset missing tissues, invalid motion cannot fabricate geometry', () => {
  assert.deepEqual(bodyTransform(undefined), [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]);
  assert.throws(()=>bodyTransform({translation_m:[1,2,3]}));
  assert.throws(()=>bodyTransform({deformation_gradient:[[0,0,0],[0,0,0],[0,0,0]]},[0,0,0]));
  assert.equal(validBodyTrajectory({frames:[]}), false);
  assert.equal(validBodyTrajectory({frames:[{time_s:1,entities:{}},{time_s:0,entities:{}}],centroids_m:{}}), false);
  assert.equal(validBodyTrajectory({frames:[{time_s:0,entities:{x:{translation_m:[1,0,0]}}}],centroids_m:{}}), false);
  assert.equal(validBodyTrajectory({frames:[{time_s:0,entities:{x:{translation_m:[1,0,0]}}}],centroids_m:{x:[0,0,0]}}), true);
});

test('computed skin displacement tapers at the thorax boundary and never moves the back or accumulates across frames', async () => {
  const { deformSkinVertices } = await import('../src/state.js');
  const field = {center_m:[0,0,0], bounds_m:{min:[-2,-2,-2],max:[2,2,2]},reference_radii_m:[1,1],thorax_y_offsets_m:[-1,1],displacement_m:[.1,.2,.3]};
  const reference = new Float32Array([1,0,1,-1,0,-1,.5,1.5,0,2,0,1,1,2,1]);
  const output = new Float32Array(reference.length);
  deformSkinVertices(reference,field,output);
  const expected = [1.1,0,1.2,-1.1,0,-1,.525,1.5,.05,2,0,1,1,2,1];
  output.forEach((value,i) => assert.ok(Math.abs(value-expected[i])<1e-6));
  const first = output.slice();
  deformSkinVertices(reference,field,output);
  assert.deepEqual(output,first);
  deformSkinVertices(reference,null,output);
  assert.deepEqual(output,reference);
  assert.deepEqual(Array.from(reference),[1,0,1,-1,0,-1,.5,1.5,0,2,0,1,1,2,1]);
});

test('skin field is translation invariant and invalid field metadata cannot enter playback', async () => {
  const { deformSkinVertices } = await import('../src/state.js');
  const field={entity_ids:['skin'],center_m:[10,20,30],bounds_m:{min:[8,18,28],max:[12,22,32]},reference_radii_m:[1,1],thorax_y_offsets_m:[-1,1],displacement_m:[.1,.2,.3]};
  const output=deformSkinVertices([11,20,31],field);
  assert.ok(Math.abs(output[0]-11.1)<1e-5);
  assert.ok(Math.abs(output[2]-31.2)<1e-5);
  const trajectory={centroids_m:{},frames:[{time_s:0,entities:{},respiration:{skin_field:field}}]};
  assert.equal(validBodyTrajectory(trajectory),true);
  for (const broken of [{reference_radii_m:[0,1]},{displacement_m:[NaN,0,0]},{thorax_y_offsets_m:[-3,3]}]) {
    assert.throws(()=>deformSkinVertices([11,20,31],{...field,...broken}));
    assert.equal(validBodyTrajectory({...trajectory,frames:[{...trajectory.frames[0],respiration:{skin_field:{...field,...broken}}}]}),false);
  }
});
