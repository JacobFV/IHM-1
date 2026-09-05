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
