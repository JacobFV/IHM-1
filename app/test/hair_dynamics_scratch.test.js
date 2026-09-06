import test from 'node:test';
import assert from 'node:assert/strict';
import {ElasticHairSystem} from '../src/hair_dynamics.js';

test('substeps avoid allocating a full coordinate copy and preserve checkpoint replay',()=>{
 const system=new ElasticHairSystem({
  centerlines_m:[0,0,0,.005,0,0,.01,0,0,.015,0,0],
  strand_offsets:[0,4],radius_m:[40e-6],
  tensile_modulus_pa:7.11e9,bending_modulus_pa:5.7e9,density_kg_m3:1312,
 });
 system.positions.slice=()=>{throw Error('Substep allocated a full coordinate copy');};
 system.step(.1);
 assert.ok(system.positions[10]<0,'Gravity must still bend the free tip');
 const checkpoint=system.checkpoint();
 const options={roots_m:[.001,.002,0],root_tangents:[1,.01,0]};
 system.step(.1,options);
 const expected=system.checkpoint();
 system.restore(checkpoint);
 system.step(.1,options);
 assert.deepEqual(system.checkpoint(),expected);
});
