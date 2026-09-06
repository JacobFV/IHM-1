import test from 'node:test';
import assert from 'node:assert/strict';
import {bodyCommand,skinPressureCapability} from '../src/embodied-live.js';
const owner=(whole_skin,regional_ids=[])=>({input_capabilities:{skin_pressure:{whole_skin,regional_ids}}});
const whole=owner(true),regional=owner(false,['region_a','region_b','residual']);
test('whole-owner explicit pressure and release remain whole commands',()=>{
 for(const pressure of [0,250,5000])assert.equal(bodyCommand('embodied',1,[],{skin_compression_pa:pressure},whole).skin_compression_pa,pressure);
});
test('regional default zero produces no whole-skin topology command',()=>{
 const result=bodyCommand('embodied',1,[],{skin_compression_pa:0},regional);
 assert.ok(!('skin_compression_pa' in result));assert.ok(!('regional_skin_pressures' in result));
 assert.throws(()=>bodyCommand('embodied',1,[],{skin_compression_pa:1},regional),/Whole-skin pressure is unavailable/);
});
test('regional boundaries preserve exact owner IDs and cannot combine topologies',()=>{
 const inputs={skin_compression_pa:0,regional_skin_pressures:{region_b:250,residual:0}};
 assert.deepEqual(bodyCommand('embodied',1,[],inputs,regional).regional_skin_pressures,inputs.regional_skin_pressures);
 assert.throws(()=>bodyCommand('embodied',1,[],{...inputs,skin_compression_pa:1},regional),/Whole-skin pressure is unavailable/);
 for(const pressures of [{invented:1},{region_a:NaN},{region_a:-1},{region_a:5001},[]])assert.throws(()=>bodyCommand('embodied',1,[],{regional_skin_pressures:pressures},regional));
 assert.throws(()=>bodyCommand('embodied',1,[],inputs,whole),/Regional skin pressure is unavailable/);
});
test('missing and malformed capabilities never guess from observation or UI defaults',()=>{
 for(const frame of [undefined,{},owner('true'),owner(false,['x','x']),owner(true,['x']),{tissue_exchange:{regional_skin:{available:true}}}]){
  assert.equal(skinPressureCapability(frame).topology,'unavailable');
  assert.ok(!('skin_compression_pa' in bodyCommand('embodied',1,[],{},frame)));
  assert.throws(()=>bodyCommand('embodied',1,[],{skin_compression_pa:30},frame),/unavailable/);
 }
 assert.match(skinPressureCapability(regional).message,/regional/);
});
