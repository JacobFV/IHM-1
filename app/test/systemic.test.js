import {test} from 'node:test';
import assert from 'node:assert/strict';
import {prepareSystemic,systemicLabel,systemicSnapshot} from '../src/systemic.js';
import {spectralSeries} from '../src/state.js';
const fixture=()=>({model_id:'ihm-body',has_body_projection:false,configuration:{protocol:'meal',sample_interval_s:30},clock:{sample_interval_s:30},centroids_m:{},
  fields:{glucose:{unit:'mg/dL'},missing:{unit:'g'}},frames:[{time_s:0,entities:{},systemic_values:{glucose:90,missing:null}},{time_s:30,entities:{},systemic_values:{glucose:110,missing:null}}],
  spectra:{frequency_hz:[0,.01],sigma_per_s:[0,.1],variables:[{name:'glucose',removed_sample_mean:100,psd:{frequency_hz:[0,.01],psd:[0,2]},laplace_real:[[0,3],[0,6]],laplace_imag:[[0,4],[0,8]]}]}});
test('systemic signal and spectra preserve native units, missingness and distinct frequency axes',()=>{
  const data=fixture(),result=prepareSystemic(data,'meal');
  assert.deepEqual(result.physiology.values.glucose,[90,110]);assert.deepEqual(result.physiology.values.missing,[null,null]);
  assert.deepEqual(spectralSeries(result.spectralRun,'glucose','laplace',1).y,[0,10]);
  assert.equal(spectralSeries(result.spectralRun,'glucose','psd').unit,'(mg/dL)²/Hz');
  assert.equal(result.spectralRun.laplace.variables[0].unit,'mg/dL · s');
  assert.equal(systemicSnapshot(data,0,'missing').value,null);
  assert.equal(systemicSnapshot(data,1,'glucose').value,110);
  assert.equal(systemicSnapshot(data,20,'glucose').value,null);
  assert.match(systemicLabel('Aorta.Glucose.concentration_mg_per_dl'),/Aorta/);
});
test('sparse systemic records cannot smuggle uncomputed body motion or irregular clocks',()=>{
  const moving=fixture();moving.frames[1].entities={skin:{translation_m:[0,0,.1]}};
  assert.throws(()=>prepareSystemic(moving,'meal'),/projection/);
  const bad=fixture();bad.frames[1].time_s=0;assert.throws(()=>prepareSystemic(bad,'meal'),/clock/);
  const missing=fixture();missing.frames[1].systemic_values.glucose='110';assert.throws(()=>prepareSystemic(missing,'meal'),/numeric/);
  const flag=fixture();flag.has_body_projection='false';assert.throws(()=>prepareSystemic(flag,'meal'),/projection/);
});
