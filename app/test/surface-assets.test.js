import {test} from 'node:test';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {SurfaceAssetCache} from '../src/surface-assets.js';
const segment={id:'bone',bone_id:'bone-id',reference_centroid_m:[0,0,0],bounds_min_m:[-1,-1,-1],bounds_max_m:[1,1,1]};
const payload={schema:'ihm.continuous-surface-binding.v1',binding_identity:'identity',segments:[segment],reference_positions_m:[[0,0,0]],weights:[[1]]};
function record(){const raw=JSON.stringify(payload),digest=createHash('sha256').update(raw).digest('hex');return {raw,binding:{binding_identity:'identity',segments:[segment],weights_sha256:digest,weights_url:'/api/surface-binding/'+digest}};}
test('surface sidecar is fetched once and verified against frozen bytes and support identity',async()=>{
 const {raw,binding}=record();let calls=0,ready=0;const cache=new SurfaceAssetCache({fetcher:async()=>{calls++;return new Response(raw);},onReady:()=>ready++});
 assert.equal(cache.get(binding),null);await cache.ensure(binding);assert.ok(cache.get(binding));await cache.ensure(binding);assert.equal(calls,1);assert.equal(ready,1);
 assert.throws(()=>cache.get({...binding,binding_identity:'other'}),/Conflicting/);
});
test('sidecar integrity failures remain explicit and arbitrary URLs are never fetched',async()=>{
 const {binding}=record();let error,calls=0;const cache=new SurfaceAssetCache({fetcher:async()=>{calls++;return new Response('{}');},onError:e=>error=e});
 await cache.ensure(binding);assert.match(error.message,/integrity/);assert.throws(()=>cache.get(binding),/integrity/);assert.equal(calls,1);
 const bad=new SurfaceAssetCache({fetcher:()=>{throw Error('must not fetch');}});await bad.ensure({...binding,weights_url:'/etc/passwd'});assert.throws(()=>bad.get({...binding,weights_url:'/etc/passwd'}),/address/);
});
