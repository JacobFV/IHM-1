import {prepareContinuousAsset} from './continuous-surface.js';
const LIMIT=64*1024*1024;
export class SurfaceAssetCache {
 constructor({fetcher=(...args)=>fetch(...args),onReady=()=>{},onError=()=>{}}={}){this.fetcher=fetcher;this.onReady=onReady;this.onError=onError;this.records=new Map();}
 async read(binding){
  const digest=binding.weights_sha256,url=binding.weights_url;
  if(!/^[0-9a-f]{64}$/.test(digest)||url!==`/api/surface-binding/${digest}`)throw Error('Invalid retained surface asset address');
  const response=await this.fetcher(url,{cache:'force-cache'});
  if(!response.ok)throw Error(`Surface binding asset unavailable (${response.status})`);
  const length=Number(response.headers.get('Content-Length'));
  if(length>LIMIT)throw Error('Surface binding asset exceeds byte envelope');
  const reader=response.body.getReader(),chunks=[];let count=0;
  try{while(true){const {value,done}=await reader.read();if(done)break;count+=value.byteLength;if(count>LIMIT)throw Error('Surface binding asset exceeds byte envelope');chunks.push(value);}}
  catch(error){await reader.cancel().catch(()=>{});throw error;}
  const bytes=new Uint8Array(count);let offset=0;for(const chunk of chunks){bytes.set(chunk,offset);offset+=chunk.length;}
  const hash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes))).map(v=>v.toString(16).padStart(2,'0')).join('');
  if(hash!==digest)throw Error('Surface binding asset integrity mismatch');
  const payload=JSON.parse(new TextDecoder().decode(bytes));
  if(payload.binding_identity!==binding.binding_identity)throw Error('Surface binding identity changed');
  if(JSON.stringify(payload.segments)!==JSON.stringify(binding.segments))throw Error('Surface binding segment support changed');
  return prepareContinuousAsset(payload);
 }
 ensure(binding){
  const key=binding.weights_sha256;let record=this.records.get(key);
  if(!record){record={identity:binding.binding_identity};this.records.set(key,record);record.promise=this.read(binding).then(asset=>{record.asset=asset;try{this.onReady();}catch(error){record.error=error;this.onError(error);}return asset;},error=>{record.error=error;this.onError(error);return null;});}
  if(record.identity!==binding.binding_identity)throw Error('Conflicting retained surface identity');
  return record.promise;
 }
 get(binding){this.ensure(binding);const record=this.records.get(binding.weights_sha256);if(record.error)throw record.error;return record.asset||null;}
}
export const surfaceAssets=new SurfaceAssetCache();
export function configureSurfaceAssets({onReady,onError}){surfaceAssets.onReady=onReady;surfaceAssets.onError=onError;}
