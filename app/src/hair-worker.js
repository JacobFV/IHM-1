import {createHairController} from './hair_dynamics.js';

// Also exported for a tiny fixture test without browser assets.
export function createHairWorkerHandler(send) {
 const controllers=new Map();
 return ({data})=>{
  const {id,sequence,geometry,request,pausedSinceLast,dispose}=data;
  if(dispose){controllers.delete(id);return;}
  try {
   if(!controllers.has(id)){
    if(!geometry)throw Error('Missing hair worker geometry');
    controllers.set(id,{controller:createHairController(geometry),lastTime:null,unsimulatedTotal:0});
   }
   const state=controllers.get(id),{controller}=state;
   if(pausedSinceLast)controller.update({...request,active:false});
   const result=controller.update(request);
   const elapsed=state.lastTime===null?0:Math.max(0,request.simulationTime-state.lastTime);
   const reset=result.diagnostics.updated&&!!result.diagnostics.reset_reason;
   const integrated=result.diagnostics.updated&&!reset?elapsed:0;
   const unsimulated=reset?elapsed:0;
   state.unsimulatedTotal+=unsimulated;
   if(result.diagnostics.updated)state.lastTime=request.simulationTime;
   result.diagnostics={...result.diagnostics,requested_time_s:request.simulationTime,
    integrated_interval_s:integrated,unsimulated_interval_s:unsimulated,
    unsimulated_total_s:state.unsimulatedTotal,
    physical_interval_status:reset?'not_integrated':result.diagnostics.updated?'integrated':'held',
    input_sampling:'interpolated_roots_endpoint_forces',coupling:'one_way_no_reaction'};
   // Transfer a copy; the controller owns and reuses its cached tube buffer.
   const positions=result.diagnostics.updated?result.positions.slice():null;
   send({id,sequence,positions,diagnostics:result.diagnostics},positions?[positions.buffer]:[]);
  } catch(error) {
   send({id,sequence,error:String(error?.message||error)});
  }
 };
}

if(typeof self!=='undefined'&&typeof document==='undefined')self.onmessage=createHairWorkerHandler((message,transfer)=>self.postMessage(message,transfer));
