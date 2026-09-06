// One worker and one running solve across all regions. Each region holds at
// most one replacement request. Backpressure never invents mechanical time:
// the existing controller explicitly reports any gap it cannot simulate.
export function createHairWorkerPool(workerFactory=()=>new Worker(new URL('./hair-worker.js',import.meta.url),{type:'module'})) {
 let worker=null,inFlight=null,nextId=1,nextSequence=1,failure=null;
 const handles=new Map(),queue=new Set();
 const fail=(message)=>{
  failure=String(message);worker?.terminate();worker=null;inFlight=null;queue.clear();
  for(const handle of handles.values()){
   handle.pending=null;
   handle.onResult({positions:null,diagnostics:{updated:false,worker_pending:false,worker_error:failure,collision_model:'none'}});
  }
 };
 const pump=()=>{
  if(inFlight||failure||!queue.size)return;
  const id=queue.values().next().value;queue.delete(id);
  const handle=handles.get(id);if(!handle?.pending){pump();return;}
  const request=handle.pending;handle.pending=null;
  const sequence=nextSequence++;
  inFlight={id,sequence,request};
  const message={id,sequence,request,pausedSinceLast:handle.pausedSinceLast};
  handle.pausedSinceLast=false;
  if(!handle.initialized){message.geometry=handle.geometry;handle.initialized=true;}
  try {worker.postMessage(message);}catch(error){fail(error.message);}
 };
 const ensureWorker=()=>{
  if(worker||failure)return;
  try {
   worker=workerFactory();
   worker.onmessage=({data})=>{
    if(!inFlight||data.id!==inFlight.id||data.sequence!==inFlight.sequence)return;
    const completed=inFlight;inFlight=null;
    const handle=handles.get(data.id);
    if(handle){
     if(data.error){handle.failure=data.error;handle.pending=null;queue.delete(data.id);}
     const superseded=handle.pending&&(handle.pending.recordKey!==completed.request.recordKey||handle.pending.simulationTime<completed.request.simulationTime);
     if(!handle.hidden&&!superseded)handle.onResult({positions:data.positions||null,diagnostics:{...data.diagnostics,updated:!!data.diagnostics?.updated,worker_error:data.error||null,worker_pending:!!handle.pending,coalesced_updates:handle.coalesced}});
    }
    pump();
   };
   worker.onerror=(event)=>{event.preventDefault?.();fail(event.message||'Hair worker failed');};
   worker.onmessageerror=()=>fail('Hair worker message could not be decoded');
  }catch(error){fail(error.message);}
 };
 return {
  createController(geometry,onResult){
   const {strands,attachment,render_strands,render_attachment,guide_interpolation}=geometry;
   geometry={strands,attachment,render_strands,render_attachment,guide_interpolation};
   const id=nextId++,handle={geometry,onResult,pending:null,initialized:false,pausedSinceLast:false,hidden:false,coalesced:0};
   handles.set(id,handle);ensureWorker();
   return {
    update(request){
     if(!handles.has(id))return;
     if(failure||handle.failure){onResult({positions:null,diagnostics:{updated:false,worker_pending:false,worker_error:failure||handle.failure,collision_model:'none'}});return;}
     handle.hidden=!request.active||!request.visible;
     if(handle.hidden){handle.pausedSinceLast=true;handle.pending=null;queue.delete(id);return;}
     if(handle.pending)handle.coalesced++;
     handle.pending=request;queue.add(id);pump();
    },
    dispose(){
     handles.delete(id);queue.delete(id);
     if(worker)worker.postMessage({id,dispose:true});
     if(!handles.size){worker?.terminate();worker=null;inFlight=null;failure=null;}
    },
   };
  },
 };
}

let sharedPool;
export function createAsyncHairController(geometry,onResult){
 sharedPool??=createHairWorkerPool();
 return sharedPool.createController(geometry,onResult);
}
