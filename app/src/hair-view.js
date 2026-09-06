import {createAsyncHairController} from './hair-worker-client.js';
import {bodyTransform} from './state.js';

export function attachElasticHair(object,geometry) {
  if(geometry.attachment?.kind!=='ElasticHairMaterialPoint')return false;
  object.userData.elasticHair?.dispose?.();
  delete object.userData.elasticHair;delete object.userData.hairWorkerToken;
  object.userData.elasticHairGeometry=geometry;
  object.userData.hairStaticReference=object.geometry.getAttribute('position').array.slice();
  return true;
}

// The controller returns canonical coordinates with the skin transform already
// applied. Applying the ordinary body matrix here would move every strand twice.
export function updateElasticHair(object,{frame,referenceCentroids={},recordKey,visible=true,enabled=false,
    gravity_m_s2=[0,-9.81,0],controllerFactory=createAsyncHairController,onDiagnostics=null}={}) {
  const geometry=object.userData.elasticHairGeometry;
  if(!geometry)return false;
  const positions=object.geometry.getAttribute('position');
  if(!enabled){
    if(object.userData.elasticHair){
      object.userData.elasticHair.dispose?.();
      delete object.userData.hairWorkerToken;
      positions.array.set(object.userData.hairStaticReference);positions.needsUpdate=true;
      object.geometry.computeVertexNormals();object.geometry.computeBoundingSphere();
      delete object.userData.elasticHair;
    }
    const id=geometry.attachment.skin_entity_id;
    object.matrixAutoUpdate=false;
    object.matrix.set(...bodyTransform(frame?.entities?.[id],referenceCentroids?.[id]));
    object.matrixWorldNeedsUpdate=true;
    object.userData.hairDiagnostics={mode:'static',physical_radius_multiplier:1,updated:false,
      local_skin_deformation:false,elastic_solver_enabled:false};
    return true;
  }
  object.userData.hairDiagnosticsListener=onDiagnostics;
  const apply=(result)=>{
    object.userData.hairDiagnostics={...object.userData.hairDiagnostics,...result.diagnostics,mode:'dynamic',elastic_solver_enabled:true};
    if(result.diagnostics.updated&&result.positions){
      positions.array.set(result.positions);positions.needsUpdate=true;
      object.geometry.computeVertexNormals();object.geometry.computeBoundingSphere();
      object.matrixAutoUpdate=false;object.matrix.identity();object.matrixWorldNeedsUpdate=true;
    }
    object.userData.hairDiagnosticsListener?.(object);
  };
  if(!object.userData.elasticHair){
    const token={};object.userData.hairWorkerToken=token;
    object.userData.elasticHair=controllerFactory(geometry,result=>{
      if(object.userData.hairWorkerToken===token)apply(result);
    });
  }
  object.userData.hairDiagnostics={...object.userData.hairDiagnostics,mode:'dynamic',elastic_solver_enabled:true,
    worker_pending:!!frame&&visible,paused:!frame||!visible};
  const id=geometry.attachment.skin_entity_id;
  // Send only skin inputs, not the complete whole-body frame on every update.
  const hairFrame=frame?{time_s:frame.time_s,entities:{[id]:frame.entities?.[id]},respiration:{skin_field:frame.respiration?.skin_field}}:null;
  const result=object.userData.elasticHair.update({frame:hairFrame,referenceCentroids:{[id]:referenceCentroids[id]},simulationTime:frame?.time_s??0,
    recordKey,visible,active:!!frame,gravity_m_s2});
  // Explicit controller injection supports deterministic source-only tests.
  if(result)apply(result);
  return true;
}
