import {createHairController} from './hair_dynamics.js';
import {bodyTransform} from './state.js';

export function attachElasticHair(object,geometry) {
  if(geometry.attachment?.kind!=='ElasticHairMaterialPoint')return false;
  object.userData.elasticHairGeometry=geometry;
  object.userData.hairStaticReference=object.geometry.getAttribute('position').array.slice();
  return true;
}

// The controller returns canonical coordinates with the skin transform already
// applied. Applying the ordinary body matrix here would move every strand twice.
export function updateElasticHair(object,{frame,referenceCentroids={},recordKey,visible=true,enabled=false,
    gravity_m_s2=[0,-9.81,0]}={}) {
  const geometry=object.userData.elasticHairGeometry;
  if(!geometry)return false;
  const positions=object.geometry.getAttribute('position');
  if(!enabled){
    if(object.userData.elasticHair){
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
  const controller=object.userData.elasticHair??=createHairController(geometry);
  const result=controller.update({frame,referenceCentroids,simulationTime:frame?.time_s??0,
    recordKey,visible,active:!!frame,gravity_m_s2},positions.array);
  object.userData.hairDiagnostics={...result.diagnostics,mode:'dynamic',elastic_solver_enabled:true};
  if(result.diagnostics.updated){
    positions.needsUpdate=true;
    object.geometry.computeVertexNormals();object.geometry.computeBoundingSphere();
  }
  object.matrixAutoUpdate=false;object.matrix.identity();object.matrixWorldNeedsUpdate=true;
  return true;
}
