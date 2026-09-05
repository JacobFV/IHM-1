import {createHairController} from './hair_dynamics.js';

export function attachElasticHair(object,geometry) {
  if(geometry.attachment?.kind!=='ElasticHairMaterialPoint')return false;
  object.userData.elasticHair=createHairController(geometry);
  return true;
}

// The controller returns canonical coordinates with the skin transform already
// applied. Applying the ordinary body matrix here would move every strand twice.
export function updateElasticHair(object,{frame,referenceCentroids={},recordKey,visible=true,
    gravity_m_s2=[0,-9.81,0]}={}) {
  const controller=object.userData.elasticHair;
  if(!controller)return false;
  const positions=object.geometry.getAttribute('position');
  const result=controller.update({frame,referenceCentroids,simulationTime:frame?.time_s??0,
    recordKey,visible,active:!!frame,gravity_m_s2},positions.array);
  object.userData.hairDiagnostics=result.diagnostics;
  if(result.diagnostics.updated){
    positions.needsUpdate=true;
    object.geometry.computeVertexNormals();object.geometry.computeBoundingSphere();
  }
  object.matrixAutoUpdate=false;object.matrix.identity();object.matrixWorldNeedsUpdate=true;
  return true;
}
