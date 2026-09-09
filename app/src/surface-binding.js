import {surfaceAssets} from './surface-assets.js';
import {bindContinuousVertices,poseContinuousVertices,continuousMaterialAnchor} from './continuous-surface.js';
// Same hard rest-envelope support as the mechanics contact skin. This is an
// inferred segment attachment, with discontinuities at joints; no FEM or cloth.
export function surfaceMetadata(frame){
 const binding=frame?.mechanics?.surface_binding||frame?.surface_binding;
 if(!binding)return null;
 const transforms=frame?.mechanics?.surface_transforms||frame?.surface_transforms;
 const continuous=binding.schema==='ihm.continuous-surface-binding.v1'&&binding.rule==='graph_regularized_linear_blend';
 if(!continuous&&(binding.schema!=='ihm.segment-surface-binding.v1'||binding.rule!=='nearest_named_bone_envelope'||binding.tie_break!=='lexicographic_native_segment_id')||binding.reference_coordinate_frame!=='canonical_rest'||binding.coordinate_frame!=='canonical_current_world'||!binding.binding_identity||!Array.isArray(binding.segments)||!binding.segments.length||!transforms)throw Error('Invalid segment surface binding');
 return {binding,transforms,continuous};
}
function vector(v){return Array.isArray(v)&&v.length===3&&v.every(Number.isFinite);}
export function bindSurfaceVertices(reference,binding){
 if(reference.length%3||!Array.from(reference).every(Number.isFinite))throw Error('Invalid surface reference coordinates');
 const segments=[...binding.segments].sort((a,b)=>a.id<b.id?-1:a.id>b.id?1:0);
 if(new Set(segments.map(s=>s.id)).size!==segments.length||segments.some(s=>typeof s.id!=='string'||typeof s.bone_id!=='string'||![s.bounds_min_m,s.bounds_max_m,s.reference_centroid_m].every(vector)||s.bounds_min_m.some((v,i)=>v>s.bounds_max_m[i])))throw Error('Invalid surface segment envelopes');
 const owners=new Uint16Array(reference.length/3);
 for(let i=0;i<owners.length;i++){
  let best=Infinity;
  for(let j=0;j<segments.length;j++){
   const s=segments[j];let distance=0;
   for(let k=0;k<3;k++){const p=reference[3*i+k],d=Math.max(s.bounds_min_m[k]-p,p-s.bounds_max_m[k],0);distance+=d*d;}
   if(distance<best){best=distance;owners[i]=j;}
  }
 }
 return {identity:binding.binding_identity,segments,owners};
}
export function poseSurfaceVertices(reference,bound,transforms,output){
 if(reference===output||reference.length!==output.length||bound.owners.length*3!==reference.length)throw Error('Invalid surface output buffer');
 const current=bound.segments.map(s=>{const t=transforms[s.id];if(!t||!vector(t.centroid_m)||!Array.isArray(t.rotation_matrix)||t.rotation_matrix.length!==3||!t.rotation_matrix.every(vector))throw Error('Missing surface segment transform');return t;});
 for(let i=0;i<bound.owners.length;i++){
  const j=bound.owners[i],s=bound.segments[j],t=current[j],x=reference[3*i]-s.reference_centroid_m[0],y=reference[3*i+1]-s.reference_centroid_m[1],z=reference[3*i+2]-s.reference_centroid_m[2];
  for(let k=0;k<3;k++)output[3*i+k]=t.centroid_m[k]+t.rotation_matrix[k][0]*x+t.rotation_matrix[k][1]*y+t.rotation_matrix[k][2]*z;
 }
 return output;
}
export function updateSegmentSurface(mesh,frame,entityId,{garment=false}={}){
 const metadata=surfaceMetadata(frame),data=mesh.userData,positions=mesh.geometry?.attributes?.position;
 if(!metadata||!positions||!garment&&!metadata.binding.surface_entity_ids?.includes(entityId)){
  if(data.surfaceBindingPending){mesh.visible=data.surfaceBindingVisible;delete data.surfaceBindingPending;}
  if(data.segmentSurfaceActive&&positions){positions.array.set(data.segmentSurfaceReference);positions.needsUpdate=true;mesh.geometry.computeVertexNormals();mesh.geometry.computeBoundingSphere();data.segmentSurfaceActive=false;}
  return false;
 }
 data.segmentSurfaceReference ||= (data.reference||data.skinReference||positions.array).slice();
 if(metadata.continuous){
  const asset=surfaceAssets.get(metadata.binding);
  if(!asset){if(!data.surfaceBindingPending)data.surfaceBindingVisible=mesh.visible;data.surfaceBindingPending=true;mesh.visible=false;return true;}
  if(data.surfaceBindingPending){mesh.visible=data.surfaceBindingVisible;delete data.surfaceBindingPending;}
  if(data.segmentSurfaceBinding?.identity!==metadata.binding.binding_identity)data.segmentSurfaceBinding={...bindContinuousVertices(data.segmentSurfaceReference,asset,{skin:!garment}),continuous:true};
  poseContinuousVertices(data.segmentSurfaceReference,data.segmentSurfaceBinding,metadata.transforms,positions.array);
 }else{
  if(data.segmentSurfaceBinding?.identity!==metadata.binding.binding_identity)data.segmentSurfaceBinding=bindSurfaceVertices(data.segmentSurfaceReference,metadata.binding);
  poseSurfaceVertices(data.segmentSurfaceReference,data.segmentSurfaceBinding,metadata.transforms,positions.array);
 }
 positions.needsUpdate=true;mesh.geometry.computeVertexNormals();mesh.geometry.computeBoundingSphere();
 mesh.matrixAutoUpdate=false;mesh.matrix.identity();mesh.matrixWorldNeedsUpdate=true;data.segmentSurfaceActive=true;
 return true;
}
// A triangle can straddle hard owners. Pick the nearest displayed vertex's
// support, retain that bone material frame, and disclose the approximation.
export function surfacePickBone(mesh,face,point){
 const bound=mesh.userData.segmentSurfaceActive&&mesh.userData.segmentSurfaceBinding;
 if(!bound||bound.continuous||!face)return null;
 const p=mesh.geometry.attributes.position.array;let best=Infinity,owner;
 for(const i of [face.a,face.b,face.c]){const d=[0,1,2].reduce((sum,k)=>sum+(p[3*i+k]-point[k])**2,0);if(d<best){best=d;owner=bound.segments[bound.owners[i]].bone_id;}}
 return owner;
}

export function surfacePickAnchor(mesh,face,point){
 const bound=mesh.userData.segmentSurfaceActive&&mesh.userData.segmentSurfaceBinding;
 if(!bound?.continuous||!face)return null;
 const p=mesh.geometry.attributes.position.array;let best=Infinity,index;
 for(const i of [face.a,face.b,face.c]){const d=[0,1,2].reduce((sum,k)=>sum+(p[3*i+k]-point[k])**2,0);if(d<best){best=d;index=i;}}
 return continuousMaterialAnchor(bound,mesh.userData.segmentSurfaceReference,index);
}
