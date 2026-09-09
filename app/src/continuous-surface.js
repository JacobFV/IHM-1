// The same rest-space linear blend used by native contact and force scatter.
// This is geometric attachment, not a deformable skin constitutive model.
const vector = v => Array.isArray(v) && v.length === 3 && v.every(Number.isFinite);
function coordinates(reference) {
  if (!reference || !Number.isInteger(reference.length) || reference.length % 3) throw Error('Invalid surface reference coordinates');
  for (const value of reference) if (!Number.isFinite(value)) throw Error('Invalid surface reference coordinates');
}
export function prepareContinuousAsset(payload) {
  if (payload?.schema !== 'ihm.continuous-surface-binding.v1' || typeof payload.binding_identity !== 'string' || !payload.binding_identity || !Array.isArray(payload.segments) || !payload.segments.length || !Array.isArray(payload.reference_positions_m) || !payload.reference_positions_m.length || !Array.isArray(payload.weights)) throw Error('Invalid continuous surface asset');
  const segments = payload.segments.map(s => {
    if (!s || typeof s.id !== 'string' || !s.id || typeof s.bone_id !== 'string' || !s.bone_id || !vector(s.reference_centroid_m)) throw Error('Invalid continuous surface segment');
    return {...s, reference_centroid_m: [...s.reference_centroid_m]};
  });
  if (new Set(segments.map(s => s.id)).size !== segments.length) throw Error('Duplicate continuous surface segment');
  const vertexCount = payload.reference_positions_m.length, width = segments.length;
  if (payload.weights.length !== vertexCount) throw Error('Invalid continuous surface weights');
  const reference = new Float64Array(vertexCount * 3), weights = new Float64Array(vertexCount * width);
  for (let i = 0; i < vertexCount; i++) {
    const p = payload.reference_positions_m[i], row = payload.weights[i];
    if (!vector(p)) throw Error('Invalid continuous surface reference positions');
    if (!Array.isArray(row) || row.length !== width || row.some(w => !Number.isFinite(w) || w < 0) || Math.abs(row.reduce((a,b) => a+b,0)-1) > 1e-12) throw Error('Invalid continuous surface weights');
    reference.set(p, i*3); weights.set(row, i*width);
  }
  return {identity: payload.binding_identity, segments, reference, weights, vertexCount};
}
// Median partition gives a balanced tree without sorting every subtree.
function buildTree(asset) {
  const indices = Uint32Array.from({length: asset.vertexCount}, (_,i) => i), p = asset.reference;
  const compare = (a,b,axis) => p[3*a+axis]-p[3*b+axis] || a-b;
  function build(lo,hi,depth) {
    if (lo >= hi) return null;
    const axis = depth%3, mid = (lo+hi)>>>1;
    let left=lo,right=hi-1;
    while (left<right) {
      const pivot=indices[(left+right)>>>1]; let i=left,j=right;
      while(i<=j) {
        while(compare(indices[i],pivot,axis)<0)i++;
        while(compare(indices[j],pivot,axis)>0)j--;
        if(i<=j){const temp=indices[i];indices[i]=indices[j];indices[j]=temp;i++;j--;}
      }
      if(mid<=j)right=j; else if(mid>=i)left=i; else break;
    }
    return {index:indices[mid],axis,left:build(lo,mid,depth+1),right:build(mid+1,hi,depth+1)};
  }
  return build(0,indices.length,0);
}
function nearest(asset,x,y,z) {
  asset.tree ||= buildTree(asset);
  const p=asset.reference,q=[x,y,z];let best=Infinity,index=Infinity;
  function visit(node) {
    if(!node)return;
    const i=node.index,dx=p[3*i]-x,dy=p[3*i+1]-y,dz=p[3*i+2]-z,d=Math.sqrt(dx*dx+dy*dy+dz*dz);
    if(d<best || d===best && i<index){best=d;index=i;}
    const delta=q[node.axis]-p[3*i+node.axis];
    visit(delta<=0?node.left:node.right);
    // Include equal-distance candidates so ties choose the first source vertex.
    if(Math.abs(delta)<=best)visit(delta<=0?node.right:node.left);
  }
  visit(asset.tree);return index;
}
function float32Matches(a,b) {return a===b || a===Math.fround(b) || Math.abs(a-b)<=Math.max(Math.abs(b)*2**-24,2**-149);}
export function bindContinuousVertices(reference,asset,{skin=false}={}) {
  coordinates(reference);
  const vertexCount=reference.length/3,width=asset.segments.length;
  if(skin && (vertexCount!==asset.vertexCount || reference.some((v,i)=>!float32Matches(v,asset.reference[i])))) throw Error('Skin reference vertices do not match continuous surface sidecar');
  const sourceIndices=new Uint32Array(vertexCount),weights=new Float64Array(vertexCount*width);
  for(let i=0;i<vertexCount;i++) {
    const source=skin?i:nearest(asset,reference[3*i],reference[3*i+1],reference[3*i+2]);
    sourceIndices[i]=source;weights.set(asset.weights.subarray(source*width,(source+1)*width),i*width);
  }
  // Use authoritative unrounded skin rest coordinates; garments retain their rest.
  return {identity:asset.identity,segments:asset.segments,weights,sourceIndices,vertexCount,reference:skin?asset.reference.slice():Float64Array.from(reference),inputReference:Float64Array.from(reference)};
}
function assertReference(reference,bound) {
  if(reference.length!==bound.inputReference.length)throw Error('Surface reference changed after binding');
  for(let i=0;i<reference.length;i++)if(reference[i]!==bound.inputReference[i])throw Error('Surface reference changed after binding');
}
function currentTransforms(segments,transforms) {
  return segments.map(s=>{
    const t=transforms?.[s.id];
    if(!t || !vector(t.centroid_m) || !Array.isArray(t.rotation_matrix) || t.rotation_matrix.length!==3 || !t.rotation_matrix.every(vector))throw Error('Missing continuous surface segment transform');
    return t;
  });
}
export function poseContinuousVertices(reference,bound,transforms,output) {
  assertReference(reference,bound);
  if(!output || output.length!==reference.length || output===reference || output===bound.reference || output===bound.inputReference || output.buffer && (output.buffer===reference.buffer || output.buffer===bound.reference.buffer || output.buffer===bound.inputReference.buffer))throw Error('Invalid continuous surface output buffer');
  const current=currentTransforms(bound.segments,transforms),width=bound.segments.length,p=bound.reference;
  for(let i=0;i<bound.vertexCount;i++) {
    let a=0,b=0,c=0;
    for(let j=0;j<width;j++) {
      const w=bound.weights[i*width+j];if(!w)continue;
      const s=bound.segments[j],t=current[j],r=t.rotation_matrix,x=p[3*i]-s.reference_centroid_m[0],y=p[3*i+1]-s.reference_centroid_m[1],z=p[3*i+2]-s.reference_centroid_m[2];
      a+=w*(t.centroid_m[0]+r[0][0]*x+r[0][1]*y+r[0][2]*z);
      b+=w*(t.centroid_m[1]+r[1][0]*x+r[1][1]*y+r[1][2]*z);
      c+=w*(t.centroid_m[2]+r[2][0]*x+r[2][1]*y+r[2][2]*z);
    }
    output[3*i]=a;output[3*i+1]=b;output[3*i+2]=c;
  }
  return output;
}
export function continuousMaterialAnchor(bound,reference,vertexIndex) {
  assertReference(reference,bound);
  if(!Number.isInteger(vertexIndex)||vertexIndex<0||vertexIndex>=bound.vertexCount)throw Error('Invalid continuous surface vertex');
  const width=bound.segments.length;
  return {identity:bound.identity,segments:bound.segments,reference_m:Array.from(bound.reference.subarray(3*vertexIndex,3*vertexIndex+3)),weights:bound.weights.slice(width*vertexIndex,width*(vertexIndex+1))};
}
export function continuousMaterialPorts(anchor,transforms,force) {
  if(!vector(force))throw Error('Invalid continuous surface force');
  const current=currentTransforms(anchor.segments,transforms),point_m=[0,0,0],ports=[];
  for(let j=0;j<anchor.segments.length;j++) {
    const w=anchor.weights[j];if(!w)continue;
    const s=anchor.segments[j],t=current[j],offset=anchor.reference_m.map((v,k)=>v-s.reference_centroid_m[k]);
    const station=t.centroid_m.map((v,k)=>v+t.rotation_matrix[k].reduce((sum,r,l)=>sum+r*offset[l],0));
    for(let k=0;k<3;k++)point_m[k]+=w*station[k];
    if(force.some(v=>v!==0))ports.push({id:s.bone_id,point_m:station,force_n:force.map(v=>w*v)});
  }
  return {point_m,ports};
}
