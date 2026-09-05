import { deformSkinVertices } from './state.js';

// Transport each shaft as a rigid local material frame. Skin stretch changes
// follicle spacing, while this kinematic projection preserves shaft length.
export function attachedHairPositions(reference, attachment, field, target=null) {
  if(!field) {const output=target||new Float32Array(reference.length);output.set(reference);return output;}
  const triangles=attachment.reference_triangles_m,bary=attachment.barycentric;
  const count=bary?.length,stride=attachment.vertices_per_shaft;
  if(!Number.isInteger(count)||stride!==8||triangles?.length!==count*9||reference.length!==count*stride*3)throw Error('Invalid material hair attachment');
  const moved=deformSkinVertices(triangles,field, new Float64Array(triangles.length));
  const out=target||new Float32Array(reference.length);
  function frame(v,offset) {
    const a=v.slice(offset,offset+3),b=v.slice(offset+3,offset+6),c=v.slice(offset+6,offset+9);
    const e=b.map((x,i)=>x-a[i]),w=c.map((x,i)=>x-a[i]);
    const n=[e[1]*w[2]-e[2]*w[1],e[2]*w[0]-e[0]*w[2],e[0]*w[1]-e[1]*w[0]];
    const le=Math.hypot(...e),ln=Math.hypot(...n);
    if(le<=0||ln<=0)throw Error('Degenerate follicle material face');
    for(let j=0;j<3;j++){e[j]/=le;n[j]/=ln;}
    return [e,[n[1]*e[2]-n[2]*e[1],n[2]*e[0]-n[0]*e[2],n[0]*e[1]-n[1]*e[0]],n];
  }
  for(let i=0;i<count;i++) {
    const old=frame(triangles,9*i),now=frame(moved,9*i),root=[0,0,0],next=[0,0,0];
    for(let j=0;j<3;j++)for(let k=0;k<3;k++){root[k]+=bary[i][j]*triangles[i*9+j*3+k];next[k]+=bary[i][j]*moved[i*9+j*3+k];}
    for(let j=0;j<stride;j++) {
      const off=(i*stride+j)*3,r=[reference[off]-root[0],reference[off+1]-root[1],reference[off+2]-root[2]];
      const local=old.map(axis=>axis[0]*r[0]+axis[1]*r[1]+axis[2]*r[2]);
      for(let k=0;k<3;k++)out[off+k]=next[k]+now[0][k]*local[0]+now[1][k]*local[1]+now[2][k]*local[2];
    }
  }
  return out;
}
