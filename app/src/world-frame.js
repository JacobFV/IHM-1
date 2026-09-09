// Environment geometry stays in its declared gravity frame. Interaction ports
// and body geometry use canonical coordinates; convert only at that boundary.
const identity=()=>[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
export function environmentWorldMatrix(state){
 if(!state?.object_coordinate_frame)return identity();
 if(state.object_coordinate_frame!=='gravity_aligned_world')throw Error('Unknown environment coordinate frame');
 const m=state.world_frame?.world_to_canonical;
 if(!Array.isArray(m)||m.length!==4||m.some(row=>!Array.isArray(row)||row.length!==4||row.some(v=>!Number.isFinite(v))))throw Error('Missing or invalid environment world transform');
 if(m[3].some((v,i)=>Math.abs(v-(i===3?1:0))>1e-8))throw Error('Non-affine environment world transform');
 for(let i=0;i<3;i++)for(let j=0;j<3;j++)if(Math.abs(m[i].slice(0,3).reduce((s,v,k)=>s+v*m[j][k],0)-(i===j?1:0))>1e-7)throw Error('Non-rigid environment world transform');
 const determinant=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
 if(Math.abs(determinant-1)>1e-7)throw Error('Reflected environment world transform');
 return m;
}
export function canonicalEnvironmentObject(state,object){
 if(!object)return object;
 const m=environmentWorldMatrix(state),vector=p=>m.slice(0,3).map(row=>row.slice(0,3).reduce((s,v,k)=>s+v*p[k],0)),point=p=>vector(p).map((v,i)=>v+m[i][3]);
 const out={...object};
 if(object.position_m)out.position_m=point(object.position_m);
 if(object.rotation_matrix)out.rotation_matrix=m.slice(0,3).map(row=>[0,1,2].map(j=>row.slice(0,3).reduce((s,v,k)=>s+v*object.rotation_matrix[k][j],0)));
 if(object.positions){out.positions=[];for(let i=0;i<object.positions.length;i+=3)out.positions.push(...point(object.positions.slice(i,i+3)));}
 return out;
}
