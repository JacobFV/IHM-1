const finiteVector=v=>v?.length===3&&v.every(Number.isFinite);

// R maps material coordinates to the scene frame; its transpose is the inverse.
export function rotateOffset(offset,rotation,inverse=false) {
  if(!rotation)return [...offset];
  return offset.map((_,i)=>offset.reduce((sum,v,j)=>sum+v*(inverse?rotation[j][i]:rotation[i][j]),0));
}

export async function advanceScene(request,path,command) {
  try{return {frame:await request(path+'/step',command),recovered:false};}
  catch(error){
    // A lost response may follow a committed force. Read the sequence; never
    // replay an uncertain command and accidentally apply its impulse twice.
    try{return {frame:await request(path),recovered:true,error};}
    catch{throw error;}
  }
}

export function cursorSpring(point,target,stiffness=80,limit=100) {
  if(!finiteVector(point)||!finiteVector(target)||!(stiffness>0)||!(limit>0))throw Error('Invalid cursor force');
  const f=point.map((v,i)=>(target[i]-v)*stiffness),norm=Math.hypot(...f);
  return norm>limit?f.map(v=>v*limit/norm):f;
}
