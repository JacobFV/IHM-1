// Implicit small-deflection, untwisted discrete elastic beams in SI units.
// Axial U=EA/(2h)(delta longitudinal displacement)^2; transverse
// bending U=EI/(2h^3)|second difference|^2. No radius amplification.
import {deformSkinVertices,bodyTransform} from './state.js';
const finite=(a,n,label)=>{if(a?.length!==n||!Array.from(a).every(Number.isFinite))throw Error(`Invalid ${label}`);};
export function materialHairRoots(attachment,field=null,transform=null){
 const b=attachment.barycentric,t=attachment.reference_triangles_m;finite(t,b.length*9,'root triangles');
 let moved=field?deformSkinVertices(t,field,new Float64Array(t.length)):t;
 if(transform){finite(transform,16,'root transform');const p=new Float64Array(t.length);for(let i=0;i<t.length;i+=3)for(let k=0;k<3;k++)p[i+k]=transform[k*4]*moved[i]+transform[k*4+1]*moved[i+1]+transform[k*4+2]*moved[i+2]+transform[k*4+3];moved=p;}
 const roots=new Float64Array(b.length*3),directions=new Float64Array(b.length*3);
 for(let s=0;s<b.length;s++){
  finite(b[s],3,'barycentric root');if(b[s].some(x=>x<0||x>1)||Math.abs(b[s].reduce((a,x)=>a+x,0)-1)>1e-8)throw Error('Invalid barycentric sum');
  const q=9*s,e=[0,0,0],f=[0,0,0];
  for(let k=0;k<3;k++){e[k]=moved[q+3+k]-moved[q+k];f[k]=moved[q+6+k]-moved[q+k];for(let j=0;j<3;j++)roots[3*s+k]+=b[s][j]*moved[q+3*j+k];}
  const n=[e[1]*f[2]-e[2]*f[1],e[2]*f[0]-e[0]*f[2],e[0]*f[1]-e[1]*f[0]],length=Math.hypot(...n);if(length<=1e-20)throw Error('Degenerate follicle triangle');
  for(let k=0;k<3;k++)directions[3*s+k]=n[k]/length;
 }
 return {roots_m:roots,root_tangents:directions};
}
export class ElasticHairSystem {
 constructor(data,{maxSubstep=1/240}={}){
  if(!Number.isFinite(maxSubstep)||maxSubstep<=0||maxSubstep>.02)throw Error('Invalid solver settings');
  const n=data.centerlines_m?.length/3;finite(data.centerlines_m,n*3,'centerlines');
  this.offsets=Array.from(data.strand_offsets||[]);this.radii=Float64Array.from(data.radius_m||[]);this.count=this.radii.length;
  if(!Number.isInteger(n)||!this.count||this.offsets.length!==this.count+1||this.offsets[0]!==0||this.offsets.at(-1)!==n)throw Error('Invalid strand offsets');
  const E=data.tensile_modulus_pa,B=data.bending_modulus_pa,rho=data.density_kg_m3;
  if(![E,B,rho].every(x=>Number.isFinite(x)&&x>0)||!Array.from(this.radii).every(x=>Number.isFinite(x)&&x>0))throw Error('Invalid hair material');
  this.positions=Float64Array.from(data.centerlines_m);this.velocities=new Float64Array(n*3);this.mass=new Float64Array(n);this.invMass=new Float64Array(n);this.root=new Float64Array(this.count*3);this.tangents=new Float64Array(this.count*3);this.spacing=new Float64Array(this.count);this.stretch=[];this.bend=[];this.rods=[];this.maxSubstep=maxSubstep;this.time=0;this.gravity=[0,-9.81,0];
  for(let s=0;s<this.count;s++){
   const first=this.offsets[s],end=this.offsets[s+1];if(!Number.isInteger(first)||!Number.isInteger(end)||end-first<3)throw Error('Hair needs at least three nodes');
   const r=this.radii[s],A=Math.PI*r*r,I=Math.PI*r**4/4,d=[0,1,2].map(k=>this.positions[3*(first+1)+k]-this.positions[3*first+k]),h=Math.hypot(...d);if(h<=0)throw Error('Zero hair segment');this.spacing[s]=h;this.rods.push({first,end,h,axial:E*A/h,bending:B*I/h**3,cache:new Map()});
   for(let k=0;k<3;k++){this.root[3*s+k]=this.positions[3*first+k];this.tangents[3*s+k]=d[k]/h;}
   for(let j=first;j<end;j++){
    this.mass[j]=rho*A*h*((j===first||j===end-1)?.5:1);this.invMass[j]=j<first+2?0:1/this.mass[j];
    if(j<end-1){const e=[0,1,2].map(k=>this.positions[3*(j+1)+k]-this.positions[3*j+k]);if(Math.abs(Math.hypot(...e)-h)>h*1e-6||e.some((v,k)=>Math.abs(v-d[k])>h*1e-6))throw Error('Only initially straight uniform rods are supported');this.stretch.push({a:j,b:j+1,h,k:E*A/h,lambda:0});}
    if(j>first&&j<end-1)this.bend.push({a:j-1,b:j,c:j+1,k:B*I/h**3,lambda:[0,0,0]});
   }
  }
 }
 totalMassKg(){return this.mass.reduce((a,x)=>a+x,0);}
 maxStrain(){let e=0;for(const c of this.stretch){let l=0;for(let k=0;k<3;k++)l+=(this.positions[3*c.b+k]-this.positions[3*c.a+k])**2;e=Math.max(e,Math.abs(Math.sqrt(l)/c.h-1));}return e;}
 step(dt,{gravity_m_s2=this.gravity,roots_m=this.root,root_tangents=this.tangents,external_forces_n=null}={}){
  if(!Number.isFinite(dt)||dt<0||dt>.25)throw Error('Hair timestep must be 0..0.25 s');finite(gravity_m_s2,3,'gravity');finite(roots_m,this.count*3,'roots');finite(root_tangents,this.count*3,'root tangents');if(external_forces_n)finite(external_forces_n,this.positions.length,'external forces');
  const directions=Float64Array.from(root_tangents);for(let s=0;s<this.count;s++){const l=Math.hypot(...directions.slice(3*s,3*s+3));if(l<=0)throw Error('Zero root tangent');for(let k=0;k<3;k++)directions[3*s+k]/=l;}
  if(dt===0)return;const startRoots=this.root.slice(),startDirections=this.tangents.slice(),targetRoots=Float64Array.from(roots_m);this.gravity=Array.from(gravity_m_s2);
  const steps=Math.ceil(dt/this.maxSubstep-1e-10),h=dt/steps,x=this.positions,w=this.invMass;
  for(let sub=1;sub<=steps;sub++){
   const before=x.slice();
   for(let s=0;s<this.count;s++){
    const rod=this.rods[s],{first,end}=rod,f=sub/steps;
    const dir=[0,1,2].map(k=>startDirections[3*s+k]*(1-f)+directions[3*s+k]*f),length=Math.hypot(...dir);
    if(length<1e-12)for(let k=0;k<3;k++)dir[k]=directions[3*s+k];else for(let k=0;k<3;k++)dir[k]/=length;
    const root=[startRoots[3*s]*(1-f)+targetRoots[3*s]*f,startRoots[3*s+1]*(1-f)+targetRoots[3*s+1]*f,startRoots[3*s+2]*(1-f)+targetRoots[3*s+2]*f],N=end-first,F=N-2;
    let factors=rod.cache.get(h);
    if(!factors){
     factors=[];
     for(let axis=0;axis<2;axis++){
      const A=Array.from({length:F},()=>new Float64Array(F));for(let i=0;i<F;i++)A[i][i]=this.mass[first+i+2]/(h*h);
      const stencil=axis===0?[-1,1]:[1,-2,1],stiffness=axis===0?rod.axial:rod.bending;
      for(let a=0;a<=N-stencil.length;a++)for(let i=0;i<stencil.length;i++)for(let j=0;j<stencil.length;j++)if(a+i>=2&&a+j>=2)A[a+i-2][a+j-2]+=stiffness*stencil[i]*stencil[j];
      const L=Array.from({length:F},()=>new Float64Array(F));for(let i=0;i<F;i++)for(let j=0;j<=i;j++){let a=A[i][j];for(let k=0;k<j;k++)a-=L[i][k]*L[j][k];L[i][j]=i===j?Math.sqrt(a):a/L[j][j];if(!Number.isFinite(L[i][j]))throw Error('Invalid hair factorization');}const R=Array.from({length:F},()=>new Float64Array(F));for(let col=0;col<F;col++){const q=new Float64Array(F);q[col]=this.mass[first+col+2]/(h*h);for(let i=0;i<F;i++){for(let j=0;j<i;j++)q[i]-=L[i][j]*q[j];q[i]/=L[i][i];}for(let i=F-1;i>=0;i--){for(let j=i+1;j<F;j++)q[i]-=L[j][i]*q[j];q[i]/=L[i][i];}for(let i=0;i<F;i++)R[i][col]=q[i];}factors.push(R);
     }
     rod.cache.set(h,factors);if(rod.cache.size>8)rod.cache.delete(rod.cache.keys().next().value);
    }
    const predicted=rod.predicted||(rod.predicted=new Float64Array(F*3)),axial=rod.axialPredicted||(rod.axialPredicted=new Float64Array(F));
    for(let j=0;j<F;j++){const i=first+j+2;let along=0;for(let k=0;k<3;k++){const d=before[3*i+k]+h*this.velocities[3*i+k]+h*h*(this.gravity[k]+(external_forces_n?.[3*i+k]||0)*w[i])-root[k]-dir[k]*(j+2)*rod.h;predicted[3*j+k]=d;along+=dir[k]*d;}axial[j]=along;}
    for(let j=0;j<N;j++)for(let k=0;k<3;k++)x[3*(first+j)+k]=root[k]+dir[k]*j*rod.h;
    for(let i=0;i<F;i++){let dx=0,dy=0,dz=0,along=0;for(let j=0;j<F;j++){const bend=factors[1][i][j];dx+=bend*predicted[3*j];dy+=bend*predicted[3*j+1];dz+=bend*predicted[3*j+2];along+=(factors[0][i][j]-bend)*axial[j];}const q=3*(first+i+2);x[q]+=dx+dir[0]*along;x[q+1]+=dy+dir[1]*along;x[q+2]+=dz+dir[2]*along;}
   }
   for(let q=0;q<x.length;q++)this.velocities[q]=(x[q]-before[q])/h;
  }
  this.root.set(targetRoots);this.tangents.set(directions);this.time+=dt;
 }
 energy(){let U=0;const x=this.positions;for(let j=0;j<this.mass.length;j++)for(let k=0;k<3;k++)U+=this.mass[j]*(.5*this.velocities[3*j+k]**2-this.gravity[k]*x[3*j+k]);
  for(let s=0;s<this.count;s++){const rod=this.rods[s],d=this.tangents.slice(3*s,3*s+3);for(let j=rod.first;j<rod.end-1;j++){let axial=0;for(let k=0;k<3;k++)axial+=d[k]*(x[3*(j+1)+k]-x[3*j+k]);U+=.5*rod.axial*(axial-rod.h)**2;}for(let j=rod.first+1;j<rod.end-1;j++){const c=[0,1,2].map(k=>x[3*(j-1)+k]-2*x[3*j+k]+x[3*(j+1)+k]),parallel=c.reduce((sum,z,k)=>sum+z*d[k],0);U+=.5*rod.bending*(c.reduce((sum,z)=>sum+z*z,0)-parallel*parallel);}}
  return U;
 }
 maxSlope(){let maximum=0;for(let s=0;s<this.count;s++){const r=this.rods[s],d=this.tangents.slice(3*s,3*s+3);for(let i=r.first;i<r.end-1;i++){let parallel=0,length=0;for(let k=0;k<3;k++){const e=this.positions[3*(i+1)+k]-this.positions[3*i+k];parallel+=d[k]*e;length+=e*e;}maximum=Math.max(maximum,Math.sqrt(Math.max(0,length-parallel*parallel))/r.h);}}return maximum;}
 resetPose(roots,tangents){finite(roots,this.root.length,'reset roots');finite(tangents,this.tangents.length,'reset directions');for(let s=0;s<this.count;s++){const d=tangents.slice(3*s,3*s+3),l=Math.hypot(...d);if(l<=0)throw Error('Zero root tangent');for(let k=0;k<3;k++){this.root[3*s+k]=roots[3*s+k];this.tangents[3*s+k]=d[k]/l;for(let i=this.offsets[s];i<this.offsets[s+1];i++)this.positions[3*i+k]=roots[3*s+k]+d[k]/l*(i-this.offsets[s])*this.spacing[s];}}this.velocities.fill(0);}
 checkpoint(){return {positions:Array.from(this.positions),velocities:Array.from(this.velocities),roots:Array.from(this.root),tangents:Array.from(this.tangents),gravity:this.gravity.slice(),time:this.time};}
 restore(s){finite(s.positions,this.positions.length,'checkpoint positions');finite(s.velocities,this.velocities.length,'checkpoint velocities');finite(s.roots,this.root.length,'checkpoint roots');finite(s.tangents,this.tangents.length,'checkpoint tangents');finite(s.gravity,3,'checkpoint gravity');if(!Number.isFinite(s.time)||s.time<0)throw Error('Invalid checkpoint clock');this.positions.set(s.positions);this.velocities.set(s.velocities);this.root.set(s.roots);this.tangents.set(s.tangents);this.gravity=s.gravity.slice();this.time=s.time;}
 writeTubePositions(target=null,sides=4){if(!Number.isInteger(sides)||sides<3)throw Error('Invalid tube sides');const out=target||new Float32Array(this.positions.length*sides);if(out.length!==this.positions.length*sides)throw Error('Wrong tube buffer');const x=this.positions;
  for(let s=0;s<this.count;s++)for(let i=this.offsets[s];i<this.offsets[s+1];i++){
   const a=3*Math.max(i-1,this.offsets[s]),b=3*Math.min(i+1,this.offsets[s+1]-1);let nx=x[b]-x[a],ny=x[b+1]-x[a+1],nz=x[b+2]-x[a+2];const norm=Math.hypot(nx,ny,nz);nx/=norm;ny/=norm;nz/=norm;
   let ux,uy,uz;if(Math.abs(ny)>.9){ux=0;uy=nz;uz=-ny;}else{ux=-nz;uy=0;uz=nx;}const un=Math.hypot(ux,uy,uz);ux/=un;uy/=un;uz/=un;
   const vx=ny*uz-nz*uy,vy=nz*ux-nx*uz,vz=nx*uy-ny*ux,r=this.radii[s],q=3*i;
   for(let j=0;j<sides;j++){const c=sides===4?(j===0?1:j===2?-1:0):Math.cos(2*Math.PI*j/sides),sn=sides===4?(j===1?1:j===3?-1:0):Math.sin(2*Math.PI*j/sides),o=(i*sides+j)*3;out[o]=x[q]+r*(c*ux+sn*vx);out[o+1]=x[q+1]+r*(c*uy+sn*vy);out[o+2]=x[q+2]+r*(c*uz+sn*vz);}
  }return out;
 }
}

// All output coordinates include the skin entity body transform: use identity
// object matrix, otherwise body motion would be applied twice.
export function createHairController(geometry,{maxCatchupS=.25,maxUpdateHz=10,...solverOptions}={}){
 if(!Number.isFinite(maxCatchupS)||maxCatchupS<=0||maxCatchupS>.25||!Number.isFinite(maxUpdateHz)||maxUpdateHz<=0||maxUpdateHz>20)throw Error('Invalid hair catchup/update budget');
 const system=new ElasticHairSystem(geometry.strands,solverOptions),attachment=geometry.attachment;
 if(system.count>512)throw Error('Interactive hair guide budget exceeded');
 const render=geometry.render_strands?new ElasticHairSystem(geometry.render_strands,solverOptions):system;
 if(render.count>4096)throw Error('Interactive render strand budget exceeded');
 const map=geometry.guide_interpolation;
 if(render!==system){if(map?.length!==render.count)throw Error('Missing guide interpolation');for(const row of map){if(!row.length||row.some(([i,w])=>!Number.isInteger(i)||i<0||i>=system.count||!Number.isFinite(w)||w<0)||Math.abs(row.reduce((sum,[,w])=>sum+w,0)-1)>1e-8)throw Error('Invalid guide interpolation');}}
 let lastTime=null,lastRecord=null,paused=false,cached=null,diagnostics=null;
 const response=(target,extra={})=>{if(target&&cached){if(target.length!==cached.length)throw Error('Wrong tube buffer');target.set(cached);}return {positions:target||cached,diagnostics:{...diagnostics,...extra}};};
 return {system,update({frame=null,referenceCentroids={},simulationTime=frame?.time_s,recordKey=null,gravity_m_s2=[0,-9.81,0],external_forces_n=null,active=true,visible=true}={},target=null){
  if(!Number.isFinite(simulationTime)||simulationTime<0)throw Error('Hair requires an explicit finite mechanical clock');finite(gravity_m_s2,3,'gravity');if(external_forces_n)finite(external_forces_n,system.positions.length,'external forces');
  if(!active||!visible){paused=true;return response(target,{paused:true,updated:false});}
  const dt=lastTime===null?0:simulationTime-lastTime;
  const resetReason=lastTime===null?'initial':recordKey!==lastRecord?'record_changed':dt<0?'rewind':paused?'resume_without_catchup':dt>maxCatchupS?'time_gap_not_simulated':null;
  if(!resetReason&&dt<1/maxUpdateHz-1e-10)return response(target,{paused:false,updated:false});
  const field=frame?.respiration?.skin_field,id=attachment.skin_entity_id,skinField=field?.entity_ids?.includes(id)?field:null,transform=bodyTransform(frame?.entities?.[id],referenceCentroids?.[id]);
  const roots=materialHairRoots(attachment,skinField,transform);
  if(resetReason){system.resetPose(roots.roots_m,roots.root_tangents);system.time=simulationTime;system.gravity=Array.from(gravity_m_s2);}
  else system.step(dt,{...roots,gravity_m_s2,external_forces_n});
  if(render!==system){
   const r=materialHairRoots(geometry.render_attachment,skinField,transform);render.resetPose(r.roots_m,r.root_tangents);
   for(let s=0;s<render.count;s++){
    const a=render.offsets[s],N=render.offsets[s+1]-a;
    for(let j=2;j<N;j++)for(const [g,w] of map[s]){
     const first=system.offsets[g],G=system.offsets[g+1]-first,q=j/(N-1)*(G-1),lo=Math.floor(q),hi=Math.min(G-1,lo+1),f=q-lo;
     for(let k=0;k<3;k++){
      const displacement=(1-f)*system.positions[3*(first+lo)+k]+f*system.positions[3*(first+hi)+k]-system.root[3*g+k]-system.tangents[3*g+k]*q*system.spacing[g];
      render.positions[3*(a+j)+k]+=w*displacement*(render.spacing[s]*(N-1)/(system.spacing[g]*(G-1)));
     }
    }
   }
  }
  lastTime=simulationTime;lastRecord=recordKey;paused=false;
  cached=render.writeTubePositions(cached);const slope=system.maxSlope(),renderSlope=render.maxSlope();diagnostics={time_s:simulationTime,reset_reason:resetReason,paused:false,updated:true,maximum_transverse_slope:slope,render_maximum_transverse_slope:renderSlope,within_small_deflection:Math.max(slope,renderSlope)<=.3,physical_radius_multiplier:1,collision_model:'none',mechanical_energy_j:system.energy(),simulated_guides:system.count,rendered_fibers:render.count,render_motion:render===system?'direct':'three_nearest_guide_displacement',maximum_update_hz:maxUpdateHz};return response(target);
 }};
}
