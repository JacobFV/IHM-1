import * as THREE from 'three';
import { environmentMaterial, materialUVs } from './environment-materials.js';
import {deformSkinVertices,bodyTransform,geometryArrays} from './state.js';

// Garments are synthesized around the retained skin's transverse sections.
// This is a geometric pattern/ease prior, not a cloth contact solve.
function hull(points) {
  const p=points.map(v=>[v[0],v[1]]).sort((a,b)=>a[0]-b[0]||a[1]-b[1]);
  const cross=(a,b,c)=>(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
  const half=q=>{const h=[];for(const v of q){while(h.length>1&&cross(h.at(-2),h.at(-1),v)<=0)h.pop();h.push(v);}return h;};
  return [...half(p).slice(0,-1),...half([...p].reverse()).slice(0,-1)];
}

function sections(positions,indices,y) {
  const nodes=new Map(),key=v=>v.map(x=>Math.round(x*1e7)).join(',');
  for(let f=0;f<indices.length;f+=3) {
    const hit=[];
    for(let e=0;e<3;e++) {
      const a=3*indices[f+e],b=3*indices[f+(e+1)%3],ya=positions[a+1],yb=positions[b+1];
      if((ya<y)===(yb<y))continue;
      const t=(y-ya)/(yb-ya);hit.push([positions[a]+t*(positions[b]-positions[a]),positions[a+2]+t*(positions[b+2]-positions[a+2])]);
    }
    if(hit.length!==2)continue;
    const keys=hit.map(key);if(keys[0]===keys[1])continue;
    for(let j=0;j<2;j++){if(!nodes.has(keys[j]))nodes.set(keys[j],{point:hit[j],neighbors:[]});nodes.get(keys[j]).neighbors.push(keys[1-j]);}
  }
  const seen=new Set(),groups=[];
  for(const id of nodes.keys()) {
    if(seen.has(id))continue;
    const q=[id],points=[];
    while(q.length){const i=q.pop();if(seen.has(i))continue;seen.add(i);const node=nodes.get(i);points.push(node.point);q.push(...node.neighbors);}
    if(points.length>=5)groups.push(points);
  }
  return groups;
}

function envelope(groups,kind,scale,y) {
  const center=p=>p.reduce((s,v)=>s+v[0],0)/p.length;
  let chosen;
  if(kind==='pelvis')chosen=groups.filter(p=>Math.abs(center(p))<.23*scale).flat();
  else {
    const target=(kind==='left'?.085:kind==='right'?-.085:0)*scale;
    chosen=[...groups].sort((a,b)=>Math.abs(center(a)-target)-Math.abs(center(b)-target))[0];
  }
  if(!chosen?.length)throw Error('Skin has no usable section for the garment pattern');
  // BP skin joins the arms to the trunk above the axilla and contains hand
  // sections close to the hip. Regional pattern masks prevent those unrelated
  // surfaces from enlarging a torso/hip contour. Bounds are stated priors in
  // this generic reference frame, not recovered anatomical boundaries.
  const halfWidth=(kind==='torso'&&y>.33&&y<.57?.175:.19)*scale;
  chosen=chosen.filter(v=>Math.abs(v[0])<=halfWidth);
  if(chosen.length<3)throw Error('Garment regional mask has insufficient skin support');
  const outline=hull(chosen);
  const xs=outline.map(p=>p[0]),zs=outline.map(p=>p[1]);
  const origin=[(Math.min(...xs)+Math.max(...xs))/2,(Math.min(...zs)+Math.max(...zs))/2];
  return {outline,origin};
}

function radialPoint({outline,origin},angle,ease) {
  const d=[Math.cos(angle),Math.sin(angle)];let distance=Infinity;
  // Every convex half-plane limits the ray. This preserves the measured section
  // envelope without copying every skin wrinkle into the garment.
  for(let i=0;i<outline.length;i++) {
    const a=outline[i],b=outline[(i+1)%outline.length],edge=[b[0]-a[0],b[1]-a[1]];
    const cross=(u,v)=>u[0]*v[1]-u[1]*v[0];
    const denom=cross(edge,d);
    if(denom< -1e-14)distance=Math.min(distance,cross(edge,[a[0]-origin[0],a[1]-origin[1]])/denom);
  }
  if(!Number.isFinite(distance)||distance<=0)throw Error('Invalid garment section envelope');
  return [origin[0]+(distance+ease)*d[0],origin[1]+(distance+ease)*d[1]];
}

export function buildGarments(skin,source={}) {
  const p=Array.from(skin.positions.flat?.()||skin.positions),ind=Array.from(skin.indices.flat?.()||skin.indices);
  if(p.length<9||p.length%3||!p.every(Number.isFinite)||ind.length%3||!ind.every(i=>Number.isInteger(i)&&i>=0&&3*i<p.length))
    throw Error('Clothing needs finite indexed canonical skin geometry');
  let ymin=Infinity,ymax=-Infinity;for(let i=1;i<p.length;i+=3){ymin=Math.min(ymin,p[i]);ymax=Math.max(ymax,p[i]);}
  const scale=(ymax-ymin)/1.7194712,center=(ymin+ymax)/2+.005135*scale;
  const yAt=y=>center+scale*y,N=64;
  const section=(y,kind)=>envelope(sections(p,ind,yAt(y)),kind,scale,y);
  const make=(id,name)=>({id,name,positions:[],indices:[],source:{id:source.id,geometry_sha256:source.geometry_sha256??null},
    physical_contact_solved:false,evidence_kind:'engineered_garment_pattern',
    assumptions:['Canonical skin transverse convex envelope','Pattern boundaries and ease are engineering priors','Torso/hip transverse regional mask: ±190 mm, reduced to ±175 mm above the axilla; scaled with reference height','Kinematic thoracic following; no friction or genital-tissue solve in this display']});
  const add=(g,v)=>{g.positions.push(...v);return g.positions.length/3-1;};
  const join=(g,a,b,keep=()=>true)=>{for(let j=0;j<N;j++){const k=(j+1)%N;if(keep(j)){g.indices.push(a[j],a[k],b[j],a[k],b[k],b[j]);}}};
  const ring=(g,y,kind,ease,angle0=-Math.PI/2)=>{
    const env=section(y,kind);return Array.from({length:N},(_,j)=>{const [x,z]=radialPoint(env,angle0+j*2*Math.PI/N,ease*scale);return add(g,[x,yAt(y),z]);});
  };
  const shirt=make('shirt','Sleeveless shirt');
  let prior;
  for(let r=0;r<=36;r++){
    const y=.035+r*(.605-.035)/36,current=ring(shirt,y,'torso',.012);
    if(prior)join(shirt,prior,current,j=>{
      const angle=-Math.PI/2+(j+.5)*2*Math.PI/N;
      return !(y>.35&&y<.59&&Math.abs(Math.cos(angle))>.84);
    });
    prior=current;
  }
  const shorts=make('shorts','Shorts');prior=null;
  for(let r=0;r<=14;r++) {
    const current=ring(shorts,.105-r*.235/14,'pelvis',.018);
    if(prior)join(shorts,prior,current);prior=current;
  }
  // Shared crotch seam: both leg tubes use the same vertex indices. There are
  // exactly three boundary loops (waist, left hem, right hem), no hidden caps.
  const outer=prior,front=shorts.positions.slice(3*outer[N/2],3*outer[N/2]+3),back=shorts.positions.slice(3*outer[0],3*outer[0]+3);
  const seam=[];
  for(let j=1;j<N/2;j++){const t=j/(N/2);seam.push(add(shorts,[front[0]*(1-t)+back[0]*t,yAt(-.13-.035*Math.sin(Math.PI*t)),front[2]*(1-t)+back[2]*t]));}
  for(const [kind,start] of [['left',0],['right',N/2]]) {
    const top=Array.from({length:N/2+1},(_,j)=>outer[(start+j)%N]).concat(kind==='left'?seam:[...seam].reverse());
    prior=top;
    const topPositions=top.map(i=>shorts.positions.slice(3*i,3*i+3));
    for(let r=1;r<=12;r++){
      const t=r/12,y=-.13-.20*t,env=section(y,kind),current=[];
      for(let j=0;j<N;j++){
        const angle=-Math.PI/2+(start+j)*2*Math.PI/N,[x,z]=radialPoint(env,angle,.014*scale);
        // Smooth the branch into the measured leg envelope across the first
        // 6 cm below the crotch; the residual is retained as pattern synthesis.
        const blend=Math.min(1,t/.3),reference=topPositions[j];
        current.push(add(shorts,[reference[0]*(1-blend)+x*blend,reference[1]*(1-t)+yAt(-.33)*t,reference[2]*(1-blend)+z*blend]));
      }
      join(shorts,prior,current);prior=current;
    }
  }
  // Arm openings remove faces. Compact their orphaned vertices so every node
  // in a downstream mechanical materialization has positive lumped mass.
  for(const garment of [shirt,shorts]) {
    const used=[...new Set(garment.indices)].sort((a,b)=>a-b),lookup=new Map(used.map((old,i)=>[old,i]));
    garment.positions=used.flatMap(i=>garment.positions.slice(3*i,3*i+3));
    garment.indices=garment.indices.map(i=>lookup.get(i));
  }
  return [shirt,shorts];
}

export class ClothingView {
  constructor(parent) {
    this.group=new THREE.Group();this.group.name='clothing';parent.add(this.group);
    this.meshes=new Map();this.enabled=new Map([['shirt',true],['shorts',true]]);
  }
  fit(skin,source) {
    this.clear();this.source=source;
    this.materializations=buildGarments(skin,source);
    for(const garment of this.materializations) {
      const geometry=new THREE.BufferGeometry();geometry.setAttribute('position',new THREE.Float32BufferAttribute(garment.positions,3));
      geometry.setIndex(garment.indices);geometry.computeVertexNormals();
      // Cloth, with the weave the fabric family already carries. A garment is an
      // open shell, so it keeps both faces.
      const material=environmentMaterial('fabric',garment.id==='shirt'?'#759b9a':'#273e54',{doubleSided:true});
      materialUVs(geometry,'fabric');
      const mesh=new THREE.Mesh(geometry,material);mesh.name=garment.id;mesh.visible=this.enabled.get(garment.id);
      mesh.userData.reference=geometry.attributes.position.array.slice();mesh.userData.garment=garment;
      this.meshes.set(garment.id,mesh);this.group.add(mesh);
    }
  }
  setEnabled(id,enabled) {
    if(!this.enabled.has(id)||typeof enabled!=='boolean')throw Error('Unknown garment or invalid enabled state');
    this.enabled.set(id,enabled);if(this.meshes.has(id))this.meshes.get(id).visible=enabled;
  }
  update(frame,centroids={}) {
    const field=frame?.respiration?.skin_field;
    for(const [id,mesh] of this.meshes) {
      const positions=mesh.geometry.attributes.position;
      deformSkinVertices(mesh.userData.reference,id==='shirt'&&field?.entity_ids?.includes(this.source?.id)?field:null,positions.array);
      positions.needsUpdate=true;mesh.geometry.computeVertexNormals();mesh.geometry.computeBoundingSphere();
      mesh.matrixAutoUpdate=false;mesh.matrix.set(...bodyTransform(frame?.entities?.[this.source?.id],centroids[this.source?.id]));mesh.matrixWorldNeedsUpdate=true;
    }
  }
  setClipping(planes){for(const mesh of this.meshes.values())mesh.material.clippingPlanes=planes;}
  clear(){for(const mesh of this.meshes.values()){mesh.geometry.dispose();mesh.material.dispose();mesh.removeFromParent();}this.meshes.clear();}
  dispose(){this.clear();this.group.removeFromParent();}
}

// ---------------------------------------------------------------------------
// The wardrobe: 33 CC0 MakeHuman garments registered onto this body's outer
// envelope and cloth-simulated against it, served by /api/clothing. Nothing is
// synthesized here — every vertex comes from the run, in the same
// bodyparts3d-display-m frame as the skin, so a garment needs no fitting step
// and is placed by the body's own transform alone.
export class WardrobeView {
  constructor(parent,{fetchGeometry}) {
    this.group=new THREE.Group();this.group.name='wardrobe';parent.add(this.group);
    this.fetchGeometry=fetchGeometry;this.meshes=new Map();this.catalog=new Map();
    this.drawOrder=[];this.active=new Set();this.pending=new Map();this.generation=0;
    this.source=null;this.clipping=null;
  }
  setCatalog(catalog) {
    this.catalog=new Map((catalog.garments||[]).map(g=>[g.id,g]));
    this.drawOrder=catalog.draw_order||[];
    for(const [id,mesh] of this.meshes)if(!this.catalog.has(id)){mesh.geometry.dispose();mesh.material.dispose();mesh.removeFromParent();this.meshes.delete(id);}
  }
  // The garment follows the skin entity it was registered against, so the body
  // transform that moves the skin moves the clothes with it.
  bind(source){this.source=source;}
  order(garment) {
    const index=this.drawOrder.indexOf((garment.slots||[])[0]);
    return index<0?this.drawOrder.length:index;
  }
  async build(id) {
    if(this.meshes.has(id))return this.meshes.get(id);
    if(this.pending.has(id))return this.pending.get(id);
    const garment=this.catalog.get(id);
    if(!garment?.geometry_url)throw Error('No served geometry for garment '+id);
    const task=(async()=>{
      const data=await this.fetchGeometry(garment.geometry_url);
      const {positions,indices}=geometryArrays(data);
      if(!indices)throw Error('Garment geometry must be indexed');
      const geometry=new THREE.BufferGeometry();
      geometry.setAttribute('position',new THREE.BufferAttribute(positions,3));
      geometry.setIndex(new THREE.BufferAttribute(indices,1));
      geometry.computeVertexNormals();
      materialUVs(geometry,'fabric');
      const mesh=new THREE.Mesh(geometry,environmentMaterial('fabric',garment.colour||'#9aa7ad',{doubleSided:true}));
      mesh.name=id;mesh.renderOrder=this.order(garment);
      mesh.userData.reference=positions.slice();mesh.userData.garment=garment;
      if(this.clipping)mesh.material.clippingPlanes=this.clipping;
      mesh.visible=this.active.has(id);
      this.meshes.set(id,mesh);this.group.add(mesh);
      return mesh;
    })().finally(()=>this.pending.delete(id));
    this.pending.set(id,task);
    return task;
  }
  // Selection is the catalog's, not this view's: whatever the tiles report is
  // worn is loaded and shown, and everything else is hidden rather than removed
  // so re-selecting a garment costs no request.
  async setActive(ids) {
    const current=++this.generation;
    this.active=new Set(ids.filter(id=>this.catalog.has(id)));
    for(const [id,mesh] of this.meshes)mesh.visible=this.active.has(id);
    const failures=[];
    await Promise.all([...this.active].map(async id=>{
      try{const mesh=await this.build(id);if(current===this.generation)mesh.visible=this.active.has(id);}
      catch(error){failures.push(id+': '+error.message);}
    }));
    return failures;
  }
  update(frame,centroids={}) {
    const field=frame?.respiration?.skin_field;
    const applied=field?.entity_ids?.includes(this.source?.id)?field:null;
    const transform=bodyTransform(frame?.entities?.[this.source?.id],centroids[this.source?.id]);
    for(const mesh of this.meshes.values()) {
      if(!mesh.visible)continue;
      const positions=mesh.geometry.attributes.position;
      deformSkinVertices(mesh.userData.reference,applied,positions.array);
      positions.needsUpdate=true;mesh.geometry.computeVertexNormals();mesh.geometry.computeBoundingSphere();
      mesh.matrixAutoUpdate=false;mesh.matrix.set(...transform);mesh.matrixWorldNeedsUpdate=true;
    }
  }
  setClipping(planes){this.clipping=planes;for(const mesh of this.meshes.values())mesh.material.clippingPlanes=planes;}
  clear(){for(const mesh of this.meshes.values()){mesh.geometry.dispose();mesh.material.dispose();mesh.removeFromParent();}this.meshes.clear();this.pending.clear();}
  dispose(){this.generation++;this.clear();this.group.removeFromParent();}
}
