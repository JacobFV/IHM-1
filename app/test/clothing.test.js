import {test} from 'node:test';
import assert from 'node:assert/strict';
import * as THREE from 'three';
import fs from 'node:fs';
import zlib from 'node:zlib';
import {buildGarments, ClothingView} from '../src/clothing.js';

function bodyGeometry() {
  const geometry=new THREE.CylinderGeometry(.16,.16,1.73,64,80,true);
  return {positions:Array.from(geometry.attributes.position.array),indices:Array.from(geometry.index.array)};
}
test('garments have finite, nondegenerate triangles, three shorts boundaries, and source provenance',()=>{
  const garments=buildGarments(bodyGeometry(),{id:'test-skin',geometry_sha256:'test-hash'});
  assert.deepEqual(garments.map(g=>g.id),['shirt','shorts']);
  for(const garment of garments) {
    assert.equal(garment.source.geometry_sha256,'test-hash');
    assert.equal(garment.physical_contact_solved,false);
    assert.ok(garment.positions.every(Number.isFinite));
    assert.equal(new Set(garment.indices).size,garment.positions.length/3,'every cloth node carries triangle mass');
    const counts=new Map();
    for(let k=0;k<garment.indices.length;k+=3) {
      const ids=garment.indices.slice(k,k+3),p=ids.map(i=>new THREE.Vector3(...garment.positions.slice(3*i,3*i+3)));
      assert.ok(new THREE.Vector3().subVectors(p[1],p[0]).cross(new THREE.Vector3().subVectors(p[2],p[0])).length()>1e-10);
      for(let j=0;j<3;j++){const e=[ids[j],ids[(j+1)%3]].sort((a,b)=>a-b).join(',');counts.set(e,(counts.get(e)||0)+1);}
    }
    assert.ok([...counts.values()].every(c=>c<=2),'manifold sewn mesh');
    {
      const adj=new Map();for(const [e,c] of counts)if(c===1){const [a,b]=e.split(',').map(Number);for(const [i,j] of [[a,b],[b,a]]){if(!adj.has(i))adj.set(i,[]);adj.get(i).push(j);}}
      assert.ok([...adj.values()].every(x=>x.length===2));
      let loops=0;const seen=new Set();for(const i of adj.keys())if(!seen.has(i)){loops++;const q=[i];while(q.length){const k=q.pop();if(seen.has(k))continue;seen.add(k);q.push(...adj.get(k));}}
      assert.equal(loops,garment.id==='shorts'?3:4,'waist and two leg hems; shirt hem, neck and armholes');
    }
  }
});
test('clothes are opaque and independently switchable, without mutating body geometry',()=>{
  const group=new THREE.Group(), source=bodyGeometry(), original=[...source.positions];
  const clothing=new ClothingView(group);clothing.fit(source,{id:'skin'});
  assert.equal(clothing.group.children.length,2);
  assert.ok(clothing.group.children.every(o=>o.visible && o.material.opacity===1 && !o.material.transparent));
  clothing.setEnabled('shorts',false);
  assert.equal(clothing.meshes.get('shorts').visible,false);assert.equal(clothing.meshes.get('shirt').visible,true);
  assert.deepEqual(source.positions,original);
  assert.throws(()=>clothing.setEnabled('unknown',true));
  clothing.dispose();assert.equal(group.children.length,0);
});
const heldSkin=new URL('../../data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz',import.meta.url);
test('held skin arms and hand sections do not become garment shoulder or hip flares',{skip:!fs.existsSync(heldSkin)},()=>{
  const garments=buildGarments(JSON.parse(zlib.gunzipSync(fs.readFileSync(heldSkin))));
  for(const garment of garments)for(let i=0;i<garment.positions.length;i+=3)
    assert.ok(Math.abs(garment.positions[i])<.215,`${garment.id} includes an arm or hand section`);
});
