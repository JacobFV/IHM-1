"""Retained native-state geometry and work-conjugacy checks; never starts native."""
from pathlib import Path
import copy,json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.articulated import CanonicalRegistration
ROOT=Path(__file__).resolve().parents[1]

def main():
    payload=json.loads((ROOT/'data/derived/canonical/mechanics.json').read_bytes())
    record=json.loads((ROOT/'data/derived/native-stream-smoke-n6e0pvqi/supine/smoke.json').read_bytes())
    native=record['initial'];reg=CanonicalRegistration(payload,native);projected=reg.project(native)
    # Not a literal count. It was 2408, then 2403 after the duplicate-surface collapse,
    # and 4000 once the display promotion landed; a literal here only ever gets edited to
    # match whatever the build now says. What this actually needs is that the projection
    # covers the mechanics entity set exactly, and that mechanics agrees with the anatomy
    # it was built from -- which a dropped or unprojected entity still fails.
    anatomy=json.loads((ROOT/'data/derived/canonical/anatomy.json').read_bytes())
    assert len(projected)==len(payload['entities'])==payload['counts']['entities']
    assert {e['id'] for e in payload['entities']}=={e['id'] for e in anatomy['entities']}
    assert max(np.linalg.norm(e['translation_m']) for e in projected.values())<1e-14
    assert abs(np.linalg.det(reg.basis)-1)<1e-14
    assert all(np.array_equal(c,reg.global_map) for c in reg.maps.values())
    # One source point must have one world location across segment associations.
    p=np.array([.15,.3,-.04,1]);assert max(np.linalg.norm(c@p-reg.global_map@p) for c in reg.maps.values())==0
    # Force and twist virtual power conserved in the common frame.
    native=record['step'];contacts=reg.contact_wrenches(native)
    for raw,c in zip(native['contacts'],contacts):
        b=native['bodies'][raw['body_frame']];v=np.array(b['origin_velocity_m_s']);w=np.array(b['angular_velocity_rad_s'])
        p0=np.dot(raw['force_n'],v)+np.dot(raw['moment_nm'],w)
        p1=np.dot(c['force_n'],reg.basis@v)+np.dot(c['moment_nm'],reg.basis@w)
        assert abs(p0-p1)<1e-12
    ident=next(iter(reg.named));point=np.array(reg.specs[ident]['centroid_m']);force=np.array([2.,-3.,7.]);mapped=reg.force(ident,point,force,native)
    assert np.allclose(reg.basis@mapped['point_m']+reg.global_map[:3,3],point)
    assert np.allclose(reg.basis@mapped['force_n'],force)
    # Shared rigid native motion must give the same transform for every mesh owner.
    changed=copy.deepcopy(record['initial']);theta=.02;r=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1.]])
    rigid=np.eye(4);rigid[:3,:3]=r;rigid[:3,3]=[.01,-.02,.03]
    for b in changed['bodies'].values():b['transform_ground']=(rigid@b['transform_ground']).tolist()
    expected=reg.global_map@rigid@np.linalg.inv(reg.global_map)
    assert all(np.allclose(t,expected,atol=1e-14) for t in reg.transforms(changed).values())
    print(json.dumps({'passed':True,'canonical_entities':len(projected),'native_segments':len(reg.bodies),'global_fit':reg.global_fit,'scope':'Retained native-state projection only; no new native launch or settled-support claim'},indent=2))
if __name__=='__main__':main()
