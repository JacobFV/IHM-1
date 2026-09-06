"""Read-only membership audit of frozen posterior quadrature; no regeneration."""
import gzip,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def audit():
    manifest_path=ROOT/'data/derived/supine-surface-contact-exmzq9pq/manifest.json';manifest=json.loads(manifest_path.read_text())
    topology_path=ROOT/'data/research/engineered_skin_territories/materialization.json';topology=json.loads(topology_path.read_text())
    arrays_path=ROOT/manifest['arrays_path'];geometry_path=ROOT/'data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz'
    if sha(arrays_path)!=manifest['arrays_sha256']:raise ValueError('Frozen quadrature hash changed')
    if sha(geometry_path)!=manifest['source_files'][str(geometry_path.relative_to(ROOT))]:raise ValueError('Frozen geometry hash changed')
    receipt=next(r for r in topology['source_receipts'] if r['path']==str(geometry_path.relative_to(ROOT)))
    if receipt['sha256']!=sha(geometry_path):raise ValueError('Topology and contact geometry mismatch')
    array=np.load(arrays_path);faces=array['face_indices'];component=np.array(topology['surface_diagnostic']['face_component_ids']);selected=topology['surface_diagnostic']['selected_component_id'];eligible=np.array(topology['contact_eligible_triangle_ids'])
    if not np.array_equal(eligible,np.where(component==selected)[0]):raise ValueError('Topology eligible mask mismatch')
    geometry=json.loads(gzip.decompress(geometry_path.read_bytes()));xyz=np.array(geometry['positions']).reshape(-1,3);triangles=np.array(geometry['indices']).reshape(-1,3)
    if len(component)!=len(triangles):raise ValueError('Topology face count mismatch')
    inv=np.linalg.inv(np.array(manifest['registration']['source_to_canonical_ground']));source=xyz@inv[:3,:3].T+inv[:3,3]
    exterior=source[triangles[eligible]];bad=np.where(component[faces]!=selected)[0];rows=[]
    for k in bad:
        yz=array['reference_points_source_m'][k,1:];a=exterior[:,1,1:]-exterior[:,0,1:];b=exterior[:,2,1:]-exterior[:,0,1:];d=yz-exterior[:,0,1:];det=a[:,0]*b[:,1]-a[:,1]*b[:,0]
        valid=np.abs(det)>=1e-16;u=np.zeros(len(det));v=u.copy();u[valid]=(d[valid,0]*b[valid,1]-d[valid,1]*b[valid,0])/det[valid];v[valid]=(a[valid,0]*d[valid,1]-a[valid,1]*d[valid,0])/det[valid]
        valid&=(u>=-1e-10)&(v>=-1e-10)&(u+v<=1+1e-10);x=exterior[:,0,0]+u*(exterior[:,1,0]-exterior[:,0,0])+v*(exterior[:,2,0]-exterior[:,0,0]);ids=np.where(valid)[0];candidate=None
        if len(ids):
            winner=ids[np.argmin(x[ids])];candidate=dict(face_index=int(eligible[winner]),source_x_m=float(x[winner]),delta_x_m=float(x[winner]-array['reference_points_source_m'][k,0]))
        rows.append(dict(quadrature_index=int(k),face_index=int(faces[k]),component=int(component[faces[k]]),body=manifest['bodies'][array['body_indices'][k]],projected_area_m2=float(array['area_m2'][k]),source_point_m=array['reference_points_source_m'][k].tolist(),exterior_only_same_ray=candidate))
    return dict(scope='Frozen reference-pose membership and same-grid ray audit only; no native run, remeshing or canonical regeneration',
        topology_basis=topology['surface_diagnostic']['selection_evidence'],exterior_proxy_not_closed_surface=True,
        sample_component_counts={str(c):int(np.sum(component[faces]==c)) for c in np.unique(component[faces])},
        total_samples=len(faces),excluded_samples=rows,excluded_projected_area_m2=float(array['area_m2'][bad].sum()),
        full_minimum_source_x_m=float(source[:,0].min()),exterior_minimum_source_x_m=float(exterior[:,:,0].min()),
        conclusion='No nested-inner samples, but one seam sample means exterior-only physics is not exactly identical. Any replacement changes contact geometry/source identity and requires a fresh validated reference/cache boundary.',
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (manifest_path,topology_path,arrays_path,geometry_path,Path(__file__))})

if __name__=='__main__':print(json.dumps(audit(),indent=2,allow_nan=False))
