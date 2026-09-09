"""Geometry/reciprocity receipts for shared continuous surface attachment."""
import gzip,json
from pathlib import Path
import numpy as np
from ihm.assembly.continuous_surface_binding import ContinuousSurfaceBinding,DEFAULT_ASSET
from ihm.assembly.surface_binding import SegmentSurfaceBinding,SKIN_ASSET
from types import SimpleNamespace

root=Path(__file__).resolve().parents[1]
payload=json.loads(gzip.decompress((root/DEFAULT_ASSET).read_bytes()));binding=ContinuousSurfaceBinding(payload)
record_path=next((root/'data/derived/browser-integration/drag-common-frame').rglob('native-frames.json'))
record=json.loads(record_path.read_text());frame=record['last'];base=root/'data/derived/embodied-sessions'/frame['id']/'runtime/mechanics'
canonical=json.loads((base/'canonical_mechanics.json').read_text());manifest=json.loads((base/'registration.json').read_text());reg=SimpleNamespace(groups=manifest['groups'],specs={e['id']:e for e in canonical['entities']})
hard=SegmentSurfaceBinding(reg);geometry=json.loads(gzip.decompress((root/SKIN_ASSET).read_bytes()));x=binding.rest;tri=np.asarray(geometry['indices']).reshape(-1,3);edges=np.unique(np.sort(np.concatenate([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]]),axis=1),axis=0);length=np.linalg.norm(x[edges[:,0]]-x[edges[:,1]],axis=1);nonzero=length>1e-10
full=binding.sampled_binding(np.arange(len(x)));rows=[]
for name,frame in record.items():
    poses=binding.frame(frame['mechanics']['entities']);posed=binding.project(full,poses);old=hard.project(hard.bind(x),hard.frame(frame['mechanics']['entities']))
    def metrics(y):
        lens=np.linalg.norm(y[edges[:,0]]-y[edges[:,1]],axis=1);ratio=lens[nonzero]/length[nonzero]
        return dict(maximum_edge_stretch_ratio=float(ratio.max()),edge_stretch_ratio_p99=float(np.quantile(ratio,.99)),minimum_edge_stretch_ratio=float(ratio.min()),maximum_edge_m=float(lens.max()))
    rows.append(dict(time_s=frame['time_s'],hard=metrics(old),continuous=metrics(posed)))
    assert metrics(posed)['maximum_edge_stretch_ratio']<metrics(old)['maximum_edge_stretch_ratio']
# Random independent native segment twists test exact work-conjugate scatter.
rng=np.random.default_rng(52);indices=rng.choice(len(x),24,replace=False);sample=binding.sampled_binding(indices);forces=rng.normal(size=(24,3));poses=binding.frame(record['last']['mechanics']['entities']);stations=binding.stations(sample,poses);surface=binding.project(sample,poses);ports=binding.scatter_forces(sample,poses,forces)
vel=rng.normal(size=(len(binding.segments),3));omega=rng.normal(size=(len(binding.segments),3));centers=np.array([poses[s['id']]['centroid_m'] for s in binding.segments]);point_vel=vel[None,:,:]+np.cross(omega[None,:,:],stations-centers[None,:,:]);skin_vel=np.sum(sample[1][:,:,None]*point_vel,axis=1)
expected_power=float(np.sum(forces*skin_vel));by_bone={s['bone_id']:j for j,s in enumerate(binding.segments)};actual_power=0
for port in ports:
    j=by_bone[port['id']];v=vel[j]+np.cross(omega[j],np.asarray(port['point_m'])-centers[j]);actual_power+=float(np.dot(port['force_n'],v))
resultant=np.sum([p['force_n'] for p in ports],axis=0);moment=np.sum([np.cross(p['point_m'],p['force_n']) for p in ports],axis=0)
force_error=float(np.max(np.abs(resultant-forces.sum(axis=0))));moment_error=float(np.max(np.abs(moment-np.cross(surface,forces).sum(axis=0))));power_error=abs(actual_power-expected_power)
assert max(force_error,moment_error,power_error)<1e-12
# A common rigid transform must reproduce that transform of every rest vertex.
a=.31;r=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1.]])
t=np.array([.1,-.2,.3]);common={s['id']:dict(centroid_m=(r@np.asarray(s['reference_centroid_m'])+t).tolist(),rotation_matrix=r.tolist()) for s in binding.segments}
rigid_error=float(np.max(np.abs(binding.project(sample,common)-(sample[0]@r.T+t))))
assert rigid_error<1e-12
fixture=dict(surface_binding=binding.manifest(),sidecar=payload,surface_transforms=poses,source_vertex_indices=indices.tolist(),expected_canonical_positions=surface.tolist())
# Full weights remain in the frozen asset; parity fixture only includes sampled rows.
fixture['sidecar']={**{k:v for k,v in payload.items() if k not in ('weights','reference_positions_m')},'reference_positions_m':sample[0].tolist(),'weights':sample[1].tolist()}
(root/'data/derived/continuous-surface-binding-fixture.json').write_text(json.dumps(fixture,indent=2)+'\n')
report=dict(passed=True,binding_identity=binding.identity,algorithm=payload['algorithm'],recorded_poses=rows,force_resultant_error_n=force_error,moment_error_nm=moment_error,virtual_power_error_w=power_error,common_rigid_motion_error_m=rigid_error,scope=payload['scope'])
(root/'data/derived/continuous-surface-binding-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
