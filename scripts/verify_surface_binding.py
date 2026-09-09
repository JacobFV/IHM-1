"""Check full skin binding against retained native browser poses, without owners."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from ihm.assembly.surface_binding import SegmentSurfaceBinding, SKIN_ASSET
from ihm.assembly.environment_dynamics import EnvironmentDynamics


def verify(root, frames_path):
    recorded = json.loads(frames_path.read_text())
    frames = list(recorded.values())
    ident = frames[0]['id']
    source = root / 'data/derived/embodied-sessions' / ident / 'runtime/mechanics'
    payload = json.loads((source / 'canonical_mechanics.json').read_text())
    manifest = json.loads((source / 'registration.json').read_text())
    specs = {e['id']: e for e in payload['entities']}
    registration = SimpleNamespace(groups=manifest['groups'], specs=specs, manifest=lambda: manifest)
    def ranking(point):
        return sorted((float(np.linalg.norm(np.maximum(np.maximum(np.asarray(g['bounds_min_m'])-point,
                        point-np.asarray(g['bounds_max_m'])),0))), body) for body,g in registration.groups.items())
    registration._ranking = ranking
    binding = SegmentSurfaceBinding.from_root(root, registration)
    points = np.asarray(json.loads(gzip.decompress((root/SKIN_ASSET).read_bytes()))['positions']).reshape(-1,3)
    bound = binding.bind(points)
    _, sample = np.unique(np.floor(points/.045).astype(int),axis=0,return_index=True)
    world = EnvironmentDynamics.__new__(EnvironmentDynamics)
    world.ids=[];world.offsets=[]
    for p in points[sample]:
        owner=registration._ranking(p)[0][1]
        bone=registration.groups[owner]['canonical_bones'][0]
        world.ids.append(bone);world.offsets.append(p-np.asarray(specs[bone]['centroid_m']))
    world.offsets=np.asarray(world.offsets)
    expected_owners=[registration._ranking(p)[0][1] for p in points]
    assert [binding.manifest()['segments'][i]['id'] for i in bound[0]] == expected_owners
    rows=[]
    for frame in frames:
        entities=frame['mechanics']['entities']
        transformed=binding.project(bound,binding.frame(entities))
        world.world_frame=None
        actual=world.skin_points(entities)
        assert np.array_equal(transformed[sample], actual), 'Canonical material skin points differ from contact support'
        w=frame['environment_state']['world_frame']
        matrix=np.asarray(w['canonical_to_world'])
        class WorldFrame:
            def points_to_world(self, p): return p @ matrix[:3,:3].T + matrix[:3,3]
        world.world_frame=WorldFrame()
        rendered_world=world.world_frame.points_to_world(transformed[sample])
        actual_world=world.skin_points(entities)
        assert np.array_equal(rendered_world,actual_world)
        rows.append(dict(time_s=frame['time_s'],canonical_max_error_m=float(np.max(np.abs(transformed[sample]-actual))),world_max_error_m=float(np.max(np.abs(rendered_world-actual_world)))))
    fixture_indices=sample[np.linspace(0,len(sample)-1,24,dtype=int)]
    fixture=dict(surface_binding=binding.manifest(), rest_vertices=points[fixture_indices].tolist(),
                 surface_transforms=binding.frame(frames[-1]['mechanics']['entities']),
                 expected_canonical_positions=binding.project(bound,binding.frame(frames[-1]['mechanics']['entities']))[fixture_indices].tolist(),
                 source_vertex_indices=fixture_indices.tolist())
    (root/'data/derived/surface-binding-fixture.json').write_text(json.dumps(fixture,indent=2)+'\n')
    # A tie chooses native name, irrespective of manifest insertion order.
    g={'z':dict(canonical_bones=['z'],bounds_min_m=[0,0,0],bounds_max_m=[1,1,1]),'a':dict(canonical_bones=['a'],bounds_min_m=[0,0,0],bounds_max_m=[1,1,1])}
    tie=SegmentSurfaceBinding(SimpleNamespace(groups=g,specs={i:{'centroid_m':[0,0,0]} for i in g}))
    assert tie.manifest()['segments'][tie.bind([[.5,.5,.5]])[0][0]]['id']=='a'
    try: binding.frame({})
    except KeyError: pass
    else: raise AssertionError('Missing bone pose must fail rather than silently use rest geometry')
    return dict(passed=True, full_skin_vertices=len(points), physics_samples=len(sample), frames=rows,
                source_frames=str(frames_path.relative_to(root)), source_frames_sha256=hashlib.sha256(frames_path.read_bytes()).hexdigest(),
                binding_identity=binding.manifest()['binding_identity'], segments=len(binding.manifest()['segments']),
                scope='Exact hard material attachment equality for all retained contact samples in two recorded actual native poses; full-resolution ownership verified against original ranking. Not a skin FEM or anatomical accuracy validation.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--frames',type=Path);parser.add_argument('--output',type=Path,default=Path('data/derived/surface-binding-verification.json'));args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    frames=args.frames or next((root/'data/derived/browser-integration/drag-common-frame').rglob('native-frames.json'))
    report=verify(root,frames);args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
