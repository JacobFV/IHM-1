"""Projection preserves exact frame clocks and all renderer inputs."""
from pathlib import Path
import json
import sys
import hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.app.projection import display_trajectory

if __name__=='__main__':
    raw=(ROOT/'data/derived/canonical/trajectory.json').read_bytes()
    data=json.loads(raw);projected=display_trajectory(data,source_sha256=hashlib.sha256(raw).hexdigest())
    assert projected['centroids_m']==data['centroids_m']
    assert len(projected['frames'])==len(data['frames'])
    for source,view in zip(data['frames'],projected['frames']):
        assert source['time_s']==view['time_s'] and source['physiology']==view['physiology']
        assert source['respiration']['skin_field']==view['respiration']['skin_field']
        assert source['entities'].keys()==view['entities'].keys()
        for id,state in source['entities'].items():
            for key in ('translation_m','rotation_matrix','deformation_gradient'):
                default=[0,0,0] if key=='translation_m' else [[1,0,0],[0,1,0],[0,0,1]]
                assert state.get(key,default)==view['entities'][id].get(key,default)
    size=len(json.dumps(projected,separators=(',',':')).encode())
    print(json.dumps({'frames':len(projected['frames']),'full_bytes':len(raw),'display_bytes':size,'exact_projection':True}))
