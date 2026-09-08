"""Compare native articulation with/without server-owned bedroom contact loads."""
import json
from pathlib import Path
import tempfile
import time
import numpy as np
from ihm.assembly.articulated import ArticulatedBodyPlant
from ihm.assembly.environment_dynamics import EnvironmentDynamics

root=Path(__file__).resolve().parents[1]
output=Path(tempfile.mkdtemp(prefix='environment-native-',dir=root/'data/derived'))
plant=None;start=time.monotonic()
try:
    plant=ArticulatedBodyPlant(root,output/'plant',environment='supine')
    initial=plant.snapshot();checkpoint=plant.checkpoint()
    baseline=plant.advance(.02)
    plant.restore(checkpoint)
    owner=EnvironmentDynamics(root,'supine',{'scene':'bedroom','objects':[]},plant.registration)
    ports=owner.advance(.02,initial['entities'])
    loaded=plant.advance(.02,forces=ports)
    delta=max(np.linalg.norm(np.array(loaded['entities'][k]['centroid_m'])-baseline['entities'][k]['centroid_m']) for k in loaded['entities'])
    assert ports,'Bedroom must physically load the body'
    assert delta>1e-10,'Environment must change native body trajectory'
    for _ in range(9):
        ports=owner.advance(.02,loaded['entities']);loaded=plant.advance(.02,forces=ports)
    plant.release(checkpoint)
    report={'passed':True,'native_body_max_displacement_difference_m':float(delta),'duration_s':loaded['time_s'],
        'final_contact_ports':len(ports),'body_impulse_ns':owner.last_impulse.tolist(),'wall_s':time.monotonic()-start,
        'scope':'Native OpenSim mechanical causal comparison, 0.2 s; not full physiological or long-horizon validation'}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (output/'environment-frame.json').write_text(json.dumps(owner.frame()))
    print(json.dumps({'output':str(output),**report},indent=2))
finally:
    if plant is not None:plant.close()
