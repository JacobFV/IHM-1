"""Compare actual native onset motion with an explicit solved initial pose."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import tempfile
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.assembly.articulated import CanonicalRegistration
from ihm.assembly.environment_dynamics import EnvironmentDynamics
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG


def main(pose_path):
    artifact=json.loads(pose_path.read_text())
    output=Path(tempfile.mkdtemp(prefix='supine-initial-motion-',dir=ROOT/'data/derived'))
    canonical=json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())
    reports={}
    for label,pose in [('default',None),('candidate',artifact['coordinates'])]:
        native=NativeMechanicalStream(ROOT,output/label,environment='supine',target_mass_kg=MECHANICAL_TARGET_MASS_KG,
            augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json',initial_pose=pose)
        try:
            assert hashlib.sha256((native.output/'inputs/subject_walk_scaled.osim').read_bytes()).hexdigest()==artifact['model_sha256']
            initial=native.snapshot()
            registration=CanonicalRegistration(canonical,{**initial,'bodies':initial['registration_reference_bodies']})
            observer=EnvironmentDynamics(ROOT,'supine',{'scene':'garden-patio','objects':[]},registration)
            points=observer.skin_points(registration.project(initial))
            rows=[]
            for _ in range(5):
                state=native.advance(.02)
                current=observer.skin_points(registration.project(state))
                rows.append({'time_s':state['time_s'],'maximum_skin_sample_step_m':float(np.max(np.linalg.norm(current-points,axis=1))),
                             'kinetic_energy_j':state['kinetic_energy_j'],'contact_force_n':state['contact_force_n'],
                             'momentum_residual_n':float(np.linalg.norm(state['momentum_balance_residual_n']))})
                points=current
            reports[label]={'initial_contact_force_n':initial['contact_force_n'],'rows':rows}
        finally:native.close()
    ratio=reports['candidate']['rows'][0]['maximum_skin_sample_step_m']/reports['default']['rows'][0]['maximum_skin_sample_step_m']
    report={'passed':ratio<.5,'first_step_motion_ratio':ratio,'initial_pose_artifact':str(pose_path),
            'initial_pose_artifact_sha256':hashlib.sha256(pose_path.read_bytes()).hexdigest(),
            'scope':'Actual92muscle native motion with one unchanged source-default registration reference and registered skin observation; no environment feedback or static equilibrium claim','arms':reports}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report|{'output':str(output)},indent=2))
    assert report['passed'],'Candidate did not halve initial registered skin motion'


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('pose',type=Path)
    main(parser.parse_args().pose)
