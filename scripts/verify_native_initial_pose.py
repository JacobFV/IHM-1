"""Explicit native initialization must match the isolated static evaluation."""
from pathlib import Path
import json
import sys
import tempfile
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream


def main():
    output=Path(tempfile.mkdtemp(prefix='native-initial-pose-',dir=ROOT/'data/derived'))
    kwargs=dict(environment='supine',target_mass_kg=77.6122029,
                augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json')
    pose={'pelvis_tx':-.005,'elbow_flex_r':.3,'elbow_flex_l':.3}
    reference=NativeMechanicalStream(ROOT,output/'reference',**kwargs)
    try:
        reference_state=reference.snapshot()
        command=['evaluate_static_pose',str(len(pose))]
        for name,value in pose.items():command.extend([name,str(value)])
        evaluated=reference._request(' '.join(command))
        assert reference.snapshot()==reference_state
    finally:reference.close()
    initialized=NativeMechanicalStream(ROOT,output/'initialized',initial_pose=pose,**kwargs)
    try:
        state=initialized.snapshot()
        assert state['initial_pose_applied'] and state['time_s']==0 and state['kinetic_energy_j']==0
        assert state['support_plane_source_x_m']==reference_state['support_plane_source_x_m']
        assert state['registration_reference_bodies']==reference_state['registration_reference_bodies'], 'Initialization changed material registration reference'
        assert state['bodies']['pelvis']['transform_ground']!=state['registration_reference_bodies']['pelvis']['transform_ground'], 'Initialized posture must move relative to the source registration'
        for name,coordinate in state['coordinates'].items():
            assert abs(coordinate['value']-evaluated['coordinates'][name]['value'])<1e-9
        force=sum((np.array(contact['force_n']) for contact in evaluated['contacts']),start=np.zeros(3))
        assert np.linalg.norm(force-state['contact_force_n'])<1e-8
        assert abs(state['metabolic_reference']['M0_w']-state['total_muscle_metabolic_w'])<1e-10
        assert all(state[name]==0 for name in ['external_work_j','muscle_metabolic_energy_j','signed_active_fiber_work_j','muscle_heat_energy_j'])
        report={'passed':True,'initial_pose':pose,'initial_force_n':state['contact_force_n'],
                'scope':'Explicit initial coordinate state matches copied-state static evaluation; no equilibrium claim'}
        (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report|{'output':str(output)},indent=2))
    finally:initialized.close()


if __name__=='__main__':main()
