#!/usr/bin/env python3
"""Record actual initial geometry versus the retained forward-swing footprint."""
import argparse,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--registration',required=True);a=p.parse_args();folder=(ROOT/a.registration).parent
    s=json.loads((folder/'initial_snapshot.json').read_text());ref=json.loads((ROOT/'data/derived/mechanics/patient_rightswing_anchored_stance98/initial_snapshot.json').read_text())
    def foot(st,side):return [c for c in st['contacts'] if c['name'].startswith('contact') and c['name'].endswith('_'+side)]
    def center(st,side):return np.mean([c['center_m'] for c in foot(st,side)],axis=0)
    forces={side:sum(c['force_n'][1] for c in foot(s,side)) for side in ('r','l')}
    data={'schema':'ihm.native-static-contact-approach-geometry.v1','registration':a.registration,'swing_footprint_reference':'data/derived/mechanics/patient_rightswing_anchored_stance98/initial_snapshot.json','right_centroid_displacement_from_swing_m':(center(s,'r')-center(ref,'r')).tolist(),'left_horizontal_anchor_error_from_swing_m':(center(s,'l')-center(ref,'l'))[[0,2]].tolist(),'horizontal_foot_span_error_from_swing_m':((center(s,'r')-center(s,'l'))-(center(ref,'r')-center(ref,'l')))[[0,2]].tolist(),'right_signed_clearances_m':{c['name']:c['center_m'][1]-c['radius_m'] for c in foot(s,'r')},'right_minimum_signed_clearance_m':min(c['center_m'][1]-c['radius_m'] for c in foot(s,'r')),'foot_normal_forces_n':forces,'actual_right_load_fraction':forces['r']/sum(forces.values()),'force_origin':'Actual native compliant contact law at solved pose; no prescribed GRF','physical_time_advanced_s':0,'actual_dynamic_contact_approach_demonstrated':False,'measurement_basis':'Untouched native initialized contact geometry and force snapshot; copied equilibrium is a separate artifact.'}
    out=folder/'static_contact_approach_geometry.json'
    if out.exists():raise ValueError('Refusing overwrite')
    out.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data,indent=2))
if __name__=='__main__':main()
