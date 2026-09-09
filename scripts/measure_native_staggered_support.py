#!/usr/bin/env python3
"""Measure a retained staggered-support endpoint; no motion is performed."""
import argparse,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--registration',required=True);p.add_argument('--output',required=True);a=p.parse_args();folder=(ROOT/a.registration).parent
    read=lambda path:json.loads(path.read_text())
    s=read(folder/'initial_snapshot.json');landing=read(ROOT/'data/derived/mechanics/patient_rightlanding_anchored_stance98/initial_snapshot.json');initial=read(ROOT/'data/derived/mechanics/patient_stance98/initial_snapshot.json')
    def foot(st,side):return [c for c in st['contacts'] if c['name'].startswith('contact') and c['name'].endswith('_'+side)]
    def centroid(st,side):return np.mean([c['center_m'] for c in foot(st,side)],axis=0)
    def com(st):
        values=[];masses=[]
        for body in st['bodies'].values():
            mass=body['mass_kg'];point=(np.array(body['transform_ground'])@np.r_[body['mass_center_local_m'],1])[:3];values.append(point*mass);masses.append(mass)
        if abs(sum(masses)-st['mass_kg'])>1e-8:raise ValueError('Body masses do not sum to native mass')
        return sum(values)/sum(masses)
    forces={side:sum(c['force_n'][1] for c in foot(s,side)) for side in ('r','l')}
    report={'schema':'ihm.native-static-staggered-support-geometry.v1','registration':a.registration,'actual_com_m':com(s).tolist(),'com_displacement_from_initial_stance_m':(com(s)-com(initial)).tolist(),'com_displacement_from_landing_m':(com(s)-com(landing)).tolist(),'right_centroid_displacement_from_initial_stance_m':(centroid(s,'r')-centroid(initial,'r')).tolist(),'right_centroid_displacement_from_landing_m':(centroid(s,'r')-centroid(landing,'r')).tolist(),'left_centroid_horizontal_displacement_from_initial_stance_m':(centroid(s,'l')-centroid(initial,'l'))[[0,2]].tolist(),'foot_span_horizontal_error_from_landing_m':((centroid(s,'r')-centroid(s,'l'))-(centroid(landing,'r')-centroid(landing,'l')))[[0,2]].tolist(),'foot_normal_forces_n':forces,'actual_right_load_fraction':forces['r']/sum(forces.values()),'physical_time_advanced_s':0,'actual_dynamic_weight_acceptance_demonstrated':False,'measurement_basis':'Untouched native initialization snapshots: actual body mass centers transformed into world frame and actual contact geometry/forces.'}
    out=ROOT/a.output
    if out.exists():raise ValueError('Refusing overwrite')
    out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
