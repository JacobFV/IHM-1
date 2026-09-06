"""Bounded native copied-state IK for an exact hand surface near a scalp shaft."""
from pathlib import Path
import argparse,json,signal,sys,tempfile,time,xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.materialize_hair_residual_native import materialize
from scripts.verify_hair_coupled_native import CandidateStream
from ihm.assembly.source_skin import file_sha256


def run(build):
    def expired(sig,frame):raise TimeoutError('20s bounded pose query budget exceeded')
    signal.signal(signal.SIGALRM,expired);signal.alarm(20);start=time.monotonic();output=Path(tempfile.mkdtemp(prefix='hair-hand-pose-',dir=ROOT/'data/derived'));identity=materialize(output/'inputs');native=CandidateStream(output/'native',identity,build_path=build)
    try:
        initial=native.snapshot();base=json.loads((ROOT/'data/research/hair_source_factory.json').read_bytes());partition=json.loads((ROOT/'data/research/hair_native_partition.json').read_bytes());scan=json.loads((ROOT/'data/research/hair_collision_face_scan.json').read_bytes());bodies={b['owner']:b for b in base['actual_frozen_native_inertia']['bodies']}
        face=scan['engineered_pose_hand_faces'][0];vertices=np.asarray(face['positions_m']);embedding=np.asarray(bodies['hand_r']['canonical_reference_to_body_local']);local=vertices@embedding[:3,:3].T+embedding[:3,3];localcenter=local.mean(0);localnormal=np.cross(local[1]-local[0],local[2]-local[0]);localnormal/=np.linalg.norm(localnormal)
        strands=partition['groups']['scalp']['prepared']['strands'];x=np.asarray(strands['centerlines_m']).reshape(-1,3);offsets=strands['strand_offsets'];torso=np.asarray(initial['bodies']['torso']['transform_ground'])@np.asarray(bodies['torso']['canonical_reference_to_body_local']);hand=np.asarray(initial['bodies']['hand_r']['transform_ground']);current=hand[:3,:3]@localcenter+hand[:3,3]
        candidates=[]
        for i,(a,b) in enumerate(zip(offsets[:-1],offsets[1:])):
            point=torso[:3,:3]@x[a+2]+torso[:3,3];candidates.append((np.linalg.norm(point-current),i,a+2,point))
        _,guide,node,point=min(candidates,key=lambda row:row[0]);tangent=torso[:3,:3]@(x[node]-x[node-1]);tangent/=np.linalg.norm(tangent);normal=np.cross(tangent,[0,1,0]);normal/=np.linalg.norm(normal);target=point-1e-5*normal
        names=['arm_flex_r','arm_add_r','arm_rot_r','elbow_flex_r','pro_sup_r'];engineered=np.array([[-2.5,2.5],[-2.5,2.5],[-3.1,3.1],[0,2.618],[0,2.09]]);xml=ET.parse(output/'inputs/subject_walk_scaled.osim').getroot();source=np.array([list(map(float,xml.find(f'.//Coordinate[@name="{n}"]/range').text.split())) for n in names]);bounds=np.c_[np.maximum(engineered[:,0],source[:,0]),np.minimum(engineered[:,1],source[:,1])];evaluations=[]
        def evaluate(q):
            if len(evaluations)>=120:raise RuntimeError('120native pose evaluation budget reached')
            response=native._request(' '.join(['hair_pose','hand_r',str(len(names)),*[item for name,value in zip(names,q) for item in (name,str(value))]]));transform=np.asarray(response['transform_ground']);center=transform[:3,:3]@localcenter+transform[:3,3];n=transform[:3,:3]@localnormal;residual=np.r_[(center-target)/.1,n-normal]
            evaluations.append({'q':q.tolist(),'center_native_ground_m':center.tolist(),'normal_native_ground':n.tolist(),'residual':residual.tolist()});return residual
        result=least_squares(evaluate,[1.8,.4,0,1.5,1.],bounds=(bounds[:,0],bounds[:,1]),max_nfev=18,ftol=1e-8,xtol=1e-8,gtol=1e-8);evaluate(result.x);last=evaluations[-1];position_error=float(np.linalg.norm(np.asarray(last['center_native_ground_m'])-target));normal_error=float(np.linalg.norm(np.asarray(last['normal_native_ground'])-normal))
        observed=native._request('observe');assert observed['time_s']==initial['time_s'] and observed['bodies']==initial['bodies'] and observed['coordinates']==initial['coordinates']
        report={'schema':'ihm.hair-hand-pose-feasibility.v1','kinematic_feasible':position_error<1e-5 and normal_error<1e-3,'native_time_advanced':False,'collision_exercised':False,'output':str(output.relative_to(ROOT)),'face':face,'scalp_guide_index':guide,'scalp_sample_id':partition['groups']['scalp']['prepared']['sample_ids'][guide],'free_node_offset':2,'target_hair_point_native_ground_m':point.tolist(),'target_face_center_native_ground_m':target.tolist(),'target_face_normal_native_ground':normal.tolist(),'pose':dict(zip(names,result.x.tolist())),'source_bounds_rad':source.tolist(),'engineered_search_bounds_rad':bounds.tolist(),'position_error_m':position_error,'normal_error':normal_error,'evaluations':evaluations,'wall_s':time.monotonic()-start,'binary_sha256':native.binary_sha256,'source_sha256':{str(p):file_sha256(p) for p in [Path(__file__),ROOT/'data/research/hair_collision_face_scan.json']},'limits':['Source shoulder ROM is broad engineering mapping; search bounds are explicit synthesis, not calibrated anatomical ROM','No other body self-contact or muscle strain acceptance; copied-state geometry feasibility only','Roots and all residual inertias unchanged; face is exact exterior skin with fixed inferred hand_r owner']}
        (output/'receipt.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k] for k in ('output','kinematic_feasible','position_error_m','normal_error','wall_s')}))
    finally:native.close();signal.alarm(0)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--build',type=Path,required=True);parser.add_argument('--run-native',action='store_true');args=parser.parse_args()
    if not args.run_native:raise SystemExit('Requires coordinated20s native pose slot')
    run(args.build)
