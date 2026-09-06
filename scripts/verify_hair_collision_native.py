"""Exact-source moving hand face collision; isolated engineered initial pose."""
from pathlib import Path
import argparse,json,signal,sys,tempfile,time,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.source_skin import file_sha256
from ihm.assembly.hair_coupled_candidate import coupled_candidate
from scripts.materialize_hair_residual_native import materialize
from scripts.verify_hair_coupled_native import CandidateStream,totals
from scripts.verify_hair_residual_native import validate
from scripts.effective_passive_energy_observation import PassiveEnergyObservation
from scripts.audit_hair_collision_failure import source_joint_forces


def initial_condition(output,pose,*,enforce_force_budget=True):
    identity=materialize(output/'inputs');q=np.array(list(pose['pose'].values()));center=np.array(pose['evaluations'][-1]['center_native_ground_m']);J=np.zeros((3,5));seen=set()
    for row in pose['evaluations'][-7:]:
        d=np.array(row['q'])-q;nonzero=np.flatnonzero(abs(d)>1e-12)
        if len(nonzero)==1:
            i=int(nonzero[0]);J[:,i]=(np.asarray(row['center_native_ground_m'])-center)/d[i];seen.add(i)
    if len(seen)!=5 or np.linalg.matrix_rank(J)!=3:raise ValueError('Complete retained native local pose Jacobian required')
    desired=.2*np.asarray(pose['target_face_normal_native_ground']);speed=np.linalg.lstsq(J,desired,rcond=None)[0]
    if np.linalg.norm(J@speed-desired)>1e-8 or np.max(abs(speed))>5:raise ValueError('Bounded approach velocity infeasible')
    path=output/'inputs/subject_walk_scaled.osim';original=path.read_bytes();xml=ET.fromstring(original)
    for (name,value),velocity in zip(pose['pose'].items(),speed):
        coordinate=xml.find(f'.//Coordinate[@name="{name}"]');coordinate.find('default_value').text=format(value,'.17g');coordinate.find('default_speed_value').text=format(velocity,'.17g')
    path.write_bytes(ET.tostring(xml,encoding='utf-8',xml_declaration=True));identity['files'][path.name]=file_sha256(path);identity['collision_initial_condition']={'pose':pose['pose'],'coordinate_speeds_rad_s':dict(zip(pose['pose'],speed.tolist())),'hand_center_velocity_jacobian_m_per_rad':J.tolist(),'desired_hand_relative_approach_m_s':desired.tolist(),'pre_pose_model_sha256':file_sha256(original),'basis':'Declared engineered articulated initial condition; no continuing-state teleport or hair root reattachment'}
    coordinates={c.get('name'):{'value':float(c.findtext('default_value')),'speed':float(c.findtext('default_speed_value'))} for c in xml.iter('Coordinate')}
    forces=source_joint_forces(output/'inputs/subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml',coordinates);maximum=max(abs(f['generalized_force_nm']) for f in forces)
    identity['source_passive_force_screen']={'maximum_abs_generalized_force_nm':maximum,'engineering_acceptance_budget_nm':100.,'passed':maximum<=100.,'basis':'Bounded experiment qualification only; original source force laws unchanged; not anatomical ROM calibration'}
    (output/'inputs/hair_residual_identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    if enforce_force_budget and maximum>100:raise ValueError('Initial pose violates source passive-force experiment budget')
    return identity


def passive(port,snapshot):
    return port.observe({'muscles':[{'name':name,'type':m['muscle_type'],'passive_energy_j':m['passive_energy_j'],'fiber_length_m':m['fiber_length_m']} for name,m in snapshot['muscles'].items()]})


def main(build,run):
    pose_path=ROOT/'data/research/hair_hand_pose_acceptance.json';pose=json.loads(pose_path.read_bytes());scan_path=ROOT/'data/research/hair_collision_face_scan.json'
    if file_sha256(scan_path)!=pose['source_sha256'][str(scan_path)] or not pose['kinematic_feasible']:raise ValueError('Exact source pose dependency mismatch')
    if not run:
        with tempfile.TemporaryDirectory(prefix='hair-collision-light-',dir=ROOT/'data/derived') as tmp:
            identity=initial_condition(Path(tmp),pose,enforce_force_budget=False);assert identity['target_native_mass_kg']==77.61218724193209
        base=json.loads((ROOT/'data/research/hair_source_factory.json').read_bytes());partition=json.loads((ROOT/'data/research/hair_native_partition.json').read_bytes());registration=json.loads((ROOT/base['actual_frozen_native_inertia']['registration_path']).read_bytes());old=json.loads((ROOT/'data/research/hair_residual_native_acceptance.json').read_bytes());snapshot=json.loads((ROOT/old['output']/'snapshot.json').read_bytes())
        baseline=coupled_candidate(partition,base,registration,snapshot);extra=coupled_candidate(partition,base,registration,snapshot,additional_contact_faces=[pose['face']])
        assert np.array_equal(baseline.hair.position_m,extra.hair.position_m) and np.array_equal(baseline.hair.mass_kg,extra.hair.mass_kg) and len(extra.triangles)==95
        assert not identity['source_passive_force_screen']['passed']
        print('PASS known invalid pose rejected by source force screen; source-bound exact handface pose, fullrank native Jacobian, bounded approach velocities, residual inventory retained; no native run');return
    def expired(sig,frame):raise TimeoutError('45s bounded collision acceptance')
    signal.signal(signal.SIGALRM,expired);signal.alarm(45);output=Path(tempfile.mkdtemp(prefix='hair-collision-native-',dir=ROOT/'data/derived'));identity=initial_condition(output,pose);native=CandidateStream(output/'native',identity,build_path=build);rows=[];started=time.monotonic()
    try:
        base=json.loads((ROOT/'data/research/hair_source_factory.json').read_bytes());partition=json.loads((ROOT/'data/research/hair_native_partition.json').read_bytes());registration=json.loads((ROOT/base['actual_frozen_native_inertia']['registration_path']).read_bytes());first=native.snapshot();(output/'initial.json').write_text(json.dumps(first,indent=2)+'\n');validate(first,identity,base,partition)
        owner=coupled_candidate(partition,base,registration,first,additional_contact_faces=[pose['face']]);port=PassiveEnergyObservation(output/'inputs/subject_walk_scaled.osim');initial=totals(first,owner);initial_passive=passive(port,first)
        index=owner.identity['sample_references'].index(('scalp',pose['scalp_sample_id']));node=owner.hair.offsets[index]+2
        # Exact added face is last only if its sourceface ID sorts last; locate
        # by the retained coordinates instead of making that ordering assumption.
        reg=owner.registration;body='hand_r';transform=reg.transforms(first)[body];vertices=np.asarray(pose['face']['positions_m'])@transform[:3,:3].T+transform[:3,3];normal=np.cross(vertices[1]-vertices[0],vertices[2]-vertices[0]);normal/=np.linalg.norm(normal);gap=float((owner.hair.position_m[node]-vertices.mean(0))@normal)
        if abs(gap-1e-5)>1e-7:raise ValueError('Actual initialized collision gap differs from declared pose')
        for interval in range(2):
            before=native.snapshot();result,audit=owner.advance(native,1e-4,[],{},tolerance_m=1e-10,max_iterations=8);energy=passive(port,result);current=totals(result,owner)
            raw=current['energy']-initial['energy'];corrected=raw+energy['source_passive_correction_j']-initial_passive['source_passive_correction_j']
            assert abs(owner.hair.time_s-result['time_s'])<1e-14 and audit['within_small_deflection']
            assert np.linalg.norm(audit['hair_body_impulse_residual_ns'])<1e-11
            assert abs(result['external_work_j']-before['external_work_j']-audit['native_hair_work_j'])<1e-11
            rows.append({'native_time_s':result['time_s'],'hair_time_s':owner.hair.time_s,'raw_total_mechanical_energy_change_j':raw,'corrected_total_mechanical_energy_change_j':corrected,'native_signed_active_fiber_work_j':result['signed_active_fiber_work_j'],'unexplained_corrected_mechanical_minus_signed_active_j':corrected-result['signed_active_fiber_work_j'],'linear_momentum_change':(current['momentum']-initial['momentum']).tolist(),'angular_momentum_change':(current['angular']-initial['angular']).tolist(),'passive_energy_observation':energy,'hair_audit':owner.frame()['audit']})
        contacts=sum(r['hair_audit']['contact_count'] for r in rows);impulse=sum(np.linalg.norm(r['hair_audit']['surface_contact_impulse_ns']) for r in rows)
        report={'schema':'ihm.hair-exact-hand-collision.v1','collision_exercised':contacts>0 and impulse>0,'initial_gap_m':gap,'output':str(output.relative_to(ROOT)),'intervals':rows,'initial_passive_energy_observation':initial_passive,'input_identity':identity,'wall_s':time.monotonic()-started,'binary_sha256':native.binary_sha256,'live_enabled':False,'twenty_ms_ready':False,'limitations':['Engineered arm-near-scalp pose; other-body gross contact and muscle strain not certified','Sparse95sourcefaces only; discrete centerline faceinterior contact, noedges/CCD/radius/selfcontact','Thelen source-corrected passive observer leaves raw nativePE unchanged; unexplained energy includes native damping/path virtualwork/numerical terms, not declared closed']}
        (output/'receipt.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');(output/'initial.json').write_text(json.dumps(first,indent=2)+'\n');(output/'final.json').write_text(json.dumps(native.snapshot(),indent=2)+'\n')
        print(json.dumps({'output':str(output),'collision_exercised':report['collision_exercised'],'contacts':contacts,'impulse_ns':impulse,'wall_s':report['wall_s']}))
    except BaseException as error:
        (output/'failure.json').write_text(json.dumps({'error_type':type(error).__name__,'error':str(error),'accepted_intervals':rows,'native_time_s':native.snapshot()['time_s'],'status':'Incomplete acceptance; no collision or closure claim'},indent=2,allow_nan=False)+'\n');raise
    finally:native.close();signal.alarm(0)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--build',type=Path,default=ROOT/'data/derived/hair-observer-build-v1');parser.add_argument('--run-native',action='store_true');a=parser.parse_args();main(a.build,a.run_native)
