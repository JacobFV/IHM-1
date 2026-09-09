"""Measure raw static/short-forward residuals of a frozen supported initialization."""
import json,sys,tempfile,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream

def main(variant):
    directory=(ROOT/variant).resolve();pose=json.loads((directory/'initial_pose.json').read_text());support=pose['support_geometry'];output=Path(tempfile.mkdtemp(prefix='supine-equilibrium-initial-',dir=ROOT/'data/derived'));(output/'protocol.py').write_bytes(Path(__file__).read_bytes())
    def write(name,x):(output/name).write_text(json.dumps(x,indent=2)+'\n')
    stream=NativeMechanicalStream(ROOT,output/'native',environment='supine',target_mass_kg=pose['target_mass_kg'],augmented_registration=str((directory/'registration.json').relative_to(ROOT)),initial_pose=dict(sorted(pose['coordinates'].items())),surface_contact_manifest=support['manifest_path'],bed_material=support['bed_material'])
    try:
        before=stream.snapshot();write('initialized.json',before)
        q={n:before['coordinates'][n]['value'] for n in sorted(pose['coordinates'])};args=['evaluate_static_pose',str(len(q))]
        for n,v in q.items():args.extend([n,str(v)])
        static=stream._request(' '.join(args));write('static_actual_initialized_coordinates.json',static)
        after=stream.advance(1e-6,actuation={n:muscle['activation'] for n,muscle in before['muscles'].items()});write('advanced_1us.json',after)
        acc={n:(after['coordinates'][n]['speed']-row['speed'])/1e-6 for n,row in before['coordinates'].items()};raw_static_max=max(map(abs,static['udot']));raw_force_max=max(map(abs,static['constrained_zero_acceleration_residual_mobility_force']));tones={n:before['muscles'][n]['activation'] for n in pose['activations']}
        source=json.loads((ROOT/'data/derived/supine-lumbar-default-receipt-t435gyzt/initialized.json').read_text());source_pose=json.loads((ROOT/'data/derived/supine-supported-tones-0uq4bcb3/initial_pose.json').read_text());displacements={n:before['coordinates'][n]['value']-v for n,v in source_pose['coordinates'].items()};rot=[abs(v) for n,v in displacements.items() if n not in ('pelvis_tx','pelvis_ty','pelvis_tz')];translation=[abs(v) for n,v in displacements.items() if n in ('pelvis_tx','pelvis_ty','pelvis_tz')]
        report=dict(variant=str(directory.relative_to(ROOT)),pose_sha256=hashlib.sha256((directory/'initial_pose.json').read_bytes()).hexdigest(),support_geometry=support,raw_static_max_abs_udot=raw_static_max,raw_static_max_abs_generalized_force=raw_force_max,finite_difference_1us_accelerations=acc,max_abs_1us_acceleration=max(map(abs,acc.values())),actual26_activations=tones,max_selected_activation_difference=max(abs(tones[n]-v) for n,v in pose['activations'].items()),initial_energy_and_metabolic={k:before[k] for k in ('kinetic_energy_j','potential_energy_j','total_muscle_metabolic_w','muscle_heat_w','signed_active_fiber_power_w','metabolic_reference','metabolic_analysis_mass_kg')},reference_geometry_matches_previous_frozen_lp6=before['registration_reference_bodies']==source['registration_reference_bodies'],actual_q_displacement_from_lp6=displacements,max_rotation_displacement_rad=max(rot),mean_abs_rotation_displacement_rad=sum(rot)/len(rot),max_translation_displacement_m=max(translation),mean_abs_translation_displacement_m=sum(translation)/len(translation),static_thresholds_met=raw_static_max<.1 and raw_force_max<.1,scope='Unloaded no-cloth zero-velocity static initialization and1us forward consistency only; no sustained stability or learned controller claim. Dynamic relaxed seed velocities/energy discarded, then true-wrench/q/tone optimization performed.')
        write('report.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('finite_difference_1us_accelerations','actual26_activations','actual_q_displacement_from_lp6')},indent=2));print(output)
    finally:stream.close()
if __name__=='__main__':main(sys.argv[1])
