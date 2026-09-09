#!/usr/bin/env python3
"""Collect the landing/load-transfer experiment into one auditable receipt.

Reads only artifacts the runs and designers already wrote, restates their
sha256s, and carries every arm's own scope string through unchanged.
"""
import argparse, hashlib, json, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

DESIGNS = {
    'right_lift': 'data/research/locomotion_control/linearization_8_mijxgp',
    'right_swing': 'data/research/locomotion_control/linearization_j0b9dci6',
    'approach_10mm': 'data/research/locomotion_control/linearization_rnwkf3bb',
    'approach_1mm': 'data/research/locomotion_control/linearization_5j_icido',
    'landing_10percent': 'data/research/locomotion_control/linearization_ufj51r6y',
    'landing_25percent': 'data/research/locomotion_control/linearization_landing25',
    'landing_45percent': 'data/research/locomotion_control/linearization_landing45',
    'landing_10percent_uncrouched': 'data/research/locomotion_control/linearization_g_x24sqn',
}
SCRIPTS = ['scripts/design_native_deflated_lqr.py', 'scripts/run_native_movement_chain.py',
           'scripts/summarize_native_movement_chain.py', 'scripts/evaluate_native_movement_program.py',
           'scripts/discretize_native_stance_lqr.py', 'scripts/linearize_native_stance.py',
           'scripts/solve_native_patient_footprint_crouched_landing98.py',
           'scripts/anchor_native_stance_foot_v2.py', 'scripts/export_native_copied_equilibrium.py',
           'scripts/export_native_copied_lqr_target.py', __file__]


def design_row(directory):
    row = {}
    for name in ('discrete_margin', 'finite_horizon_margin', 'deflated_margin'):
        path = ROOT / directory / name / 'report.json'
        if not path.exists():
            continue
        report = json.loads(path.read_text())
        row[name] = {
            'linearization_sha256': sha(path.parent / 'linearization.npz'),
            'dare_residual_relative': report.get('dare_residual_relative'),
            'reduced_dare_residual_relative': report.get('reduced_dare_residual_relative'),
            'discrete_gain_spectral_radius': report.get('discrete_gain_spectral_radius'),
            'reduced_closed_loop_spectral_radius': report.get('reduced_closed_loop_spectral_radius'),
            'full_closed_loop_radius_excluding_deflated': report.get('full_closed_loop_radius_excluding_deflated'),
            'deflated_dimension': report.get('deflated_dimension'),
            'deflated_input_reachability': report.get('deflated_input_reachability'),
            'discrete_gain_max': report.get('discrete_gain_max'),
        }
    return row


def arm_row(directory):
    path = ROOT / directory / 'chain_report.json'
    if not path.exists():
        return {'directory': directory, 'present': False}
    chain = json.loads(path.read_text())
    summary_path = ROOT / directory / 'summary.json'
    row = {
        'directory': directory, 'present': True,
        'chain_report_sha256': sha(path),
        'all_stages_completed_horizon': chain['all_stages_completed_horizon'],
        'single_native_process': chain['single_native_process'],
        'native_stream_identity': chain['native_stream_identity'],
        'chain_final_time_s': chain['chain_final_time_s'],
        'wall_seconds': chain['wall_seconds'],
        'scope': chain['scope'],
        'stages': [{'stage': s['stage'], 'program': s['program'], 'program_sha256': s['program_sha256'],
                    'elapsed_s': s['elapsed_s'], 'completed_horizon': s['completed_horizon'],
                    'error': s['error'], 'max_clipped_count': s['max_clipped_count'],
                    'final_com_speed_m_s': s['final_com_speed_m_s'],
                    'final_per_foot_clearance_m': s['final_per_foot_clearance_m'],
                    'source_sha256': s['source_sha256']} for s in chain['stages']],
    }
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        row['summary_sha256'] = sha(summary_path)
        row['final_com_forward_travel_m'] = summary['final_com_forward_travel_m']
        row['segments'] = summary['segments']
    return row


def single_run_row(directory):
    """The motivating run is a single program, not a chain, so it has no chain_report.json.

    Its report field `final_per_foot_clearance_m` records the LAST frame, which for
    this run is a frame of a fall, not an airborne foot held clear of the ground.
    Reading the retained trajectory instead, the right foot reached 5.2 mm at
    t=22.6 s and then rose again. The row keeps that correction beside the raw
    number so the receipt cannot be read the way the raw field was.
    """
    path = ROOT / directory / 'report.json'
    if not path.exists():
        return {'directory': directory, 'present': False}
    report = json.loads(path.read_text())
    return {
        'directory': directory, 'present': True, 'kind': 'single_program',
        'report_sha256': sha(path),
        'elapsed_s': report['elapsed_s'], 'completed_horizon': report['completed_horizon'],
        'horizon_s': report['horizon_s'], 'error': report['error'],
        'max_clipped_count': report['max_clipped_count'],
        'final_com_speed_m_s': report['final_com_speed_m_s'],
        'final_per_foot_clearance_m': report['final_per_foot_clearance_m'],
        'source_sha256': report['source_sha256'],
        'interpretation': (
            'Terminated at 24.13 s of a 30 s horizon in a fall. The 0.242 m right-foot '
            'clearance in this report is a fall frame, not a foot held airborne; the '
            'retained trajectory reaches 5.2 mm at t=22.6 s before rising. The runaway '
            'state is hip_adduction_r moving opposite its own reference from blend ~0.16. '
            'Clipping was a symptom, not the trigger.'
        ),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', default='data/research/locomotion_control/landing_experiment_receipt.json')
    a = p.parse_args()
    receipt = {
        'schema': 'ihm.native-landing-experiment.v1',
        'generated_unix_s': int(time.time()),
        'question': 'Why did the native 98-muscle right-foot landing program fail, and what makes it complete?',
        'scripts_sha256': {s: sha(ROOT / s) if s != __file__ else sha(s) for s in SCRIPTS},
        'gain_designs': {name: design_row(d) for name, d in DESIGNS.items()},
        'load_waypoints': {
            name: {
                'registration_sha256': sha(ROOT / f'data/derived/mechanics/{name}/registration.json'),
                'acceptance': json.loads((ROOT / f'data/derived/mechanics/{name}/native_initial_acceptance.json').read_text()),
            } for name in ('patient_footprint_landing10_anchored_stance98',
                           'patient_footprint_landing25_anchored_stance98',
                           'patient_footprint_landing45_anchored_stance98')
            if (ROOT / f'data/derived/mechanics/{name}/native_initial_acceptance.json').exists()
        },
        'arms': [arm_row(d) for d in (
            'data/research/locomotion_control/step_deflated_v1',
            'data/research/locomotion_control/step_dare_control_v1',
            'data/research/locomotion_control/step_no_approach_control_v1',
            'data/research/locomotion_control/step_transfer_v1',
        )],
        'prior_failing_run': single_run_row('data/research/locomotion_control/patient_right_landing_program'),
        'walking_demonstrated': False,
        'sustained_forward_locomotion_demonstrated': False,
        'perturbation_recovery_demonstrated_for_stepping': False,
        'scope': (
            'Single right-foot step attempt in the native 98-muscle upright plant, driven only by muscle '
            'excitation from blended static-equilibrium LQR waypoints. Shows: a settled double-support '
            'touchdown with the right foot ~16 cm ahead of the left; a matched control arm that falls when '
            'the staged contact approach is removed; and a valid Riccati design for the loaded-landing and '
            'load-transfer waypoints after deflating the exact ground-plane translation/yaw symmetry that no '
            'muscle can reach. Does NOT show: walking, sustained forward locomotion, a second step, a left '
            'swing, load beyond 45 percent on the leading foot, or perturbation recovery during stepping.'
        ),
    }
    out = ROOT / a.output
    out.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'output': a.output, 'arms': [
        {'directory': r['directory'], 'present': r['present'],
         'all_stages_completed_horizon': r.get('all_stages_completed_horizon'),
         'chain_final_time_s': r.get('chain_final_time_s'),
         'final_com_forward_travel_m': r.get('final_com_forward_travel_m')} for r in receipt['arms']]}, indent=2))


if __name__ == '__main__':
    main()
