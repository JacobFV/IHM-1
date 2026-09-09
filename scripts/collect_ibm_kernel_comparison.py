"""Assemble the kernel-comparison receipt from the retained native stance reports.

Usage: collect_ibm_kernel_comparison.py <output.json>

One row per arm of every matched ten-second push trial, with the report digest,
so the comparison in docs/IBM_CURRICULUM16_KERNEL.md can be checked against the
runs it was written from rather than retyped.
"""
import hashlib, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ('elapsed_s', 'completed_horizon', 'peak_com_displacement_m',
          'final_com_displacement_m', 'final_com_speed_m_s', 'pelvis_height_m', 'error')

RUNS = {
    'fused_kernel_positive': 'data/runtime/motor-learning/cortical-stance-native-patient-small-signal-positive-20260908/report.json',
    'fused_kernel_negative': 'data/runtime/motor-learning/cortical-stance-native-patient-small-signal-negative10-20260908/report.json',
    'curriculum16_positive': 'data/runtime/motor-learning/cortical-stance-native-curriculum16-positive-20260908/report.json',
    'curriculum16_negative': 'data/runtime/motor-learning/cortical-stance-native-curriculum16-negative-20260908/report.json',
    'shuffled_control_positive': 'data/runtime/motor-learning/cortical-stance-native-shuffled16-positive-20260908/report.json',
}

out = {'schema': 'ihm.ibm-kernel-comparison.v1',
       'question': 'Does the 16-objective IBM kernel change what the embodied stance policy can do, '
                   'and is the kernel content itself load-bearing?',
       'method': 'One stance pipeline, three kernels. Same teacher, trajectories, site count, encoder, '
                 'epochs, learning rate, readout fit and small-signal scaling; only dyn.embed differs. '
                 'Every arm is matched against a severed control that retains the same tonic baseline.',
       'runs': {}}
for name, path in RUNS.items():
    p = ROOT / path
    if not p.is_file():
        out['runs'][name] = {'status': 'not run'}
        continue
    d = json.loads(p.read_text())
    out['runs'][name] = {'report': path, 'report_sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
                         'artifact_sha256': d['artifact_sha256'],
                         'perturbation_force_n': d['perturbation_force_n'],
                         'seconds_requested': d['seconds_requested'],
                         'arms': {a: {k: r.get(k) for k in FIELDS} for a, r in d['arms'].items()}}
print(json.dumps(out, indent=2))
Path(sys.argv[1]).write_text(json.dumps(out, indent=2) + '\n')
