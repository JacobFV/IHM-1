#!/usr/bin/env python3
"""Write the audit manifest for a retained cortical stance bundle.

The bundle directory is produced by `scale_cortical_stance_regime.py`. This step
adds nothing to what the runtime loads -- the runtime pins everything through the
artifact's own `provenance` block -- it records, in one readable file, which files
are in the bundle, which IBM kernel it descends from, and which matched native
receipts were accepted for it. Named receipts are copied in and hashed so a later
reader cannot silently be shown a different run.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

FIELDS = ('elapsed_s', 'completed_horizon', 'peak_com_displacement_m',
          'final_com_displacement_m', 'final_com_speed_m_s', 'pelvis_height_m', 'error')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    import torch
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bundle', required=True)
    ap.add_argument('--kernel-bundle', required=True,
                    help='data/models/... directory holding the IBM kernel this was trained from')
    ap.add_argument('--receipt', action='append', default=[], metavar='NAME=PATH',
                    help='Named native evaluation report to copy in and summarise')
    ap.add_argument('--scope', action='append', default=[])
    args = ap.parse_args()

    bundle = Path(args.bundle)
    if not bundle.is_absolute():
        bundle = ROOT / bundle
    artifact = bundle / 'cortical_stance.pt'
    loaded = torch.load(artifact, map_location='cpu', weights_only=True)
    provenance = loaded['provenance']
    kernel_dir = Path(args.kernel_bundle)
    if not kernel_dir.is_absolute():
        kernel_dir = ROOT / kernel_dir
    kernel_manifest = json.loads((kernel_dir / 'manifest.json').read_text())
    if kernel_manifest['kernel_sha256'] != provenance['brain_checkpoint_sha256']:
        raise ValueError('Named kernel bundle is not the one this artifact was trained from')

    validation = {}
    for item in args.receipt:
        name, _, path = item.partition('=')
        source = Path(path)
        if not source.is_absolute():
            source = ROOT / source
        target = bundle / f'{name}.json'
        shutil.copyfile(source, target)
        report = json.loads(target.read_text())
        validation[name] = {'source': str(source.relative_to(ROOT)), 'sha256': sha(target),
                            'perturbation_force_n': report.get('perturbation_force_n'),
                            'arms': {arm: {k: row.get(k) for k in FIELDS}
                                     for arm, row in report.get('arms', {}).items()}}

    state = loaded['state_dict']
    manifest = {
        'schema': 'ihm.ibm-cortical-stance-bundle.v1',
        'native_model_sha256': provenance['model_sha256'],
        'target_mass_kg': provenance['target_mass_kg'],
        'dt_s': provenance['dt_s'],
        'sites': int(state['dyn.embed'].shape[0]),
        'state_width': len(provenance['state_names']),
        'muscles': len(provenance['muscle_names']),
        'sensory_sites': int(state['sensory_sites'].numel()),
        'motor_sites': int(state['motor_sites'].numel()),
        'computation_dtype': provenance.get('computation_dtype'),
        'artifact_sha256': sha(artifact),
        'brain_checkpoint_sha256': provenance['brain_checkpoint_sha256'],
        'kernel_bundle': str(kernel_dir.relative_to(ROOT)),
        'kernel_donor': {k: kernel_manifest.get(k) for k in
                         ('donor_repository', 'donor_commit', 'donor_path', 'donor_sha256',
                          'consolidation_step', 'objectives', 'is_ablation_control',
                          'ablation_basis', 'site_permutation_seed')},
        'source_sha256': {p.name: sha(p) for p in sorted(bundle.iterdir()) if p.is_file()
                          and p.name != 'manifest.json'},
        'validation': validation,
        'scope': args.scope,
    }
    (bundle / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
