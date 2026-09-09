#!/usr/bin/env python3
"""Retain an IBM-1 association kernel in IHM-1 with pinned donor provenance.

The donor checkpoint in the IBM-1 repository is never modified. Only the shared
`dyn.embed` cortico-cortical kernel is copied; task heads are deliberately left
behind because they are sensory retrieval readouts and carry no motor policy.
The retained file has the same `ibm1/implicit-v1` field shape the existing IHM
adapters already load, so the kernel identity is the only thing that changes.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SCHEMA = 'ibm1/implicit-v1'
KERNEL_FIELDS = ('sites', 'embed', 'schema', 'step', 'sources', 'weights', 'aligned')


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def donor_commit(donor):
    result = subprocess.run(['git', '-C', str(donor), 'rev-parse', 'HEAD'],
                            capture_output=True, text=True)
    if result.returncode:
        raise ValueError('Donor repository commit unavailable')
    status = subprocess.run(['git', '-C', str(donor), 'status', '--porcelain'],
                            capture_output=True, text=True, check=True)
    return result.stdout.strip(), bool(status.stdout.strip())


def main():
    import torch
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ibm-root', default=str(Path.home() / 'Documents/IBM-1'))
    ap.add_argument('--checkpoint', default='ckpt/ibm1_curriculum16.pt',
                    help='Donor checkpoint path relative to --ibm-root')
    ap.add_argument('--output', default='data/models/ibm_curriculum16_kernel_v1')
    ap.add_argument('--scope', default=(
        'Shared cortico-cortical association kernel only. Its training objectives are '
        'sensory retrieval terms; it contains no motor policy and no motor objective. '
        'Engineering sensorimotor simulation, not a validated human brain.'))
    args = ap.parse_args()

    donor = Path(args.ibm_root).resolve()
    source = (donor / args.checkpoint).resolve()
    if not source.is_file():
        raise ValueError(f'Donor checkpoint missing: {source}')
    equations = donor / 'scripts/pretrain_video_loop.py'
    commit, dirty = donor_commit(donor)

    # The donor may still be training and rewriting this file every consolidation.
    # Hash and deserialize ONE in-memory copy so the recorded digest provably
    # describes the bytes the retained kernel was taken from.
    raw = source.read_bytes()
    donor_sha = sha(raw)
    checkpoint = torch.load(io.BytesIO(raw), map_location='cpu', weights_only=False)
    if sha(source.read_bytes()) != donor_sha:
        raise ValueError('Donor checkpoint changed while being read; retry when the donor run is quiescent')
    if checkpoint.get('schema') != SCHEMA:
        raise ValueError(f'Donor is not an implicit kernel (schema {checkpoint.get("schema")})')
    embed = checkpoint['dyn.embed']
    if embed.ndim != 2 or not torch.isfinite(embed).all():
        raise ValueError('Donor kernel is not a finite two-dimensional embedding')
    if int(checkpoint.get('sites', embed.shape[0])) != embed.shape[0]:
        raise ValueError('Donor site count disagrees with the embedding')

    kernel = {'dyn.embed': embed.detach().clone().contiguous()}
    for name in KERNEL_FIELDS:
        if name in checkpoint:
            kernel[name] = checkpoint[name]
    heads = sorted(checkpoint.get('heads', {}))

    out = Path(args.output)
    if not out.is_absolute():
        out = ROOT / out
    out.mkdir(parents=True, exist_ok=False)
    torch.save(kernel, out / 'kernel.pt')
    kernel_bytes = (out / 'kernel.pt').read_bytes()
    # The adapters load with weights_only=True; prove the retained file survives it.
    reloaded = torch.load(out / 'kernel.pt', map_location='cpu', weights_only=True)
    if not torch.equal(reloaded['dyn.embed'], embed) or reloaded.get('schema') != SCHEMA:
        raise ValueError('Retained kernel does not round-trip under the runtime loader')
    equation_bytes = equations.read_bytes()
    (out / 'pretrain_video_loop.py').write_bytes(equation_bytes)
    (out / Path(__file__).name).write_bytes(Path(__file__).read_bytes())

    manifest = {
        'schema': 'ihm.ibm-kernel-bundle.v1',
        'donor_repository': 'IBM-1',
        'donor_commit': commit,
        'donor_worktree_dirty': dirty,
        'donor_path': args.checkpoint,
        'donor_sha256': donor_sha,
        'donor_bytes': len(raw),
        'kernel_sha256': sha(kernel_bytes),
        'kernel_bytes': len(kernel_bytes),
        'cortical_source_sha256': sha(equation_bytes),
        'kernel_schema': SCHEMA,
        'sites': int(embed.shape[0]),
        'embed': int(embed.shape[1]),
        'consolidation_step': checkpoint.get('step'),
        'objectives': list(checkpoint.get('sources', [])),
        'objective_weights': list(checkpoint.get('weights', [])),
        'materialization_heads_in_donor': heads,
        'heads_retained': False,
        'heads_omitted_because': 'Task heads are sensory retrieval readouts; none is a motor head',
        'donor_checkpoint_modified': False,
        'donor_run_state': ('Donor curriculum run was still training when this snapshot was taken; '
                            'the digest above pins the exact bytes copied, not the final run'),
        'scope': args.scope,
    }
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
