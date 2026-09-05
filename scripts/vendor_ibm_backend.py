#!/usr/bin/env python3
"""Snapshot actual IBM Python source bytes without modifying the donor checkout."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def vendor(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    def git(*args):
        return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()
    paths = sorted(p for p in (source / 'ibm').rglob('*')
                   if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc')
    if not paths or not (source / 'ibm/runtime/step.py').is_file():
        raise ValueError('source does not contain the IBM runtime package')
    files = {}
    for path in paths:
        rel = path.relative_to(source).as_posix()
        data = path.read_bytes()
        target = destination / 'source' / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        files[rel] = hashlib.sha256(data).hexdigest()
    identity = dict(schema='ihm.ibm-backend.v1', donor_path=str(source), donor_commit=git('rev-parse', 'HEAD'),
                    donor_status=git('status', '--porcelain', '--untracked-files=no'), files=files,
                    package_sha256=hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(),
                    scope='actual complete Python package; no datasets or learned weights',
                    donor_mutated=False)
    (destination / 'manifest.json').write_text(json.dumps(identity, indent=2) + '\n')
    return identity


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT.parent / 'IBM-1')
    parser.add_argument('--output', type=Path, default=ROOT / 'data/derived/canonical/ibm-backend')
    args = parser.parse_args()
    identity = vendor(args.source, args.output)
    print(json.dumps({k: v for k, v in identity.items() if k != 'files'}, indent=2))
