"""Joint trajectories replayed on the anatomy through the segment binding.

These are not native canonical runs and are deliberately not offered as such.
A canonical trajectory is the body's own materialized physics: physiology,
internal state, per-entity deformation.  One of these is a *kinematic replay*
-- a controller's joint trajectory pushed through OpenSim forward kinematics
and applied to each anatomical surface as one rigid segment transform.  It
carries no physiology and no deformation, and the reader is told so in the
`basis` every entry hands back.

Written by scripts/render_anatomical_motion.py --emit.
"""
import json
from pathlib import Path

DIRECTORY = 'data/derived/anatomy-segment-binding'
PREFIX = 'trajectory-'
SCHEMA = 'ihm.body-trajectory.v1'


def _entries(root):
    directory = Path(root) / DIRECTORY
    if not directory.is_dir():
        return []
    return sorted(p for p in directory.glob(PREFIX + '*.json') if p.is_file())


def _identity(path):
    return path.stem[len(PREFIX):]


def list_segment_bound(root):
    """What is on disk, with the size the reader is about to pay for."""
    rows = []
    for path in _entries(root):
        ident = _identity(path)
        if not ident:
            continue
        try:
            head = json.loads(path.read_text())
        except (ValueError, OSError):
            continue
        if head.get('schema') != SCHEMA:
            continue
        projection = head.get('projection', {})
        source = projection.get('source_trajectory', {})
        rows.append({
            'id': ident,
            'label': 'Segment-bound · ' + ident.replace('-', ' '),
            'frames': len(head.get('frames', [])),
            'entities': len(head.get('centroids_m', {})),
            'bytes': path.stat().st_size,
            'source_trajectory': source.get('path'),
            'source_sha256': source.get('sha256'),
            'binding': projection.get('binding', {}).get('path'),
            'basis': projection.get('basis'),
        })
    return {
        'trajectories': rows,
        'scope': 'Kinematic replay of a joint trajectory on the anatomical '
                 'surfaces. Not a native canonical run: no physiology, no '
                 'internal state, and every entity moves rigidly with the one '
                 'segment it is bound to.',
    }


def read_segment_bound(root, ident):
    for path in _entries(root):
        if _identity(path) == ident:
            payload = json.loads(path.read_text())
            if payload.get('schema') != SCHEMA:
                return None
            return payload
    return None
