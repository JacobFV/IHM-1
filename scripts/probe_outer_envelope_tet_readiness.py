"""Does the already-repaired whole-body skin envelope tetrahedralise as it stands?

The non-muscular tet-readiness lane deliberately excludes the skin entity body-bp3d-FJ2810 because a
watertight outer envelope was already derived at data/derived/outer-envelope. That envelope's own
receipt lists `triangle_self_intersection_not_exhaustively_tested` as unverified, so reuse is only
safe if the envelope is actually a valid PLC. This probe answers that question and nothing else: it
does NOT re-derive the envelope. It diagnoses the stored envelope, hands it to TetGen, then applies
the muscular repair order to that same envelope and hands the result to TetGen again.

Run with the isolated libigl environment, not the physiological runtime:
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/probe_outer_envelope_tet_readiness.py --self-test
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/probe_outer_envelope_tet_readiness.py \
  --output data/derived/entity-tet-ready-skin-probe-v1
"""
from pathlib import Path
import argparse
import gzip
import importlib.metadata
import json
import shutil
import sys
import time

import numpy as np
import igl.copyleft.cgal as cgal
from igl.copyleft import tetgen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_muscle_tet_ready_surfaces import sha, write, diagnose, signed_volume, repair  # noqa: E402
from verify_muscle_tet_ready_surfaces import tetrahedralize  # noqa: E402
from build_entity_tet_ready_surfaces import analyse_patches, shape_metrics  # noqa: E402

ENVELOPE = ROOT / 'data/derived/outer-envelope/outer-envelope.json.gz'


def load():
    payload = json.loads(gzip.decompress(ENVELOPE.read_bytes()))
    v = np.ascontiguousarray(np.asarray(payload['positions'], float).reshape(-1, 3))
    f = np.ascontiguousarray(np.asarray(payload['indices'], np.int64).reshape(-1, 3))
    return payload, v, f


def self_test():
    v = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [0, 0, 1.]])
    f = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]], np.int64)
    assert tetrahedralize(v, f, 'pYq1.414')['succeeded']
    assert diagnose(v, f)['closed'] and analyse_patches(v, f)['patches'] == 1
    assert ENVELOPE.is_file(), 'the outer envelope this probe reads must exist'
    print('PASS harness and envelope-presence checks')


def probe(output, flags):
    out = Path(output).resolve()
    if out.exists():
        raise ValueError('Choose a fresh output directory')
    packages = {n: importlib.metadata.version(n) for n in ('libigl', 'numpy')}
    if packages['libigl'] != '2.6.2':
        raise ValueError('This probe pins libigl 2.6.2')
    started = time.monotonic()
    payload, v, f = load()
    stored = diagnose(v, f)
    stored_trial = tetrahedralize(v, f, flags)
    rv, rf, log = repair(v, f)
    repaired = diagnose(rv, rf)
    repaired_trial = tetrahedralize(rv, rf, flags)
    out.mkdir(parents=True)
    (out / 'inputs').mkdir()
    shutil.copyfile(__file__, out / 'inputs/probe_outer_envelope_tet_readiness.py')
    shutil.copyfile(ROOT / 'data/derived/outer-envelope/manifest.json', out / 'inputs/outer-envelope-manifest.json')
    geometry = None
    if repaired_trial['succeeded']:
        body = dict(schema='ihm.entity-tet-ready-surface.v1', entity_id='body-bp3d-FJ2810',
                    name='skin', role='skin', system='integumentary', units='m',
                    frame=payload.get('frame'), representation='triangular_surface',
                    source_path='data/derived/outer-envelope/outer-envelope.json.gz',
                    source_sha256=sha(ENVELOPE),
                    positions=[float(x) for x in rv.ravel()], indices=[int(i) for i in rf.ravel()])
        geometry = out / 'geometry' / 'body-bp3d-FJ2810.json.gz'
        geometry.parent.mkdir()
        geometry.write_bytes(gzip.compress(json.dumps(body, allow_nan=False).encode(), mtime=0))
    volume = signed_volume(rv, rf)
    summary = dict(
        schema='ihm.entity-tet-ready-skin-probe.v1', tetgen_flags=flags,
        source='data/derived/outer-envelope/outer-envelope.json.gz', source_sha256=sha(ENVELOPE),
        envelope_not_re_derived=True,
        stored=dict(diagnostics=stored, shape=shape_metrics(v, f),
                    patch_analysis=analyse_patches(v, f), tetgen=stored_trial),
        repaired=dict(diagnostics=repaired, shape=shape_metrics(rv, rf),
                      patch_analysis=analyse_patches(rv, rf), tetgen=repaired_trial,
                      operations=log,
                      signed_volume_drift_m3=volume - stored['signed_volume_m3'],
                      signed_volume_relative_drift=(volume - stored['signed_volume_m3'])
                      / stored['signed_volume_m3'] if stored['signed_volume_m3'] else None,
                      tet_volume_vs_surface_relative_error=(repaired_trial['tet_volume_m3'] - volume) / volume
                      if repaired_trial.get('tet_volume_m3') and volume else None,
                      output_path=str(geometry.relative_to(out)) if geometry else None,
                      output_sha256=sha(geometry) if geometry else None),
        elapsed_s=time.monotonic() - started)
    write(out / 'summary.json', summary)
    artifacts = {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
    write(out / 'manifest.json', dict(
        schema='ihm.entity-tet-ready-skin-probe.v1', packages=packages, python=sys.version,
        outer_envelope_manifest_sha256=sha(ROOT / 'data/derived/outer-envelope/manifest.json'),
        native_cgal_sha256=sha(Path(cgal.pyigl_copyleft_cgal.__file__)),
        native_tetgen_sha256=sha(Path(tetgen.pyigl_copyleft_tetgen.__file__)),
        artifacts_sha256=artifacts, limitations=[
            'This probe reads the existing outer envelope and does not re-derive it; the aperture '
            'capping, sheet assignment and voxel oracle of data/derived/outer-envelope are unchanged.',
            'A repaired envelope written here is a tet-readiness candidate for the skin lane to '
            'accept or reject, not an adopted replacement for the published envelope.']))
    print(json.dumps(dict(
        stored_closed=stored['closed'],
        stored_self_intersecting_face_pairs=stored['self_intersecting_face_pairs'],
        stored_tetgen_succeeded=stored_trial['succeeded'],
        stored_tetgen_message=stored_trial.get('message'),
        repaired_self_intersecting_face_pairs=repaired['self_intersecting_face_pairs'],
        repaired_closed=repaired['closed'],
        repaired_tetgen_succeeded=repaired_trial['succeeded'],
        repaired_tets=repaired_trial.get('tets'),
        signed_volume_relative_drift=summary['repaired']['signed_volume_relative_drift'],
        tet_volume_vs_surface_relative_error=summary['repaired']['tet_volume_vs_surface_relative_error'],
        elapsed_s=summary['elapsed_s']), indent=2))
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--output', type=Path)
    p.add_argument('--flags', default='pYq1.414')
    a = p.parse_args()
    if a.self_test:
        self_test()
    if a.output:
        probe(a.output, a.flags)
    if not a.self_test and not a.output:
        p.error('Select --self-test or --output')
