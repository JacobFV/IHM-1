"""Independent proof that the repaired non-muscular surfaces tetrahedralise.

The forked, stdout-capturing TetGen harness is imported verbatim from
scripts/verify_muscle_tet_ready_surfaces.py; neither that module nor the muscular builder is
modified. Every recorded diagnostic is recomputed here from the stored geometry rather than
trusted, then each surface is handed to TetGen twice: once as the losslessly welded canonical
surface (the control) and once as the repaired surface.

Run with the isolated libigl environment, not the physiological runtime:
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/verify_entity_tet_ready_surfaces.py --self-test
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/verify_entity_tet_ready_surfaces.py \
  --build data/derived/entity-tet-ready-v1 --output data/derived/entity-tet-ready-verification-v1 \
  --all --workers 10

Added over the muscular verification, because this set is not homogeneous:
  * role_strata()  a stratified sample that covers every role and every face-count decile within
                   each role, and force-includes every entity the build flagged
                   cavity_convention_pending, hole_filling_required or self-union-repaired.
  * per-role and per-decile success tables against the welded control.
  * radius_proxy_m carried into every trial so thin-tube failures are separable from blob failures.
"""
from pathlib import Path
import argparse
import gzip
import importlib.metadata
import json
import multiprocessing as mp
import shutil
import sys
import time

import numpy as np
import igl
import igl.copyleft.cgal as cgal
from igl.copyleft import tetgen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_muscle_tet_ready_surfaces import sha, write, diagnose, signed_volume  # noqa: E402
from verify_muscle_tet_ready_surfaces import load_surface, welded_control, tetrahedralize  # noqa: E402
from build_entity_tet_ready_surfaces import shape_metrics  # noqa: E402

_STATE = {}


def role_strata(records, per_decile):
    """Every role, every face-count decile inside that role, plus every flagged entity."""
    picked = {}
    by_role = {}
    for r in records:
        by_role.setdefault(r['role'], []).append(r)
    for role, group in sorted(by_role.items()):
        group = sorted(group, key=lambda r: (r['after']['faces'], r['entity_id']))
        n = len(group)
        for d in range(10):
            band = group[n * d // 10:max(n * d // 10 + 1, n * (d + 1) // 10)]
            if not band:
                continue
            step = max(1, len(band) // per_decile)
            for r in band[::step][:per_decile]:
                picked.setdefault(r['entity_id'], r)
    for r in records:
        if (r['cavity_convention_pending'] or r['hole_filling_required']
                or not r['vertex_links_manifold']
                or any(s['stage'] == 'extract_outer_manifold' for s in r['operations'])):
            picked.setdefault(r['entity_id'], r)
    return sorted(picked.values(), key=lambda r: (r['role'], r['after']['faces'], r['entity_id']))


def _trial(record):
    build, flags = _STATE['build'], _STATE['flags']
    payload, rv, rf = load_surface(build / record['output_path'])
    source = json.loads(gzip.decompress((ROOT / record['source_path']).read_bytes()))
    cv, cf = welded_control(np.asarray(source['positions'], float).reshape(-1, 3),
                            np.asarray(source['indices'], np.int64).reshape(-1, 3))
    control = tetrahedralize(cv, cf, flags)
    repaired = tetrahedralize(rv, rf, flags)
    volume = signed_volume(rv, rf)
    return dict(
        entity_id=record['entity_id'], name=record['name'], role=record['role'], system=record['system'],
        faces_welded_control=int(len(cf)), faces_repaired=int(len(rf)),
        prior_blocking_reasons=record['prior_blocking_reasons'],
        prior_topological_candidate=record['prior_topological_candidate'],
        self_intersecting_pairs_before=record['before']['self_intersecting_face_pairs'],
        nonmanifold_vertices_after=record['after']['nonmanifold_vertices'],
        nonmanifold_edges_after=record['after']['nonmanifold_edges'],
        boundary_edges_after=record['after']['boundary_edges'],
        face_components_after=record['after']['face_components'],
        radius_proxy_m=record['shape_after']['radius_proxy_m'],
        min_positive_face_area_m2=record['shape_after']['min_positive_face_area_m2'],
        cavity_convention_pending=record['cavity_convention_pending'],
        void_patches_after=record['after_patch_analysis'].get('void_patches'),
        enclosed_void_volume_m3=record['after_patch_analysis'].get('enclosed_volume_m3'),
        signed_volume_before_m3=record['before']['signed_volume_m3'],
        signed_volume_after_m3=volume,
        signed_volume_relative_drift=record['signed_volume_relative_drift'],
        tet_volume_vs_surface_relative_error=(repaired['tet_volume_m3'] - volume) / volume
        if repaired.get('tet_volume_m3') and volume else None,
        welded_control=control, repaired=repaired)


def _rediagnose(record):
    _, v, f = load_surface(_STATE['build'] / record['output_path'])
    current = diagnose(v, f)
    shape = shape_metrics(v, f)
    differing = {k: (record['after'][k], current[k]) for k in current
                 if (abs(current[k] - record['after'][k]) > 1e-18 if k == 'signed_volume_m3'
                     else current[k] != record['after'][k])}
    return dict(entity_id=record['entity_id'], role=record['role'], **current), differing, shape


def _init(build, flags):
    _STATE['build'] = Path(build)
    _STATE['flags'] = flags


def self_test():
    v = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [0, 0, 1.]])
    f = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]], np.int64)
    good = tetrahedralize(v, f, 'pYq1.414')
    assert good['succeeded'] and abs(good['tet_volume_m3'] - signed_volume(v, f)) < 1e-12
    bad = tetrahedralize(v, f[:3], 'pYq1.414')
    assert not bad['succeeded'] and bad['tetgen_stdout']
    records = []
    for role, count in (('vascular', 40), ('nerve', 7), ('fluid_cavity', 2)):
        for i in range(count):
            records.append(dict(entity_id='%s-%03d' % (role, i), role=role,
                                after=dict(faces=10 * (i + 1)), operations=[],
                                cavity_convention_pending=(role == 'fluid_cavity' and i == 1),
                                hole_filling_required=False, vertex_links_manifold=True))
    chosen = role_strata(records, 1)
    roles = {r['role'] for r in chosen}
    assert roles == {'vascular', 'nerve', 'fluid_cavity'}, roles
    assert len({r['entity_id'] for r in chosen}) == len(chosen), 'no duplicates'
    assert sum(r['role'] == 'vascular' for r in chosen) == 10, 'ten deciles of the big role'
    assert sum(r['role'] == 'nerve' for r in chosen) == 7, 'a small role is covered exhaustively'
    assert any(r['cavity_convention_pending'] for r in chosen), 'flagged entities are force-included'
    print('PASS TetGen success/failure capture and role-and-decile stratification checks')


def verify(build, output, workers, flags, run_all, per_decile):
    build = Path(build).resolve()
    out = Path(output).resolve()
    if out.exists():
        raise ValueError('Choose a fresh verification directory')
    manifest = json.loads((build / 'manifest.json').read_text())
    for name, digest in manifest['artifacts_sha256'].items():
        path = (build / name).resolve()
        if not path.is_relative_to(build) or sha(path) != digest:
            raise ValueError('Build artifact changed: ' + name)
    for name, digest in manifest['source_geometry_sha256'].items():
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or sha(path) != digest:
            raise ValueError('Canonical source changed: ' + name)
    if sha(ROOT / 'data/derived/canonical/anatomy.json') != manifest['anatomy_sha256']:
        raise ValueError('Canonical anatomical identity changed')
    packages = {n: importlib.metadata.version(n) for n in ('libigl', 'numpy')}
    if packages['libigl'] != '2.6.2':
        raise ValueError('This verification pins libigl 2.6.2')
    records = [json.loads(l) for l in (build / 'entities.jsonl').read_text().splitlines()]
    records = [r for r in records if not r.get('failed')]
    out.mkdir(parents=True)
    (out / 'inputs').mkdir()
    shutil.copyfile(__file__, out / 'inputs/verify_entity_tet_ready_surfaces.py')
    shutil.copyfile(ROOT / 'scripts/verify_muscle_tet_ready_surfaces.py',
                    out / 'inputs/verify_muscle_tet_ready_surfaces.py')
    shutil.copyfile(build / 'manifest.json', out / 'inputs/build-manifest.json')
    shutil.copyfile(build / 'summary.json', out / 'inputs/build-summary.json')
    started = time.monotonic()
    context = mp.get_context('fork')
    with context.Pool(workers, initializer=_init, initargs=(str(build), flags)) as pool:
        rediagnosed = []
        mismatches = []
        shapes = {}
        for current, differing, shape in pool.imap(_rediagnose, records, chunksize=8):
            rediagnosed.append(current)
            shapes[current['entity_id']] = shape
            if differing:
                mismatches.append(dict(entity_id=current['entity_id'], differing=differing))
        with (out / 'rediagnosed.jsonl').open('w') as handle:
            for r in rediagnosed:
                handle.write(json.dumps(r, allow_nan=False) + '\n')
        if mismatches:
            write(out / 'diagnostic-mismatches.json', mismatches)
        print('Rediagnosed', len(rediagnosed), 'mismatches', len(mismatches), flush=True)
        chosen = records if run_all else role_strata(records, per_decile)
        trials = []
        for trial in pool.imap(_trial, chosen, chunksize=2):
            trials.append(trial)
            if len(trials) % 100 == 0:
                print('TetGen', len(trials), 'of', len(chosen), '%.0fs' % (time.monotonic() - started),
                      flush=True)
    order = {r['entity_id']: i for i, r in enumerate(chosen)}
    trials.sort(key=lambda t: order[t['entity_id']])
    with (out / 'tetgen-trials.jsonl').open('w') as handle:
        for t in trials:
            handle.write(json.dumps(t, allow_nan=False) + '\n')
    write(out / 'summary.json', summarise(build, records, trials, rediagnosed, mismatches, flags,
                                          run_all, time.monotonic() - started))
    artifacts = {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
    write(out / 'manifest.json', dict(
        schema='ihm.entity-tet-ready-verification.v1', packages=packages, python=sys.version,
        build_manifest_sha256=sha(build / 'manifest.json'), anatomy_sha256=manifest['anatomy_sha256'],
        muscle_verifier_sha256=sha(ROOT / 'scripts/verify_muscle_tet_ready_surfaces.py'),
        native_cgal_sha256=sha(Path(cgal.pyigl_copyleft_cgal.__file__)),
        native_tetgen_sha256=sha(Path(tetgen.pyigl_copyleft_tetgen.__file__)),
        artifacts_sha256=artifacts, limitations=[
            'A successful tetrahedralisation proves the surface is a valid closed PLC for TetGen; it '
            'does not establish anatomical correctness, inter-entity disjointness or material ownership.',
            'The welded control is the losslessly welded canonical surface with repeated-index and '
            'exactly collinear faces dropped, so it isolates the contribution of the later repair steps.',
            'Entities flagged cavity_convention_pending are tetrahedralised in the solid-inclusion '
            'reading the build wrote. A cavity-preserving reading would need the recorded seed points '
            'passed to TetGen as hole markers and is NOT proved here.',
            'Tet volume is compared against the repaired surface integral, not against any tissue mass '
            'ledger.']))
    summary = json.loads((out / 'summary.json').read_text())
    print(json.dumps({k: v for k, v in summary.items() if k != 'failures'}, indent=2))
    return summary


def summarise(build, records, trials, rediagnosed, mismatches, flags, run_all, elapsed):
    successes = [t for t in trials if t['repaired']['succeeded']]
    errors = [abs(t['tet_volume_vs_surface_relative_error']) for t in successes
              if t['tet_volume_vs_surface_relative_error'] is not None]
    drifts = [abs(r['signed_volume_relative_drift']) for r in records
              if r['signed_volume_relative_drift'] is not None]

    def table(key):
        out = {}
        for t in trials:
            bucket = out.setdefault(t[key], dict(trials=0, repaired_ok=0, control_ok=0,
                                                 control_aborts=0, repaired_aborts=0))
            bucket['trials'] += 1
            bucket['repaired_ok'] += bool(t['repaired']['succeeded'])
            bucket['control_ok'] += bool(t['welded_control']['succeeded'])
            bucket['control_aborts'] += t['welded_control'].get('exception') == 'ProcessAborted'
            bucket['repaired_aborts'] += t['repaired'].get('exception') == 'ProcessAborted'
        return dict(sorted(out.items()))

    faces = sorted(t['faces_repaired'] for t in trials)
    cuts = [faces[min(len(faces) - 1, len(faces) * i // 10)] for i in range(1, 10)]
    for t in trials:
        t['_decile'] = sum(t['faces_repaired'] > c for c in cuts)
    decile = table('_decile')
    thin = [t for t in trials if (t['radius_proxy_m'] or 1.0) < 1e-3]
    for t in trials:
        t.pop('_decile')
    return dict(
        schema='ihm.entity-tet-ready-verification.v1', build=str(build),
        entities_rediagnosed=len(rediagnosed), diagnostic_mismatches=len(mismatches),
        recorded_diagnostics_reproduced=not mismatches,
        tetgen_flags=flags, exhaustive=run_all, tetgen_trials=len(trials),
        repaired_tetrahedralized=len(successes),
        welded_control_tetrahedralized=sum(t['welded_control']['succeeded'] for t in trials),
        repaired_success_rate=len(successes) / len(trials) if trials else 0.0,
        welded_control_success_rate=sum(t['welded_control']['succeeded'] for t in trials) / len(trials)
        if trials else 0.0,
        by_role=table('role'), by_face_count_decile=decile, face_count_decile_cuts=cuts,
        thin_structure_trials=len(thin),
        thin_structure_tetrahedralized=sum(t['repaired']['succeeded'] for t in thin),
        thin_structure_control_tetrahedralized=sum(t['welded_control']['succeeded'] for t in thin),
        cavity_flagged_trials=sum(t['cavity_convention_pending'] for t in trials),
        cavity_flagged_tetrahedralized=sum(t['repaired']['succeeded'] and t['cavity_convention_pending']
                                           for t in trials),
        trials_with_prior_blockers=sum(bool(t['prior_blocking_reasons']) for t in trials),
        prior_blocked_tetrahedralized=sum(t['repaired']['succeeded'] and bool(t['prior_blocking_reasons'])
                                          for t in trials),
        trials_with_nonmanifold_vertices=sum(t['nonmanifold_vertices_after'] > 0 for t in trials),
        nonmanifold_vertex_trials_tetrahedralized=sum(t['repaired']['succeeded'] and
                                                      t['nonmanifold_vertices_after'] > 0 for t in trials),
        total_tets=sum(t['repaired'].get('tets', 0) for t in successes),
        median_tets=float(np.median([t['repaired']['tets'] for t in successes])) if successes else 0.0,
        max_absolute_tet_vs_surface_volume_error=max(errors, default=0.0),
        median_absolute_tet_vs_surface_volume_error=float(np.median(errors)) if errors else 0.0,
        max_absolute_repair_volume_drift=max(drifts, default=0.0),
        median_absolute_repair_volume_drift=float(np.median(drifts)) if drifts else 0.0,
        entities_with_volume_drift_over_1pct=[dict(entity_id=r['entity_id'], role=r['role'],
                                                   relative_drift=r['signed_volume_relative_drift'])
                                              for r in records
                                              if r['signed_volume_relative_drift'] is not None
                                              and abs(r['signed_volume_relative_drift']) > .01],
        repaired_process_aborts=sum(t['repaired'].get('exception') == 'ProcessAborted' for t in trials),
        welded_control_process_aborts=sum(t['welded_control'].get('exception') == 'ProcessAborted'
                                          for t in trials),
        failures=[dict(entity_id=t['entity_id'], name=t['name'], role=t['role'],
                       faces=t['faces_repaired'], radius_proxy_m=t['radius_proxy_m'],
                       nonmanifold_vertices=t['nonmanifold_vertices_after'],
                       nonmanifold_edges_after=t['nonmanifold_edges_after'],
                       boundary_edges_after=t['boundary_edges_after'],
                       exception=t['repaired'].get('exception'), message=t['repaired'].get('message'),
                       tetgen_stdout=t['repaired']['tetgen_stdout'][-12:])
                  for t in trials if not t['repaired']['succeeded']],
        elapsed_s=elapsed)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--build', type=Path, default=ROOT / 'data/derived/entity-tet-ready-v1')
    p.add_argument('--output', type=Path)
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--all', action='store_true')
    p.add_argument('--per-decile', type=int, default=2)
    p.add_argument('--flags', default='pYq1.414')
    a = p.parse_args()
    if a.self_test:
        self_test()
    if a.output:
        verify(a.build, a.output, a.workers, a.flags, a.all, a.per_decile)
    if not a.self_test and not a.output:
        p.error('Select --self-test or --output')
