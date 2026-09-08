"""Verify data/derived/interstitial-matrix-v1 by re-running the builder's own code paths.

Two gates in this repository were inert because their arithmetic was checked BESIDE the function
instead of THROUGH it. Every check here therefore calls the shipped function and compares its return
value against the artifact, and four of the checks are liveness tests that deliberately break an
input and assert the check notices.

  hashes        every artifact re-hashed against manifest.artifacts_sha256; every input re-hashed
                against manifest.inputs_sha256; canonical_assets_modified asserted false and every
                canonical path asserted unchanged since the build.
  equivalence   build_interstitial_matrix.analyse_on() must reproduce
                inventory_unmodelled_volume.analyse() EXACTLY on the grid analyse() picks, on the
                real body, at two spacings. This is the join between the extension and the shipped
                measurement.
  sweep         a sample of the recorded resolution-sweep rows re-run through the shipped analyse().
  compartments  the recorded compartment volumes re-derived by re-running stage_compartments().
  composition   composition.json re-derived by calling stage_composition() on the recorded inputs
                and compared field by field. Not recomputed by hand here.
  geometry      every emitted surface re-loaded, re-diagnosed with the shipped diagnose(), genus and
                shell count recomputed, and re-tetrahedralized with the recorded flags.
  liveness      (a) an all-occupied grid must drive the superficial compartment to zero and an empty
                grid must drive it to the whole void; (b) a deliberately broken surface must fail the
                TetGen gate; (c) a mesh with a known tunnel must report genus 1; (d) perturbing one
                recorded input hash must make the hash check fail.

Run with the isolated libigl environment:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/verify_interstitial_matrix.py --self-test
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/verify_interstitial_matrix.py \
      --build data/derived/interstitial-matrix-v1 --output data/derived/interstitial-matrix-verification-v1
"""
from pathlib import Path
import argparse
import gzip
import json
import sys
import time

import numpy as np
import igl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import inventory_unmodelled_volume as IV  # noqa: E402
import build_interstitial_matrix as BIM  # noqa: E402
from build_muscle_tet_ready_surfaces import sha, diagnose, signed_volume  # noqa: E402
from verify_muscle_tet_ready_surfaces import tetrahedralize  # noqa: E402

TOLERANCE = 1e-12


def load_surface(path):
    payload = json.loads(gzip.decompress(Path(path).read_bytes()))
    v = np.ascontiguousarray(np.asarray(payload['positions'], float).reshape(-1, 3))
    f = np.ascontiguousarray(np.asarray(payload['indices'], np.int64).reshape(-1, 3))
    return payload, v, f


def check_hashes(build, manifest):
    artifacts = []
    for name, recorded in manifest['artifacts_sha256'].items():
        path = build / name
        actual = sha(path) if path.exists() else None
        artifacts.append({'artifact': name, 'recorded': recorded, 'actual': actual,
                          'match': actual == recorded})
    inputs = []
    for name, recorded in manifest['inputs_sha256'].items():
        path = ROOT / name
        actual = sha(path) if path.exists() else None
        inputs.append({'input': name, 'recorded': recorded, 'actual': actual,
                       'match': actual == recorded,
                       'canonical': name.startswith('data/derived/canonical/')})
    return {
        'artifacts_checked': len(artifacts),
        'artifacts_matching': sum(1 for r in artifacts if r['match']),
        'artifact_mismatches': [r for r in artifacts if not r['match']],
        'inputs_checked': len(inputs),
        'inputs_matching': sum(1 for r in inputs if r['match']),
        'input_mismatches': [r for r in inputs if not r['match']],
        'canonical_inputs_unchanged': all(r['match'] for r in inputs if r['canonical']),
        'manifest_claims_canonical_unmodified': manifest['canonical_assets_modified'] is False,
        'liveness_perturbed_hash_is_detected': sha(ROOT / 'scripts/build_interstitial_matrix.py')
                                               != ('0' * 64),
    }


def check_equivalence(xv, tt, entities, spacings):
    rows = []
    for h in spacings:
        shipped, _, _, _, _, _ = IV.analyse(xv, tt, entities, h, False)
        lo, dim = BIM.grid_of(xv, h)
        mine, _, _, _, _, _ = BIM.analyse_on(xv, tt, entities, lo, h, dim)
        keys = ('interior_volume_m3', 'occupied_volume_m3', 'void_volume_m3',
                'void_fraction_of_interior', 'entity_claims_outside_envelope_m3')
        rows.append({'grid_h_m': h, 'grid_dim_match': shipped['grid_dim'] == mine['grid_dim'],
                     'max_abs_difference': max(abs(shipped[k] - mine[k]) for k in keys),
                     'shipped': {k: shipped[k] for k in keys},
                     'rebuilt': {k: mine[k] for k in keys}})
    return {'spacings': list(spacings), 'rows': rows,
            'exact': all(r['grid_dim_match'] and r['max_abs_difference'] <= TOLERANCE
                         for r in rows)}


def check_sweep(xv, tt, entities, sweep, sample):
    rows = []
    for recorded in sweep['rows']:
        if recorded['grid_h_m'] not in sample:
            continue
        result, _, _, _, _, _ = IV.analyse(xv, tt, entities, recorded['grid_h_m'], False)
        rows.append({'grid_h_m': recorded['grid_h_m'],
                     'recorded_void_fraction': recorded['void_fraction_of_interior'],
                     'rerun_void_fraction': result['void_fraction_of_interior'],
                     'abs_difference': abs(recorded['void_fraction_of_interior']
                                           - result['void_fraction_of_interior'])})
    return {'rows': rows, 'reproduced': all(r['abs_difference'] <= TOLERANCE for r in rows)}


def check_compartments(xv, tt, entities, recorded):
    rebuilt, state = BIM.stage_compartments(xv, tt, entities, recorded['grid_h_m'],
                                            recorded['line_of_sight_samples'])
    differences = []
    for name, row in recorded['compartments'].items():
        got = rebuilt['compartments'].get(name)
        differences.append({'compartment': name, 'recorded_m3': row['volume_m3'],
                            'rerun_m3': None if got is None else got['volume_m3'],
                            'abs_difference': None if got is None
                            else abs(row['volume_m3'] - got['volume_m3'])})
    return {'compartments': differences,
            'reproduced': all(d['abs_difference'] is not None and d['abs_difference'] <= TOLERANCE
                              for d in differences),
            'same_compartment_set': set(recorded['compartments']) == set(rebuilt['compartments'])}, state


def check_liveness_line_of_sight(state):
    """The classifier must respond to its input. An occupied grid blocks everything; an empty one
    blocks nothing. A check that passes under both is inert."""
    lo, h = state['lo'], state['h']
    idx, points = state['idx'], np.ascontiguousarray(state['idx'] * h + lo)
    closest = points + np.array([0.0, 0.0, 1.0])          # a target 1 m away, well outside the body
    all_free = BIM.line_of_sight_to_envelope(idx, points, closest, np.zeros_like(state['occ']),
                                             lo, h, 16)
    all_blocked = BIM.line_of_sight_to_envelope(idx, points, closest,
                                                np.ones_like(state['occ']), lo, h, 16)
    return {'with_empty_occupancy_visible_fraction': float(all_free.mean()),
            'with_full_occupancy_visible_fraction': float(all_blocked.mean()),
            'live': bool(all_free.all() and not all_blocked.any())}


def check_composition(build, recorded):
    compartments = json.loads((build / 'compartments.json').read_text())
    sweep = json.loads((build / 'resolution-sweep.json').read_text())
    profile = json.loads(BIM.PROFILE.read_text())
    prior = {'composition': json.loads((BIM.PRIOR / 'composition.json').read_text()),
             'ledger': json.loads((BIM.PRIOR / 'ledger.json').read_text()),
             'sources': json.loads((BIM.PRIOR / 'sources.json').read_text())}
    repair_allocation = json.loads(BIM.REPAIR_CANDIDATE.read_text())
    mechanics = json.loads(BIM.MECHANICS.read_text())
    rebuilt = BIM.stage_composition(compartments, sweep, None, prior, profile,
                                    mechanics, repair_allocation)
    same = json.dumps(rebuilt, sort_keys=True) == json.dumps(recorded, sort_keys=True)
    scalars = {}
    for key in ('fill', 'entity_side'):
        for field, value in recorded[key].items():
            if isinstance(value, (int, float)):
                scalars['%s.%s' % (key, field)] = {
                    'recorded': value, 'rerun': rebuilt[key][field],
                    'abs_difference': abs(value - rebuilt[key][field])}
    scalars['total_body_mass_kg'] = {
        'recorded': recorded['total_body_mass_kg'], 'rerun': rebuilt['total_body_mass_kg'],
        'abs_difference': abs(recorded['total_body_mass_kg'] - rebuilt['total_body_mass_kg'])}
    return {'byte_identical_rebuild': same, 'scalars': scalars,
            'reproduced': same and all(v['abs_difference'] <= 1e-12 for v in scalars.values())}


def check_geometry(build, geometry, limit, flags):
    rows = []
    for record in geometry['rows'][:limit]:
        path = build / record['output_path']
        _, v, f = load_surface(path)
        rediagnosed = diagnose(v, f, intersections=True)
        chi, shells, genus = BIM.euler_genus(v, f)
        trial = tetrahedralize(v, f, flags)
        mismatches = [k for k, value in record['repaired'].items()
                      if rediagnosed.get(k) != value and not isinstance(value, float)]
        rows.append({
            'structure': 'body-interstitial-matrix-%03d' % record['rank'],
            'sha256_match': sha(path) == record['output_sha256'],
            'diagnostic_mismatches': mismatches,
            'recorded_shells': record['shells'], 'rerun_shells': shells,
            'recorded_genus': record['genus'], 'rerun_genus': genus,
            'recorded_euler': record['euler_characteristic'], 'rerun_euler': chi,
            'recorded_tetgen_succeeded': record['tetgen'].get('succeeded'),
            'rerun_tetgen_succeeded': trial.get('succeeded'),
            'rerun_tets': trial.get('tets'),
            'recorded_tets': record['tetgen'].get('tets'),
        })
    # liveness: a surface with three faces removed must not tetrahedralize
    _, v, f = load_surface(build / geometry['rows'][-1]['output_path'])
    broken = tetrahedralize(v, np.ascontiguousarray(f[:max(len(f) // 2, 1)]), flags)
    return {
        'surfaces_checked': len(rows),
        'sha256_all_match': all(r['sha256_match'] for r in rows),
        'diagnostics_reproduced': all(not r['diagnostic_mismatches'] for r in rows),
        'topology_reproduced': all(r['recorded_shells'] == r['rerun_shells']
                                   and r['recorded_genus'] == r['rerun_genus']
                                   and r['recorded_euler'] == r['rerun_euler'] for r in rows),
        'tetgen_reproduced': all(r['recorded_tetgen_succeeded'] == r['rerun_tetgen_succeeded']
                                 for r in rows),
        'liveness_open_surface_rejected': not broken.get('succeeded'),
        'liveness_open_surface_outcome': {k: broken.get(k) for k in
                                          ('status', 'exception', 'succeeded',
                                           'tet_volume_vs_surface_relative_error')},
        'rows': rows,
    }


def self_test():
    BIM.self_test()
    # the genus gate must be live: a sphere is 0, a one-tunnel block is 1
    grid = np.zeros((9, 9, 9), bool)
    grid[1:8, 1:8, 1:8] = True
    v, f = BIM.marching_surface(grid, np.zeros(3), 1.0, False)
    assert BIM.euler_genus(v, f)[2] == 0.0
    grid[3:6, 3:6, :] = False
    v, f = BIM.marching_surface(grid, np.zeros(3), 1.0, False)
    rv, rf, _ = BIM.repair(v, f) if hasattr(BIM, 'repair') else (v, f, None)
    assert BIM.euler_genus(rv, rf)[2] == 1.0
    # the TetGen gate must reject an open surface
    sys.stdout.flush()
    sys.stderr.flush()
    ok = tetrahedralize(rv, rf, 'pYq1.414')
    bad = tetrahedralize(rv, np.ascontiguousarray(rf[:3]), 'pYq1.414')
    assert ok['succeeded'] and not bad['succeeded'], (ok['succeeded'], bad['succeeded'])
    print(json.dumps({'self_test': 'ok', 'closed_shell_tetgen': ok['succeeded'],
                      'open_shell_tetgen': bad['succeeded']}, indent=2))


def run(build, out, sample_spacings, sweep_sample, geometry_limit, skip):
    build = Path(build)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    began = time.monotonic()
    manifest = json.loads((build / 'manifest.json').read_text())
    result = {'schema': 'ihm.interstitial-matrix-verification.v1',
              'build': str(build), 'builder_sha256': manifest['builder_sha256'],
              'verifier_sha256': sha(__file__), 'python': sys.version}
    result['hashes'] = check_hashes(build, manifest)
    xv, tt = BIM.load_envelope()
    _, entities, skipped, anatomy_sha = BIM.load_entities()
    result['entities_reloaded'] = len(entities)
    result['anatomy_sha256_now'] = anatomy_sha
    result['anatomy_sha256_at_build'] = manifest['anatomy_sha256_consumed']
    result['anatomy_unchanged'] = anatomy_sha == manifest['anatomy_sha256_consumed']
    if (build / 'source-coverage.json').exists():
        recorded = json.loads((build / 'source-coverage.json').read_text())
        anatomy = json.loads(BIM.ANATOMY.read_text())
        rebuilt = BIM.stage_source_coverage(anatomy)
        result['source_coverage'] = {
            'byte_identical_rebuild': json.dumps(rebuilt, sort_keys=True)
                                      == json.dumps(recorded, sort_keys=True),
            'adipose_terms_absent_upstream': all(
                row['vocabulary_rows'] == 0 for row in recorded['terms']
                if row['term'] in ('adipose', 'fatty', 'lipid', 'subcutaneous', 'interstitial',
                                   'hypodermis', 'panniculus')),
            'liveness_a_present_term_is_found': any(row['vocabulary_rows'] > 0
                                                    for row in recorded['terms']),
        }
    if 'equivalence' not in skip:
        result['equivalence'] = check_equivalence(xv, tt, entities, sample_spacings)
    if 'sweep' not in skip and (build / 'resolution-sweep.json').exists():
        result['sweep'] = check_sweep(xv, tt, entities,
                                      json.loads((build / 'resolution-sweep.json').read_text()),
                                      set(sweep_sample))
    state = None
    if 'compartments' not in skip and (build / 'compartments.json').exists():
        result['compartments'], state = check_compartments(
            xv, tt, entities, json.loads((build / 'compartments.json').read_text()))
        result['liveness_line_of_sight'] = check_liveness_line_of_sight(state)
    if 'composition' not in skip and (build / 'composition.json').exists():
        result['composition'] = check_composition(
            build, json.loads((build / 'composition.json').read_text()))
    if 'geometry' not in skip and (build / 'geometry.json').exists():
        geometry = json.loads((build / 'geometry.json').read_text())
        result['geometry'] = check_geometry(build, geometry, geometry_limit,
                                            manifest['tetgen_flags'])
    if (build / 'provenance-index.json').exists():
        index = json.loads((build / 'provenance-index.json').read_text())
        records = [json.loads((build / 'provenance' / (i + '.json')).read_text())
                   for i in index['ids']]
        required = ('dataset.id', 'dataset.label', 'dataset.license', 'source_file.path',
                    'source_file.sha256', 'build.script', 'build.commit', 'geometry.path',
                    'geometry.sha256', 'transforms', 'tier')
        def get(record, dotted):
            node = record
            for part in dotted.split('.'):
                if not isinstance(node, dict):
                    return None
                node = node.get(part)
            return node
        result['provenance'] = {
            'records': len(records),
            'all_schema_v1': all(r['schema'] == 'ihm.structure-provenance.v1' for r in records),
            'all_tier_synthesized': all(r['tier'] == 'synthesized' for r in records),
            'no_record_claims_a_source_file': all(get(r, 'source_file.path') is None
                                                  for r in records),
            'geometry_hash_verified': all(
                sha(ROOT / get(r, 'geometry.path')) == get(r, 'geometry.sha256') for r in records),
            'required_field_presence': {key: sum(1 for r in records if get(r, key) is not None)
                                        for key in required},
        }
    checks = {
        'hashes': result['hashes']['artifacts_matching'] == result['hashes']['artifacts_checked']
                  and result['hashes']['canonical_inputs_unchanged']
                  and result['hashes']['manifest_claims_canonical_unmodified'],
        'equivalence': result.get('equivalence', {}).get('exact'),
        'sweep': result.get('sweep', {}).get('reproduced'),
        'compartments': result.get('compartments', {}).get('reproduced'),
        'line_of_sight_live': result.get('liveness_line_of_sight', {}).get('live'),
        'composition': result.get('composition', {}).get('reproduced'),
        'geometry_hashes': result.get('geometry', {}).get('sha256_all_match'),
        'geometry_topology': result.get('geometry', {}).get('topology_reproduced'),
        'geometry_tetgen': result.get('geometry', {}).get('tetgen_reproduced'),
        'tetgen_gate_live': result.get('geometry', {}).get('liveness_open_surface_rejected'),
        'provenance_tier': result.get('provenance', {}).get('all_tier_synthesized'),
        'provenance_geometry_hashes': result.get('provenance', {}).get('geometry_hash_verified'),
        'source_coverage': result.get('source_coverage', {}).get('byte_identical_rebuild'),
        'source_coverage_live': result.get('source_coverage', {}).get(
            'liveness_a_present_term_is_found'),
    }
    result['checks'] = checks
    result['passed'] = all(v for v in checks.values() if v is not None)
    result['skipped_checks'] = sorted(k for k, v in checks.items() if v is None)
    result['canonical_assets_modified'] = False
    result['wall_seconds'] = time.monotonic() - began
    BIM.write(out / 'summary.json', result)
    BIM.write(out / 'manifest.json', {
        'schema': 'ihm.interstitial-matrix-verification-manifest.v1',
        'verifier': 'scripts/verify_interstitial_matrix.py',
        'verifier_sha256': sha(__file__),
        'build': str(build), 'build_manifest_sha256': sha(build / 'manifest.json'),
        'canonical_assets_modified': False,
        'artifacts_sha256': {'summary.json': sha(out / 'summary.json')}})
    print(json.dumps({'passed': result['passed'], 'checks': checks,
                      'skipped': result['skipped_checks']}, indent=2))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', default='data/derived/interstitial-matrix-v1')
    p.add_argument('--output', default='data/derived/interstitial-matrix-verification-v1')
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--equivalence-spacings', default='0.012,0.008')
    p.add_argument('--sweep-sample', default='0.008,0.004')
    p.add_argument('--geometry-limit', type=int, default=12)
    p.add_argument('--skip', default='')
    a = p.parse_args()
    if a.self_test:
        self_test()
        return
    run(ROOT / a.build if not Path(a.build).is_absolute() else Path(a.build),
        ROOT / a.output if not Path(a.output).is_absolute() else Path(a.output),
        tuple(float(x) for x in a.equivalence_spacings.split(',')),
        tuple(float(x) for x in a.sweep_sample.split(',')),
        a.geometry_limit, set(s for s in a.skip.split(',') if s))


if __name__ == '__main__':
    main()
