"""Conservative repair of the non-muscular canonical surfaces into tetrahedralisation-ready shells.

The repair order is imported verbatim from scripts/build_muscle_tet_ready_surfaces.py so this lane
cannot drift from the muscular lane that it extends. Nothing in the muscular module is modified.

Run with the isolated libigl environment, not the physiological runtime:
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_entity_tet_ready_surfaces.py --self-test
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_entity_tet_ready_surfaces.py \
  --output data/derived/entity-tet-ready-v1 --workers 10

What this adds over the muscular lane, and why:
  * scope             every canonical entity whose role is not `muscle`. Entities whose reference
                      geometry is not a triangular surface (3 skin_layer quadrature shells, 1
                      lymphatic structural graph) and the whole-body skin envelope (already repaired
                      at data/derived/outer-envelope) are excluded by name with a recorded reason.
  * nesting analysis  step 5 of the muscular recipe flips EVERY orientable patch to positive signed
                      volume. On muscle that only ever affected slivers. On organs a patch nested
                      inside another patch can be a genuine anatomical cavity, and flipping it turns
                      a void into solid material. analyse_patches() is therefore run on the welded
                      input and on the repaired output; it labels each patch with its containment
                      depth (winding number of the candidate container evaluated at vertices of the
                      candidate content) and the signed volume the SOURCE authored for it. Odd depth
                      plus authored-negative volume is an authored cavity. Nothing is decided here:
                      the geometry is written with the muscular convention (all patches outward) and
                      every affected entity is flagged cavity_convention_pending, with a seed point
                      inside each void recorded so a mesher can be handed a hole marker later.
  * thin-structure    vascular and nerve surfaces are tubes. radius_proxy_m = 2*volume/area and the
                      minimum positive face area are recorded so near-degenerate tubes are visible.
No step fills a hole or moves an existing coordinate.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_muscle_tet_ready_surfaces import (  # noqa: E402  the muscular recipe is the contract
    sha, write, signed_volume, edge_incidence, weld_exact, diagnose, normalize, repair,
    self_intersecting_pairs)

# Repaired elsewhere and deliberately not re-derived; see data/derived/outer-envelope.
SKIN_ENVELOPE = 'body-bp3d-FJ2810'
EXCLUSION_REASONS = {
    'not_a_triangular_surface': 'reference_geometry.representation is not triangular_surface',
    'outer_envelope_already_repaired': 'whole-body skin slab; watertight envelope exists at '
                                       'data/derived/outer-envelope',
}


def surface_area(v, f):
    if not len(f):
        return 0.0
    n = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    return float(np.linalg.norm(n, axis=1).sum() / 2.0)


def shape_metrics(v, f):
    """Scale descriptors that separate tubes and sheets from blobs."""
    if not len(f):
        return dict(bbox_diagonal_m=0.0, surface_area_m2=0.0, radius_proxy_m=None,
                    min_edge_length_m=None, median_edge_length_m=None,
                    min_positive_face_area_m2=None)
    area = surface_area(v, f)
    volume = signed_volume(v, f)
    e = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    lengths = np.linalg.norm(v[e[:, 0]] - v[e[:, 1]], axis=1)
    cross = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    areas = np.linalg.norm(cross, axis=1) / 2.0
    positive = areas[areas > 0]
    return dict(
        bbox_diagonal_m=float(np.linalg.norm(v.max(axis=0) - v.min(axis=0))),
        surface_area_m2=area,
        # For a long cylinder of radius r, 2*V/A -> r. For a slab of thickness t, 2*V/A -> t.
        radius_proxy_m=float(2.0 * volume / area) if area > 0 else None,
        min_edge_length_m=float(lengths.min()),
        median_edge_length_m=float(np.median(lengths)),
        min_positive_face_area_m2=float(positive.min()) if len(positive) else None)


def patch_labels(f):
    """Connected patches under bfs_orient, with each patch consistently wound."""
    oriented, patch = igl.bfs_orient(f)
    return np.ascontiguousarray(np.asarray(oriented, np.int64)), np.asarray(patch).ravel()


def _inside(container_v, container_f, points):
    """Fraction of `points` strictly enclosed by the closed patch (container_v, container_f)."""
    w = np.asarray(igl.winding_number(np.ascontiguousarray(container_v),
                                      np.ascontiguousarray(container_f),
                                      np.ascontiguousarray(points)), float)
    return float((np.abs(w) > 0.5).mean())


def cavity_seed_point(v, f, samples=12):
    """A point strictly inside the closed patch, for use as a TetGen hole marker.

    A coarse grid over the patch bounding box is tested with the winding number and the enclosed
    candidate furthest from the surface is returned. None when the patch encloses nothing sampled.
    """
    lo, hi = v.min(axis=0), v.max(axis=0)
    grid = np.stack(np.meshgrid(*[np.linspace(a, b, samples)[1:-1] for a, b in zip(lo, hi)],
                                indexing='ij'), axis=-1).reshape(-1, 3)
    if not len(grid):
        return None
    w = np.asarray(igl.winding_number(v, f, np.ascontiguousarray(grid)), float)
    inside = grid[np.abs(w) > 0.5]
    if not len(inside):
        return None
    d = np.abs(np.asarray(igl.signed_distance(np.ascontiguousarray(inside), v, f)[0], float))
    return [float(x) for x in inside[int(np.argmax(d))]]


def analyse_patches(v, f, seed_points=False, max_patches=400, sample_points=9):
    """Per-patch signed volume and containment depth on a surface that is NOT reoriented here.

    Containment is decided by evaluating the winding number of the candidate container at up to
    `sample_points` vertices of the candidate content; the patches are disjoint after the repair,
    so any vertex of one lies strictly inside or strictly outside the other. Bounding-box
    containment prefilters the pairs. Depth is the number of patches that contain a patch, so an
    odd depth is a void under the standard even-odd rule.
    """
    if not len(f):
        return dict(patches=0, analysed=False, reason='empty')
    oriented, label = patch_labels(f)
    count = int(label.max()) + 1
    records = []
    for i in range(count):
        pf = oriented[label == i]
        used = np.unique(pf)
        boundary, nonmanifold = edge_incidence(pf)
        records.append(dict(patch=i, faces=int(len(pf)), vertices=int(len(used)),
                            signed_volume_m3=signed_volume(v, pf), closed=boundary == 0 and nonmanifold == 0,
                            bbox_min=v[used].min(axis=0), bbox_max=v[used].max(axis=0),
                            _faces=pf, _used=used))
    if count > max_patches:
        for r in records:
            r.pop('_faces'), r.pop('_used'), r.pop('bbox_min'), r.pop('bbox_max')
        return dict(patches=count, analysed=False, reason='patch_count_over_%d' % max_patches,
                    patch_signed_volumes=[r['signed_volume_m3'] for r in records])
    containment = {i: [] for i in range(count)}
    for a in range(count):
        ra = records[a]
        if not ra['closed'] or not ra['faces']:
            continue
        pts = v[ra['_used']]
        if len(pts) > sample_points:
            pts = pts[np.linspace(0, len(pts) - 1, sample_points).astype(int)]
        for b in range(count):
            if a == b or not records[b]['closed']:
                continue
            rb = records[b]
            if not (np.all(ra['bbox_min'] >= rb['bbox_min']) and np.all(ra['bbox_max'] <= rb['bbox_max'])):
                continue
            if _inside(v, rb['_faces'], pts) > 0.5:
                containment[a].append(b)
    out = []
    for i, r in enumerate(records):
        depth = len(containment[i])
        entry = dict(patch=i, faces=r['faces'], signed_volume_m3=r['signed_volume_m3'],
                     closed=r['closed'], contained_in=containment[i], depth=depth,
                     authored_negative=r['signed_volume_m3'] < 0,
                     void_by_parity=depth % 2 == 1)
        if seed_points and depth % 2 == 1 and r['closed']:
            entry['seed_point_m'] = cavity_seed_point(v, r['_faces'])
        out.append(entry)
    nested = [r for r in out if r['depth'] > 0]
    return dict(patches=count, analysed=True,
                nested_patches=len(nested),
                void_patches=[r['patch'] for r in out if r['void_by_parity']],
                authored_cavities=[r['patch'] for r in out if r['void_by_parity'] and r['authored_negative']],
                authored_negative_patches=[r['patch'] for r in out if r['authored_negative']],
                enclosed_volume_m3=float(sum(abs(r['signed_volume_m3']) for r in out if r['void_by_parity'])),
                detail=out)


def load_entity(entity):
    reference = entity['reference_geometry']
    path = ROOT / reference['path']
    digest = sha(path)
    if digest != reference['sha256']:
        raise ValueError('Canonical geometry changed: ' + entity['id'])
    raw = json.loads(gzip.decompress(path.read_bytes()))
    v = np.asarray(raw['positions'], float).reshape(-1, 3)
    f = np.ascontiguousarray(np.asarray(raw['indices'], np.int64).reshape(-1, 3))
    if raw.get('display_decimation'):
        raise ValueError('Decimated display geometry rejected: ' + entity['id'])
    if not np.isfinite(v).all() or not len(f) or f.min() < 0 or f.max() >= len(v):
        raise ValueError('Invalid source geometry: ' + entity['id'])
    return digest, v, f


def process(job):
    entity, prior, out = job
    ident = entity['id']
    try:
        digest, v, f = load_entity(entity)
        raw = diagnose(v, f, intersections=False)
        wv, wf = weld_exact(v, f)
        before = diagnose(wv, wf)
        before_patches = analyse_patches(wv, wf)
        rv, rf, log = repair(v, f)
        after = diagnose(rv, rf)
        after_patches = analyse_patches(rv, rf, seed_points=True)
        drift = after['signed_volume_m3'] - before['signed_volume_m3']
        payload = dict(schema='ihm.entity-tet-ready-surface.v1', entity_id=ident, name=entity['name'],
                       role=entity['role'], system=entity['system'], units='m',
                       frame=entity['reference_geometry']['frame'],
                       representation='triangular_surface',
                       source_path=entity['reference_geometry']['path'], source_sha256=digest,
                       positions=[float(x) for x in rv.ravel()], indices=[int(i) for i in rf.ravel()])
        geometry = Path(out) / 'geometry' / (ident + '.json.gz')
        geometry.write_bytes(gzip.compress(json.dumps(payload, allow_nan=False).encode(), mtime=0))
        cavity_suspect = bool(after_patches.get('void_patches') or before_patches.get('authored_cavities')
                              or (not after_patches.get('analysed', True)))
        record = dict(entity_id=ident, name=entity['name'], role=entity['role'], system=entity['system'],
                      source_path=entity['reference_geometry']['path'], source_sha256=digest,
                      prior_blocking_reasons=prior.get('blocking_reasons'),
                      prior_topological_candidate=prior.get('topological_candidate'),
                      raw=raw, before=before, after=after, operations=log,
                      before_patch_analysis=before_patches, after_patch_analysis=after_patches,
                      shape_before=shape_metrics(wv, wf), shape_after=shape_metrics(rv, rf),
                      output_path=str(geometry.relative_to(Path(out))), output_sha256=sha(geometry),
                      signed_volume_drift_m3=drift,
                      signed_volume_relative_drift=drift / before['signed_volume_m3']
                      if before['signed_volume_m3'] else None,
                      hole_filling_required=after['boundary_edges'] > 0,
                      vertex_links_manifold=after['nonmanifold_vertices'] == 0,
                      cavity_convention_pending=cavity_suspect,
                      tet_ready=after['closed'] and after['self_intersection_free'] and
                      after['orientation_consistent'] and after['signed_volume_m3'] > 0)
        return record
    except BaseException as error:  # a single pathological surface must not end the sweep
        return dict(entity_id=ident, name=entity['name'], role=entity['role'], system=entity['system'],
                    failed=True, error_type=type(error).__name__, error=str(error))


def self_test():
    tet_v = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [0, 0, 1.]])
    tet_f = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]], np.int64)
    solo = analyse_patches(tet_v, tet_f)
    assert solo['patches'] == 1 and solo['nested_patches'] == 0 and not solo['void_patches']
    # A cube with an inward-wound cube inside it: the classic authored cavity.
    def cube(scale, offset):
        c = np.array([[0., 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                      [0, 0, 1.], [1, 0, 1], [1, 1, 1], [0, 1, 1]]) * scale + offset
        q = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4],
                      [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]], np.int64)
        return c, q
    ov, of = cube(1.0, np.zeros(3))
    iv, if_ = cube(0.4, np.full(3, 0.3))
    assert signed_volume(ov, of) > 0
    hollow_v = np.r_[ov, iv]
    hollow_f = np.r_[of, (if_ + len(ov))[:, [0, 2, 1]]]  # inner shell wound inward = a real void
    analysis = analyse_patches(hollow_v, hollow_f, seed_points=True)
    assert analysis['patches'] == 2, analysis
    assert analysis['nested_patches'] == 1 and analysis['void_patches'] == [1], analysis
    assert analysis['authored_cavities'] == [1], analysis
    inner = analysis['detail'][1]
    assert inner['contained_in'] == [0] and inner['authored_negative']
    seed = inner['seed_point_m']
    assert seed is not None and all(0.3 < c < 0.7 for c in seed), seed
    assert abs(analysis['enclosed_volume_m3'] - 0.4 ** 3) < 1e-12
    # The muscular recipe flips the void solid; the analysis of its output must say so.
    rv, rf, _ = repair(hollow_v, hollow_f)
    after = analyse_patches(rv, rf)
    assert after['void_patches'] and not after['authored_cavities'], \
        'a flipped void must still be reported as a parity void with no surviving authored sign'
    # Disjoint, non-nested patches must not be reported as cavities.
    apart = analyse_patches(np.r_[ov, iv + 5], np.r_[of, if_ + len(ov)])
    assert apart['patches'] == 2 and apart['nested_patches'] == 0 and not apart['void_patches']
    # Shape metrics on a thin tube-like slab.
    m = shape_metrics(ov, of)
    assert abs(m['surface_area_m2'] - 6.0) < 1e-12 and abs(m['radius_proxy_m'] - 2.0 / 6.0) < 1e-12
    assert m['min_positive_face_area_m2'] > 0 and m['min_edge_length_m'] > 0
    assert cavity_seed_point(ov, of) is not None
    print('PASS nesting, authored-cavity, seed-point, flip-detection and shape-metric checks')


def build(output, roles, limit, workers, order):
    out = Path(output).resolve()
    if out.exists():
        raise ValueError('Choose a fresh output directory')
    anatomy_path = ROOT / 'data/derived/canonical/anatomy.json'
    readiness = ROOT / 'data/derived/surface-volume-readiness-v1'
    anatomy = json.loads(anatomy_path.read_text())
    prior = {r['entity_id']: r for r in map(json.loads, (readiness / 'entities.jsonl').read_text().splitlines())}
    candidates = [e for e in anatomy['entities'] if e['role'] != 'muscle']
    if roles:
        candidates = [e for e in candidates if e['role'] in roles]
    excluded = []
    entities = []
    for e in candidates:
        if e['reference_geometry']['representation'] != 'triangular_surface':
            excluded.append(dict(entity_id=e['id'], name=e['name'], role=e['role'],
                                 representation=e['reference_geometry']['representation'],
                                 reason='not_a_triangular_surface',
                                 detail=EXCLUSION_REASONS['not_a_triangular_surface']))
        elif e['id'] == SKIN_ENVELOPE:
            excluded.append(dict(entity_id=e['id'], name=e['name'], role=e['role'],
                                 representation=e['reference_geometry']['representation'],
                                 reason='outer_envelope_already_repaired',
                                 detail=EXCLUSION_REASONS['outer_envelope_already_repaired']))
        else:
            entities.append(e)
    if not entities:
        raise ValueError('No entities selected')
    if order == 'faces':
        entities.sort(key=lambda e: e['source_face_count'])
    if limit:
        entities = entities[:limit]
    packages = {n: importlib.metadata.version(n) for n in ('libigl', 'numpy')}
    if packages['libigl'] != '2.6.2':
        raise ValueError('This build pins libigl 2.6.2')
    out.mkdir(parents=True)
    (out / 'geometry').mkdir()
    (out / 'inputs').mkdir()
    shutil.copyfile(__file__, out / 'inputs/build_entity_tet_ready_surfaces.py')
    shutil.copyfile(ROOT / 'scripts/build_muscle_tet_ready_surfaces.py',
                    out / 'inputs/build_muscle_tet_ready_surfaces.py')
    shutil.copyfile(readiness / 'manifest.json', out / 'inputs/surface-volume-readiness-manifest.json')
    native = Path(cgal.pyigl_copyleft_cgal.__file__)
    shutil.copyfile(native, out / 'inputs' / native.name)
    write(out / 'excluded.json', excluded)
    jobs = [(e, prior.get(e['id'], {}), str(out)) for e in entities]
    started = time.monotonic()
    records = []
    try:
        if workers > 1:
            with mp.get_context('fork').Pool(workers, maxtasksperchild=40) as pool:
                for record in pool.imap_unordered(process, jobs, chunksize=4):
                    records.append(record)
                    if len(records) % 100 == 0:
                        print('Repaired', len(records), 'of', len(jobs),
                              '%.0fs' % (time.monotonic() - started), flush=True)
        else:
            for job in jobs:
                records.append(process(job))
        rank = {e['id']: i for i, e in enumerate(entities)}
        records.sort(key=lambda r: rank[r['entity_id']])
        with (out / 'entities.jsonl').open('w') as handle:
            for record in records:
                handle.write(json.dumps(record, allow_nan=False) + '\n')
        write(out / 'summary.json', summarise(records, excluded, time.monotonic() - started))
        sources = {r['source_path']: r['source_sha256'] for r in records if not r.get('failed')}
        artifacts = {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
        write(out / 'manifest.json', dict(
            schema='ihm.entity-tet-ready-build.v1', packages=packages, python=sys.version,
            native_module_sha256=sha(native), anatomy_sha256=sha(anatomy_path),
            muscle_builder_sha256=sha(ROOT / 'scripts/build_muscle_tet_ready_surfaces.py'),
            readiness_manifest_sha256=sha(readiness / 'manifest.json'),
            source_geometry_sha256=sources, artifacts_sha256=artifacts,
            scope='every canonical entity whose role is not muscle, minus the entities listed in '
                  'excluded.json',
            method='repair order imported unchanged from scripts/build_muscle_tet_ready_surfaces.py: '
                   'exact weld; drop repeated-index and exactly-collinear faces; collapse coincident '
                   'faces (opposite-winding pairs removed as zero-volume fins); bfs_orient patches '
                   'flipped to positive signed volume; CGAL remesh_self_intersections(stitch_all) '
                   'arrangement; CGAL self-union outer shell only where the arrangement leaves '
                   'nonmanifold edges. Added here: per-patch containment analysis before and after, '
                   'and thin-structure scale metrics.',
            canonical_geometry_modified=False, limitations=[
                'Step 5 of the imported recipe orients every patch outward. Entities whose patch '
                'analysis reports a parity void or an authored cavity are flagged '
                'cavity_convention_pending: their written geometry is the solid-inclusion reading, '
                'and a seed point inside each void is recorded so the opposite reading remains '
                'reachable without re-deriving the surface. No cavity convention is imposed here.',
                'The self-union step discards material enclosed twice by overlapping lobes of one '
                'surface; its volume change is recorded per entity.',
                'No hole is filled and no existing coordinate is moved; entities still carrying '
                'boundary edges are flagged hole_filling_required and are not certified.',
                'Containment is tested with the winding number at sampled vertices of the contained '
                'patch, and is not attempted for entities with more than 400 patches; those are '
                'reported as analysed=false rather than assumed cavity-free.',
                'Tetrahedralisation itself is proved separately by '
                'scripts/verify_entity_tet_ready_surfaces.py.']))
        print(json.dumps({k: v for k, v in json.loads((out / 'summary.json').read_text()).items()
                          if not isinstance(v, list)}, indent=2))
    except BaseException as error:
        write(out / 'failure.json', dict(type=type(error).__name__, message=str(error),
                                         completed=len(records)))
        raise


def summarise(records, excluded, elapsed):
    ok = [r for r in records if not r.get('failed')]
    failed = [r for r in records if r.get('failed')]

    def tally(key, side, subset=None):
        return sum(bool(r[side][key]) for r in (subset if subset is not None else ok))

    def by_role(predicate):
        roles = {}
        for r in ok:
            roles.setdefault(r['role'], [0, 0])
            roles[r['role']][1] += 1
            roles[r['role']][0] += bool(predicate(r))
        return {k: dict(matching=v[0], entities=v[1]) for k, v in sorted(roles.items())}

    cavity = [r for r in ok if r['cavity_convention_pending']]
    return dict(
        schema='ihm.entity-tet-ready-build.v1', entities=len(records), repaired=len(ok),
        build_failures=[dict(entity_id=r['entity_id'], role=r['role'], error_type=r['error_type'],
                             error=r['error']) for r in failed],
        excluded=excluded,
        stage_meaning=dict(raw='canonical rendering topology, seam-duplicated vertices, no intersection test',
                           before='after lossless exact-coordinate welding only; this is the defect baseline',
                           after='after the full documented repair order'),
        raw=dict(vertices=sum(r['raw']['vertices'] for r in ok), closed=sum(r['raw']['closed'] for r in ok),
                 total_boundary_edges=sum(r['raw']['boundary_edges'] for r in ok),
                 total_nonmanifold_edges=sum(r['raw']['nonmanifold_edges'] for r in ok)),
        faces_in=sum(r['before']['faces'] for r in ok), faces_out=sum(r['after']['faces'] for r in ok),
        before=dict(
            closed=tally('closed', 'before'), self_intersection_free=tally('self_intersection_free', 'before'),
            orientation_consistent=tally('orientation_consistent', 'before'),
            with_duplicate_faces=sum(r['before']['duplicate_face_copies'] > 0 for r in ok),
            with_zero_area_faces=sum(r['before']['zero_area_faces'] > 0 for r in ok),
            with_repeated_index_faces=sum(r['before']['repeated_index_faces'] > 0 for r in ok),
            with_nonmanifold_edges=sum(r['before']['nonmanifold_edges'] > 0 for r in ok),
            with_nonmanifold_vertices=sum(r['before']['nonmanifold_vertices'] > 0 for r in ok),
            with_boundary_edges=sum(r['before']['boundary_edges'] > 0 for r in ok),
            multi_component=sum(r['before']['face_components'] > 1 for r in ok),
            total_duplicate_face_copies=sum(r['before']['duplicate_face_copies'] for r in ok),
            total_zero_area_faces=sum(r['before']['zero_area_faces'] for r in ok),
            total_repeated_index_faces=sum(r['before']['repeated_index_faces'] for r in ok),
            total_nonmanifold_edges=sum(r['before']['nonmanifold_edges'] for r in ok),
            total_boundary_edges=sum(r['before']['boundary_edges'] for r in ok),
            total_self_intersecting_face_pairs=sum(r['before']['self_intersecting_face_pairs'] for r in ok)),
        after=dict(
            closed=tally('closed', 'after'), self_intersection_free=tally('self_intersection_free', 'after'),
            orientation_consistent=tally('orientation_consistent', 'after'),
            with_duplicate_faces=sum(r['after']['duplicate_face_copies'] > 0 for r in ok),
            with_zero_area_faces=sum(r['after']['zero_area_faces'] > 0 for r in ok),
            with_repeated_index_faces=sum(r['after']['repeated_index_faces'] > 0 for r in ok),
            with_nonmanifold_edges=sum(r['after']['nonmanifold_edges'] > 0 for r in ok),
            with_nonmanifold_vertices=sum(r['after']['nonmanifold_vertices'] > 0 for r in ok),
            with_boundary_edges=sum(r['after']['boundary_edges'] > 0 for r in ok),
            multi_component=sum(r['after']['face_components'] > 1 for r in ok),
            total_self_intersecting_face_pairs=sum(r['after']['self_intersecting_face_pairs'] for r in ok)),
        tet_ready=sum(r['tet_ready'] for r in ok),
        tet_ready_with_manifold_vertex_links=sum(r['tet_ready'] and r['vertex_links_manifold'] for r in ok),
        tet_ready_by_role=by_role(lambda r: r['tet_ready']),
        self_intersecting_before_by_role=by_role(lambda r: not r['before']['self_intersection_free']),
        fin_faces_before_by_role=by_role(lambda r: r['before']['duplicate_face_copies'] > 0),
        pinched_vertex_links=[r['entity_id'] for r in ok if not r['vertex_links_manifold']],
        hole_filling_required=[r['entity_id'] for r in ok if r['hole_filling_required']],
        self_union_applied=[r['entity_id'] for r in ok
                            if any(s['stage'] == 'extract_outer_manifold' for s in r['operations'])],
        arrangement_applied=sum(any(s['stage'] == 'resolve_self_intersections' for s in r['operations'])
                                for r in ok),
        cavity_convention_pending=[dict(entity_id=r['entity_id'], name=r['name'], role=r['role'],
                                        patches_after=r['after_patch_analysis'].get('patches'),
                                        void_patches_after=r['after_patch_analysis'].get('void_patches'),
                                        authored_cavities_before=r['before_patch_analysis'].get('authored_cavities'),
                                        enclosed_void_volume_m3=r['after_patch_analysis'].get('enclosed_volume_m3'),
                                        surface_volume_m3=r['after']['signed_volume_m3'])
                                   for r in cavity],
        cavity_convention_pending_count=len(cavity),
        authored_cavities_before=[r['entity_id'] for r in ok
                                  if r['before_patch_analysis'].get('authored_cavities')],
        patch_analysis_skipped=[r['entity_id'] for r in ok if not r['after_patch_analysis'].get('analysed', True)],
        thin_structures_radius_proxy_under_1mm=[r['entity_id'] for r in ok
                                                if (r['shape_after']['radius_proxy_m'] or 1.0) < 1e-3],
        max_absolute_relative_volume_drift=max((abs(r['signed_volume_relative_drift']) for r in ok
                                                if r['signed_volume_relative_drift'] is not None), default=0.0),
        elapsed_s=elapsed)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--output', type=Path)
    p.add_argument('--roles', nargs='*')
    p.add_argument('--limit', type=int, default=0)
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--order', choices=('anatomy', 'faces'), default='anatomy')
    a = p.parse_args()
    if a.self_test:
        self_test()
    if a.output:
        build(a.output, a.roles, a.limit, a.workers, a.order)
    if not a.self_test and not a.output:
        p.error('Select --self-test or --output')
