"""Independently accept the conforming-domain display: geometry, determinism and fail-closed reads.

Four properties, each measured rather than asserted from the build log:
  1. every owner's exported boundary is closed and correctly oriented -- the divergence volume of
     the exported triangles equals the summed signed volume of that owner's tetrahedra;
  2. no exported triangle is invented -- every one is a face of a tet the source volume owns;
  3. rebuilding from the same sources reproduces the committed display byte for byte;
  4. read_experiment refuses a changed source volume and a changed display, in an isolated root.

Property 4 copies the artifacts into a tempfile.TemporaryDirectory and tampers with the copies.
Nothing here writes anywhere under the repository except the temporary build directory it is given.

  .venv/bin/python scripts/verify_conforming_domain_display.py
"""
from pathlib import Path
import hashlib, gzip, json, shutil, sys, tempfile
import numpy as np

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
DISPLAY = BASE / 'data/derived/conforming-domain-display-v1'
FACES = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def signed_volume(points, triangles):
    a, b, c = (points[triangles[:, i]] for i in range(3))
    return float(np.einsum('ij,ij->i', np.cross(b - a, c - a), a).sum() / 6.0)


def check(label, ok, detail=''):
    print(f"{'PASS' if ok else 'FAIL'}  {label}{'  ' + detail if detail else ''}")
    return bool(ok)


def geometry(domain_id, tets, labels, vertices, label_of):
    """Property 1 and 2 for one exported domain."""
    payload = json.loads(gzip.decompress((DISPLAY / f'{domain_id}.json.gz').read_bytes()))
    points = np.asarray(payload['geometry']['positions_m'], np.float64).reshape(-1, 3)
    faces = {}
    for group in payload['groups']:
        index = np.asarray(group['triangles'], np.int64).reshape(-1, 3)
        for entry in group['structures']:
            faces[entry['entity_id']] = index[entry['triangle_start']:entry['triangle_start'] + entry['triangle_count']]
    volume_ok, member_ok, worst = True, True, 0.0
    source_faces = None
    for entity_id, triangles in faces.items():
        label = label_of[entity_id]
        block = tets[labels == label]
        tet_volume = float(np.einsum('ij,ij->i', np.cross(
            vertices[block[:, 1]] - vertices[block[:, 0]], vertices[block[:, 2]] - vertices[block[:, 0]]),
            vertices[block[:, 3]] - vertices[block[:, 0]]).sum() / 6.0)
        exported = signed_volume(points, triangles)
        error = abs(exported - tet_volume) / max(abs(tet_volume), 1e-15)
        worst = max(worst, error)
        if error > 1e-4:
            volume_ok = False
    # Property 2: the exported triangles are a subset of the source volume's own faces.
    if source_faces is None:
        keys = np.sort(tets[:, FACES].reshape(-1, 3), axis=1)
        source_faces = set(map(tuple, np.unique(keys, axis=0)))
    sample = np.concatenate([t for t in faces.values()])
    # Positions were deduplicated, so match on coordinates back into the source vertex array.
    quantized = np.asarray(vertices, np.float32)
    lookup = {row.tobytes(): i for i, row in enumerate(quantized)}
    mapped = np.array([lookup[row.tobytes()] for row in np.asarray(points, np.float32)], np.int64)
    for triangle in map(tuple, np.sort(mapped[sample], axis=1)):
        if triangle not in source_faces:
            member_ok = False
            break
    ok = check(f'{domain_id}: every owner boundary is closed and correctly oriented', volume_ok,
               f'worst relative divergence-volume error {worst:.2e}')
    return ok and check(f'{domain_id}: every exported triangle is a face of a source tetrahedron', member_ok,
                        f'{len(sample):,} triangles checked')


def main():
    from ihm.app.experiments import read_experiment
    ok = True
    index = json.loads((DISPLAY / 'index.json').read_bytes())

    hand = BASE / 'data/derived/hand-reflexive-grip-domain-v1'
    volume = np.load(hand / 'tetmesh-coarse.npz', allow_pickle=False)
    members = [json.loads(line) for line in (hand / 'members.jsonl').read_text().splitlines()]
    label_of = {member['entity_id']: i for i, member in enumerate(members)}
    label_of['unclaimed-complement'] = -1
    ok &= geometry('hand-reflexive-grip-coarse', volume['TT'], volume['owner'], volume['TV'], label_of)

    partition = BASE / 'data/derived/material-domains/whole-body-0.01m'
    grid = np.load(partition / 'whole-body-domain.npz', allow_pickle=False)
    ok &= geometry('whole-body-0.01m', grid['tetrahedra'], grid['material_index'], grid['vertices_m'],
                   {str(entity): i for i, entity in enumerate(grid['source_ids'])})

    # Property 3: determinism.
    from scripts.build_conforming_domain_display import build
    with tempfile.TemporaryDirectory() as scratch:
        rebuilt = build(scratch, [entry['id'] for entry in index['domains']])
        same = all(sha(Path(scratch) / entry['display_path']) == entry['display_sha256'] for entry in rebuilt['domains'])
        matches = all(a['display_sha256'] == b['display_sha256'] for a, b in zip(index['domains'], rebuilt['domains']))
        ok &= check('rebuilding from the same sources reproduces the committed display', same and matches,
                    f"{len(rebuilt['domains'])} domains")

    # Property 4: fail-closed reads in an isolated root.
    with tempfile.TemporaryDirectory() as scratch:
        root = Path(scratch)
        entry = next(e for e in index['domains'] if e['id'] == 'hand-reflexive-grip-coarse')
        for source in entry['source_hashes']:
            (root / source).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(BASE / source, root / source)
        target = root / 'data/derived/conforming-domain-display-v1'
        target.mkdir(parents=True, exist_ok=True)
        shutil.copy2(DISPLAY / entry['display_path'], target / entry['display_path'])
        (target / 'index.json').write_text(json.dumps({**index, 'domains': [entry]}, indent=2) + '\n')
        kind = 'conforming-domain-' + entry['id']
        payload = read_experiment(root, kind)
        ok &= check('isolated root reads the display', payload['id'] == entry['id'],
                    f"{payload['geometry']['triangle_count']:,} triangles")
        for name, path in (('source volume', root / 'data/derived/hand-reflexive-grip-domain-v1/members.jsonl'),
                           ('display export', target / entry['display_path'])):
            original = path.read_bytes()
            path.write_bytes(original + b'\n')
            try:
                read_experiment(root, kind)
                refused = False
            except ValueError:
                refused = True
            path.write_bytes(original)
            ok &= check(f'a changed {name} is refused', refused)
        ok &= check('the restored root reads again', read_experiment(root, kind)['id'] == entry['id'])
    print('OK' if ok else 'FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
