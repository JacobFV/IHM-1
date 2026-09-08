"""Verify the colour-palette endpoint reads what the build wrote and fails closed on everything else.

  .venv/bin/python scripts/verify_tissue_colour_palettes.py

The builder's own --self-test checks the palettes. This checks the reader: that ihm.app.experiments
.read_palette serves the index and each palette, and that it refuses a tampered artifact, a changed
input, a missing input, a build whose self-test did not pass, and an unknown palette id.

Every tamper is done inside a tempfile.TemporaryDirectory() copy of the two directories the reader
touches. Nothing is ever written under the repository root: a fixture that writes at a
caller-supplied root must never be handed the real tree.
"""
from pathlib import Path
import json
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.app.experiments import read_palette, PALETTES  # noqa: E402

MANIFEST_INPUT = 'data/derived/app/manifest.json'


def staged(tmp):
    """A minimal workspace holding only what read_palette reads. Never the repository root."""
    root = Path(tmp)/'workspace'
    (root/'data/derived/app').mkdir(parents=True)
    shutil.copytree(ROOT/PALETTES, root/PALETTES)
    shutil.copy(ROOT/MANIFEST_INPUT, root/MANIFEST_INPUT)
    return root


def refuses(root, expect, ident=None):
    try:
        read_palette(root, ident)
    except ValueError as error:
        assert expect in str(error), f'expected {expect!r}, got {error!r}'
        return str(error)
    raise AssertionError('read_palette accepted what it should have refused: '+expect)


def main():
    checks = []

    def check(name, detail):
        checks.append((name, detail))
        print('PASS '+name+' :: '+str(detail))

    with tempfile.TemporaryDirectory() as tmp:
        root = staged(tmp)
        assert not Path(root).resolve().is_relative_to(ROOT), 'fixture root must be outside the repository'

        index = read_palette(root)
        ids = [p['id'] for p in index['palettes']]
        assert index['default_palette'] == 'didactic', 'the app default must not change'
        check('index served and the default is unchanged', f'{ids}, default={index["default_palette"]}')

        served = {}
        for ident in ids:
            palette = read_palette(root, ident)
            assert palette['id'] == ident and len(palette['colours']) == 7390
            served[ident] = len(palette['colours'])
        check('every declared palette is served in full', served)

        assert index['sources'] and index['mucosa'], 'index must carry sources and the mucosal mapping'
        unresolved = [k for k, v in index['sources']['sources'].items() if not v.get('resolution_verified')]
        measured_sources = set()
        realistic = read_palette(root, 'realistic')
        for role in realistic['roles'].values():
            if role['tier'] == 'measured':
                measured_sources.update(role['sources'])
        assert not (measured_sources & set(unresolved)), 'an unresolved identifier backs a measured colour'
        check('no unresolved identifier backs a measured colour',
              f'{len(unresolved)} unresolved sources, {len(measured_sources)} cited by measured colours')

        check('unknown palette id refused', refuses(root, 'Unknown palette', 'no-such-palette'))

    with tempfile.TemporaryDirectory() as tmp:
        root = staged(tmp)
        target = root/PALETTES/'palettes/realistic.json'
        target.write_bytes(target.read_bytes().replace(b'"realistic"', b'"realistIc"', 1))
        check('tampered palette artifact refused', refuses(root, 'Palette artifact changed', 'realistic'))

    with tempfile.TemporaryDirectory() as tmp:
        root = staged(tmp)
        target = root/PALETTES/'index.json'
        target.write_bytes(target.read_bytes()+b'\n')
        check('tampered index refused', refuses(root, 'Palette artifact changed'))

    with tempfile.TemporaryDirectory() as tmp:
        root = staged(tmp)
        target = root/MANIFEST_INPUT
        target.write_bytes(target.read_bytes()+b' ')
        check('changed display manifest refused', refuses(root, 'Palette input changed'))

    with tempfile.TemporaryDirectory() as tmp:
        root = staged(tmp)
        (root/MANIFEST_INPUT).unlink()
        check('missing display manifest refused', refuses(root, 'Palette input is missing'))

    with tempfile.TemporaryDirectory() as tmp:
        root = staged(tmp)
        target = root/PALETTES/'manifest.json'
        data = json.loads(target.read_bytes())
        data['self_test']['passed'] = False
        target.write_bytes(json.dumps(data).encode())
        check('build that failed its own self-test refused', refuses(root, 'did not pass its own self-test'))

    with tempfile.TemporaryDirectory() as tmp:
        root = staged(tmp)
        target = root/PALETTES/'manifest.json'
        data = json.loads(target.read_bytes())
        data['outputs_sha256']['manifest.json'] = '0'*64
        target.write_bytes(json.dumps(data).encode())
        check('manifest that lists itself refused', refuses(root, 'lists itself'))

    with tempfile.TemporaryDirectory() as tmp:
        root = staged(tmp)
        target = root/PALETTES/'manifest.json'
        data = json.loads(target.read_bytes())
        data['schema'] = 'ihm.something-else.v1'
        target.write_bytes(json.dumps(data).encode())
        check('unexpected manifest schema refused', refuses(root, 'schema is not the one'))

    print(f'verified: {len(checks)} checks; palette endpoint serves and fails closed')


if __name__ == '__main__':
    main()
