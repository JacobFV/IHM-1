#!/usr/bin/env python
"""Retain the acquired CC0 MakeHuman garment meshes under data/raw with hashes.

Downloads (or reuses a staging copy of) the MakeHuman system asset pack and the
CC0 community clothing packs, extracts only the meshes named by the wardrobe
catalogue, and writes a retrieval receipt carrying every archive URL, archive
sha256, per-file sha256, per-asset author and per-asset licence as declared by
the upstream pack manifest.

  python scripts/retain_makehuman_garment_sources.py [--staging DIR] [--self-test]
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.garment_wardrobe import CATALOGUE, sha256_file  # noqa: E402

RAW = ROOT / 'data/raw/clothing/makehuman'
FILES_BASE = 'https://files2.makehumancommunity.org/asset_packs'
SITE_BASE = ('https://raw.githubusercontent.com/makehumancommunity/'
             'makehuman-static-website/master/content/Assets/AssetPacks')
REPO_BASE = 'https://raw.githubusercontent.com/makehumancommunity/makehuman/master'
SYSTEM_PACK_URL = f'{FILES_BASE}/makehuman_system_assets/makehuman_system_assets_cc0.zip'
SYSTEM_PACK_PAGE = 'https://static.makehumancommunity.org/assets/assetpacks/makehuman_system_assets.html'
BASE_MESH_URL = f'{REPO_BASE}/makehuman/data/3dobjs/base.obj'
LICENCE_URLS = {'LICENSE.md': f'{REPO_BASE}/LICENSE.md', 'LICENSE.ASSETS.md': f'{REPO_BASE}/LICENSE.ASSETS.md'}

DATASET = {
    'id': 'makehuman-cc0-clothes',
    'label': 'MakeHuman CC0 clothing assets',
    'version': 'system asset pack + CC0 community asset packs, retrieved 2026-09-07',
    'url': 'https://static.makehumancommunity.org/assets/assetpacks.html',
    'specimen': 'garment meshes authored against the MakeHuman base mesh; no physical garment measured',
    'units': 'decimetre (MakeHuman authoring unit)',
    'frame': 'makehuman-base-mesh',
    'attribution': 'MakeHuman project and named asset authors; CC0 waives the attribution requirement',
    'acquisition_status': 'downloaded from the upstream asset distribution',
    'license': 'CC0-1.0',
    'license_url': 'https://creativecommons.org/publicdomain/zero/1.0/',
}
ROW = re.compile(r'^\|\s*clothes\s*\|[^|]*\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|\s*\[[^\]]*\]\(([^)]+)\)\s*\|\s*([^|]+?)\s*\|')


def fetch(url, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=600) as response, open(destination, 'wb') as handle:
        shutil.copyfileobj(response, handle)
    return destination


def pack_licences(text):
    out = {}
    for line in text.splitlines():
        match = ROW.match(line)
        if match:
            out[match.group(1)] = {'author': match.group(2), 'source_url': match.group(3), 'license': match.group(4)}
    return out


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--staging', type=Path, default=None,
                        help='directory already holding the downloaded archives and extracts')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args(argv)
    if args.self_test:
        text = ('| Asset type | Thumbnail | Asset name | Author | Source | License |\n'
                '| clothes | ![a.png](a.png) | toigo_gloves_long | MargaretToigo | [asset repo](http://x/1) | CC0 |')
        parsed = pack_licences(text)
        assert parsed == {'toigo_gloves_long': {'author': 'MargaretToigo', 'source_url': 'http://x/1', 'license': 'CC0'}}, parsed
        assert {e['pack'] for e in CATALOGUE} - {'system'}, 'catalogue lost its community packs'
        assert all(e['obj'].endswith('.obj') for e in CATALOGUE)
        print('retain_makehuman_garment_sources self-test: pass')
        return 0

    staging = args.staging or (ROOT / 'data/derived/.staging-makehuman')
    staging.mkdir(parents=True, exist_ok=True)
    packs = sorted({e['pack'] for e in CATALOGUE if e['pack'] != 'system'})
    retrieved = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()

    archives = {}
    for pack in packs + ['makehuman_system_assets']:
        url = SYSTEM_PACK_URL if pack == 'makehuman_system_assets' else f'{FILES_BASE}/{pack}/{pack}_cc0.zip'
        archive = staging / f'{pack}_cc0.zip'
        if not archive.exists():
            print(f'downloading {url}')
            fetch(url, archive)
        extract = staging / ('mhassets' if pack == 'makehuman_system_assets' else f'packext/{pack}')
        if not extract.exists():
            with zipfile.ZipFile(archive) as zf:
                zf.extractall(extract)
        archives[pack] = {'url': url, 'sha256': sha256_file(archive), 'bytes': archive.stat().st_size,
                          'retrieved_at': retrieved, 'retained': False,
                          'retained_reason': 'archive is 6-272 MB and carries textures this build does not use; '
                                             'its sha256 and URL are recorded and only the used meshes are retained'}

    RAW.mkdir(parents=True, exist_ok=True)
    licences = {}
    for pack in packs:
        page = RAW / 'packs' / f'{pack}.md'
        if not page.exists():
            fetch(f'{SITE_BASE}/{pack}.md', page)
        licences[pack] = pack_licences(page.read_text())
    system_page = RAW / 'packs' / 'makehuman_system_assets.txt'
    if not system_page.exists():
        system_page.parent.mkdir(parents=True, exist_ok=True)
        system_page.write_text(
            'MakeHuman system asset pack\n'
            f'page: {SYSTEM_PACK_PAGE}\n'
            f'archive: {SYSTEM_PACK_URL}\n'
            'licence as stated on the pack page: CC0\n'
            'licence of every MHCLO-based asset per the project LICENSE.md section C: CC0 1.0 Universal\n')
    for name, url in LICENCE_URLS.items():
        if not (RAW / name).exists():
            fetch(url, RAW / name)
    base = RAW / 'base/base.obj'
    if not base.exists():
        fetch(BASE_MESH_URL, base)

    files = {}
    for entry in CATALOGUE:
        if entry['pack'] == 'system':
            source = staging / 'mhassets/clothes' / entry['asset'] / entry['obj']
            licence = {'author': 'MakeHuman project', 'license': 'CC0',
                       'source_url': SYSTEM_PACK_PAGE}
        else:
            source = staging / 'packext' / entry['pack'] / 'clothes' / entry['asset'] / entry['obj']
            licence = licences[entry['pack']].get(entry['asset'])
            if licence is None:
                raise SystemExit(f'{entry["id"]}: no upstream licence row for asset {entry["asset"]!r}')
        if licence['license'].strip().upper() not in {'CC0', 'CC0-1.0'}:
            raise SystemExit(f'{entry["id"]}: upstream licence {licence["license"]!r} is not CC0; refusing to retain')
        target = RAW / 'clothes' / entry['pack'] / entry['asset'] / entry['obj']
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            shutil.copyfile(source, target)
        files[entry['id']] = {
            'path': str(target.relative_to(ROOT)), 'sha256': sha256_file(target), 'bytes': target.stat().st_size,
            'pack': entry['pack'], 'asset': entry['asset'], 'author': licence['author'],
            'license': licence['license'], 'asset_page': licence['source_url'],
            'archive_url': archives[entry['pack'] if entry['pack'] != 'system' else 'makehuman_system_assets']['url'],
            'retrieved_at': retrieved,
        }

    receipt = {
        'schema': 'ihm.raw-retrieval.v1',
        'dataset': DATASET,
        'licence_evidence': {name: {'path': str((RAW / name).relative_to(ROOT)), 'sha256': sha256_file(RAW / name)}
                             for name in LICENCE_URLS},
        'licence_rule': 'LICENSE.md section C: "Clothes (any MHCLO-based asset) ... released under CC0 1.0 Universal". '
                        'Every retained asset additionally carries a per-asset CC0 row in its pack manifest.',
        'source_body': {'path': str(base.relative_to(ROOT)), 'sha256': sha256_file(base),
                        'bytes': base.stat().st_size, 'url': BASE_MESH_URL, 'retrieved_at': retrieved,
                        'license': 'CC0-1.0', 'license_note': 'the file header states the CC0 release and its 2020 '
                                                              'copyright holders verbatim'},
        'pack_manifests': {pack: {'path': str((RAW / 'packs' / f'{pack}.md').relative_to(ROOT)),
                                  'sha256': sha256_file(RAW / 'packs' / f'{pack}.md'),
                                  'url': f'{SITE_BASE}/{pack}.md'} for pack in packs},
        'archives': archives,
        'files': files,
    }
    out = RAW / 'retrieval.json'
    out.write_text(json.dumps(receipt, indent=1, sort_keys=True) + '\n')
    print(f'retained {len(files)} garment meshes -> {out.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
