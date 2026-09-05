#!/usr/bin/env python3
"""Acquire pinned official Z-Anatomy archive and retain source provenance."""
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw/anatomy/extended'
REVISION = '702db27cf1395bf067158e86c7659b6e37a7be2c'
REPOSITORY = 'https://github.com/Z-Anatomy/Models-of-human-anatomy'
EXPECTED_SHA256 = {
    'Z-Anatomy.zip': 'e029688545627bd0214b269e1063143abb580aad72b2c2445d6d8a9a0d9da736',
    'Readme.md': '2439e8ccea759c644bbaf2f6fd48551fe054bfda30fd3d97d8384f75273466b3',
    'License.txt': '196b66b56551a862e59872f7cdb70e6d9a6ad84e9105962fca8ec28c14e97520',
}

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def collect():
    RAW.mkdir(parents=True, exist_ok=True)
    records = []
    for name in ['Z-Anatomy.zip', 'Readme.md', 'License.txt']:
        url = f'https://raw.githubusercontent.com/Z-Anatomy/Models-of-human-anatomy/{REVISION}/{name}'
        path = RAW / name
        if not path.exists():
            temporary = path.with_suffix(path.suffix + '.partial')
            urllib.request.urlretrieve(url, temporary)
            temporary.replace(path)
        if sha256(path) != EXPECTED_SHA256[name]: raise ValueError(f'Checksum mismatch: {path}')
        records.append({'path': str(path.relative_to(ROOT)), 'url': url, 'bytes': path.stat().st_size, 'sha256': sha256(path)})
    archive = RAW / 'Z-Anatomy.zip'
    with zipfile.ZipFile(archive) as z:
        for member in z.infolist():
            target = RAW / 'extracted' / member.filename
            if not target.resolve().is_relative_to((RAW / 'extracted').resolve()): raise ValueError('Unsafe archive member')
            if member.is_dir(): continue
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists() or target.stat().st_size != member.file_size: target.write_bytes(z.read(member))
    blend = RAW / 'extracted/Z-Anatomy/Startup.blend'
    if sha256(blend) != '9f08a17ea0115fed80b2a73ecdf0a1bc2ab2f6956f37c593ce23d513ea35afcd':
        raise ValueError('Extracted blend differs from the pinned archive; remove extracted file and recollect')
    provenance = {'schema_version': 1, 'repository': REPOSITORY, 'revision': REVISION,
        'files': records, 'blend_path': str(blend.relative_to(ROOT)), 'blend_sha256': sha256(blend),
        'source_category': 'authored anatomical reference atlas derived in part from BodyParts3D',
        'independent_subject_count': 0, 'population_sample': False,
        'declared_license': 'CC-BY-SA-4.0', 'license_url': 'https://creativecommons.org/licenses/by-sa/4.0/',
        'attribution': 'Z-Anatomy — The libre 3D atlas of anatomy; Gauthier Kervyn and contributors. BodyParts3D — The Database Center for Life Science — CC-BY-SA 2.1 Japan (as credited by Z-Anatomy).',
        'third_party_license_caveat': 'README additionally credits Inner Ear (University of Dundee, CC-BY-NC-SA-4.0) and Kidney (Lissie Cowley, CC-BY-NC-4.0); blanket project license does not resolve those component restrictions. Full raw archive retains these credits; display export excludes kidney-named objects and nervous/sense-organ collection.'}
    (RAW / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance, indent=2))

if __name__ == '__main__': collect()
