#!/usr/bin/env python3
"""Archive the publisher's anatomy-derived lymph graph, never a generated graph."""
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw/lymphatic'
BASE = 'https://mdpi-res.com/d_attachment/'
SOURCES = {
    'mathematics-08-02236-s001.zip': ('mathematics/mathematics-08-02236/article_deploy/', '7c7b05eb6ae5dd7682ed0ed1fc89e8f08eb0c33d8bc3ac8711f97a5c07a0d390'),
    'mathematics-08-02236.xml': ('mathematics/mathematics-08-02236/article_deploy/', '13aa578a04f1cd1e535c05c5cce0ddfa00451f3b94b7dbb19d8f5628f8ca98f8'),
    'mathematics-08-02236.pdf': ('mathematics/mathematics-08-02236/article_deploy/', '3b5c01c352ec975a793f2ea5e9e147e8bba3dc7b347d114842b91ce1ff6355e8'),
    'computation-06-00001.pdf': ('computation/computation-06-00001/article_deploy/', '94024a082e8f318b4e8390cb15d86d3df7d4cc49bfd85372b64894d9961466e7'),
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def collect():
    RAW.mkdir(parents=True, exist_ok=True)
    records = []
    for name, (directory, expected) in SOURCES.items():
        path = RAW / name
        url = BASE + directory + name
        if not path.exists():
            temporary = path.with_suffix(path.suffix + '.partial')
            subprocess.run(['curl', '--fail', '--location', '--silent', '--show-error', '--max-time', '120', url, '-o', str(temporary)], check=True)
            if sha(temporary) != expected:
                raise ValueError(f'Publisher content changed: {name}; inspect before accepting new hash')
            temporary.replace(path)
        if sha(path) != expected:
            raise ValueError(f'Source checksum mismatch: {name}')
        records.append({'path': str(path.relative_to(ROOT)), 'url': url, 'sha256': expected, 'bytes': path.stat().st_size})
    archive = RAW / 'mathematics-08-02236-s001.zip'
    members = []
    with zipfile.ZipFile(archive) as z:
        for item in z.infolist():
            if item.is_dir():
                continue
            target = RAW / 'extracted' / item.filename
            if not target.resolve().is_relative_to((RAW / 'extracted').resolve()):
                raise ValueError('Unsafe archive path')
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(z.read(item))
            members.append({'path': str(target.relative_to(ROOT)), 'archive_member': item.filename, 'sha256': sha(target), 'bytes': target.stat().st_size})
    provenance = {
        'schema_version': 1, 'retrieved_date': '2026-09-05',
        'source_doi': '10.3390/math8122236', 'basis_doi': '10.3390/computation6010001',
        'attribution': 'Savinkov, Grebennikov, Puchkova, Chereshnev, Sazonov and Bocharov (2020). Graph Theory for Modeling and Analysis of the Human Lymphatic System. Mathematics 8(12), 2236. Basis: Tretyakova, Savinkov, Lobov and Bocharov (2017), Computation 6(1), 1; PlasticBoy project.',
        'source_revision': 'publisher supplement ZIP, member timestamps 2020-12-14; SHA-256 pinned',
        'license': 'CC-BY-4.0 (publisher article declaration; supplement has no separate license file)',
        'license_url': 'https://creativecommons.org/licenses/by/4.0/',
        'license_evidence': 'mathematics-08-02236.xml article/front/article-meta/permissions; authors publicly provide graph as article supplementary material',
        'third_party_basis': 'PlasticBoy authored anatomy; original commercial polygon assets are not acquired or redistributed here. Article license is not asserted for original PlasticBoy assets.',
        'evidence_kind': 'published_anatomy_derived_computational_graph',
        'measured_dataset': False, 'independent_subject_count': 0,
        'coordinate_units': 'mm', 'units_evidence': '2020 paper section 2.2: coordinates scaled to basal 1750 mm height; graph edge header Length(mm)',
        'wolfram_cross_reference': 'https://datarepository.wolframcloud.com/resources/63e705c2-bcc6-4458-959d-bd5f8baafd30/',
        'wolfram_acquisition': 'Resource HTTP 404 and cloud scheduled-upgrade page on 2026-09-05; publisher supplementary files used directly, not Wolfram-exported data.',
        'files': records, 'archive_members': members,
    }
    (RAW / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps({'files': len(records), 'archive_members': len(members), 'archive_sha256': sha(archive)}))
    return provenance

if __name__ == '__main__':
    collect()
