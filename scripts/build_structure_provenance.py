#!/usr/bin/env python3
"""Uniform provenance for every structure, whatever its origin. Candidate only; canonical assets are not written.

  python scripts/build_structure_provenance.py [--self-test] [--no-verify-hashes]

Emits data/derived/structure-provenance-candidate-v1/:
  schema.json   the one record shape every structure carries
  index.json    dataset registry, counts, per-record sha256, audit summary
  records/*.json one record per display structure and per canonical entity
  audit.json    field-by-field coverage before and after backfill, per model/system/role
  manifest.json inputs and outputs hashed
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.anatomy import physical_role
OUT = ROOT/'data/derived/structure-provenance-candidate-v1'
PROMOTED = ROOT/'data/derived/structure-provenance-v1'
SCHEMA_ID = 'ihm.structure-provenance.v1'

TIERS = {
    'measured': 'Geometry acquired from a dataset that describes a physical specimen; no fitted cross-source transform was applied to place it.',
    'transferred': 'Geometry from a separate source placed into this body by a fitted transform, so its position here is inferred. The fit residual is recorded where it was measured.',
    'derived': 'Geometry produced from other geometry by a recorded authoring or computation step, with no independent acquisition of its own.',
    'synthesized': 'Geometry or extent constructed from explicit priors with no source geometry of its own.'}

# Files the repository already retains that evidence a licence. No licence is asserted without one.
LICENCE_EVIDENCE = {
    'bodyparts3d': None,  # declared in data/derived/anatomy/bodyparts3d_index.json, hashed as an input
    'z-anatomy': 'data/raw/anatomy/extended/License.txt',
    'vascular': 'data/raw/vascular/vmr/LICENSE',
    'betse-tissue': 'data/raw/physiology/betse/LICENSE'}

# Structures whose provenance lives in a build receipt rather than in an upstream atlas file.
RECEIPT_BUILDS = {
    'body-detail-hair': ('scripts/build_body_details.py', 'data/derived/canonical/details.json'),
    'body-detail-microvascular': ('scripts/build_body_details.py', 'data/derived/canonical/details.json')}

# Where to look for a retained file when only its hash was recorded. Bounded, per dataset.
HASH_SEARCH_ROOTS = {
    'vascular-aorta': ['data/raw/vascular/vmr/extracted/0001_H_AO_SVD'],
    'vascular-cerebral': ['data/raw/vascular/vmr/extracted/0050_H_CERE_H'],
    'vascular-pulmonary': ['data/raw/vascular/vmr/extracted/0077_H_PULM_H'],
    'published-lymphatic-network': ['data/raw/lymphatic'],
    'betse-tissue': ['data/derived/bioelectric', 'data/raw/physiology/betse'],
    'opensim-rajagopal': ['data/raw/anatomy/opensim-models/source/Models/Rajagopal']}

REQUIRED = ['dataset.id', 'dataset.label', 'dataset.license', 'source_file.path', 'source_file.sha256',
            'build.script', 'build.commit', 'geometry.path', 'geometry.sha256', 'transforms', 'tier']


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def sha256_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_bytes())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=1, sort_keys=False))
    return sha256_file(path)


def git(*args):
    try:
        return subprocess.run(['git', *args], cwd=ROOT, capture_output=True, text=True, timeout=30).stdout.strip() or None
    except Exception:
        return None


_WORKING_TREE = {}


def working_tree():
    """What the repository actually looked like when this build ran.

    A build whose entire purpose is traceability must not answer 'which commit produced
    this' with a bare null. It usually cannot answer with a commit either: this tree is
    routinely dirty, and a commit id would then name content that is not what ran. So it
    answers with all three of the things that are true -- the commit HEAD is on, the tree
    that commit points at, and an explicit digest over every path git reports as modified,
    taken from the bytes on disk. Together those pin the exact state; the commit alone
    does not, and null says nothing at all.

    Read-only. Nothing here writes an object, an index entry or a ref.
    """
    if _WORKING_TREE:
        return dict(_WORKING_TREE)
    head = git('rev-parse', 'HEAD')
    porcelain = git('status', '--porcelain') or ''
    dirty = sorted({line[3:].split(' -> ')[-1].strip().strip('"')
                    for line in porcelain.splitlines() if line.strip()})
    entries = []
    for relative in dirty:
        path = ROOT/relative
        entries.append(relative+':'+(sha256_file(path) if path.is_file() else 'absent'))
    _WORKING_TREE.update({
        'head_commit': head,
        'head_tree': git('rev-parse', 'HEAD^{tree}'),
        'head_committed_at': git('log', '-1', '--format=%cI'),
        'clean': not dirty,
        'dirty_paths': dirty,
        'dirty_path_count': len(dirty),
        'dirty_digest': hashlib.sha256('\n'.join(entries).encode()).hexdigest() if entries else None,
        'dirty_digest_method': 'sha256 over sorted "<path>:<sha256 of the bytes on disk>" lines for '
                               'every path git status reports, "absent" for a deleted one',
        'guarantees': 'head_commit plus head_tree plus dirty_digest identify the exact bytes this '
                      'build read. head_commit alone does not while dirty_path_count is nonzero.',
        'does_not_guarantee': 'the dirty content is not stored anywhere, so this state can be '
                              'recognised again but not reconstructed from the repository. Only a '
                              'build on a clean tree is reproducible from its commit alone.'})
    return dict(_WORKING_TREE)


class Builds:
    """Every retained derived build that names a geometry file, inverted path->build and sha256->build."""
    FIELDS = ('source_geometry_sha256', 'per_entity_input_sha256', 'inputs_sha256', 'source_hashes', 'artifacts_sha256')

    def __init__(self, verify=True):
        self.by_path, self.by_hash, self.manifests = defaultdict(list), defaultdict(list), []
        for name in ('manifest.json', 'manifest_fragment.json', 'validation.json', 'inventory.json', 'index.json'):
            for path in sorted((ROOT/'data/derived').glob('*/'+name)) + sorted((ROOT/'data/derived').glob('*/*/'+name)):
                if path.stat().st_size > 40 << 20 or 'structure-provenance-candidate' in str(path):
                    continue
                try:
                    data = read_json(path)
                except Exception:
                    continue
                if not isinstance(data, dict):
                    continue
                self._absorb(path, data, verify)

    def _absorb(self, path, data, verify):
        relative = str(path.relative_to(ROOT))
        found = False
        for field in self.FIELDS:
            entries = data.get(field)
            if not isinstance(entries, dict):
                continue
            for key, value in entries.items():
                if not isinstance(value, str) or len(value) != 64:
                    continue
                found = True
                reference = {'build': str(path.parent.relative_to(ROOT)), 'manifest': relative, 'field': field, 'key': key}
                if key.startswith('data/'):
                    self.by_path[key].append(reference)
                self.by_hash[value].append(reference)
        single = data.get('source_geometry_path')
        if isinstance(single, str) and isinstance(data.get('source_geometry_sha256'), str):
            found = True
            reference = {'build': str(path.parent.relative_to(ROOT)), 'manifest': relative, 'field': 'source_geometry_sha256', 'key': single}
            self.by_path[single].append(reference)
            self.by_hash[data['source_geometry_sha256']].append(reference)
        if found:
            self.manifests.append({'path': relative, 'sha256': sha256_file(path) if verify else None,
                                   'schema': data.get('schema') or data.get('schema_version')})

    def lookup(self, geometry_path, geometry_sha):
        seen, out = set(), []
        for reference in self.by_path.get(geometry_path, []):
            key = (reference['build'], reference['field'], reference['key'])
            seen.add(key)
            out.append({**reference, 'matched_by': 'path'})
        for reference in self.by_hash.get(geometry_sha or '', []):
            key = (reference['build'], reference['field'], reference['key'])
            if key in seen:
                for item in out:
                    if (item['build'], item['field'], item['key']) == key:
                        item['matched_by'] = 'path+sha256'
                continue
            out.append({**reference, 'matched_by': 'sha256'})
        return sorted(out, key=lambda r: (r['build'], r['field'], r['key']))


class Provenance:
    def __init__(self, verify=True, out=OUT):
        self.verify = verify
        self.out = Path(out)
        self.manifest = read_json(ROOT/'data/derived/app/manifest.json')
        self.anatomy = read_json(ROOT/'data/derived/canonical/anatomy.json')
        self.bp_index = {m['element_id']: m for m in read_json(ROOT/'data/derived/anatomy/bodyparts3d_index.json')['meshes']}
        self.bp_head = read_json(ROOT/'data/derived/anatomy/bodyparts3d_index.json')
        self.za_index = {m['id']: m for m in read_json(ROOT/'data/derived/anatomy/extended/source_index.json')['meshes']}
        self.models = {m['id']: m for m in self.manifest['models']}
        self.entities = {e['id']: e for e in self.anatomy['entities']}
        self.registrations = self.anatomy['registrations']
        self.builds = Builds(verify)
        self.scripts = {}
        self.hash_cache = {}
        self.resolved = self._load_resolution()
        self.searched = set()
        self.receipts = dict(RECEIPT_BUILDS)
        self._discover_receipts()
        self.datasets = self._datasets()

    def _discover_receipts(self):
        """A structure whose geometry sits beside a build receipt is owned by that build, not by whatever entity it hangs off."""
        seen = set()
        for structure in self.manifest['structures']:
            path = structure.get('geometry_path')
            if not path or structure['id'] in self.receipts:
                continue
            directory = str(Path(path).parent)
            if directory in ('data/derived/canonical/geometry', 'data/derived/anatomy/extended/geometry', 'data/derived/app/geometry'):
                continue
            for name in ('manifest_fragment.json', 'manifest.json'):
                receipt = f'{directory}/{name}'
                if (directory, name) in seen or not (ROOT/receipt).is_file():
                    continue
                seen.add((directory, name))
                data = read_json(ROOT/receipt)
                listed = {s['id'] for s in data.get('structures', []) if isinstance(s, dict)}
                scripts = [k for k in (data.get('source_hashes') or {}) if k.startswith('scripts/') and k.endswith('.py')]
                if not listed or not scripts:
                    continue
                for identity in listed:
                    self.receipts.setdefault(identity, (scripts[0], receipt))

    # --- shared helpers -------------------------------------------------
    def _load_resolution(self):
        path = self.out/'hash-resolution.json'
        if not path.is_file():
            return {}
        kept = {}
        for digest, relative in read_json(path).get('resolved', {}).items():
            if relative and (ROOT/relative).is_file() and self.file_hash(relative) == digest:
                kept[digest] = relative
        return kept

    def resolve_hash(self, digest, model_id):
        """Find the retained file whose content matches a recorded hash. Bounded to that dataset's roots."""
        if not digest:
            return None
        if digest in self.resolved:
            return self.resolved[digest]
        roots = HASH_SEARCH_ROOTS.get(model_id, [])
        for relative in roots:
            if relative in self.searched:
                continue
            self.searched.add(relative)
            base = ROOT/relative
            if not base.exists():
                continue
            for current, _, names in os.walk(base):
                for name in names:
                    file = Path(current)/name
                    try:
                        if file.stat().st_size > (512 << 20):
                            continue
                        self.resolved.setdefault(sha256_file(file), str(file.relative_to(ROOT)))
                    except OSError:
                        continue
        return self.resolved.get(digest)

    def file_hash(self, relative):
        if relative not in self.hash_cache:
            path = ROOT/relative
            self.hash_cache[relative] = sha256_file(path) if path.is_file() else None
        return self.hash_cache[relative]

    def script(self, relative):
        if relative not in self.scripts:
            path = ROOT/relative
            commit = git('log', '-1', '--format=%H', '--', relative)
            dirty = git('status', '--porcelain', '--', relative)
            tracked = git('ls-files', '--', relative)
            if commit:
                reason = None
            elif not path.is_file():
                reason = 'the script named by this record is not present in the working tree'
            elif not tracked:
                reason = 'the script is present but not tracked in git, so no commit contains it'
            else:
                reason = 'git is unavailable or this path has no commit history'
            self.scripts[relative] = {'script': relative, 'script_sha256': sha256_file(path) if path.is_file() else None,
                                      'commit': commit, 'commit_absent_reason': reason,
                                      'commit_covers_working_tree': (dirty == '' or dirty is None) and commit is not None,
                                      'uncommitted_changes': bool(dirty),
                                      'working_tree': {k: working_tree()[k] for k in
                                                       ('head_commit', 'head_tree', 'clean', 'dirty_path_count', 'dirty_digest')}}
        return dict(self.scripts[relative])

    def licence(self, dataset_key, declared, url=None, declared_in=None):
        evidence = LICENCE_EVIDENCE.get(dataset_key)
        record = {'license': declared, 'license_url': url, 'license_evidence': None}
        if evidence and (ROOT/evidence).is_file():
            record['license_evidence'] = {'path': evidence, 'sha256': self.file_hash(evidence)}
        elif declared:
            declared_in = declared_in or 'data/derived/app/manifest.json'
            record['license_evidence'] = {'path': declared_in, 'sha256': self.file_hash(declared_in),
                                          'note': 'declared in the retained source metadata; no separate licence file retained'}
        else:
            record['license_evidence'] = None
            record['license_absent_reason'] = 'no licence declaration retained for this dataset in this repository'
        return record

    def _datasets(self):
        bp = self.bp_head
        za = self.models['z-anatomy']
        out = {
            'bodyparts3d': {'id': 'bodyparts3d', 'label': 'BodyParts3D 4.0', 'version': '4.0', 'revision': None,
                            'url': 'https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html',
                            'specimen': 'adult male reference atlas', 'units': 'mm',
                            'frame': 'bodyparts3d-mm-left-posterior-superior', 'attribution': bp.get('attribution'),
                            'acquisition_status': 'acquired', **self.licence('bodyparts3d', bp.get('license'), 'https://creativecommons.org/licenses/by/4.0/',
                                                          'data/derived/anatomy/bodyparts3d_index.json')},
            'z-anatomy': {'id': 'z-anatomy', 'label': 'Z-Anatomy', 'version': None, 'revision': za.get('source_revision'),
                          'url': 'https://github.com/Z-Anatomy/Models-of-human-anatomy',
                          'specimen': 'authored reference atlas; partly derived from BodyParts3D, not an additional measured human',
                          'units': za.get('source_units'), 'frame': za.get('source_frame'), 'attribution': za.get('attribution'),
                          'acquisition_status': 'authored', **self.licence('z-anatomy', za.get('license'), za.get('license_url'))},
            'opensim-rajagopal': {'id': 'opensim-rajagopal', 'label': 'Rajagopal2016 · published OpenSim model', 'version': '2016',
                                  'revision': self.models['opensim-rajagopal']['source']['source_revision'],
                                  'url': self.models['opensim-rajagopal']['source']['url'],
                                  'specimen': 'Rajagopal2016 source model; not the BodyParts3D subject', 'units': 'm',
                                  'frame': 'opensim-ground-m', 'attribution': self.models['opensim-rajagopal'].get('attribution'),
                                  'acquisition_status': 'published model default pose', **self.licence('opensim-rajagopal', None)},
            'published-lymphatic-network': {'id': 'published-lymphatic-network', 'label': 'Savinkov et al. 2020 · published lymphatic graph',
                                            'version': '2020', 'revision': self.models['published-lymphatic-network']['source']['source_revision'],
                                            'url': 'https://doi.org/10.3390/math8122236',
                                            'specimen': 'PlasticBoy-derived anatomical estimate; no measured subject', 'units': 'mm',
                                            'frame': 'savinkov2020-xyz-mm', 'attribution': self.models['published-lymphatic-network'].get('attribution'),
                                            'acquisition_status': 'published structural model; uncalibrated',
                                            **self.licence('published-lymphatic-network', self.models['published-lymphatic-network']['source'].get('license'))},
            'betse-tissue': {'id': 'betse-tissue', 'label': 'BETSE 1.5.1 generic tissue simulation', 'version': '1.5.1',
                      'revision': '4bba0624092f02ef52711d715d16ae7007d49aea', 'url': 'https://github.com/betsee/betse',
                      'specimen': 'generic computational tissue; no human subject', 'units': 'm',
                      'frame': 'betse-generated-planar-tissue', 'attribution': None,
                      'acquisition_status': 'executed source simulation',
                      **self.licence('betse-tissue', 'BSD 2-Clause text retained verbatim; see license_evidence',
                                     'https://github.com/betsee/betse/blob/master/LICENSE')},
            'ihm-synthesis': {'id': 'ihm-synthesis', 'label': 'IHM · built in this repository from retained inputs', 'version': None, 'revision': None,
                              'url': None, 'specimen': 'no source specimen', 'units': 'm', 'frame': 'bodyparts3d-display-m',
                              'attribution': None, 'acquisition_status': 'synthesized in this repository',
                              **self.licence('ihm-synthesis', None)}}
        for key in ('vascular-aorta', 'vascular-cerebral', 'vascular-pulmonary'):
            model = self.models[key]
            meta = model.get('source_metadata', {})
            source = model.get('source') or next((s['source'] for s in self.manifest['structures']
                                                  if s['model_id'] == key and s.get('source')), {})
            out[key] = {'id': key, 'label': source.get('label'), 'version': meta.get('Image Number'),
                        'revision': meta.get('SDR DOI') or meta.get('DOI'), 'url': 'https://www.vascularmodel.com/',
                        'specimen': f"{meta.get('Species','human')} {meta.get('Sex','?')} age {meta.get('Age','?')}; {meta.get('Disease','?')}"
                                    if meta else source.get('specimen'),
                        'units': source.get('units'), 'frame': source.get('frame'),
                        'attribution': meta.get('Citation'), 'acquisition_status': source.get('status'),
                        **self.licence('vascular', 'Stanford/UC/OSMSC data licence, research and development use with attribution; see license_evidence',
                                      'https://www.vascularmodel.com/')}
        return out

    # --- per-structure resolution --------------------------------------
    def source_file(self, structure, entity):
        """The specific upstream file this geometry was built from, with its hash."""
        identity = structure['id']
        source = structure.get('source') or {}
        if entity and entity.get('provenance', {}).get('files'):
            record = entity['provenance']['files'][0]
            return self._checked(record['path'], record['sha256'])
        if structure['model_id'] == 'bodyparts3d':
            raw = self.bp_index.get(identity.removeprefix('bp3d-'))
            if raw:
                return self._checked(raw['source_path'], raw['sha256'], bytes_=raw.get('bytes'))
        if structure['model_id'] == 'z-anatomy':
            raw = self.za_index.get(identity)
            if raw:
                return self._checked(raw['source_geometry_path'], raw['source_geometry_sha256'],
                                     note='Blender-evaluated viewport surface extracted from '
                                          + str(read_json(ROOT/'data/derived/anatomy/extended/source_index.json').get('source_blend_sha256', ''))[:12])
        if structure.get('geometry_source') and structure.get('geometry_source_sha256'):
            return self._checked(structure['geometry_source'], structure['geometry_source_sha256'])
        if identity in self.receipts:
            receipt = self.receipts[identity][1]
            return self._checked(receipt, self.file_hash(receipt),
                                 note='build receipt; it carries the hashes of the upstream geometry and code it consumed')
        if source.get('path'):
            return self._checked(source['path'], source.get('sha256'))
        if source.get('sha256'):
            located = self.resolve_hash(source['sha256'], structure['model_id'])
            if located:
                return self._checked(located, source['sha256'],
                                     note='located by content hash; the retained metadata recorded the hash without a path')
            return {'path': None, 'sha256': source['sha256'], 'sha256_verified': None, 'bytes': None,
                    'absent_reason': 'source hash retained; no file with that content is present in this working tree'}
        return {'path': None, 'sha256': None, 'sha256_verified': None, 'bytes': None,
                'absent_reason': 'no upstream file recorded for this structure'}

    def _checked(self, path, expected, bytes_=None, note=None):
        record = {'path': path, 'sha256': expected, 'bytes': bytes_, 'sha256_verified': None}
        if note:
            record['note'] = note
        if self.verify and expected:
            actual = self.file_hash(path)
            record['sha256_verified'] = actual == expected
            if actual is None:
                record['sha256_verified'] = None
                record['absent_reason'] = 'recorded source file is not present in this working tree'
        return record

    def geometry(self, structure, entity):
        path = structure.get('geometry_path') or (f"data/derived/app/geometry/{structure['id']}.json.gz")
        expected = structure.get('geometry_sha256')
        record = {'path': path, 'sha256': expected, 'sha256_verified': None,
                  'representation': (entity or {}).get('reference_geometry', {}).get('representation') or structure.get('kind'),
                  'frame': self.models[structure['model_id']].get('frame'), 'units': self.models[structure['model_id']].get('display_units'),
                  'source_vertex_count': (entity or {}).get('source_vertex_count') or (structure.get('display_geometry') or {}).get('source_vertices'),
                  'source_face_count': (entity or {}).get('source_face_count') or (structure.get('display_geometry') or {}).get('source_faces')}
        if self.verify and expected:
            actual = self.file_hash(path)
            record['sha256_verified'] = None if actual is None else actual == expected
            if actual is None:
                record['absent_reason'] = 'recorded geometry file is not present in this working tree'
        return record

    def transforms(self, structure, entity):
        model = self.models[structure['model_id']]
        chain = []
        display = model.get('display_transform')
        provenance = (entity or {}).get('provenance', {})
        if structure['id'] in self.receipts:
            chain.append({'kind': 'constructed_directly_in_canonical_frame', 'from': model['frame'], 'to': model['frame'],
                          'residual': None, 'residual_reason': 'built in the canonical frame from canonical inputs; no placement transform applied'})
        elif entity and provenance.get('source_to_canonical'):
            chain.append({'kind': 'rigid_rotation_scale_translation', 'from': self.datasets_for(structure)['frame'],
                          'to': model['frame'], **provenance['source_to_canonical'],
                          'residual': None, 'residual_reason': 'exact affine recorded at build time; no fit was performed'})
        elif display and not entity:
            chain.append({'kind': 'rigid_rotation_scale_translation', 'from': model.get('source_frame'), 'to': model['frame'],
                          **{k: v for k, v in display.items()}, 'residual': None,
                          'residual_reason': 'display placement transform recorded at build time; no fit was performed'})
        registration_id = provenance.get('registration_id') or (entity or {}).get('uncertainty', {}).get('registration_or_synthesis', {}).get('registration_id')
        report = self.registrations.get(registration_id) if registration_id else None
        if report and registration_id == 'z_anatomy':
            chain.append({'kind': 'affine_plus_thin_plate_spline', 'registration_id': registration_id,
                          'from': self.datasets['z-anatomy']['frame'], 'to': model['frame'],
                          'landmark_count': len(report.get('landmarks', [])),
                          'residual': {'metric': 'held_out_rms_m', 'value': report.get('held_out_rms_m'),
                                       'max_m': report.get('held_out_max_m'), 'method': report.get('held_out_method'),
                                       'fitted_rms_m': report.get('rms_m') or report.get('residual_rms_m')}})
        if report and registration_id == 'lymphatic_pose':
            chain.append({'kind': 'body_extent_alignment_plus_arm_centerline_pose_map', 'registration_id': registration_id,
                          'from': self.datasets['published-lymphatic-network']['frame'], 'to': model['frame'],
                          'body_scale_xyz': report.get('body_scale_xyz'), 'body_translation_m': report.get('body_translation_m'),
                          'residual': {'metric': 'node_group_association_rms_m', 'value': report.get('node_group_association_rms_m'),
                                       'method': report.get('node_group_association_threshold_m') and 'nearest lymph-node-group distance within 60 mm',
                                       'note': 'diagnostic distance, not a landmark-validated fit residual; anatomically_labeled_source_landmarks = '
                                               + str(report.get('anatomically_labeled_source_landmarks'))}})
        if not chain and entity:
            parent = next((c for c in entity.get('connections', []) if c.get('relation') == 'layer_of'), None)
            if parent:
                chain.append({'kind': 'inherited_from_parent_surface', 'parent_entity_id': parent['entity_id'],
                              'from': model['frame'], 'to': model['frame'], 'residual': None,
                              'residual_reason': 'shell quadrature layer over the parent surface; no transform of its own'})
        if not chain and structure['id'] in self.receipts:
            chain.append({'kind': 'constructed_directly_in_canonical_frame', 'from': model['frame'], 'to': model['frame'],
                          'residual': None, 'residual_reason': 'built in the canonical frame from canonical inputs; no placement transform applied'})
        if not chain:
            return [], 'no placement transform is recorded for this structure'
        return chain, None

    def datasets_for(self, structure):
        identity, model_id = structure['id'], structure['model_id']
        if model_id == 'ihm-body':
            if identity.startswith('body-skin-') or identity in self.receipts:
                return self.datasets['ihm-synthesis']
            if identity.startswith('body-za-') or identity.startswith('body-') and identity[5:].startswith('za-'):
                return self.datasets['z-anatomy']
            if 'lymphatic' in identity:
                return self.datasets['published-lymphatic-network']
            return self.datasets['bodyparts3d']
        return self.datasets.get(model_id, self.datasets['ihm-synthesis'])

    def build_of(self, structure):
        model_id = structure['model_id']
        if structure['id'] in self.receipts:
            return self.script(self.receipts[structure['id']][0])
        if model_id == 'ihm-body':
            return self.script('scripts/build_canonical_anatomy.py')
        if model_id == 'z-anatomy':
            return self.script('scripts/build_extended_anatomy.py')
        if model_id == 'opensim-rajagopal':
            return self.script('scripts/build_opensim_display.py')
        if model_id == 'published-lymphatic-network':
            return self.script('scripts/build_lymph_network.py')
        if model_id == 'betse-tissue':
            return self.script('scripts/export_betse_tissue.py')
        return self.script('scripts/build_spatial_atlas.py')

    def tier(self, structure, entity, chain):
        identity, model_id = structure['id'], structure['model_id']
        dataset = self.datasets_for(structure)
        evidence_kind = structure.get('evidence_kind') or (entity or {}).get('evidence_kind')
        quoted = {'source.specimen': (structure.get('source') or {}).get('specimen'),
                  'source.status': (structure.get('source') or {}).get('status'),
                  'evidence_kind': evidence_kind,
                  'model.calibration_status': self.models[model_id].get('calibration_status')}
        fitted = any(t.get('residual') for t in chain)
        if model_id == 'betse-tissue':
            return 'synthesized', 'Solver-generated planar tissue; the source records that its initial random geometry was not seeded reproducibly and is retained by hash.', quoted
        if identity in self.receipts:
            return 'derived', 'Built in this repository from retained inputs; the build receipt hashes every geometry, measurement and code file it consumed.', quoted
        if evidence_kind in ('synthesized_layer',) or dataset['id'] == 'ihm-synthesis':
            return 'synthesized', 'No source geometry of its own; extent follows an explicit thickness prior over a parent surface.', quoted

        if fitted:
            return 'transferred', 'A fitted cross-source transform places this geometry in the canonical body; its residual is recorded above.', quoted
        if dataset['id'] == 'bodyparts3d':
            return 'measured', "Acquired BodyParts3D mesh for the atlas specimen, status 'acquired'; only the recorded exact rigid/scale display transform is applied. The upstream imaging protocol is not retained here.", quoted
        if dataset['id'].startswith('vascular-'):
            return 'measured', 'Case-specific vascular surface segmented from a separate imaging subject; retained case metadata records species, sex, age and disease. Not registered to this body.', quoted
        if dataset['id'] == 'z-anatomy':
            return 'derived', 'Authored atlas surface that the source itself records as partly derived from BodyParts3D and not an additional measured human.', quoted
        if dataset['id'] == 'opensim-rajagopal':
            return 'derived', 'Published musculoskeletal model geometry in its source default pose; no acquisition of a specimen is retained in this repository.', quoted
        if dataset['id'] == 'published-lymphatic-network':
            return 'derived', 'Published PlasticBoy-derived structural graph; the source declares it is not measured anatomy.', quoted
        return 'derived', 'Produced from other retained geometry by a recorded build step.', quoted

    def role(self, structure, entity):
        if entity and entity.get('role'):
            return {'role': entity['role'], 'role_source': 'canonical entity record'}
        try:
            inferred = physical_role(structure.get('name') or '', structure.get('system') or '')
        except Exception:
            inferred = None
        return {'role': inferred, 'role_source': 'inferred by ihm.assembly.anatomy.physical_role from name and system'
                if inferred else 'no role recorded and none inferable'}

    def frame_relation(self, structure, chain):
        if structure['model_id'] == 'ihm-body':
            return 'canonical' if not any(t.get('residual') for t in chain) else 'registered_into_canonical'
        return 'source_frame_display_placement'

    def record(self, structure, display_present=True):
        entity = self.entities.get(structure.get('canonical_entity_id') or structure['id']) if structure['model_id'] == 'ihm-body' else None
        if structure['id'] in self.receipts:
            entity = None
        dataset = self.datasets_for(structure)
        chain, chain_note = self.transforms(structure, entity)
        geometry = self.geometry(structure, entity)
        tier, basis, quoted = self.tier(structure, entity, chain)
        record = {'schema': SCHEMA_ID, 'structure_id': structure['id'], 'model_id': structure['model_id'],
                  'canonical_entity_id': structure.get('canonical_entity_id'), 'name': structure.get('name'),
                  'system': structure.get('system'), **self.role(structure, entity),
                  'evidence_kind': structure.get('evidence_kind') or (entity or {}).get('evidence_kind'),
                  'dataset': dataset, 'source_file': self.source_file(structure, entity), 'build': self.build_of(structure),
                  'geometry': geometry, 'transforms': chain, 'transform_note': chain_note,
                  'frame_relation': self.frame_relation(structure, chain),
                  'tier': tier, 'tier_basis': basis, 'tier_evidence': {k: v for k, v in quoted.items() if v},
                  'assumptions': [{'id': a, 'statement': next((x['statement'] for x in self.anatomy['assumption_ledger'] if x['id'] == a), None)}
                                  for a in (structure.get('assumptions') or (entity or {}).get('assumptions') or [])],
                  'derived_artifacts': self.builds.lookup(geometry['path'], geometry['sha256']),
                  'display_present': display_present}
        record['completeness'] = completeness(record)
        return record


def dotted(record, key):
    value = record
    for part in key.split('.'):
        if not isinstance(value, dict):
            return None
        value = value.get(part)
    return value


def completeness(record):
    present, missing = [], []
    for key in REQUIRED:
        value = dotted(record, key)
        (present if value not in (None, '', [], {}) else missing).append(key)
    hashes = [record['source_file'].get('sha256_verified'), record['geometry'].get('sha256_verified')]
    return {'required_present': len(present), 'required_total': len(REQUIRED), 'missing': missing,
            'answerable': not ({'dataset.id', 'source_file.sha256', 'build.script', 'geometry.sha256'} & set(missing)),
            'hashes_verified': sum(1 for h in hashes if h is True), 'hashes_failed': sum(1 for h in hashes if h is False),
            'hashes_unchecked': sum(1 for h in hashes if h is None)}


# --- audit of the pre-existing state ------------------------------------
BEFORE_FIELDS = {
    'source_dataset': lambda s, e: bool((s.get('source') or {}).get('label')) or bool((e or {}).get('provenance', {}).get('source')),
    'source_file_path': lambda s, e: bool((s.get('source') or {}).get('path') or s.get('geometry_source')
                                          or ((e or {}).get('provenance', {}).get('files') or [{}])[0].get('path')),
    'source_file_hash': lambda s, e: bool((s.get('source') or {}).get('sha256') or s.get('geometry_source_sha256')
                                          or ((e or {}).get('provenance', {}).get('files') or [{}])[0].get('sha256')),
    'build_script': lambda s, e: False,
    'build_commit': lambda s, e: False,
    'transform_chain': lambda s, e: bool((e or {}).get('provenance', {}).get('source_to_canonical')
                                         or (e or {}).get('provenance', {}).get('registration_id')),
    'transform_residual': lambda s, e: bool((e or {}).get('provenance', {}).get('registration_id')),
    'license': lambda s, e: bool((s.get('source') or {}).get('license') or (e or {}).get('provenance', {}).get('license')),
    'evidence_tier': lambda s, e: bool(s.get('evidence_kind') or (e or {}).get('evidence_kind'))}

AFTER_FIELDS = {
    'source_dataset': lambda r: bool(r['dataset']['id']),
    'source_file_path': lambda r: bool(r['source_file']['path']),
    'source_file_hash': lambda r: bool(r['source_file']['sha256']),
    'build_script': lambda r: bool(r['build']['script']),
    'build_commit': lambda r: bool(r['build']['commit']),
    'transform_chain': lambda r: bool(r['transforms']),
    'transform_residual': lambda r: any(t.get('residual') for t in r['transforms']),
    'license': lambda r: bool(r['dataset']['license']),
    'evidence_tier': lambda r: bool(r['tier'])}


def audit(provenance, records):
    structures = provenance.manifest['structures']
    entities = provenance.entities
    before, after = Counter(), Counter()
    per_model = defaultdict(lambda: {'count': 0, 'before': Counter(), 'after': Counter()})
    per_system, per_role = defaultdict(lambda: {'count': 0, 'unanswerable_before': 0, 'unanswerable_after': 0}), \
                           defaultdict(lambda: {'count': 0, 'unanswerable_before': 0, 'unanswerable_after': 0})
    for structure in structures:
        entity = entities.get(structure.get('canonical_entity_id') or structure['id']) if structure['model_id'] == 'ihm-body' else None
        record = records[structure['id']]
        model = per_model[structure['model_id']]
        model['count'] += 1
        answered_before = True
        for field, test in BEFORE_FIELDS.items():
            ok = bool(test(structure, entity))
            before[field] += ok
            model['before'][field] += ok
            if field in ('source_dataset', 'source_file_hash', 'build_script'):
                answered_before &= ok
        for field, test in AFTER_FIELDS.items():
            ok = bool(test(record))
            after[field] += ok
            model['after'][field] += ok
        system, role = structure.get('system') or 'unassigned', record['role'] or 'unassigned'
        for bucket, key in ((per_system, system), (per_role, role)):
            bucket[key]['count'] += 1
            bucket[key]['unanswerable_before'] += not answered_before
            bucket[key]['unanswerable_after'] += not record['completeness']['answerable']
    total = len(structures)
    return {'schema': 'ihm.structure-provenance-audit.v1', 'structures_audited': total,
            'canonical_entities': len(entities),
            'display_structures_not_canonical': total - sum(1 for s in structures if s['model_id'] == 'ihm-body'),
            'field_coverage': {field: {'before': before[field], 'after': after[field], 'total': total,
                                       'before_percent': round(100*before[field]/total, 2), 'after_percent': round(100*after[field]/total, 2)}
                               for field in BEFORE_FIELDS},
            'per_model': {model: {'count': data['count'],
                                  'before': {f: data['before'][f] for f in BEFORE_FIELDS},
                                  'after': {f: data['after'][f] for f in AFTER_FIELDS}}
                          for model, data in sorted(per_model.items())},
            'per_system': dict(sorted(per_system.items())), 'per_role': dict(sorted(per_role.items())),
            'unanswerable_before': sum(v['unanswerable_before'] for v in per_system.values()),
            'unanswerable_after': sum(v['unanswerable_after'] for v in per_system.values()),
            'tier_counts': dict(Counter(r['tier'] for r in records.values())),
            'frame_relation_counts': dict(Counter(r['frame_relation'] for r in records.values())),
            'hash_verification': {'source_file_verified': sum(1 for r in records.values() if r['source_file'].get('sha256_verified') is True),
                                  'source_file_failed': sum(1 for r in records.values() if r['source_file'].get('sha256_verified') is False),
                                  'source_file_unchecked': sum(1 for r in records.values() if r['source_file'].get('sha256_verified') is None),
                                  'geometry_verified': sum(1 for r in records.values() if r['geometry'].get('sha256_verified') is True),
                                  'geometry_failed': sum(1 for r in records.values() if r['geometry'].get('sha256_verified') is False),
                                  'geometry_unchecked': sum(1 for r in records.values() if r['geometry'].get('sha256_verified') is None)},
            'derived_artifact_coverage': {'structures_with_at_least_one': sum(1 for r in records.values() if r['derived_artifacts']),
                                          'total_references': sum(len(r['derived_artifacts']) for r in records.values())},
            'unfilled': {'license_absent': sorted({r['dataset']['id'] for r in records.values() if not r['dataset']['license']}),
                         'source_file_path_absent': sum(1 for r in records.values() if not r['source_file']['path']),
                         'transform_absent': sum(1 for r in records.values() if not r['transforms']),
                         'residual_absent': sum(1 for r in records.values() if not any(t.get('residual') for t in r['transforms']))}}


SCHEMA = {'schema': SCHEMA_ID,
          'intent': 'One provenance record shape for every structure regardless of origin. No structure is gated, ranked or visually '
                    'discriminated by this record; the tier states what kind of evidence the geometry rests on, not whether it may be shown.',
          'tiers': TIERS,
          'fields': {
              'structure_id': 'display structure id; equals the canonical entity id for ihm-body structures',
              'canonical_entity_id': 'canonical entity this display structure resolves to, or null',
              'dataset': 'originating dataset: id, label, version, revision, url, specimen, units, frame, attribution, acquisition_status, '
                         'license, license_url, license_evidence (retained file path and hash) or license_absent_reason',
              'source_file': 'the specific upstream file: path, sha256, bytes, sha256_verified (recomputed against the working tree)',
              'build': 'script, script_sha256, commit, commit_covers_working_tree, uncommitted_changes',
              'geometry': 'the current geometry file: path, sha256, sha256_verified, representation, frame, units, vertex/face counts',
              'transforms': 'ordered chain applied to reach the current geometry; each entry carries residual {metric,value,method} '
                            'where a fit was measured, or residual null with residual_reason where the transform was exact',
              'tier': 'measured | transferred | derived | synthesized',
              'tier_basis': 'the rule that assigned the tier',
              'tier_evidence': 'the verbatim retained strings that the rule read',
              'frame_relation': 'canonical | registered_into_canonical | source_frame_display_placement',
              'assumptions': 'ids and statements from the canonical assumption ledger',
              'derived_artifacts': 'retained builds that consumed this geometry, matched by path and/or geometry sha256',
              'completeness': 'required-field count, missing list, answerable flag, hash verification tallies'},
          'required_fields': REQUIRED,
          'honesty_rule': 'An empty field is emitted where the repository holds no evidence. No licence, hash, commit or residual is inferred.'}


def build(verify=True, out=OUT, promoted=False):
    out = Path(out)
    out = out if out.is_absolute() else (ROOT/out)
    relative_out = str(out.resolve().relative_to(ROOT))
    provenance = Provenance(verify, out)
    records = {}
    for structure in provenance.manifest['structures']:
        records[structure['id']] = provenance.record(structure)
    missing_entities = [i for i in provenance.entities if i not in records]
    for identity in missing_entities:
        entity = provenance.entities[identity]
        records[identity] = provenance.record(
            {'id': identity, 'model_id': 'ihm-body', 'name': entity['name'], 'system': entity['system'], 'kind': 'mesh',
             'geometry_path': entity['reference_geometry']['path'], 'geometry_sha256': entity['reference_geometry']['sha256'],
             'evidence_kind': entity['evidence_kind'], 'canonical_entity_id': identity,
             'assumptions': entity.get('assumptions', []), 'source': {},
             'display_geometry': {'source_vertices': entity.get('source_vertex_count'), 'source_faces': entity.get('source_face_count')}},
            display_present=False)
    out.mkdir(parents=True, exist_ok=True)
    (out/'records').mkdir(exist_ok=True)
    for stale in (out/'records').glob('*.json'):
        if stale.stem not in records:
            stale.unlink()
    record_hashes = {}
    for identity, record in records.items():
        record_hashes[identity] = write_json(out/'records'/f'{identity}.json', record)
    report = audit(provenance, records)
    schema_sha = write_json(out/'schema.json', SCHEMA)
    audit_sha = write_json(out/'audit.json', report)
    index = {'schema': SCHEMA_ID, 'schema_path': relative_out+'/schema.json',
             'schema_sha256': schema_sha, 'records_directory': relative_out+'/records',
             'record_count': len(records), 'canonical_entities_without_display_structure': missing_entities,
             'display_structures': len(provenance.manifest['structures']),
             'datasets': provenance.datasets, 'tiers': TIERS,
             'counts': {'per_model': dict(Counter(r['model_id'] for r in records.values())),
                        'per_tier': report['tier_counts'], 'per_frame_relation': report['frame_relation_counts']},
             'audit_summary': {k: report[k] for k in ('field_coverage', 'unanswerable_before', 'unanswerable_after',
                                                      'hash_verification', 'derived_artifact_coverage', 'unfilled')},
             'record_sha256': record_hashes,
             'canonical_assets_modified': promoted,
             'promotion_performed': promoted,
             'display_gating': 'none; this record does not select, hide, colour or order any structure'}
    index_sha = write_json(out/'index.json', index)
    resolution_sha = write_json(out/'hash-resolution.json', {
        'purpose': 'files located by content hash where the retained metadata recorded a hash but no path',
        'search_roots': HASH_SEARCH_ROOTS,
        'resolved': {digest: relative for digest, relative in sorted(provenance.resolved.items())
                     if digest in {r['source_file'].get('sha256') for r in records.values()}}})
    inputs = {p: provenance.file_hash(p) for p in
              ['data/derived/app/manifest.json', 'data/derived/canonical/anatomy.json',
               'data/derived/anatomy/bodyparts3d_index.json', 'data/derived/anatomy/extended/source_index.json']}
    manifest = {'schema': 'ihm.structure-provenance-candidate.v1', 'python': sys.version,
                'builder': provenance.script(str(Path(__file__).resolve().relative_to(ROOT))),
                'inputs_sha256': inputs, 'consumed_build_manifests': provenance.builds.manifests,
                'outputs_sha256': {'schema.json': schema_sha, 'index.json': index_sha, 'audit.json': audit_sha,
                                   'hash-resolution.json': resolution_sha,
                                   'records_sha256_of_index_field': sha256_json(record_hashes)},
                'working_tree': working_tree(),
                'record_count': len(records), 'hash_verification_enabled': verify,
                'canonical_assets_modified': promoted, 'promotion_performed': promoted,
                'limitations': [
                    'A candidate for review. No canonical file is written and no structure is promoted.'
                    if not promoted else
                    'Promoted. This index is the traceability record the canonical body manifest points at; '
                    'it indexes the canonical entity set as it stands and is rebuilt whenever that set changes.',
                    'Licence is emitted only where this repository retains a declaration or a licence file; opensim-rajagopal and betse have none and are left empty.',
                    'Residuals are copied from the registration reports that measured them. Structures placed by an exact recorded transform carry residual null with a stated reason, not a zero.',
                    'The tier states the kind of evidence behind the geometry. It does not rank structures and nothing in the viewer gates on it.',
                    'derived_artifacts is a retained-build cross-reference matched by geometry path and hash; it does not assert that a build is current.']}
    manifest_sha = write_json(out/'manifest.json', manifest)
    print(json.dumps({'records': len(records), 'index_sha256': index_sha, 'manifest_sha256': manifest_sha,
                      'unanswerable_before': report['unanswerable_before'], 'unanswerable_after': report['unanswerable_after'],
                      'tiers': report['tier_counts'], 'hash_verification': report['hash_verification']}, indent=1))
    return index, report


def self_test(out=OUT):
    """Exercise the shape rules on synthetic inputs inside a temporary directory only."""
    failures = []
    with tempfile.TemporaryDirectory() as temporary:
        probe = Path(temporary)/'probe.bin'
        probe.write_bytes(b'ihm')
        if sha256_file(probe) != hashlib.sha256(b'ihm').hexdigest():
            failures.append('sha256_file disagrees with hashlib')
        if write_json(Path(temporary)/'a/b.json', {'x': 1}) != sha256_file(Path(temporary)/'a/b.json'):
            failures.append('write_json returned a hash that does not match the file')
    record = {'dataset': {'id': 'd', 'label': 'L', 'license': 'CC'}, 'source_file': {'path': 'p', 'sha256': 'h', 'sha256_verified': True},
              'build': {'script': 's', 'commit': 'c'}, 'geometry': {'path': 'g', 'sha256': 'k', 'sha256_verified': True},
              'transforms': [{'kind': 'x', 'residual': None}], 'tier': 'measured'}
    full = completeness(record)
    if full['missing'] or not full['answerable'] or full['hashes_verified'] != 2:
        failures.append(f'complete record misjudged: {full}')
    empty = completeness({'dataset': {}, 'source_file': {}, 'build': {}, 'geometry': {}, 'transforms': [], 'tier': None})
    if empty['answerable'] or len(empty['missing']) != len(REQUIRED):
        failures.append(f'empty record misjudged: {empty}')
    if set(TIERS) != {'measured', 'transferred', 'derived', 'synthesized'}:
        failures.append('tier vocabulary drifted from the four declared tiers')
    if not (out/'index.json').exists():
        print(json.dumps({'self_test': 'shape checks only; no build present yet', 'failures': failures}, indent=1))
        return 1 if failures else 0
    index = read_json(out/'index.json')
    if index['schema'] != SCHEMA_ID:
        failures.append('index schema id drifted')
    checked = 0
    for identity, expected in list(index['record_sha256'].items())[::250]:
        path = out/'records'/f'{identity}.json'
        if not path.exists() or sha256_file(path) != expected:
            failures.append('record hash mismatch: '+identity)
        else:
            record = read_json(path)
            if record['schema'] != SCHEMA_ID or record['structure_id'] != identity:
                failures.append('record identity mismatch: '+identity)
            checked += 1
    manifest = read_json(out/'manifest.json')
    if manifest['outputs_sha256']['index.json'] != sha256_file(out/'index.json'):
        failures.append('manifest does not bind the emitted index')
    for path, expected in manifest['inputs_sha256'].items():
        if expected and sha256_file(ROOT/path) != expected:
            failures.append('input changed since the build: '+path)
    print(json.dumps({'self_test': 'ok' if not failures else 'failed', 'records_sampled': checked,
                      'record_count': index['record_count'], 'failures': failures}, indent=1))
    return 1 if failures else 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--no-verify-hashes', action='store_true')
    parser.add_argument('--output', type=Path, default=OUT,
                        help='where to write. Defaults to the candidate directory; pass '
                             'data/derived/structure-provenance-v1 for the promoted index, which '
                             'leaves the candidate directory intact as its own provenance record.')
    parser.add_argument('--promoted', action='store_true',
                        help='record this build as the promoted index rather than a candidate')
    arguments = parser.parse_args()
    if arguments.self_test:
        raise SystemExit(self_test(arguments.output))
    build(verify=not arguments.no_verify_hashes, out=arguments.output, promoted=arguments.promoted)
