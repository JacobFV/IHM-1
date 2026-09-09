#!/usr/bin/env python3
"""Every LIVE provenance pin still matches the bytes it names.

A live pin is one that production code re-validates when it loads, so a stale one
does not fail here -- it fails much later as an unrelated ValueError deep inside
some other verify module ("Canonical source hash mismatch", "Changed engineered
materialization", "lineage manifest changed"). This module makes that whole class
of rot fail loudly, immediately, and with the name of the pin that moved.

Deliberately NOT covered: point-in-time run receipts under data/derived/<run-id>/
and data/research/<experiment>/, and the provenance/sources/ copies inside
data/models/<bundle>/. Those record the inputs as they were when the run happened;
they are historical evidence and are supposed to go on naming superseded bytes.
"""
import ast
import hashlib
import importlib
import json
from pathlib import Path
import re
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

VARIANTS = ROOT/'data/runtime/physiology/variants'
TERRITORY = 'data/research/engineered_skin_territories/materialization.json'
PREPARED = ROOT/'data/research/configured_regional_skin/prepared_v1'
SUPERSEDED_LUMBAR = 'data/derived/lumbar-muscle-native-lb45uirs'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def document_digest(document):
    return hashlib.sha256(json.dumps(document,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def constant(module_path, name):
    """Read a module-level str constant without importing the module."""
    tree = ast.parse((ROOT/module_path).read_text())
    for node in tree.body:
        if isinstance(node,(ast.Assign,ast.AnnAssign)):
            targets = node.targets if isinstance(node,ast.Assign) else [node.target]
            if any(isinstance(t,ast.Name) and t.id==name for t in targets) and isinstance(node.value,ast.Constant):
                return node.value.value
    raise AssertionError('No module-level constant '+name+' in '+module_path)


class LivePins(unittest.TestCase):
    def assertPinned(self, relative_path, pinned, owner):
        path = ROOT/relative_path
        self.assertTrue(path.is_file(), owner+' pins a file that is gone: '+relative_path)
        self.assertEqual(digest(path), pinned,
                         owner+' pins stale bytes for '+relative_path+'; rebuild it or repoint the pin')

    def test_canonical_body_pins_still_name_the_loaded_runtime_and_sources(self):
        """ihm/assembly/body.py re-hashes each of these and raises on mismatch."""
        body = json.loads((ROOT/'data/derived/canonical/body.json').read_text())
        self.assertEqual(len(body['runtime_sources']), 14)
        for name,record in body['runtime_sources'].items():
            self.assertPinned(record['path'],record['sha256'],'canonical body.json runtime_sources.'+name)
        self.assertEqual(set(body['sources']),{'anatomy','profile','mechanics','brain','respiration','peripheral'})
        for name,record in body['sources'].items():
            self.assertPinned(record['path'],record['sha256'],'canonical body.json sources.'+name)
        self.assertPinned(body['profile']['native_patient_path'],body['profile']['native_patient_sha256'],
                          'canonical body.json profile.native_patient')

    def test_frozen_skin_territory_identity_agrees_across_every_pinner(self):
        """The materialization sha256 is a frozen identity with three live owners."""
        actual = digest(ROOT/TERRITORY)
        owners = {'ihm/assembly/regional_skin_configuration.py':'MATERIALIZATION_SHA',
                  'ihm/assembly/skin_microvascular_patch.py':'TERRITORY_SHA'}
        for module_path,name in owners.items():
            self.assertEqual(constant(module_path,name),actual,
                             module_path+':'+name+' no longer names '+TERRITORY)
        self.assertEqual(constant('ihm/assembly/skin_microvascular_patch.py','TERRITORY_PATH'),TERRITORY)
        # The frozen artifact records the canonical anatomy CONTAINER it was cut from,
        # which a canonical rebuild legitimately supersedes. What must never drift is
        # the per-entity geometry it actually consumes: the live anatomy has to still
        # declare exactly those files at exactly those hashes.
        material = json.loads((ROOT/TERRITORY).read_text())
        anatomy = json.loads((ROOT/'data/derived/canonical/anatomy.json').read_bytes())
        entities = {e['id']:e for e in anatomy['entities']}
        self.assertEqual(len(material['source_receipts']),5)
        for receipt in material['source_receipts']:
            reference = entities[receipt['entity_id']]['reference_geometry']
            self.assertEqual((reference['path'],reference['sha256'],reference['units'],reference['frame']),
                             (receipt['path'],receipt['sha256'],'m',anatomy['frame']['id']),
                             'live anatomy.json moved the geometry frozen for '+receipt['entity_id'])
            self.assertPinned(receipt['path'],receipt['sha256'],TERRITORY+' source_receipts')

    def test_configured_regional_skin_identity_chain_is_unbroken(self):
        """Configuration digest is compiled into the native header the variant pins."""
        configuration = json.loads((PREPARED/'configuration.json').read_text())
        body = {k:v for k,v in configuration.items() if k!='configuration_sha256'}
        self.assertEqual(document_digest(body),configuration['configuration_sha256'],
                         'prepared_v1/configuration.json no longer matches its own configuration_sha256')
        header = (PREPARED/'native_configured_regional_skin.h').read_text()
        compiled = re.search(r'configuration_sha256="([0-9a-f]{64})"',header)
        self.assertIsNotNone(compiled,'native_configured_regional_skin.h has no configuration_sha256 constant')
        self.assertEqual(compiled.group(1),configuration['configuration_sha256'],
                         'the compiled native header names a different configuration than prepared_v1/configuration.json')
        self.assertEqual(configuration['materialization_sha256'],digest(ROOT/TERRITORY),
                         'prepared_v1/configuration.json was cut from a different materialization')

    def test_skin_microvascular_primary_receipts_match_disk(self):
        """materialize_skin_patch re-hashes each receipt before it reads anything."""
        priors = json.loads((ROOT/'data/sources/skin_microvascular_priors.json').read_text())
        self.assertTrue(priors['primary_receipts'])
        for receipt in priors['primary_receipts']:
            self.assertPinned(receipt['path'],receipt['sha256'],'skin_microvascular_priors.json primary_receipts')
            self.assertEqual((ROOT/receipt['path']).stat().st_size,receipt['bytes'])

    def test_native_variant_libraries_and_lineage_match_their_manifests(self):
        """EmbodiedRuntime.from_workspace refuses to start on any drift here."""
        manifests = {d.name:json.loads((d/'manifest.json').read_text())
                     for d in sorted(VARIANTS.iterdir()) if (d/'manifest.json').is_file()}
        self.assertGreaterEqual(len(manifests),28)
        libraries = 0
        for name,manifest in manifests.items():
            library = VARIANTS/name/'libbiogears.so.8.0.0'
            if 'library_sha256' in manifest and library.is_file():
                libraries += 1
                self.assertEqual(digest(library),manifest['library_sha256'],
                                 'built library differs from its manifest for variant '+name)
            parent = manifest.get('parent_variant')
            if parent is None: continue
            self.assertIn(parent,manifests,'variant '+name+' names a parent that is gone: '+str(parent))
            self.assertEqual(digest(VARIANTS/parent/'manifest.json'),manifest['parent_manifest_sha256'],
                             'lineage manifest changed under variant '+name)
            self.assertEqual(manifests[parent].get('library_sha256'),manifest['parent_library_sha256'],
                             'lineage library changed under variant '+name)
        self.assertGreaterEqual(libraries,28)

    def test_mechanical_plant_owner_is_the_module_the_fixtures_patch(self):
        """Guards the ArticulatedBodyPlant -> SelectiveProjectionPlant class of rot.

        A factory fixture that patches a plant the factory no longer imports runs a
        real plant instead of the fake and fails somewhere unrecognizable.
        """
        source = (ROOT/'ihm/assembly/embodied.py').read_text()
        self.assertIn('from .selective_projection import SelectiveProjectionPlant',source)
        self.assertNotIn('from .articulated import ArticulatedBodyPlant',source)
        fixture = (ROOT/'scripts/verify_regional_embodied_factory.py').read_text()
        self.assertIn("patch('ihm.assembly.selective_projection.SelectiveProjectionPlant'",fixture)
        self.assertNotIn("patch('ihm.assembly.articulated.ArticulatedBodyPlant'",fixture)
        from ihm.assembly.selective_projection import SelectiveProjectionPlant
        from ihm.assembly.articulated import ArticulatedBodyPlant
        self.assertTrue(issubclass(SelectiveProjectionPlant,ArticulatedBodyPlant))

    def test_every_patch_target_named_in_a_verify_script_still_resolves(self):
        """A renamed or moved owner has to fail here, not as a silent no-op patch."""
        pattern = re.compile(r"patch\(\s*'((?:ihm|scripts)[A-Za-z0-9_.]*)'")
        targets = {}
        for path in sorted(ROOT.glob('scripts/verify_*.py'))+sorted(ROOT.glob('scripts/test_*.py')):
            for match in pattern.finditer(path.read_text()):
                targets.setdefault(match.group(1),set()).add(path.name)
        self.assertGreaterEqual(len(targets),23)
        for target,users in sorted(targets.items()):
            parts = target.split('.')
            owner = None
            for i in range(len(parts),0,-1):
                try:
                    owner = importlib.import_module('.'.join(parts[:i]));rest = parts[i:];break
                except ImportError: continue
            self.assertIsNotNone(owner,'no importable module in patch target '+target+' used by '+str(sorted(users)))
            for attribute in rest:
                self.assertTrue(hasattr(owner,attribute),
                                'patch target '+target+' has no '+attribute+'; used by '+str(sorted(users)))
                owner = getattr(owner,attribute)

    def test_superseded_lumbar_variant_is_not_a_live_verify_dependency(self):
        """The rematerialized directory replaced it under a byte-identity receipt."""
        for path in sorted(ROOT.glob('scripts/verify_*.py')):
            if path.name==Path(__file__).name: continue
            if SUPERSEDED_LUMBAR in path.read_text():
                self.fail(path.name+' still reads the superseded '+SUPERSEDED_LUMBAR
                          +'; repoint it at data/derived/mechanics/whole_body_lumbar_current/')
        current = ROOT/'data/derived/mechanics/whole_body_lumbar_current'
        receipt = json.loads((current/'rematerialization.json').read_text())
        self.assertEqual(receipt['schema'],'ihm.lumbar-source-rematerialization.v1')
        self.assertEqual(receipt['previous_registration'],SUPERSEDED_LUMBAR+'/variant/registration.json')
        self.assertTrue(receipt['donor_fragment_byte_identical'])
        self.assertTrue(receipt['model_catalog_and_insertion_byte_identical'])
        self.assertTrue(receipt['current_audit_recomputed'])
        self.assertFalse(receipt['new_dynamic_acceptance_claimed'])
        # The receipt exists precisely because these canonical containers were rebuilt;
        # its "current" side has to still be the canonical file on disk today.
        for path,change in receipt['changed_coverage_sources'].items():
            self.assertPinned(path,change['current_sha256'],'lumbar rematerialization.json current_sha256')
            self.assertNotEqual(change['previous_sha256'],change['current_sha256'])
        for name in ('registration.json','catalog.json','audit.json','subject_with_lumbar.osim'):
            self.assertTrue((current/name).is_file(),'rematerialized lumbar variant is missing '+name)


if __name__=='__main__':unittest.main()
