"""Fail-closed static peripheral port audit; never executes a native or neural model.

A receptor ID here identifies an observation port, not a measured receptor/axon.
Coverage means implemented wiring capability, not an active session or validation.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from .sensorimotor_catalog import whole_body_effector_catalog

DEFAULT_REGISTRATION = 'data/derived/mechanics/whole_body_arm26_v2/registration.json'
BRAIN_PATH = 'data/derived/canonical/brain.json'
REFLEX_PATH = 'data/raw/sensorimotor/geyer_herr_2010.pdf'
REFLEX_SHA256 = '8a60efec9ceba5120c1679610f676ab297e8dce887ed4ae6f0890f7e17e13d94'
REFLEX_EFFECTORS = tuple(f'{name}_{side}' for side in ('r', 'l')
                        for name in ('tibant', 'soleus', 'gasmed', 'gaslat'))


def _verified(root, relative, expected=None):
    if not isinstance(relative, str) or Path(relative).is_absolute():
        raise ValueError('Audit sources must be root-relative paths')
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError('Audit source escapes root')
    try:
        with path.open('rb') as source:
            digest = hashlib.file_digest(source, 'sha256').hexdigest()
    except OSError as exc:
        raise ValueError(f'Missing audit source: {relative}') from exc
    if expected is not None and expected != digest:
        raise ValueError(f'Audit source hash mismatch: {relative}')
    return path, digest


def build_peripheral_coverage(root, registration=DEFAULT_REGISTRATION, *, controller_catalog=None):
    """Return JSON-safe, exact-plant coverage; optionally check a runtime catalog.

    ``registration`` is a root-relative manifest path or a manifest dictionary.
    No ODE, native plant, receptor simulation, or large geometry is loaded.
    """
    root = Path(root).resolve()
    provenance = {}
    if isinstance(registration, dict):
        manifest = deepcopy(registration)
    else:
        path, digest = _verified(root, str(registration))
        provenance[str(registration)] = digest
        manifest = json.loads(path.read_text())
    # Check paths before the catalog reader can attempt filesystem access.
    for name in ('model', 'catalog'):
        relative = manifest[f'{name}_path']
        _, provenance[relative] = _verified(root, relative, manifest[f'{name}_sha256'])
    for relative, expected in manifest['sources'].items():
        _, provenance[relative] = _verified(root, relative, expected)
    catalog = whole_body_effector_catalog(root, manifest)
    model = ET.parse(root / manifest['model_path']).getroot()
    bodies = {b.get('name') for b in model.findall('.//BodySet/objects/*')}
    muscles = {m.get('name'): m for m in model.findall('.//ForceSet/objects/*') if 'Muscle' in m.tag}
    brain_path, provenance[BRAIN_PATH] = _verified(root, BRAIN_PATH)
    brain = json.loads(brain_path.read_text())
    regions = [n['id'] for n in brain['nodes']]
    if len(set(regions)) != len(regions):
        raise ValueError('Duplicate brain region ID')
    for edge in brain['edges']:
        if edge['source'] not in regions or edge['target'] not in regions:
            raise ValueError('Brain edge references nonexistent population')
    # Stream hashes of preserved sources; never load donor meshes into memory.
    for source in brain['sources']:
        relative = source['preserved_path']
        _, provenance[relative] = _verified(root, relative, source['sha256'])
    if not any(s.get('role') == 'executed_neural_law' for s in brain['sources']):
        raise ValueError('Missing pinned IBM neural law')
    _, provenance[REFLEX_PATH] = _verified(root, REFLEX_PATH, REFLEX_SHA256)
    for relative in ('ihm/assembly/sensorimotor.py', 'ihm/assembly/sensorimotor_catalog.py',
                     'ihm/assembly/brain.py', 'ihm/assembly/peripheral_coverage.py'):
        _, provenance[relative] = _verified(root, relative)
    expected = {r['id']: r for r in catalog}
    if not set(REFLEX_EFFECTORS) <= expected.keys():
        raise ValueError('Registered plant lacks required source ankle effectors')
    if controller_catalog is not None:
        if not isinstance(controller_catalog, (list, tuple)):
            raise ValueError('Controller catalog must be a list of registered rows')
        if len(controller_catalog) != len(catalog) or any(not isinstance(r, dict) for r in controller_catalog):
            raise ValueError('Controller catalog differs from registered plant')
        supplied = {r.get('id'): r for r in controller_catalog}
        if supplied != expected:
            raise ValueError('Controller catalog identity, assignment or provenance differs from registered plant')
    effectors, receptors = [], []
    for row in catalog:
        key, side = row['id'], row['side']
        if side not in ('r', 'l') or not key.endswith('_' + side):
            raise ValueError(f'Invalid effector laterality: {key}')
        hemisphere = 'lh' if side == 'r' else 'rh'
        for field, suffix in [('motor_region', 'precentral'), ('sensory_region', 'postcentral')]:
            if row[field] not in regions or row[field] != f'brain-{hemisphere}-{suffix}':
                raise ValueError(f'Invalid cortical assignment: {key}')
        attachments = {p.text.rsplit('/', 1)[-1] for p in muscles[key].findall('.//PathPointSet/objects/*/socket_parent_frame')}
        if not attachments or attachments != set(row['attachment_bodies']) or attachments - bodies:
            raise ValueError(f'Attachment bodies differ from actual plant: {key}')
        receptor_ids = []
        for kind, fields, normalization in [
            ('length', ['fiber_length_m', 'fiber_length_proxy_m'], 'optimal_fiber_length_m'),
            ('force', ['tendon_force_n'], 'max_isometric_force_n')]:
            receptor_id = f'muscle-observation/{key}/{kind}'
            receptor_ids.append(receptor_id)
            receptors.append({'id': receptor_id, 'effector_id': key, 'kind': kind,
                'status': 'engineered_observation_port', 'input_fields': fields,
                'normalization_field': normalization, 'requires_sensor_basis': True,
                'measured_receptor_identity': False,
                'basis': 'native mechanical observation; fiber proxy must remain explicitly labeled'})
        reflex = {'status': 'unsupported', 'input_effectors': [], 'reason': 'No source spinal reflex primitive implemented for this effector'}
        if key in REFLEX_EFFECTORS:
            inputs = [key, 'soleus_' + side] if key.startswith('tibant') else (
                ['gasmed_' + side, 'gaslat_' + side] if key.startswith('gas') else [key])
            reflex = {'status': 'source_derived_primitive', 'input_effectors': inputs,
                'contact_receptor_id': f'foot-contact/{side}', 'source_path': REFLEX_PATH,
                'source_sha256': REFLEX_SHA256,
                'basis': 'Geyer/Herr ankle primitive transfer; GAS heads explicitly lumped; not identified axonal wiring'}
        effectors.append({'id': key, 'side': side, 'body_group': row['body_group'],
            'attachment_bodies': sorted(attachments), 'source_path': row['source_path'],
            'source_sha256': row['source_sha256'], 'receptor_ids': receptor_ids,
            'afferent': {'status': 'engineered', 'source_receptor_ids': receptor_ids,
                'target_region_id': row['sensory_region'], 'basis': row['assignment_basis'],
                'implementation': 'SensorimotorController._step', 'output_unit': 'Hz'},
            'efferent': {'status': 'engineered', 'source_region_id': row['motor_region'],
                'target_effector_id': key, 'basis': 'regional-rate gain and requested-effector-gated drive; no identified motor recruitment',
                'output_field': 'motor_excitations', 'output_unit': 'dimensionless excitation',
                'activation_owner': 'mechanical_plant'}, 'spinal_reflex': reflex})
    for side in ('r', 'l'):
        receptors.append({'id': f'foot-contact/{side}', 'side': side,
            'kind': 'contact', 'status': 'engineered_observation_port',
            'input_field': f'foot_contact_force_n.{side}', 'unit': 'N',
            'measured_receptor_identity': False, 'basis': 'actual foot normal load used for stance gate; no cutaneous receptor population'})
    gaps = [
        ('spindle_and_tendon_afferents', 'No explicit spindle, Golgi tendon organ, Ia/II/Ib axon or motor-unit populations; length/force ports are proxies.'),
        ('peripheral_nerve_anatomy', 'No measured nerve trunks, roots, dermatomes, axonal conduction paths or neuromuscular junctions.'),
        ('cutaneous_to_cortex', 'Regional rapid/slow IBM touch exists separately; signed donor response has no identified conversion to cortical Hz and is not a whole-body receptor census.'),
        ('non_ankle_reflexes', 'No spinal reflex primitives for the other 84 registered effectors.'),
        ('cranial_and_special_senses', 'No audited visual, auditory, vestibular, olfactory, gustatory, ocular, facial or bulbar peripheral routes.'),
        ('visceral_and_autonomic', 'Native systemic autonomic ownership is separate; no audited IBM visceral afferent or autonomic efferent wiring.'),
        ('absent_effectors', 'Hand/finger, neck, facial, respiratory and independent scapular muscle effectors are not added by this 92-muscle catalog.'),
        ('autonomous_motor_policy', 'No identified recruitment, calibrated descending policy or autonomous walking controller.')]
    report = {'schema': 'ihm.peripheral-coverage.v1', 'plant_model_path': manifest['model_path'],
        'plant_model_sha256': manifest['model_sha256'], 'catalog_sha256': manifest['catalog_sha256'],
        'biological_validation': False, 'runtime_active': False,
        'coverage_semantics': 'Static implemented-port capability for exact registered plant; does not establish native execution or biological coverage',
        'registration_native_verified': manifest.get('native_verified', False),
        'brain_region_ids': regions, 'effectors': effectors, 'receptors': receptors,
        'unsupported_paths': [{'id': key, 'status': 'unsupported', 'reason': reason} for key, reason in gaps],
        'source_hashes': provenance,
        'summary': {'effector_count': len(effectors), 'receptor_port_count': len(receptors),
            'engineered_afferent_effector_count': len(effectors), 'engineered_efferent_effector_count': len(effectors),
            'source_reflex_effector_count': len(REFLEX_EFFECTORS),
            'anatomically_measured_connection_count': 0}}
    report['audit_sha256'] = hashlib.sha256(json.dumps(report, sort_keys=True, allow_nan=False).encode()).hexdigest()
    return report


def write_peripheral_coverage(root, output_path, registration=DEFAULT_REGISTRATION, *, controller_catalog=None):
    """Validate completely before writing a deterministic JSON artifact."""
    report = build_peripheral_coverage(root, registration, controller_catalog=controller_catalog)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return report


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--registration', default=DEFAULT_REGISTRATION)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = write_peripheral_coverage(args.root, args.output, args.registration)
    print(json.dumps(result['summary'], sort_keys=True))
