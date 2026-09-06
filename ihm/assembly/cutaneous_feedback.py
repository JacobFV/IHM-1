"""Causal cutaneous stimuli through pinned IBM operators and explicit cortical priors.

Sites, contact areas, foundation stiffness and cortical recruitment are supplied
by the caller. This adapter owns receptor state only, never the shared brain clock.
"""
from ihm.brain.active_source import DEFAULT_SOURCE,resolve_source
from copy import deepcopy
import hashlib
import inspect
import json
import math
from pathlib import Path

import numpy as np
from ihm.brain.causal import CausalIBM
from ihm.brain.ibm_backend import IBMBackend
from ihm.brain.port_mapping import pressure_to_indentation_um
from .brain import verify_sources


def _number(value, label, minimum=-1e12, maximum=1e12):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not minimum <= value <= maximum:
        raise ValueError(f'Invalid {label}')
    return float(value)


def _vector(value, label):
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f'Expected three-vector {label}')
    return [_number(v, label) for v in value]


class _ThermalBackend:
    """Select actual pinned thermal transfer/prior medians for existing ZOH adapter."""
    def __init__(self, backend):
        from ibm.forge.priors import median_of
        self.identity = backend.identity
        item = next(i for i in backend.registry.implementations.values()
                    if i.name == 'thermoreceptor_static_dynamic')
        self.fn = item.transfer
        accepted = inspect.signature(self.fn).parameters
        self.theta = {k: median_of(v) for k, v in item.params.items() if k in accepted}

    def _transfer_spec(self, kind):
        if kind != 'slow':
            raise ValueError('Thermal adapter requires static/dynamic transfer')
        return self.fn, self.theta, 'thermoreceptor_static_dynamic'

    def transfer(self, kind, basis, parameters=None):
        self._transfer_spec(kind)
        return self.fn(basis, **{**self.theta, **(parameters or {})})


class CutaneousFeedback:
    """Material-site force and skin-temperature samples to causal sensory input.

    Unknown area disables force-to-pressure transduction and preserves raw force.
    Direct native indentation bypasses foundation stiffness and requires a stable
    material manifest/quadrature/triangle receipt plus explicit deformation basis.
    A temperature sample is a caller-supplied native skin temperature in degC,
    applied to each explicit site as a declared spatially uniform transfer prior.
    """
    def __init__(self, root, *, sites, recruitment_hz_per_response, delay_s=0., source_pin=DEFAULT_SOURCE):
        self.root = Path(root).resolve()
        self.recruitment = _number(recruitment_hz_per_response, 'recruitment gain', 0., 1e6)
        self.delay_s = _number(delay_s, 'conduction delay', 0., 10.)
        brain_path = self.root / 'data/derived/canonical/brain.json'
        brain = json.loads(brain_path.read_text())
        verify_sources(brain, self.root)
        regions = {n['id'] for n in brain['nodes']}
        if not isinstance(sites, (list, tuple)) or not 1 <= len(sites) <= 64:
            raise ValueError('Supply 1 to 64 explicitly registered material sites')
        self.sites = {}
        for raw in sites:
            if not isinstance(raw, dict):
                raise ValueError('Invalid material site')
            site = deepcopy(raw)
            key = site.get('id')
            if not isinstance(key, str) or not key or key in self.sites:
                raise ValueError('Duplicate or invalid site ID')
            if not isinstance(site.get('support_basis'), str) or not site['support_basis'].strip():
                raise ValueError('Explicit material registration/support basis required')
            site['position_m'] = _vector(site.get('position_m'), 'site position')
            site['normal'] = _vector(site.get('normal'), 'outward surface normal')
            if not math.isclose(math.sqrt(sum(v*v for v in site['normal'])), 1., rel_tol=1e-6, abs_tol=1e-6):
                raise ValueError('Surface normal must be an explicit unit outward vector')
            if site.get('sensory_region') not in regions:
                raise ValueError('Sensory region absent from canonical brain')
            area = site.get('contact_area_m2')
            site['contact_area_m2'] = None if area is None else _number(area, 'contact area', 1e-12, 10.)
            mode = site.setdefault('mechanical_input', 'force_foundation')
            if mode not in ('force_foundation', 'native_indentation'):
                raise ValueError('Unknown mechanical receptor input mode')
            if mode == 'force_foundation':
                site['stiffness_pa_per_m'] = _number(site.get('stiffness_pa_per_m'), 'foundation stiffness Pa/m', 1e-12, 1e15)
            else:
                if 'stiffness_pa_per_m' in site:
                    raise ValueError('Direct native indentation must not impose another foundation stiffness')
                identity = site.get('material_identity')
                if not isinstance(identity, dict) or set(identity) != {'manifest_sha256', 'quadrature_index', 'triangle_index'}:
                    raise ValueError('Native indentation requires stable manifest/quadrature/triangle identity')
                digest = identity['manifest_sha256']
                if not isinstance(digest, str) or len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
                    raise ValueError('Invalid native material manifest hash')
                if any(type(identity[k]) is not int or identity[k] < 0 for k in ('quadrature_index', 'triangle_index')):
                    raise ValueError('Invalid native quadrature/triangle index')
                for field in ('indentation_basis', 'area_basis'):
                    if not isinstance(site.get(field), str) or not site[field].strip():
                        raise ValueError(f'Explicit native {field} required')
            site['reference_temperature_C'] = _number(site.get('reference_temperature_C'), 'reference temperature', -100., 100.)
            self.sites[key] = site
        source_pin=resolve_source(self.root,source_pin)
        backend = IBMBackend(self.root / 'data/derived/canonical/ibm-backend',source_pin=None) if source_pin is None else IBMBackend(source_pin=source_pin)
        thermal = _ThermalBackend(backend)
        self.models, self.initial = {}, {}
        for key, site in self.sites.items():
            models = {kind: CausalIBM(backend, sites_m=[site['position_m']], kind=kind, delay_s=self.delay_s)
                      for kind in ('rapid', 'slow')}
            model = CausalIBM(thermal, sites_m=[site['position_m']], kind='slow', delay_s=self.delay_s)
            # The existing slow adapter's algebra is exact for the source thermal
            # transfer, but its default mechanoreceptor unit label must change.
            model.audit.update(input_units='degC change from explicit reference',
                               output_units='signed donor thermal response; uncalibrated linear transient',
                               thermal_reference_basis='caller prior; no nonlinear warm/cold tuning')
            model._identity = hashlib.sha256(json.dumps(model.audit, sort_keys=True).encode()).hexdigest()
            models['thermal'] = model
            self.models[key] = models
            self.initial[key] = {kind: m.checkpoint() for kind, m in models.items()}
        first = next(iter(self.models.values()))
        self.audit = {'schema': 'ihm.cutaneous-materialization.v1',
            'source_pin': None if source_pin is None else source_pin.to_dict(),
            'package_sha256': backend.identity['package_sha256'],
            'sites': list(deepcopy(self.sites).values()),
            'receptors': {kind: deepcopy(m.audit) for kind, m in first.items()},
            'receptor_audit_sites': 'representative first site; all sites use same verified transfer and priors',
            'moving_support_basis': 'native contacts may carry current point/outward normal together; source LTI support remains registered reference spatial prior, no spatial receptor rematerialization',
            'recruitment_hz_per_response': self.recruitment,
            'recruitment_basis': 'engineered absolute signed-response magnitude pooled into caller-assigned cortical populations; not measured recruitment',
            'mechanical_basis': 'compression=max(0,-force dot outward normal); p=F/A; force_foundation uses explicit p/stiffness prior; native_indentation uses supplied modeled deformation directly with material receipt, no stiffness inversion',
            'thermal_basis': 'native scalar skin temperature transferred uniformly to supplied sites; reference subtraction is a caller prior; source LTI only, no absolute temperature tuning',
            'block_basis': 'engineering ablation: immediately reset blocked receptor state and purge delayed queue; not a physiological local anesthetic model',
            'refractory_model': 'none; source LTI adaptation only, no spike/refractory process',
            'missing_temperature_basis': 'engineering dropout policy: reset thermal state/queue and emit zero; absence is not a reference-temperature sample',
            'source_hashes': {str(p.relative_to(self.root)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in [brain_path, self.root / 'ihm/assembly/cutaneous_feedback.py',
                          self.root / 'ihm/brain/causal.py', self.root / 'ihm/brain/port_mapping.py']}}
        self.identity = hashlib.sha256(json.dumps(self.audit, sort_keys=True, allow_nan=False).encode()).hexdigest()
        self.time_s = 0.

    def _observe(self, observation):
        if not isinstance(observation, dict):
            raise ValueError('Cutaneous observation required')
        time = _number(observation.get('time_s'), 'observation time', 0., 1e9)
        if abs(time - self.time_s) > 1e-9:
            raise ValueError('Cutaneous and physical clocks differ')
        temperature = observation.get('skin_temperature_C')
        if temperature is not None:
            temperature = _number(temperature, 'native skin temperature', -100., 100.)
        contacts = observation.get('contacts')
        if not isinstance(contacts, list):
            raise ValueError('Explicit contacts list required; empty means released')
        forces = {key: [0., 0., 0.] for key in self.sites}
        seen = set()
        geometry = {key: {'point_m': list(site['position_m']), 'normal': list(site['normal']),
                          'moving': False} for key, site in self.sites.items()}
        native_indentations = {key: 0. for key, site in self.sites.items() if site['mechanical_input'] == 'native_indentation'}
        for contact in contacts:
            if not isinstance(contact, dict) or contact.get('id') not in self.sites:
                raise ValueError('Unknown contact material-site ID')
            key = contact['id']
            if key in seen:
                raise ValueError('Duplicate contact ID; aggregate physical forces explicitly')
            seen.add(key)
            forces[key] = _vector(contact.get('force_n'), 'contact force N')
            site = self.sites[key]
            if site['mechanical_input'] == 'native_indentation':
                receipt = contact.get('material_identity')
                if (not isinstance(receipt, dict) or receipt != site['material_identity']
                        or any(type(receipt.get(k)) is not int for k in ('quadrature_index', 'triangle_index'))):
                    raise ValueError('Native contact material identity differs from registered site')
                if contact.get('indentation_basis') != site['indentation_basis']:
                    raise ValueError('Native indentation provenance differs from registration')
                native_indentations[key] = _number(contact.get('indentation_m'), 'native indentation m', 0., 1.)
                if ('point_m' in contact) != ('normal' in contact):
                    raise ValueError('Current native point and outward normal must be supplied together')
                if 'point_m' in contact:
                    point = _vector(contact['point_m'], 'current native point')
                    normal = _vector(contact['normal'], 'current outward surface normal')
                    if not math.isclose(math.sqrt(sum(v*v for v in normal)), 1., rel_tol=1e-6, abs_tol=1e-6):
                        raise ValueError('Current native outward normal must be a unit vector')
                    geometry[key] = {'point_m': point, 'normal': normal, 'moving': True}
            elif 'indentation_m' in contact or 'material_identity' in contact:
                raise ValueError('Native deformation supplied to force-foundation site')
            if site['mechanical_input'] == 'force_foundation' and 'normal' in contact:
                raise ValueError('Static foundation does not accept a moving normal')
            if site['mechanical_input'] == 'force_foundation' and 'point_m' in contact:
                point = _vector(contact['point_m'], 'contact application point')
                if not np.allclose(point, self.sites[key]['position_m'], atol=1e-9, rtol=0.):
                    raise ValueError('Contact point differs from registered material site')
        return forces, temperature, native_indentations, geometry

    def step(self, dt_s, observation, *, sensory_blocks=()):
        dt = _number(dt_s, 'cutaneous interval', 1e-9, .1)
        forces, temperature, native_indentations, geometry = self._observe(observation)
        if not isinstance(sensory_blocks, (list, tuple, set)) or any(not isinstance(k, str) or k not in self.sites for k in sensory_blocks):
            raise ValueError('Unknown sensory block site')
        blocks = set(sensory_blocks)
        saved = self.checkpoint()
        try:
            rows, inputs = [], {}
            for key, site in self.sites.items():
                force = forces[key]
                current = geometry[key]
                compression = max(0., -sum(f*n for f, n in zip(force, current['normal'])))
                area = site['contact_area_m2']
                pressure = None if area is None else compression / area
                native = site['mechanical_input'] == 'native_indentation'
                indentation = native_indentations[key] * 1e6 if native else (
                    None if pressure is None else float(pressure_to_indentation_um(pressure, stiffness_pa_per_m=site['stiffness_pa_per_m'])))
                delta = None if temperature is None else temperature - site['reference_temperature_C']
                stimuli = {'rapid': indentation or 0., 'slow': indentation or 0., 'thermal': delta or 0.}
                responses = {}
                for kind, model in self.models[key].items():
                    if key in blocks or kind == 'thermal' and temperature is None:
                        reset = deepcopy(self.initial[key][kind])
                        reset['time_s'] = self.time_s + dt
                        model.restore(reset)
                        responses[kind] = 0.
                    else:
                        responses[kind] = float(model.advance([stimuli[kind]], dt)[0])
                rate = min(1000., self.recruitment * sum(abs(v) for v in responses.values()))
                if rate:
                    region = site['sensory_region']
                    inputs[region] = min(1000., inputs.get(region, 0.) + rate)
                rows.append({'id': key, 'sample_time_s': self.time_s,
                    'response_time_s': self.time_s + dt, 'delay_s': self.delay_s,
                    'force_n': force, 'position_m': current['point_m'], 'outward_normal': current['normal'],
                    'reference_position_m': list(site['position_m']),
                    'reference_outward_normal': list(site['normal']),
                    'position_basis': 'current native material point; source receptor support remains reference prior' if current['moving'] else 'registered reference material point',
                    'normal_force_n': compression, 'contact_area_m2': area,
                    'pressure_pa': pressure, 'indentation_um': indentation,
                    'skin_temperature_C': temperature, 'temperature_change_C': delta,
                    'mechanical_status': 'native_modeled_indentation' if native else (
                        'unknown_area_no_transduction' if area is None else 'explicit_foundation_prior'),
                    'material_identity': deepcopy(site.get('material_identity')),
                    'indentation_basis': site.get('indentation_basis', 'explicit linear foundation prior'),
                    'area_basis': site.get('area_basis', 'caller supplied contact area'),
                    'thermal_status': 'unavailable' if temperature is None else 'source_LTI_transient_with_reference_prior',
                    'rapid_response': responses['rapid'], 'slow_response': responses['slow'],
                    'thermal_response': responses['thermal'], 'sensory_input_hz': rate,
                    'sensory_region': site['sensory_region'], 'blocked': key in blocks,
                    'support_basis': site['support_basis']})
            self.time_s += dt
            return {'schema': 'ihm.cutaneous-feedback.v1', 'time_s': self.time_s,
                'sites': rows, 'sensory_inputs_hz': inputs, 'sensory_blocks': sorted(blocks),
                'model_sha256': self.identity, 'biological_validation': False,
                'brain_advanced': False, 'timing': 'physical sample held during interval; receptor endpoint feeds NEXT shared-brain exchange',
                'recruitment_basis': self.audit['recruitment_basis']}
        except Exception:
            self.restore(saved)
            raise

    def checkpoint(self):
        return {'identity': self.identity, 'time_s': self.time_s,
                'sites': {key: {kind: model.checkpoint() for kind, model in models.items()}
                          for key, models in self.models.items()}}

    def restore(self, checkpoint):
        if not isinstance(checkpoint, dict) or checkpoint.get('identity') != self.identity:
            raise ValueError('Cutaneous checkpoint identity mismatch')
        time = _number(checkpoint.get('time_s'), 'checkpoint time', 0., 1e9)
        values = checkpoint.get('sites')
        if not isinstance(values, dict) or set(values) != set(self.models):
            raise ValueError('Invalid checkpoint sites')
        for key, models in self.models.items():
            if not isinstance(values[key], dict) or set(values[key]) != set(models):
                raise ValueError('Invalid checkpoint receptor kinds')
            for kind in models:
                if not isinstance(values[key][kind], dict) or values[key][kind].get('time_s') != time:
                    raise ValueError('Receptor checkpoint clock mismatch')
        before = self.checkpoint()
        try:
            for key, models in self.models.items():
                for kind, model in models.items():
                    model.restore(values[key][kind])
            self.time_s = time
        except Exception:
            for key, models in self.models.items():
                for kind, model in models.items():
                    model.restore(before['sites'][key][kind])
            raise
