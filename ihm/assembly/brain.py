"""Canonical regional brain executing the preserved IBM adaptive population law.

Spatial identity comes from BP3D tissue, transferred DK labels from fsaverage.
Physiology-to-neural drive and autonomic readout are explicit generic priors.
"""
from __future__ import annotations

import ast
import hashlib
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
INPUT_BASELINES = {'mean_arterial_pressure_mmHg': 90., 'oxygen_saturation': .98,
                   'core_temperature_C': 37.}


def verify_sources(data, root=None):
    """Verify all preserved source bytes before source-derived code executes."""
    base = Path(root or ROOT)
    verified = []
    for source in data['sources']:
        path = (base / source['preserved_path']).resolve()
        if not path.is_relative_to(base.resolve()):
            raise ValueError('Source snapshot must be inside the configured root')
        if hashlib.sha256(path.read_bytes()).hexdigest() != source['sha256']:
            raise ValueError(f'Source hash mismatch: {source["preserved_path"]}')
        verified.append(source['preserved_path'])
    return verified


def _load_rate_law(data, root):
    verify_sources(data, root)
    source = next(s for s in data['sources'] if s.get('role') == 'executed_neural_law')
    path = Path(root) / source['preserved_path']
    original = ast.parse(path.read_text())
    names = {'_sigmoid', 'wilson_cowan_excitatory', 'shunting_inhibition_rate'}
    selected = [n for n in original.body if isinstance(n, ast.FunctionDef) and n.name in names]
    if {n.name for n in selected} != names:
        raise ValueError('Preserved IBM source does not contain required neural laws')
    # Exact function ASTs avoid importing IBM's unrelated registry and NN dependencies.
    module = ast.Module(body=selected, type_ignores=[])
    namespace = {'np': np}
    exec(compile(module, str(path), 'exec'), namespace)
    return namespace['wilson_cowan_excitatory'], namespace['shunting_inhibition_rate']


class BodyBrain:
    """Reusable regional ODE with independent state and explicit input/output ports."""

    @classmethod
    def from_dict(cls, data, *, root=None, max_step_s=.001):
        return cls(data, root=root, max_step_s=max_step_s)

    def __init__(self, data, *, root=None, max_step_s=.001):
        if not math.isfinite(max_step_s) or not 0 < max_step_s <= .002:
            raise ValueError('max_step_s must be in (0, 0.002] for the 5 ms source rate relaxation')
        self.data = data
        self.max_step_s = max_step_s
        self.source_rate_law, self.source_inhibition_law = _load_rate_law(data, Path(root or ROOT))
        self.theta = dict(data['parameters']['ibm_wilson_cowan'])
        self.coupling = dict(data['parameters']['body_transfer_priors'])
        n = len(data['nodes'])
        if n == 0:
            raise ValueError('Brain has no regional tissue supports')
        self.ids = [node['id'] for node in data['nodes']]
        indices = {node_id: i for i, node_id in enumerate(self.ids)}
        if len(indices) != n:
            raise ValueError('Duplicate brain region identity')
        self.weights = np.zeros((n, n))
        for edge in data['edges']:
            i, j = indices[edge['source']], indices[edge['target']]
            w = float(edge['weight'])
            if not math.isfinite(w) or w < 0:
                raise ValueError('Neural edge weights must be finite nonnegative priors')
            self.weights[j, i] += w
        self.weights /= np.maximum(self.weights.sum(axis=1, keepdims=True), 1e-12)
        self.state = np.zeros((n, 3))
        self.state[:, 0] = self.theta['v_rest_mv']
        self.state[:, 1] = self.theta['r_max_hz'] / (1 + np.exp(
            -(self.theta['v_rest_mv'] - self.theta['v_half_mv']) / self.theta['slope_mv']))
        self.time_s = 0.

    def _rhs(self, state, availability, temperature_factor):
        v, r, a = state.T
        gain = self.coupling['recurrent_drive_nS_per_Hz']
        drive = availability * (self.coupling['tonic_drive_nS'] + gain * (self.weights @ r))
        x = {'neural.exc.potential': v, 'neural.exc.activity': r,
             'neural.exc.adaptation': a, 'neural.exc.ampa': drive * .7,
             'neural.exc.nmda': drive * .3,
             'neural.inh.gaba_a': self.coupling['inhibitory_conductance_per_Hz'] * r}
        rates = self.source_rate_law(x, self.theta)
        inhibition = self.source_inhibition_law(x, self.theta)
        return temperature_factor * np.column_stack((
            rates['neural.exc.potential'] + inhibition['neural.exc.potential'],
            rates['neural.exc.activity'], rates['neural.exc.adaptation']))

    def step(self, dt_s, physiology=None):
        """Advance seconds; MAP mmHg, saturation fraction, temperature Celsius.

        Missing inputs use labeled priors. Returns commands without applying them
        to an external body engine. RK4 substeps resolve IBM's 5 ms rate filter.
        """
        dt_s = float(dt_s)
        if not math.isfinite(dt_s) or not 0 <= dt_s <= 60.:
            raise ValueError('dt_s must be finite in [0, 60] seconds')
        supplied = physiology or {}
        values = {k: float(supplied.get(k, default)) for k, default in INPUT_BASELINES.items()}
        if not all(math.isfinite(v) for v in values.values()):
            raise ValueError('Brain physiology inputs must be finite')
        pressure, oxygen, temperature = (values[k] for k in INPUT_BASELINES)
        if not (0 <= pressure <= 300 and 0 <= oxygen <= 1 and 20 <= temperature <= 45):
            raise ValueError('Input outside supported units/range: MAP [0,300] mmHg, O2 [0,1], temperature [20,45] C')
        availability = float(np.clip(pressure / self.coupling['map_reference_mmHg'], 0., 1.)
                             * np.clip(oxygen / .98, 0., 1.))
        temperature_factor = self.coupling['temperature_Q10'] ** ((temperature - 37.) / 10.)
        steps = max(1, math.ceil(dt_s / self.max_step_s))
        h = dt_s / steps
        y = self.state.copy()
        for _ in range(steps):
            k1 = self._rhs(y, availability, temperature_factor)
            k2 = self._rhs(y + h * k1 / 2, availability, temperature_factor)
            k3 = self._rhs(y + h * k2 / 2, availability, temperature_factor)
            k4 = self._rhs(y + h * k3, availability, temperature_factor)
            y += h * (k1 + 2*k2 + 2*k3 + k4) / 6
        if not np.isfinite(y).all():
            raise FloatingPointError('Non-finite brain integration; original state retained')
        self.state = y
        self.time_s += dt_s
        medulla = [i for i,n in enumerate(self.data['nodes']) if n.get('autonomic_role') == 'medulla']
        medulla_rate = float(np.mean(y[medulla, 1])) if medulla else float(np.mean(y[:, 1]))
        stress = 1. - availability
        sympathetic = float(np.clip(.2 + .7 * stress + .02 * max(0., temperature - 37.)
                                     + .05 * (1. - medulla_rate / 20.), 0., 1.))
        return {'time_s': self.time_s, 'model_id': self.data['id'],
                'regional_state': {'node_ids': self.ids, 'potential_mV': y[:,0].tolist(),
                                   'activity_hz': y[:,1].tolist(), 'adaptation_mV': y[:,2].tolist()},
                'autonomic_commands': {'sympathetic_fraction': sympathetic,
                    'parasympathetic_fraction': float(np.clip(.6 - .5 * stress, 0., 1.)),
                    'applied_to_body': False, 'calibration': 'uncalibrated generic transfer prior',
                    'ports': self.data['ports']['efferent']},
                'physiology_inputs': values,
                'input_provenance': {k: 'caller_supplied' if k in supplied else 'assumed_baseline' for k in values},
                'solver': {'method': 'RK4', 'substeps': steps, 'step_s': h},
                'oxygen_perfusion_availability': availability,
                'coupling_scope': 'body-driven neural ODE; autonomic commands available, external feedback not applied'}
