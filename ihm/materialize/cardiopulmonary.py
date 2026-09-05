"""nonlinear explicit view of registered cardiopulmonary cycle processes."""
from dataclasses import asdict
from pathlib import Path
import json
import math
import numpy as np
from scipy.integrate import solve_ivp
from ihm.processes.cardiopulmonary import Parameters, COMPONENTS, initial_state, rhs, signals


class CardiopulmonaryModel:
    def __init__(self, parameters=None, subject='supine-example', state=None, time=0.):
        self.parameters = parameters or Parameters(); self.subject = subject
        self.components = COMPONENTS; self.time = float(time)
        self.state = initial_state(self.parameters) if state is None else np.asarray(state, float).copy()
        if not subject or not math.isfinite(time) or time < 0 or self.state.shape != (11,) or not np.isfinite(self.state).all() or np.any(self.state[:9] <= 0):
            raise ValueError('invalid cardiopulmonary state or subject')
        from ihm.body import body
        r = body(cardiopulmonary=self.parameters)
        self.provenance = {'representation': 'nonlinear_cardiopulmonary_lumped', 'validated_biology': False,
            'posture': 'supine', 'parameters': asdict(self.parameters),
            'components': [asdict(r.components[c]) for c in COMPONENTS],
            'processes': [asdict(v) for v in r.processes.values() if v.form == 'cardiopulmonary_lumped'],
            'topology': asdict(r.topologies['cardiopulmonary_cycle']),
            'limitations': ['deterministic trajectory, no uncertainty calibration or assimilation in this view',
                'fixed-rate phase pacemakers; no ECG, chemoreflex or baroreflex model',
                'one-way coupling: breathing changes circulation; no circulatory feedback into breathing',
                'no gas exchange, oxygen consumption or carbon dioxide chemistry in this view',
                'ideal valves and lumped chambers/vessels; no posture change or movement',
                'illustrative parameters and startup transient; not a fitted resting human']}

    def forecast(self, times):
        times = np.asarray(times, float)
        if times.ndim != 1 or not len(times) or not np.isfinite(times).all() or np.any(times < self.time) or np.any(np.diff(times) < 0):
            raise ValueError('nonempty sorted finite nonpast forecast times required')
        if times[-1]-self.time > 600 or len(times) > 100000:
            raise ValueError('forecast exceeds 600-second / 100000-output-point budget')
        unique, inverse = np.unique(times, return_inverse=True)
        if unique[-1] == self.time:
            states = np.repeat(self.state[None, :], len(times), axis=0)
        else:
            solved = solve_ivp(lambda t, x: rhs(t, x, self.parameters), (self.time, unique[-1]),
                               self.state, t_eval=unique, rtol=1e-7, atol=1e-9,
                               max_step=min(.01, 60/self.parameters.heart_rate_bpm/100))
            if not solved.success:
                raise ValueError('cardiopulmonary integration failed: '+solved.message)
            states = solved.y.T[inverse]
        if not np.isfinite(states).all() or np.any(states[:, :9] <= 0):
            raise ValueError('simulation left the positive-volume operating regime')
        all_signals = [signals(x, self.parameters) for x in states]
        return {'time': times, 'components': COMPONENTS, 'state': states,
                'signals': {key: np.asarray([s[key] for s in all_signals]) for key in all_signals[0]},
                'posture': 'supine', 'validated_biology': False,
                'uncertainty': 'not quantified; deterministic mechanistic trajectory'}

    def advance(self, time):
        result = self.forecast([time]); self.state = result['state'][-1].copy(); self.time = float(time)
        return self

    def save(self, path):
        path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
        meta = {'version': 1, 'model_type': 'cardiopulmonary', 'subject': self.subject, 'time': self.time,
                'parameters': asdict(self.parameters), 'provenance': self.provenance}
        with path.open('wb') as stream:
            np.savez_compressed(stream, state=self.state, metadata=json.dumps(meta, allow_nan=False))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as data:
            meta = json.loads(str(data['metadata']))
            if meta.get('version') != 1 or meta.get('model_type') != 'cardiopulmonary':
                raise ValueError('unsupported cardiopulmonary artifact')
            m = cls(Parameters(**meta['parameters']), meta['subject'], data['state'], meta['time'])
            m.provenance = meta['provenance']
            return m


def cardiopulmonary_model(subject='supine-example', parameters=None):
    return CardiopulmonaryModel(parameters, subject)
