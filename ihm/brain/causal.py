"""Bounded causal realization of selected pinned IBM LTI transfers.

This is an IHM state-space adapter, not IBM's periodic runtime or its full graph.
Inputs are held constant on each advance interval; outputs are signed endpoints.
"""
import copy
import hashlib
import inspect
import json
from types import SimpleNamespace

import numpy as np
from scipy import signal


class CausalIBM:
    def __init__(self, backend, *, sites_m, kind='rapid', parameters=None,
                 delay_s=0., semantics='direct'):
        sites = np.asarray(sites_m, dtype=float)
        if sites.ndim != 2 or sites.shape[1] != 3 or not len(sites) or not np.isfinite(sites).all():
            raise ValueError('explicit nonempty finite sites_m (sites,3) required')
        if not np.isfinite(delay_s) or delay_s < 0 or semantics not in ('direct', 'drift'):
            raise ValueError('invalid delay or semantics')
        fn, medians, implementation = backend._transfer_spec(kind)
        theta = {k: p.default for k, p in inspect.signature(fn).parameters.items()
                 if k != 'basis' and p.default is not inspect.Parameter.empty}
        theta.update(medians)
        supplied = dict(parameters or {})
        if set(supplied) - set(theta):
            raise ValueError('unknown donor parameter')
        theta.update(supplied)
        theta = {k: float(v) for k, v in theta.items()}
        if not all(np.isfinite(v) for v in theta.values()):
            raise ValueError('nonfinite donor parameter')
        times = [v for k, v in theta.items() if k.endswith('_s')]
        if not times or min(times) < 1e-5 or max(times) > 1e3 or max(times)/min(times) > 1e6:
            raise ValueError('time constants outside validated numerical domain [1e-5,1e3] s')
        if not 0 <= theta.get('static_fraction', .4) <= 1:
            raise ValueError('static_fraction must be in [0,1]')
        if abs(theta.get('gain', 1.)) > 1e6:
            raise ValueError('gain outside numerical domain')
        gain = theta.get('gain', 1.)
        if kind in ('rapid', 'slow'):
            ta, tm = theta['tau_adapt_s'], theta['tau_membrane_s']
            num = np.array([gain*ta, gain*(theta['static_fraction'] if kind == 'slow' else 0.)])
            factors = [ta, tm]
        elif kind in ('twitch', 'activation'):
            num = np.array([gain])
            factors = [theta['contraction_time_s']]*2
            if kind == 'activation':
                factors += [theta['tau_rise_s'], theta['tau_fall_s']]
        else:
            raise ValueError('unsupported donor transfer')
        if semantics == 'drift':
            factors.append(1.)
        den = np.array([1.])
        for tau in factors:
            den = np.convolve(den, [tau, 1.])
        # A zero gain still has well-defined donor poles; avoid scipy's leading-zero warning.
        A, B, C, D = signal.tf2ss(num if gain else [1.], den)
        if gain == 0:
            C *= 0
            D *= 0
        omega = np.r_[0., np.geomspace(1e-4, 1e5, 801)]
        basis = SimpleNamespace(omega=omega, k=len(omega))
        donor = np.asarray(backend.transfer(kind, basis, theta))
        if semantics == 'drift':
            donor = donor/(1+1j*omega)
        realized = np.array([(C @ np.linalg.solve(1j*w*np.eye(len(A))-A, B)+D).item()
                             for w in omega])
        error = float(np.max(np.abs(realized-donor))/max(1e-30, float(np.max(np.abs(donor)))))
        if not np.isfinite(realized).all() or error > 1e-9:
            raise ValueError('state-space donor parity failed; outside validated numerical domain')
        self._A, self._B, self._C, self._D = A, B, C, D
        self._delay = float(delay_s)
        self._x = np.zeros((len(A), len(sites)))
        self._held = np.zeros(len(sites))
        self._queue = []
        self.time_s = 0.
        self.audit = {
            'backend': 'IHM scipy exact-ZOH realization of pinned IBM transfer',
            'package_sha256': backend.identity['package_sha256'],
            'implementation': implementation, 'parameters': theta,
            'semantics': semantics, 'transfer': 'H(s)' if semantics == 'direct' else 'H(s)/(s+1)',
            'delay_s': self._delay, 'delay_origin': 'explicit caller conduction delay',
            'input_units': 'um indentation' if kind in ('rapid', 'slow') else 'unit drive',
            'output_units': 'donor response; no calibrated receptor mV or muscle N',
            'sites_m': sites.tolist(), 'state_order': len(A),
            'poles_real_per_s': [-1/t for t in factors],
            'donor_relative_error': error, 'parity_omega_rad_s': [0., 1e-4, 1e5],
            'parity_samples': len(omega), 'approximation': 'exact rational coefficients; floating point only',
            'full_ibm_graph_streaming': False,
        }
        self._identity = hashlib.sha256(json.dumps(self.audit, sort_keys=True).encode()).hexdigest()
        self._zoh_cache = {}

    def _integrate(self, duration):
        if duration <= 0:
            return
        if duration not in self._zoh_cache:
            ad, bd, _, _, _ = signal.cont2discrete(
                (self._A, self._B, self._C, self._D), duration, method='zoh')
            if len(self._zoh_cache) >= 64:
                self._zoh_cache.clear()
            self._zoh_cache[duration] = ad, bd
        ad, bd = self._zoh_cache[duration]
        self._x = ad @ self._x + bd @ self._held[None, :]

    def advance(self, values, dt_s):
        """Hold values over this interval; return response at its right endpoint."""
        values = np.asarray(values, dtype=float)
        if values.shape != self._held.shape or not np.isfinite(values).all():
            raise ValueError('one finite input value per support site required')
        if not np.isfinite(dt_s) or dt_s <= 0 or not np.isfinite(self.time_s + dt_s):
            raise ValueError('positive finite dt_s required')
        end = self.time_s + float(dt_s)
        self._queue.append((self.time_s + self._delay, values.copy()))
        while self._queue and self._queue[0][0] <= end:
            stamp, arriving = self._queue.pop(0)
            self._integrate(stamp-self.time_s)
            self.time_s = stamp
            self._held = arriving
        self._integrate(end-self.time_s)
        self.time_s = end
        return (self._C @ self._x + self._D @ self._held[None, :]).reshape(-1).copy()

    def checkpoint(self):
        return {'identity': self._identity, 'time_s': self.time_s,
                'state': self._x.tolist(), 'held': self._held.tolist(),
                'queue': [(t, v.tolist()) for t, v in self._queue]}

    def restore(self, checkpoint):
        """Atomic restore; different source, support, parameters or delay is rejected."""
        data = copy.deepcopy(checkpoint)
        if data.get('identity') != self._identity:
            raise ValueError('checkpoint materialization identity differs')
        state, held = np.asarray(data['state'], float), np.asarray(data['held'], float)
        time = float(data['time_s'])
        queue = [(float(t), np.asarray(v, float)) for t, v in data['queue']]
        if (state.shape != self._x.shape or held.shape != self._held.shape
                or not np.isfinite(state).all() or not np.isfinite(held).all()
                or not np.isfinite(time) or time < 0
                or any(not np.isfinite(t) or t < time or t > time + self._delay or v.shape != held.shape
                       or not np.isfinite(v).all() for t, v in queue)
                or any(a[0] > b[0] for a, b in zip(queue, queue[1:]))):
            raise ValueError('invalid causal checkpoint state')
        self._x, self._held, self.time_s, self._queue = state, held, time, queue
