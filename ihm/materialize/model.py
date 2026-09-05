"""portable explicit predictor, retaining evidence and process provenance."""
import copy
import json
import math
from pathlib import Path
from dataclasses import asdict
import numpy as np
from ihm.runtime.gaussian import transition, covariance


class Model:
    def __init__(self, components, targets, subject, mean, cov, a, b, q, provenance, time=0., seen=()):
        self.components = tuple(components); self.targets = tuple(targets); self.subject = subject
        self.mean = np.asarray(mean, float).copy(); n = len(self.components)
        self.cov = covariance(cov, n).copy(); self.a = np.asarray(a, float).copy()
        self.b = np.asarray(b, float).copy(); self.q = covariance(q, n).copy()
        if self.mean.shape != (n,) or self.a.shape != (n, n) or self.b.shape != (n,):
            raise ValueError('invalid model dimensions')
        if not all(np.isfinite(v).all() for v in (self.mean, self.a, self.b)):
            raise ValueError('nonfinite model parameters')
        if not math.isfinite(time) or time < 0 or not subject or not n or len(set(components)) != n or not set(targets) <= set(components):
            raise ValueError('invalid model metadata')
        self.time = float(time); self.provenance = copy.deepcopy(provenance)
        self.seen = set(map(tuple, seen))

    def advance(self, time, clamps=None):
        if not math.isfinite(time) or time < self.time:
            raise ValueError('time must be finite and nondecreasing')
        clamps = clamps or {}
        if any(c not in self.components or not math.isfinite(v) for c, v in clamps.items()):
            raise ValueError('invalid clamp')
        a = self.a.copy(); b = self.b.copy(); q = self.q.copy()
        mean = self.mean.copy(); cov = self.cov.copy()
        for c, value in clamps.items():
            i = self.components.index(c)
            mean[i] = value; cov[i, :] = 0; cov[:, i] = 0
            a[i, :] = 0; b[i] = 0; q[i, :] = 0; q[:, i] = 0
        f, offset, noise = transition(a, b, q, time-self.time)
        mean = f @ mean + offset; cov = f @ cov @ f.T + noise
        if not np.isfinite(mean).all():
            raise ValueError('nonfinite prediction; revise dynamics or horizon')
        covariance(cov, len(mean))
        self.mean = mean; self.cov = (cov+cov.T)/2
        if clamps:
            self.provenance.setdefault('interventions', []).append({'start': self.time, 'end': time, 'clamps': clamps.copy()})
        self.time = float(time)
        return self

    def assimilate(self, evidence, error_covariance=None):
        evidence = list(evidence)
        if not evidence:
            return self
        times = {e.time for e in evidence}; keys = [e.key for e in evidence]
        if len(times) != 1 or min(times) < self.time:
            raise ValueError('a conditioning batch must share one nonpast time')
        if any(e.subject != self.subject or e.component not in self.components for e in evidence):
            raise ValueError('subject mismatch or component outside materialization')
        if len(set(keys)) != len(keys) or self.seen.intersection(keys):
            raise ValueError('duplicate evidence')
        sources = copy.deepcopy(self.provenance.get('sources', {}))
        for e in evidence:
            manifest = {'reference': e.reference, 'species': e.species, 'kind': e.kind,
                        'bindings': list(e.bindings)}
            if e.source in sources and sources[e.source]['card'] != manifest:
                raise ValueError('same source id refers to conflicting source metadata')
            entry = sources.setdefault(e.source, {'card': manifest, 'data_sha256': []})
            if e.data_sha256 and e.data_sha256 not in entry['data_sha256']:
                entry['data_sha256'].append(e.data_sha256)
        n, k = len(self.components), len(evidence)
        r = np.diag([e.variance for e in evidence]) if error_covariance is None else covariance(error_covariance, k, positive=True)
        if not np.allclose(np.diag(r), [e.variance for e in evidence], rtol=1e-8, atol=1e-15):
            raise ValueError('error covariance diagonal must match evidence variance')
        work = copy.deepcopy(self); work.advance(evidence[0].time)
        h = np.zeros((k, n))
        for j, e in enumerate(evidence):
            h[j, self.components.index(e.component)] = 1
        gain = np.linalg.solve(h @ work.cov @ h.T + r, h @ work.cov).T
        mean = work.mean + gain @ (np.array([e.value for e in evidence]) - h @ work.mean)
        residual = np.eye(n) - gain @ h
        cov = residual @ work.cov @ residual.T + gain @ r @ gain.T
        covariance(cov, n)
        self.mean = mean; self.cov = (cov+cov.T)/2; self.time = work.time
        self.provenance['sources'] = sources
        self.seen.update(keys)
        self.provenance['evidence'].extend(asdict(e) for e in evidence)
        self.provenance.setdefault('conditioning', []).append({'time': self.time, 'error_covariance': r.tolist(),
            'assumption': 'batch errors independent of previous batches; caller responsible for cross-batch/source dependence'})
        return self

    def forecast(self, times, clamps=None):
        times = list(times)
        if any(not math.isfinite(t) or t < self.time for t in times) or any(b < a for a, b in zip(times, times[1:])):
            raise ValueError('forecast times must be sorted, finite and nonpast')
        work = copy.deepcopy(self); idx = [self.components.index(c) for c in self.targets]
        means, stds = [], []
        for time in times:
            work.advance(time, clamps)
            means.append(work.mean[idx].copy()); stds.append(np.sqrt(np.maximum(0, np.diag(work.cov)[idx])))
        return {'time': np.asarray(times), 'components': self.targets,
                'mean': np.asarray(means).reshape(len(times), len(idx)),
                'std': np.asarray(stds).reshape(len(times), len(idx)),
                'validated_biology': False}

    def save(self, path):
        metadata = {'version': 1, 'components': self.components, 'targets': self.targets,
                    'subject': self.subject, 'time': self.time, 'seen': sorted(self.seen),
                    'provenance': self.provenance}
        path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('wb') as stream:
            np.savez_compressed(stream, mean=self.mean, cov=self.cov, a=self.a, b=self.b, q=self.q,
                                metadata=json.dumps(metadata, allow_nan=False))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as d:
            meta = json.loads(str(d['metadata']))
            if meta.pop('version') != 1:
                raise ValueError('unsupported model format')
            return cls(**meta, **{k: d[k] for k in ('mean', 'cov', 'a', 'b', 'q')})
