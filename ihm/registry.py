"""one checked namespace for the four primitives and their supports."""
import math
import re
from types import MappingProxyType
from ihm.fields import Component, Support
from ihm.anatomy import Partition
from ihm.topologies import Topology
from ihm.processes import Process


class Registry:
    def __init__(self):
        self.components = {}; self.supports = {}; self.anatomy = {}
        self.topologies = {}; self.processes = {}; self._ids = set()
        self._sealed = False

    def add(self, item):
        if self._sealed:
            raise ValueError('registry is sealed')
        if not re.fullmatch(r'[a-z][a-z0-9_]*(\.[a-z][a-z0-9_]*)*', item.id):
            raise ValueError(f'invalid id: {item.id}')
        if item.id in self._ids:
            raise ValueError(f'duplicate id: {item.id}')
        table = {Component: self.components, Support: self.supports,
                 Partition: self.anatomy, Topology: self.topologies, Process: self.processes}[type(item)]
        table[item.id] = item; self._ids.add(item.id)
        return item

    def seal(self):
        regions = {r for s in self.supports.values() for r in s.regions}
        for c in self.components.values():
            if c.support not in self.supports or c.region not in self.supports[c.support].regions:
                raise ValueError(f'unknown support/region: {c.id}')
            if not c.unit or not math.isfinite(c.mean) or not math.isfinite(c.std) or c.std <= 0:
                raise ValueError(f'invalid prior: {c.id}')
        for a in self.anatomy.values():
            for region, weights in a.memberships.items():
                if region not in regions or any(not math.isfinite(w) or w < 0 for w in weights.values()) or sum(weights.values()) > 1 + 1e-12:
                    raise ValueError(f'invalid membership: {a.id}/{region}')
        for t in self.topologies.values():
            if any(a not in regions or b not in regions for a, b in t.edges):
                raise ValueError(f'unknown topology region: {t.id}')
        outputs = set()
        for p in self.processes.values():
            if p.topology not in self.topologies or len(p.inputs) != len(p.weights):
                raise ValueError(f'invalid process: {p.id}')
            if p.form not in ('affine', 'cardiopulmonary_lumped') or len(set(p.inputs)) != len(p.inputs):
                raise ValueError(f'unsupported/duplicate process inputs: {p.id}')
            if any(c not in self.components for c in (*p.inputs, p.output)):
                raise ValueError(f'dangling component: {p.id}')
            if not all(math.isfinite(v) for v in (*p.weights, p.bias, p.noise)) or p.noise < 0:
                raise ValueError(f'invalid parameters: {p.id}')
            target = self.components[p.output].region
            edges = self.topologies[p.topology].edges
            if any((self.components[c].region, target) not in edges for c in p.inputs):
                raise ValueError(f'process outside topology: {p.id}')
            outputs.add(p.output)
        if set(self.components) - outputs:
            raise ValueError(f'components without dynamics: {set(self.components) - outputs}')
        for name in ('components', 'supports', 'anatomy', 'topologies', 'processes'):
            setattr(self, name, MappingProxyType(getattr(self, name)))
        self._sealed = True
        return self

    def describe(self):
        return [{'id': c.id, 'unit': c.unit, 'support': c.support, 'region': c.region,
                 'prior': c.provenance} for c in self.components.values()]
