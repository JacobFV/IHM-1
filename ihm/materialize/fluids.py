"""conservative compartment fluid exchange on heterogeneous interface graphs.

Each edge uses Q = G(P_a - P_b + head), P = P0 + (V-V0)/compliance.
Fixed head can represent a pump or a linearized oncotic difference. This is
not a full Starling/glycocalyx model. Valves are fixed-regime resistances;
negative predicted flow is reported, never silently clipped.
"""
from dataclasses import dataclass, asdict
import math
import re
import numpy as np
from ihm.fields import Support, Component
from ihm.anatomy import Partition
from ihm.topologies import Topology
from ihm.processes import Process
from ihm.materialize import Request, materialize


@dataclass(frozen=True)
class FluidCompartment:
    id: str
    system: str
    volume_mL: float
    pressure_mmHg: float
    compliance_mL_per_mmHg: float
    std_mL: float


@dataclass(frozen=True)
class FluidEdge:
    a: str
    b: str
    interface: str
    conductance_mL_per_s_mmHg: float
    head_mmHg: float = 0


@dataclass(frozen=True)
class FluidNetwork:
    compartments: tuple[FluidCompartment, ...]
    edges: tuple[FluidEdge, ...]

    def __post_init__(self):
        ids = [c.id for c in self.compartments]
        if not ids or len(ids) > 512 or len(set(ids)) != len(ids):
            raise ValueError('unique compartments, at most 512, required')
        for c in self.compartments:
            if not re.fullmatch('[a-z][a-z0-9_]*', c.id) or not re.fullmatch('[a-z][a-z0-9_]*', c.system):
                raise ValueError('invalid fluid id or system')
            if any(not math.isfinite(v) or v <= 0 for v in (c.volume_mL, c.compliance_mL_per_mmHg, c.std_mL)) or not math.isfinite(c.pressure_mmHg):
                raise ValueError('invalid fluid compartment parameters')
        seen = set()
        for e in self.edges:
            if e.a not in ids or e.b not in ids or e.a == e.b or not re.fullmatch('[a-z][a-z0-9_]*', e.interface):
                raise ValueError('invalid fluid edge')
            if not math.isfinite(e.conductance_mL_per_s_mmHg) or e.conductance_mL_per_s_mmHg <= 0 or not math.isfinite(e.head_mmHg):
                raise ValueError('invalid fluid conductance/head')
            key = (tuple(sorted((e.a, e.b))), e.interface)
            if key in seen:
                raise ValueError('duplicate fluid interface')
            seen.add(key)

    @classmethod
    def example(cls):
        return cls((
            FluidCompartment('capillary', 'vascular', 50, 20, 5, 5),
            FluidCompartment('interstitium', 'interstitial', 200, 0, 50, 20),
            FluidCompartment('initial_lymph', 'lymphatic', 5, 0, 2, 1),
            FluidCompartment('collector', 'lymphatic', 10, 1, 2, 2),
            FluidCompartment('node', 'lymphoid', 5, 1, 2, 1),
            FluidCompartment('vein', 'venous', 100, 5, 20, 10)), (
            FluidEdge('capillary', 'interstitium', 'capillary_filtration', .001, -18),
            FluidEdge('interstitium', 'initial_lymph', 'lymphatic_uptake', .005, 1),
            FluidEdge('initial_lymph', 'collector', 'lymphatic_collecting', .003, 2),
            FluidEdge('collector', 'node', 'lymph_node_transit', .003, 1),
            FluidEdge('node', 'vein', 'lymphovenous_return', .002, 5),
            FluidEdge('capillary', 'vein', 'vascular_return', .002)))


def register_fluids(r, network):
    by_id = {c.id: c for c in network.compartments}
    for system in sorted({c.system for c in network.compartments}):
        regions = tuple('fluid_'+c.id for c in network.compartments if c.system == system)
        r.add(Support('fluid_'+system, 'named_fluid_compartments', regions))
    regions = tuple('fluid_'+c.id for c in network.compartments)
    r.add(Topology('fluid_local', tuple((x, x) for x in regions)))
    r.add(Partition('fluid_systems', {'fluid_'+c.id: {c.system: 1.} for c in network.compartments}))
    for c in network.compartments:
        ident = f'fluid.{c.id}.volume'
        r.add(Component(ident, 'fluid_'+c.system, 'fluid_'+c.id, 'mL', c.volume_mL, c.std_mL))
        r.add(Process('fluid.storage.'+c.id, (ident,), ident, 'fluid_local', (0.,), 0., 0.,
                      'conserved_storage;illustrative_compliance'))
    for interface in sorted({e.interface for e in network.edges}):
        pairs = set()
        for e in network.edges:
            if e.interface == interface:
                a, b = 'fluid_'+e.a, 'fluid_'+e.b
                pairs.update(((a, b), (b, a), (a, a), (b, b)))
        r.add(Topology('fluid_'+interface, tuple(sorted(pairs))))
    for j, e in enumerate(network.edges):
        a, b = by_id[e.a], by_id[e.b]; g = e.conductance_mL_per_s_mmHg
        offset = g*(a.pressure_mmHg - a.volume_mL/a.compliance_mL_per_mmHg
                    -b.pressure_mmHg + b.volume_mL/b.compliance_mL_per_mmHg + e.head_mmHg)
        weights = (g/a.compliance_mL_per_mmHg, -g/b.compliance_mL_per_mmHg)
        inputs = (f'fluid.{a.id}.volume', f'fluid.{b.id}.volume')
        for dst, sign in ((a, -1), (b, 1)):
            r.add(Process(f'fluid.edge_{j}.to_{dst.id}', inputs, f'fluid.{dst.id}.volume',
                          'fluid_'+e.interface, tuple(sign*w for w in weights), sign*offset, 0.,
                          'conservative_linear_hydraulics;illustrative_parameters;fixed_valve_regime'))


def fluid_model(network, subject='synthetic-fluid'):
    from ihm.body import body
    targets = tuple(f'fluid.{c.id}.volume' for c in network.compartments)
    m = materialize(body(fluids=network), Request(targets, subject))
    m.provenance.update({'representation': 'conservative_fluid_network', 'network': asdict(network),
                         'limitations': ['linear compliance and fixed valve regime',
                           'no active valve gating, solute transport or protein mass balance',
                           'negative volume/flow marks extrapolation outside regime',
                           'illustrative parameter values; no human predictive validation']})
    return m


def flows(model, network):
    by_id = {c.id: c for c in network.compartments}
    h = np.zeros((len(network.edges), len(model.components))); offset = []
    for j, e in enumerate(network.edges):
        a, b = by_id[e.a], by_id[e.b]; g = e.conductance_mL_per_s_mmHg
        h[j, model.components.index(f'fluid.{a.id}.volume')] = g/a.compliance_mL_per_mmHg
        h[j, model.components.index(f'fluid.{b.id}.volume')] = -g/b.compliance_mL_per_mmHg
        offset.append(g*(a.pressure_mmHg-a.volume_mL/a.compliance_mL_per_mmHg
                        -b.pressure_mmHg+b.volume_mL/b.compliance_mL_per_mmHg+e.head_mmHg))
    value = h @ model.mean + offset
    return {'mL_per_s': value, 'std_mL_per_s': np.sqrt(np.maximum(0, np.diag(h @ model.cov @ h.T))),
            'reverse_flow_edges': np.flatnonzero(value < 0).tolist()}
