"""heterogeneous integumentary electrical networks on explicit tissue supports.

The apical extracellular RC sheet and basal cell membrane network are distinct.
Basal extracellular potential is grounded. This is a prescribed-bath approximation,
not a bidomain electrodiffusion solution. Gap junctions connect only declared
cell contacts; the extracellular graph follows separate geometric contacts.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np
from ihm.fields import Component, Support
from ihm.anatomy import Partition
from ihm.topologies import Topology
from ihm.processes import Process
from ihm.materialize import Request, materialize


@dataclass(frozen=True)
class Cell:
    position_m: tuple[float, float, float]
    phenotype: str = 'keratinocyte'
    capacitance_f: float = 1e-11
    sodium_conductance_s: float = 1e-11
    potassium_conductance_s: float = 1e-10
    chloride_conductance_s: float = 2e-11
    sodium_inside_mM: float = 15
    potassium_inside_mM: float = 140
    chloride_inside_mM: float = 30
    pump_current_a: float = 1e-12


@dataclass(frozen=True)
class SurfaceSite:
    position_m: tuple[float, float, float]
    capacitance_f: float = 1e-6
    barrier_conductance_s: float = 1e-6
    pump_current_a: float = -3e-8
    wounded: bool = False


@dataclass(frozen=True)
class SkinPatch:
    cells: tuple[Cell, ...]
    surface: tuple[SurfaceSite, ...]
    gap_edges: tuple[tuple[int, int, float], ...]
    surface_edges: tuple[tuple[int, int, float], ...]
    temperature_k: float = 310.15
    sodium_outside_mM: float = 140
    potassium_outside_mM: float = 4
    chloride_outside_mM: float = 110

    def __post_init__(self):
        if not self.cells or not self.surface or len(self.cells)*4+len(self.surface) > 512:
            raise ValueError('skin patch requires cells/surface and at most 512 electrical/ionic states')
        if any(not math.isfinite(v) or v <= 0 for v in (self.temperature_k, self.sodium_outside_mM, self.potassium_outside_mM, self.chloride_outside_mM)):
            raise ValueError('positive bath concentrations and temperature required')
        for c in self.cells:
            values = (c.capacitance_f, c.sodium_conductance_s, c.potassium_conductance_s,
                      c.chloride_conductance_s, c.sodium_inside_mM, c.potassium_inside_mM, c.chloride_inside_mM)
            if any(not math.isfinite(v) or v <= 0 for v in values) or not math.isfinite(c.pump_current_a):
                raise ValueError('invalid membrane parameters')
            if c.phenotype not in ('keratinocyte', 'fibroblast', 'melanocyte'):
                raise ValueError('declare a supported cell phenotype')
        for s in self.surface:
            if not all(math.isfinite(v) for v in (s.capacitance_f, s.barrier_conductance_s, s.pump_current_a)) or min(s.capacitance_f, s.barrier_conductance_s) <= 0:
                raise ValueError('invalid surface parameters')
        for elements in (self.cells, self.surface):
            for element in elements:
                if len(element.position_m) != 3 or not np.isfinite(element.position_m).all():
                    raise ValueError('positions must be finite 3-vectors in meters')
        for edges, elements in ((self.gap_edges, self.cells), (self.surface_edges, self.surface)):
            seen = set()
            for a, b, g in edges:
                if not isinstance(a, int) or not isinstance(b, int) or a == b or not (0 <= a < len(elements) and 0 <= b < len(elements)) or not math.isfinite(g) or g <= 0:
                    raise ValueError('invalid conductance edge')
                pair = tuple(sorted((a, b)))
                if pair in seen or np.linalg.norm(np.subtract(elements[a].position_m, elements[b].position_m)) == 0:
                    raise ValueError('duplicate edge or coincident connected sites')
                seen.add(pair)

    @classmethod
    def line(cls, n=9, spacing_m=2e-4, wound=True):
        """small synthetic cross-section; arbitrary 3D graphs can be supplied directly."""
        if not isinstance(n, int) or n < 3 or not math.isfinite(spacing_m) or spacing_m <= 0:
            raise ValueError('line requires n >= 3 and positive spacing')
        surface = tuple(SurfaceSite((i*spacing_m, 0., 0.),
                        barrier_conductance_s=1e-4 if wound and i == n//2 else 1e-6,
                        pump_current_a=0. if wound and i == n//2 else -3e-8,
                        wounded=bool(wound and i == n//2)) for i in range(n))
        cells = tuple(Cell((i*spacing_m, 0., -1e-4),
                           potassium_conductance_s=1e-10*(1+.1*(i % 3)))
                      for i in range(n) if not (wound and i == n//2))
        gaps = tuple((i, i+1, 1e-11) for i in range(len(cells)-1)
                     if cells[i+1].position_m[0]-cells[i].position_m[0] < 1.5*spacing_m)
        return cls(cells, surface, gaps, tuple((i, i+1, 1e-6) for i in range(n-1)))


def register_skin(r, patch):
    cell_regions = tuple(f'skin_cell_{i}' for i in range(len(patch.cells)))
    surface_regions = tuple(f'skin_surface_{i}' for i in range(len(patch.surface)))
    r.add(Support('skin_cell_membranes', 'skin_patch_m', cell_regions))
    r.add(Support('skin_apical_surface', 'skin_patch_m', surface_regions))
    r.add(Partition('skin_cell_types', {region: {c.phenotype: 1.} for region, c in zip(cell_regions, patch.cells)}))
    r.add(Partition('skin_wound_partition', {region: {'wound' if s.wounded else 'intact': 1.}
                                          for region, s in zip(surface_regions, patch.surface)}))
    r.add(Topology('skin_membrane_local', tuple((x, x) for x in cell_regions)))
    r.add(Topology('skin_barrier_local', tuple((x, x) for x in surface_regions)))
    def network_topology(ident, edges, regions):
        pairs = {(regions[i], regions[i]) for edge in edges for i in edge[:2]}
        for a, b, _ in edges:
            pairs.update(((regions[a], regions[b]), (regions[b], regions[a])))
        r.add(Topology(ident, tuple(sorted(pairs))))
    network_topology('skin_gap_junction', patch.gap_edges, cell_regions)
    network_topology('skin_extracellular', patch.surface_edges, surface_regions)
    thermal_voltage = 8.314462618*patch.temperature_k/96485.33212
    mechanism = 'physical_RC_form;illustrative_parameters;fixed_bath;linearized_Nernst'
    for i, (region, c) in enumerate(zip(cell_regions, patch.cells)):
        prefix = f'skin.cell_{i}'
        ions = (('sodium', 1, c.sodium_inside_mM, patch.sodium_outside_mM, c.sodium_conductance_s),
                ('potassium', 1, c.potassium_inside_mM, patch.potassium_outside_mM, c.potassium_conductance_s),
                ('chloride', -1, c.chloride_inside_mM, patch.chloride_outside_mM, c.chloride_conductance_s))
        total_g = sum(ion[-1] for ion in ions)
        reversal = sum(g*thermal_voltage/z*math.log(outside/inside) for _, z, inside, outside, g in ions)
        resting = (reversal-c.pump_current_a)/total_g
        vm = prefix+'.membrane_voltage'
        r.add(Component(vm, 'skin_cell_membranes', region, 'V', resting, .03))
        inputs = [vm]; weights = [-total_g/c.capacitance_f]
        bias = (reversal-c.pump_current_a)/c.capacitance_f
        for ion, z, inside, outside, g in ions:
            ident = prefix+'.'+ion+'_inside'
            r.add(Component(ident, 'skin_cell_membranes', region, 'mmol/L', inside, inside*.2))
            # A reservoir approximation: concentrations are uncertain but constant.
            # No unsupported ion mass-balance kinetics are invented here.
            r.add(Process('skin.reservoir.'+ident, (ident,), ident, 'skin_membrane_local', (0.,), 0., 0., 'fixed_intracellular_reservoir_approximation'))
            derivative = -g*thermal_voltage/(z*inside*c.capacitance_f)
            inputs.append(ident); weights.append(derivative); bias -= derivative*inside
        r.add(Process('skin.membrane.cell_'+str(i), tuple(inputs), vm, 'skin_membrane_local',
                      tuple(weights), bias, 2*(.02**2)*total_g/c.capacitance_f, mechanism))
    for i, (region, s) in enumerate(zip(surface_regions, patch.surface)):
        ident = f'skin.surface_{i}.potential'
        r.add(Component(ident, 'skin_apical_surface', region, 'V', s.pump_current_a/s.barrier_conductance_s, .02))
        rate = s.barrier_conductance_s/s.capacitance_f
        r.add(Process('skin.barrier.site_'+str(i), (ident,), ident, 'skin_barrier_local',
                      (-rate,), s.pump_current_a/s.capacitance_f, 2*.01**2*rate,
                      'physical_RC_form;illustrative_barrier_pump_and_shunt_parameters'))
    for name, edges, elements, component in (
        ('gap_junction', patch.gap_edges, patch.cells, lambda i: f'skin.cell_{i}.membrane_voltage'),
        ('extracellular', patch.surface_edges, patch.surface, lambda i: f'skin.surface_{i}.potential')):
        for edge_id, (a, b, g) in enumerate(edges):
            for src, dst in ((a, b), (b, a)):
                weight = g/elements[dst].capacitance_f
                r.add(Process(f'skin.{name}.edge_{edge_id}.to_{dst}', (component(src), component(dst)),
                              component(dst), 'skin_'+name, (weight, -weight), 0., 0.,
                              'Ohmic_conservative_pair;illustrative_conductance'))
    return r


def skin_model(patch, subject='synthetic-skin', view='combined'):
    from ihm.body import body
    targets = []
    if view not in ('combined', 'surface', 'membrane'):
        raise ValueError('skin view must be combined, surface or membrane')
    if view in ('combined', 'surface'):
        targets.extend(f'skin.surface_{i}.potential' for i in range(len(patch.surface)))
    if view in ('combined', 'membrane'):
        targets.extend(f'skin.cell_{i}.membrane_voltage' for i in range(len(patch.cells)))
    m = materialize(body(skin=patch), Request(tuple(targets), subject))
    m.provenance.update({'representation': 'heterogeneous_skin_RC_network', 'geometry': asdict(patch),
                         'voltage_reference': 'basal extracellular bath = 0 V; apical potential is negative of TEP',
                         'limitations': ['fixed extracellular bath and intracellular ion reservoirs',
                           'Nernst linearization around declared concentrations',
                           'no cell migration, differentiation, healing or morphogenetic outcome model',
                           'electrical and gap-junction graphs independent under the prescribed-bath approximation',
                           'conductance and pump parameters uncalibrated; no human predictive validation']})
    return m


def electric_field(model, patch):
    """signed edge field along a->b, and its propagated standard deviation."""
    h = np.zeros((len(patch.surface_edges), len(model.components)))
    vectors = []
    for j, (a, b, _) in enumerate(patch.surface_edges):
        displacement = np.subtract(patch.surface[b].position_m, patch.surface[a].position_m)
        length = np.linalg.norm(displacement); vectors.append(displacement/length)
        h[j, model.components.index(f'skin.surface_{a}.potential')] = 1/length
        h[j, model.components.index(f'skin.surface_{b}.potential')] = -1/length
    return {'V_per_m': h @ model.mean, 'std_V_per_m': np.sqrt(np.maximum(0, np.diag(h @ model.cov @ h.T))),
            'direction': np.asarray(vectors), 'reference': 'edge direction a->b; E = -grad(apical potential)'}
