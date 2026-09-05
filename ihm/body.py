"""organ-scale reference body; all numeric defaults are synthetic weak priors.

No coefficient here claims an empirically identified human physiological law.
Each process is a local affine approximation, replaceable by a fitted form.
"""
from ihm.registry import Registry
from ihm.fields import Component, Support
from ihm.anatomy import Partition
from ihm.topologies import Topology
from ihm.processes import Process

# id, support, region, units, center, uncertainty, relaxation time in seconds
DECLARATIONS = (
 ('neural.autonomic', 'peripheral_nerves', 'autonomic', '1', 0, 1, 30),
 ('cardiac.rate', 'heart_tissue', 'heart', 'bpm', 70, 25, 20),
 ('cardiac.stroke_volume', 'heart_tissue', 'heart', 'mL', 70, 25, 30),
 ('blood.pressure', 'vascular_tree', 'systemic', 'mmHg', 90, 30, 15),
 ('blood.flow', 'vascular_tree', 'systemic', 'L/min', 5, 2, 20),
 ('blood.oxygenation', 'vascular_tree', 'systemic', '1', .96, .08, 30),
 ('blood.volume', 'vascular_tree', 'systemic', 'L', 5, 1, 3600),
 ('respiratory.ventilation', 'airway_tree', 'lung', 'L/min', 6, 3, 20),
 ('respiratory.carbon_dioxide', 'airway_tree', 'lung', 'mmHg', 40, 15, 60),
 ('metabolic.glucose', 'interstitial', 'systemic', 'mmol/L', 5, 2, 1800),
 ('metabolic.lactate', 'interstitial', 'systemic', 'mmol/L', 1, 1, 300),
 ('metabolic.oxygen_consumption', 'organ_tissue', 'muscle', 'mL/min', 250, 150, 60),
 ('endocrine.insulin', 'vascular_tree', 'systemic', 'mU/L', 8, 8, 600),
 ('endocrine.cortisol', 'vascular_tree', 'systemic', 'nmol/L', 300, 200, 3600),
 ('renal.filtration', 'nephron', 'kidney', 'mL/min', 100, 50, 120),
 ('renal.urine_flow', 'nephron', 'kidney', 'mL/min', 1, 1, 300),
 ('extracellular.sodium', 'interstitial', 'systemic', 'mmol/L', 140, 10, 3600),
 ('extracellular.potassium', 'interstitial', 'systemic', 'mmol/L', 4, 1, 1800),
 ('hepatic.clearance', 'organ_tissue', 'liver', 'L/min', 1, .8, 600),
 ('digestive.glucose_delivery', 'gut_lumen', 'gut', 'mmol/min', 0, 3, 900),
 ('immune.crp', 'organ_tissue', 'systemic', 'mg/L', 2, 10, 14400),
 ('thermal.core', 'organ_tissue', 'core', 'degC', 37, 1, 900),
 ('thermal.skin', 'body_surface', 'skin', 'degC', 33, 3, 300),
 ('effector.activation', 'motor_units', 'muscle', '1', .05, .2, 2),
 ('effector.fatigue', 'motor_units', 'muscle', '1', .05, .2, 300),
 ('mechanical.force', 'musculoskeletal', 'muscle', 'N', 10, 100, 1),
 ('mechanical.velocity', 'musculoskeletal', 'limb', 'm/s', 0, 1, 1),
 ('structural.muscle_mass', 'organ_tissue', 'muscle', 'kg', 25, 10, 1e7),
 ('device.heart_rate', 'sensor_array', 'wrist', 'bpm', 70, 30, 1),
 ('device.glucose', 'sensor_array', 'skin', 'mmol/L', 5, 2, 300),
 ('device.temperature', 'sensor_array', 'skin', 'degC', 33, 3, 5),
)

# Input -> output sensitivity at the illustrative reference state.
COUPLINGS = (
 ('neural.autonomic', 'cardiac.rate', 15, 'autonomic_pathway'),
 ('cardiac.rate', 'blood.flow', .05, 'vascular_transport'),
 ('cardiac.stroke_volume', 'blood.flow', .05, 'vascular_transport'),
 ('blood.flow', 'blood.pressure', 8, 'vascular_transport'),
 ('blood.volume', 'blood.pressure', 5, 'vascular_transport'),
 ('blood.pressure', 'renal.filtration', .3, 'vascular_transport'),
 ('renal.filtration', 'renal.urine_flow', .01, 'renal_exchange'),
 ('renal.urine_flow', 'blood.volume', -.02, 'renal_exchange'),
 ('respiratory.ventilation', 'blood.oxygenation', .005, 'gas_exchange'),
 ('respiratory.ventilation', 'respiratory.carbon_dioxide', -2, 'gas_exchange'),
 ('metabolic.oxygen_consumption', 'respiratory.ventilation', .005, 'gas_exchange'),
 ('blood.oxygenation', 'metabolic.lactate', -3, 'metabolic_exchange'),
 ('digestive.glucose_delivery', 'metabolic.glucose', .2, 'metabolic_exchange'),
 ('endocrine.insulin', 'metabolic.glucose', -.03, 'endocrine_transport'),
 ('metabolic.glucose', 'endocrine.insulin', 2, 'endocrine_transport'),
 ('endocrine.cortisol', 'metabolic.glucose', .001, 'endocrine_transport'),
 ('renal.filtration', 'extracellular.potassium', -.002, 'renal_exchange'),
 ('renal.urine_flow', 'extracellular.sodium', .5, 'renal_exchange'),
 ('effector.activation', 'metabolic.oxygen_consumption', 300, 'metabolic_exchange'),
 ('effector.activation', 'effector.fatigue', .2, 'mechanical_link'),
 ('effector.activation', 'mechanical.force', 300, 'mechanical_link'),
 ('effector.fatigue', 'mechanical.force', -20, 'mechanical_link'),
 ('structural.muscle_mass', 'mechanical.force', 1, 'mechanical_link'),
 ('mechanical.force', 'mechanical.velocity', .001, 'mechanical_link'),
 ('metabolic.oxygen_consumption', 'thermal.core', .001, 'thermal_exchange'),
 ('immune.crp', 'thermal.core', .01, 'thermal_exchange'),
 ('thermal.core', 'thermal.skin', .5, 'thermal_exchange'),
 ('cardiac.rate', 'device.heart_rate', 1, 'device_coupling'),
 ('metabolic.glucose', 'device.glucose', 1, 'device_coupling'),
 ('thermal.skin', 'device.temperature', 1, 'device_coupling'),
)


from ihm.fields.systems import declarations, COUPLINGS as SYSTEM_COUPLINGS

DECLARATIONS = DECLARATIONS + declarations()
COUPLINGS = COUPLINGS + SYSTEM_COUPLINGS


def body(overrides=(), skin=None, fluids=None, cardiopulmonary=None):
    r = Registry()
    for support in sorted({d[1] for d in DECLARATIONS}):
        regions = tuple(sorted({d[2] for d in DECLARATIONS if d[1] == support}))
        r.add(Support(support, 'named_compartments', regions))
    for ident, support, region, unit, mean, std, tau in DECLARATIONS:
        r.add(Component(ident, support, region, unit, mean, std))
    regions = sorted({d[2] for d in DECLARATIONS})
    r.add(Partition('organ_partition', {x: {x: 1.} for x in regions}))
    r.add(Partition('functional_partition', {
        'heart': {'cardiovascular': 1.}, 'lung': {'respiratory': 1.},
        'muscle': {'locomotor': .6, 'metabolic': .4},
        'kidney': {'excretory': .5, 'fluid_balance': .5}}))
    edges = {'local': {(x, x) for x in regions}}
    for a, b, _, t in COUPLINGS:
        edges.setdefault(t, set()).add((r.components[a].region, r.components[b].region))
    for t, pairs in edges.items():
        r.add(Topology(t, tuple(sorted(pairs))))
    times = {d[0]: d[-1] for d in DECLARATIONS}
    processes = {}
    for c in r.components.values():
        tau = times[c.id]
        p = Process('relax.' + c.id, (c.id,), c.id, 'local', (-1/tau,), c.mean/tau,
                    2*c.std**2/tau)
        processes[p.id] = p
    for a, b, sensitivity, topology in COUPLINGS:
        weight = sensitivity / times[b]
        p = Process('couple.' + a + '.to.' + b, (a,), b, topology, (weight,),
                    -weight*r.components[a].mean, 0.)
        processes[p.id] = p
    seen = set()
    for p in overrides:
        if p.id not in processes or p.id in seen:
            raise ValueError(f'unknown or duplicate process override: {p.id}')
        old = processes[p.id]
        if (old.inputs, old.output, old.topology) != (p.inputs, p.output, p.topology):
            raise ValueError('a fitted form must preserve process ontology')
        processes[p.id] = p; seen.add(p.id)
    for p in processes.values():
        r.add(p)
    if skin is not None:
        from ihm.materialize.skin import register_skin
        register_skin(r, skin)
    if fluids is not None:
        from ihm.materialize.fluids import register_fluids
        register_fluids(r, fluids)
    if cardiopulmonary is not None:
        from ihm.processes.cardiopulmonary import register_cardiopulmonary
        register_cardiopulmonary(r, cardiopulmonary)
    return r.seal()
