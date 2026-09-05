"""trace symbolic targets before allocating any state."""
from dataclasses import dataclass, asdict
import numpy as np
from ihm.materialize.model import Model

LIBRARY = {
    'vascular': ('vascular.radius', 'vascular.compliance', 'capillary.filtration', 'venous.pressure'),
    'immune': ('immune.crp', 'lymph_node.activated_t_cells', 'immune.neutrophils'),
    'hematopoietic': ('hematopoietic.erythrocyte_production', 'blood.hematocrit', 'blood.platelets'),
    'digestive': ('digestive.glucose_delivery', 'digestive.water_absorption', 'digestive.motility'),
    'hepatic': ('hepatic.clearance', 'hepatic.glucose_production', 'hepatic.bile_flow'),
    'endocrine': ('endocrine.insulin', 'endocrine.cortisol', 'endocrine.thyroxine', 'endocrine.aldosterone'),
    'nervous': ('neural.autonomic', 'neural.sensory_activity', 'csf.pressure'),
    'skeletal': ('bone.mineral_density', 'joint.cartilage_stress', 'tendon.tension'),
    'reproductive': ('reproductive.estradiol', 'reproductive.testosterone', 'uterine.endometrial_thickness'),
    'cardiovascular': ('blood.pressure', 'blood.flow', 'device.heart_rate'),
    'respiratory': ('blood.oxygenation', 'respiratory.carbon_dioxide'),
    'glucose': ('metabolic.glucose', 'device.glucose'),
    'renal': ('renal.filtration', 'renal.urine_flow', 'extracellular.potassium'),
    'thermal': ('thermal.core', 'device.temperature'),
    'lymphatic': ('lymphatic.flow', 'interstitial.volume', 'lymph_node.activated_t_cells'),
    'blood': ('blood.hematocrit', 'blood.viscosity', 'blood.oxygenation'),
    'integumentary': ('integumentary.transepithelial_potential', 'integumentary.barrier_integrity', 'sweat.flow'),
    'movement': ('mechanical.force', 'mechanical.velocity'),
}


@dataclass(frozen=True)
class Request:
    targets: tuple[str, ...]
    subject: str
    max_states: int = 512

    def __post_init__(self):
        if not self.targets or not self.subject or not isinstance(self.max_states, int) or not 0 < self.max_states <= 2048:
            raise ValueError('targets, subject and state budget (1..2048) required')
        if len(set(self.targets)) != len(self.targets):
            raise ValueError('duplicate targets')


def materialize(registry, request, population_prior=None):
    if not registry._sealed:
        raise ValueError('seal the ontology before materializing')
    required = set(request.targets)
    if required - registry.components.keys():
        raise ValueError(f'unknown targets: {required - registry.components.keys()}')
    selected = {}
    while True:
        previous = len(required)
        for p in registry.processes.values():
            if p.output in required:
                selected[p.id] = p; required.update(p.inputs)
        if len(required) > request.max_states:
            raise ValueError('materialization exceeds state budget')
        if len(required) == previous:
            break
    components = tuple(sorted(required)); idx = {c: i for i, c in enumerate(components)}
    n = len(idx); a = np.zeros((n, n)); b = np.zeros(n); q = np.zeros((n, n))
    for p in selected.values():
        if p.form != 'affine':
            raise ValueError('nonlinear cycle processes require cardiopulmonary_model()')
        i = idx[p.output]
        for c, weight in zip(p.inputs, p.weights):
            a[i, idx[c]] += weight
        b[i] += p.bias; q[i, i] += p.noise
    cs = [registry.components[c] for c in components]
    provenance = {'representation': 'compartment_scalar_gaussian',
                  'validated_biology': False,
                  'limitation': 'illustrative affine dynamics; no parameter-posterior uncertainty',
                  'components': [asdict(c) for c in cs],
                  'processes': [asdict(p) for p in selected.values()],
                  'topologies': [asdict(registry.topologies[t]) for t in sorted({p.topology for p in selected.values()})],
                  'anatomy': [{'id': v.id, 'memberships': {k: dict(w) for k, w in v.memberships.items()}} for v in registry.anatomy.values()], 'evidence': []}
    model = Model(components, request.targets, request.subject,
                 np.array([c.mean for c in cs]), np.diag([c.std**2 for c in cs]),
                 a, b, q, provenance)

    if population_prior is not None:
        from ihm.forge.population import initialize_population
        initialize_population(model,population_prior)
    return model
