"""Mechanical/world subcycling inside a slower neural/physiological exchange."""
from copy import deepcopy
import math
import numpy as np


def advance_body_world(plant, world, dt, state, forces, actuation, object_forces=(), max_step=.005,
                       respiratory_projector=None, final_full=True, observation_entity_ids=()):
    """Return endpoint mechanics and impulse-equivalent load quadrature ports.

User body forces follow their initial material stations. Object force targets
already contain local rigid stations or cloth node IDs. Contact forces are
recomputed against each accepted native pose; no extra body mass is introduced.
"""
    count = max(1, math.ceil(dt / max_step))
    h = dt / count
    # World samples attach to a small set of named bones. Unknown world owners
    # retain the full projection rather than assuming their observation needs.
    required = set(getattr(world, 'ids', state['entities'])) | set(observation_entity_ids)
    required.update(force['id'] for force in forces)
    anchors = []
    for force in forces:
        entity = state['entities'].get(force['id'])
        local = None
        if entity is not None:
            rotation = np.asarray(entity.get('rotation_matrix', np.eye(3)))
            local = rotation.T @ (np.asarray(force['point_m']) - entity['centroid_m'])
        anchors.append((force, local))
    loads = []
    respiratory = None
    positive_work = 0.
    for index in range(count):
        applied = []
        for force, local in anchors:
            port = deepcopy(force)
            if local is not None:
                entity = state['entities'][force['id']]
                port['point_m'] = (np.asarray(entity['centroid_m']) +
                    np.asarray(entity.get('rotation_matrix', np.eye(3))) @ local).tolist()
            applied.append(port)
        applied += world.advance(h, state['entities'],
                                 **({'object_forces': object_forces} if object_forces else {}))
        if respiratory_projector is not None:
            # These stations and forces belong to the accepted pose at the
            # beginning of this substep, not to the final exchange pose.
            sample = respiratory_projector(applied, state['entities'])
            if respiratory is None:
                respiratory = deepcopy(sample)
                respiratory.update(external_pressure_pa=0., force_ports=[], ignored_nonrespiratory_ids=[])
                respiratory['quadrature'] = {
                    'substeps': count, 'interval_s': h,
                    'basis': 'Mean generalized load at each force-generating pose; native lung volume held over exchange'}
            respiratory['external_pressure_pa'] += sample['external_pressure_pa'] / count
            for port in sample['force_ports']:
                weighted = deepcopy(port)
                for key in ('force_n', 'moment_nm'):
                    if key in weighted:
                        weighted[key] = (np.asarray(weighted[key]) / count).tolist()
                weighted['generalized_force_pa'] /= count
                weighted['quadrature_weight'] = 1. / count
                weighted['sample_time_s'] = state['time_s']
                respiratory['force_ports'].append(weighted)
            respiratory['ignored_nonrespiratory_ids'].extend(sample['ignored_nonrespiratory_ids'])
        if hasattr(plant, 'advance_observation') and not (final_full and index == count-1):
            state = plant.advance_observation(h, forces=applied, actuation=actuation, entity_ids=required)
        else:
            state = plant.advance(h, forces=applied, actuation=actuation)
        positive_work += state['positive_muscle_work_j']
        for force in applied:
            loads.append({**force, 'force_n': (np.asarray(force['force_n']) / count).tolist()})
    state['positive_muscle_work_j'] = positive_work
    audit = {'interval_s': h, 'substeps': count,
        'load_basis': 'Force/position quadrature over mechanical substeps; impulses preserved',
        'biological_validation': False}
    if respiratory is not None:
        respiratory['ignored_nonrespiratory_ids'] = list(dict.fromkeys(respiratory['ignored_nonrespiratory_ids']))
        audit['respiratory_load'] = respiratory
    return state, loads, audit
