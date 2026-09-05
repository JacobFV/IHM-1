"""Three-mode thoracic mechanics driven by native gas volume or pressure.

The native engine owns gas storage. This constrained viscoelastic substructure
owns rib, anterior chest and diaphragm motion; it is not a contact/FEM lung.
All dynamics use SI units and the anatomical X-left/Y-superior/Z-anterior frame.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import numpy as np


def build(root):
    """Derive geometry and bindings; retain uncalibrated parameter provenance."""
    root = Path(root)
    path = root/'data/derived/canonical/anatomy.json'
    anatomy = json.loads(path.read_text())
    lungs = [e for e in anatomy['entities'] if 'lobe' in e['name'].lower() and 'lung' in e['name'].lower()]
    if len(lungs) != 5:
        raise ValueError('Expected the five canonical lung lobes')
    low = np.min([e['bounds_m']['min'] for e in lungs], axis=0)
    high = np.max([e['bounds_m']['max'] for e in lungs], axis=0)
    center = (low+high)/2
    width, height, depth = high-low
    a, b = width/2, depth/2
    # Elliptical-cylinder linear volume Jacobian. q = lateral radius increment,
    # anterior depth increment, caudal diaphragm excursion; posterior wall fixed.
    areas = np.array([math.pi*b*height, math.pi*a*height/2, math.pi*a*b])
    fractions = np.array([.18, .17, .65])
    compliance = .0002/98.0665
    stiffness = np.diag(areas**2/(compliance*fractions))
    # A positive coupling energy penalizes differences from the generic mode
    # partition without changing total static compliance or its equilibrium.
    for i, j in [(0, 1), (1, 2), (0, 2)]:
        row = np.zeros(3)
        row[i], row[j] = areas[i]*fractions[j], -areas[j]*fractions[i]
        stiffness += 2./compliance * np.outer(row, row)
    bindings = []
    for entity in anatomy['entities']:
        name = entity['name'].lower()
        kind = None
        if entity in lungs:
            kind = 'lung'
        elif name == 'diaphragm':
            kind = 'diaphragm'
        elif name.endswith(' rib') and entity['role'] == 'rigid_bone':
            kind = 'rib'
        elif name in ['body of sternum', 'manubrium', 'xiphoid process']:
            kind = 'sternum'
        elif name.endswith('thoracic vertebra') and 'intervertebral' not in name and entity['role'] == 'rigid_bone':
            kind = 'posterior_support'
        if kind:
            position = np.asarray(entity['centroid_m'])
            basis = np.zeros((3, 3))
            if kind == 'rib':
                # Rigid centroid translation only. No implied rib articulation.
                basis[0, 0] = np.clip((position[0]-center[0])/a, -1, 1)
                basis[2, 1] = np.clip((position[2]-low[2])/depth, 0, 1)
                basis[1, 1] = .25 * np.clip((position[2]-low[2])/depth, 0, 1)
            elif kind == 'sternum':
                basis[2, 1], basis[1, 1] = 1., .25
            elif kind == 'diaphragm':
                basis[1, 2] = -1.
            elif kind == 'lung':
                basis[0, 0] = (position[0]-center[0])/a
                basis[2, 1] = (position[2]-low[2])/depth
                basis[1, 2] = -np.clip((high[1]-position[1])/height, 0, 1)
            bindings.append({'entity_id': entity['id'], 'kind': kind, 'centroid_m': entity['centroid_m'], 'translation_basis': basis.tolist()})
    skin_ids = [e['id'] for e in anatomy['entities'] if e['role'] == 'skin']
    skin = {'entity_ids': skin_ids, 'center_m': center.tolist(), 'bounds_m': {'min': [low[0]-.04, low[1]-.14, low[2]-.04], 'max': [high[0]+.04, high[1]+.07, high[2]+.04]}, 'thorax_y_m': [float(low[1]), float(high[1])], 'reference_radii_m': [float(a), float(b)], 'posterior_z_m': float(low[2]), 'reference_height_m': float(height), 'basis': 'thoracic-smoothstep-v1'}
    return {'schema_version': 1, 'model_id': 'ihm-reduced-thoracic-mechanics-v1', 'anatomy_sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'frame': anatomy['frame'], 'bindings': bindings, 'skin_field': skin, 'volume_jacobian_m2': areas.tolist(), 'stiffness_n_m': stiffness.tolist(), 'damping_ns_m': (.08*stiffness).tolist(), 'effective_mass_kg': [1.5, 1., .5], 'reference_dimensions_m': {'transverse': float(width), 'anteroposterior': float(depth), 'height': float(height)}, 'priors': {'chest_wall_compliance_m3_pa': {'value': compliance, 'range_l_cmH2O': [.1, .3], 'evidence_kind': 'generic_prior', 'calibrated': False, 'note': 'Local linear resting compliance; range is an exploration bracket, not a fitted confidence interval.'}, 'mode_volume_fractions': {'value': fractions.tolist(), 'evidence_kind': 'synthesized_prior', 'calibrated': False}, 'viscoelastic_relaxation_s': {'value': .08, 'range': [.02, .2], 'evidence_kind': 'numerical_biomechanical_prior', 'calibrated': False}, 'effective_mass_kg': {'value': [1.5, 1., .5], 'evidence_kind': 'synthesized_reduced_mass', 'calibrated': False}}, 'sources': [{'title': 'Marconi and De Lazzari (2020), In silico study of airway/lung mechanics in normal human breathing', 'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7239037/', 'doi': '10.1016/j.matcom.2020.05.014', 'supports': 'Pressure-volume compliance; lung and chest-wall recoil add; local linear approximation only, not their full nonlinear model.'}, {'title': 'Aliverti et al. (1999), Chest wall kinematics and respiratory muscle action in walking healthy humans', 'url': 'https://pubmed.ncbi.nlm.nih.gov/10484561/', 'supports': 'Regional rib-cage/abdominal volume accounting; no fitted transfer to this anatomy.'}, {'title': 'Boussuges et al. (2009), Diaphragmatic motion studied by M-mode ultrasonography', 'url': 'https://journal.chestnet.org/article/S0012-3692%2809%2960127-6/abstract', 'supports': 'Observed quiet-breathing diaphragm excursion provides an order-of-magnitude comparison, not calibration.'}], 'limitations': ['Three displacement modes, linearized cavity volume, rigid translated ribs with constrained orientations; no costovertebral joint solve.', 'The regional skin field and anisotropic lung shape are kinematic reductions of the solved thoracic modes, not volumetric FEM or pleural contact.', 'Native gas volume is prescribed; the constraint pressure is a mechanical reaction, not a measured pleural pressure or a second gas state.', 'Generic compliance, effective mass, damping and regional partition are uncalibrated priors.', 'Thoracic vertebrae are prescribed reference supports; support reaction and boundary work must be audited by whole-body mechanics.', 'Lung/skin/self contact, tissue sliding, abdominal organ contact, and bidirectional respiratory physiology feedback remain unresolved.']}


def deform_skin(points, field):
    """Apply the documented regional basis to reference vertices without mutation."""
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 3 or not np.isfinite(p).all():
        raise ValueError('Expected finite Nx3 skin vertices')
    c = np.asarray(field['center_m'])
    low, high = np.asarray(field['bounds_m']['min']), np.asarray(field['bounds_m']['max'])
    a, b = field['reference_radii_m']
    # Coordinates relative to center make the renderer field translation invariant.
    y0, y1 = field['thorax_y_m']
    # Stored thoracic offsets remain invariant if the complete field is translated.
    if 'thorax_y_offsets_m' in field:
        y0, y1 = np.array(field['thorax_y_offsets_m']) + c[1]
    ramp = lambda x: np.clip(x, 0, 1)**2*(3-2*np.clip(x, 0, 1))
    vertical = ramp((p[:, 1]-low[1])/(y0-low[1])) * ramp((high[1]-p[:, 1])/(high[1]-y1))
    lateral = 1-ramp((np.abs(p[:, 0]-c[0])-a)/(max(high[0]-c[0], c[0]-low[0])-a))
    weight = vertical*lateral
    q = np.asarray(field['displacement_m'])
    out = p.copy()
    out[:, 0] += weight*q[0]*np.clip((p[:, 0]-c[0])/a, -1, 1)
    anterior = np.clip((p[:, 2]-c[2]+b)/(2*b), 0, 1)
    out[:, 2] += weight*q[1]*anterior
    return out


class BodyRespiration:
    """Backward Euler constrained damped mechanics with discrete energy audit.

    M qdd + D qd + K q = f + A*p, and A.q = native deltaV in volume mode.
    Pressure mode is provided for model checks and explicit alternative forcing.
    """
    def __init__(self, payload):
        self.payload = payload
        self.a = np.asarray(payload['volume_jacobian_m2'], float)
        self.k = np.asarray(payload['stiffness_n_m'], float)
        self.d = np.asarray(payload['damping_ns_m'], float)
        self.m = np.diag(payload['effective_mass_kg'])
        for matrix in [self.k, self.d, self.m]:
            if matrix.shape != (3, 3) or not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.T) or np.linalg.eigvalsh(matrix).min() <= 0:
                raise ValueError('Positive symmetric thoracic operators required')
        if self.a.shape != (3,) or not np.isfinite(self.a).all() or np.min(self.a) <= 0:
            raise ValueError('Positive finite volume Jacobian required')
        self.q, self.v = np.zeros(3), np.zeros(3)
        self.time = self.work = self.dissipation = self.numerical_dissipation = 0.
        self.max_substep = .01
        self._operators = {}
        self.lung_ids = {b['entity_id'] for b in payload['bindings'] if b['kind'] == 'lung'}

    @classmethod
    def from_dict(cls, payload):
        return cls(payload)

    def step(self, dt, drivers=None):
        drivers = {} if drivers is None else drivers
        if not isinstance(drivers, dict) or set(drivers)-{'pressure_pa', 'volume_change_m3', 'generalized_forces_n', 'lung_volume_ratios'}:
            raise ValueError('Unknown respiratory driver')
        if isinstance(dt, bool) or not math.isfinite(float(dt)) or not 0 < float(dt) <= 1:
            raise ValueError('dt must be finite in (0,1] seconds')
        dt = float(dt)
        if 'pressure_pa' in drivers and 'volume_change_m3' in drivers:
            raise ValueError('Prescribe pressure or volume, not both')
        for key in ['pressure_pa', 'volume_change_m3']:
            if key in drivers and (isinstance(drivers[key], bool) or not math.isfinite(float(drivers[key]))):
                raise ValueError('Finite scalar respiratory boundary required')
        f = np.asarray(drivers.get('generalized_forces_n', [0., 0., 0.]), float)
        if f.shape != (3,) or not np.isfinite(f).all():
            raise ValueError('Expected finite three-mode generalized force')
        ratios = drivers.get('lung_volume_ratios', {})
        if not isinstance(ratios, dict) or set(ratios)-self.lung_ids:
            raise ValueError('Unknown lung volume binding')
        if any(isinstance(r, bool) or not math.isfinite(float(r)) or not .25 <= float(r) <= 4 for r in ratios.values()):
            raise ValueError('Positive lung ratio in .25..4 required')
        volume_mode = 'volume_change_m3' in drivers
        initial_volume = float(self.a@self.q)
        target = float(drivers.get('volume_change_m3', initial_volume))
        pressure = float(drivers.get('pressure_pa', 0.))
        n = max(1, math.ceil(dt/self.max_substep))
        h = dt/n
        key = round(h, 12)
        if key not in self._operators:
            inverse = np.linalg.inv(self.m+h*self.d+h*h*self.k)
            self._operators[key] = (inverse, inverse@self.a)
        inverse, compliance_direction = self._operators[key]
        residual = 0.
        q, v = self.q.copy(), self.v.copy()
        work, dissipation, numerical_dissipation = self.work, self.dissipation, self.numerical_dissipation
        for i in range(n):
            old_q, old_v = q.copy(), v.copy()
            rhs = self.m@old_v+h*(f-self.k@old_q)
            velocity_free = inverse@rhs
            if volume_mode:
                wanted = initial_volume+(target-initial_volume)*(i+1)/n
                pressure = (wanted-self.a@(old_q+h*velocity_free))/(h*h*(self.a@compliance_direction))
            v = velocity_free+h*compliance_direction*pressure
            q = old_q+h*v
            dq, dv = q-old_q, v-old_v
            work += float((f+self.a*pressure)@dq)
            dissipation += float(h*v@self.d@v)
            numerical_dissipation += float(.5*dv@self.m@dv+.5*dq@self.k@dq)
            residual = max(residual, float(np.linalg.norm(self.m@dv/h+self.d@v+self.k@q-f-self.a*pressure)))
        dims = self.payload['reference_dimensions_m']
        if not np.isfinite(q).all() or np.any(np.abs(q) > .25*np.array([dims['transverse']/2, dims['anteroposterior'], dims['height']])):
            raise ValueError('Respiratory forcing exceeded the small-deformation domain')
        # Commit only after the candidate solution and derived shapes are safe.
        self.q, self.v = q, v
        self.work, self.dissipation, self.numerical_dissipation = work, dissipation, numerical_dissipation
        self.time += dt
        elastic, kinetic = float(.5*q@self.k@q), float(.5*v@self.m@v)
        field = dict(self.payload['skin_field'])
        # Preserve relative vertical anchors for origin-invariant field evaluation.
        # Body geometry stores offsets explicitly at build time in later versions;
        # the centered reference height defines the same thoracic band here.
        field['thorax_y_offsets_m'] = [-dims['height']/2, dims['height']/2]
        field['displacement_m'] = self.q.tolist()
        field.update(lateral_expansion_m=float(self.q[0]), anterior_expansion_m=float(self.q[1]), diaphragm_descent_m=float(self.q[2]))
        translations, deformations, entities = {}, {}, {}
        for binding in self.payload['bindings']:
            entity_id, kind = binding['entity_id'], binding['kind']
            translation = np.asarray(binding['translation_basis'])@self.q
            state = {'translation_m': translation.tolist()}
            if kind != 'lung' or entity_id in ratios:
                translations[entity_id] = state['translation_m']
            if kind == 'diaphragm':
                # Tissue volume preserved; dome descent is its centroid DOF.
                shorten = math.exp(-self.q[2]/dims['height'])
                deformation = np.diag([shorten, shorten**-2, shorten])
                deformations[entity_id] = state['deformation_gradient'] = deformation.tolist()
            elif kind == 'lung' and entity_id in ratios:
                shape = np.array([1+2*self.q[0]/dims['transverse'], 1+self.q[2]/dims['height'], 1+self.q[1]/dims['anteroposterior']])
                if np.min(shape) <= 0:
                    raise ValueError('Thoracic deformation exceeded positive shape range')
                shape *= (float(ratios[entity_id])/float(np.prod(shape)))**(1/3)
                deformations[entity_id] = state['deformation_gradient'] = np.diag(shape).tolist()
            entities[entity_id] = state
        volume = float(self.a@self.q)
        return {'schema_version': 1, 'model_id': self.payload['model_id'], 'time_s': self.time, 'displacement_m': self.q.tolist(), 'velocity_m_s': self.v.tolist(), 'volume_change_m3': volume, 'constraint_pressure_pa': pressure, 'diaphragm_descent_m': float(self.q[2]), 'chest_dimensions_m': {'transverse': dims['transverse']+2*float(self.q[0]), 'anteroposterior': dims['anteroposterior']+float(self.q[1]), 'height': dims['height']+float(self.q[2])}, 'entities': entities, 'mechanics_drivers': {'prescribed_translations_m': translations, 'deformation_gradients': deformations}, 'skin_field': field, 'audit': {'force_balance_residual_n': residual, 'volume_constraint_residual_m3': volume-target if volume_mode else None, 'elastic_energy_j': elastic, 'kinetic_energy_j': kinetic, 'accumulated_boundary_work_j': self.work, 'accumulated_dissipation_j': self.dissipation, 'implicit_numerical_dissipation_j': self.numerical_dissipation, 'energy_balance_residual_j': elastic+kinetic+self.dissipation+self.numerical_dissipation-self.work, 'volume_scope': 'exact linearized modal cavity constraint; not a closed surface mesh integral', 'pressure_scope': 'reaction incremental to reference, not observed pleural pressure', 'biological_validation': False}}
