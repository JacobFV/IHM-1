"""Resolve mechanically admissible traction fields to audited native boundaries.

Registered cells are non-overlapping source surface elements, not fluid stores.
Only complete uniform whole-Skin normal loading has a native pressure port.
"""
from copy import deepcopy
import math
from .body_exchange import ORGANS, number


def vector(value, label):
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError('Expected three components: ' + label)
    return [number(v, label) for v in value]


class RegionalDrainageBoundary:
    def __init__(self, cells, owner_scope_receipts=None):
        self.cells = {}
        source_elements = set()
        for raw in cells:
            row = deepcopy(raw)
            key = row['id']
            if key in self.cells:
                raise ValueError('Duplicate source surface cell')
            if row['native_owner'] not in [o + '.extracellular' for o in ORGANS]:
                raise ValueError('Unknown extracellular owner')
            identity = row.get('source_identity', {})
            mesh, triangle = identity.get('mesh_sha256'), identity.get('triangle')
            if not isinstance(mesh, str) or not mesh or isinstance(triangle, bool) or not isinstance(triangle, int) or triangle < 0:
                raise ValueError('Source mesh hash and triangle identity required')
            area = number(row['area_m2'], 'registered source area', True)
            if area <= 0:
                raise ValueError('Positive source cell area required')
            normal = vector(row['outward_normal'], 'outward normal')
            if not math.isclose(math.sqrt(math.fsum(x*x for x in normal)), 1., abs_tol=1e-12, rel_tol=0):
                raise ValueError('Unit outward normal required')
            # A source element may only occur once, even under another display ID.
            if (mesh, triangle) in source_elements:
                raise ValueError('Duplicate source element ownership')
            source_elements.add((mesh, triangle))
            self.cells[key] = row
        self.scopes = deepcopy(owner_scope_receipts or {})

    def map(self, tractions, native_snapshot):
        """Return command candidates without sending commands or changing state.

        Every traction has cell_id, area_m2 and traction_pa (outward convention).
        Missing cells are unknown, not implicitly unloaded. Whole-owner scope is
        an explicit caller receipt, never inferred from surface-area weights.
        """
        rows = []
        seen = set()
        for load in tractions:
            key = load['cell_id']
            if key not in self.cells or key in seen:
                raise ValueError('Unknown or duplicate loaded cell')
            seen.add(key)
            cell = self.cells[key]
            area = number(load['area_m2'], 'loaded area', True)
            if not math.isclose(area, cell['area_m2'], abs_tol=0., rel_tol=1e-12):
                raise ValueError('Loaded area differs from registered source area')
            traction = vector(load['traction_pa'], 'traction Pa')
            normal = cell['outward_normal']
            pressure = -math.fsum(t*n for t, n in zip(traction, normal))
            tangent = [t + pressure*n for t, n in zip(traction, normal)]
            rows.append({'cell_id': key, 'native_owner': cell['native_owner'],
                         'source_identity': deepcopy(cell['source_identity']),
                         'area_m2': area, 'outward_normal': list(normal),
                         'traction_pa': traction, 'compressive_pressure_pa': pressure,
                         'tangential_traction_pa': tangent,
                         'force_n': [area*t for t in traction],
                         'normal_load_n': area*pressure, 'independent_store': False})
        commands, unresolved, audits = [], [], []
        for owner in sorted({r['native_owner'] for r in rows}):
            selected = [r for r in rows if r['native_owner'] == owner]
            registered = [r for r in self.cells.values() if r['native_owner'] == owner]
            area = math.fsum(r['area_m2'] for r in selected)
            normal_load = math.fsum(r['normal_load_n'] for r in selected)
            pressures = [r['compressive_pressure_pa'] for r in selected]
            reasons = []
            if len(selected) != len(registered):
                reasons.append('incomplete_registered_surface_loading')
            scope = self.scopes.get(owner, {})
            if scope.get('complete_native_owner_surface') is not True or not scope.get('source_identity'):
                reasons.append('whole_native_owner_surface_scope_unresolved')
            if any(p < 0 for p in pressures):
                reasons.append('tensile_pressure_unsupported')
            if any(p > 5000 for p in pressures):
                reasons.append('pressure_exceeds_native_command_limit_5000_pa')
            if any(math.sqrt(math.fsum(t*t for t in r['tangential_traction_pa'])) > 1e-10 for r in selected):
                reasons.append('tangential_traction_has_no_native_fluid_port')
            if any(not math.isclose(p, pressures[0], rel_tol=1e-12, abs_tol=1e-10) for p in pressures):
                reasons.append('nonuniform_regional_pressure_has_no_native_port')
            if owner != 'Skin.extracellular':
                reasons.append('owner_has_no_audited_external_pressure_adapter')
            values = native_snapshot.get('values', {})
            installed = all(isinstance(values.get(k), (int, float)) and not isinstance(values.get(k), bool) and math.isfinite(values[k]) for k in (
                'tissue.compression.Skin.requested_pa', 'tissue.compression.Skin.applied_pa'))
            if owner == 'Skin.extracellular' and not installed:
                reasons.append('native_skin_boundary_installation_unobserved')
            audits.append({'native_owner': owner, 'loaded_area_m2': area,
                           'registered_area_m2': math.fsum(r['area_m2'] for r in registered),
                           'normal_load_n': normal_load,
                           'force_n': [math.fsum(r['force_n'][i] for r in selected) for i in range(3)],
                           'owner_scope_receipt': deepcopy(scope)})
            if reasons:
                unresolved.append({'native_owner': owner, 'reasons': reasons,
                                   'required_extension': 'native_regional_interstitial_pressure_boundaries'})
            else:
                # This is uniformity-checked quadrature, not homogenization of
                # nonuniform regional loading into a fictitious regional solve.
                pressure = normal_load / area
                commands.append({'native_owner': owner, 'method': 'skin_compression',
                                 'pressure_pa': pressure, 'scope': 'uniform_whole_native_owner',
                                 'boundary_node': 'IHMSkinExternalPressure',
                                 'drive_path': 'IHMGroundToSkinExternalPressure',
                                 'compliance_path': 'IHMSkinE3ToExternalPressure',
                                 'normal_load_residual_n': pressure*area-normal_load,
                                 'regional_native_state_created': False})
        return {'schema': 'regional_drainage_boundary_v1',
                'native_time_s': number(native_snapshot.get('time_s'), 'native time', True),
                'regional_tractions': rows, 'owner_audits': audits,
                'native_command_candidates': commands, 'unresolved': unresolved,
                'unloaded_cell_ids': sorted(set(self.cells)-seen),
                'native_commands_sent': False,
                'regional_drainage_predicted': False}
