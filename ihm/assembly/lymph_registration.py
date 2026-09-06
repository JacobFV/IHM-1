"""Source-only exclusive skin inventories; no anatomical or native assignment inference."""
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np


class SkinTerritoryInventory:
    """One triangle has one material owner, including an explicit unresolved owner.

    Named masks are submitted annotations with a provenance receipt, not inferred
    drainage labels. This validates geometry/accounting only, never their anatomy.
    """
    def __init__(self, geometry_path, expected_sha256, regions=()):
        path = Path(geometry_path)
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != expected_sha256:
            raise ValueError('Changed skin geometry')
        geometry = json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)
        xyz = np.asarray(geometry['positions'], dtype=float).reshape(-1, 3)
        triangles = np.asarray(geometry['indices'])
        if triangles.dtype.kind not in 'iu' or triangles.size % 3:
            raise ValueError('Invalid triangle indices')
        triangles = triangles.reshape(-1, 3)
        if not len(triangles) or triangles.min() < 0 or triangles.max() >= len(xyz) or not np.isfinite(xyz).all():
            raise ValueError('Invalid surface')
        self.area = np.linalg.norm(np.cross(xyz[triangles[:, 1]] - xyz[triangles[:, 0]],
                                           xyz[triangles[:, 2]] - xyz[triangles[:, 0]]), axis=1) / 2
        if np.any(self.area <= 0):
            raise ValueError('Degenerate skin triangle')
        self.geometry_sha256 = expected_sha256
        self.owners = np.full(len(triangles), 'unresolved', dtype=object)
        self.ids = ['unresolved']
        self.annotations = []
        occupied = set()
        for region in regions:
            name = region['id']
            if not isinstance(name, str) or not name or name in self.ids:
                raise ValueError('Duplicate or invalid territory')
            if region.get('native_owner') is not None:
                raise ValueError('Source registration cannot bind an installed engineering native owner')
            if not region.get('annotation_receipt') or region.get('geometry_sha256') != expected_sha256:
                raise ValueError('Missing annotation or geometry receipt')
            ids = region['triangle_ids']
            if not ids or any(isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < len(triangles) for i in ids):
                raise ValueError('Invalid triangle membership')
            if len(set(ids)) != len(ids) or occupied.intersection(ids):
                raise ValueError('Duplicate material ownership')
            occupied.update(ids)
            self.owners[ids] = name
            self.ids.append(name)
            self.annotations.append({'id': name, 'annotation_receipt': region['annotation_receipt'],
                                     'anatomical_validity': 'not_established_by_inventory_validation'})
        self.areas = {name: float(self.area[self.owners == name].sum()) for name in self.ids}

    def report(self):
        total = float(self.area.sum())
        return {'schema': 'skin_territory_inventory_v1', 'geometry_sha256': self.geometry_sha256,
                'triangle_count': len(self.area), 'area_m2': total,
                'regions': [{'id': name, 'triangle_count': int(np.sum(self.owners == name)),
                             'area_m2': area, 'area_fraction': area / total,
                             'native_owner': None} for name, area in self.areas.items()],
                'exclusive_complete_material_inventory': True,
                'anatomical_drainage_registration_validated': False,
                'annotations': self.annotations, 'native_commands': []}

    def reduce_contacts(self, contacts):
        """Conserve source-facet force and area; pressure is an engineering reduction.

        normal_pressure_pa is a supplied compressive scalar; force_n preserves the
        full supplied vector. This does not infer orientation or a depth pressure.
        Each contact represents a disjoint clipped facet patch; same-triangle
        multiple patches require a future geometric overlap proof and fail closed.
        """
        seen = set()
        rows = {name: {'contact_area_m2': 0., 'normal_force_n': 0., 'force_n': np.zeros(3)}
                for name in self.ids}
        for contact in contacts:
            i = contact['triangle_id']
            if isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < len(self.area) or i in seen:
                raise ValueError('Duplicate or invalid contact triangle')
            seen.add(i)
            area, pressure = contact['area_m2'], contact['normal_pressure_pa']
            if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not np.isfinite(x) for x in [area, pressure]):
                raise ValueError('Invalid contact scalar')
            if area <= 0 or area > self.area[i] * (1 + 1e-12) or pressure < 0:
                raise ValueError('Contact exceeds source area or tensile pressure')
            force = np.asarray(contact['force_n'], dtype=float)
            if force.shape != (3,) or not np.isfinite(force).all():
                raise ValueError('Invalid contact force')
            row = rows[self.owners[i]]
            row['contact_area_m2'] += area
            row['normal_force_n'] += area * pressure
            row['force_n'] += force
        return {'regions': {name: {**row, 'force_n': row['force_n'].tolist(),
                                  'whole_territory_mean_pressure_pa': row['normal_force_n'] / self.areas[name]
                                  if self.areas[name] else None} for name, row in rows.items()},
                'reduction_evidence': 'engineering_force_preserving_surface_average',
                'interstitial_pressure_identified': False, 'native_commands': []}
