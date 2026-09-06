"""Auditable regional views of native interstitial and lymph observations.

This module allocates observations, never integrates fluid or protein. Its
regional rows are an alternative view of the native owners, not extra stores.
"""
from copy import deepcopy
import math
from pathlib import Path
import json

from .body_exchange import ORGANS, digest, number
from .body_microstructure import anatomical_supports


class NativeRegionalLymphCoverage:
    def __init__(self, anatomy, organs=ORGANS, native_identity=None,
                 source_hashes=None, topology=None):
        self.organs = tuple(organs)
        if not self.organs or len(set(self.organs)) != len(self.organs) or any(o not in ORGANS for o in self.organs):
            raise ValueError('Unknown or duplicate native organ')
        ids = [e['id'] for e in anatomy['entities']]
        if len(ids) != len(set(ids)):
            raise ValueError('Duplicate anatomical entity identity')
        self.supports = anatomical_supports(deepcopy(anatomy), self.organs)
        self.native_identity = deepcopy(native_identity or {})
        self.source_hashes = dict(source_hashes or {})
        self.topology = None if topology is None else {
            'status': 'source_graph_present_without_native_owner_correspondence',
            'model_id': topology.get('model_id'),
            'nodes': len(topology.get('nodes', [])),
            'edges': len(topology.get('edges', [])),
            'source_graph': deepcopy(topology.get('source_graph', topology.get('provenance'))),
            'directed': topology.get('directed'),
        }

    @classmethod
    def from_workspace(cls, root, native_identity, organs=ORGANS):
        root = Path(root).resolve()
        if not native_identity:
            raise ValueError('Explicit native identity receipt is required')
        anatomy_path = root / 'data/derived/canonical/anatomy.json'
        anatomy = json.loads(anatomy_path.read_text())
        files = [anatomy_path, root / 'scripts/native_tissue_ports.h',
                 root / 'ihm/assembly/body_microstructure.py',
                 root / 'ihm/assembly/lymph_coverage.py']
        topology_path = root / 'data/derived/canonical/lymphatic_graph.json'
        topology = json.loads(topology_path.read_text()) if topology_path.exists() else None
        if topology is not None:
            files.append(topology_path)
            receipt = topology.get('source_graph')
            if receipt is not None:
                path, expected = receipt.get('path'), receipt.get('sha256')
                if not path or not expected:
                    raise ValueError('Missing lymph source graph receipt')
                source = (root / path).resolve()
                if not source.is_relative_to(root) or digest(source) != expected:
                    raise ValueError('Changed lymph source graph')
                files.append(source)
        hashes = {str(p.relative_to(root)): digest(p) for p in files}
        model = cls(anatomy, organs, native_identity, hashes, topology)
        for rows in model.supports.values():
            for row in rows:
                geometry = row.get('source_geometry')
                if geometry is None:
                    continue
                path, expected = geometry.get('path'), geometry.get('sha256')
                if not path or not expected:
                    raise ValueError('Missing anatomical geometry receipt')
                source = (root / path).resolve()
                if not source.is_relative_to(root) or digest(source) != expected:
                    raise ValueError('Changed anatomical support geometry')
                model.source_hashes[str(source.relative_to(root))] = expected
        return model

    def observe(self, snapshot):
        values = snapshot.get('values', {})
        unknowns = []

        def quantity(key, unit, nonnegative=False):
            raw = values.get(key)
            # Native unavailable ports can arrive as null or IEEE NaN. Neither
            # is a zero observation; reject invalid finite values and infinity.
            missing = raw is None or (isinstance(raw, float) and math.isnan(raw))
            value = None if missing else number(raw, key, nonnegative)
            if missing:
                unknowns.append(key)
            return {'source_key': key, 'unit': unit, 'value': value,
                    'status': 'unobserved' if missing else 'observed_native'}

        owners = []
        regions = []
        conservation = []
        for organ in (*self.organs, 'Lymph'):
            owner = 'Lymph' if organ == 'Lymph' else organ + '.extracellular'
            prefix = 'tissue.' + owner
            quantities = {
                'volume_ml': quantity(prefix + '.volume_ml', 'mL', True),
                'albumin_mass_g': quantity(prefix + '.Albumin.mass_g', 'g', True),
                'pressure_mmhg': quantity(prefix + '.pressure_mmhg', 'mmHg'),
                'albumin_concentration_g_per_l': quantity(prefix + '.Albumin.concentration_g_per_l', 'g/L', True),
            }
            owners.append({'native_owner': owner,
                           'native_compartment': 'Lymph' if organ == 'Lymph' else organ + 'TissueExtracellular',
                           'classification': 'native_lumped_state',
                           'accounting_owner': True, 'quantities': quantities})
            if organ == 'Lymph':
                continue
            used = {'volume_ml': [], 'albumin_mass_g': []}
            supports = self.supports[organ]
            for i, support in enumerate(supports):
                allocated = {}
                for field in used:
                    q = quantities[field]
                    total = q['value']
                    value = None if total is None else (total - math.fsum(used[field]) if i == len(supports) - 1 else total * support['fraction'])
                    if value is not None:
                        used[field].append(value)
                    allocated[field] = {**q, 'value': value,
                                        'status': 'unobserved' if value is None else 'allocated_native_observation'}
                regions.append({**deepcopy(support), 'id': support['id'] + '.interstitial',
                                'native_owner': owner, 'classification': 'regional_allocated_output' if support['material_attachment'] else 'unlocalized_native_view',
                                'accounting_owner': False, 'independent_store': False,
                                'quantities': allocated,
                                'regional_pressure_mmhg': None,
                                'regional_protein_flux_g_per_s': None})
            conservation.append({'native_owner': owner, **{
                field + '_residual': None if quantities[field]['value'] is None else math.fsum(used[field]) - quantities[field]['value']
                for field in used}})

        pathways = []
        for organ in self.organs:
            for suffix, stage in ((organ + 'E3To' + organ + 'L1', 'interstitial_drainage'),
                                  (organ + 'L1To' + organ + 'L2', 'serial_lymph_path'),
                                  (organ + 'ToLymphValve', 'lymph_return_valve')):
                pathways.append({'native_path': suffix, 'organ': organ, 'stage': stage,
                                 'flow': quantity('tissue.path.' + suffix + '.flow_ml_per_s', 'mL/s'),
                                 'regional_topology': 'unresolved', 'protein_flux_g_per_s': None})
        pathways.append({'native_path': 'LymphToVenaCava', 'organ': None,
                         'stage': 'lumped_venous_return',
                         'flow': quantity('tissue.path.LymphToVenaCava.flow_ml_per_s', 'mL/s'),
                         'regional_topology': 'unresolved', 'protein_flux_g_per_s': None})
        topology = deepcopy(self.topology) if self.topology is not None else {'status': 'absent_source_graph'}
        topology.update(native_regional_mapping='absent', regional_lymph_reservoirs='absent',
                        regional_flow_solution='absent', allocated_lymph_fluid_or_protein=False)
        return {'schema': 'native_regional_lymph_coverage_v1',
                'time_s': number(snapshot.get('time_s'), 'time', True),
                'native_identity': deepcopy(self.native_identity),
                'source_hashes': dict(self.source_hashes),
                'native_owners': owners, 'regional_views': regions,
                'pathways': pathways, 'allocation_conservation': conservation,
                'unobserved_source_keys': sorted(set(unknowns)),
                'unlocalized_native_owners': ['Lymph'] + [o + '.extracellular' for o, rows in self.supports.items() if any(r['material_attachment'] is None for r in rows)],
                'topology': topology, 'whole_body_mass_closure_claimed': False,
                'limitations': [
                    'Sum native owners OR their regional views; never add both. These views also overlap body_exchange allocations of the same owners.',
                    'Surface-area allocation is an anatomical display prior, not measured interstitial volume or drainage territory.',
                    'Lymph remains one unlocalized native owner; no graph node receives a second fluid or protein store.',
                    'Albumin is the exported protein species, not total protein. Missing observations remain unknown.',
                    'Serial path flows are distinct observations, not additive drainage; no concentration-times-flow protein flux or integrated closure is inferred.',
                ]}
