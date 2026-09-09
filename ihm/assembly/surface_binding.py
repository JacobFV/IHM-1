"""Distributed material attachment matching EnvironmentDynamics skin contacts.

Hard assignment to native segment envelopes is an explicit geometric prior, not
skin FEM or joint-continuity reconstruction. Coordinates remain canonical here;
the shared world/view rigid transform must be applied only after attachment.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import numpy as np

SKIN_ASSET = 'data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz'


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


class SegmentSurfaceBinding:
    """Freeze registration support, then emit only 22 segment poses per frame."""
    def __init__(self, registration, *, source_identity=None,
                 surface_entity_ids=('body-bp3d-FJ2810',)):
        segments = []
        for body in sorted(registration.groups):
            group = registration.groups[body]
            bone = group['canonical_bones'][0]
            segments.append(dict(id=body, bone_id=bone,
                bounds_min_m=np.asarray(group['bounds_min_m'], float).tolist(),
                bounds_max_m=np.asarray(group['bounds_max_m'], float).tolist(),
                reference_centroid_m=np.asarray(registration.specs[bone]['centroid_m'], float).tolist()))
        if not segments:
            raise ValueError('Surface attachment needs retained native segment supports')
        manifest = dict(schema='ihm.segment-surface-binding.v1',
            coordinate_frame='canonical_current_world',
            reference_coordinate_frame='canonical_rest',
            rule='nearest_named_bone_envelope',
            tie_break='lexicographic_native_segment_id',
            segments=segments, surface_entity_ids=list(surface_entity_ids),
            source_identity=deepcopy(source_identity or {}),
            equation='current = bone_centroid + bone_rotation @ (rest_vertex - reference_bone_centroid)',
            scope='Hard per-rest-vertex native segment attachment, exactly matching sampled world skin support. No FEM, skin continuity, respiration displacement, or added inertial owner. Independent cloth positions override this attachment.')
        manifest['binding_identity'] = _digest(manifest)
        self._manifest = manifest
        self._segments = segments

    @classmethod
    def from_root(cls, root, registration):
        root = Path(root)
        sources = {p: hashlib.sha256((root / p).read_bytes()).hexdigest()
                   for p in ('data/derived/canonical/mechanics.json', SKIN_ASSET)}
        sources['registration_manifest_sha256'] = _digest(registration.manifest())
        sources['binding_implementation_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        return cls(registration, source_identity=sources)

    def manifest(self):
        return deepcopy(self._manifest)

    def frame(self, entities):
        result = {}
        for segment in self._segments:
            entity = entities[segment['bone_id']]
            centroid = np.asarray(entity['centroid_m'], float)
            rotation = np.asarray(entity['rotation_matrix'], float)
            if centroid.shape != (3,) or rotation.shape != (3, 3) or not np.isfinite(centroid).all() or not np.isfinite(rotation).all():
                raise ValueError('Surface segment pose must be finite canonical centroid and rotation')
            result[segment['id']] = dict(centroid_m=centroid.tolist(), rotation_matrix=rotation.tolist())
        return result

    def bind(self, rest_vertices):
        points = np.asarray(rest_vertices, float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError('Surface vertices must be finite N by 3 canonical rest coordinates')
        # Sorted segment IDs make argmin's first tie identical to _ranking.
        distances = np.stack([np.linalg.norm(np.maximum(np.maximum(
            np.asarray(s['bounds_min_m']) - points,
            points - np.asarray(s['bounds_max_m'])), 0), axis=1)
            for s in self._segments], axis=1)
        owners = np.argmin(distances, axis=1)
        centers = np.asarray([s['reference_centroid_m'] for s in self._segments])
        return owners, points - centers[owners]

    def project(self, binding, transforms):
        owners, offsets = binding
        result = np.empty_like(offsets)
        for index, segment in enumerate(self._segments):
            selected = owners == index
            pose = transforms[segment['id']]
            # Match the contact implementation's matrix-vector operation exactly.
            result[selected] = np.asarray(pose['centroid_m']) + np.asarray([
                np.asarray(pose['rotation_matrix']) @ offset for offset in offsets[selected]]).reshape(-1, 3)
        return result
