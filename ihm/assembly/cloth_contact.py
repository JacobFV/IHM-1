"""Approximate discrete vertex self-contact for the spring cloth.

This is a particle thickness guard, not triangle/edge collision detection: coarse
faces can intersect between vertices and fast motion can tunnel between steps.
Position projection preserves the free-pair centre of mass but can change spring
potential energy. Its correction is deliberately NOT converted into velocity.
The separate perfectly inelastic normal impulse cannot add kinetic energy.
"""
import numpy as np
from scipy.spatial import cKDTree


class ClothSelfContact:
    """Immutable topology cache, reusable across mesh checkpoints/restores.

    ``separation_m`` is the minimum centre-to-centre vertex spacing. Spring
    neighbors, and their neighbors, are excluded so local bend/shear degrees of
    freedom do not collide with their own material neighborhood.
    """

    def __init__(self, mesh, separation_m=.012, iterations=3):
        if not np.isfinite(separation_m) or separation_m <= 0:
            raise ValueError('separation_m must be finite and positive')
        if not isinstance(iterations, int) or iterations < 1:
            raise ValueError('iterations must be a positive integer')
        self.separation_m = float(separation_m)
        self.iterations = iterations
        self.count = len(mesh.x)
        neighbors = [set() for _ in range(self.count)]
        for a, b in zip(mesh.a, mesh.b):
            neighbors[int(a)].add(int(b))
            neighbors[int(b)].add(int(a))
        self.excluded = set()
        for a, immediate in enumerate(neighbors):
            local = set(immediate)
            for b in immediate:
                local.update(neighbors[b])
            self.excluded.update((min(a, b), max(a, b)) for b in local)

    def resolve(self, mesh):
        """Project overlaps and remove approaching normal relative velocity.

        Mutates ``mesh.x`` and ``mesh.v``. Free pairs conserve linear momentum;
        fixed vertices act as external supports. Returned ``support_impulse_ns``
        is the impulse imparted to those supports (opposite the mobile impulse).
        No body contact force port should receive internal free-pair impulses.
        """
        if len(mesh.x) != self.count:
            raise ValueError('Mesh topology differs from the cached topology')
        mass = np.broadcast_to(np.asarray(mesh.mass, dtype=float), (self.count,))
        if np.any(~np.isfinite(mass)) or np.any(mass <= 0):
            raise ValueError('Vertex masses must be finite and positive')
        inverse = np.where(mesh.fixed, 0., 1. / mass)
        before = float(.5 * np.sum(mass[:, None] * mesh.v ** 2))
        support = np.zeros(3)
        contacts = set()
        max_correction = 0.
        for _ in range(self.iterations):
            pairs = sorted(cKDTree(mesh.x).query_pairs(self.separation_m))
            for a, b in pairs:
                if (a, b) in self.excluded:
                    continue
                weight = inverse[a] + inverse[b]
                if weight == 0:
                    continue
                delta = mesh.x[b] - mesh.x[a]
                distance = float(np.linalg.norm(delta))
                depth = self.separation_m - distance
                if depth <= 0:
                    continue
                if distance > 1e-12:
                    normal = delta / distance
                else:
                    # Exact coincidence has no geometric normal. Use material
                    # rest separation; a deterministic axis is the final fallback.
                    delta = mesh.rest[b] - mesh.rest[a]
                    length = np.linalg.norm(delta)
                    normal = delta / length if length > 1e-12 else np.array([1., 0., 0.])
                correction = depth * normal / weight
                mesh.x[a] -= inverse[a] * correction
                mesh.x[b] += inverse[b] * correction
                max_correction = max(max_correction, depth)
                contacts.add((a, b))
                relative = float(np.dot(mesh.v[b] - mesh.v[a], normal))
                if relative < 0:
                    impulse = -relative * normal / weight
                    mesh.v[a] -= inverse[a] * impulse
                    mesh.v[b] += inverse[b] * impulse
                    if mesh.fixed[a]:
                        support -= impulse
                    if mesh.fixed[b]:
                        support += impulse
        after = float(.5 * np.sum(mass[:, None] * mesh.v ** 2))
        return {'pair_count': len(contacts), 'max_overlap_m': max_correction,
                'kinetic_energy_change_j': after - before,
                'support_impulse_ns': support.tolist()}
