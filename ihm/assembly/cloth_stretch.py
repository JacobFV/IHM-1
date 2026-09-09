"""Engineering unilateral yarn-extension constraints for the spring cloth.

Material edges can compress/fold freely; stretch/shear edges have a declared
extension ceiling. This is discrete strain limiting, not calibrated fabric or
an implicit continuum solve. Position projection changes potential energy and
is separately audited; it is never converted into a velocity kick.
"""
import numpy as np


class ClothStretchConstraint:
    def __init__(self, mesh, max_extension_ratio=1.12, iterations=24, tolerance=1e-4):
        if not np.isfinite(max_extension_ratio) or max_extension_ratio < 1:
            raise ValueError('Extension ratio must be finite and at least one')
        if type(iterations) is not int or iterations < 1:
            raise ValueError('Positive integer iterations required')
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('Positive finite tolerance required')
        # Two-edge bend springs are not material yarn lengths.
        keep = np.asarray(mesh.edges)[:, 2] >= .99
        self.a = np.asarray(mesh.a, int)[keep].copy()
        self.b = np.asarray(mesh.b, int)[keep].copy()
        self.rest = np.asarray(mesh.length, float)[keep].copy()
        if np.any(self.rest <= 0) or not np.isfinite(self.rest).all():
            raise ValueError('Finite positive material rest lengths required')
        self.ratio = float(max_extension_ratio)
        self.limit = self.ratio*self.rest
        self.iterations = iterations
        self.tolerance = float(tolerance)
        self.count = len(mesh.x)
        # Disjoint edges in each color can be projected simultaneously without
        # conflicting scatter updates. Colors remain deterministic on replay.
        occupied = [set() for _ in range(self.count)]
        colors = []
        for edge, (a, b) in enumerate(zip(self.a, self.b)):
            color = 0
            while color in occupied[a] or color in occupied[b]:
                color += 1
            while color >= len(colors):
                colors.append([])
            colors[color].append(edge)
            occupied[a].add(color); occupied[b].add(color)
        self.colors = [np.array(c, int) for c in colors]

    def resolve(self, mesh, gravity=None):
        """Mutate positions/velocities, returning honest convergence/work data.

        Alternate this with body/static/self-contact projection; a final stretch
        solve alone can pull cloth through colliders. Incompatible anchors or
        collision constraints require changed physical boundaries, not new rest
        lengths. ``converged`` reports only yarn constraints at this instant.
        """
        if len(mesh.x) != self.count:
            raise ValueError('Cached topology does not match mesh')
        mass = np.broadcast_to(np.asarray(mesh.mass, float), (self.count,))
        if not np.isfinite(mass).all() or np.any(mass <= 0):
            raise ValueError('Finite positive vertex masses required')
        inverse = np.where(mesh.fixed, 0., 1./mass)
        gravity = np.zeros(3) if gravity is None else np.asarray(gravity, float)
        if gravity.shape != (3,) or not np.isfinite(gravity).all():
            raise ValueError('Gravity requires three finite components')
        x0 = mesh.x.copy()
        kinetic0 = float(.5*np.sum(mass[:, None]*mesh.v**2))

        def spring_energy():
            lengths = np.linalg.norm(mesh.x[mesh.b]-mesh.x[mesh.a], axis=1)
            return float(.5*np.sum(mesh.stiffness*(lengths-mesh.length)**2))

        potential0 = spring_energy()
        support = np.zeros(3)
        projected = set()
        for sweep in range(self.iterations):
            for color in self.colors:
                a, b = self.a[color], self.b[color]
                delta = mesh.x[b]-mesh.x[a]
                length = np.linalg.norm(delta, axis=1)
                normal = delta/np.maximum(length[:, None], 1e-12)
                weights = inverse[a]+inverse[b]
                depth = np.maximum(0., length-self.limit[color])
                active = (depth > 0) & (weights > 0)
                if not np.any(active):
                    continue
                projected.update(color[active].tolist())
                correction = np.where(active, depth/np.maximum(weights, 1e-30), 0.)[:, None]*normal
                mesh.x[a] += inverse[a, None]*correction
                mesh.x[b] -= inverse[b, None]*correction
                speed = np.sum((mesh.v[b]-mesh.v[a])*normal, axis=1)
                magnitude = np.where(active, np.maximum(0., speed)/np.maximum(weights, 1e-30), 0.)
                impulse = magnitude[:, None]*normal
                mesh.v[a] += inverse[a, None]*impulse
                mesh.v[b] -= inverse[b, None]*impulse
                support += impulse[mesh.fixed[a]].sum(axis=0)
                support -= impulse[mesh.fixed[b]].sum(axis=0)
            ratios = np.linalg.norm(mesh.x[self.b]-mesh.x[self.a], axis=1)/self.rest
            if np.all(ratios <= self.ratio+self.tolerance):
                break
        ratios = np.linalg.norm(mesh.x[self.b]-mesh.x[self.a], axis=1)/self.rest
        displacement = mesh.x-x0
        result = {
            'max_extension_ratio': float(ratios.max(initial=0.)),
            'extension_limit_ratio': self.ratio,
            'converged': bool(np.all(ratios <= self.ratio+self.tolerance)),
            'violating_edges': int(np.count_nonzero(ratios > self.ratio+self.tolerance)),
            'fixed_fixed_violations': int(np.count_nonzero((ratios > self.ratio+self.tolerance) & mesh.fixed[self.a] & mesh.fixed[self.b])),
            'projected_edges': len(projected), 'iterations': sweep+1,
            'max_projection_m': float(np.linalg.norm(displacement, axis=1).max(initial=0.)),
            'spring_energy_change_j': spring_energy()-potential0,
            'gravity_projection_energy_change_j': float(-np.sum(mass[:, None]*displacement*gravity)),
            'kinetic_energy_change_j': float(.5*np.sum(mass[:, None]*mesh.v**2))-kinetic0,
            'support_impulse_ns': support.tolist(),
            'scope': 'Uncalibrated unilateral stretch/shear edge limit; no area constraint, continuum fabric, or collider compatibility guarantee. Projection work is nonconservative and separately reported.'}
        return result
