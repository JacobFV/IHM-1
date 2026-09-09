"""One-sided sampled upper envelope for a loose supine blanket.

This engineering contact geometry treats the blanket as a top cover: it cannot
wrap beneath the body or represent upright/arbitrary cloth. The upper union of
sample spheres has gaps and sharp seams; it is not a continuous anatomical skin.
"""
import numpy as np
from scipy.spatial import cKDTree


class SupineClothCover:
    """Cache the skin XY tree once per prescribed body exchange.

    Reconstruct when body points or velocities change. Query returns vertical
    penetration depth, outward upper-sphere normal, and source sample owner.
    Uncovered vertices have depth=0, normal=(0,0,1), owner=-1.
    """
    def __init__(self, points, velocity=None, radius_m=.048):
        self.points = np.array(points, float, copy=True)
        if self.points.ndim != 2 or self.points.shape[1] != 3 or not np.isfinite(self.points).all():
            raise ValueError('Skin points must be finite N by 3')
        self.velocity = np.zeros_like(self.points) if velocity is None else np.array(velocity, float, copy=True)
        if self.velocity.shape != self.points.shape or not np.isfinite(self.velocity).all():
            raise ValueError('Finite skin velocity must match points')
        if not np.isfinite(radius_m) or radius_m <= 0:
            raise ValueError('Positive finite radius required')
        self.radius = float(radius_m)
        self.tree = cKDTree(self.points[:, :2])

    def query(self, x):
        x = np.asarray(x, float)
        if x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all():
            raise ValueError('Cloth vertices must be finite N by 3')
        pairs = self.tree.sparse_distance_matrix(cKDTree(x[:, :2]), self.radius, output_type='coo_matrix')
        height = np.full(len(x), -np.inf)
        owner = np.full(len(x), -1, int)
        normal = np.zeros_like(x); normal[:, 2] = 1.
        if len(pairs.data):
            dz = np.sqrt(np.maximum(0., self.radius**2-pairs.data**2))
            tops = self.points[pairs.row, 2]+dz
            np.maximum.at(height, pairs.col, tops)
            # Deterministic lowest sample index resolves equal-height seams.
            candidates = tops == height[pairs.col]
            chosen = np.full(len(x), len(self.points), int)
            np.minimum.at(chosen, pairs.col[candidates], pairs.row[candidates])
            covered = chosen < len(self.points)
            owner[covered] = chosen[covered]
            delta_xy = x[covered, :2]-self.points[owner[covered], :2]
            delta_z = height[covered]-self.points[owner[covered], 2]
            normal[covered, :2] = delta_xy/self.radius
            normal[covered, 2] = delta_z/self.radius
            # At the hemisphere rim the normal is horizontal and the height
            # derivative is singular; its analytic unit normal remains finite.
            normal[covered] /= np.maximum(np.linalg.norm(normal[covered], axis=1)[:, None], 1e-15)
        return np.maximum(0., height-x[:, 2]), normal, owner

    def project(self, mesh, h, friction=.45):
        """Project onto the envelope along its normal and remove approach only.

        ``query`` returns a VERTICAL penetration depth against a height field.
        Repairing that depth vertically displaces the vertex by 1/normal_z times
        the true geometric overlap, and normal_z goes to zero at every sample
        sphere's rim, so a vertex near the body silhouette used to be thrown
        metres of spring extension away from its neighbours. The minimum-norm
        repair of the constraint z >= height(x, y) is along the analytic surface
        normal, by depth*normal_z; that is the correction applied here. Curvature
        can leave a small residual, which the caller's alternating passes absorb.

        Normal and Coulomb tangential impulses use an uncalibrated engineering
        friction coefficient. Position correction does not become velocity.
        Boundary work is estimated
        as impulse dotted with prescribed sample velocity; native body work can
        differ under explicit coupling. Gravity/spring projection work is separate.
        """
        if not np.isfinite(h) or h <= 0:
            raise ValueError('Positive finite substep required')
        if not np.isfinite(friction) or friction < 0:
            raise ValueError('Finite nonnegative friction required')
        mass = np.broadcast_to(np.asarray(mesh.mass, float), (len(mesh.x),))
        if not np.isfinite(mass).all() or np.any(mass <= 0):
            raise ValueError('Finite positive vertex masses required')
        depth, normal, owner = self.query(mesh.x)
        depth[mesh.fixed] = 0.
        active = (depth > 0) & (owner >= 0)
        impulses = np.zeros_like(self.points)
        before_kinetic = float(.5*np.sum(mass[:, None]*mesh.v**2))

        def spring_energy():
            if not all(hasattr(mesh, key) for key in ('a','b','length','stiffness')):
                return None
            lengths = np.linalg.norm(mesh.x[mesh.b]-mesh.x[mesh.a], axis=1)
            return float(.5*np.sum(mesh.stiffness*(lengths-mesh.length)**2))

        before_spring = spring_energy()
        repair = (depth*normal[:, 2])[:, None]*normal
        vertical = repair[:, 2].copy()
        mesh.x += repair
        boundary_work = 0.
        if np.any(active):
            skin_velocity = self.velocity[owner[active]]
            speed = np.sum((mesh.v[active]-skin_velocity)*normal[active], axis=1)
            normal_dv = np.maximum(0., -speed)
            dv = normal_dv[:, None]*normal[active]
            relative = mesh.v[active]+dv-skin_velocity
            tangent = relative-np.sum(relative*normal[active], axis=1)[:, None]*normal[active]
            tangent_speed = np.linalg.norm(tangent, axis=1)
            # Coulomb impulse uses this actual normal impulse as its budget.
            # No invented resting normal load or persistent friction is added.
            reduction = np.minimum(tangent_speed, friction*normal_dv)
            dv -= reduction[:, None]*tangent/np.maximum(tangent_speed[:, None], 1e-15)
            impulse = mass[active, None]*dv
            mesh.v[active] += dv
            np.add.at(impulses, owner[active], -impulse)
            boundary_work = float(np.sum(impulse*skin_velocity))
        kinetic_change = float(.5*np.sum(mass[:, None]*mesh.v**2))-before_kinetic
        after_spring = spring_energy()
        mesh.skin_barrier_state = {
            'friction_coefficient': float(friction),
            'vertex_count': int(np.count_nonzero(active)),
            'covered_vertex_count': int(np.count_nonzero(owner >= 0)),
            'max_projection_m': float(np.linalg.norm(repair, axis=1).max(initial=0.)),
            'sum_projection_m': float(np.linalg.norm(repair, axis=1).sum()),
            'max_vertical_penetration_m': float(depth.max(initial=0.)),
            'spring_energy_change_j': None if before_spring is None else after_spring-before_spring,
            'gravity_projection_energy_change_j': float(np.sum(mass*vertical)*9.81),
            'kinetic_energy_change_j': kinetic_change,
            'prescribed_boundary_work_j': boundary_work,
            'contact_impulse_dissipation_j': boundary_work-kinetic_change,
            'scope': 'Supine loose top cover over upper sample-sphere envelope only; no underside wrapping or upright cloth. Minimum-norm normal-direction repair; its work is nonconservative; gravity audit assumes canonical -Z gravity at 9.81 m/s2. Prescribed boundary work is not measured native work.'}
        return impulses
