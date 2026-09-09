"""Discrete hinge bending energy for the triangulated cloth sheet.

The existing sheet resists only edge extension. A crease that keeps every edge
length is therefore a zero-energy deformation, so folds cost nothing and the
sharp creases observed in the blanket are an unpenalised mesh mode rather than
a rendering artifact. This module supplies the missing shell term.

The energy is the quadratic (small-angle) bending model of Bergou, Wardetzky,
Harmon, Zorin and Grinspun, "A Quadratic Bending Model for Inextensible
Surfaces" (SGP 2006), whose per-hinge form agrees with the discrete-shells
hinge energy 3*k*|edge|^2*theta^2/(A0+A1) to leading order in the dihedral
angle. It is exact only for a FLAT rest state and only isometric deformation;
at large fold angles it understates the restoring torque, and the rigidity is
an engineering value, not a measured Kawabata textile bending modulus.

Its Hessian is constant and positive semidefinite, so it can be handed to an
implicit solve unchanged and cannot itself create energy.
"""
import numpy as np
from scipy import sparse

SCOPE = ('Quadratic flat-rest hinge bending (Bergou et al. 2006); constant PSD Hessian; '
         'exact for isometric deformation of a planar rest sheet, understates large-fold torque; '
         'uncalibrated engineering rigidity, not a measured textile bending modulus.')


def hinges_from_triangles(triangles):
    """Interior edges shared by exactly two triangles, as (p, q, r, s) flaps.

    ``p, q`` are the hinge edge; ``r, s`` are the opposite vertices. Edges used
    by one triangle are boundaries and edges used by more than two make the
    surface nonmanifold: both are excluded rather than silently averaged.
    """
    triangles = np.asarray(triangles, int)
    if triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError('Triangles require an M by 3 index array')
    shared = {}
    for face in triangles:
        for i in range(3):
            p, q = int(face[i]), int(face[(i+1) % 3])
            if p == q:
                raise ValueError('Degenerate triangle with a repeated vertex')
            shared.setdefault((min(p, q), max(p, q)), []).append(int(face[(i+2) % 3]))
    flaps = []
    for (p, q), opposite in sorted(shared.items()):
        if len(opposite) == 2:
            flaps.append((p, q, opposite[0], opposite[1]))
    return np.asarray(flaps, int).reshape(-1, 4)


class ClothBending:
    """Immutable hinge stencils built once from a planar rest configuration.

    ``rigidity_n_m`` is the bending rigidity B in N*m. The fabric bending
    length (B / (areal weight))**(1/3) is reported so the caller can see how
    the choice compares with the mesh spacing; a bending length below the
    spacing means the mesh cannot resolve its own drape.
    """

    def __init__(self, rest_positions, triangles, rigidity_n_m=1e-3, planarity_tolerance_m=1e-9):
        rest = np.asarray(rest_positions, float)
        if rest.ndim != 2 or rest.shape[1] != 3 or not np.isfinite(rest).all():
            raise ValueError('Rest positions require finite N by 3 values')
        if not np.isfinite(rigidity_n_m) or rigidity_n_m < 0:
            raise ValueError('Bending rigidity must be finite and nonnegative')
        if not np.isfinite(planarity_tolerance_m) or planarity_tolerance_m <= 0:
            raise ValueError('Positive finite planarity tolerance required')
        self.count = len(rest)
        self.rigidity = float(rigidity_n_m)
        self.flaps = hinges_from_triangles(triangles)
        if len(self.flaps) and self.flaps.max() >= self.count:
            raise ValueError('Triangle indices exceed the vertex count')
        p, q, r, s = (rest[self.flaps[:, i]] for i in range(4)) if len(self.flaps) else (np.zeros((0, 3)),)*4
        edge = q-p
        length = np.linalg.norm(edge, axis=1)
        if np.any(length <= 0):
            raise ValueError('Hinge edges require positive rest length')
        unit = edge/length[:, None]
        # Signed positions of the opposite vertices in the hinge frame.
        along_r = np.sum((r-p)*unit, axis=1)
        along_s = np.sum((s-p)*unit, axis=1)
        perp_r = (r-p)-along_r[:, None]*unit
        perp_s = (s-p)-along_s[:, None]*unit
        height_r = np.linalg.norm(perp_r, axis=1)
        height_s = np.linalg.norm(perp_s, axis=1)
        if np.any(height_r <= 0) or np.any(height_s <= 0):
            raise ValueError('Degenerate zero-area rest triangle in a hinge flap')
        offset = np.abs(np.sum(np.cross(perp_r, perp_s), axis=1))/np.maximum(height_r*height_s, 1e-300)
        if np.any(offset > planarity_tolerance_m):
            raise ValueError('This quadratic bending model requires a planar rest flap')
        # K is the affine dependency of the four coplanar rest points, scaled so
        # that |K . x| equals |edge| * dihedral angle for a small fold. Then
        # c * |K . x|^2 reproduces the discrete-shells hinge energy exactly to
        # leading order, with the cotangent weights implied rather than fitted.
        weight_r = -length/height_r
        weight_s = -length/height_s
        fraction_r = along_r/length
        fraction_s = along_s/length
        weight_q = -weight_r*fraction_r-weight_s*fraction_s
        weight_p = -(weight_q+weight_r+weight_s)
        self.stencil = np.column_stack([weight_p, weight_q, weight_r, weight_s])
        area = .5*length*(height_r+height_s)
        self.coefficient = 3.*self.rigidity/np.maximum(area, 1e-300)
        self.rest_length = length
        residual = np.einsum('ei,eij->ej', self.stencil, np.stack([p, q, r, s], axis=1))
        self.rest_stencil_residual_m = float(np.linalg.norm(residual, axis=1).max(initial=0.))
        if self.rest_stencil_residual_m > 1e-6*max(1., float(length.max(initial=1.))):
            raise ValueError('Hinge stencil does not annihilate its own rest flap')

    def _combination(self, x):
        return np.einsum('ei,eij->ej', self.stencil, x[self.flaps])

    def energy(self, positions):
        """Bending energy in joules. Exactly zero for any affine image of rest."""
        x = np.asarray(positions, float)
        if x.shape != (self.count, 3):
            raise ValueError('Positions must match the cached vertex count')
        if not len(self.flaps):
            return 0.
        return float(np.sum(self.coefficient*np.sum(self._combination(x)**2, axis=1)))

    def forces(self, positions):
        """Conservative forces, minus the energy gradient. Sum to zero exactly."""
        x = np.asarray(positions, float)
        if x.shape != (self.count, 3):
            raise ValueError('Positions must match the cached vertex count')
        force = np.zeros_like(x)
        if not len(self.flaps):
            return force
        combination = self._combination(x)
        contribution = -2.*self.coefficient[:, None, None]*self.stencil[:, :, None]*combination[:, None, :]
        for corner in range(4):
            np.add.at(force, self.flaps[:, corner], contribution[:, corner])
        return force

    def hessian(self):
        """Constant symmetric positive semidefinite 3N by 3N energy Hessian."""
        size = 3*self.count
        if not len(self.flaps):
            return sparse.csr_matrix((size, size))
        rows, cols, data = [], [], []
        axis = np.arange(3)
        for i in range(4):
            for j in range(4):
                block = 2.*self.coefficient*self.stencil[:, i]*self.stencil[:, j]
                rows.append(np.repeat(3*self.flaps[:, i][:, None]+axis, 1, axis=0).ravel())
                cols.append((3*self.flaps[:, j][:, None]+axis).ravel())
                data.append(np.repeat(block, 3))
        return sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                                 shape=(size, size)).tocsr()

    def dihedral_angles_deg(self, positions):
        """Actual hinge fold angles, for geometric acceptance rather than looks."""
        x = np.asarray(positions, float)
        if x.shape != (self.count, 3):
            raise ValueError('Positions must match the cached vertex count')
        if not len(self.flaps):
            return np.zeros(0)
        p, q, r, s = (x[self.flaps[:, i]] for i in range(4))
        edge = q-p
        length = np.maximum(np.linalg.norm(edge, axis=1), 1e-300)
        unit = edge/length[:, None]
        left = np.cross(unit, r-p)
        right = np.cross(s-p, unit)
        left /= np.maximum(np.linalg.norm(left, axis=1), 1e-300)[:, None]
        right /= np.maximum(np.linalg.norm(right, axis=1), 1e-300)[:, None]
        cosine = np.clip(np.sum(left*right, axis=1), -1., 1.)
        return np.degrees(np.arccos(cosine))

    def state(self, positions, areal_density_kg_m2=None):
        angles = self.dihedral_angles_deg(positions)
        record = {'hinge_count': int(len(self.flaps)), 'rigidity_n_m': self.rigidity,
                  'bending_energy_j': self.energy(positions),
                  'max_fold_angle_deg': float(angles.max(initial=0.)),
                  'p99_fold_angle_deg': float(np.percentile(angles, 99)) if len(angles) else 0.,
                  'mean_fold_angle_deg': float(angles.mean()) if len(angles) else 0.,
                  'folds_over_60_deg': int(np.count_nonzero(angles > 60.)),
                  'rest_stencil_residual_m': self.rest_stencil_residual_m,
                  'mean_hinge_edge_m': float(self.rest_length.mean()) if len(self.flaps) else 0.,
                  'scope': SCOPE}
        if areal_density_kg_m2:
            record['bending_length_m'] = float((self.rigidity/(areal_density_kg_m2*9.81))**(1/3))
        return record
