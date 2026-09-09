"""Triangle-level and continuous self-collision for the cloth sheet.

``cloth_contact.ClothSelfContact`` guards only centre-to-centre vertex spacing.
Two coarse faces can therefore pass through each other between their vertices,
and a fast vertex can cross a face inside one step without ever being close to
one of its corners. This module adds the two primitive pairs that a cloth sheet
actually needs: vertex against triangle, and edge against edge, both as a
discrete proximity test and as a continuous swept test over the step.

Momentum properties are exact rather than approximate. Every response acts
along the closest-point separation direction with barycentric weights that sum
to zero, so the impulse conserves linear momentum, and because the weighted
material points are separated only along that same direction, it conserves
angular momentum as well. The normal impulse is perfectly inelastic and can
only remove kinetic energy. Position repair is applied separately and is never
converted into velocity, matching the rest of this cloth stack.

Limits: proximity is enforced at a declared sheet thickness, not a measured
fabric one; there is no friction here; the continuous test assumes vertices
move on straight lines within the step; and repeated projection is a Gauss
Seidel sweep with an iteration cap, so a dense pile can end a step with a
residual overlap. The residual is reported rather than hidden.
"""
import numpy as np
from scipy.spatial import cKDTree

SCOPE = ('Vertex/triangle and edge/edge proximity plus straight-line continuous self collision; '
         'exact linear and angular momentum, perfectly inelastic normal impulse, position repair '
         'never converted to velocity; declared sheet thickness, no friction, capped Gauss-Seidel sweeps.')


def _closest_point_triangle(point, a, b, c):
    """Vectorized closest point on a triangle with barycentric weights."""
    ab = b-a; ac = c-a; ap = point-a
    d1 = np.sum(ab*ap, axis=1); d2 = np.sum(ac*ap, axis=1)
    bp = point-b; d3 = np.sum(ab*bp, axis=1); d4 = np.sum(ac*bp, axis=1)
    cp = point-c; d5 = np.sum(ab*cp, axis=1); d6 = np.sum(ac*cp, axis=1)
    vc = d1*d4-d3*d2; vb = d5*d2-d1*d6; va = d3*d6-d5*d4
    denominator = va+vb+vc
    u = np.zeros(len(point)); v = np.zeros(len(point)); w = np.zeros(len(point))
    interior = (va > 0) & (vb > 0) & (vc > 0) & (denominator > 0)
    with np.errstate(invalid='ignore', divide='ignore'):
        inverse = np.where(interior, 1./np.where(denominator == 0, 1., denominator), 0.)
    u[interior] = (va*inverse)[interior]; v[interior] = (vb*inverse)[interior]; w[interior] = (vc*inverse)[interior]
    corner_a = (d1 <= 0) & (d2 <= 0) & ~interior
    corner_b = (d3 >= 0) & (d4 <= d3) & ~interior & ~corner_a
    corner_c = (d6 >= 0) & (d5 <= d6) & ~interior & ~corner_a & ~corner_b
    u[corner_a] = 1.; v[corner_b] = 1.; w[corner_c] = 1.
    done = interior | corner_a | corner_b | corner_c
    edge_ab = (vc <= 0) & (d1 >= 0) & (d3 <= 0) & ~done
    edge_ac = (vb <= 0) & (d2 >= 0) & (d6 <= 0) & ~done & ~edge_ab
    edge_bc = ~done & ~edge_ab & ~edge_ac
    def parameter(numerator, denominator_):
        return np.clip(np.where(np.abs(denominator_) < 1e-300, 0., numerator/np.where(denominator_ == 0, 1., denominator_)), 0., 1.)
    t = parameter(d1, d1-d3); u[edge_ab] = (1-t)[edge_ab]; v[edge_ab] = t[edge_ab]
    t = parameter(d2, d2-d6); u[edge_ac] = (1-t)[edge_ac]; w[edge_ac] = t[edge_ac]
    t = parameter(d4-d3, (d4-d3)+(d5-d6)); v[edge_bc] = (1-t)[edge_bc]; w[edge_bc] = t[edge_bc]
    return np.column_stack([u, v, w])


def _closest_segment_segment(p0, p1, q0, q1):
    """Vectorized clamped closest-point parameters on two segments."""
    d1 = p1-p0; d2 = q1-q0; r = p0-q0
    a = np.sum(d1*d1, axis=1); e = np.sum(d2*d2, axis=1); f = np.sum(d2*r, axis=1)
    c = np.sum(d1*r, axis=1); b = np.sum(d1*d2, axis=1)
    denominator = a*e-b*b
    safe = np.where(np.abs(denominator) < 1e-300, 1., denominator)
    s = np.where(np.abs(denominator) < 1e-300, 0., np.clip((b*f-c*e)/safe, 0., 1.))
    t = (b*s+f)/np.where(e < 1e-300, 1., e)
    t = np.clip(t, 0., 1.)
    s = np.clip(np.where(a < 1e-300, 0., (b*t-c)/np.where(a < 1e-300, 1., a)), 0., 1.)
    return s, t


def _cubic_roots(coefficients):
    """Real roots in [0, 1] of a possibly degenerate cubic, ascending."""
    roots = []
    trimmed = np.trim_zeros(np.asarray(coefficients, float), 'f')
    if not len(trimmed) or not np.isfinite(trimmed).all():
        return roots
    if len(trimmed) == 1:
        return roots
    for root in np.roots(trimmed):
        if abs(root.imag) <= 1e-9*max(1., abs(root.real)) and -1e-9 <= root.real <= 1.+1e-9:
            roots.append(min(1., max(0., float(root.real))))
    return sorted(roots)


class ClothSelfCollision:
    """Immutable primitive topology built once from the triangle list.

    ``thickness_m`` is the declared sheet thickness enforced between
    non-adjacent primitives. ``iterations`` caps the Gauss-Seidel sweeps.
    """

    def __init__(self, triangles, vertex_count, thickness_m=.006, iterations=4):
        triangles = np.asarray(triangles, int)
        if triangles.ndim != 2 or triangles.shape[1] != 3:
            raise ValueError('Triangles require an M by 3 index array')
        if type(vertex_count) is not int or vertex_count < 1:
            raise ValueError('Positive integer vertex count required')
        if triangles.size and (triangles.min() < 0 or triangles.max() >= vertex_count):
            raise ValueError('Triangle indices outside the vertex range')
        if not np.isfinite(thickness_m) or thickness_m <= 0:
            raise ValueError('Positive finite thickness required')
        if type(iterations) is not int or iterations < 1:
            raise ValueError('Positive integer iteration count required')
        self.triangles = triangles.copy()
        self.count = vertex_count
        self.thickness = float(thickness_m)
        self.iterations = iterations
        edges = set()
        for face in triangles:
            for i in range(3):
                p, q = int(face[i]), int(face[(i+1) % 3])
                edges.add((min(p, q), max(p, q)))
        self.edges = np.asarray(sorted(edges), int).reshape(-1, 2)

    def _pairs(self, start, end):
        """Broad phase over swept primitive bounds; returns candidate pairs."""
        margin = self.thickness
        vertex_low = np.minimum(start, end)-margin
        vertex_high = np.maximum(start, end)+margin
        face_low = np.minimum(start[self.triangles], end[self.triangles]).min(axis=1)-margin
        face_high = np.maximum(start[self.triangles], end[self.triangles]).max(axis=1)+margin
        edge_low = np.minimum(start[self.edges], end[self.edges]).min(axis=1)-margin
        edge_high = np.maximum(start[self.edges], end[self.edges]).max(axis=1)+margin

        def overlaps(low_a, high_a, low_b, high_b, radius_scale=1.):
            centre_a = .5*(low_a+high_a); centre_b = .5*(low_b+high_b)
            reach = .5*np.max(np.linalg.norm(high_a-low_a, axis=1), initial=0.)
            reach += .5*np.max(np.linalg.norm(high_b-low_b, axis=1), initial=0.)
            if not len(centre_a) or not len(centre_b):
                return np.zeros((0, 2), int)
            found = cKDTree(centre_a).query_ball_tree(cKDTree(centre_b), reach*radius_scale+1e-12)
            rows = np.concatenate([np.full(len(f), i, int) for i, f in enumerate(found)]) if any(found) else np.zeros(0, int)
            cols = np.concatenate([np.asarray(f, int) for f in found if f]) if any(found) else np.zeros(0, int)
            if not len(rows):
                return np.zeros((0, 2), int)
            keep = np.all((low_a[rows] <= high_b[cols]) & (high_a[rows] >= low_b[cols]), axis=1)
            return np.column_stack([rows[keep], cols[keep]])

        vertex_face = overlaps(vertex_low, vertex_high, face_low, face_high)
        if len(vertex_face):
            faces = self.triangles[vertex_face[:, 1]]
            vertex_face = vertex_face[~np.any(faces == vertex_face[:, 0][:, None], axis=1)]
        edge_edge = overlaps(edge_low, edge_high, edge_low, edge_high)
        if len(edge_edge):
            first = self.edges[edge_edge[:, 0]]; second = self.edges[edge_edge[:, 1]]
            distinct = edge_edge[:, 0] < edge_edge[:, 1]
            shared = np.any(first[:, :, None] == second[:, None, :], axis=(1, 2))
            edge_edge = edge_edge[distinct & ~shared]
        return vertex_face, edge_edge

    def _stencils(self, x, vertex_face, edge_edge):
        """Contact stencils: indices, signed weights, normal and current gap."""
        indices = []; weights = []; normals = []; gaps = []
        if len(vertex_face):
            faces = self.triangles[vertex_face[:, 1]]
            point = x[vertex_face[:, 0]]
            bary = _closest_point_triangle(point, x[faces[:, 0]], x[faces[:, 1]], x[faces[:, 2]])
            closest = np.einsum('ij,ijk->ik', bary, x[faces])
            delta = point-closest
            distance = np.linalg.norm(delta, axis=1)
            indices.append(np.column_stack([vertex_face[:, 0], faces]))
            weights.append(np.column_stack([np.ones(len(faces)), -bary]))
            normals.append(delta/np.maximum(distance[:, None], 1e-300))
            gaps.append(distance)
        if len(edge_edge):
            first = self.edges[edge_edge[:, 0]]; second = self.edges[edge_edge[:, 1]]
            s, t = _closest_segment_segment(x[first[:, 0]], x[first[:, 1]], x[second[:, 0]], x[second[:, 1]])
            point_a = x[first[:, 0]]+s[:, None]*(x[first[:, 1]]-x[first[:, 0]])
            point_b = x[second[:, 0]]+t[:, None]*(x[second[:, 1]]-x[second[:, 0]])
            delta = point_a-point_b
            distance = np.linalg.norm(delta, axis=1)
            indices.append(np.column_stack([first, second]))
            weights.append(np.column_stack([1-s, s, -(1-t), -t]))
            normals.append(delta/np.maximum(distance[:, None], 1e-300))
            gaps.append(distance)
        if not indices:
            return (np.zeros((0, 4), int), np.zeros((0, 4)), np.zeros((0, 3)), np.zeros(0),
                    np.zeros(0, bool))
        kind = np.concatenate([np.full(len(i), n == 0, bool) for n, i in enumerate(indices)])
        return (np.concatenate(indices), np.concatenate(weights), np.concatenate(normals),
                np.concatenate(gaps), kind)

    def _advancing(self, start, travel, vertex_face, edge_edge):
        """Drop candidate pairs that provably cannot reach contact this step."""
        if len(vertex_face):
            faces = self.triangles[vertex_face[:, 1]]
            bary = _closest_point_triangle(start[vertex_face[:, 0]], start[faces[:, 0]],
                                           start[faces[:, 1]], start[faces[:, 2]])
            gap = np.linalg.norm(start[vertex_face[:, 0]]-np.einsum('ij,ijk->ik', bary, start[faces]), axis=1)
            reach = travel[vertex_face[:, 0]]+travel[faces].max(axis=1)
            vertex_face = vertex_face[gap-reach <= self.thickness]
        if len(edge_edge):
            first = self.edges[edge_edge[:, 0]]; second = self.edges[edge_edge[:, 1]]
            s_, t_ = _closest_segment_segment(start[first[:, 0]], start[first[:, 1]],
                                              start[second[:, 0]], start[second[:, 1]])
            point_a = start[first[:, 0]]+s_[:, None]*(start[first[:, 1]]-start[first[:, 0]])
            point_b = start[second[:, 0]]+t_[:, None]*(start[second[:, 1]]-start[second[:, 0]])
            gap = np.linalg.norm(point_a-point_b, axis=1)
            reach = travel[first].max(axis=1)+travel[second].max(axis=1)
            edge_edge = edge_edge[gap-reach <= self.thickness]
        return vertex_face, edge_edge

    def _continuous(self, start, end, vertex_face, edge_edge):
        """Straight-line swept coplanarity roots inside the step, if any.

        A conservative advancement filter runs first: two primitives whose start
        separation exceeds the thickness plus the total distance their vertices
        travel cannot touch during the step, whatever their paths. That is an
        exact exclusion, not a heuristic cutoff, and it leaves only the pairs
        that actually need a cubic root.
        """
        hits = []
        velocity = end-start
        travel = np.linalg.norm(velocity, axis=1)
        vertex_face, edge_edge = self._advancing(start, travel, vertex_face, edge_edge)
        for pair, is_vertex_face in ((vertex_face, True), (edge_edge, False)):
            for row in pair:
                if is_vertex_face:
                    face = self.triangles[row[1]]
                    nodes = np.array([row[0], face[0], face[1], face[2]], int)
                else:
                    nodes = np.concatenate([self.edges[row[0]], self.edges[row[1]]])
                p = start[nodes]; d = velocity[nodes]
                if is_vertex_face:
                    order = (1, 2, 3, 0)
                else:
                    order = (0, 1, 2, 3)
                base, first, second, probe = order
                # ((b - a) x (c - a)) . (p - a) is a cubic in the step fraction.
                a0 = p[first]-p[base]; a1 = d[first]-d[base]
                b0 = p[second]-p[base]; b1 = d[second]-d[base]
                c0 = p[probe]-p[base]; c1 = d[probe]-d[base]
                cross0 = np.cross(a0, b0)
                cross1 = np.cross(a0, b1)+np.cross(a1, b0)
                cross2 = np.cross(a1, b1)
                coefficients = [float(cross2@c1), float(cross2@c0+cross1@c1),
                                float(cross1@c0+cross0@c1), float(cross0@c0)]
                for root in _cubic_roots(coefficients):
                    position = start+root*velocity
                    if is_vertex_face:
                        bary = _closest_point_triangle(position[nodes[0]][None], position[nodes[1]][None],
                                                       position[nodes[2]][None], position[nodes[3]][None])[0]
                        closest = bary@position[nodes[1:]]
                        gap = float(np.linalg.norm(position[nodes[0]]-closest))
                        weight = np.concatenate([[1.], -bary])
                    else:
                        s, t = _closest_segment_segment(position[nodes[0]][None], position[nodes[1]][None],
                                                        position[nodes[2]][None], position[nodes[3]][None])
                        point_a = position[nodes[0]]+s[0]*(position[nodes[1]]-position[nodes[0]])
                        point_b = position[nodes[2]]+t[0]*(position[nodes[3]]-position[nodes[2]])
                        gap = float(np.linalg.norm(point_a-point_b))
                        weight = np.array([1-s[0], s[0], -(1-t[0]), -t[0]])
                    if gap <= self.thickness:
                        # The side the primitives started on decides which way
                        # a tunnelled vertex must be pushed back. A repair that
                        # separated them on the wrong side would leave the sheet
                        # inside out while reporting no overlap.
                        if is_vertex_face:
                            direction = np.cross(p[2]-p[1], p[3]-p[1])
                        else:
                            direction = np.cross(p[1]-p[0], p[3]-p[2])
                        length = float(np.linalg.norm(direction))
                        if length < 1e-12:
                            offset = np.einsum('i,ij->j', weight, p)
                            length = float(np.linalg.norm(offset))
                            if length < 1e-12:
                                continue
                            direction = offset
                        direction = direction/length
                        if is_vertex_face:
                            side = float((p[0]-p[1])@direction)
                        else:
                            side = float(np.einsum('i,ij->j', weight, p)@direction)
                        if side < 0:
                            direction = -direction
                        # Prefer the contact-time weighted separation: an impulse
                        # along it acts at one coincident material point and so
                        # exerts exactly zero torque. The oriented primitive
                        # normal is only the degenerate fallback.
                        offset = np.einsum('i,ij->j', weight, position[nodes])
                        span = float(np.linalg.norm(offset))
                        if span > 1e-12 and float(offset@direction) > 0:
                            direction = offset/span
                        hits.append((root, nodes, weight, is_vertex_face, direction,
                                     position[nodes].copy()))
                        break
        return hits

    def resolve(self, mesh, previous_positions=None):
        """Apply impulses and position repair; return an auditable record.

        ``previous_positions`` enables the continuous test over the step that
        produced ``mesh.x``. Without it only the discrete proximity state is
        enforced and fast tunnelling is not detected: the record says which.
        """
        x = mesh.x
        if x.shape != (self.count, 3) or not np.isfinite(x).all():
            raise ValueError('Finite mesh positions must match the cached topology')
        if mesh.v.shape != x.shape or not np.isfinite(mesh.v).all():
            raise ValueError('Finite mesh velocities must match positions')
        mass = np.broadcast_to(np.asarray(mesh.mass, float), (self.count,))
        if not np.isfinite(mass).all() or np.any(mass <= 0):
            raise ValueError('Finite positive vertex masses required')
        inverse = np.where(mesh.fixed, 0., 1./mass)
        before_kinetic = float(.5*np.sum(mass[:, None]*mesh.v**2))
        momentum_before = np.sum(mass[:, None]*mesh.v, axis=0)
        angular_before = np.sum(np.cross(x, mass[:, None]*mesh.v), axis=0)
        support = np.zeros(3)
        start = x.copy() if previous_positions is None else np.asarray(previous_positions, float)
        if start.shape != x.shape or not np.isfinite(start).all():
            raise ValueError('Previous positions must be finite and match the mesh')
        vertex_face, edge_edge = self._pairs(start, x)
        continuous = self._continuous(start, x, vertex_face, edge_edge) if previous_positions is not None else []
        swept_impulses = 0
        swept_repair = 0.
        torque = np.zeros(3)
        for _, nodes, weight, _kind, direction, contact_positions in continuous:
            swept_impulses += self._impulse(mesh, nodes, weight, inverse, support, start,
                                            normal=direction, reference=contact_positions,
                                            torque=torque)
            # Separate the crossed primitives on their original side. This is a
            # position repair only; it is never turned into velocity.
            scale = float(np.sum(weight**2*inverse[nodes]))
            separation = float(np.einsum('i,ij->j', weight, x[nodes])@direction)
            if scale > 0 and separation < self.thickness:
                correction = (self.thickness-separation)/scale*direction
                x[nodes] += (weight*inverse[nodes])[:, None]*correction
                swept_repair = max(swept_repair, float(np.abs((self.thickness-separation)/scale*weight*inverse[nodes]).max()))
        contacts = 0
        maximum_overlap = 0.
        position_repair = 0.
        clearance = float('inf')
        for sweep in range(self.iterations):
            vertex_face, edge_edge = self._pairs(x, x)
            nodes, weights, normals, gaps, _kind = self._stencils(x, vertex_face, edge_edge)
            clearance = float(gaps.min(initial=float('inf')))
            overlapping = np.flatnonzero(gaps < self.thickness)
            if not len(overlapping):
                break
            for row in overlapping:
                index = nodes[row]; weight = weights[row]
                # Recompute the separation from the CURRENT positions with the
                # sweep's barycentric weights. Reusing the sweep-start normal
                # after earlier corrections have moved these vertices would put
                # the impulse off the closest-point line and leak torque.
                offset = np.einsum('i,ij->j', weight, x[index])
                distance = float(np.linalg.norm(offset))
                if distance >= self.thickness or distance < 1e-12:
                    continue
                normal = offset/distance
                depth = self.thickness-distance
                maximum_overlap = max(maximum_overlap, depth)
                scale = float(np.sum(weight**2*inverse[index]))
                if scale <= 0:
                    continue
                correction = depth/scale*normal
                reference = x[index].copy()
                x[index] += (weight*inverse[index])[:, None]*correction
                position_repair = max(position_repair, float(np.abs(depth/scale*weight*inverse[index]).max()))
                contacts += 1
                self._impulse(mesh, index, weight, inverse, support, None, normal=normal,
                              reference=reference, torque=torque)
        after_kinetic = float(.5*np.sum(mass[:, None]*mesh.v**2))
        return {'vertex_triangle_and_edge_edge_contacts': int(contacts),
                'swept_contacts': int(swept_impulses),
                'max_swept_repair_m': float(swept_repair),
                'continuous_test_applied': previous_positions is not None,
                'max_overlap_m': float(maximum_overlap),
                'max_position_repair_m': float(position_repair),
                'residual_min_distance_m': clearance if not contacts else self.minimum_distance(x),
                'sweeps': int(sweep+1),
                'thickness_m': self.thickness,
                'kinetic_energy_change_j': after_kinetic-before_kinetic,
                'linear_momentum_change_ns': (np.sum(mass[:, None]*mesh.v, axis=0)-momentum_before).tolist(),
                'angular_momentum_change_nms': (np.sum(np.cross(x, mass[:, None]*mesh.v), axis=0)-angular_before).tolist(),
                'impulse_torque_residual_nms': float(np.linalg.norm(torque)),
                'fixed_support_impulse_ns': support.tolist(),
                'scope': SCOPE}

    def _impulse(self, mesh, nodes, weight, inverse, support, start, normal=None,
                 reference=None, torque=None):
        """Perfectly inelastic normal impulse; zero-sum weights conserve momentum.

        The impulse acts along the closest-point separation, so the two weighted
        material points coincide in that direction and the net torque about any
        origin is zero. The residual is accumulated rather than assumed.
        """
        if reference is None:
            reference = (mesh.x if start is None else start)[nodes]
        if normal is None:
            offset = np.einsum('i,ij->j', weight, reference)
            length = float(np.linalg.norm(offset))
            if length < 1e-12:
                return 0
            normal = offset/length
        approach = float(np.sum(weight*(mesh.v[nodes]@normal)))
        if approach >= 0:
            return 0
        scale = float(np.sum(weight**2*inverse[nodes]))
        if scale <= 0:
            return 0
        magnitude = -approach/scale
        mesh.v[nodes] += (weight*inverse[nodes])[:, None]*(magnitude*normal)
        if torque is not None:
            torque += np.cross(np.einsum('i,ij->j', weight, reference), magnitude*normal)
        # The mobile primitives receive (sum of their weights) * j * n, so the
        # anchors receive the opposite, which is (sum of the fixed weights) * j
        # * n because every stencil's weights sum to zero.
        fixed = inverse[nodes] == 0
        if np.any(fixed):
            support += np.sum(weight[fixed])*magnitude*normal
        return 1

    def minimum_distance(self, positions):
        """Smallest non-adjacent primitive distance, for penetration acceptance."""
        x = np.asarray(positions, float)
        vertex_face, edge_edge = self._pairs(x, x)
        _, _, _, gaps, _ = self._stencils(x, vertex_face, edge_edge)
        return float(gaps.min(initial=float('inf')))
