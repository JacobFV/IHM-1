"""Experimental collision geometry adapter for the passive cloth solver.

All sampled spheres participate, including continuous straight-line relative
node/sphere sweeps. A detected sphere contributes its tangent at the ORIGINAL
node position: this halfspace contains the original feasible node and excludes
the entire sphere at both endpoints and along the translating-plane sweep.
Every outer solve starts from the unchanged original mesh. Failure is atomic.

These are vertex contacts, not triangle collision or friction. Sphere motion is
prescribed constant translation over one step. Bounded boxes are TOP supports
only: their side/underside geometry is absent, and entry below the top is rejected
if an original-state-feasible top constraint cannot be imposed. Constraints are
conservative and may reject feasible curved paths; there is no position repair.
"""
from copy import copy
import numpy as np
from .cloth_passive_step import PassiveStepRejected


def _vector(value, name):
    x = np.asarray(value, dtype=float)
    if x.shape != (3,) or not np.isfinite(x).all():
        raise ValueError(name + ' requires finite three-vector')
    return x.copy()


def _box_interval(start, end, lower, upper, axes):
    """Closed line/rectangle interval in normalized step time, or None."""
    lo, hi = 0., 1.
    for axis in axes:
        delta = end[axis] - start[axis]
        if abs(delta) < 1e-15:
            if start[axis] < lower[axis] or start[axis] > upper[axis]:
                return None
        else:
            a, b = sorted(((lower[axis]-start[axis])/delta,
                           (upper[axis]-start[axis])/delta))
            lo, hi = max(lo, a), min(hi, b)
            if lo > hi:
                return None
    return lo, hi


class ClothPassiveGeometry:
    def __init__(self, stepper, max_outer_iterations=8, tolerance_m=None, up_axis=2):
        if type(max_outer_iterations) is not int or max_outer_iterations < 1:
            raise ValueError('positive integer outer iteration limit required')
        if type(up_axis) is not int or up_axis not in (0, 1, 2):
            raise ValueError('up_axis must be 0, 1 or 2')
        self.stepper = stepper
        self.max_outer_iterations = max_outer_iterations
        self.tolerance = stepper.tolerance if tolerance_m is None else float(tolerance_m)
        if not np.isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError('positive finite geometry tolerance required')
        self.up_axis = up_axis

    def step(self, mesh, dt_s, *, spheres=(), floor_height_m=None, boxes=(),
             gravity=(0., 0., -9.81), forces_n=None):
        h = float(dt_s)
        if not np.isfinite(h) or h <= 0:
            raise ValueError('positive finite step required')
        x0 = np.asarray(mesh.x, float).copy()
        v0 = np.asarray(mesh.v, float).copy()
        if x0.ndim != 2 or x0.shape[1] != 3 or v0.shape != x0.shape or not np.isfinite(x0).all() or not np.isfinite(v0).all():
            raise ValueError('finite N by 3 mesh state required')
        g = _vector(gravity, 'gravity')
        up = np.eye(3)[self.up_axis]
        axes = [i for i in range(3) if i != self.up_axis]
        sphere_data = []
        for i, sphere in enumerate(spheres):
            center = _vector(sphere['center_m'], 'sphere center')
            velocity = _vector(sphere.get('velocity_m_s', [0.,0.,0.]), 'sphere velocity')
            radius = float(sphere['radius_m'])
            if not np.isfinite(radius) or radius <= 0:
                raise ValueError('sphere radius must be finite and positive')
            relative = x0-center
            distance = np.linalg.norm(relative, axis=1)
            if np.any(distance < radius-self.tolerance) or np.any(distance == 0):
                raise PassiveStepRejected('Initial sampled sphere geometry infeasible', {'sphere':i})
            sphere_data.append((center, velocity, radius, sphere.get('owner', i), relative/distance[:,None]))
        box_data = []
        for i, box in enumerate(boxes):
            lower = _vector(box['min_m'], 'box min')
            upper = _vector(box['max_m'], 'box max')
            if np.any(upper <= lower):
                raise ValueError('box bounds must have positive extent')
            inside = np.all((x0[:,axes] >= lower[axes]) & (x0[:,axes] <= upper[axes]), axis=1)
            if np.any(inside & (x0[:,self.up_axis] < upper[self.up_axis]-self.tolerance)):
                raise PassiveStepRejected('Initial top-cover box geometry infeasible', {'box':i})
            box_data.append((lower, upper, box.get('owner', 'box:'+str(i))))
        planes = {}
        if floor_height_m is not None:
            floor = float(floor_height_m)
            if not np.isfinite(floor):
                raise ValueError('finite floor height required')
            if np.any(x0[:,self.up_axis] < floor-self.tolerance):
                raise PassiveStepRejected('Initial floor geometry infeasible')
            for node in range(len(x0)):
                planes[('floor',node)] = dict(node=node, normal=up.tolist(), offset_m=floor, owner='floor')

        def violations(end):
            keys = []
            minimum_gap = float('inf')
            for i, (center, velocity, radius, owner, normals) in enumerate(sphere_data):
                start_relative = x0-center
                delta = end-x0-h*velocity
                length2 = np.sum(delta*delta, axis=1)
                fraction = np.clip(-np.sum(start_relative*delta, axis=1)/np.maximum(length2, 1e-30), 0., 1.)
                gap = np.linalg.norm(start_relative+fraction[:,None]*delta, axis=1)-radius
                minimum_gap = min(minimum_gap, float(np.min(gap, initial=float('inf'))))
                keys.extend(('sphere',i,int(node)) for node in np.flatnonzero(gap < -self.tolerance))
            for i, (lower, upper, owner) in enumerate(box_data):
                for node, (start, finish) in enumerate(zip(x0, end)):
                    interval = _box_interval(start, finish, lower, upper, axes)
                    if interval is None:
                        continue
                    heights = [start[self.up_axis]+t*(finish[self.up_axis]-start[self.up_axis]) for t in interval]
                    if min(heights) < upper[self.up_axis]-self.tolerance:
                        keys.append(('box',i,node))
            return keys, minimum_gap

        def add_constraints(keys):
            added = 0
            for key in keys:
                if key in planes:
                    continue
                kind, index, node = key
                if kind == 'sphere':
                    center, velocity, radius, owner, normals = sphere_data[index]
                    normal = normals[node]
                    planes[key] = dict(node=node, normal=normal.tolist(),
                        offset_m=float(normal@center+radius), velocity_m_s=velocity.tolist(), owner=owner)
                else:
                    lower, upper, owner = box_data[index]
                    if x0[node,self.up_axis] < upper[self.up_axis]-self.tolerance:
                        raise PassiveStepRejected('Top-only box entry needs an infeasible original-state tangent; reduce step or use side collision geometry', {'box':index,'node':node})
                    planes[key] = dict(node=node, normal=up.tolist(), offset_m=float(upper[self.up_axis]), owner=owner)
                added += 1
            return added

        # Seed with a ballistic sweep, then validate the actual coupled solve.
        predicted = x0+h*v0+h*h*g
        add_constraints(violations(predicted)[0])
        for iteration in range(1, self.max_outer_iterations+1):
            trial = copy(mesh)
            trial.x = x0.copy(); trial.v = v0.copy()
            diagnostics = self.stepper.step(trial, h, gravity=g, planes=list(planes.values()), forces_n=forces_n)
            missed, min_gap = violations(trial.x)
            if not missed:
                # The solver has already checked floor halfspaces and energy.
                # The geometric sweep verifies the curved primitives separately.
                mesh.x[:] = trial.x; mesh.v[:] = trial.v
                diagnostics['geometry'] = dict(outer_iterations=iteration,
                    sphere_count=len(sphere_data), box_count=len(box_data),
                    contact_constraints=len(planes), minimum_swept_sphere_gap_m=None if not np.isfinite(min_gap) else min_gap,
                    sphere_motion='prescribed constant translation; start-state tangents',
                    scope='all-sphere vertex sweeps; static floor and bounded box tops; no triangle, side, underside, friction or self collision')
                return diagnostics
            if not add_constraints(missed):
                raise PassiveStepRejected('Constrained sphere/box sweep remains infeasible; no state changed', {'missed_contacts':missed})
        raise PassiveStepRejected('Geometry contact-set iteration limit; no state changed', {'outer_iterations':self.max_outer_iterations})

    def advance(self, mesh, dt_s, *, min_dt_s=.000625, spheres=(), floor_height_m=None,
                boxes=(), gravity=(0., 0., -9.81), forces_n=None):
        """Advance one interval atomically with rejection-driven bisection.

        A rejected candidate retries from its original state at two half steps.
        Prescribed sphere centers advance by elapsed time, never by an attempted
        step's endpoint. A failure anywhere leaves the caller's entire interval
        unchanged. Velocities satisfy each accepted INTERNAL step's kinematics;
        endpoint velocity is not the outer interval's average displacement rate.
        """
        h = float(dt_s); minimum = float(min_dt_s)
        if not np.isfinite(h) or h <= 0 or not np.isfinite(minimum) or minimum <= 0:
            raise ValueError('positive finite interval and minimum step required')
        trial = copy(mesh)
        trial.x = np.array(mesh.x, dtype=float, copy=True)
        trial.v = np.array(mesh.v, dtype=float, copy=True)
        # Materialize once: generators and caller-owned records must remain
        # stable across retries and accepted internal steps.
        originals = []
        for sphere in spheres:
            record = dict(sphere)
            record['center_m'] = _vector(record['center_m'], 'sphere center')
            record['velocity_m_s'] = _vector(record.get('velocity_m_s', [0.,0.,0.]), 'sphere velocity')
            originals.append(record)
        box_records = [dict(box) for box in boxes]
        accepted = []
        rejected = 0

        def segment(elapsed, interval):
            nonlocal rejected
            moving = [{**s, 'center_m':s['center_m']+elapsed*s['velocity_m_s']} for s in originals]
            try:
                result = self.step(trial, interval, spheres=moving,
                    floor_height_m=floor_height_m, boxes=box_records,
                    gravity=gravity, forces_n=forces_n)
            except PassiveStepRejected as error:
                rejected += 1
                half = interval*.5
                if half < minimum*(1.-1e-12):
                    raise PassiveStepRejected('Adaptive geometry reached minimum step; whole interval unchanged',
                        {'failed_elapsed_s':elapsed, 'failed_dt_s':interval,
                         'minimum_dt_s':minimum, 'rejected_attempts':rejected,
                         'cause':str(error), 'cause_diagnostics':error.diagnostics}) from error
                segment(elapsed, half)
                segment(elapsed+half, half)
            else:
                accepted.append((elapsed, interval, result))

        segment(0., h)
        first = accepted[0][2]; last = accepted[-1][2]
        boundary_work = sum(r['prescribed_boundary_work_j'] for _, _, r in accepted)
        external_work = sum(r['external_force_work_j'] for _, _, r in accepted)
        energy_excess = last['energy_after_j']-first['energy_before_j']-boundary_work-external_work
        if not np.isfinite(energy_excess) or energy_excess > self.stepper.energy_tolerance:
            raise PassiveStepRejected('Adaptive interval energy audit failed; whole interval unchanged',
                {'energy_excess_j':float(energy_excess)})
        reactions = []
        for elapsed, interval, result in accepted:
            for reaction in result['plane_reactions']:
                reactions.append({**reaction, 'sample_elapsed_s':elapsed, 'interval_s':interval})
        residuals = np.asarray([r['momentum_residual_ns'] for _, _, r in accepted])
        support_impulse = np.sum([r['impulse_ns'] for r in reactions],axis=0) if reactions else np.zeros(3)
        result = dict(accepted=True, solver_success=True,
            solver_message='All internal steps accepted with unchanged physical tolerances',
            iterations=sum(r['iterations'] for _, _, r in accepted),
            constraint_violation_m=max(r['constraint_violation_m'] for _, _, r in accepted),
            stationarity_residual=max(r['stationarity_residual'] for _, _, r in accepted),
            energy_before_j=first['energy_before_j'], energy_after_j=last['energy_after_j'],
            prescribed_boundary_work_j=float(boundary_work), external_force_work_j=float(external_work),
            energy_excess_j=float(energy_excess), momentum_residual_ns=residuals.sum(axis=0).tolist(),
            max_inner_momentum_residual_ns=float(np.max(np.abs(residuals),initial=0.)),
            max_extension_ratio=max(r['max_extension_ratio'] for _, _, r in accepted),
            plane_reactions=reactions, total_support_impulse_ns=np.asarray(support_impulse).tolist(),
            accepted_substeps=len(accepted), finest_dt_s=min(interval for _, interval, _ in accepted),
            rejected_attempts=rejected,
            substeps=[dict(elapsed_s=elapsed, dt_s=interval, geometry=r['geometry']) for elapsed, interval, r in accepted],
            scope='Atomic adaptive conservative sphere/top geometry; exact summed work and impulses; kinematics apply to each internal step, not outer displacement')
        mesh.x[:] = trial.x; mesh.v[:] = trial.v
        return result
