# Conservative free thorax equations

`ihm/assembly/thoracic_free_dynamics.py` supplies a loadable source-only right-hand
side for the existing thoracic material mechanism. It includes the moving
parent, all 48 moving material partitions, diaphragm descent, and the two
unresolved rib locks. It introduces no material mass, stiffness, damping,
activation, gravity, or native physiological connection. Loads must be supplied
explicitly as generalized forces in the declared velocity convention.

The subsystem contains **20.235757397 kg**. The actual current native torso is
**27.654676965 kg**: applying this subsystem requires first composing the
prospective **7.418919568 kg** cervical debit and its exact first/second moments.
The remaining thoracic core and moving source lamina masses conserve that
cervical-reduced ledger. Nothing here changes the current 22-body native model.
The source maps and uncertainty remain those of
[the executable thoracic mechanism](THORACIC_EXECUTABLE_MECHANISM.md).

## Mixed velocity coordinates and equations

Let world position be `X=t+R x(q)`, where `R` maps torso-frame vectors to world,
and `J=dx/dq`. The 32-component velocity vector is
`u=(v,omega,s)=(R^T tdot,omega_body,qdot)`. Two internal rib coordinates and
speeds remain locked at zero; the 30 active mass directions must be positive
definite. Angular velocity is a quasi-velocity. It is never interpreted as
Euler-angle derivative.

For each material point:

```
A = [I, -[x]cross, J]
Xdot = R A u
Xddot = R (A udot + c)
c = omega cross v + omega cross (omega cross x)
    + 2 omega cross (J s) + H[s,s]
M udot + b = Q + Qconstraint
b = integral A^T c dm + core gyroscopic contribution
```

Each nonlinear material term is an independent rib rotation with fixed
interpolation weights, so `H[s,s]=sum_j (axis_j cross J_j) s_j^2`.
There are no mixed internal Hessians for this recipe. Linear sternum translation
and diaphragm descent have zero second derivative; their mass and couplings
still enter `M` and the rotational convective terms. The rigid residual core
uses its full COM offset and tensor, including `omega cross (I_core omega)`.
The integration uses the same fixed reference triangle masses and positive
three-point quadrature as the kinetic module. This exactly integrates products
of the piecewise linear material fields within the declared geometric prior.

`evaluate(q,u,Q)` computes the mass matrix, bias, and constrained acceleration.
`solve(evaluated,Q)` reuses that same evaluated state for another explicit load.
It does not reevaluate changed coordinates or velocities. Locked-coordinate
loads produce only formal balancing loads on eliminated columns and do zero
work. Because those material Jacobian columns are zero, neither source point
loads nor inertial loads reconstruct physical locked-rib hinge wrenches. The
API explicitly withholds those physical reactions (`None`).
`pose_rates(R,u)` returns `tdot=R v` and `Rdot=R[omega]cross`. The world origin
acceleration is `R(vdot+omega cross v)`; it is not simply `R vdot`.
The bounded integrator in `thoracic_trajectory.py` preserves proper rotations
using a Rodrigues exponential, advances internal coordinates with `s`, and
enforces locks and the engineering geometry domain. Its second-order explicit
Lie midpoint update is not exactly energy conserving or symplectic. The
initial RHS gives midpoint body speeds and coordinates; the midpoint RHS
updates endpoint speeds, and the midpoint angular speed determines the
full-step SO(3) increment. Cavity facets are checked at midpoint and endpoint;
this does not prove absence of all material self-intersections.

## Verified differential conservation

The retained probe is
`data/research/thoracic_mechanism/free_dynamics_probe_v1/probe.json`. It binds
source files and the actual mechanism manifest. At one nonzero pose/velocity,
it checks the differential free-energy balance using an independent centered
finite difference of the complete mass matrix along `qdot`:

```
dT/dt = u^T M udot + 0.5 u^T Mdot u = u^T Q
```

With body momentum `p=(M u)linear` and angular momentum about the moving parent
origin `h=(M u)angular`, zero external load also satisfies
`pdot+omega cross p=0` and `hdot+omega cross h+v cross p=0`. The latter transport
term is necessary for angular momentum about a fixed world origin. Tests check
these independently of the acceleration residual, all 48 material directional
Hessians against finite differences, SO(3) pose rates, and energy/work balance
under geometric cavity pressure plus an actual source-node force.

This is local differential evidence. The additional short trajectory check
below does not establish long-duration stability. The
cavity remains an engineered lung-envelope geometry with explicit offset
limitations, not measured gas volume. The native chest compliances and ideal
pressure source retain their existing ownership. Recoil and pressure/flow signs
must be explicitly transferred and validated before native coupling; no source
work is assigned to muscle chemistry.

Reproduce under 1 GiB without compilation or native processes:

```
OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 nice -n 10 .venv/bin/python -m scripts.verify_thoracic_free_dynamics
OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 nice -n 10 .venv/bin/python -m scripts.probe_thoracic_free_dynamics --output <fresh directory>
```


## Bounded finite trajectory acceptance

`data/research/thoracic_mechanism/trajectory_probe_v1/probe.json` retains the
actual starting state, midpoint, endpoints, velocities, source receipts and
invariants for 50 ms of free motion with one 50 ms step versus two 25 ms steps.
Parent translation and rotation are nonzero, with angular speed
`[0.6,-0.4,0.5] rad/s`; the same 48 material partitions and two locks are used.

| Invariant absolute error | 50 ms step | Two 25 ms steps |
|---|---:|---:|
| Kinetic energy (J) | 1.42985e-7 | 1.98808e-8 |
| World linear momentum (kg m/s) | 3.25335e-5 | 8.13342e-6 |
| World angular momentum about fixed world origin (kg m²/s) | 1.21834e-5 | 3.01924e-6 |

Momentum errors reduce to about one quarter on refinement, consistent with the
second-order update. Energy error decreases further in this case; no universal
order beyond the method's second order is inferred. Proper rotations, locked
coordinates/speeds, unchanged input ownership and unchanged subsystem mass
also pass. Two trajectory tests ran in 25.9 seconds under 1 GiB without native
processes. The [native composition plan](THORACIC_NATIVE_COMPOSITION_PLAN.md)
retains the missing cervical mass and common-parent coupling requirements.
