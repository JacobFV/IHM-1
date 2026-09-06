# Stationary BIClong surface path and wrapping work

The isolated four-wrap candidate restores sampled native length/force work
consistency without adding corrective torques. The existing exact BRA cylinder
construction is retained; BIClong now solves a stationary path on its original
ellipsoid surface. All four changes use private wrapping-object copies and
muscle-specific references. Original94 other muscle paths, force parameters,
body masses/inertias and joints remain unchanged. The raw model and production
native library are not modified.

## Stationary surface equations

In the retained wrap frame, let `D=diag(1/a_i²)` and parameterize the surface
curve by arc length. The geodesic equation is

```
x' = v
v' = −(vᵀ D v)/(xᵀ D² x) D x
```

The acceleration is normal to the ellipsoid; surface-tangential curvature is
zero. Initial unit tangent is the projected incoming straight segment. Four
shooting unknowns are two contact-chart coordinates, surface arc length and
outgoing straight-segment length. Four residuals enforce initial surface
tangency plus three-dimensional closure of the outgoing tangent ray onto the
second route endpoint. Thus both endpoint tangencies and the intervening surface
stationarity are solved together. This removes the heuristic plane constraint,
not merely its chord-sum error.

The high-accuracy Python reference uses DOP853 and bounded least squares. The
isolated C++ class uses bounded Newton shooting with centered numerical Jacobian,
backtracking,128 RK4 steps and a256-step refinement check. Its residual threshold
is1e−11 in dimensionless/scaled coordinates; refinement endpoint displacement
must stay below1nm and the refined residual below1e−8. These are numerical solver
gates, not anatomical compliance or force parameters.

Native ellipsoid axes remain(0.0275592293673,0.0220473834938,0.0220473834938)m and
quadrant remains−y. The original native contact/wrap result seeds the branch.
The solver requires exterior endpoints, a native arc below40mm, a local contact
chart, corrected contact shifts no greater than2mm, and sampled surface/speed/
quadrant checks. The line search also restricts normalized arc below1.81,
approximately39.906mm; this is an explicit numerical exploration bound. It
rejects nonconvergence, departure from this neighborhood, unsupported no-wrap
cases that intersect the ellipsoid, or geometry-domain violations. It does not
silently return the old path. These restrictions identify a bounded local
stationary branch; they do not prove global shortest-route uniqueness or a
complete physiological range.

`c1` and `sv` fields recorded for the corrected ellipsoid remain the normalized
**native seed-plane** quantities. They are not a plane of the corrected geodesic;
using the earlier planar-section analyzer on those corrected points would be
incorrect. The corrected surface points and contacts are recorded separately
in the same result.

## Actual native checks

The fixture compiled and completed104 raw plus104 candidate poses in12.760s
under4GiB, one thread and nice10. No physical time advanced. Its native result
is `data/derived/stationary-wrap-candidate-4xo9l1c6/report.json`.

- 400 accepted corrected-wrap observations pass the original1e−6m virtual-work
 gate, with maximum moment-arm versus length-derivative error2.776e−10m.
- All9,776 observations of the other94 muscle lengths/statuses remain unchanged
 within1e−12m.
- At diagnostic coordinate speed0.07rad/s, maximum native lengthening-speed
 versus length-derivative discrepancy is1.943e−11m/s. Maximum unit-tension body
 power plus lengthening speed is1.246e−17W.
- Maximum unit-tension net force is1.570e−16N and net moment5.663e−16Nm, preserving
 action/reaction on the wrapping body and route endpoints.
- All224 local and32 held-out corrected-wrap observations are accepted. Sixteen
 extreme BIClong source-range samples reject explicitly: eight unsupported
 no-wrap branches, six combined quadrant/surface/speed guards, and two contact
 neighborhood departures. Those cases are not counted as accepted work tests.
- The executable directly exercises four cylinder rejections: inside endpoint,
 finite-cap contact, intersecting unsupported no-wrap state and incompatible
 extra winding. All four reject in both raw and candidate runs.

Independent DOP853 reference comparisons at the two retained BIClong seeds
agree with native contacts/arc within3.555e−13m. Eight repeated seed cases per
side, interspersed with other sampled configurations, agree in entry contact
within2.75e−15m. This is sampled repeat/branch-continuity evidence, not a global
history-independence or trajectory-integrator proof. Independent endpoint
finite differences of the reference length agree with unit endpoint tangent
forces within6.74e−11. Native entry contacts shift79µm right and181µm left from
the heuristic plane construction.

The native probe uses unit tension and diagnostic velocities; it does not
choose activation to settle a support pose. Force-law XML remains exact, but
full98-muscle fiber equilibrium, metabolic/chemical ownership and effective
potential acceptance remain separate gates. A consistent wrapped geometry does
not itself correct the independently identified Thelen passive-energy primitive
or certify contact-energy gradients.

## Reproduction and promotion boundary

Committed receipts are in `data/research/stationary_wrap_candidate/`, including
all16 rejection messages and hashes of complete native observations/models.
Light checks require no native process:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_stationary_wrap_materializer.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_ellipsoid_geodesic_reference.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_stationary_wrap_acceptance.py
```

`scripts/run_stationary_wrap_probe.py --prepare` makes an isolated copied model
and build recipe. `--run <directory>` requires the shared heavy slot and uses
archived attested libraries with bounded process-group timeout/kill/reap.
Loading the candidate XML requires its experimental class registrations; this
is not an installed production plant variant. Independent review, explicit
model/build identity and complete effective-gradient acceptance are required
before coupling it to the support solver or brain-driven runtime. Rejected
source-domain cases must remain rejected unless a separate physical route or
expanded stationary branch is justified and verified.
