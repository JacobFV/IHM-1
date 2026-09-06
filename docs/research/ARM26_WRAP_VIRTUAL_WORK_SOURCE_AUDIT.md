# Arm26 wrap geometry and virtual work

The retained 50-call native receipt (`native-wrap-work-ovst3fef`) shows persistent, sided length–moment-arm discrepancies in bilateral BRA and BIClong. The source audit identifies concrete geometric reasons the stored length and force mapping need not be work-conjugate. It does not change the native library, raw model, tension laws or acceptance threshold. Exact source hashes, four wrap definitions, and the native receipt hash are in `data/research/arm26_wrap_geometry/audit.json`.

## Different definitions in the current implementation

`GeometryPath.cpp:1202` sums ordinary straight-segment distances but substitutes `PathWrapPoint::getWrapLength()` between two points on the same wrapping object. That scalar is supplied separately by each wrap algorithm. `GeometryPath.cpp:409` applies unit-direction tension along segments connecting different mobilized bodies, and `MomentArmSolver.cpp` projects the resulting spatial forces through the system Jacobian and coordinate coupling vector. It does not differentiate the stored wrapping arc length.

`PathWrapPoint.cpp:86` returns zero `getdPointdQ()`. Its velocity at line98 is only the wrapping body's station velocity, excluding sliding of the contact location relative to that body. `GeometryPath::computeLengtheningSpeed()` sums station-to-station speeds; this is another force-consistent convention that need not equal the time derivative of its separately stored approximate length. Explicit moving-point work compensation in `produceForces()` applies to `MovingPathPoint`, not `PathWrapPoint`.

For a truly stationary frictionless wrapped path, contact sliding and surface-path variation cancel in the first variation, so endpoint tangent forces can legitimately produce `Q = −T dL/dq` without explicit contact-point derivatives. This cancellation cannot simply be assumed for nonstationary geometric constructions or an unrelated chord-sum length.

## BRA: axial tangency bypass, not outer iteration tolerance

Both BRA paths contain exactly one WrapCylinder on the corresponding humerus. Its retained radius is 0.0176379067951 m, length 0.0551184587346 m, quadrant `all`; bilateral rigid transforms are exported without alteration. The `hybrid` XML label is not an ellipsoid algorithm selection for this cylinder.

`GeometryPath.cpp:1067` sets `singleWrap = (wrapSetSize==1)`. `WrapCylinder.cpp:469–578` estimates axial contact locations by geometric plane/apex construction. Lines642–647 then store the helical arc length `sqrt((r theta)^2 + (z2−z1)^2)`. Crucially, the axial tangent correction at lines676–693 requires `!singleWrap`. Thus these BRA paths never run that correction, regardless of visualization. With visualization disabled, they additionally retain only two display wrap points; the scalar length remains helical rather than that chord.

The source threshold at line50 is **0.1 degree**, despite nearby comments saying one degree. It applies only to the skipped correction; tightening it alone cannot fix BRA. The 0.002 m segment sampling likewise belongs to the alternate branch. The GeometryPath outer 0.0005 m stop is immaterial because single wraps receive one outer iteration.

A conservative replacement for an accepted, fixed cylinder winding branch has an exact closed form. In wrap coordinates, let endpoint radial distances be rho1/rho2, tangent-plane straight lengths `a=sqrt(rho1²−r²)`, `b=sqrt(rho2²−r²)`, selected projected arc `c=r theta`, and axial endpoint displacement dz. With `H=a+c+b`:

```
L = sqrt(H²+dz²)
z_contact1 = z_endpoint1 + dz*a/H
z_contact2 = z_endpoint1 + dz*(a+c)/H
```

This follows from unrolling the cylinder and allocating the common axial slope over both straight segments and the surface arc. It supplies exact endpoint tangency, length and tension directions from one construction. `cylinder_branch()` in the new audit script implements this reference calculation and requires explicit winding sense; it does not select an anatomical branch.

Its usable domain is strictly exterior endpoints, a continuous nondegenerate selected winding branch, source-compatible quadrant and finite-cylinder contacts. The reference helper checks exterior/nondegenerate inputs; finite length, anatomical admissibility, branch continuity and whole-model registration remain caller preconditions and are **not certified**. A point entering the cylinder, end-cap contact or branch transition must reject or enter an explicitly modeled alternative. The source's rule returning no-wrap only when both contacts lie outside displayed cylinder length must not silently become a new physical end-cap model.

## BIClong: heuristic plane and discretized length

Each BIClong has one humeral WrapEllipsoid, range2–3, quadrant `−y`, with semiaxes 0.0275592293673, 0.0220473834938, 0.0220473834938 m. It is an axisymmetric ellipsoid, not a sphere. The source hybrid construction (`WrapEllipsoid.cpp:334–480`) blends a major-axis method and a sampled fan; blend thresholds are0.7073/0.9 and fan count300. This selected plane changes with the endpoints.

Lines535–550 find tangencies **within that plane**, then trace its intersection with the ellipsoid. `calcTangentPoint()` uses a dimensionless scaled residual-squared threshold1e−8 with outer cap50 and inner cap1000; it is not a native moment-arm tolerance. `CalcDistanceOnEllipsoid()` at822 onward chooses `floor(contact_chord_length/0.001)` segments, capped499, and sums their chord lengths at968. For contact separation below1mm it returns the contact chord directly. Segment-count changes can introduce additional derivative discontinuities.

Even a perfectly integrated generic plane section of this nonspherical ellipsoid need not be a surface geodesic. Consequently, merely refining the 1mm polyline or tangent tolerance is not sufficient to guarantee cancellation of path/plane variation. The small fixture uses the **actual retained semiaxes** with an illustrative oblique section and obtains nonzero geodesic curvature10.6284 m⁻¹; a principal-plane control is zero to floating precision. This demonstrates the structural limitation, not the actual BIClong plane or a quantitative attribution of its observed error. Separating actual chord integration, tangent residual and plane-selection contributions requires native contact/plane/wrap-point observation, which the 50-call receipt does not contain.

The conservative correction is a constrained stationary surface geodesic with free contact endpoints on the retained ellipsoid, using its existing transform, axes, route range and quadrant. The hybrid plane may initialize a branch but cannot constrain the final solution. Enforce endpoint tangency and zero surface-tangential curvature; integrate arc length to a converged error bound. Derive length, lengthening speed, contact forces and moment arms from that same solution. Contact reaction on the humerus must be retained. An alternative explicitly constrained planar guide would require additional physical guide reactions and anatomical evidence; it cannot be introduced as an invisible derivative correction.

## Bounded next acceptance and materialization

First add an isolated read-only native observer for the four paths: actual route endpoints, wrap-frame transform, chosen contact locations, wrap status, surface/polyline length, number of points, and ellipsoid plane. Sample the already retained seed and ±1e−4/±1e−5 rad cases. This is sufficient to test cylinder endpoint stationarity and separate BIClong chord error from nonstationary-plane effects without a new support solve. No native build/run was performed in this audit.

Then implement an **opt-in geometry variant**, preserving the raw source model and all tension/fiber parameters. For cylinder branches, use the exact solution above. For ellipsoid branches, use a converged stationary surface route with declared branch/quadrant/end-point domain. Keep the original94 unaffected muscle paths unchanged, preserve all inertias, and add no guide mass or ad hoc generalized torque. Confirm actual native `−dL/dq`, speed/work and body-wrench mapping agree under both sided finite differences and small closed loops, including nonzero speeds. Reject branch boundaries; do not hide them by interpolation or relaxing the existing1e−6 m/rad gate. Only then rerun the complete98-muscle effective-gradient prerequisite. The earlier Thelen passive-energy normalization issue remains a separate numerical-merit correction.

Reproduce the retained source inventory and analytic fixtures:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit_arm26_wrap_geometry.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_arm26_wrap_geometry.py
```

The fixture passes36 independent endpoint finite-difference checks over two selected winding senses and three axial displacements. Maximum exact-cylinder endpoint-gradient error is1.433×10⁻¹⁰. This verifies the analytic correction's local geometry, not native Arm26 correction acceptance. All inspected OpenSim source files carry their original Apache2.0 notices; retained model provenance remains unchanged.
