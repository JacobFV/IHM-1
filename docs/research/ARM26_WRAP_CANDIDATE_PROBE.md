# Isolated Arm26 wrap observation and cylinder candidate

The source-only fixture adds three registered experimental wrap classes in its
own executable. `ObservedArm26Cylinder` and `ObservedArm26Ellipsoid` call the
unchanged native wrapping routines and record their actual local endpoints,
contacts, status, wrap transform, quadrant, arc and polyline. Ellipsoid c1/sv
are normalized with the returned factor, because native code leaves those
fields internally scaled. No production OpenSim class registration is replaced.

`ExactArm26Cylinder` calls the native cylinder routine first to retain its
observed XY contact/winding choice. For an accepted branch it only replaces
axial contact positions and surface arc with the exact unrolled-cylinder
construction. It requires exterior endpoints, native XY tangencies, consistent
incoming/outgoing winding, a nondegenerate arc less than one revolution,
retained `all` quadrant, and both corrected contacts strictly within finite
length. A native no-wrap result is accepted only when the projected straight
segment clears the infinite cylinder. Cap, inside, near-surface and incompatible
branch cases reject explicitly. This deliberately does not invent a cap law,
change winding, disable a muscle or patch torque after force evaluation.

The materializer changes four observer class tags, adds two private BRA cylinder
clones, and changes only the two BRA wrap references in a copy of the frozen
98-muscle model. Reversing these explicit changes recovers the original XML bytes.
The original cylinders are shared with triceps paths, so changing their force
geometry directly is forbidden; the native scope gate caught that first attempt.
Only two BRA geometries are candidates; BIClong remains observed. All muscle
parameters, original remaining96 muscle paths, body inertias and joints remain
unchanged. The standalone probe performs no muscle equilibration, controller
execution or physical-time advancement: unit-tension geometry, force and moment
arm tests isolate wrapping from the muscle constitutive law.

The prepared schedule has104 cases per variant: both signs at1e−3,1e−4,1e−5 rad,
the seed, source-bound endpoints and quartiles for eight shoulder/elbow
coordinates, and eight held-out combined poses. XML shoulder ranges extend far
beyond a calibrated physiological domain; testing their samples establishes
explicit acceptance/rejection coverage, not physiological ROM or proof over the
whole interval. Each accepted BRA observation compares native moment arm with
centered path-length derivative (second-order inward differences at source
bounds) and checks net body force/moment under unit
tension. The other96 muscle lengths/statuses are compared between variants.
Native assembled coordinates and rejected path messages are retained.

Offline BIClong analysis reports exact observed contact plane/surface residuals,
endpoint tangency, planar-section geodesic curvature, and32/64-point Gauss arc
quadrature against the stored chord sum. Endpoints are projected to the exact
section for quadrature; that displacement is reported separately, so tangent
solver residual cannot masquerade as arc discretization. BIClong geometry is
not replaced in this stage.

Light checks:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_arm26_wrap_materializer.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_arm26_contact_analysis.py
```

Prepare copies and exact identities with
`python scripts/run_arm26_wrap_probe.py --prepare`; run the returned directory
with `--run` only when the shared native slot is assigned. The runner uses the
attested archived loader/libraries, one O0 translation unit,4GiB address-space
limits, one thread, nice10 and a90s combined deadline with process-group kill
and reap. It never updates the production build pointer. Preparation and
analytic tests alone do not certify native candidate acceptance.

## Actual native result

`data/derived/arm26-wrap-candidate-po_0h6uq/report.json` passed the sampled
candidate and untouched-path gates in7.989s, including compilation and both
static runs. All104 cases produced208 accepted BRA observations, with maximum
moment-arm/path-derivative discrepancy2.776e−10m. Unit-tension resultant force
was exactly zero and maximum resultant moment2.498e−16Nm. All9,984 observations
of the other96 muscle lengths/statuses matched within1e−12m. Maximum BRA length
change over the sampled domain was26.077µm. No sampled branch was rejected;
unsupported-domain rejection guards exist, but this result does not claim their
boundary cases were exercised or certify the complete continuous source range.

The actual local cylinder axial tangency residual decreases from0.0144831 to
1.638e−15. No torque correction is applied: changed contacts, exact helical
length and the unchanged native tension-to-body-wrench mapping produce the
improvement. This is an executable isolated candidate. Loading its copied XML
requires the experimental class registrations in this probe; it is not yet a
production adapter/factory variant or a dynamic force-speed acceptance.

Actual BIClong contacts confirm that arc discretization alone is insufficient.
Across the local sided observations, the chosen planar section has geodesic
curvature0.948–1.819m⁻¹. Refining its stored chord arc adds1.688–1.945µm, with
projected contact displacement0.079–0.120µm and32/64-point quadrature agreement
better than7e−18m. Surface-tangency dot residuals remain approximately
3.94e−5–4.55e−5. Correcting only the arc primitive changes right adduction error
from0.000107194 to0.000109818m and left adduction from0.000233557 to0.000236497m;
these still fail the original1e−6m gate. The observed nonstationary surface route
therefore requires a geometry correction, not just more chord samples. No
ellipsoid replacement is implemented or certified here.

Committed evidence in `data/research/arm26_wrap_candidate/acceptance.json` pins
the complete native observations and retains these metrics and six-direction
arc-refinement attribution. `failed_shared_scope.json` preserves the first
native attempt's explicit failure: changing a shared cylinder perturbed18
triceps observations by up to18µm. The final private-clone fix preserves those
paths. An earlier observer-only compile failure (Rotation matrix indexing) is
retained at `arm26-wrap-candidate-rhpz5osx/compile.log`.

Next promotion gates remain: independently review the exact candidate's domain
and XML ownership, run controlled nonzero-speed native length/force work checks,
and integrate only through an explicit build/model identity. For BIClong,
construct a stationary ellipsoid surface route retaining source axes/transform/
quadrant and branch evidence, then repeat the same untouched-path and work
checks. The complete98-muscle energy-gradient prerequisite remains failed until
BIClong and the separately identified numerical-merit issues are resolved.
