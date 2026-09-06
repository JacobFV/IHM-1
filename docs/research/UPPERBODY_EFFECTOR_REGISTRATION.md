# Source-derived upper-body actuator extension

The materialized `data/derived/mechanics/whole_body_arm26_v1/subject_with_arms.osim` contains 92 muscles: the retained 80 lower-limb muscles and 12 bilateral upper-arm actuators. This is an XML source construction with checked registration and controller ports. **Native force, moment arms, range of motion and stable movement are not yet verified.** The native mechanical agent owns those acceptance tests.

## Evidence and what was missing

The canonical mechanics file declares 443 muscle actuators: 80 transferred OpenSim lower-limb muscles and 363 geometry-derived muscles. Its 24 additional tendon/ligament paths are not muscles. The geometry-derived attachment rule projects principal surface extrema onto nearby bones, and its strength estimates use a generic PCSA/tension prior. These are useful synthesis proposals but do not establish real upper-body muscle attachment or strength calibration.

The held [Arm26 source](https://github.com/opensim-org/opensim-models/blob/master/Models/Arm26/arm26.osim) provides six upper-arm Thelen2003 muscles, anatomical path points, wrapping surfaces and muscle parameters. It credits the OpenSim adaptation team and Holzbaur, Murray and Delp, with CC BY 3.0 terms retained in the XML. The [original upper-extremity paper](https://nmbl.stanford.edu/publications/pdf/Holzbaur2005.pdf) describes a richer model; Arm26 is its educational reduction, not that complete model. [Thelen's muscle-mechanics paper](https://uwnmbl.engr.wisc.edu/pubs/jbme03.pdf) and the retained OpenSim implementation identify the force/activation model. Source implementation hashes, XML hashes, credits, primary URLs and remaining sources are retained in `data/sources/sensorimotor/upperbody_sources.json`. No new bulk acquisition or native build was needed.

The new names are `arm26_{TRIlong,TRIlat,TRImed,BIClong,BICshort,BRA}_{r,l}`. Source Fmax values are 798.52, 624.30, 624.30, 624.30, 435.56 and 987.26 N respectively. These remain **source-model parameters**, not person-specific measurements. Activation/deactivation constants, contraction/force-length/force-velocity parameters and pennation are copied from the donor. Fiber and tendon slack lengths scale with the registration length ratio. No generic torque actuator substitutes for these muscles.

## Explicit registration

The source and target shoulder/elbow joint frames provide the registration landmarks. Each side uses uniform scale from humeral landmark distance and a proper rotation aligning the humeral axis. Source base points map about the donor shoulder to the native torso's acromial point. Humeral points map into the existing humerus frame. The source fuses ulna, radius and hand; triceps/brachialis insertions map to native ulna, while biceps insertions map to radius using the native neutral radioulnar frame offset. This biceps mapping is an explicit named-anatomy assignment beyond the donor's fused body label.

The left side is a sagittal reflection prior. It is not an independently measured left arm. Four wrapping surfaces per side (two cylinders and two ellipsoids) remain present. Sizes and centers follow the same similarity transform as their bodies. Orientation is transformed using the source's body-fixed XYZ convention; mirrored local geometry restores a proper rotation and reverses the local z quadrant where required. The source `WrapObject.cpp` explicitly uses `setRotationToBodyFixedXYZ`, which anchors this convention. Wrapping references are renamed consistently with the new objects.

All native bodies, joint XML, body masses and inertias remain unchanged. This prevents duplicating donor arm mass. The registration is anchored at zero arm coordinates; parent/child native offset orientations are checked to cancel there. Matching two landmarks does **not** validate wrapping or moment arms across a three-dimensional shoulder motion. Arm26 omits rotator-cuff/deltoid/pectoral/latissimus recruitment and independent scapular motion; the twelve actuators cannot establish full shoulder control or all upper-body motion.

## Runtime and controller integration

`registration.json` binds source XML, generated model, generated catalog and builder identity. `catalog.json` provides source/target path points, source muscle law, Fmax/length parameters, attachment bodies and explicit decoder/registration assumptions.

```python
manifest = json.loads((root / 'data/derived/mechanics/whole_body_arm26_v1/registration.json').read_text())
rows = whole_body_effector_catalog(root, manifest)
controller = SensorimotorController.from_root(root, muscle_catalog=rows)
```

The catalog reader checks exact model/catalog/source hashes and compares every effector identity, force normalization and optimal length against the generated XML. It does not modify the default 80-muscle plant. The caller must load the matching 92-muscle model. The native adapter must replace only the original 80 fitted lower-limb paths; the added 12 retain actual `GeometryPath` and wrapping objects.

All twelve new muscles have individually addressed, delayed, IBM-mediated descending excitation ports. They receive zero command unless requested; ankle reflex laws are not broadcast onto the arms. The existing decoder remains an engineered uncalibrated mapping. Native mechanics owns activation and force; there is no second controller activation law.

## Neck and wrist: concrete next construction

The held `Neck3dof_point_constraint.osim` contains only five Schutte muscle elements, 26 bodies, 27 named coordinates, 21 coordinate couplers and a point constraint. Its credits/publications are placeholders. It is an engine test fixture, not a complete calibrated neck atlas. Seven cervical masses sum to 1.63 kg; skull mass is 4.7 kg. Two 1 kg auxiliary bodies belong to the fixture construction and must not be imported as anatomical tissue. Its many coordinates are coupled; they do not justify claiming 27 independent cervical DOFs.

A defensible next native neck extension is to preserve the source's three master neck rotations and cervical coupling functions, map T1 onto the torso, retain seven physical cervical segments and a skull, and bind the five source paths to their exact named segments. The head/neck mass and inertia must be removed from the current lumped torso with a positive residual inertia check and parallel-axis accounting. Auxiliary frames should remain kinematic supports, not extra material owners. The unilateral incomplete muscle set requires additional source muscles or explicit mirrored priors before agonist/antagonist coverage; passive joint laws and equilibrium need independent verification. Placeholder source attribution must be resolved before treating its numerical coefficients as published anatomical calibration.

The held [WristModel](https://github.com/opensim-org/opensim-models/blob/master/Models/WristModel/wrist.osim) has 25 deprecated Schutte actuators, 28 bodies and 14 named coordinates. It includes both `ECU_pre-surgery` and `ECU_post-surgery`: these are alternative routes for one muscle, not simultaneous additive strength. It credits Gonzalez, Buchanan and Delp's 1997 wrist-architecture/moment-arm model. The native hand is currently lumped, so the next construction needs the donor carpal/metacarpal/thumb/index joint tree, its couplers and paths, with the existing hand's mass reassigned rather than duplicated. Finger muscles connected to welded digits do not establish independently articulated fingers. Resolve the surgical ECU alternative and deprecated-law compatibility before creating active catalog entries.

The richer MoBL-ARMS dynamic model has been identified at the primary [SimTK project](https://simtk.org/projects/upexdyn); the web request returned 403 and no model bytes were acquired. It is recorded as unacquired. Respiratory muscle force remains under the coordinated respiratory-native interface; neither the Arm26 model nor the canonical inferred diaphragm path identifies distributed diaphragm/intercostal actuation.

## Bounded checks and retained results

`verify_upperbody_effectors` verifies source bytes remain unchanged, the 92-muscle/8-wrap count, exact unchanged joint tree/mass, machine-precision elbow landmark alignment, proper registration scale, individual right-arm command without contralateral/antagonist recruitment, and rejection of a changed materialized model. The preceding sensorimotor tests remain applicable.

```
OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 -- nice -n 10 .venv/bin/python -m scripts.verify_upperbody_effectors
```

The source build took 0.19 s and 71,964 KiB peak RSS; the bounded test used 80,660 KiB peak RSS. Both ran under a 1 GiB address-space limit. No native compilation, optimizer, browser, full-body simulation or live IBM access occurred. Native mechanical acceptance must record actual force/moment arms and bilateral registration behavior before these XML checks can support a movement claim.
