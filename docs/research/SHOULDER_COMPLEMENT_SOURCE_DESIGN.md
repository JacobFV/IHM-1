# Source shoulder complement and registration design

The isolated package supplies **15 missing right-shoulder muscle compartments**
from an acquired MoBL-ARMS 4.1 donor, with complete muscle laws, attachments,
function coefficients, wrap definitions and source kinematic dependencies. It
installs nothing. The current plant's bilateral TRIlong/BIClong/BICshort already
cross the shoulder; this package does not rename those as deltoid or cuff muscles,
and does not add another copy of their force.

## Actual source acquisition and reuse status

Previously held Arm26 contains only six arm muscles. The 38 canonical shoulder
and trunk candidates in `lumbar_shoulder_coverage/audit.json` use inferred
extrema/nearest-bone attachments and engineering force estimates; they cannot
supply measured routes. Held cervical sources establish separate cervical/shoulder
anchor dependencies, not a complete dedicated shoulder force complement.

The newly retained public [CEINMS-RT model repository](https://github.com/CEINMS-RT/UpperLimbModel)
is pinned to commit `459e2ebbf47acb72b0ecbc59950a2d0a983d28db`. Its unscaled
`MOBL_ARMS_41.osim` is **868,277 bytes**, not the repository's separate scaled model.
Repository metadata, commit, file listing, README and license are retained.
The exact donor XML names Saul, Murray, Goehler, Daly, Vidt and Crouch. This is a
maintained public redistribution; byte equivalence to the original SimTK release
has **not** been established.

`data/research/shoulder_complement/acquisition.json` records **1,270,521 acquired
source bytes**, URL/length/SHA-256 for every retained response. No geometry meshes,
ZIPs or native libraries were downloaded. Official SimTK listings show the
unimanual 4.1 tutorial ZIP as 63 MB and bimanual ZIP as 31 MB. The selected download
confirmation required login; its response is retained. Neither ZIP is counted as
held data. The JSON-LD combined package size is 95 MB, not a downloaded byte count.

The current [SimTK primary project download page](https://simtk.org/frs/?group_id=657)
states noncommercial-use conditions alongside BSD-3-Clause wording and required
citations. The CEINMS repository has an Apache-2.0 license file; the older
[OpenSim model catalog](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53090607/Musculoskeletal%2BModelssrc%3Dbreadcrumbs-parent)
listed MIT for this model. These conflicting records are retained, not resolved
by choosing the least restrictive label. This local research package does not
assert that the repository license overrides upstream conditions, or authorize
public distribution/deployment. The catalog is a verified page, not a newly held
catalog file; the authoritative current project page and mirror license bytes are
held.

Primary [Saul et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4282829/)
describes a seven-DOF dynamic model with 50 actuators, based on a 50th-percentile
male. The study derived peak force using volumes and measured joint moments from
five healthy young men, fitting specific tension of 50.8 N/cm² for 32 measured
muscles. Routes combine fixed/moving points and wrapping surfaces. These are
population/model priors, not force measurements in IHM's individual anatomy.
The complete public article HTML is held; no subject covariance or raw attachment
sample distribution is supplied here. A [McFarland 2019 primary record](https://pubmed.ncbi.nlm.nih.gov/30835272/)
is also retained because the model repository cites that subsequent development.
The exact acquired XML, rather than an assumed publication version, determines
all extracted parameter values.

## Concrete model inputs

All 15 use `Millard2012EquilibriumMuscle`. Complete subtrees retain the source's
activation/fiber settings and active/passive/tendon curve definitions; the table
is only a compact view. Lengths are metres, force newtons and pennation radians.

| Compartment | Fmax (N) | Optimal fiber (m) | Tendon slack (m) | Pennation (rad) |
| --- | ---: | ---: | ---: | ---: |
| DELT1 | 1218.9 | 0.0976 | 0.0930 | 0.38397244 |
| DELT2 | 1103.5 | 0.1078 | 0.1095 | 0.26179939 |
| DELT3 | 201.6 | 0.1367 | 0.0380 | 0.31415927 |
| SUPSP | 499.2 | 0.0682 | 0.0395 | 0.12217305 |
| INFSP | 1075.8 | 0.0755 | 0.0308 | 0.33161260 |
| SUBSC | 1306.9 | 0.0873 | 0.0330 | 0.34906585 |
| TMIN | 269.5 | 0.0741 | 0.0713 | 0.41887902 |
| TMAJ | 144.0 | 0.1624 | 0.0200 | 0.27925268 |
| PECM1 | 444.3 | 0.1442 | 0.0028 | 0.29670597 |
| PECM2 | 658.3 | 0.1385 | 0.0890 | 0.45378560 |
| PECM3 | 498.1 | 0.1385 | 0.1320 | 0.43633231 |
| LAT1 | 290.5 | 0.2540 | 0.1200 | 0.43633231 |
| LAT2 | 317.5 | 0.2324 | 0.1765 | 0.33161256 |
| LAT3 | 189.0 | 0.2789 | 0.1403 | 0.36651914 |
| CORB | 208.2 | 0.0932 | 0.0970 | 0.47123890 |

The complement contains 44 fixed, 13 moving and 8 conditional path points on
`thorax`, `clavicle`, `scapula` and `humerus`; 23 distinct wrap objects are referenced.
DELT2 moves with shoulder rotation. Pectoral moving points depend on elevation;
latissimus and teres-major routes also contain conditional dependencies on
`elv_angle` and `shoulder_rot`. Source `PathWrap` order, method and index ranges
remain unchanged, including torus, cylinder, sphere and ellipsoid objects. A
straight-line endpoint import would discard material source behavior.

`shoulder_forces.xml` preserves the 15 complete source muscle subtrees under a
ForceSet. `source_dependencies.xml` preserves the full Ground, BodySet, JointSet
and ConstraintSet, including 12 bodies, 12 CustomJoints, 26 coordinates and 13
coordinate couplers. Six base coordinates and seven independent limb coordinates
are distinct from the 13 dependent coordinates. Dependencies are an inventory;
they are **not** a loadable registered IHM Model or a mass addition request.

The donor couples sternoclavicular and acromioclavicular motion to elevation,
uses inverse rotations through phantom bodies, and decomposes shoulder motion
through `shoulder0/1/2`. Default elevation is about 30 degrees and the source
thorax ground rotation is about −90 degrees around y; assuming zero defaults or
renaming these angles to current flexion/adduction/rotation is invalid.

This donor still lacks explicit trapezius, serratus anterior, rhomboid, levator
scapulae and pectoralis-minor actuators. Prescribed scapular rhythm is not active
force control of those girdle muscles. Full independent girdle muscle coverage
needs a further source, not fabricated forces.

## Coherent integration design and gates

The preferred later variant retains a moving, coupled source girdle rather than
silently welding scapular/clavicular paths to torso. It must satisfy these gates:

1. **Register anatomy in explicit frames.** Identify sternoclavicular,
   acromioclavicular, glenohumeral and humeral landmarks in both donor and canonical
   anatomy. Donor meshes are not yet held, so model landmarks alone do not certify
   anatomical alignment. Keep rigid transforms/scales, residuals and uncertainty
   separate from untouched source muscle parameters. A mirrored left model needs
   explicit reflection of vectors, joint axes and wrap quadrants plus validation;
   it is not obtained by suffixing names.
2. **Reconcile kinematics before adding force.** Preserve the donor's elevation
   couplers and moving/conditional point functions. Current three shoulder
   coordinates and source Euler decompositions differ. Define and test a mapping
   of complete thorax-to-humerus orientation and angular-velocity Jacobian across
   the source-supported range, including singular configurations. Decide whether
   the existing interface exposes mapped coordinates or the immutable variant
   explicitly changes its coordinate contract. No such map has been invented here.
3. **Partition mass once with cervical/thoracic work.** Source clavicle mass is
   0.156 kg and scapula 0.70396 kg per side. Their bilateral source-total 1.71992 kg
   is a candidate prior, not an accepted subtraction. First obtain the exact
   existing torso ledger after actor cervical/head/thoracic partitions. Transform
   proposed girdle COM and full inertia into that common frame; subtract mass,
   first moment and origin inertia from that residual owner, then recover residual
   COM/inertia and reject nonphysical tensors. Preserve total mass, COM and full
   inertia. Donor thorax mass is zero and the donor's 0.0001 kg numerical phantom
   bodies must not become extra physical tissue. Existing humerus/forearm/hand
   masses remain their current owners unless separately reconciled. The donor's
   entire 4.77882 kg body inventory must never be added wholesale.
4. **Reconcile existing shoulder paths.** The active Arm26 long biceps/triceps and
   short biceps use registered torso origins; a moving girdle requires auditing
   those origins against source anatomy. Adding the 15 new compartments does not
   justify duplicating or silently replacing the existing six arm actuators per
   side. Source Millard laws and active Arm26 Thelen laws remain distinct. The
   independent Arm26 wrap audit and source properties must be preserved.
5. **Accept geometry and physics in isolation.** First instantiate the unchanged
   donor in the exact engine version and compare source path lengths, conditional
   transitions, wrapping, moment arms and coupled-coordinate Jacobians over a
   bounded pose set. Then test the registered variant against those properties,
   pressure-free geometry and tension-feasible torque cones. Determine singular
   or unsupported ranges explicitly. No arbitrary balancing torque, force fitting
   to a desired pose, passive-law import or asserted shoulder stability follows
   from successful XML parsing. Native engine loading and torque tests have not
   been run in this task.

The immediate next request is bounded donor-model geometry acceptance in a
separate native slot, followed by registration/mass-partition design coordination
with the cervical/thoracic owner. Missing source geometry can be acquired from
pinned public files in a separately bounded batch. No production plant mutation
is authorized by these artifacts.

## Reproduction

```
.venv/bin/python scripts/materialize_shoulder_complement.py
.venv/bin/python scripts/verify_shoulder_complement.py
```

Four light tests verify source hashes, units/values, exact subtree preservation,
wrap/frame/coordinate dependency closure, moving and conditional point retention,
zero installed mass and deterministic generation. They import no native engine,
load no meshes and issue no native commands. Full original donor/source bytes are
preserved, including their formatting; extracted XML changes serialization only.
