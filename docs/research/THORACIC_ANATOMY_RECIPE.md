# Source-registered thoracic mechanical recipe

`data/research/thoracic_anatomy/v3/manifest.json` is an opt-in anatomy and inertia
recipe. It retains full geometry, exact source node/face bindings, 24 candidate
rib hinges, a sternum translation coordinate and a diaphragm material descent
mode. It supplies concrete geometry and mass inputs for a reduced mechanical
implementation. It does not yet implement a loadable coupled native thorax,
identify stiffness, or authorize native activation.

## Frozen source and frames

The source is the factory's frozen `canonical_mechanics.json` under
`data/derived/audits/cutaneous-factory-ylro2d66/body/mechanics`, already bound by
`data/research/cervical_inertia/v2/manifest.json`. Original compressed geometry
bytes and SHA256 identities are retained. The 63 structures contain 611,686
vertices and 852,148 source triangles. No display decimation was used:

| Source structures | Count | Mechanical treatment in this recipe |
| --- | ---: | --- |
| Ribs | 24 | Candidate rigid carriers with local COM/inertia and posterior hinge |
| Costal cartilages | 14 | Registered links; stiffness and deformation law absent |
| Sternum pieces | 3 | Compound anterior-translation candidate |
| Diaphragm | 1 | Source material surface and fixed-anchor descent weights |
| Intercostal muscle groups | 6 | Material shares and 66 candidate adjacent-rib attachment pairs |
| Thoracic/lumbar vertebrae | 15 | Posterior geometric supports; mass stays in parent torso |

The retained BodyParts3D geometry is attributed to DBCLS under CC BY 4.0.
Canonical meters are mapped by the **same frozen proper-rigid torso embedding**
used by the cervical recipe. There is no additional fit, scaling, snapping or
moving-pose registration. This inherits the global 62.356 mm COM/envelope proxy
registration residual; it is not a measured anatomical landmark registration.
Native model identity and prerequisite cervical partition are explicit.

## Attachment and joint materialization

Every endpoint records an entity ID, original vertex index, original triangle
index, one-hot barycentric coordinates and reference position in the native
torso frame. The nearest-node method is reproducible but does not certify the
closest continuous surface point or identify a measured enthesis. Reference
gaps are preserved rather than closed silently.

Each rib is paired by source ordinal with its thoracic vertebra. The nearest
5% of rib vertices form a posterior patch; its centroid and first principal
axis define a **candidate** hinge. Axis sign is deterministic and is not
labeled inspiration. Maximum nearest rib/vertebra gap is 0.340 mm, but proximity
does not establish the actual costovertebral/costotransverse joint axis. Two
patches have first/second singular-value ratios below 2 and are marked
`ambiguous_geometric_axis`. The ratio threshold is an engineering conditioning
screen, not an anatomical validation threshold. Floating-rib mechanics remain
explicitly uncalibrated. Each rib also has a proposed body frame and inertia
about its COM expressed locally, so a future explicit PinJoint conversion need
not infer frames again.

Costochondral, cartilage/sternum and sternum-piece candidates retain 30 point
pairs. No cartilage material prior from the canonical display model is enabled.
The diaphragm has 16 candidate lower-rib/xiphoid/lumbar connections, with a
maximum gap of 4.163 mm. These costal/sternal/crural associations are anatomical
hypotheses expressed as source-node proximity, not recovered fiber insertions.
The atlas does not identify central-tendon fibers or a calibrated abdominal
support model in this recipe.

The six intercostal surfaces are aggregate groups, not 66 identified individual
muscle fascicles. Bilateral grouping is inferred from canonical X position.
Each group is tested against adjacent ribs to produce upper/lower candidates.
Nine of 66 pairings are marked `unsupported_by_source_proximity`: an endpoint
gap exceeds half the corresponding rib-to-rib span. The maximum endpoint gap
is 50.286 mm. The other 57 are **unvalidated geometric candidates**, not accepted
muscle insertions. No fiber orientation, pennation, physiological cross section,
optimal length, force or stiffness is fabricated from these point pairs.

## Material mass and inertia

Original canonical `mass_role=material_partition_proxy` shares are retained as
explicit engineering allocations. Their open-surface volume and generic
allocation uncertainty remain visible. The original masses are not treated as
independently measured tissue masses, nor are their display stiffness or
isotropic inertia proxies promoted to mechanics.

The recipe computes exact uniform triangular-lamina COM and second moments,
normalized to each allocated mass. This avoids using a signed integral of an
unclosed volume and avoids filling curved rib concavities with a convex solid.
It is a surface-distribution prior; it does not assert bone density, shell
thickness or tissue volume. Degenerate triangles contribute zero area and are
counted. All material inertia tensors pass physical second-moment checks.

| Material allocation | Mass (kg) |
| --- | ---: |
| Ribs | 1.075515 |
| Costal cartilage | 0.120888 |
| Sternum | 0.215785 |
| Diaphragm | 0.383416 |
| Intercostal groups | 0.502449 |
| Total debit | 2.298053 |

This is a **composed prospective partition**. First apply the exact cervical
ledger, leaving torso mass 20.235757397 kg; then debit the thoracic shares,
leaving 17.937704388 kg. Do not apply the thoracic replacement directly to the
current unmodified 27.654676965 kg native torso. Posterior vertebral supports
are read for geometry but remain within the parent mass. Recombination gives
mass error zero, COM error 5.65e-17 m and inertia error 2.24e-16 kg m². The
remaining torso second-moment eigenvalues are positive: 0.339390, 0.386191 and
0.702257 kg m². No mode adds mass to the whole-body total.

## Concrete kinematic starting mechanism

The candidate coordinates are 24 independent rib hinge angles (radians), one
anterior sternum translation (meters), and one diaphragm descent (meters).
`hinge_motion` returns actual rotated source positions and their analytic
Jacobian, including the angular motion missing from the existing display modes.
No physiological admissible angle range or closed-chain constraint is claimed.

For the diaphragm, each source node has a descent weight equal to squared
nearest-anchor distance divided by the squared maximum distance. Candidate
anchor nodes have weight zero and the largest weight is one. The descent axis
is canonical inferior transformed into the same torso frame. The compressed
weight field is separately hashed. `diaphragm_motion` returns positions and
its exact material Jacobian. The map is an engineered fixed-anchor reduction;
its strain/folding domain is unvalidated and moving rib anchors must be coupled
before physical use.

Area integration of the squared, linearly interpolated nodal weights gives a
fixed-anchor diaphragm modal kinetic coefficient of **0.042538251 kg**. This is
derived from the allocated 0.383416468 kg surface material, not an added point
mass. Generalized torso cross terms, moving attachments and intercostal/
cartilage material motion must be included in a complete kinetic model.

The next native conversion can construct the rib carriers and explicit hinge
frames, combine the three sternum pieces, and implement the diaphragm material
coordinate with consistent mass/Jacobians. It must also resolve the two
ambiguous axes, reject/replace unsupported attachments, close moving rib/sternal
constraints, and define material motion for every debited share. Debiting mass
without representing its subsequent kinetic motion would be invalid.

## Respiratory ownership and remaining acceptance

This recipe follows [the mechanical-port audit](RESPIRATORY_MECHANICAL_PORT_AUDIT.md).
It does not compute an enclosed cavity from open tissue meshes or relabel a
convex envelope as native gas volume. A separate cavity construction/calibration
must provide `V(q)` and `J_V`, with pressure work `Q_p·qdot = p*Vdot` and the same
material map used by bed/cloth/contact tractions. Source-node bindings provide
a basis for that shared map; resultant torso wrenches alone are insufficient.

Native bilateral chest-wall compliance and ideal respiratory pressure drive
remain the owners until an explicit mechanical transfer/replacement is made.
The signed worker's source work observer uses pressure difference and source
flow; its receipt cannot be assigned to diaphragm muscle chemistry by label.
No added stiffness, damping, active specific tension, ATP cost or neural
recruitment law is enabled here. No native build, simulation or default-model
change was performed.

## Verification

The source builder and checks run under a 1 GiB address-space cap. Six tests
cover analytic lamina moments, hinge Jacobian/virtual work, exact source-node
binding, fixed diaphragm anchors/Jacobian, unsupported inputs, all 63 retained
geometry hashes, all 202 attachment pairs and the composed torso mass ledger.

```
OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 nice -n 10 .venv/bin/python -m scripts.verify_thoracic_anatomy
```

Materialize afresh with `python -m scripts.build_thoracic_anatomy --output
<fresh workspace directory>`. The builder refuses to overwrite an existing
artifact and does not alter native sources.
