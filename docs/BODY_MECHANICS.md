The canonical body now has an executable mechanical model in `ihm/assembly/mechanics.py`, built by `scripts/build_body_mechanics.py` from `data/derived/canonical/anatomy.json`. It combines all 2,403 canonical entities in one connected structural graph: 233 rigid bones, 2,164 affine soft solids, 443 muscle actuators representing 425 anatomical muscle entities, and 24 additional tendon/ligament paths, joined by 3,336 support links. The counts moved from 2,408 entities and 257 rigid bones when the entity-record repair landed: five BodyParts3D surfaces authored twice under two element ids collapsed to one entity each, and 23 intervertebral discs left `rigid_bone` for `cartilage`, so a disc is now the compliance between two vertebrae instead of an infinitely stiff body. See `data/derived/entity-record-repair-promotion-v1/receipt.json`. Multiple Rajagopal actuators share anatomical muscles. Exactly 80 actuator paths and their original source parameters are retained from the native OpenSim export; 363 actuator paths use explicit geometric synthesis priors.

Run from the repository root:

```bash
.venv/bin/python scripts/build_canonical_anatomy.py
.venv/bin/python scripts/build_body_mechanics.py
.venv/bin/python scripts/verify_body_mechanics.py
```

The runtime API is:

```python
import json
from ihm.assembly.mechanics import BodyMechanics
body = BodyMechanics.from_dict(json.load(open('data/derived/canonical/mechanics.json')))
state = body.step(.02, {'activation': {'recfem_r': .1}})
# Alternative independent boundary input:
state = body.step(.02, {'volume_ratios': {'body-bp3d-FJ2422': 1.05}})
```

`activation` accepts actuator IDs, native muscle names, or canonical muscle IDs. Missing activations are zero. `volume_ratios` prescribes a positive ratio in [0.25, 4]; `pressure_pa` prescribes transmural pressure on the stable hydrostatic branch. Supply one boundary type per tissue. `external_forces_n` accepts canonical entity IDs and three-component SI force vectors. Inputs apply for one call; omitting a prior boundary returns that tissue's affine deformation to its unloaded reference. These are explicit inputs, not inferred endogenous physiology.

`state['entities'][id]` contains `translation_m`, `rotation_matrix`, `deformation_gradient`, and current `centroid_m`. Apply transforms about the reference anatomical centroid:

```
x_current = centroid_reference + translation + rotation @ deformation_gradient @ (x_reference - centroid_reference)
```

The source geometry is immutable. All mechanics use meters, kilograms, seconds, newtons, pascals and joules. Native source forces remain archived evidence; they are not transplanted as equilibrated generic-body resting forces.

Bones translate rigidly. Their reference orientations are **constrained**, and the required holding moments are returned in `orientation_reaction_torques_nm`. An exploratory unconstrained rotation scheme failed a sustained activation check on tiny atlas structures; it is not the shipped numerical model. The implementation makes no prediction of articulated joint rotation, gait, posture or bone strain. The generic skeleton starts stress-free, without gravity or active resting tone. Named joints, contact, sliding fascia and physiological pretension are not yet resolved.

Soft solids use one affine compressible neo-Hookean constitutive element per anatomical surface. Its energy density and first Piola stress are:

```
W = mu/2 (F:F - 3) - mu log(J) + lambda/2 log(J)^2
P = mu (F - F^-T) + lambda log(J) F^-T, J = det(F) > 0
```

`neo_hookean` evaluates this law directly. `tetra_force_energy` also provides an actual constant-strain tetrahedral force/energy kernel, tested against finite differences and rigid rotations. The whole-body runtime uses the affine reduction, not a fitted volumetric tetra mesh. Hydrostatic pressure is solved against this law; a prescribed volume sets `J` and returns the required pressure reaction. For cavity entities, the skin parent and the structural lymph graph, the affine transform is a boundary carrier, with zero material volume. Their `constitutive` field is `affine_boundary_carrier`, distinct from an elastic solid; any returned surrogate pressure is a virtual compliance reaction, not measured cavity pressure. A native cavity-volume input is not a myocardium mass-volume measurement.

Muscle shape uses volume-preserving axial contraction with lateral expansion. Its affine elastic stress equilibrates a saturating active stress; maximum free logarithmic shortening 0.35 is an explicit engineering prior. This shape degree of freedom is a quasistatic reduction separate from the force-producing line actuator. Each actuator applies equal/opposite forces to the endpoints of every segment, including intermediate wrapping reactions. Tendon and ligament paths are tension-only elastic lines. The actuator force-length response is a Gaussian reduction with normalized width 0.45 and a passive linear stiffness prior; it is **not** a claim to reproduce Millard muscle equilibrium or native force-velocity dynamics. OpenSim slack/fiber lengths are retained and geometrically scaled, while source maximum isometric force and pennation are retained unchanged. Native tendon compliance semantics remain recorded in the source parameters; the reduced line law has no independent tendon-fiber equilibrium solve.

Native path registration first maps OpenSim ground axes into the canonical display axes, then fits each bone group's axis-aligned bounds to the corresponding canonical bone bounds. Native path points retain source-body ownership; origin/insertion points are projected to full canonical bone vertices. The artifact records transforms, path-length scaling and projection distances. Wrapped geometry is frozen at the source rest configuration and transferred; moving-body wrapping is not recomputed on canonical surfaces. Endpoint projection does not establish a measured anatomical insertion. The remaining muscles use principal surface extrema, same-side bone proximity and a distinct-bone insertion prior. These can be wrong for complex, facial or intrinsic muscles, and carry broad uncertainty rather than confidence percentages.

Effective support connectors use `EA/L` stiffness and an explicitly assumed damping ratio. A geometric minimum-distance bone tree supplies connected skeletal support; tissue centers and actuator targets add physical force exchanges. These edges are generic supports, not falsely named anatomical joints. Vector dashpots include balancing moments so internal force and torque sums cancel. Rigid translations use linearly implicit attachment stiffness/damping with 2 ms substeps and a conservative polyline stiffness bound. Reference orientation constraints do no work because angular velocity is zero.

The mass ledger matches the shared 70.7713 kg generic profile, read from `data/derived/canonical/profile.json` rather than restated in the builder. That figure is composed over this specimen's own measured interior (`data/derived/interstitial-composition-prior-v1/ledger.json`) and supersedes the inherited BioGears StandardMale 77.1107029 kg constant, which needed a 1185.8 kg/m3 void and a negative Siri fat fraction to close. Tissue volume uses closed surface integration when available, otherwise an explicitly uncertain signed integral of the open surface and finally a bounds-based ellipsoid. Vascular wall volumes use an assumed 0.3 mm wall thickness; skin layers use the anatomy artifact's shell thickness. The skin parent, lymphatic structural graph and four cardiac cavities have zero material volume, each with 1e-6 kg numerical carrier inertia debited from the material allocation. Estimated material masses are uniformly normalized to the remaining body mass; the original density priors, unscaled total, multiplier and exclusions remain visible. This is an approximate mass allocation, not a measured partition or a proof that atlas surfaces pack without overlap.

Material parameters carry value, unit, basis, source IDs, dependency group and sensitivity range or explicitly unquantified uncertainty. Ranges are sensitivity envelopes, not probability intervals. Primary evidence includes ex-vivo microindentation showing kidney/liver moduli around 0.5–3 kPa and spatial heart moduli of 1–30 kPa ([Tissue micromechanics study](https://pubmed.ncbi.nlm.nih.gov/33176223/)); the adopted tissue-scale interpretation is an uncertain transfer. Human tibialis anterior tendon testing reports a tangent modulus of 1.2 GPa at maximal isometric load ([Maganaris and Paul](https://pubmed.ncbi.nlm.nih.gov/10562354/)); applying that to generic tendons is a prior. The 0.3 MPa generic muscle specific tension and its broad envelope are a synthesis assumption informed by the measurement approach in [Maganaris et al.](https://pubmed.ncbi.nlm.nih.gov/11181594/), not a reported measurement for every muscle. Other listed tissue moduli and support coefficients are explicitly assumed. [Rajagopal2016 source files](https://github.com/opensim-org/opensim-models/tree/master/Models/Rajagopal) provide the preserved native muscle parameters.

Verification checks constitutive objectivity, tetrahedral force-energy gradients, inversion rejection, internal force and torque cancellation, zero unloaded motion, source actuator bindings, isolated activation moving rigid bones and deforming its muscle, pressure/volume response, linear momentum conservation, mass total, time refinement and sustained one-second activation. The audit reports returned-state elastic/kinetic energy, active/external work, quasistatic affine boundary work, dissipation and discrete energy residual. Implicit numerical damping and work quadrature error are separate from biological uncertainty. These tests establish executable numerical behavior; they do not establish empirically calibrated anatomical insertions, material properties, joint behavior, local strain or whole-body physiological validity.
