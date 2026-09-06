# Supine geometric feasibility: retained surface/contact audit

The retained reference pose cannot be supported by **global lowering alone** on the current rigid bed within the skin foundation's declared compression domain. This is a geometric/contact-model limitation, not a reason to relax equilibrium or strain checks. It does not prove that every articulated pose is infeasible.

## Bounded evidence

Inspected `ihm/assembly/articulated.py`, `ihm/assembly/supine_contact.py`, `scripts/build_supine_surface_contact.py`, the contact manifest and its 503,073-byte quadrature, and `data/derived/supine-support-5ma720yd/initial_native.json`. No native launch, build, full skin mesh load, or download was performed.

The outer surface is `body-bp3d-FJ2810`, retained canonical geometry `data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz` (meters; `bodyparts3d-display-m`; SHA-256 `e724b62e6fc21764a0cf79c1874e96b9852755dc420711360a02d0b7afd93d05`). The builder retains 21,381 posterior source-X ray intersections, original face IDs and 5 mm projected-area cells. It transforms through the one global registration and attaches each point to the nearest named bone envelope. These attachments are inferred, not dissected skin compartments.

The contact manifest is `data/derived/supine-surface-contact-exmzq9pq/manifest.json`, SHA-256 `6c2216146edf30f30f4b1f3135d1817c453336444e73dc94d430cb45daa17458`; quadrature SHA-256 `a976dfe58cc61b5438e349f77d0a8d06018a2b299e39c937c931c9321125a9f4`.

| Quantity | Retained value |
|---|---:|
| Total native mass / weight | 77.6122 kg / 761.376 N |
| Full skin minimum source X | -0.268803 m |
| Sampled minimum clearance above that plane | 0.0461 mm |
| Pelvis minimum clearance | 41.766 mm |
| Heel (`calcn`) minimum clearances R / L | 160.594 / 159.980 mm |
| Skin total thickness / allowed compression | 6.6 / 3.3 mm |
| Maximum normal force after lowering to sampled compression limit | 12.4378 N |
| Contacting samples / projected area | 87 / 0.002175 m² |
| Contact sample source-Y interval | 1.29252–1.60752 m |
| Weight COM projected source Y | 0.921413 m |
| Resultant contact source Y | 1.58620 m |

The 12.4378 N calculation uses `plane = min(sample_X) + 0.0033` at zero velocity. Using the full-mesh minimum is 46 micrometers stricter; the sampled result is therefore a generous screen, not a certified full-surface strain bound. The first point was identified by the surface worker as occipital/head. Within the torso-owned samples, source-Y bands 1.45–1.8 m, 1.2–1.45 m, and 0.9–1.2 m have minimum gaps 0.046, 3.200 and 24.196 mm respectively. These are geometric bands, not independently segmented anatomy.

There are two separate obstructions at this unchanged pose:

- Available reaction is only 1.63% of weight. At the allowed stretch 0.5 the retained law gives at most 14.459 kPa. Even if every supporting cell reached this pressure, at least 0.052659 m² would be needed to carry weight, versus 0.002175 m² currently active.
- All active projected contact points lie cranial to the COM projection. With a stationary friction law giving zero tangential traction at zero velocity, a nonnegative sum of these normal forces cannot balance the gravity moment. Arbitrarily increasing stiffness would not fix this support-polygon failure. Internal muscle forces cannot supply a missing net external wrench.

## Degrees of freedom and model gaps

The native plant has 22 bodies and **no separate head or cervical bodies**. Canonical skull/neck contact points map onto the same rigid torso as thoracic points. Back coordinates move that whole torso relative to pelvis; they cannot change occiput-to-thorax geometry. Independent head positioning requires a real cervical/head articulation and anatomical attachment calibration. Splitting the current torso would require partitioning its existing mass/inertia and any muscle paths, preserving total mass; adding a donor head's mass on top would double count it.

Available pelvis orientation/translation, back, hip, knee, ankle, shoulder and elbow coordinates can change the support geometry. A bounded pose search should first use actual root orientation and source-permitted joint coordinates, checking a supine morphological envelope and full generalized residuals. It must not rotate the registration and gravity together and call that a pose correction. Global translation preserves interpoint height gaps; physical body rotation and articulation need not. The held passive elbow law also favors an isolated zero-torque angle near 1.56979 rad, so a straight-arm seed is not torque neutral (see `SUPINE_INITIALIZATION.md`). This angle is only a seed, not a prescribed equilibrium.

The current quadrature is a reference posterior envelope. Significant articulation can expose omitted surfaces or put different body regions over the same bed area. A broad pose result must be checked against rebuilt/current surface quadrature, including visibility/overlap and contact moments; the reference ray sample alone is insufficient.

## Actionable physical solution

1. **Separate skin deformation from bed deformation.** The current law has a stationary rigid plane and 6.6 mm skin columns. It contains no mattress travel. Introduce an explicit finite-thickness bed layer with its own displacement, constitutive response, energy/dissipation and bottom support. Couple skin and bed in series: total approach equals skin compression plus bed compression, with equal interface traction and separate domain bounds. Maintain equal/opposite world forces and moments; keep the existing skin bound. This represents a new physical system and cannot be implemented by relabeling skin penetration as mattress compression.
2. **Determine geometry before selecting material parameters.** For unchanged orientation, engaging pelvis after first head contact requires roughly 41.8 mm relative accommodation; heels require roughly 160 mm. These are first-contact geometric differences, not mattress thickness recommendations. A flat compliant bed could let early contacts sink, but finite thickness, bottoming-out, support footprint and moment balance must be checked. Actual articulated pose adjustment may reduce the needed travel. Physical contoured supports could also change the bed surface, but their geometry must be explicit, source-bound and included in the external wrench accounting; do not create invisible contact preload.
3. **Use a bounded pose/contact feasibility stage before another native solve.** Optimize only real available coordinates and explicit bed displacement variables; enforce skin and bed limits and nonnegative reactions. Check whether the projected COM can lie inside the support region and whether the admissible reaction capacity exceeds weight. This is a necessary screen, not a substitute for constrained generalized static equilibrium, passive muscle forces, or subsequent sustained forward acceptance.
4. **Calibrate head/torso registration before interpreting pressure anatomically.** One global COM/envelope fit has 62.356 mm RMS residual, larger than the skin layer thickness. Replacing its landmarks and correcting joint morphology must preserve the common world transform and exclusive native inertial ownership (see `ANATOMICAL_REGISTRATION_RESIDUALS.md`). A plausible visual support patch is not validation of the native pivot or tissue assignment.

No mattress/foam constitutive source or calibrated bed parameter set was found in the retained `data/sources` catalog or supine implementation documents. Existing E=3 kPa is explicitly a generic skin prior. The held `data/sources/sensorimotor/upperbody_sources.json` identifies `Neck3dof_point_constraint.osim`, but it is an OpenSim test fixture with placeholder credits/publications, not a validated donor cervical anatomy model. Neither source justifies a calibrated mattress or head/neck claim. The next source acquisition should specify the intended bed and obtain thickness, load–indentation/hysteresis data and boundary conditions; until then an explicit sensitivity model may be labeled an engineering model but cannot be reported as calibrated support.

No supported state, new parameter value, or changed acceptance threshold is established by this report.
