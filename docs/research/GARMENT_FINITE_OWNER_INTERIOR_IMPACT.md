# Finite-owner garment interior impact step

`ihm/assembly/garment_vertex_face_step.py` advances explicitly supplied finite cloth and surface owners through one supported interior impact. It is separate from the prescribed skin boundary in `GarmentFeedback`: no finite skin masses are invented or copied from native tissue. The native garment path retains its strict guard until an accepted native impulse interface is available.

The operator accepts Cloth-compatible owners with position, velocity, positive nodal mass, a common clock, an elastic timestep bound and a pure `elastic_forces` method. The surface owner supplies triangles. The caller explicitly supplies candidate cloth node and surface face IDs. Those candidates are bounded (default maximum 4096 pairs); they are not a claim of complete whole-body collision coverage. Duplicate object/material identities, self-owner contact, fixed support nodes, malformed inputs, nonfinite state and clock disagreement reject before any state commits.

## Supported step

1. Apply one explicit elastic/external/gravity kick to copied states.
2. Find the first interior node/triangle event using the corrected linear swept kernel. Multiple candidate impacts, edge/vertex, degenerate, tangential, coplanar or initial impact cases reject.
3. Drift all owners to the event. Call the existing `resolve_node_triangle_contact` in `contact_dynamics.py` at its verified coincident interior point. The impulse uses inverse effective mass `1/m_node + Σ(beta_i²/m_i)`, zero normal restitution, and the existing static/kinetic Coulomb law. The surface receives `−beta_i × impulse`; the cloth receives the opposite total impulse at the same spatial point.
4. Reconsider every declared candidate for the remaining interval. The impacted pair currently supports a translating planar target with an interior convex linear contact path. A deforming postimpact face, edge exit, remaining inward normal velocity or another event rejects. This bounded translating-face case allows zero-restitution contact to continue through the rest of the step without manufacturing a separating velocity or skipping a contact check.
5. Drift the validated remainder. Reevaluate elastic energy and all ledgers; commit every owner and the common clock only after acceptance.

This is one kick plus event-subdivided drift, not continuous force integration after collision. The remaining drift is explicitly ballistic. General deforming-face sustained contact, multiple impacts, self-contact, thickness, initial containment and full candidate discovery remain unimplemented. Inputs outside the declared supported event cases fail closed; no full CCD guarantee follows from a clear candidate set.

## Work and acceptance

Normal restitution is exactly zero. Reports separate normal impact loss and friction loss, total paired impulse, momentum residual, angular impulse residual, external work, gravity work and the numerical elastic-integration energy defect. External and gravity work refer to their discrete kick, using the midpoint of pre/post-kick velocity. They are not inferred from a contact-altered endpoint displacement. No integration defect is silently credited as heat.

The caller must supply a finite nonnegative `max_energy_defect_j`. A larger defect rejects before mutation. Momentum and angular residuals use explicit absolute numerical thresholds of 1e-10 N·s and 1e-10 N·m·s. These are numerical acceptance thresholds for this bounded operator, not calibrated physical tolerances. Pure ballistic collision fixtures close their energy/momentum ledgers to floating-point precision.

## Small fixture acceptance and remaining native work

`scripts/verify_garment_vertex_face_step.py` exercises ten cases with tiny all-mobile finite owners: first event at 0.1 s and validated remainder to 0.2 s, exact paired barycentric impulse, normal/friction loss, rigid and Galilean frame invariance, multiple and edge event rejection, deforming postimpact and newly encountered event rejection, energy-budget rollback, candidate limits, duplicate-owner rejection, and external kick work. Existing garment guard and feedback tests remain separate regressions.

For the normal centroid impact of a 1 kg node at −1 m/s into three stationary 1 kg triangle nodes, effective inverse mass is 4/3 kg⁻¹. Impulse is 0.75 N·s, final common normal velocity is −0.25 m/s, and dissipated normal energy is 0.375 J. The remainder moves the common plane by −0.025 m. These values arise from the supplied masses; they are not human skin material parameters.

A real native body uses articulated inertia, not these independent skin nodal masses. Native integration needs an actual contact-point inverse effective mass and impulse/velocity response from its owning multibody state, synchronized event-time checkpoint/advance, work/energy readout and rollback. Barycentric skin mapping can transmit a wrench to that owner but cannot substitute arbitrary nodal mass for articulated inertia. No native process, full-body trajectory or browser display was run or accepted in this increment.

Independent review caught a local-versus-owner event clock error before commit. A regression starting both owners at 10 s failed when the event reported 0.1 s; the corrected receipt reports absolute `event.time_s=10.1` and explicit `elapsed_time_s=0.1`, with the committed owner clock at 10.2 s.

Independent review passed the corrected ten-test suite and additional common-gravity and static/sliding friction cases: impact loss and kick work close the energy ledger, with no further concrete defect found in that bounded review. Every positive postimpact remainder is checked; no small-duration remainder bypass is used.
