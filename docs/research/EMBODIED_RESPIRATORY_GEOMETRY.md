# Native-volume breathing within articulated geometry

`ihm/assembly/embodied_respiration.py` makes breathing visible within the continuing articulated body and supplies the force Jacobian of that same geometry. It uses the held `RespiratoryLoadPort.b` and canonical respiratory bindings. It introduces no independent gas state, passive chest/lung recoil, tissue inertia, dynamic mode integration or elastic energy. Native BioGears volume and recoil remain the owners.

## Geometry and articulation

The three displacement coordinates are `q = b * (V_native - V_initial)`, with volume in m³. The source normalization gives `a·b = 1`; therefore the linearized modal volume increment exactly matches native volume increment. The held generic stiffness only selects the static distribution b; its stiffness, damping and effective mass are not integrated again.

For reference point x, mechanical centroid c, proper rotation R, mechanical shape F and reference centroid c0, a binding maps the display point to:

`x_world = c + R F [B q + D(q,V) (x-c0)]`.

Ribs and sternum retain translation-only respiratory modes; posterior supports retain their zero source mode. The diaphragm uses the existing volume-preserving exponential shape. Lung anisotropic shape follows the source modal proportions and is normalized so `det(F D) = V_native/V_initial`. Thus a pre-existing mechanical shape can affect anisotropy but cannot add an independent lung volume. Uniform native volume ratio across the five canonical lobes is an explicit prior; their atlas mesh volumes are not equal to the native initial gas volume. No pleural/contact or rib-joint mechanics are claimed.

`geometry(current_lung_volume_ml, mechanical_entities, time_s)` returns composed entities and a `skin_field`. The skin field retains the existing `thoracic-smoothstep-v1` reference-coordinate formulation. The viewer must deform reference vertices first, then apply the current skin entity R/F/centroid transform exactly once. `coordinate_frame='canonical_reference_before_entity_transform'` and diagnostic `entity_transforms` clarify that order; the latter must not be applied a second time. Skin field deformation is not a new mass distribution.

## Work-conjugate force port

`project_load(forces, mechanical_entities, current_lung_volume_ml)` accepts current-world point forces `{id, point_m, force_n}`. It finds their material coordinates in the composed map and computes the derivative `J = ∂x_world/∂V` at fixed articulation. For deforming lungs and diaphragm, J includes both centroid translation and the off-center shape derivative. For skin, the material point is recovered through the current smooth field and its exact volume derivative is rotated through R F. A bounded local Newton solve performs only the point inversion.

The return pressure is `p_external = -Σ F·J`, so instantaneous external virtual work is `-p_external dV`. A force may move with the articulated body without silently retaining the original global direction of its volume Jacobian. `point_position(...)` exposes the forward material map for independently checking this relation.

This is an instantaneous virtual-work identity, **not exact finite-step external work**. Holding pressure over a native step approximates nonlinear lung/diaphragm shape loading; the caller must retain its temporal work quadrature and coupling latency. No recoil force is returned. Nonrespiratory force IDs are explicitly listed as ignored.

Native contacts may supply a resultant `{point_m, force_n, moment_nm}` about a body origin. Rigid rib/sternum bindings have zero virtual angular motion, so moments contribute zero in these translation-only modes. A resultant wrench does not determine distributed strain work in a deforming lung, diaphragm or skin field; those projections fail visibly and require resolved point tractions. The helper does not infer a canonical respiratory entity from an OpenSim contact body or invent a contact area.

## Verification and scope

`scripts/verify_embodied_respiration.py` uses tiny arrays only. Six tests cover source payload posterior supports, lung volume ratios, volume-preserving diaphragm shape, rigid rotation and virtual-work conjugacy, off-center lung/diaphragm/skin Jacobians against centered finite differences, wrench scope, zero-increment composition and invalid volume rejection.

The actual held payload contains 45 bindings: 24 ribs, 12 posterior supports, three sternum pieces, one diaphragm and five lung lobes. For an illustrative native volume change from 3,000 to 3,300 mL, computed q is approximately `(0.7476, 1.0443, 5.3839)` mm and each lung shape determinant is 1.1. This verifies the implemented display constraint, not a biological measurement or a native breathing experiment. No native process, build, browser or broad geometry job was used.

Unresolved boundaries include native gas-to-atlas shape calibration, regional lobe partition, diaphragm/abdominal contact, rib rotation, local bed or garment traction mapping, and distribution of body gravity across respiratory modes. The source small-deformation guard remains active. This module supplies consistent reduced geometry and load projection; it does not close those whole-body mechanical questions.
