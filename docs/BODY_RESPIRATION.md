# Respiratory mechanics and visible chest motion

The canonical respiratory substructure converts the **actual native left/right lung gas-volume change** into connected chest movement. It has three mechanical displacement coordinates: lateral rib expansion, anterior chest expansion, and caudal diaphragm excursion. Ribs, sternum, diaphragm and lung centroids become prescribed boundaries of `BodyMechanics`; attached unconstrained tissues move because the structural solver transmits their forces. A regional skin field makes the solved chest expansion visible on the original skin mesh.

This is an executable, generic reduced model. It is **not** a volumetric lung/chest FEM, a pleural contact solution, a patient-calibrated model, or bidirectional native respiratory physiology.

## Mechanical equations

For displacement `q` in metres and cavity volume Jacobian `A` in square metres:

```
M q'' + D q' + K q = f + A p
A · q = native_current_gas_volume - native_reference_gas_volume
```

`M` is an effective modal mass, `D` is positive viscous damping, and `K` is a positive symmetric coupled stiffness. The reference cavity is an elliptical cylinder derived from the bounds of the five canonical lung lobes. For lateral radius `a`, AP radius `b`, and height `h`, `A = [πbh, πah/2, πab]`. The constraint conserves the **linearized reduced cavity volume**, not an integral of the rendered mesh's enclosed volume. The unknown pressure is the constraint reaction incremental to reference. It is not a measured pleural pressure.

The pressure alternative applies `A p` without the native volume constraint. It exists for explicit model experiments and verification; the canonical native replay uses volume mode. Prescribing both pressure and volume is rejected.

Backward Euler solves the three coordinates with the volume constraint at every substep. Within each native interval, the supplied volume boundary is linearly interpolated. The audit reports elastic and kinetic energy, pressure/external work, viscous dissipation, and the positive numerical dissipation of backward Euler. Their discrete balance is checked to floating-point precision.

## Parameter evidence and uncertainty

The model uses 0.2 L/cmH₂O local chest-wall compliance, an exploratory range of 0.1–0.3 L/cmH₂O, volume fractions `[0.18, 0.17, 0.65]`, effective masses `[1.5, 1, 0.5]` kg, and damping `D = 0.08 K`. These are **uncalibrated generic/modeling priors**. The stated ranges are exploration brackets, not confidence intervals or estimates extracted from a cohort. Mode fractions, reduced masses and damping are not measured for this anatomical scaffold.

[Marconi and De Lazzari's respiratory mechanics model](https://pmc.ncbi.nlm.nih.gov/articles/PMC7239037/) supplies the pressure–volume and additive recoil relationships motivating the reduction. Its nonlinear lung/chest model is not reimplemented here. [Aliverti and colleagues' walking study](https://pubmed.ncbi.nlm.nih.gov/10484561/) supports regional chest/abdominal volume accounting. Its measured regional partition is not transferred to this body. [Boussuges and colleagues' ultrasound study](https://journal.chestnet.org/article/S0012-3692%2809%2960127-6/abstract) supplies an independent order-of-magnitude comparison for quiet-breathing diaphragm motion; agreement is not validation.

At the implemented priors, a sustained 0.5-L increment yields about 9 mm diaphragm descent. The physiology determines the waveform; there is no added sine wave or visually amplified breathing amplitude.

## Body and renderer contract

Build with `.venv/bin/python scripts/build_body_respiration.py`. It writes `data/derived/canonical/respiration.json`, retaining the anatomy SHA-256, sources, assumptions and bindings. The source anatomical frame is X left, Y superior, Z anterior, in metres.

```python
respiration = BodyRespiration.from_dict(respiration_asset)
state = respiration.step(dt, {
    "volume_change_m3": native_delta_left_plus_right_mL * 1e-6,
    "lung_volume_ratios": {lobe_id: current_side_mL / reference_side_mL},
})
mechanical = mechanics.step(dt, {
    **other_mechanical_drivers,
    **state["mechanics_drivers"],
})
```

Remove lung IDs from `other_mechanical_drivers['volume_ratios']`: each tissue accepts one volume, pressure or affine deformation boundary. The returned `mechanics_drivers` contains:

- `prescribed_translations_m`: target translations of 24 ribs, three sternal structures, diaphragm, supplied lung lobes, and zero translations of 12 thoracic vertebral supports.
- `deformation_gradients`: anisotropic lung tensors whose determinants equal the supplied native lobe ratios, plus an isochoric diaphragm tensor. BodyMechanics owns their constitutive energy once.

BodyMechanics imposes translations by Dirichlet elimination of the implicit velocity equations. Targets are interpolated through its substeps. Free neighbors receive the attachment forces; the prescribed reaction forces and accumulated boundary work are returned. Posterior supports make the reference-pose constraint explicit and prevent free whole-body translation. Rigid orientations remain constrained: ribs translate, but anatomically articulated rib rotation is unresolved.

Use **the final BodyMechanics entity transforms** in `frames.entities`; do not overlay the respiratory transforms afterward. `state['entities']` describes the respiratory boundary only. The separate `state['skin_field']` applies only to its listed skin entity IDs. Apply the regional displacement to each reference vertex before its existing body transform, never cumulatively to last frame's vertices. `deform_skin` in `ihm/assembly/respiration.py` is the reference formula. Smooth tapers bound the field laterally and vertically; head and feet are unchanged by the respiratory field. The shader and this Python implementation must agree.

## Interaction scope

Respiratory boundaries exchange forces with the entire existing connected attachment graph. A local generalized force can redistribute the three chest coordinates while the native volume constraint remains exact. Structural boundary reactions are audited; they are not currently iterated back into the respiratory mode partition or the native engine. The chest modal energy and the structural energy are separate audits of that partitioned system, not an asserted globally calibrated energy model.

Unresolved mechanics include rib joint articulation, local rib bending, pleural sliding/contact, lung fissure contact, tissue self-contact, skin's own volumetric constitutive solve, abdominal organ contact and bulging, and interactions between all anatomical surfaces. The linear cavity model rejects excursions beyond its explicit small-deformation range. These limits must remain visible in model metadata and user-facing claims.

## Verification

Run `.venv/bin/python scripts/verify_body_respiration.py` for zero-input rest; sustained pressure/compliance response; volume conservation and force balance; rib/chest/diaphragm direction; isochoric diaphragm and native lung determinants; local force redistribution; regional skin motion; timestep refinement; origin invariance; transactional rejection; and real attachment-neighbor response to prescribed chest motion.

Run `.venv/bin/python scripts/verify_body_respiration.py --native` for all 1,500 intervals of the 30-second 50-Hz native baseline. It checks shared clocks, bounded whole-body translations, the exact reduced volume constraint and prescribed-position residuals, and records the input hash and elapsed computation time. This is trajectory replay verification, not biological validation.
