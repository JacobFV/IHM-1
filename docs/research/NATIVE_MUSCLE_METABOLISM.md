# Source-native muscle metabolic power

`scripts/native_muscle_metabolism.h` wraps the held OpenSim `Umberger2010MuscleMetabolicsProbe`. It adds a measurement component to the native mechanical model, not a second muscle actuator, tissue mass or inertia. It covers activation/maintenance heat, shortening/lengthening heat and signed active-fiber mechanical work through the actual native implementation. Positive mechanical work alone omits isometric cost and is not interchangeable with the extra metabolic demand consumed by BioGears.

## Source and version qualification

The installed OpenSim header and held `.cpp` implementation, plus its upstream test fixtures, were inspected directly. The primary [Uchida et al. 2016 paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0150378) describes an updated muscle energetic model and the limitations of using positive fiber power as a metabolic surrogate. The OpenSim header explicitly says its equations may differ from the publications; the held executable implementation is the exact constitutive authority here. This wrapper does not reimplement the equations from prose. The source header also contains an inconsistent sign annotation in the shortening section; executable code uses negative fiber velocity for shortening and positive velocity for lengthening.

Frozen evidence hashes are in `data/sources/native-muscle-metabolism.json`. No new native build or numerical job was run for this wrapper during the resource incident. API/source review is complete; compilation and actual muscle-state acceptance remain gates for the serialized mechanical build owner.

## Explicit source defaults and ownership

The wrapper selects activation/maintenance, shortening and mechanical work, while switching **whole-body basal power off** because BioGears owns basal metabolic demand. It retains these held OpenSim defaults:

- Aerobic scaling 1.5, effort scaling 1, Bhargava recruitment enabled.
- Signed negative mechanical work included; source protection against negative total metabolic power enabled.
- Minimum muscle heat 1 W/kg retained.
- Slow-twitch composition 0.5 per muscle, density 1,059.7 kg/m³, specific tension 0.25 MPa. These are generic source defaults, not measurements of the canonical individual or each named muscle.

Optional per-muscle slow-twitch fraction and provided analysis mass must be explicit. Otherwise the native probe computes analysis mass as `Fmax / specific_tension * density * optimal_fiber_length`. The unmodified held 80-muscle example implies 40.5616268688 kg by this formula. This is a metabolic analysis quantity, **not additional body mass**; the body segments already own mechanical mass/inertia. Uniform segment inertia rescaling does not establish consistent muscle-specific composition or metabolic mass calibration.

Because the retained source minimum heat is nonzero at rest, merely disabling the probe's whole-body basal term does not make its raw total an incremental demand above BioGears basal metabolism. The caller must preserve raw modeled total, choose and retain an explicit reference muscle state, and report the signed difference. Adding raw muscle total to BioGears BMR double-counts an unspecified resting component. A negative difference cannot silently become zero or be supplied to an interface that only accepts nonnegative exercise demand. A fixed initial-state reference is a coupling assumption, not a clinical calibration; activation, isometric support and shortening changes must all remain visible relative to it.

## Interface

After all source and added muscles are installed, and before final connections/`initSystem`:

```cpp
IHMNativeMuscleMetabolism metabolism(model);
```

For a valid native state:

```cpp
const auto rate = metabolism.sample(state); // realizes Dynamics
// rate.total_muscle_metabolic_w
// rate.active_fiber_work_w (signed)
// rate.muscle_heat_w
// rate.analysis_mass_kg
// rate.muscles[name]: type, power/work/heat, mass, slow-twitch fraction
```

The probe returns total, basal and per-muscle power. The helper requires basal to be zero, all declared muscles to be covered, finite outputs and agreement of the per-muscle sum with total. Heat is total minus the source-consistent signed active-fiber work; it is not a separately evaluated probe with altered clamping interactions. The actual supported classes are `Millard2012EquilibriumMuscle` and `Thelen2003Muscle`; other classes fail visibly until reviewed. Thelen support is based on the source API, not a completed numerical fixture.

The adapter owns energy accumulation. If using endpoint trapezoidal power quadrature, record that numerical approximation and checkpoint cumulative metabolic energy alongside state, excitation, loads and reference power. The probe operation is `value`, so it introduces no hidden energy integrator requiring independent rollback. Keep mechanical work and metabolic energy separate in telemetry; do not add active work to the metabolic total again.

Required next native gates are zero-activation/reference behavior, positive isometric activation heat with zero fiber work, shortening and lengthening branches, per-muscle/total agreement, additional Thelen coverage, finite state, matched reference subtraction and checkpoint replay. No passing claim for those gates is made by this source-only implementation.
