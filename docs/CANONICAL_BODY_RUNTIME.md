# One generic body

IHM's primary materialization is `ihm-body`: one anatomical identity set and coordinate frame, a shared generic demographic profile, source-informed mechanics and an IBM-derived regional neural model driven by native physiological state. The source viewers remain available for evidence inspection.

This is an executable research assembly, **not a complete or validated human digital twin**. Calibrated source parameters coexist with registered geometry and explicit engineering priors. Global source participation does not make every cross-system interaction measured.

## Build and execute

```bash
.venv/bin/python scripts/build_canonical_anatomy.py --append
.venv/bin/python scripts/build_body_profile.py
.venv/bin/python scripts/build_body_brain.py
.venv/bin/python scripts/build_body_mechanics.py
.venv/bin/python scripts/build_canonical_body.py
```

The existing native adapter must include `body_compartments.csv` telemetry. Recompile it with `scripts/build_native_adapter.py` against the acquired engine. Run the native model using the generated `IHMGenericMale` patient and the existing `NativeConfig` / `run_native` interface, then materialize that recorded execution:

```bash
.venv/bin/python scripts/build_canonical_body.py \
  --native-directory data/derived/canonical/native_baseline_v1
```

```python
from ihm.human import ImplicitHuman
body = ImplicitHuman.open('.').materialize('body')
body.describe()
body.certainty('body-skin-dermis')
body.simulate('data/derived/canonical/native_baseline_v1', 'artifacts/body.json')
```

`/api/body` describes the model; `/api/body/certainty?entity=...` exposes provenance and separate uncertainty dimensions; `/api/body/trajectory` serves its recorded shared-clock state; `/api/body/spectra` exposes finite temporal transforms. `POST /api/body/scenarios` accepts native scenario configuration while requiring the canonical patient. Completed canonical jobs provide `/api/body/trajectory?run=...`.

## State ownership and transfer

BioGears owns the coupled physiological balances, including blood, heart, respiration, renal and endocrine state, thermal state, and aggregate tissue/lymph fluid. The read-only C++ adapter records actual named compartment volumes and pressures on the native clock; overlapping aggregate compartment volumes must not be summed as disjoint storage.

At each native sample interval, mechanics receives prescribed volume ratios and the IBM-derived neural model receives mean arterial pressure, oxygen saturation and core temperature. Each has its own documented numerical substeps and must arrive at the same elapsed time. Native sampling is not replaced by a rendered heartbeat sine wave.

The heart transfer projects a lumped chamber-side volume ratio onto the ventricular **cavity** surface. This explicitly does not establish atrial/ventricular partition, myocardial contraction or pressure-volume calibration. The five registered lung lobes receive ipsilateral gas-volume ratios through uniform isotropic expansion. Pleural contact, tissue/gas partition and regional ventilation remain unresolved. These maps are directly synthesized transfers and are labeled accordingly.

Mechanical reactions and neural autonomic commands are reported, but are not fed back into BioGears. Its existing nervous system remains the physiological controller. Adding another uncalibrated controller would count control twice. BETSE/wound, JOS3, CSF, reproductive and vascular CFD remain identifiable regional or alternative materializations. There is no claim that a registered vascular surface already has a 3D flow solution, or that the undirected lymph graph establishes calibrated drainage.

The neural initial state is an explicit source-law resting initialization, not a fitted equilibrium for the current physiology. Early neural transients must be interpreted accordingly.

## Precision and certainty

Every anatomical entity carries source identity, dependency group, assumptions and three separate dimensions:

* Biological variability and subject calibration. Unknown variability remains unquantified.
* Registration or synthesis uncertainty. Held-out registration errors are geometric diagnostics; engineering ranges are assumptions, not confidence intervals.
* Numerical/display approximation. Full source geometry is preserved; GPU float32 conversion and sparse transform export have explicit numerical bounds.

`ihm.assembly.certainty.fuse` accepts scalar estimates only with explicit units, justified positive variances and dependency groups. It applies equal-weight covariance intersection within related-source groups and independent Gaussian fusion between groups. Distinct group independence is a caller assertion. Derived atlases do not create extra independent subjects. This utility does not manufacture variances for the unquantified anatomy.

`body.json` pins anatomical, profile, mechanics and brain artifacts and the executable body/brain/mechanics/certainty/temporal Python modules by SHA-256. Materialization rejects mismatched identities, source hashes and patient identity, or a mismatch between physiology and compartment clocks. The API rejects trajectories and spectra based on stale canonical artifacts or changed native receipts. Queued scenarios pin their submitted source hashes and reject intervening rebuilds. Input byte snapshots provide the hashes recorded in outputs, preventing an input rewrite during integration from changing the claimed evidence. Each trajectory records its native input hashes and transfer assumptions.

## Display and temporal reduction

Explicit transforms use `c + translation + R F (vertex − c)` about the recorded reference bounding-box center. This center is not asserted to be the physical center of mass. Sparse trajectory entries omit only transforms whose translation components and matrix deviations are at most `1e-9` (meters for translation, dimensionless for matrices); absence means identity for that frame. The full final mechanical state is retained separately. Playback is a finite recorded interval, not extrapolation.

The canonical spectral export uses the full native samples, independently of display frame reduction. It records Hann Welch densities and finite Laplace transforms at `s = sigma + 2πif`, with the sample mean explicitly removed, physical time units, input hashes and finite integration interval. These finite-window descriptors are not inferred system poles or causal transfer functions. Predictive and mechanistic interpretations still require separate experiments and validation.

## Remaining scientific limits

The mechanics reduction constrains orientation and does not solve articulated joints, gait, volumetric contact, vessel-wall fluid interaction or a gravity-loaded bed contact equilibrium. Rest coordinates are stress-free reference priors. Some source surfaces overlap or remain open; the explicit mass allocation is normalized to the generic profile and does not constitute measured organ masses. Many stiffnesses, attachments, skin depths and transfer coefficients remain priors. Whole-body membrane-voltage coupling, microcirculation, tissue-resolved lymph drainage and complete cross-system uncertainty propagation remain research work.

See [canonical anatomy](CANONICAL_ANATOMY.md), [mechanics](BODY_MECHANICS.md), and [IBM brain integration](BODY_BRAIN.md) for source details and verification scope.

The first generic baseline was executed before the companion compartment receipt was added. Its explicit post-execution audit preserves `summary.before-compartment-receipt.json`, verifies the original physiological CSV hash and the companion clock, and records that the original execution-time compartment hash did not exist. Subsequent native executions produce both receipts automatically.
