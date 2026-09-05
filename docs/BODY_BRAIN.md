# Brain in the canonical body

The canonical body reuses executable code and anatomical assets already held by IBM-1. Its 80 regional populations reference existing BodyParts3D tissue entities, including 68 Desikan–Killiany cortical regions and 12 coarse subcortical populations. This is a generic regional neural model with explicit physiological inputs; it is not a subject-specific whole-brain simulation.

Build and verify from the IHM-1 repository:

```bash
.venv/bin/python scripts/build_body_brain.py --ibm-root ~/Documents/IBM-1
.venv/bin/python scripts/verify_body_brain.py
```

The builder resolves the actual MNE sample cache from IBM's `data/sources/mne-sample/raw/.location.yaml`. It reads fsaverage's full left/right pial surfaces and `aparc` annotations. `--asset-root` can explicitly select the fsaverage subject directory; `--asset-python` selects a Python containing nibabel for this ingestion step. The default is IBM's existing virtual environment. It does not install anything, alter the IBM workspace, set environment variables, or modify Python's import path to import IBM. Runtime needs only IHM's existing NumPy dependency and the preserved local source snapshots.

The output files are:

- `data/derived/canonical/brain.json`: regional tissue identities, connections, registration, source hashes, parameters, ports and assumptions.
- `brain-sources/`: all eleven original source files, including complete IBM neural/anatomy/topology code and the four full surface/annotation files. Copies are byte-preserving; source and preserved paths, sizes and SHA-256 hashes are recorded.
- `brain-registered-surfaces.npz`: all original pial vertices and triangles under the recorded affine, with original annotation indices. This is registration evidence; it is not another competing cortex in the body's organ inventory.
- `brain-state.json`: an executed baseline followed by a synthetic low-perfusion/low-oxygen perturbation, with complete regional states. These are model outputs, not recorded measurements or native-body feedback.

Generated outputs follow the repository's existing `data/derived` exclusion. Rebuild them on a new checkout with the explicitly configured source locations.

## Source model and runtime

`ihm.assembly.brain.BodyBrain` parses the complete, preserved `ibm/processes/neural.py` and compiles its exact `_sigmoid`, `wilson_cowan_excitatory` and `shunting_inhibition_rate` function ASTs. It first verifies every preserved source hash. This isolates IBM's concrete nonlinear population equations from its unrelated registry/NN imports without rewriting those equations or modifying IBM. The complete original file remains available to inspect its context and parameter provenance.

The imported excitatory law evolves membrane potential, firing rate and adaptation. IBM's default membrane time constant is 15 ms, rate relaxation is 5 ms, adaptation time is 0.5 s, rest potential is −65 mV and rate ceiling is 100 Hz. The defaults retain IBM's literature/weak-prior identity; they are not fitted participant measurements. The imported shunting law supplies activity-dependent inhibitory feedback.

IHM provides the regional wiring and integrates the coupled ODE with RK4 at at most 1 ms by default. Eight nearest regional neighbors supply 640 directed edges, weighted by the 40 mm exponential geometric prior described in IBM's association topology and normalized by incoming sum. These edges are inferred connections, not measured tractography. Synaptic conductance here is a regional instantaneous drive with a fixed AMPA/NMDA split; the full IBM receptor dynamics, spectral runtime, axonal conduction delays, tract measurements and EEG/BOLD observation model are not included.

```python
import json
from pathlib import Path
from ihm.assembly.brain import BodyBrain

brain = BodyBrain.from_dict(json.loads(
    Path('data/derived/canonical/brain.json').read_text()))
state = brain.step(0.1, {
    'mean_arterial_pressure_mmHg': 90.0,
    'oxygen_saturation': 0.98,
    'core_temperature_C': 37.0,
})
print(state['regional_state']['activity_hz'])
print(state['autonomic_commands'])
```

`from_dict(..., root=...)` roots preserved paths explicitly when the complete model directory has been relocated. Each instance has independent neural state. `step` accepts seconds in `[0, 60]`; a configurable `max_step_s` must lie in `(0, 0.002]`. Nonfinite values, invalid saturation units, implausible input ranges and source hash mismatches fail explicitly. Missing physiology inputs are marked `assumed_baseline`, while supplied values are marked `caller_supplied`; a caller must separately preserve whether supplied values were live, recorded or synthetic.

MAP and oxygen saturation scale available synaptic drive; temperature scales kinetics using an explicit Q10 prior. The medulla population's activity and physiological stress produce bounded sympathetic and parasympathetic commands. These commands are returned with `applied_to_body: false`; this component never silently feeds them into a native engine. Bidirectional physiological calibration is not claimed. The canonical body's owner can connect these explicit ports while preserving the scope of that coupling.

## Anatomical registration and uncertainty

The source is fsaverage surface tkRAS millimeters: right, anterior, superior. The destination is the existing BP3D display frame in meters: left, superior, anterior. Twenty corresponding named gyral vertex centroids fit an axis-preserving affine scale and translation after the explicit axis conversion. The build records the 4×4 transform, every landmark pair, individual residuals and the original BP display transform.

This build's landmark RMS is **15.65 mm**. That is an in-sample centroid mismatch between two different templates, not a held-out anatomical error bound. Named matching regions inherit homologous BP gyral support identities; finer unmatched DK regions map to the nearest appropriate-hemisphere coarse BP cortical support, with the support distance exposed. The transferred region's position is an affine source centroid and is not guaranteed to be on the target gyral surface. The body's BP tissue geometry remains authoritative, and the transform introduces no duplicate organ.

The twelve subcortical populations use actual BP medulla, pons, midbrain, hypothalamus, thalamic and cerebellar tissue supports. Their dynamics use coarse cortical-law proxies. Nucleus-specific cell populations and calibrated autonomic center dynamics remain explicit omissions.

The available BP nerve inventory supplies orbital/eye nerve interfaces. Optic nerve interfaces point to the coarse thalamic relay, oculomotor/trochlear interfaces to midbrain, and other orbital sensory interfaces to pons. These are anatomical interface priors, with explicit canonical entity IDs, rather than reconstructed axonal routes or simulated peripheral conduction. This source inventory does not provide systemic peripheral/autonomic nerve trajectories or spinal-cord geometry.

Verification checks source-byte identity, 68 actual cortical parcels, canonical tissue links, frame handedness and scale, independent runtime state, the hand-evaluated IBM rate law, finite responses, suppression under reduced oxygen/perfusion, autonomic direction, rejected invalid input and agreement between 1 ms and 0.5 ms integrations. Passing these checks establishes source fidelity and numerical behavior within this reduced model; it does not establish empirical physiological calibration.
