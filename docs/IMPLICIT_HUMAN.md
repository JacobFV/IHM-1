# Source-backed implicit human

`ImplicitHuman` binds acquired evidence by path and SHA-256. It exposes actual
native scalar fields and materializes predictors according to the evidence that
can support them. Source families remain separate when no measured registration,
shared subject identity, or cross-family covariance exists.

```python
from ihm import ImplicitHuman
human = ImplicitHuman.open('.')
human.save('artifacts/integrated-human.json')
print(human.describe())
print(human.fields()[0])

population = human.materialize('population')
# Observation tuple: value, independent measurement variance, exact source unit.
component, unit = population.components[0], population.units[0]
posterior = population.condition({component: (population.mean[0], 1.0, unit)})

field = human.materialize('skin-field')
print(field.predict('18-29', 'female', 'LA-1'))

circuit = human.materialize('skin-lymph')
print(circuit.step(.02))
print(circuit.response(1e-5 + 1j, 'Aorta1'))

run_id = human.describe()['temporal_runs'][0]
temporal = human.materialize('temporal', run_id=run_id)
forecast = temporal.forecast(temporal.center, 100)

native = human.materialize('native', seconds=60, sample_hz=50)
# Runs the actual locally compiled source engine; output directory must be fresh.
# native.run('artifacts/my-native-run')
```

| Materialization | Evidence and resulting behavior | Scope |
|---|---|---|
| `population` | 23 survey-weighted measured states, full covariance, Bayesian conditioning | Concurrent population prediction; no dynamics |
| `skin-field` | Six fitted phenotypic coefficients, subject-held-out human observations and clustered uncertainty | Lateral wound field above epidermis; no membrane-voltage inference |
| `skin-lymph` | Actual native R/C/source/gate snapshot, conservative descriptor integration and causal Laplace response | Frozen operating point with prescribed arterial/venous boundaries |
| `temporal` | Native or measured waveform SVD/reduced dynamics, chronological holdout | Autonomous empirical forecast; no causal intervention interpretation |
| `native` | Original BioGears coupled source equations and patient initialization | Whole-engine simulation, explicit action timeline; independent calibration remains unestablished |

`save()` verifies all evidence again. `load()` rejects changed hashes, unknown
asset keys and substituted paths. Existing in-memory materializations retain
the evidence read at their creation. Open a new substrate after rebuilding
an evidence atlas. Manifests do not bundle the several gigabytes of raw assets.

The original `body()` registry retains the IBM-style illustrative modeling
scaffold and its dependency closure. `ImplicitHuman` is the entry point for
source-backed materializations. In particular, the original `build whole-body`
command has illustrative dynamics even when its initial belief is measured.

## Integration contract

Native cardiovascular, blood chemistry, respiration, renal, endocrine, energy,
gastrointestinal, hepatic, nervous, tissue and other source systems exchange
state inside their engine. The exported circuit graph retains their shared
nodes and compartments without double counting alternative views. A selected
skin–lymph circuit can be compiled to an explicit predictor from this topology.

Cross-source exchange operators implement atomic conserved extensive fluxes,
positive implicit ion transport and physical-frame mapping. They require explicit
coefficients and registrations. Merely giving an atlas mesh and a native
compartment the same organ name does not connect their physics.

The posture hydrostatic artifact samples a hypothetical static arterial column
on physical atlas coordinates using source-native blood density and arterial
reference pressure. It distinguishes supine and upright gravity. It does not
predict organ perfusion, venous pooling, baroreflexes or a native posture change.

## Additional source models

`materialize('opensim')` returns a native mechanics runner with the independently
checked wrap-cache correction. Its explicit activation, joint perturbation,
engine variant and output hashes are recorded. `materialize('reproductive')`
runs the source gonadotropin model; prescribed ovarian inputs do not constitute
an autonomous menstrual cycle. `materialize('csf')` runs the published
intracranial model with source parameters. A `PressureInput.from_native_map(path)`
can drive it from a unit-checked native MAP recording; there is no ICP feedback.
The cellular recording is an acquired `bioelectric` asset, not an invented
whole-skin predictor. See [its native field view](BIOELECTRIC_VIEW.md).

Population conditioning uses correlation-scaled covariance factors and
noise-whitened SVD. Correlation eigenvalues within the covariance validator's
floating-point tolerance are treated as numerical null directions. This retains
small positive assay noise even for exactly aliased variables, and keeps unit
scaling from deciding the effective rank.


`materialize('thermal', profile='supine_blanket', seconds=3600, dt=30).run()`
executes the pinned JOS-3 solver with the measured total-insulation boundary.
The other source profiles are `lying_default`, `supine_mattress`, and
`supine_duvet`. The returned 85-node trajectory carries coefficients, heat
balance audits and source identity. See [thermal source and refinement](THERMAL_MODEL.md).
The workbench's saved cases use 0.9375-second steps; the default Python step
is an explicit computational choice, not a claim of timestep-independent accuracy.

The source-target audit retains all 475 upstream literature records per run,
including unmapped or ambiguous entries. Exact-quantity comparisons expose
point deviations without inventing a clinical pass tolerance. See
[native target audit](NATIVE_TARGET_AUDIT.md). The searchable catalog also
retains thermal capacities, conductances and boundary values with their units.
