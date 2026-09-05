# A source-grounded digital human: revised design

The original illustrative compartment scaffold is not the intended digital human.
The acquired data makes the replacement architecture concrete: independent 3D
anatomical supports, source-specific dynamical models, physiological state and
parameter distributions, and explicit coupling interfaces. The initial posture is
supine at rest. This document describes how the acquired assets constrain that
design, and separates working imports from missing physical integration.

## What now exists locally

| Domain | Actual imported assets | What they establish |
|---|---|---|
| Skeleton and soft-tissue anatomy | BodyParts3D OBJ atlas, ontology tables | Reference surfaces and anatomical membership; not subject registration |
| Musculoskeletal mechanics | OpenSim models, geometry, joints, muscle paths, wrapping and parameters | A coherent musculoskeletal model family with local body frames; variants overlap |
| Vascular flow | Three human VMR meshes, flow/BC files and cerebral 3D pressure/velocity time states | Case-specific discretizations and archived CFD outputs; simulated flow differs from in-vivo measurement |
| Whole-system physiology | HumMod, Physiomodel, Physiolibrary, BioGears source/model/validation assets | Existing cross-system equations and author-defined parameterizations; not a license to mix isolated constants |
| Cellular bioelectricity | BETSE source, published example configurations, locally executed full solver | Ion/electrodiffusive process implementation; a generic tissue run is not human wound calibration |
| Human wound bioelectricity | 40-person clinical table, 155 numeric observations and 16 summary targets | Layer-specific lateral field measurements; not membrane voltage or identified channel kinetics |
| Population observations | 129 NHANES tables, codebooks, survey design, 9,254 participant IDs | Cross-sectional measured physiology and demographic variation; not temporal identification |
| Anatomical/cellular semantics | 1.76 million ASCT+B records; VCCF vessels, organ crosswalks and geometry summaries | Quantitative/anatomical identity constraints at different levels; not universal adjacency |

Exact counts, revisions, hashes, paths and source-specific limitations live in
`data/derived/collection-summary.json` and the domain reports linked from README.
The searchable evidence catalog distinguishes model candidates, model parameters,
source-reported fits, validation targets, geometry summaries and empirical priors.

## 1. Separate physical identity, anatomical support and source observation

The four IBM primitives still apply. However, the body cannot be represented by
one variable per organ with one shared topology. A single physical component may
be materialized on many vessels, epithelial cells, muscle fascicles, joint
surfaces, interstitial volumes and lymphatic segments. Each support has its own
coordinate frame, units, reference configuration, anatomical annotations and
boundary interfaces.

The source evidence already forced new distinctions in the registry: serum/plasma
sodium and potassium remain blood chemistry; they are not silently observations
of tissue-specific interstitial concentrations. Peripheral pulse counts remain
arterial pulse state, not sinoatrial pacing. Timed voided urine collections remain
collection-flow state, not an instantaneous nephron output measurement.

ASCT+B identifies structures, cell types and biomarkers. It does not directly give
channel conductance, membrane capacitance, transport coefficients or dynamic
causation. VCCF branching is an anatomical relation; venous flow direction cannot
be obtained by blindly following the same edge orientation as arteries.

## 2. Keep coherent model families and explicit specimen identity

OpenSim models carry muscle definitions, wrapping geometry, joint constraints and
parameters designed to work together. The muscle's maximum force, fiber length,
tendon slack length, pennation and path geometry form one model context; copying
only one number into an unrelated body invalidates that context.

Likewise, HumMod/BioGears parameters belong to coupled equations, feedback loops,
operating conditions, units and solver assumptions. The first integrated resting
physiology backend should preserve a coherent native model's cardiovascular,
respiratory, renal, endocrine, metabolic and fluid-exchange network. The current
IHM illustrative affine coefficients are not a substitute for that backend.

That native-backend step now executes: BioGears built locally for ARM64 and ran
its resting StandardMale patient for 60 seconds after stabilization, recording
22 physiological outputs at 50 Hz. The source equations are unchanged. See
[execution evidence](NATIVE_BACKEND.md). It remains a separate native materialization;
shared-state assimilation, conservative coupling to the 3D assets and independent
calibration are still required. Rest was configured; supine posture was not
explicitly established in this backend.

BodyParts3D, OpenSim, VMR cases and NHANES participants do not represent one person.
The default target is a population-conditioned reference body with explicit
source-specific uncertainty. A subject-specific twin additionally needs matching
subject anatomy and measurements. Sex, age, body size, health status and posture
must be carried rather than averaged away without notice.

## 3. Register 3D anatomy before coupling geometric solvers

Within one OpenSim model, local attachment frames and transforms are retained.
BodyParts3D meshes retain their native coordinates, with no guessed unit scale.
VMR preserves original meshes and fields, case IDs and unconfirmed unit details.
Semantic matches between FMA/UBERON annotations are only candidate correspondences;
they are not affine/nonrigid transforms.

A valid registration artifact must store source/target frame IDs, physical length
units, transformations, landmarks or surface constraints, residuals and uncertainty.
For mechanics it must also preserve pose and reference configuration. For flow it
must preserve lumen boundaries, cap orientation and inlet/outlet identity. Until
those exist, rendering two assets together does not establish anatomical alignment.

## 4. Couple systems through physical interface contracts

The exchange interface carries extensive fluxes (volume, mass, charge, energy)
and appropriate intensive variables (pressure, concentration, temperature,
potential), with signs, units, support mapping and time windows.

- **3D vascular to systemic circulation:** pressure/flow or impedance at named
  boundaries, integrated flux conservation, consistent wall assumptions. Regional
  CFD data alone cannot reconstruct the full vascular tree.
- **Vascular to interstitial to lymphatic:** transvascular fluid/protein exchange,
  extracellular compliance, valve-dependent lymph transport and lymphovenous
  return. Whole-body paired human lymph geometry and measured flow remain missing.
- **Muscle to skeleton:** activation, tendon force, path-dependent moment arms and
  joint constraints; source geometry and force-length parameters stay together.
- **Muscle to metabolism/perfusion:** work/heat production and substrate demand,
  then feedback through blood flow and oxygen delivery. Supine rest avoids
  locomotor mechanics but does not remove basal metabolic coupling.
- **Skin to extracellular fluid/blood:** membrane ion flux, pump stoichiometry,
  extracellular electrodiffusion, vascular substrate supply and lymphatic removal.
  BETSE provides a much stronger starting mechanism than the illustrative RC sheet,
  but its generic tissue example does not identify human integumentary parameters.
- **Endocrine/neural to organ processes:** explicitly typed receptor/rate/gain
  effects with transport delays and source-specific response curves, rather than
  arbitrary correlations treated as forces.

Native model state must not be duplicated by an independent coarse component
controlling the same conserved pool. Multirate coupling needs conservation checks
at every exchange step and checks that narrowing the explicit view preserves the
boundary response needed by the target.

## 5. Calibrate at identifiable levels, then validate externally

The evidence catalog preserves equation context, units, revisions and parameter
basis. Only three currently curated BioGears values have an explicit local chain
to source-author fit statements; many more constants may be evidence based, but
that chain has not yet been recovered. None is independently refitted merely by
copying the original number.

NHANES now supplies measured joint population-state moments across 23 physical
components. These initialize state beliefs in an exported IHM model. They do not
identify conductances, clearances, delays or feedback gains. The fit/holdout split
is by participant; positive survey/subsample weights are required. Below-detection
substitutions are retained with flags and limitations. The internal holdout is not
external validation and survey point moments are not assay-error variances.

Dynamic calibration must use informative time series or perturbations and retain
measurement operators, missingness/censoring, parameter uncertainty and
identifiability. Literature targets in BioGears can guide comparison, but reused
calibration targets cannot simultaneously count as independent validation.

## 6. Current implementation boundary

Working: raw acquisition; integrity manifests; anatomy/path/geometry parsers;
vascular field decoding; broad coefficient/validation indexing; semantic
crosswalk candidates; survey-aware population-state estimation and materialization;
BETSE execution; searchable source catalog.

Not yet working: one jointly registered 3D whole-body geometry; a connected native
whole-body physiology backend inside IHM; cross-engine conservative coupling;
patient-specific parameter calibration; independently validated human wound and
lymphatic dynamics; validated whole-body predictions.

The design now follows acquired anatomy and model equations rather than adding
more invented coefficients to the original demonstration. Existing illustrative
predictors remain executable examples, clearly separated from these data assets.
