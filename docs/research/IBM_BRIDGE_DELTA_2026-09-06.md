# IBM bridge source delta and repeatable import path

Read-only audit of the developing IBM-1 checkout on 2026-09-06. No donor source,
training process, checkpoint or dataset was changed or loaded. Source ASTs and
package-file hashes were compared; one isolated inhibition-function evaluation
used NumPy. This is compatibility evidence, not whole-brain validation.

## Exact identities

| Source | Identity |
| --- | --- |
| Current clean IBM-1 HEAD | `398375cc993ac17907d4ddd9f6d53eed2fab31de` |
| Current `ibm/` package hash, same file-map algorithm as IHM vendor | `f969da521ef87ece047da189de17118538a650ebab65605a3254a6be38dffa60` |
| Existing IHM package pin | `0bcb65606a350d586c1c5fc1634913aee7f90c7de4201eef57770c6af715346f` |
| Existing donor commit receipt | `63a5668f6e456b8ca5df5848e5d839cf77357bc7` plus recorded dirty source files |

There are 29 changed and six new package files, with no deletions relative to the
pinned artifact. The old receipt includes uncommitted runtime/materializer edits;
its commit alone cannot reproduce it. The per-file hashes are the actual identity.

IHM currently has two separately pinned paths:

- `IBMBackend` verifies the complete `ibm-backend/source/ibm` package and enforces
  `PINNED_PACKAGE_SHA256`. `SnapshotLoader` executes captured verified source bytes,
  bypasses bytecode and refuses a different IBM identity in an already-imported
  process. The tactile adapter uses these actual source operators.
- `BodyBrain._load_rate_law` executes three exact function ASTs from the independent
  `brain-sources/ibm-neural.py` receipt. Its regional graph and physiological,
  sensory and motor conversions remain IHM engineering priors. Updating only the
  complete package pin does not update this regional ODE law.

## Compatibility changes that matter

**Regional neural law:** `_sigmoid` and `wilson_cowan_excitatory` have identical
ASTs. `shunting_inhibition_rate(x, theta)` retains its signature but changes its
normalization: the older law divides by an inhibition-shortened `tau_eff`, whereas
the current law divides by `tau_membrane_s * g_leak`. A bounded source-function
probe at V=-60 mV, inhibitory conductance=0.2, Erev=-70 mV, tau=0.015 s and
leak=1 gives -160 versus -133.33333333333334 mV/s. This is a real behavioral change,
not an import-only update. Re-run physiological perturbation, step-size, checkpoint
and skin-to-motor acceptance before changing the runtime pin. No trained posterior
has been imported or inferred from this correction.

**Sensorimotor:** `ibm/processes/transduction.py` is byte-identical, including the
rapid/slow and thermal static/dynamic transfers. This supports source compatibility
for the existing bounded receptor channels, but not their surrounding registry,
materializer or runtime after upgrade. Current IBM substrate documentation still
explicitly lacks peripheral receptor/motor-unit geometry; no full afferent/efferent
wiring should be inferred from the new brain source.

**Vascular:** current `vascular_tree_adjacency` accepts a diameter-dependent callable
`viscosity_pa_s` as well as a scalar. New `vascular_prior` supplies tiered measured
statistics, synthesis, support and theta-prior APIs. Vascular processes now read
`structural.lumen_radius`, `segment_length`, `branch_order` and capillary density;
these are structural inputs, not blood state to be overwritten each exchange.
An IHM import must preserve vascular-tree parent/radius units and geometry provenance.
Its native regional fluid owner remains authoritative; importing an IBM vascular
prior is not authorization to create a second blood inventory.

**Microcircuit:** the existing within-site builder signature is unchanged. New
`microcircuit_prior` adds `parameter_priors`, `laminar_gain_priors`, `theta_prior`
and `microcircuit_generative`, with tiered MICrONS/H01/Allen-derived summaries.
Mouse circuit connectivity, limited human census and population kinetics are
separate evidence scopes. These APIs do not make IHM's current 80-population
regional graph an identified microcircuit.

**Materialization/runtime:** `build(..., substrate=None)` is new;
`MaterializationRequest` adds `max_tier` and `tier_ceilings`. Use explicit ceilings
and retain emitted tier records. `advance` adds a `propagate` argument, and internal
solver handling changes despite largely stable public signatures. Whole graph
runtime compatibility and uncertainty behavior have not been executed in this audit.
Do not call `build_substrate` merely to test imports: it resolves template/atlas
assets and can materialize substantial geometry.

The new prior modules load small evidence JSON when available and otherwise use
compiled-in statistics. Their documented fallback is intentional and labeled by
`from_evidence`; source-package hashing alone does not describe which evidence
path ran. For reproducible imports, explicitly choose compiled-in statistics or
archive and hash the small measurement summaries, then record the choice. Relevant
paths are `microvasculature@1/measured.json`, `microcircuit_statistics@1/measured.json`,
`human_microcircuit_census@1/measured.json` and `population_kinetics@1/measured.json`
under the donor source-card evidence directories. No raw microscopy or learned
checkpoint is needed for this step.

## Safe repeated import, without donor changes

1. Resolve a donor commit and export its committed `ibm/` bytes into a **new**
   content-addressed IHM candidate directory, preferably using `git archive` for a
   single immutable revision. Record commit, file map and aggregate hash. If using
   `scripts/vendor_ibm_backend.py` on the working tree, require a fresh output and
   verify HEAD/status plus all source hashes before and after copying. Its current
   default overwrites the active artifact; do not use that default for updates.
2. In a fresh process, verify the complete candidate and its explicitly selected
   evidence summaries. Keep old sessions on their original source identity;
   `SnapshotLoader` intentionally disallows hot-swapping IBM in one interpreter.
3. Select the candidate with an explicit version/hash policy in the IHM adapter.
   The current hard pin correctly rejects this new hash; copying files alone is
   not a completed upgrade. Do not loosen verification to accept arbitrary source.
4. Bind the regional `BodyBrain` law to the **same candidate neural source** and
   record its hash while retaining existing canonical geometry. Do not rebuild
   anatomical assets or hand-copy the corrected inhibition formula. A versioned
   regional source manifest is sufficient for this source-law update.
5. Run bounded existing `verify_ibm_body_parity.py`, `verify_ibm_causal.py`,
   `verify_body_brain.py`, sensorimotor/cutaneous and runtime timing fixtures, plus
   the old/new inhibition regression. Only then start new body sessions on the
   candidate. Existing sessions and their receipts remain reproducible.

The next useful increment is a versioned candidate importer and a shared source
identity for package-backed receptors and regional neural laws. Vascular and
microcircuit materializations should follow as separate bounded integrations with
explicit supports, units, tier provenance and exclusive native storage ownership.
