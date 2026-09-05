# Systemic experiment viewer

The viewer keeps one canonical body and selects a recorded whole-body protocol.
Systemic records have their own active trajectory, native variable metadata and
finite-window spectra; they do not replace the saved resting body trajectory.
Pending loads have no playable stale trajectory. Selecting a regional study,
source inspection, or another protocol invalidates the previous request epoch.

Respiratory records can supply an explicitly computed body projection. Sparse
six-hour digestion records retain reference geometry and expose their actual
sample interval. No oscillation is synthesized between samples. Garments remain
enabled independently and follow only a supplied computed thoracic field.

Implementation plan: first verify the signal/spectral adapter against distinct
frequency axes, physical units, null values and rejected uncomputed motion;
then connect independent experiment selection, native mechanism readouts,
recorded time playback and PSD/finite Laplace plots. Verify request races and
mode isolation in the browser, followed by real held respiratory artifacts.

The readouts describe native model execution. Mechanism descriptions and
limitations travel with each source-pinned record. Spectra remain finite-window
descriptors, not inferred physiological poles or proof of causal coupling.

## Implemented behavior

The left panel selects a completed native protocol, while the right panel shows
selected native values and the executing mechanism descriptions supplied by the
record. Field labels retain source IDs internally and expose owner/support
metadata on hover. The plot follows the selected snapshot, including explicit
unavailable values. PSD and finite Laplace magnitudes use the record's distinct
frequency axes and physical units. Signal-source controls identify and remain
locked to the active experiment. Playback speed is explicit; no additional
physiological samples or body motion are interpolated.

The viewer retains source-fitted default clothing. A computed respiratory field
deforms the shirt kinematically. Sparse records with `has_body_projection:false`
must have empty entity transforms and no respiratory displacement field; the
adapter rejects records that violate that boundary. A projection flag must be a
boolean, and sample times must match the declared cadence.

Completing an anatomy mesh request applies its geometry at the selected body
time. It no longer resets the timeline. Similarly, a late resting-body record
updates its saved data without resetting an active regional or systemic clock.
These two races were reproduced with delayed-response browser regressions.

## Verification, 2026-09-05

The adapter/unit suite passes 25 tests. The final systemic and regional browser
suites pass six checks, including the real 1,801-frame native apnea record,
null-value handling, distinct finite-transform units, no synthesized sparse
motion, canceled requests, late anatomy geometry and late resting trajectories.
The body suite also passed four checks after the clock fix. The production app
build passes. Screenshot and verification receipt are retained under
`artifacts/verification/systemic-viewer/`.

A broader browser run initially passed 24 checks, skipped the opt-in native
scenario test, and exposed the regional clock reset. The deterministic regression
and final targeted run confirm its correction; that earlier failed run is not
reported as a fully passing suite.

The real browser acceptance currently covers short respiratory experiments.
Sparse long-duration behavior is checked with an explicit fixture while the
matched six-hour native computations finish. The index refreshes periodically
and can discover those completed records without changing source families.
