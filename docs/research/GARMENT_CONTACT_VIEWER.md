# Canonical-body garment contact viewer

The “Shorts panel · elastic tissue contact” materialization plays the accepted
local mechanics record within the existing canonical body. It replaces the
surface representation of `body-bp3d-FJ3132`, `body-bp3d-FJ3133` and
`body-bp3d-FJ3134`, plus exactly 225 participating shorts triangles. These owners
are not rendered twice. Other anatomy and the remainder of both garments retain
reference coordinates during this experiment. No resting breathing trajectory is
combined with its separate mechanical clock.

The default body stays dressed. Shirt and shorts visibility remain independently
switchable; switching shorts off also switches the computed panel off. Existing
anatomy layer visibility, tissue opacity and section planes apply to the
participating surfaces. “Focus contact region” changes the camera only. Skin is
not one of the solved contact owners; hiding it for inspection is an explicit
layer choice, not physical skin deformation or tucking.

## Recorded geometry, time and precision

The viewer reads `data/derived/garment-tissue-display-v2/display.json.gz`, SHA-256
`60d0489d9353e4518a5c0be39698c63864de940878f6dd7ee3da6c69fee26634`.
There are 49 snapshots over 0.24 s, computed with a 25 microsecond native
mechanical step and exported every 5 ms. The display retains all 1,327 tissue
boundary nodes, 2,676 tissue boundary faces, 135 panel nodes and 225 panel faces.
Its float32 position quantization has a maximum recorded error of 4.77 nm.

That small export error does not establish comparable anatomical accuracy. The
underlying mechanics uses a synthesized 4 mm voxel/tetrahedral materialization
of source surfaces. Its cell diagonal is about 6.93 mm, a discretization scale
rather than a measured error bound. The API reads these values from the
hash-bound material domain manifest and exposes the occupancy construction basis.
Full tetrahedral data and the source geometry remain separately retained.

Positions are displayed at their actual recorded timestamps with no temporal
interpolation or added motion. Slower playback affects viewing time only.
Contact-force arrays are the preceding 5 ms mean, with explicitly zero force at
the initial frame. Display resultants sum those stored nodal force vectors;
maximum displacements are computed directly from stored positions relative to
their first snapshot. Other plotted metrics retain the solver's recorded values.
Unrecorded initial Jacobian/energy metrics remain missing, not invented zero.
No spectrum is supplied for this experiment.

## Owner replacement and data integrity

Before masking any garment face, the adapter validates every panel node's
retained source coordinates against the current garment constructor output,
every retained source triangle against the current garment topology, and every
local-to-source panel face mapping (allowing its documented reversed winding).
The body source IDs must exist in the canonical manifest. Invalid clocks,
nonfinite positions, inconsistent ownership and mismatched geometry fail before
scene mutation. Closing the experiment restores the original garment index
buffer and source owner visibility.

`GET /api/body/experiments/garment-contact` checks the export's manifest hash and
all source/runtime hashes, including the source skin, garment constructor,
retained generated garment, source tissue surfaces and mechanics code. The API
also binds geometric precision metadata to the exact domain manifest. Neither
the donor files nor the retained garment constructor are modified by this viewer.

Late anatomy loading immediately reapplies owner masking. Late baseline or
experiment responses cannot reset an active local clock or reopen an experiment
after the user chooses another view. Signal selectors identify and lock to the
computed contact record instead of labeling its signals as the resting baseline.

## Scientific limits

The contact solve covers an extracted panel with prescribed far-field boundary
positions and three source-shaped tissue domains. Its friction, prestrain and
constitutive values are engineering priors, not calibrated penile skin/textile
measurements. The displayed inspector exposes these values and the solver's
impulse/work receipts and unresolved edge samples. Full shorts containment,
pelvis/thigh contact, skin/fascial shells, cloth bending and continuous/self
collision remain unvalidated or absent. Small mean glans differences in a
timestep refinement do not establish convergence of every cloth node.

## Verification

Run `npm test`, `npm run build`, and
`APP_URL=http://127.0.0.1:8765 npx playwright test` from `app/` against the running
local server. The focused browser suite is
`test/browser/garment-contact.spec.js`; it uses the actual accepted artifact,
checks recorded force/time readouts, default clothes, owner-view switching and a
deliberately delayed response. Unit tests verify exact face replacement,
restoration, source mismatch rejection, missing metrics and the held export.

From the repository root, run
`OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_garment_contact_viewer.py`.
It validates the held artifact and tests changed-source, changed-export and
escaping-path rejection using isolated temporary files. Its receipt is retained
under `artifacts/verification/garment-contact-viewer/reader.json`.

The same change updates the client scenario validator to admit every currently
registered native implementation through `whole_body_integrity_energy`, while
continuing to reject unknown implementations. A regression exposed that the old
three-entry list prevented the newer native default from reaching transport.

The completed verification run passed 29 unit tests, the production build, four
reader checks, and 28 browser tests. The existing opt-in browser case that
launches a separate native hemorrhage scenario was skipped. The full browser
run included the actual garment artifact, native apnea playback, regional
studies, native scenario transport, source anatomy and responsive panel controls.
Logs and inspected garment screenshots are retained under
`artifacts/verification/garment-contact-viewer/`.
