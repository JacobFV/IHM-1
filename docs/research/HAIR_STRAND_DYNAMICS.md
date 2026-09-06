# Rooted elastic hair, retained source populations, and display resolution

The current materialization is `data/derived/hair/elastic_v3`. It retains 92,934 scalp roots and the existing 201,528 body roots. The viewer displays 2,048 scalp fibers and 512 body fibers at their physical radius, driven by 64 and 32 independent elastic guide strands respectively. This is 2.204% and 0.254% of the retained populations. Render density is an explicit level of detail; it does not thicken fibers or claim full-population dynamics. The previous rigid `body-detail-hair` display is replaced, while its source geometry and root sample remain intact. Earlier v1/v2 artifacts are retained and not registered.

## Acquired human measurements

Primary PDFs and text are retained under `data/raw/hair`, with URL, size, SHA256 and acquisition timestamp receipts. The tracked source card `data/sources/human-hair-mechanics.json` also embeds the receipts. `data/measurements/hair/strand_reference.json` supplies the builder's numerical values:

- [Müllner et al., 2020](https://www.mdpi.com/1420-3049/25/9/2143): tensile modulus 7.11 ± 0.44 GPa (Table 2) and density 1,312 ± 43 kg/m³ (section 2.5). Commercial human hair; tensile test uses a 30-fiber bundle, 10 mm gauge, 0.01 mm/s crosshead. Density uses ten hairs from two individuals. These are reference measurements, not a canonical subject posterior.
- [Ezawa et al., 2019](https://www.mdpi.com/2079-9284/6/2/24): 5.7 GPa bending modulus is the representative example in Figure 2c, **not a cohort mean**. The study includes 156 Japanese females aged 10–70, with 30 fibers per bending measurement at 20°C and 65% RH. The model uses the paper's `B = E I` relationship. Circular section is an explicit simplification of real, often elliptical hair.
- [Alsharif and AlGhamdi, 2022](https://www.dovepress.com/article/download/80369): 120 Arabian adults, 60 men and 60 women, aged 18–60. Table 2 frontal/vertex/occipital densities are 143.9/147.1/153.6 per cm² and diameters 83.5/87/90.7 µm. Reported regional SDs are retained but not sampled as fiber-level variation. This cohort is not matched to the canonical atlas subject.

The scalp haircut is an explicit 30 mm, initially straight scenario. Body lengths 0.5–2 mm and radii 8–16 µm remain the existing vellus morphology prior. Scalp mechanical properties are transferred to body hair without empirical validation. Humidity, chemical treatment, and temperature dependence are not modeled. The acquired 2016 structure/mechanics paper supports these limitations; its parameters are not used. The retained XPBD paper was investigated, but the initial prototype failed the fixed-root energy test and was replaced.

## Mechanics and attachment

`app/src/hair_dynamics.js` implements untwisted, small-deflection discrete beams in SI units. For segment spacing h and circular area A and second moment I:

- A = πr²; I = πr⁴/4; lumped mass = ρAh, half mass at endpoints.
- Axial energy = EA/(2h) times squared longitudinal displacement difference.
- Bending energy = EI/(2h³) times squared transverse second difference.
- Backward Euler solves `(M/dt² + K) u = M/dt² u_predicted` by direct Cholesky factorization of each small strand system. Default internal step is at most 1/240 s.

Root and first segment are clamped. Subsequent particles respond to gravity and supplied forces in newtons. Root motion can do mechanical work; the energy decay test therefore holds roots fixed. The linear model reports maximum transverse slope and flags values above 0.3. This threshold is a declared small-deflection operating bound, not a failure calibration. There is no hair/body or hair/hair collision, torsion, aerodynamic drag, growth, damage, or two-way reaction force into skin.

Roots are actual source triangle barycentric coordinates, bound to the skin geometry SHA256. The scalp coverage is an atlas-frame prior: `abs(x)<0.105` and `(y>0.79 or (y>0.72 and z<0.015))`, all three corners included, outward face-normal/radial dot product above 0.45. The mask is not a measured scalp segmentation. Body roots reuse the retained source IDs, triangle indices and barycentrics; source entity/hash and reconstructed positions are checked before materialization.

Displayed fibers retain their own material roots and radii. Three nearest guides in the same retained anatomical region provide inverse-distance weighted displacement at normalized strand arc length, scaled by render/guide length. They are not independently solved beams; material-frame differences and sparse guide resolution limit fidelity. The controller reports both guide and rendered transverse slopes. Tube geometry has four sides per node and an existing buffer is reused. No per-frame Three.js TubeGeometry creation occurs.

## Viewer and public API

```js
const controller = createHairController(geometry);
const result = controller.update({
  frame, referenceCentroids, simulationTime, recordKey,
  active: true, visible: true, gravity_m_s2: [0, -9.81, 0],
  external_forces_n: null, // optional SI vector per guide particle
}, mesh.geometry.attributes.position.array);
if (result.diagnostics.updated) mesh.geometry.attributes.position.needsUpdate = true;
```

Output already includes respiratory skin deformation and the skin entity body transform. Set the object matrix to identity and bypass the old rigid attachment transform. Default update limit is **10 Hz**, with a hard configuration ceiling of 20 Hz. Time comes from the caller's mechanical simulation clock; there is no wall-time catchup. Hidden/inactive calls do not advance mechanics. Resume, rewind, record change, or a gap greater than 0.25 s resets to the current roots and declares the reason. Thus coarse 5 s archived snapshots do not masquerade as continuously integrated hair motion. Same-clock calls return cached geometry. No work is required for a stationary viewer with no new frames.

`GET /api/body/experiments/hair-strands` validates source and geometry hashes. `ImplicitHuman.materialize('hair-strands')` additionally checks its frozen manifest/module assets. Run `scripts/register_hair_strands.py` after rebuilding the app manifest to restore the two verified geometry routes. The registration stores population, guide, display density and evidence metadata in each structure. It does not change native physiology or canonical root populations.

## Verification and limits

Run with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10`:

- `node --test app/test/hair_dynamics.test.js`: ten tests cover mass/radius, rooted gravity bending, moving clamp, deterministic checkpoint, invalid input isolation, non-increasing fixed-root total energy, timestep refinement, force contrast, rate/pause/reset, and render transfer.
- `node scripts/verify_hair_strands.mjs`: verifies actual asset SHA, guide/render barycentric roots, finite geometry, slope and strain, actual gravity displacement, cached updates and pause. At 0.1 s the v3 scalp maximum guide coordinate displacement is 0.344 mm, maximum transverse slope 0.01732; body displacement is 30.36 nm, slope 0.0000350. Physical short stiff body hair should not be exaggerated into visible waving.
- The 30 mm fixture tip at 0.1 s differs by 10.53 µm between 1/240 and 1/480 s internal steps, below the declared 15 µm numerical tolerance. This is numerical evidence, not comparison against filmed human hair.

Actual asset verification and small CPU timing receipts live beside the materialization. Timings were collected under competing workload and are not a real-time guarantee. No browser or GPU acceptance was run by this worker, and no actual strand collision is claimed. Full density rendering, nonlinear curls, contacts and subject-specific hair calibration remain separate work.

The viewer integration keeps strand dynamics **opt-in and off by default** under the current competing IBM workload; static physical strands remain visible. Final optimized warm scalp timings were 37.6, 313.3 and 14.0 ms, with body timings 0.66, 42.6 and 11.2 ms. The large variation prevents a sustained 10 Hz claim. No more performance loops were run. Off-thread execution would reduce UI blocking but would not eliminate the underlying CPU cost.

`scripts/verify_hair_registration.py` passes four lightweight acceptance branches: valid public materialization, changed geometry rejection, changed module rejection, and frozen public module rejection. `scripts/acquire_hair_sources.py` rechecks all five retained PDF byte hashes and can reacquire exact source bytes if absent; publisher byte changes fail visibly.

## Bounded worker execution (2026-09-05)

The opt-in viewer now runs the same controller in one shared dedicated module
worker (`app/src/hair-worker.js`). A single solve runs at once; each hair region
holds at most one newest pending request. Only skin inputs and the strand/root/
guide data cross the worker boundary. Tube coordinates return by transferable
buffer. The controller keeps its own cache; transferred coordinates are a copy.
Disabled or hidden regions cannot apply stale worker results. Worker failures
are displayed and do not trigger a synchronous solver fallback.

Coalescing queued frames does not advance an invented animation clock. Submitted
timestamps drive the existing fixed-substep solver with endpoint root
interpolation. Intermediate input samples are discarded under backpressure;
this is an input-sampling assumption, not preservation of every forcing history.
Existing gaps over 0.25 s, resume, rewind, and record changes still reset. The
readout marks reset intervals **NOT physically integrated**, and reports the
accumulated forward `unsimulated_total_s`. The model computes no follicle/body
reaction impulses and has no body, garment, or self-contact. Adding coupled
impulse exchange would require a different queue protocol that conserves it.

Small source-only checks:

```
node --max-old-space-size=256 --test --test-concurrency=1 app/test/hair-worker.test.js app/test/hair-view.test.js app/test/hair_dynamics.test.js app/test/hair_dynamics_scratch.test.js
node --max-old-space-size=256 scripts/benchmark_hair_worker.mjs
```

One 128-guide, seven-node, 40-update fixture run used 273.33 ms for synchronous
main-thread solving versus 0.476 ms for main-thread worker dispatch; worker round
trips took 257.00 ms, excluding startup. The main event loop executed 201 timer
ticks during worker execution. Final tube coordinates matched exactly, with no
unintegrated intervals. Process RSS was 121.54 MiB. These are noisy fixture
measurements, not a full-viewer responsiveness guarantee. Normals, bounds, and
GPU uploads remain on the main thread; the off-thread solver still costs CPU.

## Opt-in authoritative mechanical hair prototype

`ihm/assembly/hair_dynamics.py` owns a separate SI beam state for explicit
mechanical experiments. It ports the straight uniform small-deflection beam
model, including physical guide masses, axial/bending stiffness, moving
two-node clamps, implicit substeps, and checkpoint/restore. It does **not**
silently replace or synchronize the JavaScript view solver.

`ihm/assembly/hair_feedback.py` connects that state to the same registered
native body load/checkpoint interface used by `GarmentFeedback`. It computes
free-node contact against `MovingSurfaceContact`, distributes equal/opposite
surface impulses through retained vertex owners, and includes both force and
moment from the moving follicle clamp. Native/body and hair intervals are
iterated from common checkpoints; nonconvergence rolls both back. Acceptance
also requires root and first-segment clamp synchronization within 1 nm, even
when the requested whole-surface convergence tolerance is looser. No root snap
or unreported state reset repairs a mismatch.

The receipt separates physical impact/friction dissipation, external gravity
work, hair energy change, numerical beam energy defect, native interface work,
and interface quadrature discrepancy. It reports linear/angular impulse
residuals and small-deflection validity. Numerical losses are not relabeled as
physical heat. The finite-body fixture includes translational mass and
rotational inertia; it exercises moving contact, force/torque feedback,
momentum accounting, energy accounting, exact replay, and rollback. Run only
the small fixture with:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_hair_feedback.py
```

The bounded check completed in 0.22 s with 77,944 KiB peak RSS. Its one-guide,
four-node moving-body interval produced two contacts; linear/angular momentum
and the coupled energy receipt agree within 1e-18 in SI units. Impact loss was
6.20e-15 J, numerical beam energy defect −4.33e-14 J, and interface work
quadrature discrepancy 4.90e-21 J. These are synthetic numerical checks, not
anatomical calibration or a production throughput claim.

This is a bounded physical prototype, **not live anatomical or garment hair
coupling**. Contact is discrete nearest face-interior centerline-node contact
with zero radial clearance, a 4 mm search band, and no shaft/edge/CCD/self
certificate. Existing Coulomb coefficients are explicit caller inputs. The
surface patch is limited to 65,536 vertices/131,072 triangles; candidate
queries retain the existing 100,000-pair batch cap. Hair is limited to 512
guides, 3–32 nodes per guide, and 240 substeps per parent interval. These caps
reject oversized requests; they do not reduce the source strand population.

Required before anatomical live enablement:

- Join attachment `sample_ids` to population `ids` and retained `face_index`;
  current view attachments contain triangle coordinates and source SHA but
  omit source face/node indices. Verify source bytes/hash and tangent rules;
  do not substitute nearest-face registration. The low-level constructor
  verifies barycentric root reconstruction and fixed owner indices, but does
  not itself verify anatomical source provenance.
- Provide source-bound material vertex owners and a physically justified
  inertial policy. Native canonical hair proxies already contribute mass;
  added physical guide masses must be reconciled explicitly. Do not subtract
  an entire proxy or multiply guide masses by rendered-population counts.
- Attach the Python owner explicitly to an authoritative native session, and
  bridge its JSON-ready state/receipt to the viewer. `frame()` currently marks
  `viewer_synchronized: false`. The JavaScript worker has no reaction sink and
  must remain separate until this ownership transfer exists.
- Follicle faces spanning different rigid owners are rejected because their
  clamp orientation derivative is not implemented. Whole-body source surfaces,
  live native execution, browser display, and deformable garment contact have
  not been validated by this fixture.
