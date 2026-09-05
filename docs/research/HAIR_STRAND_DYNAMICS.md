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
