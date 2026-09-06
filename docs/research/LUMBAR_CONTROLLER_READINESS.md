# Opt-in lumbar controller and factory readiness

The retained 98-muscle registration at `data/derived/lumbar-muscle-native-lb45uirs/variant/registration.json` appends six Gait2392 donor muscles to the unchanged 92-muscle catalog: bilateral `gait2392_ercspn`, `gait2392_intobl`, and `gait2392_extobl`. Each retains source force and fiber-length normalization and an explicit contralateral cortical engineering prior. These labels do not establish measured trunk recruitment, spinal reflex connections, or a standing policy.

`EmbodiedRuntime.from_workspace(..., augmented_registration=<workspace-relative registration>)` now permits this explicit override. Omitting the keyword preserves the original 92-muscle default. The override validates model/catalog/source hashes, exact muscle IDs and normalization before starting either native owner. Symlink redirects and outside-workspace paths are rejected; files are bounded to 32 MiB each and 64 MiB total. The runtime retains all verified bytes under `inputs`, records registration/model/catalog hashes, and requires the initialized plant catalog to equal the checked controller catalog. Existing native library/session identity, mechanical checkpoint ownership and source stability checks remain active. No API selector or automatic activation/default promotion was added.

The bounded verifier `scripts/verify_lumbar_controller_readiness.py` checks:

- Exact six donor IDs, preserved 92 rows, and source/model/catalog normalization.
- Each added muscle's synthetic CE observation becomes a normalized delayed afferent in the real IBM regional brain. Explicit descending commands drive only the named added effector; sensory/motor blocks remove arrived state and queued commands immediately, with delayed release.
- Original 92 motor outputs and brain state remain exactly equal when the six added sensory and descending ports are silent. Added nonzero afferents intentionally influence their shared cortical region.
- Reversed catalog ordering preserves outputs within `1e-12`; checkpoints nevertheless retain exact catalog/order/source identity and reject a different catalog or prior. Same-controller replay is exact.
- Mock factory forwarding, complete frozen source receipts, and rejection of altered/symlinked model files before native startup.

Validation command:

```
PYTHONPATH=.:scripts OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m unittest scripts.verify_lumbar_controller_readiness scripts.verify_regional_embodied_factory scripts.verify_candidate_embodied_factory scripts.verify_embodied_intake_mass scripts.verify_sensorimotor
```

All 26 tests passed. Controller stimuli are explicitly synthetic fixtures, not measured proprioception. These checks do not establish native coupled feedback, an equilibrated supported initial pose, calibrated trunk co-contraction, or whole-body stability. Native initialization is validated below; no feedback advance is authorized until supported-reference evidence exists.

## Actual zero-step native acceptance

`data/derived/audits/lumbar-controller-factory-miey4gjq/verification.json` records a passing 98-muscle factory initialization with `regional_skin=True`, actual native CE observation validation for the six added muscles, exact controller checkpoint restore, and cleanup. Wall time was 1.777 s, parent peak RSS 771,544 KiB and maximum child peak RSS 196,028 KiB. Physiology exited by controlled SIGTERM (-15), mechanics exited 0, and both were reaped. All clocks remained zero; no feedback advance or supported equilibrium was tested.

The run retained these SHA256 identities:

- Registration: `6d94486340828c9a06147e119d7fd8fdfa2393fd8d68fa4cb2f5e030e4a83492`
- Model: `7f6c40224fbc43359d6b2997efeb499c5ffe4d19c786d37b849312efd20d2114`
- Catalog: `106e4f0d33f8f6624f1ee9eb0660ff988a161c228c0e6a524a662417e489d1ee`
- Physiology library: `bc91cbab829c04bfa2df7490433af5bee752e7514332ac715e9f0bd11a17650a`
- Initial physiological state: `cba7ffb523728d86e8f07b2a30b215b040676c30580e8accda722d175bd687b3`
- Regional neural law: `f8c7d3172af71b399d9e6301f068d8fde8adf3363941ef0d3d3ec50ef62f663b` (legacy preserved law; no candidate selected in this run).
- Factory manifest: `b0f6626b7fecadf4f6ed029027b8c1347ad15856e810362b848977fc44f6ce26`

The verifier's `--native` mode is separately queue-gated, with a 30-second alarm and 2 GiB soft/4 GiB hard address-space limits. It must run with native threads limited to one and nice level 10. Default invocation remains a lightweight suite.
