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

All 26 tests passed. Controller stimuli are explicitly synthetic fixtures, not measured proprioception. These checks do not establish native coupled feedback, an equilibrated supported initial pose, calibrated trunk co-contraction, or whole-body stability. Native factory initialization is queued separately; no feedback advance is authorized until supported-reference evidence exists.
