# Native receptor endpoints into the shared IBM brain

This opt-in bridge reads the actual native Nervous chemoreceptor, carotid/aortic baroreceptor and pulmonary-stretch rates in Hz. It does not implement receptor equations, advance a second brain, control breathing or claim measured axonal recruitment.

`native_nervous_afferents.h` obtains standard C++ base-class member pointers through protected access, then reads the actual engine-owned `Nervous` object checked by RTTI. It constructs no observer object and performs no fictitious derived cast, offset arithmetic, access-control macro, source mutation or state write. This is explicitly a source-internal observer; the isolated build records the exact Nervous header/source, observer header and native library identities. Separate runner wrappers include the existing signed/regional adapters, replacing only the observation call. Existing runners and the composition source freeze remain untouched.

`EmbodiedRuntime.from_workspace(..., native_afferent_allocation=...)` selects the corresponding observation adapter. Allocation must explicitly provide all four channels, with fractions summing to one across existing medulla IDs. The current IDs are `brain-bp3d-FJ1769` and `brain-bp3d-FJ1831`. No default allocation is supplied. Native absolute firing is allocated as additive IBM sensory Hz—an inferred regional recruitment prior, not an identified excess-rate conversion or measured baseline. Existing single-source IBM candidate selection is supported through the shared brain's immutable source identity.

Every endpoint receipt binds native library/executable/initial-state/manifest SHA256 identities, observer header SHA256, native command sequence, elapsed 20 ms tick, absolute time/origin and all four actual rates. The receipt has its own canonical SHA256. Missing, nonfinite, out-of-range, stale, repeated, altered and wrong-owner endpoints reject. A receipt at time t feeds the sole shared-brain advance over [t,t+20 ms]; the newly accepted native endpoint is reserved for the next interval. No baseline is fabricated. Full native rates and each allocation contribution remain in the frame audit.

The optional `native_sensory_blocks` command names native channels. A block removes their input immediately; release uses the next actual receipt without creating a second adaptation state. Pre-native failures restore bridge and brain state. Missing/invalid observation after native mutation terminates the runtime. Bridge checkpoints bind native/source/brain/allocation identity and replay exactly; they do not make native serializer replay exact. No efferent respiratory or autonomic writer was added.

Light verification:

```
PYTHONPATH=.:scripts OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m unittest scripts.verify_native_afferents scripts.verify_embodied_runtime scripts.verify_regional_embodied_factory scripts.verify_candidate_embodied_factory scripts.verify_body_input_capabilities
```

All 25 tests pass, including real shared IBM causal response versus blocked control, exactly one brain advance, clock/source/tamper guards, replay, prior-endpoint runtime consumption, rollback and uncertain-native abort. Native receipt fixtures in this suite are synthetic and labeled accordingly.

`build_native_afferent_adapter.py` prepares separate signed or `--regional` runners from an accepted parent adapter and its frozen inputs; `--run` is heavy-queue gated. Current regional preparation correctly refuses the installed parent because `scripts/native_tissue_ports.h` differs from its receipt. A refreshed accepted parent adapter must precede observer compilation and native observer/control parity acceptance. New composition DSOs alone do not satisfy that requirement. No observer native run, default promotion or coupled respiratory-control acceptance is claimed.
