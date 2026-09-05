# Retained native environment inputs

`ihm/assembly/native_environment_evidence.py` retains the input bytes that the earlier source snapshot omitted: startup-identified linked libraries, the native executable, its ELF interpreter, and the explicit native runtime resource trees. It copies regular files into a closed, content-addressed archive inside each experiment. Original absolute dependency paths remain provenance metadata; they are never returned as read dependencies or used by archive resolution.

## API and integration point

```python
from ihm.assembly.native_environment_evidence import (
    freeze_native_environment, resolve_native_environment,
)

with NativeSession(config, native_directory) as native:
    # NativeSession has received ready. Capture before save/meal/step actions.
    immutable_inputs = freeze_native_environment(root, native_directory)
    # Continue the original shared native process.

# Later publication checks the archive, not mutable external libraries.
immutable_inputs = resolve_native_environment(root, native_directory)
```

Both functions return workspace-relative paths mapped to SHA-256 digests. The mapping includes the original native manifest, the archive manifest, and every archived blob. Publication should merge it into source dependencies for **each** compared experiment, including controls. No API call edits `systemic.json`, `source_receipts.json`, the native manifest, or original binaries/resources.

The archive is `native/environment-inputs/manifest.json` plus `blobs/<sha256>`. It preserves 38 startup-linked dependency files in the current native build, including system libraries and `libbiogears_cdm`, and the executable. PT_INTERP is parsed directly from the retained ELF bytes without executing a loader or shell. The legacy native manifest omitted the direct loader line in `ldd`; the interpreter therefore has an explicitly current-capture identity, not a claimed startup hash.

The resource allowlist is `patients`, `substances`, `environments`, `nutrition`, `config`, `ecg`, `xsd`, `UCEDefs.conf`, and `BioGearsConfiguration.xml`. Capture follows these explicit resource roots and records logical paths, empty directories, modes, original identities and content hashes. It does not walk the experiment directory, state outputs, logs, or simulation trajectories. A second inventory verifies content and tree membership after copying. Per-file checks compare source identity and bytes before/after the copy. Copies use independent storage; resolution rejects symlinks, hardlinks and paths outside the archive.

## Capture timing is part of the evidence

Capturing after native `ready` is still **after initialization**. Existing runs have no startup resource inventory against which to compare those bytes. The receipt always records `resources_verified_at_process_start: false`; it distinguishes `post_start_before_first_action` from `post_start_after_commands` using the command-receipt prefix at capture time. The prefix length/hash and command count are metadata; dynamic journals themselves are not archived or bound as immutable inputs.

For future runs, calling immediately after readiness prevents earlier user actions from preceding capture, but does not retroactively prove what initialization read. For a strict startup resource guarantee, a future launcher must inventory and hold resources before spawning, then validate the same versions after readiness. This module does not claim that stronger guarantee.

Linked dependency copies must match the exact hashes in the native startup manifest. If the current executable or main library has changed, the module can recover the matching version from that run's already-validated adjacent `frozen-sources.json` archive. The source of copied bytes is recorded separately from the original startup path. It never substitutes a newer binary, searches unrelated directories, or relaxes workspace containment for historical copies.

## Failure and historical retention

Capture stages into a new `.environment-capture-*` directory. Successful capture atomically renames it to `environment-inputs`. Failure leaves its copied bytes and `failure.json` in place. A subsequent attempt uses a new staging directory. An existing successful archive is validated without repinning it, even when original resources, libraries, or code have evolved.

The first real respiratory capture encountered a newer live executable and correctly failed. That failure remains retained. The next capture recovered the old executable from its prior source snapshot and completed without changing the original run. Resolution independently checks every startup dependency identity, executable identity, archive closure, and content digest. External original paths are not opened during resolution.

## Verified captures

Nine existing runs were captured: `respiratory_v3/{rest,apnea}`, `six_hour_v3/{hydration,meal}`, and `exertion_v2/{rest,hydration,meal,exercise,meal_exercise}` beneath `data/derived/systemic/`. Each has 214 verified paths: two manifests and 212 unique file blobs, including 172 runtime resource files. Retained blob bytes total approximately 84.7 MB per run. All nine correctly record capture after commands, including active long-running experiments. These are post-start environment observations, not newly rerun or retroactively revalidated physiology.

Run behavioral verification and optionally capture a real native directory:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_environment_evidence.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_environment_evidence.py --capture data/derived/systemic/respiratory_v3/rest/native
```

The verifier uses real temporary files and an ELF executable. It checks detached survival after original library/resource mutation and executable deletion; rejection of corrupted blobs, external symlinks, hardlinks, and omitted startup dependencies; exclusion of state outputs; empty-directory retention; explicit capture timing; historical executable recovery; and preservation of failed capture attempts. API-missing, missing-dependency, hardlink, and historical-executable recovery tests were observed failing before implementation of the corresponding behavior.

An archive is evidence preservation, not a complete portable runtime or a bitwise replay guarantee. The legacy manifest does not identify arbitrary later `dlopen` modules, all process environment values, kernel/CPU behavior or native serializer omissions. Resource timing remains qualified above. Original libraries retain their licenses; local copying does not grant redistribution rights.
