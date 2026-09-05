# Retained native environment inputs

`ihm/assembly/native_environment_evidence.py` retains the input bytes that the earlier source snapshot omitted: startup-identified linked libraries, the native executable, its ELF interpreter, and the explicit native runtime resource trees. It copies regular files into a closed, content-addressed archive inside each experiment. Original absolute dependency paths remain provenance metadata; they are never returned as read dependencies or used by archive resolution.

## API and integration point

```python
from ihm.assembly.native_environment_evidence import (
    freeze_native_environment, resolve_native_environment,
    materialize_native_resources, resolve_native_execution_inputs,
)

# Inside the process owner, after writing manifest.json and BEFORE Popen:
immutable_inputs = freeze_native_environment(root, native_directory, before_start=True)
resource_tree = materialize_native_resources(root, native_directory)
# The owner must select resource_tree's roots for startup before spawning.

# Legacy after-readiness capture remains available, with honest timing:
with NativeSession(config, native_directory) as native:
    # NativeSession has received ready. Capture before save/meal/step actions.
    immutable_inputs = freeze_native_environment(root, native_directory)
    # Continue the original shared native process.

# After selection, the owner writes execution-inputs.json before Popen.
# Later publication validates selected detached inputs for pre_start archives.
immutable_inputs = resolve_native_environment(root, native_directory)
```

The freeze/resolve functions return workspace-relative paths mapped to SHA-256 digests. Capture returns the original native manifest, the archive manifest, and every archived blob. Public resolution additionally requires the execution selection for `pre_start` captures and binds its receipt, patient/state input, detached materialization receipt and every selected resource file. Publication should merge the public resolver mapping into source dependencies for **each** compared experiment, including controls. Capture and materialization deliberately use an archive-only internal check: the execution selection does not exist at those earlier lifecycle stages. `materialize_native_resources` returns an absolute `Path` to `native/runtime-resources`. No API call edits `systemic.json`, `source_receipts.json`, the native manifest, or original binaries/resources.

The archive is `native/environment-inputs/manifest.json` plus `blobs/<sha256>`. It preserves 38 startup-linked dependency files in the current native build, including system libraries and `libbiogears_cdm`, and the executable. PT_INTERP is parsed directly from the retained ELF bytes without executing a loader or shell. The legacy native manifest omitted the direct loader line in `ldd`; the interpreter therefore has an explicitly current-capture identity, not a claimed startup hash.

The resource allowlist is `patients`, `substances`, `environments`, `nutrition`, `config`, `ecg`, `xsd`, `UCEDefs.conf`, and `BioGearsConfiguration.xml`. Capture follows these explicit resource roots and records logical paths, empty directories, modes, original identities and content hashes. It does not walk the experiment directory, state outputs, logs, or simulation trajectories. A second inventory verifies content and tree membership after copying. Per-file checks compare source identity and bytes before/after the copy. Copies use independent storage; resolution rejects symlinks, hardlinks and paths outside the archive.

## Capture timing is part of the evidence

Capturing after native `ready` is still **after initialization**. Existing runs have no startup resource inventory against which to compare those bytes. The receipt always records `resources_verified_at_process_start: false`; it distinguishes `post_start_before_first_action` from `post_start_after_commands` using the command-receipt prefix at capture time. The prefix length/hash and command count are metadata; dynamic journals themselves are not archived or bound as immutable inputs.

Future process owners can request `before_start=True` after writing the native manifest and before spawning. A new capture records `stage: pre_start`; nonempty receipts or an existing process-output log reject that phase claim. This remains a caller-declared lifecycle boundary, not independent process tracing. The capture alone does not attest consumption, so `resources_verified_at_process_start` remains false. The process owner must select the detached resource roots before `Popen`, validate their identity, and record the actual launch association separately. Existing successful archives are never relabeled when a later caller passes a different phase flag.

## Detached resource consumption

`materialize_native_resources` first validates the archive, then recreates its exact logical resource inventory using regular blob copies. It does not open any captured original resource paths, and does not change the native directory's existing resource-root symlinks. Each resource file is a distinct inode, owned by the current user, with link count one. Empty directories and ordinary permission bits are preserved; setuid/setgid/sticky bits are not recreated. The private top-level materialization directory is newly owned storage.

Logical paths must be relative, normalized, beneath an allowed resource root, unique, and have a recorded directory parent. Absolute paths, traversal, duplicate entries, missing roots, file parents, extra output files, indirect files/directories and shared hardlinks are rejected. An existing `runtime-resources` tree is checked in full rather than repaired or overwritten. `materialization.json` binds its file hashes to the archive manifest. Failed attempts remain in separate `.resource-materialization-*` directories.

The native process owner may replace only the resource symlink roots it created itself with links to this detached tree. It must preserve unrelated files, handle startup cleanup, and associate the consumed tree with its launch receipt. When initializing from a patient XML, it should also compare the detached patient input to `patient_identity_input_sha256`; a change between earlier patient parsing and capture must fail. Explicit saved-state inputs are outside the resource-tree helper. Libraries and the ELF loader are archived but not rebound by this helper, and process environment/dynamic `dlopen` inputs remain outside its guarantee.

## Execution selection and relocation

`resolve_native_execution_inputs(root, native_directory)` validates the separate pre-Popen `execution-inputs.json`. Public `resolve_native_environment` invokes it automatically for every `pre_start` archive. Missing selection receipts, detached trees, copied states or patient files fail resolution. Historical post-start archives still resolve their original archive mapping, with no new startup claim or requirement for artifacts that never existed.

The selection receipt must match the actual native manifest, environment manifest, and resource materialization SHA-256 values. All nine native resource-root links must select their corresponding detached roots. Saved-state mode requires `command[2]` to name the independent `input-state.xml`, with matching state, patient-identity and selection hashes. Patient initialization instead requires `command[2] == '-'`, no state configuration, and the exact configured patient XML under the detached resource inventory. Selected regular files are owned, unshared inodes. Missing, changed, symlinked or hardlinked files fail, even when an outside file has matching bytes.

Recorded absolute paths are interpreted lexically, without opening their original locations. Relocation preserves the recorded `native_directory` workspace-relative layout beneath a new workspace root. The resolver maps native-local `runtime-resources` and `input-state.xml` identities to this relocated directory, validating retained resource-link text against the old or relocated native-local target. A relocation can alternatively recreate the root links using `runtime-resources/<root>` relative targets. Traversal, unrelated trees and live runtime links are rejected. This preserves evidence verification after moving an archive; it does not make its old executable command runnable in the new location.

The receipt records selected startup inputs, not independent observation of every native file read or proof of successful initialization. Native readiness and command receipts remain separate execution evidence. Executable/library rebinding, environment reconstruction and dynamic loader tracing remain outside this resolver's guarantee.

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
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_environment_evidence.py --held data/derived/systemic/exertion_v2/rest/native
```

The verifier uses real temporary files and an ELF executable. It checks detached survival after original library/resource mutation and executable deletion; rejection of corrupted blobs, external symlinks, hardlinks, and omitted startup dependencies; exclusion of state outputs; empty-directory retention; explicit capture timing; historical executable recovery; and preservation of failed capture attempts. API-missing, missing-dependency, hardlink, and historical-executable recovery tests were observed failing before implementation of the corresponding behavior.

The added prestart suite also tests malformed logical inventories, prestart rejection after an action, unchanged legacy phase metadata, idempotent materialization, and a separate process actually opening detached bytes after the live resource changes. The held-archive check relocates a copy of the real `exertion_v2/rest` archive into an isolated temporary workspace, materializes all 172 resource files, and has a separate consumer read `UCEDefs.conf`. The original run and archive remain unchanged. Its receipt is retained under `artifacts/verification/native-environment-prestart/held-resources.json`.

Execution-selection regression first failed on acceptance of a prestart archive without any selection receipt. It now covers both saved-state and patient-XML modes, all three receipt associations, input identity, missing/changed/indirect/shared files, a missing tree, unsafe selected paths, live resource-root substitution, and complete relocation with the original native directory hidden. The real fresh native session `data/derived/audits/session-fldvfext/actions` also passes selected-input resolution and relocation alongside the old post-start archive; neither original is rewritten.

An archive is evidence preservation, not a complete portable runtime or a bitwise replay guarantee. The legacy manifest does not identify arbitrary later `dlopen` modules, all process environment values, kernel/CPU behavior or native serializer omissions. Resource timing remains qualified above. Original libraries retain their licenses; local copying does not grant redistribution rights.
