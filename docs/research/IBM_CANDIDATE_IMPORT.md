# Immutable IBM candidates and shared runtime source identity

`capture_ibm_candidate.py` now exports an exact committed `ibm/` tree through
`git archive` into a fresh directory. It accepts ordinary Python source only,
excludes working-tree edits, datasets, training scripts and checkpoints, and
refuses any existing destination or destination inside the donor checkout.
It records donor status, exact commit, per-file and package hashes, regional neural
source hash, selected function AST hashes and the delta from the active regional
source. No active artifact is overwritten or promoted.

```sh
PYTHONPATH=. .venv/bin/python scripts/capture_ibm_candidate.py \
  --revision COMMIT_SHA --output data/derived/ibm-candidates/COMMIT_SHA
```

The required explicit output is intentional. A candidate directory is write-once
through this API; every import independently verifies its manifest and source
bytes, rejecting edits or extra/missing files. It is not an OS-enforced read-only
filesystem. `SourcePin` pins the manifest, package and neural source together.

For an opt-in fresh process:

```python
from ihm.brain.candidate import SourcePin
from ihm.brain.ibm_backend import IBMBackend
from ihm.assembly.sensorimotor import SensorimotorController
from ihm.assembly.cutaneous_feedback import CutaneousFeedback

pin = SourcePin.load(candidate_directory / 'source_pin.json')
backend = IBMBackend(source_pin=pin)
controller = SensorimotorController.from_root(root, source_pin=pin)
skin = CutaneousFeedback(root, sites=explicit_sites,
                        recruitment_hz_per_response=explicit_gain,
                        source_pin=pin)
```

`BodyBrain(..., source_pin=pin)` and `BodyBrain.from_dict` also support this pin.
The candidate regional brain calls the actual functions in the same verified
`ibm.processes.neural` module used by the package runtime. No donor process formula
is copied into IHM. Existing canonical geometry and its engineering parameter
assignments remain unchanged. Brain outputs carry `source_identity`, cutaneous
materializations carry the pin, and sensorimotor checkpoint identity includes the
brain source identity. Cross-law controller checkpoint restoration is rejected.

Omitting the pin preserves the existing default package and independently retained
legacy regional law. Explicit candidates do not silently replace a live process:
`SnapshotLoader` rejects a different IBM package identity once any IBM modules are
loaded. Start a new process/session to select another candidate. The embodied
factory accepts `EmbodiedRuntime.from_workspace(..., source_pin=pin)` and forwards
the identical pin to controller and cutaneous construction. It verifies and loads
the selected package before source receipts, output creation or native startup.

## Evidence policy and bounded acceptance

Candidates currently use **compiled source priors only**. The importer does not
fetch optional measured-summary JSON or raw data. It rejects an enclosing donor
`data/sources` + `ibm` evidence root, preventing implicit evidence loading from
another checkout. Candidate tactile materialization explicitly disables the new
automatic structural substrate fallback. Tests confirm the current microcircuit,
human census, kinetics and vascular summary objects report `from_evidence=False`.
Their compiled source values retain the source's own evidence descriptions; they
are not claimed as measurements acquired by this import.

```sh
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python scripts/verify_ibm_candidate.py
```

The bounded tests cover exact commit versus dirty working tree, overwrite refusal,
changed/extra source rejection, pin mismatch, shared package/regional/receptor
identity, cross-law checkpoint rejection and hot-swap refusal. A separate fresh
process exercises actual rapid/slow source parity, 64-sample periodic receptor
windows, causal ZOH responses and one regional brain interval. No full graph,
substrate, native body or training run is performed.

Retained candidate:
`data/derived/ibm-candidates/398375cc993ac17907d4ddd9f6d53eed2fab31de/`.
Its package hash is
`f969da521ef87ece047da189de17118538a650ebab65605a3254a6be38dffa60`;
regional neural source hash is
`8a08d9d7f43999d9ae312af817e67edd39a292edcb7e575f49be43b20b326a79`.
The sibling `398375cc993ac17907d4ddd9f6d53eed2fab31de-acceptance.json` records
shared identity and the intentional inhibition-law probe delta: -160 to
-133.33333333333334 mV/s at the explicit fixture state. The candidate manifest names
`shunting_inhibition_rate` as the changed selected symbol. This is an accepted
bounded compatibility result, not biological calibration or approval to replace
existing body sessions. Broader coupled behavior and new vascular/microcircuit
materialization remain separate work.

## Body factory retention

The candidate directory must resolve inside the factory's configured workspace.
This is an explicit criterion for root-relative receipt paths; direct standalone
`IBMBackend` usage can still use a candidate elsewhere. Invalid pins, changed
manifests/source, a mismatched `source_pin.json`, or a different already-loaded IBM
identity fail before native owners are created.

Body `output/inputs` retains the exact candidate manifest, pin JSON and every
candidate Python source file at its original workspace-relative path. The final
body manifest contains `brain_source_pin`, `brain_candidate_loaded_modules` with
per-module source hashes, and the ordinary source/loaded-code receipts for the
candidate and snapshot-loader implementation. Sources are checked again after
factory initialization. Replaying the archive in a relocated workspace requires
an explicit rebased pin; original absolute paths are retained as provenance.

`verify_candidate_embodied_factory.py` uses fresh Python processes and fake native
owners to check preflight rejection, single-pin forwarding, full source retention
and preservation of optional regional-native selection. An additional lightweight
check verified 102 candidate files, 72 loaded IBM module receipts and actual IHM
loaded-code attestation. `verify_regional_embodied_factory.py` default/regional
routing remains passing. Actual native factory acceptance is recorded below. App configuration and default
promotion remain unchanged.

## Actual candidate factory acceptance

`data/derived/audits/candidate-regional-factory-o81_fuoo/verification.json` records
a successful real `regional_skin=True` factory initialization with candidate
`398375cc...` and one native skin sensor. The fixture reuses `skin-contact-0`, its
cortical assignment, gain 0.1 and reference temperature 33 C from the previously
retained cutaneous factory receipt; these remain explicitly engineering priors.
The native material manifest/quadrature/triangle identity and transformed current
coordinate frame match the receptor binding.

All clocks remain zero. Shared brain/receptor identity, 102 archived candidate
files and 72 loaded IBM module receipts were verified. Initialization, snapshot
and cleanup took 2.067 s; both native processes were reaped. The fixture uses a
30 s alarm, parent 2 GiB soft / 4 GiB hard address-space limit, nice 10 and one
BLAS/OMP thread; child launchers retain their existing 4 GiB hard ceilings.

Run only in a coordinated native slot:

```sh
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 \
  .venv/bin/python scripts/verify_candidate_embodied_factory.py --native
```

This establishes source-bound native initialization and cleanup, not a coupled
step, receptor response under native motion, equilibrium, full brain completion
or promotion of the candidate to the default runtime.

### Opt-in embodied API selection

`POST /api/embodied/sessions` accepts an optional configuration member:

```json
{"regional_skin":true,"ibm_candidate":{"commit":"398375cc993ac17907d4ddd9f6d53eed2fab31de","manifest_sha256":"397b3fd3e20b7a15da0055bbba9c83fec6cfaee04135e34c99407da8280d300b"}}
```

The configured workspace's `data/derived/ibm-candidates/<commit>` directory is the server-owned candidate registry. Installation remains a server/operator action using the capture helper. The API accepts no filesystem path, donor location, revision alias, upload, or dataset reference. Exact lowercase commit and manifest hashes select an installed artifact; all SourcePin package, neural-law, manifest and source hashes are verified before actor creation and again on its owner thread before calling the factory. The factory's existing immutable import identity checks still reject process hot-swapping before native startup. Selecting a different package in a process that already loaded IBM requires a fresh server process.

Resolution rejects redirected pin paths, symlink ancestors or tree entries, missing/changed/unmanifested source files, extra root files and even unexpected empty directories. Candidate storage must remain server-controlled during initialization; these checks are not a filesystem locking protocol against a concurrent privileged operator changing directories. The donor working tree and optional evidence datasets are never imported by selection.

Creation, listing, snapshot and command responses include `brain_source_selection` for opt-in sessions: mode, commit, manifest, package and neural source hashes. This is the verified selection, not a claim that initialization has finished; clients must also inspect session status. The runtime manifest separately retains the actual import receipt and source bytes. Omitting `ibm_candidate` retains legacy factory defaults and response behavior. No candidate is automatically promoted.

Validation: `PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m unittest scripts.verify_embodied_candidate_api scripts.verify_embodied_actor scripts.verify_embodied_http` exercises candidate forwarding/identity, malformed selectors, tampering, redirected and symlink paths, extra entries, owner-thread revalidation, default preservation, actor lifecycle and HTTP routing. This suite mocks the native factory; it makes no new native initialization or stepping claim.
