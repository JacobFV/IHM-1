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
factory has not been changed by this increment; its owner must forward one pin to
both controller and cutaneous construction and retain it in the body receipt.

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
