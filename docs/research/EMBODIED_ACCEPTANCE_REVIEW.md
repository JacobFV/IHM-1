# Embodied acceptance review, 2026-09-05

Read-only adversarial review against UI commit `441f63583d52ad9413c38e8a24b5007e5245209a`.
No native engine, browser, production build, or acquisition ran. Root and other
agents are actively integrating these files; findings describe these exact bytes:

| File | SHA-256 |
| --- | --- |
| `ihm/app/embodied.py` | `a15a95ee3b0a953f112a8e4b70bad19d920703f7b2333b20956fdb13bf159e0e` |
| `ihm/assembly/embodied.py` | `6ce4643b2380922605e9ada17c554d6b8107d8e1886d18202f87d0b9e1eb2759` |
| `ihm/app/__init__.py` | `57878d0560c6812dfc2453a5cad3174e0cde864ab724f675dd7b6954a25b1810` |

## 1. Cleanup failure releases the resource slot — reproduced, blocking

`ihm/app/embodied.py:55–59` suppresses `body.close()` failure and unconditionally
sets `actor.closed=True`. The create gate at line 100 treats that flag as proof
that all resources are released. Lines 66 and 120–121 then hide the cleanup
failure and remove the actor on a repeated close request. A cleanup exception
can therefore leave a child alive while the service accepts another native body.
This is especially material under the machine's resource constraints.

Smallest reproduction, with a simulated child and no subprocess:

```python
from pathlib import Path
from tempfile import TemporaryDirectory
from ihm.app.embodied import BodyActor, EmbodiedSessions

class Body:
    failed = False
    alive = True
    def snapshot(self):
        return {'sequence': 0, 'time_s': 0, 'entities': {}}
    def close(self):
        raise RuntimeError('child termination unconfirmed')

with TemporaryDirectory() as directory:
    body = Body()
    actor = BodyActor(lambda: body, Path(directory) / 'actor')
    actor.ready.result(2)
    try:
        actor.request_close()
    except RuntimeError:
        pass
    actor.thread.join(2)
    sessions = EmbodiedSessions(directory)
    sessions.actors['a' * 32] = actor
    print(body.alive, actor.status())
    print(not any(not a.closed for a in sessions.actors.values()))
    print(sessions.command('a' * 32, 'close', {}))
```

Observed: simulated child alive `True`; status `closed=True, error=None`;
resource gate allows create `True`; retry close reports success. Both cleanup
attempts raised. This tests the actor contract, not an actual OS kill failure.

Smallest correction: distinguish termination confirmation from worker-thread
exit. Preserve cleanup errors and retain the resource slot on unconfirmed
termination. Adapter close methods must expose actual child termination rather
than a flag set before `wait()` succeeds. A retry needs an owner-thread cleanup
path or explicit supervisor, not deletion of the failed actor.

Related source-only shutdown issue: `EmbodiedSessions.close()` (line 125) merely
requests cancellation for an initializing actor and returns. `Server.server_close`
at `ihm/app/__init__.py:114–117` does not await those workers; they are daemon
threads. Server process exit can interrupt initialization cleanup. The pending
cancellation test establishes eventual cleanup while the process stays alive,
not cleanup before server exit. Shutdown should await termination with a bounded,
explicit escalation policy. No actual process-exit experiment was run.

## 2. Default native inputs reintroduce the audited thermal mismatch — source-proven, blocking promotion

`ihm/assembly/embodied.py:47,64` selects `native_baseline_v1` and
`whole_body_integrity_depletion`. The baseline was initialized under
`saturation_bounds_heatflux`; the retained thermal audit documents its initial
storage loss and the still-defective later boundary. See
`LONG_HORIZON_THERMAL_AUDIT.md:25–39,103–107`. This is a known constitutive/state
mismatch, not simply an unspecified calibration gain.

Smallest reproduction is inspection of the factory's selected state/variant
and the referenced retained thermal receipts. No new physiological trajectory
was run. Select an explicitly paired, frozen state and accepted library identity
as a startup configuration; preserve old variants. Do not promote a final RH
library as globally accepted merely because its thermal correction passed:
separate glucose/acid–base defects and coupled-load acceptance still matter.
Expose the selected state/library identity and research limitations to the UI.

## 3. The factory does not yet select the upper-body or breathing integrations — source-proven

`ihm/assembly/embodied.py:69` omits `augmented_registration`, selecting the
80-muscle lower-limb model. The plant's optional bilateral Arm26 model is not
selected by this call. The UI correctly derives available muscle keys from the
frame and will show only 80 here. Source-native acceptance of the 92-muscle model
must precede selecting it; the first Arm26 native load failure is retained by
the mechanics agent and is not hidden by this review.

At lines 156–158 and 180–182, `entities` is the mechanical frame directly and
there is no `respiration.skin_field`. Thus the newly available respiratory
projection cannot move the live skin/chest in this revision. The renderer at
`app/src/main.js:932` consumes that field correctly when present. Integrate the
respiratory composition once, with the agreed canonical-reference skin field
before the entity transform. Preserve `centroid_m - translation_m` as the
reference centroid. The respiratory helper's author confirmed that invariant.

These are missing factory bindings, not evidence that native muscle or
respiratory dynamics have been validated or that UI fallback is acceptable.

## 4. Hidden 120-second horizon terminates a normal next tick — source-proven

The factory at `ihm/assembly/embodied.py:64` hard-codes `horizon_s=120`.
`NativeSession.step()` at `ihm/native/session.py:258–262` rejects an interval
beyond the remaining horizon. `EmbodiedRuntime.step()` sends respiratory load,
exercise, and optional compression first (lines 145–149); it has already marked
the native owner touched. At the next requested interval after 120 s, the
runtime therefore aborts the session rather than reporting a normal stop.

Smallest fixture: place native elapsed ticks at 6000 and request `.02` s; the
native step precondition raises. The ordering and resulting abort are established
from source; no 6001-step run was performed. The UI continuously requests `.02`
while running and has no horizon metadata/end status.

Validate remaining horizon before touching any owner. Include remaining time
and an explicit completed state, or deliberately select a longer bounded horizon.
A horizon limit should neither silently extend execution nor look like an
uncertain solver failure.

## 5. Exact-input archive is incomplete for the neural/mechanical execution — source-proven retention gap

`ihm/assembly/embodied.py:48–61` copies a fixed Python module list and six JSON
assets. `brain.json` references the executed pinned law at
`data/derived/canonical/brain-sources/ibm-neural.py` (SHA
`f8c7d3172af71b399d9e6301f068d8fde8adf3363941ef0d3d3ec50ef62f663b`), but those
referenced source bytes are not copied into the run. `BodyBrain` verifies and
loads the current preserved source, so this is a future reproducibility gap,
not a demonstrated wrong law in current execution. The lazily imported
`body_microstructure` module can also be absent from the initial receipt list:
its first import occurs inside `NativeTissueExchange.from_workspace()` at
`body_exchange.py:66`, after capture. Its path hash later appears in tissue
frames, but the exact executable bytes are not thereby archived.

`NativeMechanicalStream` retains the model/path/motion inputs and a build
manifest (`mechanical_stream.py:38–57`). Its mechanical shared libraries are
hashed at original runtime paths by the build helper, not detached or copied
into an execution-input archive. Transitive/system dependencies are not an
actual `ldd` closure in that manifest. Future replacement can make a recorded
run non-reconstructible. Physiology's `NativeSession` now has a stronger
before-start environment archive and detached resources; it should not be
described as missing those protections.

Smallest correction: enumerate dependencies after required imports, include
executed pinned-law bytes and exact mechanical runtime dependency closure,
bind copies to the manifest, and execute from selected immutable inputs where
possible. Keep historical source paths as provenance. Do not resolve them by
reading or modifying the live IBM checkout.

## Checks and non-findings

Executed under `OPENBLAS_NUM_THREADS=1`, `prlimit --as=1073741824`, `nice -n 10`:

```sh
.venv/bin/python -m unittest scripts.verify_embodied_actor scripts.verify_embodied_runtime
```

All five tests passed in 0.008 s. The separate fake cleanup reproduction above
confirmed finding 1. An initial attempt at that fixture omitted the OpenBLAS
thread variable and exited under the address-space cap before running; the
single-thread rerun produced the reported result. No native process was started.

Strict integer command sequence, single owner thread, recorded command before
advance, and current-snapshot recovery match the committed UI. The UI does not
resubmit an uncertain force impulse. Successful local checkpoint release and
pre-native rollback are present. Post-native uncertainty intentionally aborts;
no false whole-body rollback is claimed. Accepted frame shape and physiology
metadata match `LiveBodyHistory` and panels; unavailable units remain explicit.
These checks do not establish native 92-muscle startup, mechanical convergence,
long-horizon homeostasis, or coupled whole-body energy closure.
