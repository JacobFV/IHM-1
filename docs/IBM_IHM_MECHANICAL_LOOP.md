# IBM–IHM reduced mechanical loop

The IBM checkout's `scripts/embody.py` now loads IHM's `MechanicalPeripheral`
(`ihm/assembly/mechanical_peripheral.py`). This couples the existing reduced SI
`BodyMechanics` to `BodyPeripheral`; it is not the native articulated runtime
and is not biologically validated. No kernel weights were trained or changed.

Each interval advances mechanics using the preceding interval's peripheral
motor activations, passes its transforms, path lengths and forces to the
peripheral, and retains the new motor activations for the next interval. IBM
feeds the returned local spindle rate into its cord; the cord output becomes
the next descending command. Peripheral motor transport and activation filtering
remain active. All 249 peripheral muscle IDs must exist among the mechanics'
964 muscle/connective paths. The remaining paths have no IBM command.

`proprioceptor_rates_hz` remains the legacy mixed length/force signal. New
`spindle_rates_hz` and `tendon_rates_hz` separate those components, keyed by bare
muscle ID. IBM uses strict spindle indexing and verifies the complete key set.
These are local transduced rates, not the delayed `brain_inputs_hz` arrival
stream. IBM still uses the 64 cutaneous channels for its cortical input.
Mechanics now returns actual affine-anchor path lengths, avoiding reconstruction
from translations alone. Reference stretch below 1e-12 is roundoff-suppressed.

## Measurements

20 steps, dt=10 ms, pressure stimulation, seed 0, existing 34-source checkpoint:

| Arm | SD of command mean over time | Peak stretch arc | Peak spindle Hz |
|---|---:|---:|---:|
| Full | 0.0252852 | 0.0176576 | 6.14495 |
| No cord | 0.0000411 | 0 | 6.27681 |
| Severed kernel | 0.0254199 | 0.0159502 | 5.54192 |

Raw logs: [measurements](measurements/ibm_ihm_mechanics/). The stretch arc now
receives mechanical length feedback. Similar full/severed temporal variation
still does not establish useful cortical control; these runs are a causal wire
check, not learned movement or a task-performance comparison.

A paired perturbation/control test at dt=1 ms prescribes 5% axial deformation
of recfem_r in the mechanics, with initially zero descending command. Local
spindle response first appears at 1 ms; stretch and alpha differences first
appear at 31 ms, exactly 30 ms after that input, and muscle activation differs
at 46 ms. Thus the cord delay is 30 ms, but boundary-to-muscle latency is longer:
it includes peripheral motor transport and explicit coupling. No cortical
latency is measured by this test. A separate activation/no-activation comparison
checks that muscle commands cause deformation, translation and stretch feedback;
a force-only check ensures force cannot masquerade as spindle input.

IBM currently reports **2/249 muscles mapped to segments** with these canonical
IDs. Its per-muscle arcs operate even for unmapped IDs. This is a separate naming/
coverage issue; these results do not validate anatomical segment routing.
Reciprocal and autogenic arcs remain unfed in `embody.py`; logs expose every arc
by name. The new tendon signal is available for a future explicit Ib contract.

## Reproduce

From IHM:

```sh
.venv/bin/python scripts/verify_sensorimotor_loop.py
.venv/bin/python scripts/verify_body_peripheral.py
.venv/bin/python scripts/verify_body_mechanics.py
```

From IBM (repeat with `--no-cord` and `--sever`, distinct output paths):

```sh
PYTHONPATH=. .venv/bin/python scripts/embody.py --steps 20 --stimulate --seed 0 --out /tmp/full.json
```

No training objective or demonstration corpus is supplied by this handoff.
Training remains subsequent work, with a task metric and matched severed control
needed before attributing learned movement to the kernel. The pre-existing
peripheral/nerve edits remain uncommitted; this work does not commit them or
change the running browser/server processes.
