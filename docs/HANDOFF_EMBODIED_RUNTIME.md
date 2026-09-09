# Handoff — embodied runtime, metabolic coupling and the workbench transport

Written 2026-09-08 at the end of a debugging session, for whoever wires up full
integration next. Everything below was measured in this workspace on this
revision; where a number is quoted it came out of an actual run, and the command
that produced it is given so you can re-derive it rather than trust it.

## 1. Read this first: the live body runs for about one second

`Start body` in the workbench (or `POST /api/embodied/sessions` then `/step`)
now advances roughly 51 intervals of 20 ms and then stops. Before this session
it advanced **one**. That is an improvement, not a working body, and it is the
binding constraint on any integration milestone that needs a body to run.

Three separate failures sit on top of each other. One is fixed, two are not.

| # | Failure | Where | State |
|---|---|---|---|
| 1 | Muscle metabolic decrement exceeds the native budget on step 2 | `ihm/assembly/embodied.py` exchange | **Fixed** in this session |
| 2 | Simbody integrator cannot advance past t≈1.026 s, upright | native mechanical plant | **Open**, diagnosed |
| 3 | Budget check returns at −10.43 W vs 10.37 W allowed, supine | reference `M0_w` is one instantaneous sample | **Open**, diagnosed, fix proposed but not applied |

Do not spend time rediscovering these. Sections 3–5 have the measurements.

## 2. What changed in the working tree, and what did not

**Nothing is committed.** All of this is working-tree only.

Changed by this session:

- `ihm/assembly/embodied.py` — the metabolic exchange lag (section 3).
- `scripts/verify_embodied_runtime.py` — regression coverage for it.
- `app/src/main.js`, `app/src/style.css`, `app/test/browser/body.spec.js` — the
  transport/playback UI (section 6).

**Changed by something else, concurrently, and still moving.** The tree was
clean at session start; these appeared during it and were still being written at
17:56, so another agent or editor is live in this repo right now:

```
ihm/assembly/peripheral.py            scripts/build_body_peripheral.py
ihm/assembly/reflexes.py              scripts/build_neuromechanical_experiments.py
ihm/assembly/garment_wardrobe.py      scripts/build_garment_wardrobe.py
docs/BODY_PERIPHERAL.md               scripts/fit_garments_to_envelope.py  (+404 lines)
scripts/verify_all.py
```

plus untracked `docs/IBM_IHM_NERVE_JOIN.md`, `docs/patches/`,
`docs/research/GARMENT_RIGID_ENVELOPE_GAP.md`, `ihm/assembly/garment_remesh.py`,
`scripts/enrich_peripheral_routes.py`, `scripts/peripheral_route_catalog.py`,
`scripts/verify_garment_remesh.py`, `scripts/verify_peripheral_join.py`.

A `git commit -a` here bundles that unrelated in-progress work with the coupling
fix. Commit by explicit path, or reconcile with whoever owns it first.

## 3. The metabolic coupling failure (fixed) — what it was

### Symptom

Every embodied session died on its **second** step, in every environment, with
or without an environment selection:

```
SIGNED_NATIVE_FAILURE: signed muscle decrement exceeds native discretionary budget:
requested_W=-43.860322  allowed_decrement_W=8.686126
```

`upright` −43.860 W, `supine` −36.489 W, `free` −36.110 W, against the same
8.686126 W allowance every time. `EnvironmentDynamics` was uninvolved.

### Mechanism

Two models disagreed about the size of a muscle-energy fluctuation.

**The plant side.** `ihm/assembly/embodied.py` computed the increment as a raw
instantaneous difference against a reference frozen at t=0:

```python
metabolic_w   = (energy - previous) / dt        # dt = 0.02 s
incremental_w = metabolic_w - self.reference_metabolic_w
```

That reference is a single instantaneous sample —
`scripts/native_mechanical_stream.cpp:105`, `metabolism.sample(state)` — of an
Umberger/Uchida probe (`scripts/native_muscle_metabolism.h`; whole-body basal is
explicitly off, and the constructor asserts `power[1] == 0` to keep it off).

**The native side.** `scripts/build_signed_muscle_variant.py` apportions muscle's
share of non-brain energy by its blood-flow fraction:

```cpp
const double discretionary = nonbrainNeededEnergy_kcal * fraction - mandatory
  + 0.8 * exerciseEnergyRequested_kcal + 0.5 * otherEnergyDemandAboveBasal_kcal;
```

and `scripts/native_signed_muscle_port.h:67-69` throws when a decrement would
push consumption below the obligatory amino-acid floor. That headroom is
**8.686 W**.

**The mismatch is in the fluctuation, not the offset.** The offset cancels
correctly — that is what the increment scheme is for. But the plant's rate is
far noisier than the whole native allowance. Bare plant, upright, no actuation,
100 steps of 20 ms, reference `M0_w` = 99.103 W:

```
step   metabolic_W   signed_work_W    heat_W
   1       105.571          0.494    105.077
  11        78.873        -75.414    154.287
  51       137.427       -183.120    320.547
 100       105.135          4.277    100.858
```

It swings 77–137 W and never settles. Downward excursions of 20–40 W are
ordinary — three to five times the entire budget.

Step 1 survived by luck of sign: its increment was **+0.543 W**, and positive
increments are not bounded by that check. Step 2 was −43.860 W.

### Why the test suite missed it

`scripts/verify_embodied_runtime.py`'s stub plant returned a **constant** rate
(105 W against a 100 W reference), so its increment was always +5 W and never
negative. The guard was never exercised in the failing direction.

The real multi-step checks did catch it. `scripts/verify_environment_embodied.py`
failed identically, and its leftover scratch directories under `data/derived/`
recorded `passed: false` from **2026-09-07 16:54** onward — predating the
environment-dynamics commit `cba3d7d`.

### Statistics that justify the fix

Bare plant, 250 steps of 20 ms (5 s), increments relative to `M0_w`:

| Environment | per-step min | max | mean | steps below −8.686 W |
|---|---|---|---|---|
| supine | −7.31 W | +29.01 W | +1.16 W | 0 / 250 |
| upright | −23.68 W | +38.32 W | +3.17 W | 19 / 250 |

Averaged over a **1.0 s** window, every excursion fits the budget in both
(supine −0.34…+5.61 W, upright −0.84…+6.83 W). The signal is fast fluctuation
around a small positive mean, not a sustained decrement. That is the entire
basis for the fix.

Reproduce: the scripts are gone with the scratchpad, but the shape is
`ArticulatedBodyPlant(root, fresh_output, environment=...)`, then
`advance(.02)` in a loop differencing `muscle_metabolic_energy_j`.

### The fix

`ihm/assembly/embodied.py`, constant at line 30, state at 374, applied at 496,
exchanged at 510, committed at 529:

```python
METABOLIC_EXCHANGE_TAU_S = 2.

blend      = 1. - math.exp(-dt / self.metabolic_exchange_tau_s)
exchanged_m = filter['m_w'] + blend * (incremental_w - filter['m_w'])
exchanged_h = filter['h_w'] + blend * (delta_h       - filter['h_w'])
exchanged_w = exchanged_m - exchanged_h
```

Properties that matter, and that the tests pin:

- **The ledger identity survives.** The port validates
  `delta_m == delta_h + delta_w`. Both channels use one shared coefficient and
  work is *derived* as their difference, so the identity holds by construction
  rather than by tolerance.
- **Energy is delayed, never created or destroyed.** The lag's steady-state
  output is the mean of its input.
- **It is not a mask.** A decrement that persists still arrives in full, and is
  still rejected if it genuinely exceeds the budget — which is exactly what
  failure #3 below is.
- **Rollback-safe.** The filter commits only after the whole interval succeeds,
  so a rejected step leaves it where the last accepted one was.
- **It is an uncalibrated engineering choice**, disclosed in the frame under
  `coupling.metabolic_exchange_basis`, with the raw increments still reported
  beside the exchanged ones. τ = 2 s was chosen because 1 s windows already fit
  the budget in the table above. It is not an identified substrate kinetic and
  must not be described as one.

### Verification after the fix

- `scripts/verify_environment_embodied.py` — **passes** (`duration_s` 0.06,
  `contact_count` 121, wall 13.3 s). It had been failing since 2026-09-07.
- `python -m unittest scripts.verify_embodied_runtime` — 13/13, including two new
  tests: one drives a −44 W excursion through the exchange and asserts it lands
  inside the budget *and* that a sustained decrement still arrives in full;
  the other asserts a rejected interval leaves the lag untouched.
- `cd app && npm test` — 100/100.

## 4. Open failure: upright integrator, t≈1.026 s

```
Integrator step failed at time 1.025899526275925
Unable to advance time past 1.0259.  (Required condition 't1 > t0' was not met.)
```

Raised from `ihm/native/mechanical_stream.py` via `ArticulatedBodyPlant.advance`.
**Independent of metabolism** — the bare unactuated plant shows its violent
precursor at the same moment (step 51 above: −183 W signed fiber work, +320 W
heat, both at t = 1.02 s). The upright body has no balance control and appears
to be collapsing into a stiff configuration the integrator cannot step through.

This is what now caps every upright run at ~1 s. It is the first thing to fix if
integration needs a body that stands.

## 5. Open failure: supine budget returns at −10.43 W

After the lag, a long supine run eventually reaches:

```
requested_W=-10.425061  allowed_decrement_W=10.368900
```

Note both numbers moved: the request is a quarter of what it was, and the
allowance itself is dynamic. This is the lag working correctly and reporting a
**sustained** discrepancy rather than noise — the actuated body settles roughly
10 W below the `M0_w` sampled at t = 0.

**Proposed fix, deliberately not applied:** re-base `M0_w` on a short settling
window at initialization instead of one instantaneous sample. That would centre
the increments near zero and remove the sustained breach. It changes what the
metabolic reference *means*, so it wanted a human decision rather than mine.

Two alternatives were considered and rejected:

- *Clamping the increment to the allowance* — would require clamping all three
  channels together to preserve the identity, which means inventing a heat/work
  split. Fabrication; do not do this.
- *Raising τ* — delays the arrival but cannot fix a sustained mean. The lag
  converges to the true mean by design.

## 6. The workbench transport UI

The original complaint: two play-like controls, and no way to tell which ran the
simulation. They were `Start body` (left column → Simulation) which runs the
live body, and a bottom-centre `▶` which replays a stored trajectory from
`/api/body/trajectory` and is disabled whenever the body is live
(`app/src/main.js`, `$("play").disabled = live || frames < 2`).

Now the bottom-centre slot has one tenant at a time, via `syncTransport()`:

- Labelled `RECORDING ▶ 1x`, so it reads as a scrubber, not a run control.
- While live it is **replaced** by a `● Live · 3.42 s computed` pill (dot pulses
  when running, steady and "Paused" when paused) rather than greyed out.
- With no recording it reads `NO RECORDING` and the `▶` carries a tooltip
  saying why.
- The pane is "Recording playback" and states that the simulation starts from
  the Simulation section. The recorded clock reads `Recorded · 1.000 s`.

One latent bug fixed on the way: `#transport` and `#live-readout` both carry
`display: flex`, which overrides the `hidden` attribute. Added
`#transport[hidden], #live-readout[hidden] { display: none; }`, matching the
convention already used eight times in that stylesheet.

## 7. Operational notes that will otherwise cost you an hour

- **One native body at a time.** `POST /api/embodied/sessions` refuses while any
  session is unclosed — `One native body at a time while sharing machine
  resources`. Always `POST /api/embodied/sessions/<id>/close`, including after a
  failed or interrupted run. List with `GET /api/embodied/sessions` and look for
  `closed: false`.
- **Browser tests need the dev server, not the production one.** Run
  `cd app && npm run dev` (Vite, `:5173`, proxies `/api` to 8765) and then
  `npx playwright test`. Pointing `APP_URL` at `:8765` makes the suite crawl and
  emit spurious failures — that is a harness misconfiguration, not a regression.
- **Session scratch is large.** Each embodied session writes ~140 MB under
  `data/derived/embodied-sessions/<uuid>/`, and failed verification runs leave
  ~124 MB dirs like `data/derived/environment-embodied-*`. Both are gitignored
  and safe to delete when no session is open. 1.1 GB of them was cleared during
  this session.
- **Servers left running:** workbench on `:8765`
  (`.venv/bin/python -m ihm serve --port 8765`) and Vite on `:5173`.

## 8. The scene ball, for context

The original question was why the ball in "Garden patio" does not fall. It does
— it was frozen because the body died on step 2. Verified in isolation
(`EnvironmentDynamics`, garden-patio, stub registration): `ball-large-1` falls
from y = 0.18 m to rest at −0.851 m (plane −0.96 + radius 0.11) by t ≈ 1.2 s.

Expect a drop and settle, **not a bounce**. Ground contact is a damped penalty
force (stiffness 1500, damping 8, `ihm/assembly/environment_dynamics.py`), and
the other engine path is explicitly restitution-free —
`ihm/assembly/interactive_scene.py:194`,
`jn=max(0.,-vn*self.mass)  # zero restitution, no manufactured bounce`.
In garden-patio only `ball-large` is simulated; table, chair and trees are
display-only, because the scene has no box collider.

## 9. Open decisions

1. Re-base `M0_w` on a settling window (section 5)? Unblocks sustained supine runs.
2. Chase the upright integrator failure (section 4)? Unblocks upright runs at all.
3. Is τ = 2 s acceptable as a disclosed engineering choice, or should the
   exchange interval itself be reconsidered? The 50 Hz exchange is what forces a
   filter to exist.
