# Regional skin fluid and albumin transport

I9 now has an independently executable, body-attached regional transport model,
not a claim of coupled native forearm physiology. `RegionalSkinTransport` owns
three private storage compartments (blood proxy, interstitium, lymph) and two
finite external reservoirs. `native_blood_storage_connected=False` is enforced
on restart. Native full-skin state remains owned by BioGears.

Run `.venv/bin/python scripts/verify_body_transport.py` for checks and
`.venv/bin/python scripts/build_body_transport.py` to regenerate
`data/derived/canonical/skin-transport.json`. The artifact records the initial
parameter card, checkpoint states for actual intervention runs, exact input and
implementation hashes, fixed seed 417 (no stochastic steps), dimensional checks,
per-compartment integrated balances, and timestep refinement.

## Sources and scope of transfer

The pinned `data/derived/coupling/native-skin-circuit.json` supplies SI node
pressures, volumes, resistances and compliances. Its source native XML hash is
verified before use; the same XML supplies actual albumin concentrations and
SkinTissue mass (3.3 kg). No concentration is guessed from oncotic pressure.

The checked-out BioGears implementation is the authoritative equation source:

* `.../src/engine/Systems/Tissue.cpp`, `CalculateOncoticPressure`, lines
  2075–2123: total protein = 1.6 albumin; Landis–Pappenheimer oncotic pressure.
* `.../src/engine/Systems/Diffusion.cpp`,
  `CalculateMacromoleculeDiffusion`, lines 578–632: albumin permeability
  coefficient 0.03306 mL/min/kg, small-pore reflection 0.954, Peclet law,
  and nonnegative-only blood-to-tissue mass transfer. The unused large-pore
  reflection 0.097 is deliberately not turned into an invented second pathway.
* The following `CalculatePassiveLymphDiffusion` function: donor-concentration
  convection from tissue to lymph and lymph to vena cava, positive-only.
* `.../src/engine/Controller/BioGears.cpp`, Skin initialization around lines
  4440–4520: extracellular compliance, filtration resistance, lymph driving
  pressure, uptake resistance and one-way valve topology.

The public [BioGears tissue methodology](https://www.biogearsengine.com/documentation/_tissue_methodology.html)
documents Starling exchange and the protein oncotic approximation; its older
statement that albumin transport is future work is superseded by the pinned
source implementation above. [BioGears sepsis model paper](https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2019.01321/full)
provides primary published context for albumin/lymph transport in the engine.
These sources support the model family, not regional patient calibration.

## Geometry and parameter card

The region is the existing canonical `body-detail-microvascular-left-0` unit,
left forearm skin, with its original `body-bp3d-FJ2810` face 192812,
barycentric coordinates and 1 mm normal depth. This is an actual material
attachment to the canonical body, though the microvascular tree is synthetic
and the supplying native terminal vessel is unidentified.

The sum of cylindrical edge lumen volumes is 2.785973275e-11 m³. Dividing by
native Skin1 volume gives extensive scale `s=1.646040731e-7`. Storage, tissue
mass, compliance and permeability-surface product multiply by s; resistance
divides by s; pressures and concentrations retain native values. This is a
homogeneous extensive-scaling prior, not a measured forearm perfusion law.
In particular, the canonical tree's separate steady pressure solution is not
silently treated as the pressure solution of this storage model.

Lymph volume/compliance and its central return resistance are also scaled by s.
This whole-body-lymph allocation is an explicit uncalibrated prior. The initial
regional lymph volume exceeds the interstitial volume because those source
compartments describe different native anatomical scopes. No physiological
interpretation of that allocation is validated. A measured regional lymph
volume, compliance and transit dataset would be needed to replace it.

External reservoirs each start at 1000 times the sum of regional storage
volumes, with measured native arterial-proxy/vena-cava albumin concentrations.
They have prescribed pressures and dynamically conserved contents. This large
finite-reservoir initialization is a boundary prior, not a body inventory.

## Equations and numerics

All runtime quantities are SI: m³, Pa, seconds and kg of albumin. At each
storage, `P=P0+(V−V0)/C`. For an ordinary hydraulic connection,
`Q=(Psource−Ptarget)/R`. The filtration law is
`Qf=(Pb−Pi−πb+πi)/Rf`, with total protein T in g/dL and
`π=133.322387415*(2.1T+0.16T²+0.009T³)` Pa. Albumin in kg/m³ converts to
g/dL by division by 10. The native fluid oncotic reflection is 1; the separate
macromolecular reflection is 0.954, preserving the distinction in native code.

For albumin, `PS=0.03306*1e-6/60*tissue_mass_kg` m³/s,
`Pe=Qf*(1−0.954)/PS`, and
`J=max(0, PS*Pe*(cb−ci*exp(−Pe))/(1−exp(−Pe)))` kg/s.
The zero-Pe limit is `max(0,PS*(cb−ci))`; a Taylor expansion and `expm1`
avoid native-form cancellation/overflow at small or negative Pe. Reverse
fluid filtration is allowed but albumin back-diffusion is suppressed exactly
as in the source algorithm; this is a material model limitation.

Lymph uptake is `max(0,Pi+drive−Pl)/Ruptake`; return is
`max(0,Pl−Pv)/Rreturn`. Both carry donor albumin concentration and are
multiplied by `1−obstruction`. Thus complete obstruction blocks uptake and
return, not merely the terminal duct. This is the defined experimental control.
All flow/mass transfers enter one storage and leave another. Capacitance is
storage, not a fictitious transport path to ground.

The explicit finite-volume integrator uses at most 1 s substeps and further
reduces the substep to limit aggregate donor volume and mass loss to 10%.
Simultaneous donor accounting preserves positivity and conservation; no
post-step clipping creates or destroys matter. A reservoir-exhaustion/stiffness
guard raises instead of silently continuing an invalid experiment. Pressure
boundary work is external; there is no energy/temperature balance claim.

## Recorded verification

Tests first failed on the absent implementation, then passed after implementation.
Actual 600 s experiments at identical initial state give interstitial changes:

| Experiment | Interstitial volume change |
|---|---:|
| Baseline | −0.1136% |
| Vena cava boundary +600 Pa | +0.2025% |
| Arterial conductance/inflow fraction 0.25 | −0.4104% |
| Complete lymph uptake and return obstruction | +0.3598% |

The baseline is a transient native operating-point transfer, not a new steady
state; other-organ lymph influx and intracellular fluid exchange are excluded.
Do not interpret these changes as validated edema predictions. The obstruction
and raised-venous scenarios exceed the baseline interstitial volume; reduced
supply lowers integrated arterial inflow. One-way return is explicitly checked
under adverse pressure. A 3600 s stress run combines +5000 Pa venous pressure,
zero arterial inflow and complete obstruction and retains positive storage.

At 60 s under +600 Pa, dt = 1, 0.5, 0.25 s gives interstitial volumes
2.139608862e-10, 2.139608972e-10, 2.139609030e-10 m³; successive-difference
ratio 1.914 is consistent with first-order convergence. Restart through JSON is
bitwise identical for the same step schedule. Local volume/albumin ledger
residuals are around 1e-24 SI in the stress run; full-reservoir residuals are
4.6e-22 m³ and 7.3e-20 kg, with the larger reservoirs' floating-point summation
error reported separately. Each node's integrated balance is checked.

Not included: native-state feedback, intracellular exchange, sweating,
inflammation-dependent permeability, active pump adaptation, interstitial
mechanics, blood-cell partitioning, nonlinear regional lymph compliance,
or independently validated regional calibration.
