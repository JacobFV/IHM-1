# Executable coupled cervical/thoracic owner

The source-only owner now represents the complete **27.654676965 kg** original
native torso partition. The nine cervical bodies and 48 thoracic material
partitions share one parent twist. Their kinetic terms are assembled together;
there is no extra torso mass or additive force that pretends to change inertia.
The native target and all source recipes are bound through
`native_composition_v1/plan.json`; the original torso removal ledger is checked
against the frozen cervical source receipt before use.

## Coordinates and material derivatives

`CoupledThoracicMetric` uses 32 internal coordinates: the previous 26 thoracic
coordinates, then neck `pitch2, roll2, yaw2, pitch1, roll1, yaw1`. Its 38 speeds
are parent body-frame translation/angular velocity followed by those internal
coordinate derivatives. Two rib speeds remain eliminated, leaving 36
independent source-system speeds. Neck coordinates are restricted to an
explicit ±0.05 rad engineering domain and their expanded source ranges.

All 18 retained cervical couplers are explicit affine two-knot splines. The
implementation constructs the exact 24-by-6 expansion `q_full=N r+offset`, checks
rank six and rejects cycles, unsupported functions and unresolved joints. For
this retained recipe `Ndot=0`; this is not assumed for arbitrary donor models.
Every joint's sequential source rotations and fixed frames carry analytic first
matrix derivatives and directional second derivatives. Product rules preserve
cross derivatives along the neck chain. The zero-state body frames reproduce
the retained source recipe; nonzero-state frames and derivatives are checked
against an independent evaluator and finite differences.

For each cervical body, translational and angular Jacobians both include the
same parent twist and six neck rates. Full rotated COM inertia, parent/neck
cross terms, relative-angular acceleration, and absolute-angular gyroscopic
terms enter its metric and inertial bias. Those terms are embedded in the same
parent block as the thoracic material metric. Tests reproduce the exact native
reference 6-by-6 parent spatial mass block, then independently check kinetic
energy, differential work, and world linear/angular momentum balance.

## Atomic state and pressure ownership

`ThoracicCompositeOwner` owns position, proper rotation, all internal
coordinates/speeds, mechanical time and geometric-pressure work. Public state
reads return independent copies. `advance(dt,pressure_pa)` computes a complete
Lie midpoint candidate and validates it before publishing any state. Failure
at the initial or midpoint evaluation leaves state and revision unchanged.

`checkpoint()` returns an opaque, instance-bound token backed by a private copy.
`restore(token)` restores all mechanical state and its work ledger; a monotonic
revision advances to avoid implying that old proposals are current. Saved state
cannot be edited through the token. Foreign/released tokens fail, and at most
eight live checkpoints are retained. Callers serialize access to an owner.
These tokens roll back no native physiology, controller, or external source.

The pressure port exposes the engineered geometric cavity `V(q)`, its speed
Jacobian, `Vdot`, outward mechanical pressure, generalized force and power:
`Q=p J_V`, `u^T Q=p Vdot`. Neck and rigid-parent components of the closed-cavity
pressure force are zero. A constant-pressure step records exact geometric work
`p*(V_end-V_start)` and reports midpoint `p*Vdot*dt` separately, including its
quadrature residual. Native source-flow signs, native gas volume, chest recoil,
ideal drive transfer and muscle chemistry remain explicitly unassigned.
Physical locked-rib reaction wrenches are still withheld.

## Installed native API and executable replacement operator

`native_api_evidence_v1` retains Apache-2.0 header bytes and notices. The installed
`SimbodyMatterSubsystem.h`, `MobilizedBody_Custom.h` and `Body.h` exactly match
the retained source headers. The custom mobilizer interface supplies transforms,
H/H-transpose, Hdot and N/Ndot mappings for an attached body. It does not expose
an arbitrary distributed mass-matrix override. Attaching a custom mobilizer to
one unchanged rigid torso therefore does not implement this material metric.

The installed matter API supplies the necessary same-state operators:

- `calcM` gives the native tree mass matrix.
- `calcResidualForceIgnoringConstraints` with zero applied forces and zero
  udot gives the inertial bias. Its documented residual is
  `M udot + inertial_bias - applied_force`.
- `calcFrameJacobian` gives the world-frame torso-origin spatial Jacobian.
- `calcBiasForFrameJacobian` gives physical world spatial acceleration at
  native udot zero; `multiplyByN` supplies native qdot from u.

Simbody spatial vectors are angular-first and world-expressed. The adapter
converts them to the linear-first body-frame twist. Its body linear bias is
`R^T a_world - omega cross v`; its angular bias is `R^T alpha_world`.

`compose_native_operator` implements the full replacement, where J maps native
speeds to the parent twist, S maps the enlarged native/internal speeds to the
38 composite speeds, beta is the parent-map acceleration bias, and subscript 0
is the removed original torso:

```
M_new = diag(M_native - J^T M_0 J, 0) + S^T M_composite S
b_new = embed(b_native - J^T(M_0 beta+b_0))
        + S^T(M_composite [beta;0]+b_composite)
Q_pressure_new = S^T Q_pressure_composite
```

The operator exposes the enlarged metric and inertial bias. It does not pass
an effective force to the unchanged native solver. Native constraints,
prescribed motion and Ndot still require explicit projection in the full
mechanical integrator. Other native forces must be audited by owner; especially
mass-dependent gravity must not keep an old torso contribution while also
adding its replacement. Tests use the actual composite material metric with a
clearly declared native-remainder fixture and a nontrivial parent map/bias;
this algebraic proof is not represented as a native runtime extraction.

## Prepared scheduled native proof

`native_thoracic_operator_snapshot.cpp` is a separate, read-only extractor for
the frozen 22-body/92-muscle model. It initializes one native state, sets a small
declared velocity vector, and emits the same-state operators above. It performs
no integration and applies no force. `build_thoracic_native_snapshot.py` prepares
one translation unit with retained input/library hashes, a 2 GiB memory cap,
120-second CPU compile cap and wall-clock timeouts. Native work requires the
parent agent's scheduling grant. `verify_thoracic_native_snapshot.py` is ready
to check the actual extracted reference mass/bias block, full active rank and
new cross-inertia once such a snapshot exists. Compilation and extraction are
not yet claimed by this source-only document.
