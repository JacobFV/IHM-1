# Fixed-activation static potential: derivation and bounded solver proposal

No native run, integration, controller change or new static solve was performed for this design. The target remains the separately bound 98-muscle, seam-omitted reference. All physical equilibrium thresholds remain unchanged.

## What is conservative at zero speed

Gravity has potential −Σ m g·x_COM. Each zero-speed expression force Q_i(q_i,0) has a scalar primitive −∫Q_i dq_i. The 20 retained expressions are separable; their velocity terms vanish in this diagnostic. The native expression component does not provide a corresponding potential override, so simply querying native total potential energy would omit these primitives.

The existing skin/mattress energy is already source-defined. The series law solves equal pressure between bounded skin and mattress compression, and stores the sum of their energies. On an interior valid branch, differentiating with respect to total approach gives pressure; the two internal deformation derivatives cancel at equal pressure. The unilateral zero-contact branch is continuous, and the same source deformation limits must reject invalid trials. No energy is assigned to the omitted seam cell. Friction/damping at zero velocity contributes zero force here; this does not make a moving contact experiment conservative.

## Muscle derivation and its limits

For fixed activation a, zero fiber velocity and a smooth path length L(q), define fiber force

F_f(l_f,a) = F_iso [a f_AL(l_f/l_opt) + f_PE(l_f/l_opt)].

Let the constant pennation height be h and the fiber projection s(l_f)=sqrt(l_f²−h²), so tendon length l_t=L−s(l_f), and cos(alpha)=s/l_f. A candidate *effective* potential is

U_eff(L,l_f;a) = U_t(L−s(l_f)) + ∫ F_f(l_f,a) dl_f.

Its fiber derivative is F_f−T/cos(alpha). Thus an interior stationary fiber satisfies the native static force balance T=F_f cos(alpha). Along a differentiable stationary branch l_f*(L), the envelope derivative is dU_eff(L,l_f*(L);a)/dL=T. With source-consistent moment arm r_i=−∂L/∂q_i, the force is Q_i=−∂U_eff/∂q_i=T r_i.

This derivation requires a stable, locally single-valued internal branch, finite pennation geometry and no unaccounted minimum-fiber constraint. Positive internal second derivative is necessary for eliminating the fiber as a local energy minimum. Descending active force-length behavior can introduce nonconvexity/multiple branches; blindly taking a global fiber-energy minimum could switch away from the branch selected by native initialization. Native clamped/minimum-length warning branches need separate inequality/reaction treatment and must not silently pass the interior derivation.

The active primitive is a mathematical fixed-activation work potential, **not physiological stored energy**, chemical energy or metabolic credit. If activation changes, its parameter derivative supplies an additional term; in moving muscles, force-velocity and damping invalidate this static potential reduction. It must never be substituted into the runtime signed muscle energy ledger or used to reset M0.

Retained Thelen source initializes a bounded force-equilibrium solve and explicitly distinguishes minimum-fiber warnings. Its active curve is Gaussian, passive curve exponential above normalized length one, and tendon curve exponential/linear with zero tension below slack. The diagnostic analytically integrates those exact functions. Millard source likewise uses active + passive + damping fiber force and explicitly sets velocities to zero for static equilibrium. The same local derivation applies structurally, but the 80 actual Millard curves, lower-bound behavior and pennation have not yet passed a native full-body effective-energy gradient check.

All 80 fitted paths in the actual assembled model omit independent moment-arm functions. `FunctionBasedPath.cpp` consequently computes moment arms as the negative derivative of its length polynomial. The six lumbar geometry paths have retained native virtual-work validation. The 12 arm paths still require branch-aware wrapping/virtual-work coverage for this proposed whole-body merit function.

## Actual source-law evidence

`scripts/diagnose_fixed_activation_potential.py` reproduces the six zero-pennation lumbar Thelen branches at all 126 retained observations in `lumbar-muscle-native-lb45uirs/observations.csv`. It preserves activation 0.05 and exact model parameters. Maximum tendon-force discrepancy is 8.56965e−6 N, within each muscle's native initialization tolerance 1e−8 F_iso. The maximum effective-energy length-derivative error is 9.25511e−9 N with a central 1e−7 m difference. The minimum sampled internal stiffness is 72,861.66 N/m. These results validate these sampled branches, not every possible length or the other 92 muscles.

`verify_fixed_activation_potential.py` checks those per-muscle tolerances, work derivatives, positive sampled internal stiffness and all 20 zero-speed expression-force primitives. No runtime muscle evaluation was added. The retained native lumbar virtual-work test independently bounds moment-arm error around 1.07e−10 m.

## Concrete next solver and its prerequisite gate

The proposed solver is reduced-coordinate, source-bounded trust-region minimization of the complete *effective* potential, with a physical force gradient and a safeguarded BFGS/Newton model. Use a sensitivity-scaled coordinate metric, report its scales, and choose radius from actual versus predicted energy decrease. The effective energy replaces the poorly conditioned acceleration-squared merit; it does not replace the final acceleration test. The prototype should first remain a separate diagnostic, not a default adapter behavior.

Before any optimization, a bounded native fixture must explicitly emit/check gravity energy, skin/bed energy, the 20 joint primitives and all 98 fixed-activation muscle primitives, along with tendon force, path length, activation, fiber velocity, pennation and branch/clamp status. Check finite-difference energy gradients against actual generalized virtual work at the new reference and a few source-bounded perturbations. Small closed-loop work checks and repeated-length/branch checks must fail on inconsistent moment arms or history-dependent branch selection. Passive energy must be counted once; the active primitive and expression contributions cannot be assumed present in `calcPotentialEnergy()`. A scalar integral assembled from noisy finite differences of the whole-body force is not an adequate substitute.

True kinematic couplers are eliminated or enforced exactly. The coordinate-to-mobility work map must be explicit: pull back generalized force through the independent-coordinate velocity map, including both knee and dependent beta components, rather than dropping dependent mobility residuals. Verify the resulting gradient by virtual work. Held global translation/heading gauges remain symmetry choices whose actual residuals must still pass.

Force/pitch/roll balance should emerge from stationarity of the free root coordinates under gravity/contact energy. Adding those force balances as extra optimization equalities would produce mathematical multipliers that are not physical contact reactions; a constrained-energy stationary point with nonzero such multipliers is not equilibrium. Similarly, XML coordinate bounds are retained as trial-domain limits, but a bound KKT multiplier must not be counted as an anatomical stop torque. At any candidate, all actual accelerations, force/moment balances, gauges and native constraint errors must meet the original thresholds. A minimum resting against an unsupported coordinate bound fails that acceptance.

If the complete effective-gradient prerequisite fails or branch handling remains unresolved, retain a sensitivity-scaled support-tangent trust-region force-root method instead. Do not guess inertia weights from the emitted M·udot vector, which does not determine the full mass matrix. No physics coefficient, activation, ROM bound or equilibrium tolerance is adjusted by either numerical option.

Audit receipt SHA-256: `69e41d790e9e7e81fb80d8071cfa43bd84b5564aeb8b2e11039a4a97e9b4fdf3` (`data/derived/lumbar-supine-static-1g1q08u2/fixed_activation_potential_audit.json`).

## Prepared isolated native prerequisite probe

The source-only build manifest `data/derived/effective-potential-build-pqbqcago/manifest.json` copies the frozen engine source/headers and actual 98-muscle, seam-omitted input directory. The only engine additions are a new diagnostic header and command. It links the independently attested archived libraries; SDK headers and every copied/source input are hashed. It never changes `latest`, the old archive or physiological energy ledgers. Compilation is queued, with a 4GiB address-space limit, one thread/nice10 and a 60s process-group kill/reap cap. Python source/pullback fixtures pass; C++ compilation has not yet been attempted.

The new header evaluates separate gravity, joint, skin, mattress, native passive-muscle and active-effective primitives after the same copied-state muscle equilibration. Gravity/contact generalized forces are independently obtained from body spatial wrenches via the native system Jacobian transpose; passive expression forces are mapped into mobilities directly. It emits every muscle's native tendon force, length, moment arms, activation, fiber velocity, minimum-length margin, pennation and internal stiffness. The active Millard curve is integrated with bounded adaptive Simpson quadrature; the Thelen Gaussian uses its analytic primitive. The merit is explicitly labeled numerical only.

The native fixture is capped at 65 evaluations/45s. It plans one base, two source-bounded perturbations for each of 31 independent coordinates and one repeat (64 actual requests). It reconstructs a second-order directional gradient using the *actual assembled* coordinate-difference matrix, rejecting ill-conditioning. It independently checks muscle length/moment-arm derivatives, muscle energy/work, gravity, joint and combined skin/mattress work. It verifies native N=I instead of assuming coordinate rates equal mobilities, and explicitly pulls back both knee and dependent beta residuals. Ideal coupler reaction work must cancel; actuator forces outside the documented merit must be zero.

The preregistered diagnostic gradient gate is 1e−4 in each native generalized-force unit, with length/moment-arm error <=1e−6m, positive interior branch stiffness, fiber minimum margin >1e−10m and fiber velocity <=1e−5m/s across sampled states. These are numerical prerequisite tolerances, not relaxed physical equilibrium criteria. Actual errors and branch failures are retained regardless of gate outcome. A passed diagnostic alone does not activate an optimizer, certify global branch uniqueness, or promote an initial physiological reference. Further trust-region implementation remains conditional on actual results.
