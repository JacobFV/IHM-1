# Respiratory mechanical ownership and next interface

The unified runtime currently has a real **external-load → BioGears respiratory circuit → airflow/gas** path, followed by a native-volume-driven display deformation. It does not implement **IBM brain → phrenic activation → mechanical diaphragm/thorax → airflow**.

## Existing paths

| Component | Actual ownership and behavior |
|---|---|
| `scripts/native_body_ports.h` | Observes native lung volume, requested driver amplitude/frequency, cycle fraction, applied circuit pressure, airway flow and arterial gases. These observations do not actuate anything. |
| Held BioGears `Nervous.cpp` / `Respiratory.cpp` | Native chemoreceptor feedback updates breathing amplitude/frequency; `RespiratoryDriver()` generates the circuit pressure waveform and writes `EnvironmentToRespiratoryMuscle`. This is a lumped physiological controller, not simulated phrenic/diaphragm recruitment. |
| `ihm/assembly/brain.py`, `sensorimotor.py` | IBM neural state receives MAP, saturation and temperature in this runtime. Its autonomic commands explicitly have `applied_to_body=False`. Somatic decoding reaches the existing muscle catalog; it has no respiratory effector port. Arterial CO2 is observed from BioGears but is not passed into this IBM respiratory feedback loop. |
| `ArticulatedBodyPlant` | The retained model has 22 bodies and 92 muscles. Torso includes the chest; there are no separate rib/diaphragm joints, and no diaphragm/intercostal actuator in the catalog. Native dynamics and contacts therefore cannot solve breathing deformation. |
| `EmbodiedRespiration` | Computes `q = b*(Vnative−Vreference)`, with `b=K^-1*A/(Aᵀ*K^-1*A)`. It changes rib/sternum centroid translations, diaphragm shape, lung determinants and a smooth skin display field. No respiratory velocity, inertia, active fiber state, mechanical recoil energy or support reaction is integrated. |
| `BodyRespiration` in `respiration.py` | A separate three-mode mass/damping/stiffness integrator exists, but the unified factory does not instantiate it. Its held effective masses `[1.5,1,0.5] kg`, stiffness and damping are generic priors; activating it unchanged would add an unpartitioned 3 kg and duplicate chest recoil already owned by BioGears. |

`EmbodiedRuntime.step()` first advances articulation, then projects supplied user point forces through the respiratory shape Jacobian. Its pressure is `p_ext = −Σ F·∂x/∂V`. The native hook adds this to the generated driver **after PreProcess and before Process**, preserving native gas storage and passive recoil. The hook records `W_ext = −p_ext*(V_after−V_before)` for the held pressure interval.

That work is a circuit-boundary receipt. It is not a closed, finite-body diaphragm/torso work balance. The rigid plant sees point-force work on its own motion, while the respiratory Jacobian describes additional display displacement. These terms can form a valid decomposition only when the same physical material map/contact trajectory and external owner are shared and their summed work is audited. Current native bed quadrature and garment skin surfaces follow rigid segment transforms, not the respiratory display field. Their resultants are retained as `unresolved_contact_wrenches`; they are not fed into respiratory strain work. Thus a mattress does not yet resist displayed chest expansion, and respiratory pressure does not generate a returned finite-body diaphragm/support wrench.

## Concrete next implementation interface

Create an explicit, source-bound `ThoracicMechanicalPort` with the following contract; do not silently promote the display module:

- **Identity/ownership:** frozen anatomy and material hashes, explicit source face/node bindings, frame epoch, one named neural drive owner, native cavity/gas owner and exclusive recoil/inertia ownership. Preserve the common canonical/native rigid world. Partition any added thoracic mechanical mass/inertia from the existing torso, including COM and cross-coupling effects; do not import the standalone 3 kg prior as extra mass.
- **Mechanical state:** actual generalized chest/diaphragm coordinates and velocities, reference cavity volume, nonlinear `V(q)` and `J_V(q)`, source material-point positions/velocities, actuator state, elastic energy and dissipative history. The three existing modes can be a declared reduced starting geometry, but rib rotations, abdominal support and muscle insertions remain uncalibrated. Report that reduction explicitly.
- **Boundary inputs:** native generated respiratory demand, pleural/gas pressure, registered user/bed/cloth tractions and current torso/support pose. Use the same breathing-deformed source surface in contact and visualization. A resultant torso wrench alone cannot reconstruct respiratory virtual work.
- **Boundary outputs:** cavity volume/rate, gas-facing mechanical pressure, generalized active/contact forces, equal/opposite torso/abdominal support wrenches, source-surface motion and separate active work, passive stored energy, dissipation and pressure-volume work receipts. For pressure force use `Q_p=J_Vᵀ*p`; verify `Q_p·qdot=p*Vdot` with a declared pressure sign. This is an interface constraint, not permission to prescribe both volume and pressure independently.

The final causal implementation needs an **opt-in native respiratory boundary variant**. Retain native airways, gas transport and lung recoil. Identify the existing left/right pleural-to-respiratory-muscle chest-wall elements and transfer their mechanical role explicitly before adding a physical thorax stiffness; never retain overlapping passive recoil silently. Replace the native ideal respiratory pressure drive with the chosen mechanical actuator boundary when that actuator is enabled; otherwise the same inspiration would be powered twice. Initially the native controller may remain the clearly labeled demand owner. A later IBM/phrenic decoder requires explicit source-supported recruitment and afferent mapping; relabeling native pressure as neural spikes is not that implementation.

## Unified runtime ordering and acceptance

The current two-process runtime has no verified whole-physiology rollback. Therefore a fixed-point loop must not repeatedly call `signed_step()` while restoring mechanics alone. Stage a coupled respiratory/circuit candidate inside a bounded native transaction, or expose a specifically verified reversible circuit candidate API. Execute physiological PreProcess, meal consumption and the signed muscle energy boundary once per accepted interval. Commit the compatible thorax/cavity state once, then publish mechanics, physiology and the shared surface at the same endpoint. Any failure after irreversible physiology work remains terminal until a genuinely exact joint checkpoint exists.

Before live enablement, require zero-boundary native parity; common-clock and impulse/power identities; cavity-volume consistency; unchanged total body mass; isolated chest-contact perturbations affecting airflow through resolved tractions; and actuator work/chemical/heat accounting that does not double the existing basal or muscle budget. A display-only test cannot establish these properties.

## Review validation

Read-only source trace plus bounded retained model XML inspection confirmed 22 bodies/92 muscles and the missing respiratory joints/actuators. Six `verify_embodied_respiration.py` and three `verify_respiratory_feedback.py` tests pass. Those tests establish display/material Jacobian and external pressure conjugacy only. No native process, build, full geometry load, model change or supported-breathing claim was made.
