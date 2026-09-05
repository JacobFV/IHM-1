# Regional non-neural skin electricity

`ihm.assembly.skin_bioelectric.build_skin_electric(root)` attaches the existing
heterogeneous `SkinPatch.line` electrical laws to the canonical forearm material
point. Its pinned source mesh, barycentric coordinates, tangent axes and source
hashes accompany the artifact. Nine surface sites span 1.6 mm and nine cells sit
0.1 mm below that tangent line. This geometry is an illustrative prior.

The membrane law is `C dVm/dt = sum(g_ion (E_ion - Vm)) - I_pump + I_gap`.
Nernst valences are Na +1, K +1 and Cl -1; source concentrations remain fixed
reservoirs. The apical law is `C dphi/dt = I_pump - g_barrier phi + I_surface +
I_electrode`. Basal extracellular voltage is zero. Apical voltage is not cell
membrane voltage; TEP is the negative of apical voltage. These are separate
graphs under a prescribed bath approximation, without neural coupling.

The exact affine matrix exponential records deterministic mean trajectories.
Independent 12-point Gauss quadrature evaluates ionic, pump, barrier, electrode
and internal edge currents at each interval. Node capacitor charge change is
compared to their integrated sum. Internal Ohmic pairs cancel in the global
ledger. This audits charge exchange with fixed reservoirs and boundary supplies;
it does not assert finite ion mass conservation or solve electrodiffusion.

Four experiments share fixed indexing: equilibrium baseline, central electrical
barrier shunt (100-fold conductance and zero local pump), balanced end electrode
currents (+/-30 nA), and a +10 mV initial membrane perturbation. The shunt retains
cells and is an electrical wound surrogate, not anatomical ablation. All values
are existing illustrative source coefficients or explicitly declared priors.
Human wound-field data and its fitted phenotype support only the lateral field
observation operator, never microscopic conductances, pumps or Vm.

Run `PYTHONPATH=. .venv/bin/python scripts/verify_body_skin_electric.py` for source
generator parity, equilibrium, opposing shunt gradients, electrode effect,
membrane perturbation, per-node charge closure and halved-time-step agreement.
Run `PYTHONPATH=. .venv/bin/python scripts/build_body_skin_electric.py` to write
`data/derived/canonical/skin-electric.json`. The artifact contains `anchor`,
`geometry`, `clock`, `experiments`, `runtime_sources` and measurement source
hashes. Each experiment has time-linked `frames` containing
`membrane_voltage_V`, `apical_voltage_V` and signed `edge_field_V_m` arrays.

No whole-body skin activation, uncertainty calibration, cell migration,
regeneration, healing law, neural conversion or native blood feedback is claimed.
