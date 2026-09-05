# Native circuit and compartment atlas

Run `.venv/bin/python scripts/verify_native_circuits.py` then `.venv/bin/python scripts/build_native_circuit_atlas.py`. The default source is `data/derived/physiology/native_baseline_v2/states/native_stabilized.xml`; `--state PATH` selects another actual native saved state and `--output DIRECTORY` selects the artifact directory. `ihm.native.circuits.native_circuit_graph(state_path)` exposes the importer directly.

The compact `data/derived/native-circuits/index.json` links to `graph.json`, approximately 6.7 MB. Both preserve the native state's absolute path, SHA256, byte length and simulation time. They label parameters **native-engine initialized or stabilized, not independently calibrated**. This is extraction of actual native coefficients and physiological state, not an illustrative substitute or reimplementation of BioGears equations.

## Identity and units

Node and path IDs combine the physical family and exact native name, for example `fluid:Aorta1` and `thermal:Ambient`. A fluid Ambient and thermal Ambient remain distinct. Named circuit views reference these unique objects; 138 nodes and 247 paths are shared among views. Their counts, volumes and flows must not be multiplied by circuit membership. Conflicting duplicate object definitions, dangling endpoints or invalid circuit references fail import. Alternative equipment circuits remain present without claiming they are active.

Each scalar retains raw text, parsed value, native unit, SI value/unit, status and native read-only flag. Current, Next and Baseline properties stay separate. Explicit conversions cover actual native fluid/thermal units and supported electrical units; dimensions are checked against circuit quantity type. For example mmHg converts with 133.322387415 Pa/mmHg; mL/mmHg converts to m³/Pa; thermal kcal/degC capacitance converts to J/K using 4,184 J/kcal. Absolute node Celsius temperatures receive the Kelvin offset; temperature differences do not. Unknown circuit units, incompatible dimensions and NaN or unsupported infinite circuit values are rejected.

The source deliberately has infinite reservoir volumes at Ambient, GasSource and OxygenTank. Their nine current/next/baseline properties retain raw `INF` and status `infinite_boundary`, with JSON null numeric values. They are not replaced by fabricated finite volumes. These ideal boundary reservoirs are not finite body compartments. All serialized finite circuit properties have finite SI values. Some noncircuit biochemical/system units are preserved as `unconverted_native_unit` with raw values and no invented conversion; such values must not enter SI computations until supported.

## Native inventory

The stabilized state contains:

| Native category | Unique count |
|---|---:|
| Fluid nodes / paths | 202 / 348 |
| Thermal nodes / paths | 20 / 47 |
| Circuit views | 13 |
| Liquid / gas / thermal / tissue compartments | 117 / 33 / 20 / 13 |
| Liquid / gas / thermal compartment links | 112 / 37 / 47 |
| Circuit scalar properties | 2,509 |
| Exported System elements / scalar properties | 15 / 292 |

The path mechanism inventory includes 154 resistance-bearing paths, 73 compliances, 57 pressure sources, 21 flow sources, 15 heat sources, three temperature sources and seven thermal capacitances. A path may carry more than one kind of property, so these are not independent physical component counts. Gates, valves and switches are retained separately as native states. All 196 compartment links resolve to actual native circuit paths and all native compartment-node mappings resolve.

System elements retain their literal native type: blood chemistry, cardiovascular, drugs, endocrine, energy, gastrointestinal, hepatic, nervous, renal, respiratory, tissue, environment, anesthesia machine, ECG and inhaler. Configuration, patient scalar properties and active-substance properties are also retained. This inventories exported state and parameters; it does not claim that every internal equation or biochemical mechanism is serialized in the XML. The ECG state does not imply an electrical circuit: **this actual state contains no electrical circuit nodes or paths**. The typed importer supports that family when a future source actually exports it.

## Compartment crosswalks

Compartment IDs preserve native medium (`liquid:Aorta`, `gas:Ambient`, `tissue:BoneTissue`). Records contain their exact native anatomical label, mapped circuit-node IDs, child compartment IDs and all exported scalar properties. Repeated substance quantities carry their substance name explicitly. Links preserve medium-specific source and target compartments plus circuit path. Native compartment graphs retain their original memberships.

Examples include Aorta mapping to Aorta1/Aorta2/Aorta3, renal and nephron hierarchies, extracellular and intracellular tissue compartments, and native vascular parent/child groups. These are source semantic mappings, not atlas coordinate registrations. Parent compartments aggregate children; their volumes must not be added to child volumes. A native compartment name is the anatomical label evidence; no guessed coordinate transform or subject correspondence is introduced.

## Balance diagnostics

Every circuit view has its own incidence audit, using current signed flow (m³/s), heat transfer rate (W) or current (A). Positive net flow means inflow. Reports identify reference nodes, missing fluxes, storage-connected nodes, all-branch net inflow, transport-only net inflow and storage-branch net inflow. Unknown flux stays incomplete, not zero. The maximum absolute reported balance excludes reference nodes, which are allowed to exchange with boundaries.

Compliance/capacitance branch flux belongs in the native circuit's algebraic balance. Removing it can create transport residual because volume or energy is being stored. The importer therefore never labels transport-only residual as a leak. Nor does a single snapshot establish a derivative of volume, time-integrated mass conservation or physiological conservation across every biochemical process. Current and Next state often coincide after native advancement; dividing their difference by an invented timestep would be wrong.

At the stabilized snapshot, FullCardiovascular's maximum complete nonreference net inflow is approximately 5.30e-16 m³/s; Respiratory is 2.11e-15 m³/s. The combined Temperature view's residual is 8.17e-7 W. Partial Cardiovascular and Cerebral views have 4.16e-5 and 1.96e-5 m³/s residuals, and internal/external thermal views each have about 83.56 W residual because the complete coupled view has different boundaries. Alternative equipment views can be unsolved or incompatible with the active state. These values are retained as diagnostics without a universal pass threshold or a false conclusion of native leakage.

## JSON contract and tests

The full graph contains `nodes`, `paths`, `circuits`, `compartments`, `compartment_links`, `compartment_graphs`, `systems`, `configuration_properties`, `patient_properties`, `active_substances`, `summary`, `source`, `simulation_time`, `parameter_status` and `limitations`. Nodes and paths expose their `properties` dictionary keyed by native scalar names; circuits and crosswalks use stable IDs. The compact index keeps source, summary, names, circuit-level balance summaries and `graph_path`.

The test fixture verifies mixed pressure/flow units, independent fluid/thermal names, shared circuit deduplication, compliance storage balancing, thermal capacitance, unknown-unit rejection, dimension mismatch rejection, NaN rejection and infinite-reservoir serialization. It also imports the real stabilized state, verifies its 222/395/13 counts, all path endpoints, 15 System elements and strict finite JSON serialization. The numerical audit tests known incidence cancellation rather than assuming every partial native circuit should balance independently.
