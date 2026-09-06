# Native Skin inventory to finite epithelial patch: staged contract

`ihm/assembly/skin_epithelial_inventory.py` supplies a conservative **alternative snapshot**, not another live native compartment. `bind_inventory` subtracts each represented species from its original owner into a shadow partition and retains the exact complement. `run_passive` reassembles terminal proposed owner masses. Neither function writes native state, reserves a running owner, advances native time, or adds its masses to an existing body total. Multiple calls are alternatives; they cannot be summed as concurrent patches.

## Source support and materialization

`exact_patch_support` uses the existing byte/hash-checked exterior component contract and an explicit ordered subset of its source triangle IDs. It computes subset area from canonical metre vertices. The retained full inferred component is 1.7804602548390722 m² (`data/research/engineered_skin_territories/materialization.json`); exclusivity, self intersections and closure remain unvalidated. A support identity binds geometry, component evidence and triangle IDs. Production callers must obtain area/identity together from this helper; the low-level partition accepts explicitly labelled synthetic supports for tests.

For area A, epithelial depth h, water fraction w and cellular fluid fraction f:

* apical film volume = A × supplied film depth;
* cell fluid volume = A × h × w × f;
* basal fluid volume = A × h × w × (1−f).

These are engineering allocation priors, not measured epidermal morphology or a mapping of every native Skin cell to a keratinocyte. Each supplied species mass is partitioned by represented fluid volume / owner fluid volume. Na/K/Cl moles use explicit supplied molar masses in g/mol. Unmodelled supplied species retain their masses and do not participate in transport. Unreported species remain unresolved. A supplied unknown/NaN mass fails instead of becoming zero. There is no numerical calcium or fixed-anion inventory invented to neutralize the three ions.

`native_skin_owner` converts existing owning observations from mL to m³. Skin IC supplies cells, Skin EC or a declared regional EC owning leaf supplies basal fluid. The regional aggregate is rejected as nonowning. The apical film must be a distinct explicitly supplied external fluid owner; native Skin EC is not silently cloned into that position. The current native bulk observer supplies no apical film inventory. Geometry-to-region membership remains an explicit missing source mapping, not the engineering regional installation fractions.

## Native diffusion already owns the ions and a voltage

Source: `data/raw/physiology/biogears/projects/biogears/libBiogears/src/engine/Systems/Diffusion.cpp`.

* `Initialize` constructs vascular/EC/IC diffusion sets for tissues including Skin (lines 37–87).
* `CalculateNonLinearDiffusionMethods` invokes active ion diffusion on every set (lines 272–285).
* `CalculateActiveIonDiffusion` reads native tissue membrane potential (line 461), computes Na/K/Cl Nernst and pump terms, sets a new native tissue voltage (line 500), then debits/credits actual EC/IC masses with donor caps and mass balancing (lines 550–585).
* `SodiumPotassiumPump` specifies an empirical rate (lines 727–746). This routine does not itself debit an explicit ATP reservoir. Its existing rate must not be reused as a free energy source in the new operator.

The existing voltage accessor is `SETissueCompartment::GetMembranePotential(ElectricPotentialUnit::V)` on `SkinTissue`, guarded by `HasMembranePotential()`. `scripts/native_tissue_ports.h` currently omits that readout. An observational field `tissue.Skin.bulk_membrane_potential_v` could expose it without another state owner. It is native bulk IC-minus-EC membrane potential, not a human epidermal transepithelial observation.

For the new three-compartment operator, basal φ=0; TEP=φ_basal−φ_apical; basal membrane voltage=φ_cell−φ_basal; apical membrane voltage=φ_cell−φ_apical. Native bulk membrane potential is not automatically assigned to either epidermal membrane. Both initial patch potentials require explicit provenance.

## Charge, work and unresolved native interfaces

The initial fixed-background charge is reported explicitly as `L φ_initial − F Σ(z n_initial)` per fluid compartment. It combines unresolved other charge with the electrical initial-condition prior. It is not evidence for hidden chemical compensation. Its owner and the owner of initial capacitor field energy must be named. Capacitor charge and energy thereafter use the existing two-capacitor operator.

This binding exposes **passive transport only**. Instantaneous electrochemical forces drive the existing integrator; terminal energy decay alone is not used as a kinetic admission test. Dissipated energy is reported as a proposed heat credit to a named thermal owner; no native heat has actually been added. There is no ATP or other active work supply. Active integration requires a finite donor/work port and cannot be enabled by naming a pump.

Before live native transactions are possible, implement and verify:

1. A persistent disjoint owner partition with source-face overlap rejection and per-owner volume/mass reservations; aggregate views must never be debited.
2. Exclusion of represented IC/EC fractions from existing native active-ion diffusion and other duplicated transport. Regional leaves need actual diffusion aggregation semantics preserved for the residual native owner.
3. Atomic epoch/sequence-checked mass transfers, native substance balancing and a defined timestep order; these proposals cannot be applied twice or to newer snapshots.
4. Persistent membrane capacitor charge, fixed-background charge and electrical energy owners; a justified mapping from patch membrane states to any native bulk voltage readout.
5. Finite active-work supply and thermal-energy credits, with ATP stoichiometry only where supported. Initial field energy must be transferred or explicitly initialized, never silently manufactured on each rebind.
6. Full species inventory coverage and explicit source-area-to-region membership. Current seven-species readouts are partial observations.

Validation: `scripts/verify_skin_epithelial_inventory.py` checks exact retained source-face area, rejected invalid faces, native volume units, all-supplied-species partition identity, caller snapshot preservation, rejected aggregates/missing film/over-allocation/nonfinite priors, passive species conservation, background-charge relation, distinct voltages and preserved Albumin. These are contract tests, not whole-body RC calibration or regeneration evidence. No native source edits/builds/runs are part of this increment.
