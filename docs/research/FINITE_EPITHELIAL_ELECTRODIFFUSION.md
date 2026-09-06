# Finite epithelial ion, charge and energy component

`ihm/assembly/epithelial_electrodiffusion.py` adds an executable finite-reservoir epithelial layer. It supplies a conservative alternative preparation to the previous independent fixed-bath RC graphs. It does not change native physiology, the regional RC patch, IBM, or anatomy. No skin regeneration, calibrated whole-body electricity, or validated intact-human trajectory is claimed.

## Acquired human observation

[Abe et al. (2019), PLOS ONE](https://doi.org/10.1371/journal.pone.0219198) supply an original XLSX with a small human supplement alongside their predominantly porcine experiments. Original article XML and XLSX are retained in `data/research/skin_epithelial`, with byte sizes, SHA256, acquisition URLs, license and preparation in `source-receipt.json`. The article grants CC BY 4.0; attribution and license link are retained.

The human observation is posterior forearm in one healthy adult female, aged over 20, during sequential tape stripping. `Sheet1!C48:E51` holds four observations:

| Tape strips | Reported TEP mV | TEWL g/m²/h |
| ---: | ---: | ---: |
| 0 | −10.122183 | 8.51 |
| 1 | −9.557791 | 8.69 |
| 4 | −9.364351 | 10.80 |
| 7 | −8.764410 | 13.18 |

The protocol mentions ten strips, but the supplied worksheet has no corresponding row; no value is imputed. These are repeated states of one participant, not four donors. Other workbook groups are porcine and are excluded from the human extraction. No human current, capacitance or electrical charging transient is supplied. Tape-strip count is not elapsed time.

Main-text porcine TEP uses subepidermis relative to surface and positive values, whereas the human supplement gives negative values. Human wiring polarity was not independently resolved. The model conditions only on observed **magnitude**, requiring an explicit polarity choice; it does not silently flip the source. Its basal-positive convention remains a prior. Cultured human keratinocyte Vm (Mauro 1990, discussed in the calibration-gap audit) is a different observable and preparation and is not used to set this participant's cell Vm.

The only fitted coefficients are the intercept and slope of a descriptive `|TEP| = a + b * tape_strips` relation. Three observed states (0, 1, 4) identify two coefficients; seven strips is held aside. The held prediction exceeds observation by 0.060408 mV. This small within-participant check provides neither population uncertainty nor channel/RC identification. The fit is never used to infer conductance, capacitance, pump strength or a healing law.

## Conservative mechanism

There are separate apical bath, epithelial cell and basal bath ion inventories for Na+, K+ and Cl−. Three directed interfaces connect apical→cell, cell→basal and apical→basal; the last represents a paracellular route. The two membrane interfaces have explicit total capacitances. Basal voltage is the reference, not an infinite ionic reservoir.

Let `B` be the compartment-interface incidence matrix, `ξ` the interface transport extents in mol, `n=n0+Bξ`, and `z=(+1,+1,−1)`. Two membrane capacitances construct the graph capacitance matrix `L`. A fixed initial background/countercharge is declared as `q_fixed=L φ0−F Σ(z n0)`. Thus `q=L φ0+F Σ[z(n−n0)]`, and solving `L φ=q` with basal reference enforces the charge-voltage relation. No evolving compensating chloride is inserted. Changing an observation-conditioned initial voltage redefines initial countercharge; it is not represented as physical transport.

Electrochemical potential is `μ=RT log(c/c_ref)+zFφ`, with ideal activities and reference concentration 1 mol/m³. A passive interface transfers `J=g(μ_source−μ_destination)/(zF)²` mol/s. Each molar transfer is subtracted from one inventory and added to another. The `transport()` interface exposes individual passive/active molar fluxes and rates. `simulate()` integrates transport extents rather than independent compartment amounts, retaining species conservation by construction.

Free energy includes ideal solution mixing and membrane storage `U_electric=½ φᵀLφ`. For passive transfers, dissipation is `Σ J Δμ ≥ 0`. Prescribed active transfer has signed external power `−Σ J_active Δμ`; its integral is explicitly recorded. The ledger is `ΔG + dissipated_energy − active_work = 0`. A negative active-work term means energy returned to the ideal actuator. This is an external-work port; it is **not an ATP reservoir, pump efficiency or native metabolic expenditure**. The demonstration's 3Na outward/2K inward flux is a declared prior, without ATP kinetics or a claim of endogenous pump activity.

Outputs expose inventories, transfer extents, charge, all node voltages, apical and basal cell Vm, TEP, electrical energy, total free-energy change, dissipated energy and external active work. Positivity and finite inputs are checked; reservoir exhaustion raises instead of silently replenishing ions.

## Uncertainty and ownership

`EpithelialPatch` has no physiological coefficient defaults: callers supply volumes, concentrations, conductances, capacitances, temperature and initial voltages with explicit preparation/parameter provenance. `condition_tep_magnitude` validates source quantity and preparation and requires explicit polarity. It only initializes the measured magnitude; cell Vm and transfer parameters remain priors.

The demo uses four capacitance/conductance sensitivity scenarios (0.5× and 2× combinations), retaining every parameter array. Its min/max envelope is explicitly a prior-scenario envelope, not posterior probability or a calibrated confidence interval. The only human-conditioned electrical value is initial |TEP|. Fixed volumes, nominal temperature, bath ions, cell Vm, active transfer and all kinetics remain stated illustrative assumptions.

These inventories belong to an **independent epithelial preparation**, including its finite baths. They are not extra native Skin or blood volumes. There is no native engine object, concentration overwrite, native mass mutation or duplicate native ATP consumption. Integration into native physiology would need a single-owner exchange port with the exact opposite signed molar and work transfers applied to native stores. GI reuse requires different species/stoichiometric cycle definitions and GI-specific evidence; neutral glucose, amino acids and ATP are not represented by this three-ion primitive. Skin coefficients must not transfer to GI.

The model is one effective living epithelial layer, not a stratified epidermis or electrodiffusive spatial PDE. It omits water/osmotic mechanics, activity coefficients, channel gating, membrane area changes, differentiation, perfusion, nerves, migration and regeneration.

## Reproduction and verification

Run:

```sh
.venv/bin/python scripts/verify_epithelial_electrodiffusion.py
.venv/bin/python scripts/build_epithelial_electrodiffusion.py
```

The builder parses retained original bytes without network access and writes the extracted receipt and small source-conditioned scenario artifact. Tests verify species conservation, capacitor charge, passive free-energy decline, active-work closure, equilibrium, distinct coupled TEP/Vm, invalid-input rejection, source-preparation separation, observation-only identification and explicit prior uncertainty. No native engine or heavy job is involved. The four generated scenarios close relative energy to at most approximately 1.1×10−8, with species residuals below 4×10−33 mol and charge residuals below 7×10−28 C; these are numerical conservation checks, not physiological validation.

The next experiment remains paired voltage/current and impedance on the same area-defined human preparation, followed by ionic substitutions or well-controlled transport perturbations. Only such observations can identify RC/current parameters beyond the source-conditioned initial voltage.
