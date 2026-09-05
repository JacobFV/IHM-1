# Integumentary voltage interactions

IHM-1 distinguishes three electrical quantities:

| Quantity | Definition | Support / observation |
|---|---|---|
| Membrane voltage | intracellular minus local extracellular potential | a particular cell membrane; intracellular electrode or calibrated optical measurement |
| Transepithelial potential | basal minus apical extracellular potential | epithelial barrier; potential difference across tissue |
| Extracellular field | negative spatial gradient of extracellular potential | tissue/fluid interfaces; spatially referenced electrode measurements |

The supplied model grounds the basal bath. Its apical potential is therefore the
negative of transepithelial potential. It does not identify either quantity with
ECG voltage, neural activity, electrodermal conductance or a latent representation.

## What the literature constrains

- [Nuccitelli et al., imaging mouse and human skin wound fields](https://pmc.ncbi.nlm.nih.gov/articles/PMC3086402/): direct wound-field observations motivate a spatial extracellular model. Species and measurement geometry must accompany extracted values.
- [Human keratinocyte ENaC/galvanotaxis experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC3666251/): channel perturbation and migration experiments motivate cell-specific transduction mechanisms. They do not justify one electrotaxis coefficient for all skin cells.
- [Age and human wound electric fields](https://pmc.ncbi.nlm.nih.gov/articles/PMC3228273/): barrier recovery and visual closure differ, and age/context matter for comparison.
- [Levin, molecular bioelectricity](https://pmc.ncbi.nlm.nih.gov/articles/PMC4244194/): the conceptual motivation for non-neural membrane states, ion channels/pumps and gap-junction networks. This is a mechanistic synthesis across biological systems, not human skin calibration data.

No numeric parameter has been fitted from these articles in this repository.
Source cards are catalog-only. Regenerative and anatomical pattern-control claims
must retain species, tissue, developmental stage and intervention context.

## Executed mechanisms

For each cell, membrane current is the sum of individual Ohmic ionic currents,
an imposed pump current and currents through declared gap junctions:

```
C_i dV_i/dt = -sum_k g_ik (V_i - E_ik) - I_pump_i
             + sum_j g_ij (V_j - V_i)
E_ion = R T / (z F) log(c_out / c_in)
```

The runtime linearizes Nernst potentials around declared ion concentrations.
Intracellular concentrations are uncertain fixed reservoirs. It does not simulate
pump stoichiometry, ion depletion or nonlinear channel gating. Neighbor currents
are paired so total electrical charge transfer across each junction cancels;
voltage rates differ when cell capacitances differ.

For an apical surface site:

```
C_i dphi_i/dt = I_i - G_barrier_i phi_i + sum_j G_ij (phi_j - phi_i)
E_edge = -(phi_b - phi_a) / distance(a, b)
```

A wound removes active cells at the example's central site and makes the surface
barrier a shunt. The cell-contact graph has no junction across that gap, while the
extracellular sheet remains conductive. The default line is a small synthetic
cross-section, not a scanned skin mesh. Custom 3D positions and graph edges can
represent separately sampled cellular and extracellular supports.

## Where the next mechanisms belong

Cell migration, proliferation, differentiation, calcium signaling, ion/pump
regulation, matrix remodeling and barrier repair belong in processes over physical
state. They are not implemented by renaming electric field magnitude as a healing
score. Fibroblast and melanocyte phenotypes are supported as individually
parameterized cells, but no phenotype-specific coefficients are inferred for them.
A proper bidomain/electrodiffusion extension must let extracellular potentials
react to membrane currents and enforce ion mass balance; the present prescribed
bath model cannot make that claim.

Skin–blood–lymph coupling needs geometric interfaces and matched units: dermal
capillary filtration, interstitial protein/ion transport, lymph uptake and tissue
oxygenation. The registry admits these separate supports, but the current skin
and fluid examples have not been coupled into one calibrated wound model.
