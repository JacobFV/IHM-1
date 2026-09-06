# Skin bioelectric calibration: evidence and next experiment

Audit date: 2026-09-05. The retained regional artifact is a charge-audited illustrative RC network, not a calibrated human epidermis. The smallest useful next step is an area-defined epithelial equivalent circuit fitted to paired voltage/current/impedance observations. Human wound-field regression cannot identify the existing microscopic coefficients. No regeneration, healing outcome, or Levin-style anatomical control is demonstrated.

## Four different observables

| Quantity | Definition and present meaning | What it cannot establish |
| --- | --- | --- |
| Bulk tissue Nernst potential | `RT/(zF) log(c_out/c_in)` for native BioGears Skin pools, treating concentration ratios as activities | Cell-specific concentrations, open-channel conductance, membrane voltage, epithelial polarity |
| Keratinocyte membrane voltage | Intracellular minus local extracellular potential; code fixes basal extracellular bath to zero | TEP or a wound field; one reversal potential alone does not determine resting Vm |
| Transepithelial potential (TEP) | Basal minus apical bath potential; code `TEP = -apical` | A single cell membrane voltage or identified Na/K pump rate |
| Lateral wound field | `E_ab=(apical[a]-apical[b])/distance`, V/m | Wound current without conductivity/geometry; cell migration or repair rate |

The established human field observation is above epidermis and below stratum corneum. Its measurement plane and sign convention must accompany every comparison. mV/mm equals V/m numerically; mV does not. Native bulk intracellular Skin aggregates are not a keratinocyte layer. The Nernst calculation also assumes ideal concentration activities and equilibrium for each individual ion, not simultaneous equilibrium of all ions.

## Existing coefficients and support

Values below are in `ihm/materialize/skin.py`; interventions are in `ihm/assembly/skin_bioelectric.py`. Every numerical coefficient remains an illustrative prior. Physical RC/Nernst/Ohmic form supplies dimensional consistency, not parameter evidence.

| Parameter | Existing value | Meaning and missing calibration |
| --- | --- | --- |
| Cell capacitance | 10 pF | Total capacitance per synthetic cell; no membrane area or human measurement binding |
| Na conductance | 10 pS | Effective linear whole-cell leak, not a single-channel conductance |
| K conductance | 100, 110, 120 pS repeating | Whole-cell leak with invented spatial heterogeneity |
| Cl conductance | 20 pS | Constant leak; omits measured voltage/Ca dependence |
| Intracellular Na/K/Cl | 15/140/30 mM | Fixed reservoirs; not measured keratinocyte state |
| Extracellular Na/K/Cl | 140/4/110 mM | Fixed bath; different from native bulk tissue boundary |
| Temperature | 310.15 K | Fixed 37°C prior; not retained native skin temperature |
| Cell pump | 1 pA outward | Net electrogenic offset; no ATP turnover, ion stoichiometric ledger, saturation, or expression-to-current calibration |
| Gap edge | 10 pS | Effective contact conductance; contact area, connexin state and human pair measurements absent |
| Surface capacitance | 1 µF/site | Aggregate apical-sheet storage; not keratinocyte capacitance; site area undefined |
| Barrier conductance | 1 µS/site | Aggregate through-sheet shunt to basal reference; no TEER/area mapping |
| Surface pump | −30 nA/site | Net current entering apical capacitor is negative; yields apical −30 mV and TEP +30 mV at uniform equilibrium; not equal to the cell pump |
| Surface edge | 1 µS | Lateral conductance, requiring cross-section/length or sheet resistance calibration |
| Geometry | 9 sites, 200 µm spacing, cells 100 µm deep | Synthetic 1.6 mm line; pinned anatomical position does not calibrate thickness/area |
| Wound intervention | Central barrier 100 µS, pump zero | 100-fold shunt, cells retained in assembly experiment; not ablation or healing |
| Electrode intervention | ±30 nA endpoints | Prescribed balanced supply, without measured electrode contact area/impedance |
| Membrane perturbation | +10 mV at cell 0 | Initial condition, not a mechanotransduction law |
| Source uncertainty | Vm initial SD 30 mV; ion SD 20%; surface initial SD 20 mV; OU diffusion amplitudes derived from 20 mV cell and 10 mV surface scales | Declared Gaussian noise, not observed population or parameter uncertainty; assembly reports deterministic means only |

At this audit, `skin-electric.json` records baseline Vm −81.605 to −83.283 mV and uniform apical −30 mV, producing essentially zero baseline lateral field. Local cell `C/sum(g)` is about 0.067–0.077 s before gap contributions, and isolated surface `C/g` is 1 s. These are consequences of chosen coefficients. Equality of a simulated TEP to a reported human range does not calibrate its separate current and resistance.

The two electrical graphs are disconnected under a prescribed basal bath. Surface shunting cannot alter cell Vm; membrane perturbation cannot generate an apical field. `native_skin_ionic_boundary` only reads bulk concentrations and temperature and returns reversals. It neither drives the regional RC patch nor writes ion mass to BioGears. There is no epithelial Ca state, membrane strain dependence, channel gating, ATP balance, or feedback from native perfusion to pump activity.

## Actual retained native boundary

Read from the final unchanged receipt in `data/derived/audits/native-skin-compression-s725ktdm/unchanged/receipts.jsonl`, at 0.4 s and 304.338644 K (31.188644°C). These are simulator observations, not acquired human measurements. Exact values and source hashes are in [the evidence receipt](../../data/sources/skin-bioelectric-calibration-evidence.json).

| Ion | Bulk inside mM | Bulk outside mM | Computed Nernst mV |
| --- | ---: | ---: | ---: |
| Na+ | 15.013157 | 143.831075 | +59.262940 |
| K+ | 115.543993 | 4.469919 | −85.293918 |
| Cl− | 4.794911 | 105.363298 | −81.034261 |

The chloride concentration especially differs from the RC prior (4.795 versus 30 mM inside). Substituting it silently would change model predictions without establishing keratinocyte validity. The native compression audit is a 0.4 s, uniform whole-Skin extracellular pressure experiment; it provides no calibrated local force-to-channel law.

## Primary human evidence and transfer limits

Mauro, Pappone and Isseroff (1990) report average resting Vm **−24 mV** in cultured undifferentiated proliferative human keratinocytes. Raising bath Ca from 0.15 to 2 mM doubled a voltage-gated current and shifted activation. Their abstract describes a major voltage-independent leak and a smaller depolarization-activated chloride current. This motivates differentiation/Ca-specific calibration, but its mean is not a direct target for intact adult forearm. Sample count, uncertainty and complete clamp traces were not acquired here. [Primary study, DOI 10.1002/jcp.1041430103](https://pubmed.ncbi.nlm.nih.gov/1690740/).

Foulds and Barker directly measured transcutaneous voltage in 17 normal volunteers and observed anatomical variation. The accessible abstract supplies no site-level numerical table, so no such data are invented. [Primary human study](https://pubmed.ncbi.nlm.nih.gov/6639877/). Dubé et al. studied human tissue-engineered skin **in vitro** and porcine healing **in vivo**; TEP and Na/K-ATPase expression changed with epithelial formation. Its introductory 10–60 mV human range is background, not a newly acquired participant-level calibration cohort. Expression/TEP correlation does not identify a current per cell. [Primary tissue-engineering study](https://pubmed.ncbi.nlm.nih.gov/20486795/).

The already ingested Nuccitelli human observations support lateral fields and age/site associations. They remain a phenotype model; no channel, pump or RC identification follows from its regression coefficients. [Primary wound-field study](https://pubmed.ncbi.nlm.nih.gov/22092802/).

Morin et al. (2020) provide an acquired, small human hydration table: four volunteers, three forearm sites each. The four human rows of Table 2 are retained as [CSV](evidence/skin-bioelectric/morin-2020-human-hydration.csv), with the table XML and acquisition hash. They report fitted hydration-response rates, **not capacitor charging rates**. Human conductance, MIX and 1 kHz measurements depend on instrument and hydration; the article's four-electrode temperature experiments use pig skin. These data can benchmark hydration-response timing only. [Primary article](https://doi.org/10.1038/s41598-020-73684-y).

| Human instrument observable | Published rate s−1 | Published 95% response time min | R² |
| --- | ---: | ---: | ---: |
| HP conductance | 0.0136 | 4 | 0.991 |
| NE conductance | 0.0051 | 10 | 0.831 |
| NE MIX | 0.0095 | 5 | 0.997 |
| NE 1 kHz | 0.0062 | 8 | 0.764 |

Jemec and Serup found hydration-associated changes in human capacitance, hysteresis and elasticity; capacitance poorly predicted untreated-skin mechanics. This supports controlling hydration in mechanical experiments, not a universal capacitance-to-stiffness mapping. [Primary human mechanics study](https://medicaljournalssweden.se/actadv/article/view/16983).

Human primary keratinocytes exposed to 10% cyclic stretch at 0.5 Hz for 24 h showed Piezo1-associated Ca responses and changes in EMT-related measures. This is culture evidence of mechanical/ionic interaction; it cannot set a millisecond force-to-current coefficient for forearm indentation. [Primary stretch study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9124769/). A separate PIEZO1 study includes human agonist responses but its knockout indentation current comparisons are mouse experiments; species must remain explicit. [Primary mechanosensation study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9512397/).

No inspected primary human dataset here jointly supplies per-cell Na/K/Cl conductances, pump current, capacitance and resting Vm under the same conditions. Reported impedance capacitances of skin barriers and single-channel conductances must not be repurposed as whole-cell parameters.

## Smallest identifiable next model and experiment

Start with **one area-defined epidermal preparation and one aggregate through-layer electrical state**, retaining the basal/apical reference convention. Use a human epidermal equivalent or appropriately maintained human skin preparation; label its culture/excision status. Record area, anatomical origin, layer thickness, temperature, bathing ions/pH, hydration, differentiation and electrode configuration. Do not pool intact human forearm, porcine skin and cultured cells into one parameter likelihood.

Measure paired open-circuit TEP, a small calibrated current-step voltage response, and frequency-dependent impedance on the same preparation. Include blank bath/electrode subtraction and paired controls. Fit an aggregate source-plus-parallel-RC law only if residuals support it:

`C_A d(phi)/dt = j_active - g_A phi + j_applied`,

where `phi=V_apical−V_basal`, `C_A` is F/m², `g_A` is S/m², and currents are A/m². For measured area A, convert to existing node totals as `C=A C_A`, `g=A g_A`, `I=A j`. At equilibrium `phi0=j_active/g_A`; a current-step plateau identifies `g_A=delta_j/delta_phi`; the exponential response identifies `C_A=tau*g_A`; then the effective offset is `j_active=g_A*phi0`. A short-circuit current measurement supplies an independent check. This offset is aggregate active transport, not an identified molecular Na/K pump.

The identifiability reason is concrete: passive voltage relaxation alone only identifies `g/C` and `I/C`; multiplying all three by a common factor preserves the entire trajectory. A baseline TEP alone gives only `I/g`. Current amplitude, area and a measured transient break that ambiguity. If impedance requires a constant-phase element or multiple time constants, retain that model and measurement bandwidth rather than assign a fictitious unique capacitance.

Fit per preparation, retain donor clusters, measurement uncertainty and residuals, and reserve a second current amplitude or hydration condition for prediction. Only then extend to a small lateral sheet with measured geometry and spatial voltage observations to identify lateral conductance and wound-shunt properties. Human wound-field observations can evaluate that observation operator with matching depth and wound geometry; they do not retrospectively identify all conductances.

For a subsequent cellular model, acquire perforated-patch resting Vm, cell capacitance and I–V curves from the same human keratinocyte preparation across controlled Ca/differentiation conditions. Ionic substitution and justified pharmacological perturbations are needed to separate currents; one resting Vm is insufficient. A mechanical extension needs strain/current/time observations at matched stimulus durations with Ca or ATP readouts and controls. Until those measurements exist, keep ionic/mechanical coefficients explicitly unidentified and the native bulk boundary read-only. Neither an RC transient nor electrotaxis alone is evidence of restored anatomy or regeneration.
