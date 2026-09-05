# Human skin bioelectric measurements

The [Nuccitelli et al. 2011 clinical study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3228273/)
is downloaded with a SHA-256 receipt. `scripts/collect_skin_evidence.py` extracts
Table I into `data/derived/integumentary/`: 160 table cells (155 numeric), 40 human
participants, five author-marked outliers and 16 published group means/SEMs.
Participant IDs are only unique within age/sex groups; repeated sites and scans
are not independent subjects. Individual raw scans are unavailable.

The measured quantity is a lateral electric field above the epidermis and below
the stratum corneum. It is neither cellular membrane voltage nor the deeper
intraepidermal field. Numerically, mV/mm equals V/m. Values retain site, age, sex,
missingness, uncertainty and the authors' exclusion flags. The paper's table labels
the younger group 18–29; its methods say 18–25. This discrepancy is explicitly
preserved. Most per-person entries summarize three scans with SEM; the exact
replicate count is not known for each cell.

These observations constrain a future measurement operator and boundary-field
prediction. They do not identify individual channel conductances or validate a
whole-body wound-healing model. The source reports no correlation between field
amplitude and wound-healing rate in its mouse experiments. Species and measured
layer must remain explicit when extending the evidence.

The generic BETSE run and these human observations are separate artifacts. No
fit between them is claimed. Next identification work needs tissue geometry,
barrier properties, ion concentrations, channel/pump kinetics and perturbation
measurements with compatible experimental conditions.
