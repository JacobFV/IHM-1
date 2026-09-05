"""heterogeneous system inventory; numbers are broad illustrative initialization.

Each row is a distinct physical component, not a learned embedding. The tuple
schema matches body.DECLARATIONS. Slow anatomy is ordinary structural state.
Reproductive anatomy is an inventory, not an assertion about every subject.
"""
# field, support, region, relaxation seconds, (name, unit, center, std)
SYSTEMS = (
 ('vascular', 'arterial_wall', 'arteries', 30, (
  ('pulse_rate', 'bpm', 70, 25), ('radius', 'mm', 2, 1), ('compliance', 'mL/mmHg', 1, 1), ('smooth_muscle_tone', '1', .5, .4),
  ('endothelial_permeability', 'm/s', 1e-7, 1e-7), ('nitric_oxide', 'umol/L', .1, .1))),
 ('venous', 'venous_tree', 'veins', 30, (
  ('pressure', 'mmHg', 5, 5), ('volume', 'L', 3, 1), ('flow', 'L/min', 5, 3))),
 ('capillary', 'capillary_bed', 'capillaries', 60, (
  ('pressure', 'mmHg', 20, 10), ('filtration', 'mL/min', 2, 2),
  ('recruitment', '1', .5, .4), ('exchange_area', 'm2', 100, 80))),
 ('blood', 'vascular_tree', 'systemic', 1800, (
  ('sodium', 'mmol/L', 140, 10), ('potassium', 'mmol/L', 4, 1),
  ('chloride', 'mmol/L', 103, 10), ('total_co2', 'mmol/L', 25, 10),
  ('calcium', 'mmol/L', 2.4, .4), ('phosphate', 'mmol/L', 1.2, .5),
  ('total_protein', 'g/L', 70, 20), ('uric_acid', 'umol/L', 300, 100),
  ('cholesterol', 'mmol/L', 5, 2), ('hematocrit', '1', .4, .15), ('hemoglobin', 'g/L', 140, 40), ('viscosity', 'mPa.s', 4, 2),
  ('albumin', 'g/L', 40, 15), ('oncotic_pressure', 'mmHg', 25, 10), ('ph', '1', 7.4, .2),
  ('bicarbonate', 'mmol/L', 24, 8), ('platelets', '1e9/L', 250, 150),
  ('fibrinogen', 'g/L', 3, 2), ('thrombin', 'nmol/L', 1, 2),
  ('free_fatty_acids', 'mmol/L', .4, .4), ('triglycerides', 'mmol/L', 1, 1),
  ('bilirubin', 'umol/L', 10, 10), ('urea', 'mmol/L', 5, 3), ('creatinine', 'umol/L', 80, 40))),
 ('interstitial', 'interstitial_space', 'interstitium', 600, (
  ('pressure', 'mmHg', -1, 5), ('volume', 'L', 12, 5), ('protein', 'g/L', 15, 15),
  ('oxygen', 'mmol/L', .04, .03), ('hydraulic_conductivity', 'm/s', 1e-7, 1e-7))),
 ('lymphatic', 'lymphatic_tree', 'lymph_vessels', 120, (
  ('pressure', 'mmHg', 2, 3), ('flow', 'mL/min', 2, 2), ('volume', 'mL', 200, 150),
  ('valve_open_fraction', '1', .5, .4), ('lymphangion_tension', 'N/m', .1, .1),
  ('protein', 'g/L', 20, 15), ('lipid', 'mmol/L', 1, 1))),
 ('lymph_node', 'lymphoid_tissue', 'lymph_nodes', 3600, (
  ('antigen', 'ng/mL', 1, 3), ('dendritic_cells', 'cells/uL', 10, 10),
  ('activated_t_cells', 'cells/uL', 100, 100), ('activated_b_cells', 'cells/uL', 100, 100))),
 ('immune', 'immune_cell_pools', 'immune_tissue', 3600, (
  ('neutrophils', 'cells/uL', 4000, 3000), ('monocytes', 'cells/uL', 500, 400),
  ('macrophages', 'cells/uL', 100, 100), ('t_cells', 'cells/uL', 1000, 800),
  ('b_cells', 'cells/uL', 200, 200), ('nk_cells', 'cells/uL', 200, 200),
  ('il6', 'pg/mL', 2, 10), ('tnf', 'pg/mL', 2, 10), ('complement', 'g/L', 1, 1),
  ('immunoglobulin', 'g/L', 10, 5))),
 ('hematopoietic', 'bone_marrow', 'marrow', 86400, (
  ('erythrocyte_production', 'cells/s', 2e6, 1e6), ('leukocyte_production', 'cells/s', 1e6, 1e6),
  ('platelet_production', 'cells/s', 1e6, 1e6), ('stem_cells', 'cells/uL', 10, 10))),
 ('splenic', 'splenic_pulp', 'spleen', 3600, (
  ('blood_volume', 'mL', 100, 80), ('erythrocyte_clearance', 'cells/s', 2e6, 1e6))),
 ('alveolar', 'alveolar_surface', 'alveoli', 10, (
  ('oxygen_partial_pressure', 'mmHg', 100, 40), ('carbon_dioxide_partial_pressure', 'mmHg', 40, 20),
  ('volume', 'L', 3, 2), ('surfactant', 'mg/mL', 1, 1), ('compliance', 'L/cmH2O', .2, .15))),
 ('airway', 'airway_tree', 'lung', 60, (
  ('resistance', 'cmH2O.s/L', 2, 2), ('mucus_thickness', 'um', 10, 10), ('ciliary_frequency', 'Hz', 10, 5))),
 ('renal', 'nephron', 'kidney', 600, (
  ('sodium_reabsorption', 'mmol/min', 13, 6), ('potassium_secretion', 'mmol/min', .05, .05),
  ('urine_osmolality', 'mmol/kg', 500, 400), ('renin', 'ng/mL', 1, 1))),
 ('urinary', 'urinary_tract', 'bladder', 1800, (
  ('collection_flow', 'mL/min', 1, 1), ('volume', 'mL', 200, 200), ('pressure', 'cmH2O', 10, 10), ('outflow', 'mL/min', 0, 5))),
 ('hepatic', 'hepatic_lobule', 'liver', 1800, (
  ('glycogen', 'g', 80, 60), ('glucose_production', 'mmol/min', .8, .6),
  ('bile_flow', 'mL/min', .5, .5), ('ammonia_clearance', 'mmol/min', .1, .1))),
 ('digestive', 'gut_lumen', 'gut', 1800, (
  ('gastric_volume', 'mL', 200, 200), ('gastric_ph', '1', 2, 2), ('motility', '1/min', 3, 3),
  ('water_absorption', 'mL/min', 5, 5), ('lipid_absorption', 'mmol/min', .1, .1),
  ('amino_acid_absorption', 'mmol/min', .1, .1), ('microbial_biomass', 'g', 100, 100))),
 ('pancreatic', 'pancreatic_tissue', 'pancreas', 600, (
  ('insulin_secretion', 'mU/min', 10, 10), ('glucagon_secretion', 'ng/min', 50, 50),
  ('bicarbonate_secretion', 'mmol/min', .1, .1))),
 ('endocrine', 'vascular_tree', 'systemic', 3600, (
  ('glucagon', 'ng/L', 100, 80), ('aldosterone', 'pmol/L', 300, 250),
  ('vasopressin', 'pg/mL', 2, 2), ('thyroxine', 'pmol/L', 15, 8),
  ('tsh', 'mU/L', 2, 2), ('pth', 'pmol/L', 4, 3), ('growth_hormone', 'ng/mL', 2, 3),
  ('leptin', 'ng/mL', 10, 10), ('adrenaline', 'nmol/L', .2, .2))),
 ('adipose', 'adipose_tissue', 'adipose', 86400, (
  ('lipid_mass', 'kg', 15, 12), ('lipolysis', 'mmol/min', .1, .1), ('heat', 'W', 10, 10))),
 ('bone', 'bone_matrix', 'skeleton', 1e7, (
  ('mineral_density', 'g/cm3', 1, .5), ('osteoblast_density', 'cells/mm3', 100, 100),
  ('osteoclast_density', 'cells/mm3', 10, 10), ('strain', '1', .001, .001))),
 ('joint', 'articular_surface', 'joints', 60, (
  ('angle', 'rad', 0, 1), ('angular_velocity', 'rad/s', 0, 1),
  ('synovial_volume', 'mL', 2, 2), ('cartilage_stress', 'kPa', 100, 100))),
 ('tendon', 'tendon_fibers', 'tendons', 10, (
  ('tension', 'N', 10, 100), ('strain', '1', .01, .02))),
 ('neural', 'central_neural_tissue', 'brain_spinal_cord', 10, (
  ('central_drive', '1', 0, 1), ('sensory_activity', 'Hz', 10, 20))),
 ('sensory', 'sensory_receptors', 'receptors', 1, (
  ('nociceptor_activity', 'Hz', 1, 10), ('mechanoreceptor_activity', 'Hz', 1, 10),
  ('retinal_current', 'pA', 10, 10), ('hair_cell_potential', 'V', -.05, .03))),
 ('csf', 'csf_space', 'brain_spinal_cord', 600, (
  ('pressure', 'mmHg', 10, 8), ('flow', 'mL/min', .3, .3))),
 ('integumentary', 'epidermal_barrier', 'skin', 3600, (
  ('barrier_integrity', '1', .9, .3), ('transepithelial_potential', 'V', .03, .03),
  ('water_loss', 'g/m2/h', 10, 10), ('hydration', '1', .5, .3), ('melanin', 'mg/g', 1, 1))),
 ('dermal', 'dermal_matrix', 'skin', 86400, (
  ('collagen', 'mg/g', 100, 80), ('fibroblast_density', 'cells/mm3', 100, 100),
  ('interstitial_pressure', 'mmHg', 0, 5))),
 ('sweat', 'sweat_ducts', 'skin', 60, (
  ('flow', 'mL/min', .1, .2), ('sodium', 'mmol/L', 40, 30))),
 ('follicular', 'hair_follicles', 'skin', 86400, (
  ('growth_rate', 'mm/day', .3, .3), ('sebum_flow', 'mg/day', 1, 1))),
 ('reproductive', 'gonadal_tissue', 'gonads', 86400, (
  ('estradiol', 'pmol/L', 200, 300), ('progesterone', 'nmol/L', 5, 10),
  ('testosterone', 'nmol/L', 10, 10), ('gamete_production', 'cells/day', 1e6, 1e6))),
 ('uterine', 'uterine_tissue', 'uterus', 86400, (
  ('endometrial_thickness', 'mm', 5, 5), ('contractile_pressure', 'mmHg', 5, 10))),
 ('placental', 'placental_interface', 'placenta', 60, (
  ('blood_flow', 'mL/min', 500, 300), ('oxygen_transfer', 'mL/min', 20, 20))),
)


def declarations():
    return tuple((f'{field}.{name}', support, region, unit, center, std, tau)
                 for field, support, region, tau, components in SYSTEMS
                 for name, unit, center, std in components)

# Deliberately small, explicit couplings. Inventory is not a claim of complete
# dynamics. Interface mechanisms lacking a calibrated form remain prior-only.
COUPLINGS = (
 ('blood.pressure', 'capillary.pressure', .1, 'arterial_capillary'),
 ('venous.pressure', 'capillary.pressure', .5, 'venous_capillary'),
 ('capillary.pressure', 'capillary.filtration', .05, 'capillary_interstitial'),
 ('blood.oncotic_pressure', 'capillary.filtration', -.03, 'capillary_interstitial'),
 ('interstitial.pressure', 'capillary.filtration', -.04, 'capillary_interstitial'),
 ('capillary.filtration', 'interstitial.volume', .1, 'capillary_interstitial'),
 ('interstitial.volume', 'interstitial.pressure', .2, 'interstitial_hydraulic'),
 ('interstitial.pressure', 'lymphatic.flow', .1, 'lymphatic_uptake'),
 ('lymphatic.flow', 'interstitial.volume', -.05, 'lymphatic_uptake'),
 ('lymphatic.flow', 'venous.volume', .0001, 'lymphovenous_return'),
 ('lymphatic.protein', 'lymph_node.antigen', .01, 'lymph_node_transit'),
 ('lymph_node.antigen', 'lymph_node.activated_t_cells', 1, 'antigen_presentation'),
 ('immune.il6', 'immune.crp', .1, 'cytokine_transport'),
 ('immune.crp', 'vascular.endothelial_permeability', 1e-9, 'endothelial_exchange'),
 ('hematopoietic.erythrocyte_production', 'blood.hematocrit', 1e-8, 'marrow_release'),
 ('blood.hematocrit', 'blood.viscosity', 3, 'blood_rheology'),
 ('hepatic.glucose_production', 'metabolic.glucose', .2, 'hepatic_portal'),
 ('digestive.water_absorption', 'blood.volume', .01, 'intestinal_absorption'),
 ('pancreatic.insulin_secretion', 'endocrine.insulin', .1, 'pancreatic_secretion'),
 ('endocrine.aldosterone', 'renal.sodium_reabsorption', .001, 'renal_endocrine'),
 ('endocrine.vasopressin', 'renal.urine_flow', -.05, 'renal_endocrine'),
 ('thermal.core', 'sweat.flow', .1, 'sudomotor'),
 ('sweat.flow', 'thermal.skin', -.1, 'evaporative_exchange'),
 ('neural.central_drive', 'neural.autonomic', .2, 'central_peripheral'),
 ('neural.central_drive', 'effector.activation', .1, 'neuromuscular'),
 ('sensory.nociceptor_activity', 'neural.sensory_activity', .5, 'sensory_afferent'),
)
