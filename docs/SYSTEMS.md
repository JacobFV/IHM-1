# System coverage

The executable scaffold currently has 180 components on 39 supports. These counts are not the coverage of the imported anatomy/physiology corpus; see [the revised design](DIGITAL_HUMAN_DESIGN.md) and [collection summary](../data/derived/collection-summary.json).

| Support | Regions | Physical components |
|---|---|---|
| adipose_tissue | adipose | `adipose.lipid_mass`, `adipose.lipolysis`, `adipose.heat` |
| airway_tree | lung | `respiratory.ventilation`, `respiratory.carbon_dioxide`, `airway.resistance`, `airway.mucus_thickness`, `airway.ciliary_frequency` |
| alveolar_surface | alveoli | `alveolar.oxygen_partial_pressure`, `alveolar.carbon_dioxide_partial_pressure`, `alveolar.volume`, `alveolar.surfactant`, `alveolar.compliance` |
| arterial_wall | arteries | `vascular.pulse_rate`, `vascular.radius`, `vascular.compliance`, `vascular.smooth_muscle_tone`, `vascular.endothelial_permeability`, `vascular.nitric_oxide` |
| articular_surface | joints | `joint.angle`, `joint.angular_velocity`, `joint.synovial_volume`, `joint.cartilage_stress` |
| body_surface | skin | `thermal.skin` |
| bone_marrow | marrow | `hematopoietic.erythrocyte_production`, `hematopoietic.leukocyte_production`, `hematopoietic.platelet_production`, `hematopoietic.stem_cells` |
| bone_matrix | skeleton | `bone.mineral_density`, `bone.osteoblast_density`, `bone.osteoclast_density`, `bone.strain` |
| capillary_bed | capillaries | `capillary.pressure`, `capillary.filtration`, `capillary.recruitment`, `capillary.exchange_area` |
| central_neural_tissue | brain_spinal_cord | `neural.central_drive`, `neural.sensory_activity` |
| csf_space | brain_spinal_cord | `csf.pressure`, `csf.flow` |
| dermal_matrix | skin | `dermal.collagen`, `dermal.fibroblast_density`, `dermal.interstitial_pressure` |
| epidermal_barrier | skin | `integumentary.barrier_integrity`, `integumentary.transepithelial_potential`, `integumentary.water_loss`, `integumentary.hydration`, `integumentary.melanin` |
| gonadal_tissue | gonads | `reproductive.estradiol`, `reproductive.progesterone`, `reproductive.testosterone`, `reproductive.gamete_production` |
| gut_lumen | gut | `digestive.glucose_delivery`, `digestive.gastric_volume`, `digestive.gastric_ph`, `digestive.motility`, `digestive.water_absorption`, `digestive.lipid_absorption`, `digestive.amino_acid_absorption`, `digestive.microbial_biomass` |
| hair_follicles | skin | `follicular.growth_rate`, `follicular.sebum_flow` |
| heart_tissue | heart | `cardiac.rate`, `cardiac.stroke_volume` |
| hepatic_lobule | liver | `hepatic.glycogen`, `hepatic.glucose_production`, `hepatic.bile_flow`, `hepatic.ammonia_clearance` |
| immune_cell_pools | immune_tissue | `immune.neutrophils`, `immune.monocytes`, `immune.macrophages`, `immune.t_cells`, `immune.b_cells`, `immune.nk_cells`, `immune.il6`, `immune.tnf`, `immune.complement`, `immune.immunoglobulin` |
| interstitial | systemic | `metabolic.glucose`, `metabolic.lactate`, `extracellular.sodium`, `extracellular.potassium` |
| interstitial_space | interstitium | `interstitial.pressure`, `interstitial.volume`, `interstitial.protein`, `interstitial.oxygen`, `interstitial.hydraulic_conductivity` |
| lymphatic_tree | lymph_vessels | `lymphatic.pressure`, `lymphatic.flow`, `lymphatic.volume`, `lymphatic.valve_open_fraction`, `lymphatic.lymphangion_tension`, `lymphatic.protein`, `lymphatic.lipid` |
| lymphoid_tissue | lymph_nodes | `lymph_node.antigen`, `lymph_node.dendritic_cells`, `lymph_node.activated_t_cells`, `lymph_node.activated_b_cells` |
| motor_units | muscle | `effector.activation`, `effector.fatigue` |
| musculoskeletal | limb, muscle | `mechanical.force`, `mechanical.velocity` |
| nephron | kidney | `renal.filtration`, `renal.urine_flow`, `renal.sodium_reabsorption`, `renal.potassium_secretion`, `renal.urine_osmolality`, `renal.renin` |
| organ_tissue | core, liver, muscle, systemic | `metabolic.oxygen_consumption`, `hepatic.clearance`, `immune.crp`, `thermal.core`, `structural.muscle_mass` |
| pancreatic_tissue | pancreas | `pancreatic.insulin_secretion`, `pancreatic.glucagon_secretion`, `pancreatic.bicarbonate_secretion` |
| peripheral_nerves | autonomic | `neural.autonomic` |
| placental_interface | placenta | `placental.blood_flow`, `placental.oxygen_transfer` |
| sensor_array | skin, wrist | `device.heart_rate`, `device.glucose`, `device.temperature` |
| sensory_receptors | receptors | `sensory.nociceptor_activity`, `sensory.mechanoreceptor_activity`, `sensory.retinal_current`, `sensory.hair_cell_potential` |
| splenic_pulp | spleen | `splenic.blood_volume`, `splenic.erythrocyte_clearance` |
| sweat_ducts | skin | `sweat.flow`, `sweat.sodium` |
| tendon_fibers | tendons | `tendon.tension`, `tendon.strain` |
| urinary_tract | bladder | `urinary.collection_flow`, `urinary.volume`, `urinary.pressure`, `urinary.outflow` |
| uterine_tissue | uterus | `uterine.endometrial_thickness`, `uterine.contractile_pressure` |
| vascular_tree | systemic | `blood.pressure`, `blood.flow`, `blood.oxygenation`, `blood.volume`, `endocrine.insulin`, `endocrine.cortisol`, `blood.sodium`, `blood.potassium`, `blood.chloride`, `blood.total_co2`, `blood.calcium`, `blood.phosphate`, `blood.total_protein`, `blood.uric_acid`, `blood.cholesterol`, `blood.hematocrit`, `blood.hemoglobin`, `blood.viscosity`, `blood.albumin`, `blood.oncotic_pressure`, `blood.ph`, `blood.bicarbonate`, `blood.platelets`, `blood.fibrinogen`, `blood.thrombin`, `blood.free_fatty_acids`, `blood.triglycerides`, `blood.bilirubin`, `blood.urea`, `blood.creatinine`, `endocrine.glucagon`, `endocrine.aldosterone`, `endocrine.vasopressin`, `endocrine.thyroxine`, `endocrine.tsh`, `endocrine.pth`, `endocrine.growth_hormone`, `endocrine.leptin`, `endocrine.adrenaline` |
| venous_tree | veins | `venous.pressure`, `venous.volume`, `venous.flow` |

The optional NHANES population prior updates 23 physical components from measured data. Dynamic coefficients are unchanged. Source model parameters and geometry are indexed separately until their model context and interface contracts are implemented.
