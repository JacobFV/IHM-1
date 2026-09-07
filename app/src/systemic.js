const labels={arterial_co2_mmhg:'Arterial CO₂ pressure',arterial_o2_mmhg:'Arterial O₂ pressure',
  respiratory_request_cmh2o:'Requested respiratory pressure',respiratory_applied_cmh2o:'Applied respiratory pressure',
  respiratory_request_per_min:'Requested breathing rate',lung_volume_ml:'Lung gas volume',
  metabolic_rate_w:'Metabolic rate',stomach_carbohydrate_g:'Stomach carbohydrate',stomach_protein_g:'Stomach protein',
  stomach_fat_g:'Stomach fat',stomach_water_ml:'Stomach water',liver_glycogen_g:'Liver glycogen',
  muscle_glycogen_g:'Muscle glycogen',insulin_synthesis_pmol_per_min:'Insulin synthesis',urine_volume_ml:'Urine volume',
  'Aorta.Glucose.concentration_mg_per_dl':'Aorta · glucose concentration'};
export const systemicLabel=name=>labels[name]||String(name).replaceAll('.',' · ').replaceAll('_',' ');
