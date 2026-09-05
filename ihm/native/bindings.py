"""Observation bindings to original engine systems; no inferred anatomical registration."""
SYSTEM_VARIABLES={
 'Cardiovascular': ['ArterialPressure','HeartRate','MeanArterialPressure','SystolicArterialPressure','DiastolicArterialPressure','CardiacOutput','BloodVolume'],
 'Respiratory':['TotalLungVolume','RespirationRate','TidalVolume'],
 'BloodChemistry':['OxygenSaturation','ArterialCarbonDioxidePressure','ArterialOxygenPressure','ArterialBloodPH'],
 'Renal':['GlomerularFiltrationRate','RenalBloodFlow','UrineProductionRate'],
 'Endocrine':['InsulinSynthesisRate'],
 'Energy':['CoreTemperature','SkinTemperature','TotalMetabolicRate','SweatRate']}
def variable_bindings(columns,revision):
    result={}
    for column in columns:
        name=column.split('(')[0]
        system=next((s for s,names in SYSTEM_VARIABLES.items() if name in names),None)
        unit=column.split('(',1)[1].rstrip(')') if '(' in column else '1'
        result[column]={'quantity':name,'unit':unit,'engine_system':system,'source_path':f'projects/biogears/libBiogears/src/engine/Systems/{system}.cpp' if system else None,'source_revision':revision,'evidence_kind':'source_simulation' if system else 'engine_diagnostic','anatomical_registration':None,'uncertainty':None,'independently_calibrated':False}
    return result
