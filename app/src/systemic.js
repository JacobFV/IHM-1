const labels={arterial_co2_mmhg:'Arterial CO₂ pressure',arterial_o2_mmhg:'Arterial O₂ pressure',
  respiratory_request_cmh2o:'Requested respiratory pressure',respiratory_applied_cmh2o:'Applied respiratory pressure',
  respiratory_request_per_min:'Requested breathing rate',lung_volume_ml:'Lung gas volume',
  metabolic_rate_w:'Metabolic rate',stomach_carbohydrate_g:'Stomach carbohydrate',stomach_protein_g:'Stomach protein',
  stomach_fat_g:'Stomach fat',stomach_water_ml:'Stomach water',liver_glycogen_g:'Liver glycogen',
  muscle_glycogen_g:'Muscle glycogen',insulin_synthesis_pmol_per_min:'Insulin synthesis',urine_volume_ml:'Urine volume',
  'Aorta.Glucose.concentration_mg_per_dl':'Aorta · glucose concentration'};
export const systemicLabel=name=>labels[name]||String(name).replaceAll('.',' · ').replaceAll('_',' ');

export function prepareSystemic(data,id,label=id) {
  if(data?.model_id!=='ihm-body'||!Array.isArray(data.frames)||!data.frames.length||!data.fields||!data.centroids_m)
    throw Error('Systemic record is missing its canonical body or native fields');
  if(typeof data.has_body_projection!=='boolean')throw Error('Systemic body projection status must be an explicit boolean');
  const dt=data.clock?.sample_interval_s??data.configuration?.sample_interval_s;
  if(!Number.isFinite(dt)||dt<=0)throw Error('Systemic sample clock is invalid');
  const names=Object.keys(data.fields),times=data.frames.map(f=>f.time_s);
  for(let i=0;i<data.frames.length;i++){
    const frame=data.frames[i];
    if(!Number.isFinite(times[i])||(i&&Math.abs(times[i]-times[i-1]-dt)>Math.max(1e-8,dt*1e-7)))
      throw Error('Systemic sample clock is irregular or reversed');
    if(!frame.entities||!frame.systemic_values)throw Error('Systemic frame has no native state');
    if(!data.has_body_projection&&(Object.keys(frame.entities).length||frame.respiration))
      throw Error('Systemic record contains motion without a computed body projection');
    for(const name of names)if(frame.systemic_values[name]!==null&&!Number.isFinite(frame.systemic_values[name]))
      throw Error('Native systemic fields must be numeric or explicitly null');
  }
  const physiology={time_s:times,values:Object.fromEntries(names.map(name=>[name,data.frames.map(f=>f.systemic_values[name])])),
    units:Object.fromEntries(names.map(name=>[name,data.fields[name].unit])),
    metadata:{source_kind:`Native systemic · ${label}`,sample_interval_s:dt}};
  const source=data.spectra||{},frequency=source.frequency_hz||source.laplace?.frequency_hz||[],sigma=source.sigma_per_s||source.laplace?.sigma_per_s||[];
  const variables=source.variables||[];
  const spectralRun={id:`systemic:${id}`,source_kind:`Native systemic · ${label}`,sample_rate_hz:1/dt,
    limitations:[`Recorded ${dt} s samples · Nyquist ${(1/(2*dt)).toPrecision(4)} Hz. Finite-window model descriptors; no inferred physiological poles.`],
    variables:variables.map(v=>({id:v.name,label:systemicLabel(v.name),frequency_hz:v.psd?.frequency_hz||frequency,
      psd:Array.isArray(v.psd)?v.psd:v.psd?.psd||[],psd_unit:`(${data.fields[v.name]?.unit||'source unit'})²/Hz`})),
    laplace:{frequency_hz:frequency,sigma_per_s:sigma,variables:variables.map(v=>({id:v.name,real:v.laplace_real||[],imag:v.laplace_imag||[],
      unit:`${data.fields[v.name]?.unit||'source unit'} · s`,removed_sample_mean:v.removed_sample_mean}))}};
  return {trajectory:data,physiology,spectralRun,label,id};
}

export function systemicSnapshot(data,index,name) {
  const value=data?.frames?.[index]?.systemic_values?.[name];
  return {name,label:systemicLabel(name),value:Number.isFinite(value)?value:null,unit:data?.fields?.[name]?.unit||'',
    time_s:data?.frames?.[index]?.time_s??null};
}
