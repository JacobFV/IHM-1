"""Source circuit responses and conditional physical posture fields."""
from pathlib import Path
import json
import numpy as np
from ihm.human import digest
from ihm.coupling import hydrostatic_pressure

ROOT=Path(__file__).resolve().parents[1]
def build(root=ROOT):
 root=Path(root);derived=root/'data/derived'
 graph=json.loads((derived/'native-circuits/graph.json').read_text())
 circuit=json.loads((derived/'coupling/native-skin-circuit.json').read_text())
 anatomy=json.loads((derived/'anatomy/bodyparts3d_index.json').read_text())
 systems={s['type']:s['properties'] for s in graph['systems']}
 density=systems['BioGearsBloodChemistrySystemData']['BloodDensity']
 assert density['unit']=='kg/m^3' and density['value']>0
 pressure=systems['BioGearsCardiovascularSystemData']['MeanArterialPressure']
 assert pressure['si_unit']=='Pa'
 # Reference geometry is the complete ventricular wall bounding-box midpoint.
 heart=[m for m in anatomy['meshes'] if m['element_id']=='FJ2428']
 assert len(heart)==1
 ref=np.mean(heart[0]['bounds_in_source_coordinates'],axis=0)*.001
 points=np.array([np.mean(m['bounds_in_source_coordinates'],axis=0) for m in anatomy['meshes']])*.001
 positions={name:hydrostatic_pressure(points,pressure['si_value'],ref,density['value'],gravity)
            for name,gravity in [('supine',[0,9.80665,0]),('upright',[0,0,-9.80665])]}
 hydro=dict(kind='conditional_static_hydrostatic_field',equation='p(x) = p_ref + rho g dot (x-x_ref)',
  source_frame='bodyparts3d-mm-left-posterior-superior',position_unit='m',pressure_unit='Pa',
  reference_position_m=ref.tolist(),reference_landmark='FJ2428 ventricular wall bounding-box midpoint; geometric reference, not measured pressure catheter',
  density_kg_m3=density['value'],reference_pressure_pa=pressure['si_value'],
  coefficients_source=graph['source'],gravity_m_s2={'supine':[0,9.80665,0],'upright':[0,0,-9.80665]},
  structures=[dict(id='bp3d-'+m['element_id'],position_m=points[i].tolist(),**{name:dict(pressure_pa=float(p[i]),offset_pa=float(p[i]-pressure['si_value'])) for name,p in positions.items()}) for i,m in enumerate(anatomy['meshes'])],
  limitations=['Conditional static fluid column calculation using native blood density and arterial reference pressure on separate reference anatomy.',
    'No cross-specimen registration is implied; no regional perfusion, vessel collapse, viscous losses, venous pooling, baroreflex or native posture dynamics are predicted.',
    'All structure midpoints sample one hypothetical connected arterial pressure field; these values are not organ, venous or tissue pressures.'])
 path=derived/'coupling/hydrostatic.json';path.write_text(json.dumps(hydro,separators=(',',':'),allow_nan=False))
 poles=circuit['model']['finite_poles_per_s']
 result=dict(schema_version=1,source=graph['source'],summary=dict(nodes=len(graph['nodes']),paths=len(graph['paths']),compartments=len(graph['compartments']),
  native_coupling='Executed inside BioGears; source circuit topology and initialized coefficients retained',
  cross_source_coupling='Explicit physical operators available; geometry is not automatically coupled across specimens'),
  models=[dict(id=circuit['id'],kind=circuit['model']['kind'],path='data/derived/coupling/native-skin-circuit.json',sha256=digest(derived/'coupling/native-skin-circuit.json'),
    nodes=len(circuit['model']['node_names']),paths=len(circuit['model']['native_paths']),finite_poles_per_s=poles,
    time_constants_s=[-1/r for r,i in zip(poles['real'],poles['imag']) if r<0 and abs(i)<1e-10],
    zero_modes=sum(abs(r)+abs(i)<1e-12 for r,i in zip(poles['real'],poles['imag'])),limitations=circuit['limitations']),
    dict(id='posture-hydrostatic',kind=hydro['kind'],path='data/derived/coupling/hydrostatic.json',sha256=digest(path),structures=len(hydro['structures']),limitations=hydro['limitations'])],
  operators=['atomic conservative volume/mass/ion/charge/energy exchange','positive conservative implicit ion transport','electrodiffusive Scharfetter-Gummel flux','implicit reciprocal hydraulic network','native frozen hydraulic descriptor Laplace response'],
  limitations=['Source simulation consistency does not establish independent human calibration.', 'Field measurement above epidermis does not identify cell membrane conductances.'])
 (derived/'coupling/index.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 return result
if __name__=='__main__':print(json.dumps(build()['summary'],indent=2))
