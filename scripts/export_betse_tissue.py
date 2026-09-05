"""Export the reviewed, locally generated BETSE run to safe portable arrays.

Use .venv-physiology/bin/python. Only the hash-pinned local solver output below
may be unpickled. The web API consumes JSON, never pickle or executable objects.
"""
from pathlib import Path
import hashlib
import gzip
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'data/derived/physiology/betse_run'
EXPECTED='e278f33ba55ef1dbcb9928352cadffd4afea1a58a02b3b3eab021ecf1004b773'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
 path=RUN/'SIMS/sim_1.betse.gz'
 if sha(path)!=EXPECTED:raise ValueError('Refuse to unpickle unreviewed solver artifact')
 from betse.science.filehandling import loadSim
 sim,cells,p=loadSim(str(path))
 times=np.asarray(sim.time,float);voltage=np.asarray(sim.vm_ave_time,float)
 concentration=np.asarray(sim.cc_time,float);gating=np.asarray(sim.gjopen_time,float)
 assert voltage.shape==(len(times),212) and concentration.shape==(len(times),4,212)
 assert np.all(np.diff(times)>0) and np.isfinite(voltage).all() and np.isfinite(concentration).all()
 mem_ids=np.asarray(cells.mem_to_cells,int);counts=np.bincount(mem_ids,minlength=212)
 gate=np.array([np.bincount(mem_ids,weights=v,minlength=212)/counts for v in gating])
 assert np.isfinite(gate).all() and gate.min()>=0 and gate.max()<=1
 field_values={'Vmem':voltage,'gap_junction_open':gate}
 units={'Vmem':'V','gap_junction_open':'1'}
 for i,name in sim.ionlabel.items():field_values[name]=concentration[:,i,:];units[name]='mol/m^3'
 fields=[dict(id=name,label=name.replace('_',' '),unit=units[name],range=[float(v.min()),float(v.max())]) for name,v in field_values.items()]
 centers=np.asarray(cells.cell_centres,float);vertices=np.concatenate(cells.cell_verts).astype(float)
 center=(vertices.min(0)+vertices.max(0))/2;scale=2/np.ptp(vertices,axis=0).max()
 positions=[];indices=[];cell_ids=[];source_polygons=[]
 for i,polygon in enumerate(cells.cell_verts):
  polygon=np.asarray(polygon,float);source_polygons.append(polygon.tolist())
  points=np.vstack([centers[i],polygon]);start=len(positions)//3
  points=(points-center)*scale
  positions.extend(np.c_[points,np.zeros(len(points))].ravel().tolist());cell_ids.extend([i]*len(points))
  for j in range(len(polygon)):indices.extend([start,start+1+j,start+1+(j+1)%len(polygon)])
 frames=[dict(values={name:values[i].tolist() for name,values in field_values.items()}) for i in range(len(times))]
 metadata=json.loads((RUN/'run_summary.json').read_text())
 provenance=dict(label='BETSE 1.5.1 generic tissue simulation',url='https://github.com/betsee/betse',sha256=EXPECTED,
  frame='betse-generated-planar-tissue',specimen='generic computational tissue; no human subject',units='m',status='executed source simulation',
  source_revision=metadata['revision'],configuration_sha256=sha(RUN/'sim_config.yaml'))
 limitations=[
  'Actual planar solver cells shown in 3D coordinates; no invented depth or human skin registration.',
  'Generic computational tissue, not a human wound experiment. This 0.035-second simulation does not validate wound healing.',
  'Vmem is native cell-averaged transmembrane voltage, distinct from human extracellular lateral field in V/m.',
  'Cell concentration units follow native mol/m^3 (numerically mM); proteins/anion are source lumped species.',
  'Gap-junction open fraction is the mean of incident membrane entries; unequal cell membrane counts prevent interpreting edge-average voltage as cell-average voltage.',
  'Initial random geometry was not seeded reproducibly; the exact generated geometry, configuration and simulation are retained by hash.']
 geometry=dict(positions=positions,indices=indices,cell_ids=cell_ids,scalar_fields=fields,times=times.tolist(),frames=frames,
  source_centers_m=centers.tolist(),source_polygons_m=source_polygons,display_transform=dict(scale=float(scale),translation=(-center*scale).tolist()),
  source_membrane_to_cell=mem_ids.tolist(),source_membrane_voltages_V=np.asarray(sim.vm_time).tolist(),limitations=limitations)
 out=ROOT/'data/derived/bioelectric';out.mkdir(parents=True,exist_ok=True)
 (out/'tissue.json').write_text(json.dumps(dict(schema_version=1,source_kind='source_model_simulation',source=provenance,cells=212,membrane_edges=len(mem_ids),time_s=times.tolist(),fields=fields,frames=frames,limitations=limitations),separators=(',',':'),allow_nan=False))
 gp=ROOT/'data/derived/app/geometry/betse-tissue-cells.json.gz'
 gp.write_bytes(gzip.compress(json.dumps(geometry,separators=(',',':'),allow_nan=False).encode(),mtime=0))
 mp=ROOT/'data/derived/app/manifest.json';manifest=json.loads(mp.read_text())
 manifest['models']=[m for m in manifest['models'] if m['id']!='betse-tissue']
 manifest['structures']=[m for m in manifest['structures'] if m['model_id']!='betse-tissue']
 model=dict(id='betse-tissue',name='BETSE · cellular bioelectricity',description='212 actual planar solver cells with membrane voltage, ionic concentrations and junction gating. Generic tissue; not human skin calibration.',
  frame='betse-tissue-display',source_frame=provenance['frame'],source_units='m',display_units='normalized',
  bounds=dict(min=[-1,-1,0],max=[1,1,0]),calibration_status='source simulator execution; no human validation',limitations=limitations)
 manifest['models'].append(model)
 if not any(s['id']=='bioelectric' for s in manifest['systems']):manifest['systems'].append(dict(id='bioelectric',name='Cell bioelectricity',color='#a187d5'))
 manifest['structures'].append(dict(id='betse-tissue-cells',name='212 cells · native bioelectric fields',model_id='betse-tissue',system='bioelectric',kind='scalar_mesh',
  geometry_url='/api/geometry/betse-tissue-cells',color='#a187d5',source=provenance,geometry_sha256=sha(gp),default_visible=True,calibration_status='source simulation',limitations=limitations))
 mp.write_text(json.dumps(manifest,separators=(',',':'),allow_nan=False))
 print(json.dumps(dict(cells=212,membrane_edges=len(mem_ids),frames=len(times),cell_mean_voltage_V=float(voltage[-1].mean()),fields=fields)))
if __name__=='__main__':build()
