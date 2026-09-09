"""Serial actual native sphere, rigid-skin, and MM support comparison.

This diagnostic reports completion/failure separately; it does not accept equilibrium.
"""
import json,tempfile,time,traceback
from pathlib import Path
import numpy as np
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.assembly.articulated import CanonicalRegistration
root=Path(__file__).resolve().parents[1];output=Path(tempfile.mkdtemp(prefix='surface-support-current-',dir=root/'data/derived'));start=time.monotonic();results=[];manifest='data/derived/supine-surface-contact-5jqy1juo/manifest.json'
import argparse
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--manifest',required=True);args=parser.parse_args();manifest=args.manifest
mechanics=json.loads((root/'data/derived/canonical/mechanics.json').read_text());skinrow=next(e for e in mechanics['entities'] if e['role']=='skin')
import gzip
skin=np.array(json.loads(gzip.decompress((root/skinrow['reference_geometry']['path']).read_bytes()))['positions']).reshape(-1,3);skin=skin[::50]
for name,surface,bed in [('sphere',None,None),('surface-rigid',manifest,None),('surface-MM',manifest,'MM')]:
 stream=None;rows=[];static=None
 try:
  stream=NativeMechanicalStream(root,output/name,environment='supine',target_mass_kg=77.6122029,augmented_registration='data/derived/mechanics/whole_body_lumbar_current/registration.json',surface_contact_manifest=surface,bed_material=bed)
  initial=stream.snapshot();registration=CanonicalRegistration(mechanics,initial);owners=[registration._ranking(p)[0][1] for p in skin];bodyids=list(registration.groups);masks={b:np.array(owners)==b for b in bodyids}
  def skinpoints(state):
   transforms=registration.transforms(state);p=np.zeros_like(skin)
   for b in bodyids:
    t=transforms[b];mask=masks[b];p[mask]=skin[mask]@t[:3,:3].T+t[:3,3]
   return p
  original_skin=skinpoints(initial);previous=original_skin
  static=stream._request('evaluate_static_pose 1 pelvis_tx 0');(output/(name+'-initial-static.json')).write_text(json.dumps(static,indent=2))
  def row(frame):
   global previous
   current=skinpoints(frame);gravity=np.asarray(frame['gravity_m_s2']);normal=-gravity/np.linalg.norm(gravity);f=np.asarray(frame['contact_force_n']);foundation=frame.get('surface_foundation') or {}
   value={'time_s':frame['time_s'],'support_n':float(f@normal),'weight_n':frame['mass_kg']*np.linalg.norm(gravity),'support_fraction':float(f@normal/(frame['mass_kg']*np.linalg.norm(gravity))),'skin_20ms_jump_m':float(np.max(np.linalg.norm(current-previous,axis=1))),'skin_total_displacement_m':float(np.max(np.linalg.norm(current-original_skin,axis=1))),'kinetic_energy_j':frame['kinetic_energy_j'],'pelvis_tx_m':frame['coordinates']['pelvis_tx']['value'],'lumbar_extension_rad':frame['coordinates'].get('lumbar_extension',{}).get('value'),'maximum_skin_penetration_m':foundation.get('maximum_penetration_m'),'maximum_bed_indentation_m':foundation.get('maximum_bed_indentation_m')};previous=current;return value
  rows.append(row(initial))
  for step in range(20):
   frame=stream.advance(.005)
   if step%4==3:rows.append(row(frame))
  result={'arm':name,'completed':True,'rows':rows};results.append(result);print(json.dumps(result),flush=True)
 except BaseException as e:
  result={'arm':name,'completed':False,'error':str(e),'rows':rows};results.append(result);print(json.dumps(result),flush=True)
 finally:
  if stream:stream.close()
  report={'surface_manifest':manifest,'results':results,'wall_s':time.monotonic()-start,'scope':'Current98muscle source defaultpose, serial native sphere vsrigidskin vsMMmattress support. No world/blanket forces, no physiology/brain, no equilibrium acceptance.'};(output/'report.json').write_text(json.dumps(report,indent=2))
print(str(output),flush=True)
