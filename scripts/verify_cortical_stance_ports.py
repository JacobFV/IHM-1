#!/usr/bin/env python3
"""Bounded signed cortical port audit using actual E/I equations; no native job."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import torch
from ihm.native.cortical_stance import CorticalStancePolicy,load_cortical_stance


def main():
 torch.set_num_threads(1);torch.manual_seed(47)
 legacy_path=ROOT/'data/runtime/motor-learning/cortical-stance-patient-slow-20260908/cortical_stance.pt'
 legacy,prior=load_cortical_stance(legacy_path,dt_s=.01)
 assert legacy.dyn.n==512 and legacy.encoder_kind=='random_projection'
 raw=(legacy_path.parent/'pretrain_video_loop.py').read_bytes()
 ns={'__name__':'_port_verification_core'};exec(compile(raw,'pretrain_video_loop.py','exec'),ns)
 nstate,nmuscle=len(legacy.state_names),len(legacy.muscle_names)
 assert nstate==258 and nmuscle==98
 def construct(sensory=None,motor=None):
  dyn=ns['CorticalDynamics'](1024,legacy.dyn.embed.shape[1],legacy.dyn.idx.shape[1],'cpu')
  return CorticalStancePolicy(dyn,torch.arange(2*nstate) if sensory is None else sensory,
   torch.arange(2*nstate,2*nstate+nmuscle+1) if motor is None else motor,
   state_names=legacy.state_names,muscle_names=legacy.muscle_names,x0=legacy.x0,u0=legacy.u0,
   state_scale=legacy.state_scale,reference_normalization=True,encoder_kind='signed_identity')
 policy=construct();expected=torch.cat((torch.eye(nstate),-torch.eye(nstate)))
 assert torch.equal(policy.encoder,expected)
 assert int(torch.linalg.matrix_rank(policy.encoder))==258
 invalid=[]
 for name,sensory,motor in [('short_sensory',torch.arange(515),None),('long_sensory',torch.arange(517),torch.arange(517,616)),
  ('short_motor',None,torch.arange(516,614)),('overlapping',None,torch.arange(515,614))]:
  try:construct(sensory,motor)
  except ValueError:invalid.append(name)
  else:raise AssertionError('Invalid ports accepted: '+name)
 with torch.no_grad():
  nominal=policy.state();sever=policy.state();weights=policy.dyn.edge_weights()
  for index in range(30):
   zero,nominal=policy.advance(policy.x0[None],nominal,ticks=10,weights=weights)
   assert torch.equal(zero[0],policy.u0)
   vector=policy.x0+policy.state_scale*torch.randn(nstate)*10
   output,sever=policy.advance(vector[None],sever,ticks=10,sever=True,weights=weights)
   assert torch.equal(output[0],policy.u0)
 provenance=copy.deepcopy(prior['provenance']);provenance.update(encoder_kind='signed_identity',reference_normalization=True)
 artifact={'schema':'ihm.ibm-cortical-stance.v1','provenance':provenance,'state_dict':policy.state_dict()}
 errors=[]
 with tempfile.TemporaryDirectory(prefix='ihm-signed-port-audit-') as directory:
  directory=Path(directory);(directory/'pretrain_video_loop.py').write_bytes(raw);path=directory/'model.pt'
  torch.save(artifact,path);restored,_=load_cortical_stance(path,dt_s=.01)
  assert restored.encoder_kind=='signed_identity' and restored.dyn.n==1024 and torch.equal(restored.encoder,expected)
  assert all(torch.equal(value,restored.state_dict()[key]) for key,value in policy.state_dict().items())
  for name,mutate in [('matrix_value',lambda x:x.__setitem__((0,0),.5)),('matrix_shape',None)]:
   changed=copy.deepcopy(artifact)
   if mutate:mutate(changed['state_dict']['encoder'])
   else:changed['state_dict']['encoder']=changed['state_dict']['encoder'][:-1]
   torch.save(changed,path)
   try:load_cortical_stance(path,dt_s=.01)
   except (ValueError,RuntimeError):pass
   else:errors.append('Accepted corrupt signed encoder: '+name)
 report={'schema':'ihm.cortical-signed-ports-verification.v1','passed':not errors,'errors':errors,
  'actual_cortical_sites':1024,'state_width':nstate,'encoder_rank':258,'sensory_width':516,'motor_width':99,
  'nominal_exact_baseline_ticks':300,'arbitrary_state_sever_exact_baseline_ticks':300,
  'invalid_ports_rejected':invalid,'roundtrip_all_tensors_exact':True,'legacy_512_load_compatible':True,
  'source_sha256':{'script':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
   'policy':hashlib.sha256((ROOT/'ihm/native/cortical_stance.py').read_bytes()).hexdigest()},
  'scope':'Actual E/I equation constructor and port contract only; no training or native mechanics'}
 out=ROOT/'data/runtime/motor-learning/cortical_stance_ports_verification.json';out.write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))
 if errors:raise SystemExit(1)

if __name__=='__main__':main()
