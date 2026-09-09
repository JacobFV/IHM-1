"""Preserve local motor gain while declaring smaller cortical drive and float64."""
from pathlib import Path
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import torch,numpy as np
from ihm.native.cortical_stance import load_cortical_stance

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--artifact',required=True);p.add_argument('--output',required=True);p.add_argument('--scale',type=float,default=100);a=p.parse_args()
 if not 1<=a.scale<=1000:raise ValueError('Bounded declared gain scaling required')
 torch.set_num_threads(1);source=Path(a.artifact);policy,artifact=load_cortical_stance(source)
 teacher=ROOT/'data/research/locomotion_control/delay_design_e5s_lepu/candidate_state1.0_R1000.0.npz'
 if hashlib.sha256(teacher.read_bytes()).hexdigest()!=artifact['provenance']['teacher_sha256']:raise ValueError('Teacher identity differs')
 with np.load(teacher,allow_pickle=False) as data:
  x0=data['x0'].copy();u0=data['u0'].copy()
 policy=policy.double()
 with torch.no_grad():
  policy.state_scale.mul_(a.scale);policy.decoder.weight.mul_(a.scale)
  policy.x0.copy_(torch.from_numpy(x0));policy.u0.copy_(torch.from_numpy(u0))
  state=policy.state();weights=policy.dyn.edge_weights()
  nominal,state=policy.advance(policy.x0[None],state,ticks=1000,weights=weights)
  assert torch.equal(nominal[0],policy.u0)
  sever,_=policy.advance(policy.x0[None]+policy.state_scale[None]*.1,policy.state(),ticks=1000,sever=True,weights=weights)
  assert torch.equal(sever[0],policy.u0)
 provenance=dict(artifact['provenance']);provenance.update(parent_artifact_sha256=artifact['artifact_sha256'],computation_dtype='float64',sensory_scale_multiplier=a.scale,decoder_scale_multiplier=a.scale,
  numeric_basis='Same E/I equations and trained embeddings in float64; larger native state scales and inverse decoder scaling preserve local DC gain while reducing nonlinear cortical drive. Exact native x0/u0 restored from bound teacher.')
 report={'schema':'ihm.cortical-stance-regime-scaling.v1','scale':a.scale,'dtype':'float64','nominal_exact_baseline':True,'sever_exact_baseline':True,'parent_artifact_sha256':artifact['artifact_sha256'],'native_stance_demonstrated':False,
  'scope':'Engineering small-signal cortical operating regime; actual nonlinear native verification required'}
 out=Path(a.output);out.mkdir(parents=True,exist_ok=False)
 artifact.update(state_dict=policy.state_dict(),provenance=provenance,report=report)
 artifact.pop('artifact_sha256',None)
 torch.save(artifact,out/'cortical_stance.pt')
 (out/'pretrain_video_loop.py').write_bytes((source.parent/'pretrain_video_loop.py').read_bytes())
 (out/'teacher.npz').write_bytes(teacher.read_bytes())
 (out/'cortical_stance.py').write_bytes((ROOT/'ihm/native/cortical_stance.py').read_bytes())
 (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
 restored,_=load_cortical_stance(out/'cortical_stance.pt');assert restored.x0.dtype==torch.float64 and torch.equal(restored.x0,policy.x0)
 report['artifact_sha256']=hashlib.sha256((out/'cortical_stance.pt').read_bytes()).hexdigest()
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');(out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
