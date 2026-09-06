"""Check persisted metadata stage plus mask-equivalence rejection, no native jobs."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from stage_skin_layer_candidate import support_equivalence
from prepare_skin_layer_migration import digest, normalized
from ihm.assembly.body import CanonicalBody


def verify(epoch):
 epoch=Path(epoch);report=json.loads((epoch/'acceptance.json').read_bytes())
 for path,sha in report['outputs'].items():
  target=(epoch/path).resolve()
  if not target.is_relative_to(epoch.resolve()) or digest(target.read_bytes())!=sha:raise ValueError('Changed staged output: '+path)
 body=CanonicalBody.from_workspace(epoch/'root')
 mechanics=body.assets['mechanics'];expected=deepcopy(mechanics);normalized(expected)
 if mechanics!=expected:raise ValueError('Staged mass/inertia/damping laws differ')
 if abs(sum(e['mass_kg'] for e in mechanics['entities'])-mechanics['mass_allocation']['target_mass_kg'])>1e-9:raise ValueError('Mass allocation total differs')
 if sum(e['shell']['thickness_m'] for e in mechanics['entities'] if e['role']=='skin_layer')!=.0066:raise ValueError('Changed held shell thickness')
 differences=json.loads((epoch/'differences.json').read_bytes())
 if differences['registration_equal'] is not True or report['native_mass_migrated'] is not False:raise ValueError('Unexpected registration/native claim')
 # A historical full-anatomy hash can differ, but referenced geometry cannot.
 row={'id':'skin','reference_geometry':{'path':'g','sha256':'held'},'bounds_m':{'min':[0,0,0],'max':[1,1,1]},'centroid_m':[.5,.5,.5]}
 original={'entities':[row]};candidate=deepcopy(original);candidate['entities'][0]['volume_m3']=1
 evidence={'source_receipts':[{'entity_id':'skin','path':'g','sha256':'held'}]}
 assert support_equivalence(original,candidate,evidence)['historical_evidence_rewritten'] is False
 candidate['entities'][0]['reference_geometry']['sha256']='changed'
 try:support_equivalence(original,candidate,evidence)
 except ValueError:pass
 else:raise AssertionError('Changed geometry identity accepted')
 return {'outputs_verified':len(report['outputs']),'entities_loaded':len(body.entities),'mass_laws_verified':True,'native_executed':False}

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('epoch',type=Path);args=parser.parse_args();print(json.dumps(verify(args.epoch)))
