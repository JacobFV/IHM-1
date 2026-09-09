#!/usr/bin/env python3
"""Native adapter smoke verification, distinct from sustained balance acceptance."""
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.stance_lqr_controller import EngineeredLQRStanceController
from ihm.native.mechanical_stream import NativeMechanicalStream
ROOT=Path(__file__).resolve().parents[1]
def run(artifact,registration,output):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    reg=json.loads(Path(registration).read_bytes());catalog=json.loads((ROOT/reg['catalog_path']).read_bytes())
    with np.load(artifact) as d:mass=float(d['target_mass_kg'])
    controller=EngineeredLQRStanceController.from_root(ROOT,muscle_catalog=catalog,artifact_path=artifact,registration_path=registration)
    native=NativeMechanicalStream(ROOT,out/'native',environment='upright',target_mass_kg=mass,augmented_registration=str(Path(registration).resolve().relative_to(ROOT)))
    controller.bind_native(native);controller.retain_sources(out)
    checkpoint=native.checkpoint();neural=controller.checkpoint();records=[]
    try:
        for arm in ('feedback','sensory_block','motor_block'):
            native.restore(checkpoint);controller.restore(neural)
            for step in range(3):
                state=native.snapshot();obs={'time_s':state['time_s'],'joints':state['coordinates'],'muscles':state['muscles'],'foot_contact_force_n':state['foot_contact_force_n']}
                result=controller.step(.01,obs,sensory_blocks=[controller.muscles[0]] if arm=='sensory_block' else [],motor_blocks=[controller.muscles[0]] if arm=='motor_block' else [])
                commands=result['motor_excitations']
                if arm=='sensory_block' and any(commands.values()):raise AssertionError('Sensory block leaked privileged feedback')
                if arm=='motor_block' and commands[controller.muscles[0]]!=0.:raise AssertionError('Motor block did not zero output')
                final=native.advance(.01,actuation=commands)
                records.append({'arm':arm,'time_s':final['time_s'],'commands':commands,'diagnostics':result['lqr_stance']})
    finally:native.release(checkpoint);native.close()
    report={'schema':'ihm.engineered-lqr-adapter-native-verification.v1','controller':controller.controller_metadata,'records':records,
        'passed':True,'scope':'Three30ms matched native protocol arms; does not establish sustained balance or walking'}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'output':str(out/'report.json'),'passed':True}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--artifact',required=True);p.add_argument('--registration',required=True);p.add_argument('--output',required=True);a=p.parse_args();run(a.artifact,a.registration,a.output)
