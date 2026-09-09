"""Native muscle stretch -> named spinal arc -> native activation at1ms exchange."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.assembly.ibm_controller import IBMImplicitController
from ihm.assembly.sensorimotor_catalog import native_muscle_catalog
import torch


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    args=parser.parse_args();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(1)
    controller=IBMImplicitController.from_root(ROOT,muscle_catalog=native_muscle_catalog(ROOT),sites=128)
    ids=controller.muscles;index=ids.index('soleus_r');initial_cord=controller.checkpoint()
    plant=NativeMechanicalStream(ROOT,out/'plant',environment='free',target_mass_kg=70)
    initial=plant.snapshot();token=plant.checkpoint();records={}
    # One experimental evoked spindle channel: native soleus normalized fiber
    # length referenced to its actual initial native length. This avoids tonic
    # baseline afference obscuring the external perturbation's latency.
    base=initial['muscles']['soleus_r']['fiber_length_m']
    try:
        for arm in ('pulse','pulse_stretch_blocked','sham'):
            plant.restore(token);controller.restore(initial_cord);cord=controller.cord
            rows=[];next_command={m:0. for m in ids}
            for step in range(80):
                state=plant.snapshot();time=state['time_s']
                length=state['muscles']['soleus_r']['fiber_length_m']
                ia_hz=max(0.,(length/base-1)*200.)
                stretch=np.zeros(len(ids),np.float32)
                if arm!='pulse_stretch_blocked':stretch[index]=min(1.,ia_hz/100.)
                response=cord.step(np.zeros(len(ids),np.float32),stretch=stretch)
                command={m:float(response['alpha'][i]) for i,m in enumerate(ids)}
                # Output produced at this neural endpoint is applied on the next
                # explicit mechanics exchange, just as the unified runtime does.
                forces=[]
                if arm!='sham' and .005-1e-12<=time<.015-1e-12:
                    point=plant.body_point(body='calcn_r',station_m=[.15,0.,0.])['point_source_m']
                    forces=[{'body':'calcn_r','point_m':point,'force_n':[0.,300.,0.]}]
                endpoint=plant.advance(.001,forces=forces,actuation=next_command)
                rows.append({'sample_time_s':time,'endpoint_time_s':endpoint['time_s'],
                    'fiber_length_m':length,'evoked_ia_hz':ia_hz,
                    'stretch_contribution':float(response['stretch'][index]),
                    'renshaw_contribution':float(response['renshaw'][index]),
                    'computed_alpha':command['soleus_r'],'applied_excitation':next_command['soleus_r'],
                    'native_activation':endpoint['muscles']['soleus_r']['activation'],
                    'force_pulse_applied':bool(forces)})
                next_command=command
            records[arm]=rows
            (out/'partial.json').write_text(json.dumps(records,indent=2)+'\n')
            print(json.dumps({'arm':arm,'ia_peak_hz':max(r['evoked_ia_hz'] for r in rows),
                'stretch_peak':max(r['stretch_contribution'] for r in rows)}),flush=True)
    finally:plant.release(token);plant.close()
    pulse=records['pulse'];sham=records['sham'];blocked=records['pulse_stretch_blocked']
    def first_difference(field,left,right,tolerance=1e-9):
        return next((a['sample_time_s'] for a,b in zip(left,right) if abs(a[field]-b[field])>tolerance),None)
    ia=first_difference('evoked_ia_hz',pulse,sham)
    named=first_difference('stretch_contribution',pulse,sham)
    alpha=first_difference('computed_alpha',pulse,sham)
    excitation=first_difference('applied_excitation',pulse,sham)
    activation=next((a['endpoint_time_s'] for a,b in zip(pulse,sham) if abs(a['native_activation']-b['native_activation'])>1e-9),None)
    measured=None if ia is None or named is None else named-ia
    passed=measured is not None and abs(measured-.03)<1e-8 and abs(alpha-named)<1e-8 and excitation>=alpha+.001-1e-8
    assert not any(r['stretch_contribution'] for r in blocked)
    receipt={'schema':'ihm.native-reflex-latency.v1','passed':passed,
        'dt_s':.001,'muscle':'soleus_r','initial_native_fiber_length_m':base,
        'first_pulse_specific_ia_s':ia,'first_named_stretch_s':named,'first_alpha_s':alpha,
        'first_applied_excitation_interval_start_s':excitation,'first_native_activation_endpoint_s':activation,
        'named_stretch_delay_s':measured,'records':records,
        'source_sha256':{'verifier':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'IBM_cord':controller.identity['cord_sha256'],'native_stream':hashlib.sha256((ROOT/'ihm/native/mechanical_stream.py').read_bytes()).hexdigest()},
        'scope':['Actual native dynamic fiber length and external300N foot-force pulse; free body',
            'Evoked soleus spindle increment200Hz per reference strain; initial native fiber length is reference, no fabricated deformation',
            'Fixed zero descending request isolates named stretch; other afferent channels zero; Renshaw preserved',
            'Native mechanics plus actual IBM SegmentalCord; no1ms physiology or full cortical-control claim',
            'Matched restored-state pulse/sham differences separate pulse from passive free-body motion']}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('records','source_sha256','scope')},indent=2),flush=True)
    if not passed:raise SystemExit('Native pulse did not demonstrate required reflex latency; inspect receipt')

if __name__=='__main__':main()
