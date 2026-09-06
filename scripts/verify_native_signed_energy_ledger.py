"""Tiny actual native ledger/replay fixture; requires the coordinated heavy slot."""
from pathlib import Path
import argparse,json,math,resource,sys,tempfile,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
ROOT=Path(__file__).resolve().parents[1]
ENERGIES=('muscle_metabolic_energy_j','signed_active_fiber_work_j','muscle_heat_energy_j','positive_active_fiber_work_j')

def close(a,b):
    assert math.isfinite(a) and math.isfinite(b)
    assert abs(a-b)<=128*sys.float_info.epsilon*max(1.,abs(a),abs(b)),(a,b)

def check_interval(before,after,*,endpoint_check=True):
    dt=after['time_s']-before['time_s'];assert dt>0
    ref=before['metabolic_reference'];assert after['metabolic_reference']==ref
    m,w,h=(after[k]-before[k] for k in ENERGIES[:3])
    close(m,h+w);close(ref['M0_w'],ref['H0_w']+ref['W0_w'])
    delta_m=m/dt-ref['M0_w'];delta_w=w/dt-ref['W0_w'];delta_h=h/dt-ref['H0_w']
    close(delta_m,delta_h+delta_w)
    if endpoint_check:
        close(m,.5*dt*(before['total_muscle_metabolic_w']+after['total_muscle_metabolic_w']))
        close(w,.5*dt*(before['signed_active_fiber_power_w']+after['signed_active_fiber_power_w']))
    return {'dt_s':dt,'metabolic_energy_j':m,'signed_work_j':w,'heat_energy_j':h,'DeltaM_w':delta_m,'DeltaW_w':delta_w,'DeltaH_w':delta_h}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-native',action='store_true');args=parser.parse_args()
    if not args.run_native:raise SystemExit('Requires explicitly coordinated --run-native resource slot')
    output=Path(tempfile.mkdtemp(prefix='signed-energy-ledger-',dir=ROOT/'data/derived'));started=time.monotonic()
    plant=NativeMechanicalStream(ROOT,output/'native',environment='supine',target_mass_kg=77.6122029,augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json')
    try:
        initial=plant.snapshot();assert all(initial[k]==0 for k in ENERGIES)
        first=plant.advance(.002);intervals=[check_interval(initial,first)]
        checkpoint=plant.checkpoint();second=plant.advance(.002);intervals.append(check_interval(first,second))
        restored=plant.restore(checkpoint);assert all(restored[k]==first[k] for k in (*ENERGIES,'metabolic_reference'))
        replay=plant.advance(.002);assert all(replay[k]==second[k] for k in (*ENERGIES,'metabolic_reference','bodies'))
        # Fail after changing an excitation to exercise the native catch rollback.
        try:plant._request('advance .002 0 1 arm26_BRA_r .8 trailing')
        except ValueError as error:assert 'trailing command data' in str(error)
        else:raise AssertionError('Invalid native command accepted')
        rejected=plant._request('observe');assert all(rejected[k]==second[k] for k in (*ENERGIES,'metabolic_reference','muscles','time_s'))
        plant.release(checkpoint)
        for name,frame in [('initial',initial),('first',first),('second',second),('replay',replay),('after_rejection',rejected)]:
            (output/(name+'.json')).write_text(json.dumps(frame,allow_nan=False)+'\n')
        plant.close()
        report={'passed':True,'intervals':intervals,'reference':initial['metabolic_reference'],'positive_work_j':second['positive_active_fiber_work_j'],'signed_work_j':second['signed_active_fiber_work_j'],'wall_s':time.monotonic()-started,'maximum_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'scope':'Native endpoint quadrature and cumulative signed ledger, nonzero checkpoint replay and rejected-command rollback; no absolute chemical-to-thermal supply closure claim'}
        (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({**report,'output':str(output)},indent=2))
    finally:plant.close()
if __name__=='__main__':main()
