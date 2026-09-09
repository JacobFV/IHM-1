"""Actual native mechanics/physiology regression, with retained energy receipts."""
from pathlib import Path
import argparse,json,math,tempfile,time
from ihm.assembly.embodied import EmbodiedRuntime


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seconds',type=float,default=5.)
    parser.add_argument('--environment',choices=('supine','upright','free'),default='supine')
    args=parser.parse_args()
    if not math.isfinite(args.seconds) or not .02<=args.seconds<=120 or abs(args.seconds/.02-round(args.seconds/.02))>1e-8:
        parser.error('--seconds must be a multiple of .02 in [.02,120]')
    root=Path(__file__).resolve().parents[1]
    output=Path(tempfile.mkdtemp(prefix='sustained-metabolic-',dir=root/'data/derived'))
    started=time.monotonic();body=None;rows=[]
    report={'passed':False,'environment':args.environment,'requested_duration_s':args.seconds,
        'scope':'Actual native mechanics, selected baseline controller and native physiology; bounded runtime and energy accounting, not gait or biological calibration'}
    try:
        body=EmbodiedRuntime.from_workspace(root,output/'body',environment=args.environment)
        report['metabolic_reference']=body.metabolic_reference
        totals={key:0. for key in ('m_raw','h_raw','w_raw','m_sent','h_sent','w_sent')}
        for _ in range(round(args.seconds/.02)):
            frame=body.step({});c=frame['coupling'];dt=c['exchange_interval_s']
            for owner in (frame['mechanics']['time_s'],frame['physiology']['elapsed_s']):
                assert abs(owner-frame['time_s'])<1e-8,'Common clock diverged'
            for channel,raw,sent in (
                ('m','native_extra_metabolic_demand_w','exchanged_metabolic_increment_w'),
                ('h','muscle_heat_increment_w','exchanged_heat_increment_w'),
                ('w','signed_work_increment_w','exchanged_work_increment_w')):
                totals[channel+'_raw']+=c[raw]*dt;totals[channel+'_sent']+=c[sent]*dt
                residual=totals[channel+'_raw']-totals[channel+'_sent']-c['metabolic_pending_energy_j'][channel+'_j']
                assert abs(residual)<1e-8*(1+abs(totals[channel+'_raw'])),'Unaccounted lag energy'
            rows.append({'time_s':frame['time_s'],'raw_m_w':c['native_extra_metabolic_demand_w'],
                'exchanged_m_w':c['exchanged_metabolic_increment_w'],
                'pending_energy_j':c['metabolic_pending_energy_j']})
        report.update(passed=True,duration_s=frame['time_s'],energy_totals_j=totals,
            exchanged_range_w=[min(row['exchanged_m_w'] for row in rows),max(row['exchanged_m_w'] for row in rows)])
    except BaseException as error:
        report.update(error=str(error),duration_s=0 if not rows else rows[-1]['time_s'])
        raise
    finally:
        if body is not None:body.close()
        report['wall_s']=time.monotonic()-started
        (output/'verification.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        (output/'energy_trace.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
        print(json.dumps({'output':str(output),**report},indent=2),flush=True)


if __name__=='__main__':main()
