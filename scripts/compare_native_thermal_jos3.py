"""Execute the held JOS-3 equations under matched generic external conditions.

Independent source physiological state and setpoints are deliberately retained.
This is a model comparison, not fitting JOS-3 to BioGears or human validation.
"""
from pathlib import Path
import argparse,json,sys,tempfile
import numpy as np
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from ihm.native.thermal import load_source,thermal_clock,heat_step_audit
from scripts.audit_long_horizon_thermal import sha
from ihm.assembly.systemic_evidence import freeze_sources

def run(seconds=21600,dt=15):
    steps=thermal_clock(seconds,dt)
    package,provenance=load_source(BASE)
    profile_path=BASE/'data/derived/canonical/profile.json'
    profile=json.loads(profile_path.read_text())
    out=Path(tempfile.mkdtemp(prefix='thermal-matched-jos3-',dir=BASE/'data/derived/audits'))
    sources={str(p.relative_to(BASE)):sha(p) for p in (Path(__file__),profile_path,BASE/'ihm/native/thermal.py')}
    sources.update({'data/raw/thermal/JOS-3/'+p:h for p,h in provenance['source_sha256'].items()})
    freeze_sources(BASE,out,sources)
    for posture in ('standing','lying'):
        model=package.JOS3(height=profile['height_m'],weight=profile['mass_kg'],age=profile['age_years'],fat=100*profile['body_fat_fraction'],sex=profile['sex'],ex_output='all')
        model.posture=posture;model.Ta=22;model.Tr=22;model.Va=.1;model.RH=60;model.Icl=.5;model.PAR=1
        frames=[dict(time_s=0,central_blood_c=float(model._bodytemp[0]),mean_skin_c=float(model.TskMean))]
        audits=[];previous_state=model._bodytemp.copy();previous_profile=sys.getprofile()
        def observe(frame,event,arg):
            if frame.f_code is model._run.__func__.__code__ and event=='return' and arg is not None:
                local=frame.f_locals;capacity=model._cap;step=local['dtime']
                transfer=(local['arr_bf']+local['arr_cdt'])*capacity[:,None]/step
                boundary=local['arrB']*capacity/step;heat=local['arrQ']*capacity/step
                audits.append(heat_step_audit(capacity,transfer,boundary,local['arr_to'],heat,previous_state,model._bodytemp,step))
            if previous_profile:previous_profile(frame,event,arg)
        try:
            sys.setprofile(observe)
            for i in range(steps):
                model.simulate(1,dtime=dt);previous_state=model._bodytemp.copy();last=model._history[-1]
                frames.append(dict(time_s=(i+1)*dt,central_blood_c=float(model._bodytemp[0]),mean_skin_c=float(model.TskMean),metabolic_w=float(last['Met']),respiratory_w=float(last['RES'])))
        finally:sys.setprofile(previous_profile)
        assert len(audits)==steps
        assert all(np.isfinite(list(row.values())).all() for row in frames)
        maxima={k:max(row[k] for row in audits) for k in ('heat_balance_residual_W','node_equation_residual_max_W','internal_column_sum_max_W_K')}
        assert max(maxima.values())<1e-7
        report=dict(source=provenance,source_hashes=sources,generic_profile=profile,
            configuration=dict(posture=posture,air_temperature_c=22,radiant_temperature_c=22,air_speed_m_s=.1,relative_humidity_percent=60,uniform_local_clothing_clo=.5,physical_activity_ratio=1,seconds=seconds,dt_s=dt),
            frames=frames,heat_step_audits=audits,maximum_equation_residuals=maxima,
            limitations=['Independent JOS-3 setpoints and physiological state; central blood is a different observable from the BioGears lumped core.',
                'Uniform local clothing prior; no mattress or bedding, including in the lying comparison.',
                'Conservation and model comparison do not establish human measurement agreement.'])
        (out/(posture+'.json')).write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'posture':posture,'final':frames[-1],'maximum_equation_residuals':maxima,'output_dir':str(out)}),flush=True)
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=21600);parser.add_argument('--dt',type=float,default=15);args=parser.parse_args();run(args.seconds,args.dt)
