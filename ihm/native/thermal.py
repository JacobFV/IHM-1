"""Pinned authors' JOS-3 execution and explicit Akimoto 2025 bedding boundary.

The published solver is imported from a verified source checkout, without changing
its thermal equations. A context-managed boundary substitution implements the
paper's total-resistance inputs; it is identified separately from source defaults.
"""
from pathlib import Path
import contextlib,hashlib,importlib.util,json,subprocess,sys,threading
import numpy as np

REVISION='3c74ee2af2f79aa360093cc517e38bf465ec8c5b'
SOURCE_URL='https://github.com/TanabeLab/JOS-3'
BEDDING_DOI='10.1016/j.buildenv.2025.113074'
_LOCK=threading.RLock()
MAX_DURATION_S=604800
MIN_STEP_S=.01
MAX_STEP_S=3600
MAX_STEPS=100000

def thermal_clock(seconds,dt):
    """Bound work and reject invalid clocks before constructing the source model."""
    if any(isinstance(v,(bool,np.bool_)) or not isinstance(v,(int,float,np.integer,np.floating)) or not np.isfinite(v) for v in (seconds,dt)):
        raise ValueError('duration and timestep must be finite numbers, not booleans')
    if not 0<seconds<=MAX_DURATION_S or not MIN_STEP_S<=dt<=MAX_STEP_S:
        raise ValueError('duration must be in (0, 604800] s; timestep in [0.01, 3600] s')
    steps=round(seconds/dt)
    if not 1<=steps<=MAX_STEPS or abs(seconds/dt-steps)>1e-9:
        raise ValueError('duration must contain an integer 1..100000 timesteps')
    return steps

def total_resistance(total_insulation_clo):
    value=np.asarray(total_insulation_clo,float)
    if not np.isfinite(value).all() or (value<=0).any():raise ValueError('positive finite total insulation required')
    rt=.155*value
    return rt,rt/(.38*16.5)

def heat_step_audit(capacity,transfer,boundary,operative,heat,old,new,dt):
    """Exact backward-Euler ledger; transfer[i,j] is heat into i from j (W/K)."""
    c,w,b,to,q,t0,t1=map(lambda x:np.asarray(x,float),(capacity,transfer,boundary,operative,heat,old,new))
    if c.ndim!=1 or not c.size:raise ValueError('nonempty one-dimensional heat capacities required')
    n=len(c)
    if dt<=0 or not np.isfinite(dt) or w.shape!=(n,n) or any(x.shape!=(n,) for x in [b,to,q,t0,t1]) or not all(np.isfinite(x).all() for x in [c,w,b,to,q,t0,t1]) or (c<=0).any() or (w<0).any() or (b<0).any():raise ValueError('invalid finite heat system')
    with np.errstate(over='raise',invalid='raise',divide='raise'):
        k=np.diag(w.sum(axis=1))-w
        storage=c*(t1-t0)/dt
        boundary_heat=b*(to-t1)
        residual=storage+k@t1-boundary_heat-q
        result=dict(heat_balance_residual_W=float(abs(storage.sum()-boundary_heat.sum()-q.sum())),
            node_equation_residual_max_W=float(abs(residual).max()),internal_column_sum_max_W_K=float(abs(k.sum(axis=0)).max()),
            storage_W=float(storage.sum()),boundary_heat_into_body_W=float(boundary_heat.sum()),net_metabolic_minus_evaporative_respiratory_W=float(q.sum()))
    if not np.isfinite(list(result.values())).all():raise FloatingPointError('nonfinite heat ledger')
    return result

def load_source(root):
    root=Path(root);folder=root/'data/raw/thermal/JOS-3'
    metadata=json.loads((root/'data/raw/thermal/source.json').read_text())
    if metadata['revision']!=REVISION:raise ValueError('unexpected JOS-3 revision')
    actual=subprocess.check_output(['git','-C',str(folder),'rev-parse','HEAD'],text=True).strip()
    changed=subprocess.check_output(['git','-C',str(folder),'status','--porcelain','--untracked-files=all'],text=True).strip()
    if actual!=REVISION or changed:raise ValueError('JOS-3 checkout must be clean at the pinned revision')
    required={str(p.relative_to(folder)) for p in (folder/'src/jos3').glob('*.py')}|{'LICENSE','README.md'}
    if set(metadata['source_sha256'])!=required:raise ValueError('incomplete JOS-3 source hash inventory')
    for relative,expected in metadata['source_sha256'].items():
        if hashlib.sha256((folder/relative).read_bytes()).hexdigest()!=expected:raise ValueError('JOS-3 source changed: '+relative)
    name='jos3'
    if name in sys.modules and Path(sys.modules[name].__file__).resolve() != (folder/'src/jos3/__init__.py').resolve():raise ValueError('A different JOS-3 package is already loaded')
    if name not in sys.modules:
        spec=importlib.util.spec_from_file_location(name,folder/'src/jos3/__init__.py',submodule_search_locations=[str(folder/'src/jos3')])
        module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    # Keep historical acquisition/artifact provenance unchanged; identify this execution separately.
    provenance={**metadata,'runtime_adapter_path':'ihm/native/thermal.py',
        'runtime_adapter_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    return sys.modules[name],provenance

@contextlib.contextmanager
def bedding_boundary(package,total_insulation):
    """Paper Eq7/8 boundary, under a lock; original source functions restored."""
    with _LOCK:
        th=sys.modules[package.__name__+'.thermoregulation'];original=(th.dry_r,th.wet_r)
        rt,ret=total_resistance(total_insulation)
        th.dry_r=lambda *a,**kw:np.full(17,rt)
        th.wet_r=lambda *a,**kw:np.full(17,ret)
        try:yield
        finally:th.dry_r,th.wet_r=original

PROFILES={
 'lying_default':dict(label='JOS-3 · lying · source default environment',paper=False),
 'supine_mattress':dict(label='JOS-3 · supine · nude mattress/pillow · 22.6 °C',paper=True,total_insulation_clo=1.15,coverage_percent=23.3,sleepwear='nude',bedding='mattress and pillow'),
 'supine_blanket':dict(label='JOS-3 · supine · pajamas and blanket · 22.6 °C',paper=True,total_insulation_clo=2.78,coverage_percent=94.1,sleepwear='pajamas',bedding='blanket'),
 'supine_duvet':dict(label='JOS-3 · supine · pajamas and duvet · 22.6 °C',paper=True,total_insulation_clo=5.46,coverage_percent=94.1,sleepwear='pajamas',bedding='duvet'),
}

def run_thermal(root,profile='lying_default',seconds=3600,dt=30):
    if not isinstance(profile,str) or profile not in PROFILES:raise ValueError('unknown thermal profile')
    steps=thermal_clock(seconds,dt)
    with _LOCK:
        package,provenance=load_source(root);config=PROFILES[profile]
        model=package.JOS3(**(dict(height=1.68,weight=61.5,fat=15,age=20,sex='female',ci=2.59) if config['paper'] else {}),ex_output='all')
        model.posture='lying'
        if config['paper']:model.To=22.6;model.PAR=1.0
        initial=model._bodytemp.copy();history=[];audits=[];states=[initial.copy()];captured={}
        previous_profile=sys.getprofile()
        def observe(frame,event,arg):
            if frame.f_code is model._run.__func__.__code__ and event=='return' and arg is not None:
                loc=frame.f_locals;c=model._cap;d=loc['dtime']
                w=(loc['arr_bf']+loc['arr_cdt'])*c[:,None]/d
                b=loc['arrB']*c/d;q=loc['arrQ']*c/d
                audit=heat_step_audit(c,w,b,loc['arr_to'],q,states[-1],model._bodytemp,d)
                audits.append(audit)
                captured.update(transfer_W_K=w.copy(),boundary_W_K=b.copy(),operative_C=loc['arr_to'].copy(),heat_W=q.copy(),hc_W_m2K=loc['hc'].copy(),hr_W_m2K=loc['hr'].copy(),dry_resistance_m2K_W=loc['r_t'].copy(),evaporative_resistance_m2kPa_W=loc['r_et'].copy())
            if previous_profile:previous_profile(frame,event,arg)
        context=bedding_boundary(package,config['total_insulation_clo']) if config['paper'] else contextlib.nullcontext()
        with context:
            try:
                sys.setprofile(observe)
                for _ in range(steps):
                    model.simulate(1,dtime=dt);history.append(model._history[-1]);states.append(model._bodytemp.copy())
            finally:sys.setprofile(previous_profile)
        states=np.asarray(states);matrix=sys.modules[package.__name__+'.matrix']
        channels=[]
        def add(id,unit,values):channels.append(dict(id=id,unit=unit,values=np.asarray(values,float).tolist()))
        add('MeanSkinTemperature','degC',[np.average(row[matrix.INDEX['skin']],weights=model._bsa) for row in states])
        add('CentralBloodTemperature','degC',states[:,0])
        for region,idx in zip(matrix.BODY_NAMES,matrix.INDEX['core']):add('Core_'+region,'degC',states[:,idx])
        for region,idx in zip(matrix.BODY_NAMES,matrix.INDEX['skin']):add('Skin_'+region,'degC',states[:,idx])
        # Fluxes label first sample missing: initialization belongs to a different boundary.
        for key,unit in [('Met','W'),('RES','W'),('CO','L/h')]:add(key,unit,[np.nan]+[r[key] for r in history])
        for ch in channels:
            if ch['values'][0]!=ch['values'][0]:ch['values'][0]=None
        coefficient=dict(node_count=len(model._cap),body_regions=matrix.BODY_NAMES,node_index=matrix.IDICT,capacity_J_K=model._cap.tolist(),conductance_W_K=model._cdt.tolist(),body_surface_area_m2=model._bsa.tolist(),last_step={k:v.tolist() for k,v in captured.items()})
        audit=dict(maximum_heat_balance_residual_W=max(x['heat_balance_residual_W'] for x in audits),maximum_node_equation_residual_W=max(x['node_equation_residual_max_W'] for x in audits),maximum_internal_column_sum_W_K=max(x['internal_column_sum_max_W_K'] for x in audits),steps=audits)
        limitations=["Separate published thermoregulation model, not a matched BioGears patient or clinically calibrated bed-rest prediction.","Lying is the source supine heat-transfer posture. No lateral position, organ hydrostatics, body geometry registration or sleep-stage dynamics is inferred.","Initialization uses the authors' standing neutral setpoint procedure; at t=0 posture/environment changes. The one-hour transient is not an equilibrated sleeper.","Thermal coefficients and flows are source-model parameters; conservation and timestep refinement do not establish human measurement agreement.","Source reported sensible heat loss uses the previous temperature; the audit uses implicit-step endpoint temperature to match the solved heat equation."]
        if config['paper']:limitations += ["Whole-body measured total insulation from Table9 is applied uniformly to17regions, as one paper comparison method. Regional measurements were not supplied; local bedding/contact variation is not reconstructed.","Explicit bedding boundary adapter uses RT=0.155 IT and ReT=RT/(0.38*16.5), not clothing Icl. The paper prose calls RT reciprocal IT, but Eq2/7 and dimensions establish the conversion used here.","The paper simulated over2000minutes to steady state; this one-hour transient is not a reproduction of its final skin temperatures."]
        return dict(id=profile,label=config['label'],source_kind='source_model_simulation',time_s=(np.arange(len(states))*dt).tolist(),channels=channels,configuration={**config,'seconds':seconds,'dt_s':dt,'posture':model.posture,'air_temperature_C':model.Ta.tolist(),'radiant_temperature_C':model.Tr.tolist(),'relative_humidity_percent':model.RH.tolist(),'air_speed_m_s':model.Va.tolist(),'physical_activity_ratio':model.PAR,'height_m':model._height,'weight_kg':model._weight,'sex':model._sex,'cardiac_index_L_min_m2':model._ci},source=provenance,limitations=limitations,audit=audit,coefficients=coefficient,node_temperature_C=states.tolist())
