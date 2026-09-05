"""Ursino–Lodi 1997 two-state CSF/cerebrovascular source model, Appendix A/B."""
from pathlib import Path
import csv,hashlib,json
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root as find_root

SOURCE_URL='https://www.researchgate.net/publication/14112614_A_simple_mathematical_model_of_the_interaction_between_intracranial_pressure_and_cerebral_hemodynamics'
REVIEW_SHA256='69a989fb59cbca3e2317a6553ded045dbed340ff4558177ea6ae523d23b386ea'
MMHG_PA=133.322387415
PARAMETERS={
    'Ro':(526.3,'mmHg s/mL'),'Rpv':(1.24,'mmHg s/mL'),'Rf':(2380.,'mmHg s/mL'),
    'delta_Ca1':(.75,'mL/mmHg'),'delta_Ca2':(.075,'mL/mmHg'),'Can':(.15,'mL/mmHg'),
    'kE':(.11,'1/mL'),'kR':(49100.,'mmHg^3 s/mL'),'tau':(20.,'s'),
    'qn':(12.5,'mL/s'),'G':(1.5,'mL/mmHg per fractional CBF change'),'Pvs':(6.,'mmHg')}
LIMITATIONS=[
    'Source simulation of a separate literature model, not a calibrated or matched BioGears patient.',
    'The reduced model assumes arterial pressure exceeds ICP and ICP exceeds sinus pressure (collapsed terminal-vein Starling regime); it excludes venous compliance and rapid cardiac/respiratory pulsatility.',
    'Use exact Appendix A3 capillary conservation including CSF production, rather than the optional approximate B1 used to simplify computation.',
    'Literal Appendix A9/A10 and Table 1 imply saturation Ca=0.525/0.1125 mL/mmHg; prose describing sixfold/half basal limits is inconsistent. Equation/table implementation is explicit and not silently retuned.',
    'Rounded Table 1 initial values are approximately, not exactly, at equilibrium; the computed equilibrium and resulting baseline drift are reported.',
    'Initial ICP/compliance remain the rounded published basal state when replacing the arterial boundary; the resulting transient is not a matched native initial condition.',
    'Native MAP coupling is one-way prescribed forcing with source unit/hash and interpolation; subjects are unmatched and ICP is never fed back into the native engine.',
    'Primary source author-upload indexed text was acquired and hashed as a retrieval record; no hash of unavailable original PDF bytes is claimed.',
    'Numerical and source-target checks do not establish clinical ICP accuracy, patient identifiability, glymphatic transport or posture-specific venous mechanics.'
]


class SourceEvidenceUnavailable(RuntimeError):pass


def classify_source_response(data,expected):
    prefix=data[:20000].lower()
    if any(x in prefix for x in (b'just a moment',b'cf-chl-',b'checking your browser',b'temporarily unavailable',b'access denied')):return 'access_challenge'
    if expected=='pdf':return 'pdf_acquired_requires_review' if data.startswith(b'%PDF-') and b'%%EOF' in data[-2048:] else 'wrong_content_type'
    if expected=='html':return 'html_acquired_requires_review' if b'<html' in prefix or b'<!doctype html' in prefix else 'wrong_content_type'
    raise ValueError('expected source must be PDF or HTML')


def _source_review(root):
    path=Path(root)/'data/raw/csf/ursino_lodi_1997/indexed-author-source-review.json'
    if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest()!=REVIEW_SHA256:
        raise SourceEvidenceUnavailable('Reviewed original author-upload Table 1 and Appendices A–C retrieval is absent or changed; source-faithful execution requires verified evidence.')
    acquisition=path.parent/'acquisition.json'
    return dict(url=SOURCE_URL,doi='10.1152/jappl.1997.82.4.1256',authors=['Mauro Ursino','Carlo Alberto Lodi'],
        citation='J Appl Physiol 82:1256–1269 (1997)',source_location='Table 1 and Appendices A–C, author-upload full text',
        review_path=str(path.relative_to(root)),review_sha256=REVIEW_SHA256,
        hash_scope='indexed web-tool source retrieval, not original PDF bytes',copyright='1997 American Physiological Society; no open reuse license inferred',
        acquisition_attempts=json.loads(acquisition.read_text()) if acquisition.exists() else None)


@dataclass
class PressureInput:
    time_s: np.ndarray | None
    pressure_mmHg: np.ndarray
    provenance: dict

    @classmethod
    def constant(cls,value_mmHg=100.):
        if not np.isfinite(value_mmHg) or value_mmHg<=0:raise ValueError('positive finite arterial pressure required')
        return cls(None,np.array([float(value_mmHg)]),dict(kind='source_published_constant' if value_mmHg==100 else 'explicit_constant_override',value=value_mmHg,unit='mmHg',source='Table 1' if value_mmHg==100 else 'caller boundary override'))

    @classmethod
    def from_samples(cls,time_s,values,*,unit,provenance=None):
        t=np.asarray(time_s,float);p=np.asarray(values,float)
        if unit not in ('mmHg','Pa') or t.ndim!=1 or len(t)<2 or t.shape!=p.shape or not np.isfinite(t).all() or not np.isfinite(p).all() or np.any(np.diff(t)<=0) or np.any(p<=0):raise ValueError('finite increasing seconds and positive pressure with mmHg/Pa units required')
        p=p/MMHG_PA if unit=='Pa' else p
        return cls(t,p,dict(provenance or {},source_unit=unit,interpolation='piecewise linear; analytic segment slope in mmHg/s; no extrapolation'))

    @classmethod
    def from_native_map(cls,path):
        path=Path(path).resolve()
        with path.open() as f:
            reader=csv.reader(f);header=next(reader);data=np.array([[float(v) for v in row] for row in reader])
        if header[0]!='#Time(s)' or 'MeanArterialPressure(mmHg)' not in header:raise ValueError('native MAP requires explicit seconds and mmHg CSV headers')
        time=data[:,0];pressure=data[:,header.index('MeanArterialPressure(mmHg)')]
        return cls.from_samples(time-time[0],pressure,unit='mmHg',provenance=dict(kind='native_MAP_one_way_forcing',
            path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),source_start_s=float(time[0]),
            input_quantity='MeanArterialPressure',source_samples=len(time),coupling='one way; source patients unmatched; no ICP feedback'))

    def evaluate(self,time):
        if not np.isfinite(time):raise ValueError('finite time required')
        if self.time_s is None:return float(self.pressure_mmHg[0]),0.
        t=self.time_s
        if time<t[0]-1e-9 or time>t[-1]+1e-9:raise ValueError('pressure interpolation refuses extrapolation')
        i=int(np.clip(np.searchsorted(t,time,side='right')-1,0,len(t)-2))
        slope=(self.pressure_mmHg[i+1]-self.pressure_mmHg[i])/(t[i+1]-t[i])
        return float(self.pressure_mmHg[i]+slope*(time-t[i])),float(slope)


class CSFModel:
    def __init__(self,parameters=None):
        self.parameters={k:v[0] for k,v in PARAMETERS.items()}
        if parameters is not None:
            if not set(parameters)<=set(self.parameters):raise ValueError('unknown source parameter')
            self.parameters.update(parameters)
        if not all(np.isfinite(v) and v>0 for v in self.parameters.values()):raise ValueError('source parameters must be positive finite quantities')
        if self.parameters['Can']<=self.parameters['delta_Ca2']/2:raise ValueError('source sigmoid would admit nonpositive compliance')

    def regulated_compliance(self,x):
        p=self.parameters;delta=p['delta_Ca1'] if x<0 else p['delta_Ca2']
        # Algebraically identical to A9 with ks=delta/4; tanh avoids overflow.
        return p['Can']-delta/2*np.tanh(2*p['G']*x/delta)

    def evaluate(self,state,Pa,dPa_dt=0.,injection_ml_s=0.):
        state=np.asarray(state,float)
        if state.shape!=(2,) or not np.isfinite(state).all() or not all(np.isfinite(v) for v in (Pa,dPa_dt,injection_ml_s)):raise ValueError('finite two-state and boundary inputs required')
        Pic,Ca=state;p=self.parameters
        if not Pa>Pic>p['Pvs'] or Ca<=0:raise ValueError('outside source collapsed-vein regime: require Pa > Pic > Pvs and Ca > 0')
        Va=Ca*(Pa-Pic);Ra=p['kR']*p['Can']**2/Va**2
        # A3 exact capillary balance: arterial inflow splits into CSF + venous flow.
        parallel=1/p['Rf']+1/p['Rpv'];Pc=(Pa/Ra+Pic*parallel)/(1/Ra+parallel)
        q=(Pa-Pc)/Ra;qf=max(0.,(Pc-Pic)/p['Rf']);qo=max(0.,(Pic-p['Pvs'])/p['Ro']);qv=(Pc-Pic)/p['Rpv']
        x=(q-p['qn'])/p['qn'];dc=(self.regulated_compliance(x)-Ca)/p['tau']
        Cic=1/(p['kE']*Pic)
        dp=(Ca*dPa_dt+dc*(Pa-Pic)+qf-qo+injection_ml_s)/(Cic+Ca)
        dVa=Ca*(dPa_dt-dp)+dc*(Pa-Pic)
        return dict(Pic_mmHg=float(Pic),Ca_ml_per_mmHg=float(Ca),Pa_mmHg=float(Pa),Pc_mmHg=float(Pc),
            Va_ml=float(Va),Ra_mmHg_s_per_ml=float(Ra),Cic_ml_per_mmHg=float(Cic),q_ml_s=float(q),qf_ml_s=float(qf),qo_ml_s=float(qo),
            qv_ml_s=float(qv),dPic_mmHg_s=float(dp),dCa_ml_per_mmHg_s=float(dc),dVa_ml_s=float(dVa),
            capillary_balance_ml_s=float(q-qf-qv),intracranial_balance_ml_s=float(Cic*dp-dVa-qf+qo-injection_ml_s))

    def linearization(self,Pa=100.):
        def rhs(state):
            d=self.evaluate(state,Pa);return np.array([d['dPic_mmHg_s'],d['dCa_ml_per_mmHg_s']])
        result=find_root(rhs,[9.5,self.parameters['Can']],tol=1e-11)
        if not result.success or np.linalg.norm(rhs(result.x))>1e-8:raise RuntimeError('source equilibrium solve failed')
        steps=[1e-4,1e-7];jac=np.column_stack([(rhs(result.x+np.eye(2)[i]*h)-rhs(result.x-np.eye(2)[i]*h))/(2*h) for i,h in enumerate(steps)])
        poles=np.linalg.eigvals(jac)
        return dict(arterial_pressure_mmHg=Pa,equilibrium_state=result.x.tolist(),jacobian_native_units=jac.tolist(),eigenvalues_real=poles.real.tolist(),eigenvalues_imag=poles.imag.tolist(),
            note='finite-difference local linearization of this source model; stability is operating-point dependent')


def run_csf(root,*,pressure_input=None,duration_s=600.,samples=601,method='DOP853',parameters=None):
    source=_source_review(root)
    if not np.isfinite(duration_s) or not 0<duration_s<=3600 or not isinstance(samples,int) or not 2<=samples<=10001 or method not in ('DOP853','Radau'):raise ValueError('invalid bounded integration settings')
    port=pressure_input or PressureInput.constant();model=CSFModel(parameters);times=np.linspace(0,float(duration_s),samples)
    if port.time_s is not None and (port.time_s[0]>0 or duration_s>port.time_s[-1]+1e-9):raise ValueError('simulation requires a supplied pressure value throughout the horizon')
    def rhs(time,state):
        pa,rate=port.evaluate(time);d=model.evaluate(state,pa,rate);return [d['dPic_mmHg_s'],d['dCa_ml_per_mmHg_s']]
    # Resolve every interval of a sampled source; its derivative is discontinuous
    # at knots. The direct native MAP has a short, bounded 50 Hz time basis.
    max_step=min(.5,float(np.min(np.diff(port.time_s)))) if port.time_s is not None else .5
    sol=solve_ivp(rhs,(0,float(duration_s)),[9.5,.15],t_eval=times,method=method,rtol=1e-9,atol=1e-11,max_step=max_step)
    if not sol.success or not np.isfinite(sol.y).all():raise RuntimeError('source CSF integration failed: '+sol.message)
    observations=[model.evaluate(state,*port.evaluate(t)) for t,state in zip(times,sol.y.T)]
    metadata=[('Pic_mmHg','mmHg','Pa',MMHG_PA),('Ca_ml_per_mmHg','mL/mmHg','m^3/Pa',1e-6/MMHG_PA),
        ('Pa_mmHg','mmHg','Pa',MMHG_PA),('Pc_mmHg','mmHg','Pa',MMHG_PA),('Va_ml','mL','m^3',1e-6),
        ('q_ml_s','mL/s','m^3/s',1e-6),('qf_ml_s','mL/s','m^3/s',1e-6),('qo_ml_s','mL/s','m^3/s',1e-6)]
    channels=[dict(id=name,unit=unit,si_unit=si,si_scale=scale,values=[d[name] for d in observations]) for name,unit,si,scale in metadata]
    return dict(schema_version=1,id='ursino_lodi_1997',system='cerebral_csf',source_kind='source_model_simulation',execution_available=True,
        source=source,parameters=[dict(id=k,value=v,unit=PARAMETERS[k][1],source='Table 1',status='source basal value' if v==PARAMETERS[k][0] else 'explicit sensitivity override') for k,v in model.parameters.items()],
        time_s=times.tolist(),state_values=sol.y.tolist(),channels=channels,pressure_input=port.provenance,
        initial_state=dict(Pic_mmHg=9.5,Ca_ml_per_mmHg=.15,source='Table 1, rounded source initialization'),
        solver=dict(method=method,rtol=1e-9,atol=1e-11,max_step_s=max_step,nfev=sol.nfev),
        validation=dict(max_capillary_balance_ml_s=max(abs(d['capillary_balance_ml_s']) for d in observations),
            max_intracranial_balance_ml_s=max(abs(d['intracranial_balance_ml_s']) for d in observations)),
        linearization=model.linearization(),limitations=LIMITATIONS)


def csf_source_status(root):
    try:source=_source_review(root)
    except SourceEvidenceUnavailable:return dict(schema_version=1,id='ursino_lodi_1997',execution_available=False,parameters=[],missing_evidence=['Reviewed author-upload Table 1 and Appendices A–C retrieval'])
    return dict(schema_version=1,id='ursino_lodi_1997',execution_available=True,source=source,limitations=LIMITATIONS)


def build_csf_status(root):
    root=Path(root).resolve();out=root/'data/derived/csf';out.mkdir(parents=True,exist_ok=True)
    a=run_csf(root);b=run_csf(root,method='Radau');diff=abs(np.array(a['state_values'])-np.array(b['state_values']))
    if not np.allclose(a['state_values'],b['state_values'],rtol=2e-7,atol=2e-8):raise RuntimeError('CSF independent solver disagreement')
    a['validation']['independent_solver_max_abs_difference']=diff.max(axis=1).tolist();runs=[('baseline',a)]
    native=root/'data/derived/physiology/native_baseline_v2/native_multisystem.csv'
    if native.exists():
        port=PressureInput.from_native_map(native);runs.append(('native_map_driven',run_csf(root,pressure_input=port,duration_s=float(port.time_s[-1]),samples=601)))
    ramp=PressureInput.from_samples([0,10,20,600],[100,100,90,90],unit='mmHg',provenance=dict(kind='source_Figure_7_pressure_ramp',source='Figure 7 arterial boundary; basal Table 1 parameters'))
    runs.append(('hypotension',run_csf(root,pressure_input=ramp)))
    models=[]
    for name,run in runs:
        path=out/(name+'.json');path.write_text(json.dumps(run,separators=(',',':'),allow_nan=False))
        models.append(dict(id='ursino_lodi_1997_'+name,execution_available=True,source_kind=run['source_kind'],trajectory_path=str(path.relative_to(root)),
            source=run['source'],duration_s=run['time_s'][-1],pressure_input=run['pressure_input'],parameters=run['parameters'],
            channels=[{k:v for k,v in c.items() if k!='values'} for c in run['channels']],validation=run['validation'],limitations=LIMITATIONS))
    (out/'index.json').write_text(json.dumps(dict(schema_version=1,models=models),separators=(',',':'),allow_nan=False));return dict(execution_available=True,models=models)
