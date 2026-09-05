"""Execute hash-pinned official Physiome gonadotropin equations without rewriting math."""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import xml.etree.ElementTree as ET
import numpy as np
from scipy.integrate import solve_ivp

CELLML_SHA256='46a4cd42dda1299862ffba0086cbb910710dd9a55b1212a4b97d17341e9cb573'
CODE_SHA256='910a322c4ac848f4467915a86858d56c3e0bf8da62f6faaf40929d851370c195'
LIMITATIONS=[
    'Standalone source-model simulation, not a coupled or individually calibrated human reproductive system.',
    'Estradiol, progesterone and inhibin are prescribed time functions. This source is not an autonomous ovarian cycle generator and has no validated emergent cycle-period target.',
    'The official generated example spans days 0–10; extending its nonperiodic inputs does not generate repeated menstrual cycles.',
    'Fixed input-function shifts implement delay parameters; this is not a delay differential equation with endogenous hormonal history.',
    'Source narrative contains conflicting historical delay notes. Executed generated equations and latest metadata include fixed shifts; no unsupported correction is invented.',
    'Parameters are source-publication values and source initial conditions; no individual clinical fitting, measured holdout or cross-model patient matching is performed.',
    'Assay unit U for inhibin has no supplied conversion to mass or molarity. It is preserved rather than converted into a fabricated SI concentration.',
    'Solver agreement and synthesis/release/clearance mass accounting validate execution, not clinical outcomes or fertility prediction.'
]


def _legend(text):
    match=re.fullmatch(r'(.+) in component (.+) \((.+)\)',text)
    if not match:raise ValueError('unexpected official source legend: '+text)
    return dict(id=match[1],component=match[2],unit=match[3],legend=text)


def load_source(root):
    directory=Path(root)/'data/raw/reproductive/schlosser_selgrade_2000'
    cellml=directory/'schlosser_selgrade_2000.cellml';code=directory/'generated.py'
    for path,expected in ((cellml,CELLML_SHA256),(code,CODE_SHA256)):
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError('source hash differs from reviewed official model: '+path.name)
    # This exact downloaded module has been inspected: only math/numpy imports,
    # generated pure equations and definitions execute on import. Plot entrypoint
    # is guarded by __main__; no third-party main script is executed.
    spec=importlib.util.spec_from_file_location('ihm_schlosser_official_source',code)
    source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
    tree=ET.fromstring(cellml.read_bytes());ns={'c':'http://www.cellml.org/cellml/1.0#'}
    variables={(c.get('name'),v.get('name')):v.attrib for c in tree.findall('c:component',ns) for v in c.findall('c:variable',ns)}
    states,constants=source.initConsts();ls,la,lt,lc=source.createLegends()
    state_meta=[_legend(v) for v in ls];constant_meta=[_legend(v) for v in lc]
    for item,value in list(zip(state_meta,states))+list(zip(constant_meta,constants)):
        xml=variables[(item['component'],item['id'])]
        if xml['units']!=item['unit'] or float(xml['initial_value'])!=value:raise ValueError('codegen differs from pinned XML initial value/unit')
        item['initial_value' if item in state_meta else 'value']=float(value)
        item['parameter_status']='official CellML source value; not independently fitted'
    units={e.get('name'):[dict(u.attrib) for u in e] for e in tree.findall('c:units',ns)}
    if units['day']!=[{'multiplier':'86400.0','units':'second'}] or _legend(lt)['unit']!='day':raise ValueError('source time unit changed')
    literals=[]
    for component in tree.findall('c:component',ns):
        for e in component.iter():
            if e.tag.endswith('}cn'):
                literals.append(dict(component=component.get('name'),value=float(''.join(e.itertext()).strip()),unit=next((v for k,v in e.attrib.items() if k.endswith('units')),'dimensionless')))
    return source,dict(states_count=len(states),constants_count=len(constants),algebraic_count=source.sizeAlgebraic,time_unit='day',
        states=state_meta,parameters=constant_meta,algebraic=[_legend(v) for v in la],unit_definitions=units,math_literals=literals,
        audited_against_cellml=True,source_hashes={'cellml':CELLML_SHA256,'official_generated_python':CODE_SHA256})


def run_reproductive(root,*,duration_days=10.,samples=501,method='BDF',rtol=1e-9,atol=1e-10):
    if not np.isfinite(duration_days) or not 0<duration_days<=30 or not isinstance(samples,int) or not 2<=samples<=10001 or method not in ('BDF','DOP853') or not np.isfinite(rtol) or not np.isfinite(atol) or not 1e-12<=rtol<=1e-4 or not 1e-13<=atol<=1e-4:raise ValueError('invalid bounded source integration settings')
    source,audit=load_source(root);initial,constants=source.initConsts();t=np.linspace(0,float(duration_days),samples)
    solution=solve_ivp(lambda time,state:source.computeRates(time,state,constants),(0,float(duration_days)),initial,
        method=method,t_eval=t,rtol=rtol,atol=atol,max_step=.1)
    if not solution.success or solution.y.shape!=(4,samples) or not np.isfinite(solution.y).all() or np.min(solution.y)<-1e-7:raise RuntimeError('source integration failed or produced invalid state: '+solution.message)
    algebraic=source.computeAlgebraic(constants,solution.y,t)
    if not np.isfinite(algebraic).all():raise RuntimeError('source algebraic outputs are nonfinite')
    rates=np.array([source.computeRates(time,state,constants) for time,state in zip(t,solution.y.T)]).T
    lhs_lh=rates[0]+constants[8]*rates[1];rhs_lh=algebraic[10]-constants[8]*algebraic[0]
    lhs_fsh=rates[2]+constants[15]*rates[3];rhs_fsh=algebraic[11]-constants[15]*algebraic[2]
    prescribed={'E2','P4','Ih','E2_dE','P4_dP','Ih_dIh'}
    channel_si={'microg':(1e-9,'kg'),'microg_l':(1e-6,'kg/m^3'),
        'microg_day':(1e-9/86400,'kg/s'),'microg_l_day':(1e-6/86400,'kg/m^3/s'),
        'ng_L':(1e-9,'kg/m^3'),'nmol_L':(1e-6,'mol/m^3'),'U_L':(None,None)}
    channels=[]
    for kind,meta,values in [('state',audit['states'],solution.y),('algebraic',audit['algebraic'],algebraic)]:
        for item,v in zip(meta,values):
            channels.append(dict(**item,kind='prescribed_time_input' if item['id'] in prescribed else kind,
                si_scale=channel_si[item['unit']][0],si_unit=channel_si[item['unit']][1],values=v.tolist(),minimum=float(v.min()),maximum=float(v.max()),maximum_time_days=float(t[np.argmax(v)])))
    provenance=json.loads((Path(root)/'data/raw/reproductive/schlosser_selgrade_2000/provenance.json').read_text())
    return dict(schema_version=1,id='schlosser_selgrade_2000',source_kind='source_model_simulation',system='reproductive_endocrine',
        source=provenance,source_audit=audit,time_unit='day',time_days=t.tolist(),time_s=(t*86400).tolist(),
        state_values=solution.y.tolist(),channels=channels,parameters=audit['parameters'],
        solver=dict(method=method,rtol=rtol,atol=atol,max_step_days=.1,function_evaluations=solution.nfev,successful=True),
        validation=dict(state_nonnegative=bool(solution.y.min()>=0),source_default_window=duration_days==10,
            max_LH_release_cancellation_microg_per_day=float(max(abs(lhs_lh-rhs_lh))),
            max_FSH_release_cancellation_microg_per_day=float(max(abs(lhs_fsh-rhs_fsh))),
            cycle_period_target_days=None,cycle_period_status='not applicable: externally prescribed nonperiodic ovarian inputs'),limitations=LIMITATIONS)


def build_reproductive(root):
    root=Path(root).resolve();a=run_reproductive(root);b=run_reproductive(root,method='DOP853')
    x=np.array(a['state_values']);y=np.array(b['state_values']);difference=abs(x-y)
    if not np.allclose(x,y,rtol=2e-6,atol=2e-5):raise RuntimeError('independent solver agreement failed')
    a['validation']['independent_solver']='DOP853 versus BDF, identical source equations/initial conditions'
    a['validation']['max_solver_absolute_difference_by_state']=difference.max(axis=1).tolist()
    a['validation']['max_solver_normalized_difference']=float(np.max(difference/(1+abs(y))))
    out=root/'data/derived/reproductive';out.mkdir(parents=True,exist_ok=True)
    (out/'trajectory.json').write_text(json.dumps(a,separators=(',',':'),allow_nan=False))
    index=dict(schema_version=1,models=[dict(id=a['id'],system=a['system'],source_kind=a['source_kind'],source=a['source'],
        trajectory_path='data/derived/reproductive/trajectory.json',duration_days=10.,states=a['source_audit']['states'],
        parameters=a['parameters'],channels=[{k:v for k,v in c.items() if k!='values'} for c in a['channels']],
        validation=a['validation'],limitations=LIMITATIONS)])
    (out/'index.json').write_text(json.dumps(index,separators=(',',':'),allow_nan=False));return a
