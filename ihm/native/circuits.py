"""Import source-native circuit/compartment graphs without duplicating shared views."""
from pathlib import Path
from collections import Counter
import hashlib
import math
import xml.etree.ElementTree as ET

MMHG=133.322387415
# Explicit dimensional conversions used by native electrical/fluid/thermal circuits.
UNITS={
    'Pa':(1.,'Pa'),'mmHg':(MMHG,'Pa'),'cmH2O':(98.0665,'Pa'),'psi':(6894.757293168,'Pa'),
    'm^3':(1.,'m^3'),'L':(.001,'m^3'),'mL':(1e-6,'m^3'),
    'm^3/s':(1.,'m^3/s'),'L/s':(.001,'m^3/s'),'mL/s':(1e-6,'m^3/s'),
    'L/min':(.001/60,'m^3/s'),'mL/min':(1e-6/60,'m^3/s'),
    'Pa s/m^3':(1.,'Pa s/m^3'),'mmHg s/mL':(MMHG/1e-6,'Pa s/m^3'),
    'mmHg min/L':(MMHG*60/.001,'Pa s/m^3'),'mmHg min/mL':(MMHG*60/1e-6,'Pa s/m^3'),
    'cmH2O s/L':(98.0665/.001,'Pa s/m^3'),
    'm^3/Pa':(1.,'m^3/Pa'),'mL/mmHg':(1e-6/MMHG,'m^3/Pa'),'L/cmH2O':(.001/98.0665,'m^3/Pa'),
    'mmHg s^2/mL':(MMHG/1e-6,'Pa s^2/m^3'),'cmH2O s^2/L':(98.0665/.001,'Pa s^2/m^3'),
    'K':(1.,'K'),'degC':(1.,'K'),'J':(1.,'J'),'kcal':(4184.,'J'),
    'W':(1.,'W'),'kcal/hr':(4184./3600,'W'),'J/K':(1.,'J/K'),'kcal/degC':(4184.,'J/K'),
    'K/W':(1.,'K/W'),'degC s/kcal':(1./4184,'K/W'),
    'V':(1.,'V'),'mV':(.001,'V'),'A':(1.,'A'),'mA':(.001,'A'),'C':(1.,'C'),
    'F':(1.,'F'),'uF':(1e-6,'F'),'ohm':(1.,'ohm'),'Ohm':(1.,'ohm'),'H':(1.,'H'),
    's':(1.,'s'),'min':(60.,'s'),'hr':(3600.,'s'),'':(1.,'1')}
LIMITATIONS=[
    'Circuit coefficients and state values are native-engine initialized/stabilized, not independently calibrated human parameters.',
    'Circuit views share globally unique family/name nodes and paths; circuit memberships must not be summed as extra physical compartments.',
    'Some saved circuits are inactive equipment or alternative topologies. No activation is inferred from their presence or their stored flow values.',
    'Snapshot incidence balances include compliance/capacitance branch flux. Transport-only residual may be storage, not leakage; no dVolume/dt claim follows from one snapshot.',
    'Reference nodes and infinite reservoirs are boundary conditions, not finite body compartments. Reference heat/potential values can be algebraic, not physical tissue temperature.',
    'Compartment children and circuit-node mappings are native semantic relations, not geometric registrations. Parent and child volumes must not be added.',
    'Native System state properties inventory does not enumerate unexported source equations or identify every biochemical parameter; empty equipment states remain explicit.',
    'Anatomical labels are native compartment names; connections do not independently establish anatomical geometry or clinical validity.'
]


def _tag(e): return e.tag.rsplit('}',1)[-1]
def _text(e,name): return next(((x.text or '').strip() for x in e if _tag(x)==name),'')
def _texts(e,name): return [(x.text or '').strip() for x in e if _tag(x)==name]


def _scalar(e,strict=False):
    raw=e.get('value'); unit=e.get('unit',''); name=_tag(e)
    value=float(raw)
    if math.isnan(value): status='invalid_nan'
    elif math.isinf(value): status='infinite_boundary' if name in ('Volume','NextVolume','VolumeBaseline') and value>0 else 'nonfinite_native'
    else: status='finite'
    if strict and status not in ('finite','infinite_boundary'):raise ValueError(f'invalid nonfinite circuit quantity {name}: {raw}')
    if strict and unit not in UNITS: raise ValueError(f'unknown native circuit unit {unit!r} for {name}')
    conversion=UNITS.get(unit); si_value=None
    if conversion and math.isfinite(value):
        si_value=value*conversion[0]
        if unit=='degC' and name in ('Temperature','NextTemperature'):si_value+=273.15
        if not math.isfinite(si_value):raise ValueError('SI conversion overflow')
    return dict(value=value if math.isfinite(value) else None,raw_value=raw,unit=unit or 'dimensionless',
        si_value=si_value,si_unit=conversion[1] if conversion else None,
        status=status if conversion else 'unconverted_native_unit',read_only=e.get('readOnly'))


def _properties(e,strict=False):
    result={}
    def visit(parent,prefix=''):
        counts=Counter(_tag(c) for c in parent); seen=Counter()
        for c in parent:
            tag=_tag(c); seen[tag]+=1
            key=prefix+tag+(f'[{seen[tag]}]' if counts[tag]>1 else '')
            if c.get('value') is not None:result[key]=_scalar(c,strict)
            elif len(c):visit(c,key+'/')
    visit(e)
    return result


def _net_audit(circuit,nodes,paths):
    flux={'fluid':'Flow','thermal':'HeatTransferRate','electrical':'Current'}[circuit['family']]
    storage={'fluid':'Compliance','thermal':'Capacitance','electrical':'Capacitance'}[circuit['family']]
    net={ident:0. for ident in circuit['nodes']}; transport=dict(net); missing=set(); storage_nodes=set()
    for ident in circuit['paths']:
        p=paths[ident]; q=p['properties'].get(flux,{}).get('si_value')
        is_storage=any(name in p['properties'] for name in (storage,'Next'+storage,storage+'Baseline'))
        if is_storage:storage_nodes.update((p['source'],p['target']))
        if q is None:missing.update((p['source'],p['target']));continue
        for node,sign in ((p['source'],-1),(p['target'],1)):
            net[node]+=sign*q
            if not is_storage:transport[node]+=sign*q
    reports=[]
    for ident in circuit['nodes']:
        reports.append(dict(node=ident,net_inflow_si=net[ident] if ident not in missing else None,
            transport_net_inflow_si=transport[ident] if ident not in missing else None,
            storage_branch_net_inflow_si=net[ident]-transport[ident] if ident not in missing else None,
            touches_storage=ident in storage_nodes,is_reference=ident in circuit['reference_nodes'],complete=ident not in missing))
    interior=[r['net_inflow_si'] for r in reports if r['complete'] and not r['is_reference']]
    return dict(method='signed incidence on this circuit view only; positive is net inflow',
        flux_unit={'fluid':'m^3/s','thermal':'W','electrical':'A'}[circuit['family']],
        max_abs_net_inflow_si=max(map(abs,interior)) if interior else None,
        complete_interior_nodes=len(interior),incomplete_nodes=len(missing),nodes=reports,
        interpretation='algebraic snapshot diagnostic, not time-integrated volume/energy conservation; inactive views may retain incompatible state')


def native_circuit_graph(state_path):
    path=Path(state_path).resolve(); raw=path.read_bytes(); root=ET.fromstring(raw)
    manager=next((e for e in root if _tag(e)=='CircuitManager'),None)
    if manager is None:raise ValueError('native state has no CircuitManager')
    nodes={};paths={};circuits=[]
    for e in manager:
        tag=_tag(e); family=next((f for f in ('fluid','thermal','electrical') if tag.lower().startswith(f)),None)
        if not family:raise ValueError('unsupported native circuit element '+tag)
        name=_text(e,'Name');ident=family+':'+name
        if not name:raise ValueError('native circuit object lacks Name')
        properties=_properties(e,strict=True)
        expected={'Pressure':'Pa','PressureSource':'Pa','Volume':'m^3','Flow':'m^3/s','FlowSource':'m^3/s',
            'Compliance':'m^3/Pa','Inertance':'Pa s^2/m^3','Temperature':'K','TemperatureSource':'K',
            'Heat':'J','HeatSource':'W','HeatTransferRate':'W','Voltage':'V','VoltageSource':'V',
            'Current':'A','CurrentSource':'A','Charge':'C','Inductance':'H',
            'Resistance':{'fluid':'Pa s/m^3','thermal':'K/W','electrical':'ohm'}[family],
            'Capacitance':{'fluid':'m^3/Pa','thermal':'J/K','electrical':'F'}[family]}
        for key,quantity in properties.items():
            base=key.removeprefix('Next').removesuffix('Baseline')
            if base in expected and quantity['si_unit']!=expected[base]:
                raise ValueError('native circuit property has incompatible dimensions: '+key)

        if tag.endswith('Node'):
            record=dict(id=ident,family=family,name=name,properties=properties,circuits=[],compartments=[])
            if ident in nodes and nodes[ident]!=record:raise ValueError('conflicting duplicate node '+ident)
            nodes[ident]=record
        elif tag.endswith('Path'):
            mechanisms=[k for k in ('Resistance','Compliance','Capacitance','Inductance','Inertance','PressureSource','FlowSource','HeatSource','TemperatureSource','VoltageSource','CurrentSource') if any(n in properties for n in (k,'Next'+k,k+'Baseline'))]
            record=dict(id=ident,family=family,name=name,source=family+':'+_text(e,'SourceNode'),target=family+':'+_text(e,'TargetNode'),
                properties=properties,mechanisms=mechanisms,gate_states={_tag(x):(x.text or '').strip() for x in e if _tag(x) in ('Switch','NextSwitch','Valve','NextValve','PolarizedState','NextPolarizedState')},circuits=[])
            if ident in paths and paths[ident]!=record:raise ValueError('conflicting duplicate path '+ident)
            paths[ident]=record
        elif tag.endswith('Circuit'):
            circuits.append(dict(id=ident,family=family,name=name,nodes=list(dict.fromkeys(family+':'+v for v in _texts(e,'Node'))),
                paths=list(dict.fromkeys(family+':'+v for v in _texts(e,'Path'))),reference_nodes=[family+':'+v for v in _texts(e,'ReferenceNode')]))
        else:raise ValueError('unsupported native circuit structure '+tag)
    for p in paths.values():
        if p['source'] not in nodes or p['target'] not in nodes:raise ValueError('path has unknown endpoint '+p['id'])
    for c in circuits:
        for key,objects in [('nodes',nodes),('paths',paths)]:
            for ident in c[key]:
                if ident not in objects:raise ValueError('circuit references unknown '+ident)
                objects[ident]['circuits'].append(c['id'])
        if any(ident not in c['nodes'] for ident in c['reference_nodes']):raise ValueError('reference outside circuit')
        for ident in c['paths']:
            if paths[ident]['source'] not in c['nodes'] or paths[ident]['target'] not in c['nodes']:raise ValueError('path endpoint outside circuit')
        c['balance_audit']=_net_audit(c,nodes,paths)
    compartments=[];links=[];graphs=[];unmapped=[]
    compmanager=next((e for e in root if _tag(e)=='CompartmentManager'),[])
    for e in compmanager:
        tag=_tag(e); name=_text(e,'Name'); kind=tag.replace('Compartment','').replace('Link','').replace('Graph','').lower()
        family='thermal' if kind=='thermal' else ('electrical' if kind=='electrical' else 'fluid')
        if tag.endswith('Compartment'):
            refs=[family+':'+v for v in _texts(e,'Node')]
            children=[kind+':'+v for v in _texts(e,'Child')]
            ident=kind+':'+name
            record=dict(id=ident,kind=kind,name=name,native_anatomical_label=name,nodes=refs,children=children,
                properties=_properties(e),label_evidence='literal native CompartmentManager name',substances=_texts(e,'Substance'),
                substance_quantities=[dict(substance=_text(q,'Substance'),properties=_properties(q)) for q in e if _tag(q)=='SubstanceQuantity'])
            compartments.append(record)
            for node in refs:
                if node in nodes:nodes[node]['compartments'].append(ident)
                else:unmapped.append(dict(compartment=ident,node=node))
        elif tag.endswith('Link'):
            native_path=_text(e,'Path'); path_id=family+':'+native_path if native_path else None
            links.append(dict(id=kind+':'+name,kind=kind,name=name,path=path_id,source=kind+':'+_text(e,'SourceCompartment'),
                target=kind+':'+_text(e,'TargetCompartment'),properties=_properties(e),path_resolved=path_id in paths if path_id else False))
        elif tag.endswith('Graph'):
            graphs.append(dict(id=kind+':'+name,kind=kind,name=name,compartments=[kind+':'+v for v in _texts(e,'Compartment')],links=[kind+':'+v for v in _texts(e,'Link')]))
    systems=[]
    for e in root:
        if _tag(e)=='System':
            typename=next((v for k,v in e.attrib.items() if k.endswith('type')),'System')
            systems.append(dict(type=typename,properties=_properties(e),native_state_field_count=len(e),
                scope='exported system state, not a complete source-equation inventory'))
    nonfinite=[dict(object=obj['id'],property=k,status=v['status']) for obj in [*nodes.values(),*paths.values()] for k,v in obj['properties'].items() if v['status']!='finite']
    summary=dict(unique_nodes=len(nodes),unique_paths=len(paths),circuit_views=len(circuits),
        nodes_by_family=dict(Counter(n['family'] for n in nodes.values())),paths_by_family=dict(Counter(p['family'] for p in paths.values())),
        compartments_by_kind=dict(Counter(c['kind'] for c in compartments)),links_by_kind=dict(Counter(c['kind'] for c in links)),
        mechanisms=dict(Counter(m for p in paths.values() for m in p['mechanisms'])),system_types=[s['type'] for s in systems],
        system_scalar_fields=sum(len(s['properties']) for s in systems),
        circuit_scalar_fields=sum(len(obj['properties']) for obj in [*nodes.values(),*paths.values()]),
        shared_nodes=sum(len(n['circuits'])>1 for n in nodes.values()),shared_paths=sum(len(p['circuits'])>1 for p in paths.values()),
        nonfinite_properties=nonfinite,unmapped_compartment_nodes=unmapped,
        electrical_circuits_present=any(c['family']=='electrical' for c in circuits))
    simulation=next((e for e in root if _tag(e)=='SimulationTime'),None)
    return dict(schema_version=1,source=dict(path=str(path),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw)),
        simulation_time=_scalar(simulation) if simulation is not None else None,
        parameter_status='native_engine_initialized_or_stabilized_not_independently_calibrated',
        nodes=list(nodes.values()),paths=list(paths.values()),circuits=circuits,compartments=compartments,
        compartment_links=links,compartment_graphs=graphs,systems=systems,
        configuration_properties=next((_properties(e) for e in root if _tag(e)=='Configuration'),{}),
        patient_properties=next((_properties(e) for e in root if _tag(e)=='Patient'),{}),
        active_substances=[dict(name=_text(e,'Name'),properties=_properties(e)) for e in root if _tag(e)=='ActiveSubstance'],
        summary=summary,limitations=LIMITATIONS)
