"""Frozen native hydraulic descriptor dynamics; no illustrative coefficients."""
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.linalg import eigvals

SKIN_NODES=['Skin1','Skin2','SkinE1','SkinE2','SkinE3','SkinI','SkinL1','SkinL2','Lymph','Aorta1','VenaCava','Ground']
SKIN_BOUNDARIES=['Aorta1','VenaCava','Ground']
LIMITATIONS=[
    'Native coefficients, source pressures/flows and gate states are frozen at the saved operating point; this is not a replacement for nonlinear native advancement.',
    'Aorta, vena cava and ground are prescribed pressure boundaries; other-organ lymph inflows are frozen native snapshot fluxes.',
    'Native Closed gates conduct and Open gates block. Conducting valve reverse flow is reported as a frozen-gate validity violation, not automatically re-solved.',
    'Small-signal transfer uses zero initial perturbation. Ideal pressure constraints create infinite descriptor eigenvalues; only finite poles are physical dynamical modes of this approximation.',
    'An isolated intracellular compliance with a frozen flow source can have a zero pole and secular volume drift. Long-horizon physiological stability is not asserted.',
    'Capacitance/compliance paths represent storage, not transport to ground; volume changes are derived from storage flux at native volume-bearing nodes.',
    'Patient-level parameters remain native-engine initialized, not independently calibrated; source anatomical names are semantic supports without geometric registration.'
]


class FluidCircuit:
    @classmethod
    def from_native(cls,graph,node_names,fixed_pressure_names,*,circuit_name='FullCardiovascular'):
        self=cls();self.names=list(node_names);self.fixed_names=list(fixed_pressure_names)
        if len(set(self.names))!=len(self.names) or not self.names or not self.fixed_names or len(set(self.fixed_names))!=len(self.fixed_names) or not set(self.fixed_names)<=set(self.names):raise ValueError('distinct selected nodes and nonempty pressure boundaries required')
        self.source=graph['source'];self.circuit_name=circuit_name
        lookup={n['name']:n for n in graph['nodes'] if n['family']=='fluid'}
        if not set(self.names)<=set(lookup):raise ValueError('selected native node missing')
        self.native_nodes=[lookup[name] for name in self.names]
        self.fixed=np.array([self.names.index(n) for n in self.fixed_names],int)
        self.free=np.array([i for i in range(len(self.names)) if i not in self.fixed],int)
        if not len(self.free):raise ValueError('at least one dynamic/algebraic internal node required')
        self.pressures=np.array([self._quantity(n['properties'],'Pressure','Pa',required=True) for n in self.native_nodes])
        self.initial_pressures=self.pressures.copy();self.time_s=0.
        self.volumes={n['name']:self._quantity(n['properties'],'Volume','m^3') for n in self.native_nodes}
        self.volumes={n:v for n,v in self.volumes.items() if v is not None}
        self.initial_volumes=dict(self.volumes)
        self.native_paths=[];self.crossing_paths=[];self.external=np.zeros(len(self.names))
        selected=set(self.names)
        for path in graph['paths']:
            if path['family']!='fluid' or 'fluid:'+circuit_name not in path.get('circuits',[]):continue
            a,b=path['source'].split(':',1)[1],path['target'].split(':',1)[1]
            if a in selected and b in selected:self.native_paths.append(path)
            elif (a in selected)!=(b in selected):
                inside=a if a in selected else b
                q=self._quantity(path['properties'],'Flow','m^3/s',required=inside not in self.fixed_names)
                treatment='pressure_boundary_absorbs_crossing' if inside in self.fixed_names else 'fixed_snapshot_inflow'
                flow=(q if b in selected else -q) if q is not None else None
                self.crossing_paths.append(dict(id=path['id'],name=path['name'],source=path['source'],target=path['target'],inside_node=inside,
                    inflow_m3_s=flow,treatment=treatment,source_properties=path['properties']))
                if treatment=='fixed_snapshot_inflow':self.external[self.names.index(inside)]+=flow
        self._compile()
        return self

    @classmethod
    def from_dict(cls,data):
        if data.get('schema_version')!=1 or data.get('kind')!='native_frozen_hydraulic_descriptor':raise ValueError('unsupported circuit serialization')
        crossings=[dict(id=p['id'],name=p['name'],family='fluid',source=p['source'],target=p['target'],properties=p['source_properties'],
            gate_states={},circuits=['fluid:'+data['circuit_name']]) for p in data['crossing_paths']]
        graph=dict(source=data['source'],nodes=data['native_nodes'],paths=data['native_paths']+crossings)
        model=cls.from_native(graph,data['node_names'],data['fixed_pressure_names'],circuit_name=data['circuit_name'])
        state=np.asarray(data['pressures_pa'],float);initial=np.asarray(data['initial_pressures_pa'],float)
        flows=np.asarray(data['ideal_flows_m3_s'],float)
        if state.shape!=model.pressures.shape or initial.shape!=state.shape or flows.shape!=model.ideal_flows.shape or not all(np.isfinite(v).all() for v in (state,initial,flows)) or not np.isfinite(data['time_s']) or data['time_s']<0:raise ValueError('invalid serialized circuit state')
        for key in ('volumes_m3','initial_volumes_m3'):
            if set(data[key])!=set(model.volumes) or not all(np.isfinite(v) for v in data[key].values()):raise ValueError('invalid serialized circuit volumes')
        model.pressures=state;model.initial_pressures=initial;model.ideal_flows=flows;model.time_s=float(data['time_s'])
        model.volumes=dict(data['volumes_m3']);model.initial_volumes=dict(data['initial_volumes_m3'])
        return model

    @staticmethod
    def _quantity(props,key,unit,required=False):
        if key not in props:
            if required:raise ValueError('native property required: '+key)
            return None
        item=props[key];v=item.get('si_value')
        if item.get('si_unit')!=unit or v is None or not np.isfinite(v):raise ValueError('invalid native SI property: '+key)
        return float(v)

    def _compile(self):
        n=len(self.names);self.M=np.zeros((n,n));self.K=np.zeros((n,n));self.b=self.external.copy();self.branches=[];ideal=[]
        for p in self.native_paths:
            props=p['properties'];a=self.names.index(p['source'].split(':',1)[1]);b=self.names.index(p['target'].split(':',1)[1])
            inc=np.zeros(n);inc[a]=1;inc[b]-=1
            gates=p.get('gate_states',{});blocked=any(gates.get(k)=='Open' for k in ('Valve','Switch'))
            if any(gates.get(k) not in (None,'Open','Closed') for k in ('Valve','Switch')):raise ValueError('unknown native gate state')
            branch=dict(name=p['name'],incidence=inc,source_index=a,target_index=b,native=p,blocked=blocked)
            if blocked:branch['kind']='blocked';self.branches.append(branch);continue
            r=self._quantity(props,'Resistance','Pa s/m^3');c=self._quantity(props,'Compliance','m^3/Pa')
            q=self._quantity(props,'FlowSource','m^3/s');s=self._quantity(props,'PressureSource','Pa')
            if any(k in props for k in ('Inertance','Inductance')):raise ValueError('inertance requires a separate branch state and is not supported by this compiler')
            if sum(v is not None for v in (r,c,q))>1 or (c is not None and s is not None) or (q is not None and s is not None):raise ValueError('unsupported compound native hydraulic path')
            if r is not None:
                if r<=0:raise ValueError('native resistance must be positive')
                branch.update(kind='resistance',resistance=r,pressure_source=s or 0.)
                self.K+=np.outer(inc,inc)/r;self.b-=inc*(s or 0.)/r
            elif c is not None:
                if c<=0:raise ValueError('native compliance must be positive')
                branch.update(kind='compliance',compliance=c);self.M+=c*np.outer(inc,inc)
            elif q is not None:
                branch.update(kind='flow_source',flow_source=q);self.b-=inc*q
            else:
                branch.update(kind='ideal',pressure_source=s or 0.,ideal_index=len(ideal));ideal.append(branch)
            self.branches.append(branch)
        self.ideal=ideal;self.A=np.column_stack([x['incidence'] for x in ideal]) if ideal else np.empty((n,0))
        # Fixed-to-fixed ideal constraints cannot determine their reaction flow.
        if any(not np.any(x['incidence'][self.free]) for x in ideal):raise ValueError('ideal branch between two prescribed pressure nodes is redundant')
        f=self.free;self.pressure_scale=1e4;self.flow_scale=1e-6
        nf=len(f);ni=len(ideal);self.descriptor_M=np.zeros((nf+ni,nf+ni));self.descriptor_K=self.descriptor_M.copy()
        self.descriptor_M[:nf,:nf]=self.M[np.ix_(f,f)]*self.pressure_scale/self.flow_scale
        self.descriptor_K[:nf,:nf]=self.K[np.ix_(f,f)]*self.pressure_scale/self.flow_scale
        self.descriptor_K[:nf,nf:]=self.A[f];self.descriptor_K[nf:,:nf]=self.A[f].T
        self.ideal_flows=np.zeros(ni)
        # A nonsingular backward-Euler pencil establishes an index-1 solvable network.
        if np.linalg.matrix_rank(self.descriptor_K+self.descriptor_M)<nf+ni:raise ValueError('native selection has redundant constraints or unconstrained pressure modes')

    def _rhs(self,boundaries,boundary_rates=None):
        f=self.free;h=self.fixed;nf=len(f)
        result=np.zeros(nf+len(self.ideal));rate=np.zeros(len(h)) if boundary_rates is None else boundary_rates
        result[:nf]=(self.b[f]-self.K[np.ix_(f,h)]@boundaries-self.M[np.ix_(f,h)]@rate)/self.flow_scale
        result[nf:]=(-np.array([i['pressure_source'] for i in self.ideal])-self.A[h].T@boundaries)/self.pressure_scale
        return result

    def step(self,dt,boundary_pressures_pa=None):
        if not np.isfinite(dt) or dt<=0:raise ValueError('positive finite timestep required')
        boundary_pressures_pa={} if boundary_pressures_pa is None else boundary_pressures_pa
        if not set(boundary_pressures_pa)<=set(self.fixed_names) or not all(np.isfinite(v) for v in boundary_pressures_pa.values()):raise ValueError('only declared finite pressure boundaries may be changed')
        previous=self.pressures.copy();updated=previous.copy()
        for name,value in boundary_pressures_pa.items():updated[self.names.index(name)]=value
        oldstate=np.r_[previous[self.free]/self.pressure_scale,self.ideal_flows/self.flow_scale]
        matrix=self.descriptor_K+self.descriptor_M/dt
        rhs=self._rhs(updated[self.fixed],(updated[self.fixed]-previous[self.fixed])/dt)+self.descriptor_M@oldstate/dt
        z=np.linalg.solve(matrix,rhs)
        updated[self.free]=z[:len(self.free)]*self.pressure_scale
        self.ideal_flows=z[len(self.free):]*self.flow_scale
        rates=(updated-previous)/dt;flows={};net_out=np.zeros(len(self.names));storage=np.zeros(len(self.names));violations=[]
        for branch in self.branches:
            inc=branch['incidence'];kind=branch['kind']
            if kind=='blocked':q=0.
            elif kind=='resistance':q=(inc@updated+branch['pressure_source'])/branch['resistance']
            elif kind=='compliance':q=branch['compliance']*(inc@rates);storage+=inc*q
            elif kind=='flow_source':q=branch['flow_source']
            else:q=self.ideal_flows[branch['ideal_index']]
            flows[branch['name']]=float(q);net_out+=inc*q
            gates=branch['native'].get('gate_states',{})
            if gates.get('Valve')=='Closed' and q < -1e-15:violations.append(dict(path=branch['name'],flow_m3_s=float(q),reason='reverse flow through frozen conducting valve'))
        residual=net_out-self.external
        deltas={name:float(storage[self.names.index(name)]*dt) for name in self.volumes}
        for name,dv in deltas.items():self.volumes[name]+=dv
        if not np.isfinite(updated).all() or not all(np.isfinite(v) for v in flows.values()):raise FloatingPointError('nonfinite circuit state')
        self.pressures=updated;self.time_s+=dt
        return dict(time_s=self.time_s,pressures_pa=dict(zip(self.names,updated.tolist())),flows_m3_s=flows,volume_deltas_m3=deltas,
            volumes_m3=dict(self.volumes),balance=dict(max_abs_free_node_residual_m3_s=float(max(abs(residual[self.free]))),
                free_node_residual_m3_s={self.names[i]:float(residual[i]) for i in self.free},
                required_pressure_boundary_inflow_m3_s={self.names[i]:float(residual[i]) for i in self.fixed},
                storage_rate_m3_s={self.names[i]:float(storage[i]) for i in self.free}),
            gate_violations=violations,negative_volume_nodes=[name for name,v in self.volumes.items() if v<0],
            scaled_linear_condition=float(np.linalg.cond(matrix)))

    def response(self,s,boundary_name):
        if boundary_name not in self.fixed_names:raise ValueError('unknown pressure boundary')
        points=np.atleast_1d(np.asarray(s,complex))
        if points.ndim!=1 or not np.isfinite(points).all():raise ValueError('finite Laplace points required')
        j=self.fixed_names.index(boundary_name);u=np.zeros(len(self.fixed));u[j]=1
        pressures={name:[] for name in self.names};flows={b['name']:[] for b in self.branches}
        for point in points:
            nf=len(self.free);rhs=np.zeros(nf+len(self.ideal),complex)
            rhs[:nf]=-(self.K[np.ix_(self.free,self.fixed)]+point*self.M[np.ix_(self.free,self.fixed)])@u/self.flow_scale
            rhs[nf:]=-self.A[self.fixed].T@u/self.pressure_scale
            # At a pole the resolvent is undefined; callers must avoid that s.
            z=np.linalg.solve(self.descriptor_K+point*self.descriptor_M,rhs)
            p=np.zeros(len(self.names),complex);p[self.free]=z[:nf]*self.pressure_scale;p[self.fixed]=u
            iq=z[nf:]*self.flow_scale
            for name,value in zip(self.names,p):pressures[name].append(value)
            for b in self.branches:
                if b['kind']=='resistance':q=(b['incidence']@p)/b['resistance']
                elif b['kind']=='compliance':q=point*b['compliance']*(b['incidence']@p)
                elif b['kind']=='ideal':q=iq[b['ideal_index']]
                else:q=0j
                flows[b['name']].append(q)
        return dict(pressures_pa_per_pa={k:np.array(v) for k,v in pressures.items()},flows_m3_s_per_pa={k:np.array(v) for k,v in flows.items()})

    def finite_poles_per_s(self):
        poles=eigvals(-self.descriptor_K,self.descriptor_M)
        return poles[np.isfinite(poles)]

    def to_dict(self):
        poles=self.finite_poles_per_s();return dict(schema_version=1,kind='native_frozen_hydraulic_descriptor',source=self.source,
            circuit_name=self.circuit_name,node_names=self.names,fixed_pressure_names=self.fixed_names,time_s=self.time_s,
            pressures_pa=self.pressures.tolist(),initial_pressures_pa=self.initial_pressures.tolist(),ideal_flows_m3_s=self.ideal_flows.tolist(),volumes_m3=dict(self.volumes),
            initial_volumes_m3=self.initial_volumes,native_nodes=self.native_nodes,native_paths=self.native_paths,crossing_paths=self.crossing_paths,
            descriptor=dict(M=self.descriptor_M.tolist(),K=self.descriptor_K.tolist(),free_node_names=[self.names[i] for i in self.free],
                ideal_flow_path_names=[i['name'] for i in self.ideal],pressure_scale_pa=self.pressure_scale,flow_scale_m3_s=self.flow_scale,
                equation='M zdot + K z = b; z contains pressure/pressure_scale and ideal_flow/flow_scale'),
            finite_poles_per_s=dict(real=poles.real.tolist(),imag=poles.imag.tolist()),
            parameter_status='native engine operating-point values, frozen, not independently calibrated',limitations=LIMITATIONS)


def build_skin_circuit(root):
    root=Path(root).resolve();graph_path=root/'data/derived/native-circuits/graph.json';graph_raw=graph_path.read_bytes();graph=json.loads(graph_raw)
    source_path=Path(graph['source']['path'])
    if hashlib.sha256(source_path.read_bytes()).hexdigest()!=graph['source']['sha256']:raise ValueError('native source state hash differs from circuit graph provenance')
    model=FluidCircuit.from_native(graph,SKIN_NODES,SKIN_BOUNDARIES)
    initial=model.to_dict();baseline=model.step(.02)
    nativeflows={p['name']:p['properties'].get('Flow',{}).get('si_value') for p in model.native_paths}
    comparison={name:dict(native_snapshot_m3_s=q,frozen_next_step_m3_s=baseline['flows_m3_s'][name],difference_m3_s=baseline['flows_m3_s'][name]-q) for name,q in nativeflows.items() if q is not None}
    # Nonzero sigma avoids the genuine zero pole of frozen intracellular storage.
    frequencies=np.r_[0.,np.logspace(-6,0,121)];sigma=1e-8
    response=model.response(sigma+2j*np.pi*frequencies,'Aorta1')
    encode=lambda mapping:{k:dict(real=v.real.tolist(),imag=v.imag.tolist()) for k,v in mapping.items()}
    perturb=FluidCircuit.from_native(graph,SKIN_NODES,SKIN_BOUNDARIES)
    baseline_pressure=perturb.pressures[perturb.names.index('Aorta1')]
    samples=[]
    for i in range(300):
        step=perturb.step(1.,{'Aorta1':baseline_pressure+100.})
        if i%5==0 or i==299:samples.append(step)
    result=dict(schema_version=1,id='native_skin_lymph_circuit',graph_source=dict(path=str(graph_path.relative_to(root)),sha256=hashlib.sha256(graph_raw).hexdigest()),model=initial,baseline_step=baseline,
        native_snapshot_comparison=comparison,
        laplace=dict(boundary='Aorta1',sigma_per_s=sigma,frequency_hz=frequencies.tolist(),
            pressures_pa_per_pa=encode(response['pressures_pa_per_pa']),flows_m3_s_per_pa=encode(response['flows_m3_s_per_pa'])),
        perturbation=dict(boundary='Aorta1',delta_pressure_pa=100.,duration_s=300.,dt_s=1.,samples=samples,
            label='hypothetical +100 Pa boundary step of frozen native coefficient model; not native patient intervention validation'),limitations=LIMITATIONS)
    out=root/'data/derived/coupling';out.mkdir(parents=True,exist_ok=True)
    (out/'native-skin-circuit.json').write_text(json.dumps(result,separators=(',',':'),allow_nan=False))
    return result
