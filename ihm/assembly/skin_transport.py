"""Isolated, body-attached regional fluid/albumin transport experiment.

Native skin laws/operating-point coefficients are transferred with an explicit
geometric scale assumption. This module never writes native-engine state.
"""
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def albumin_flux_kg_s(q, blood_c, tissue_c, ps, sigma=.954):
    """Native positive-only Peclet law, with analytic zero-flow limit.

    q, ps: m3/s; concentrations: kg/m3; output: kg/s.
    """
    pe = q * (1 - sigma) / ps
    if abs(pe) < 1e-6:
        a = 1 + pe / 2 + pe * pe / 12
        b = 1 - pe / 2 + pe * pe / 12
    elif pe > 0:
        a = pe / -math.expm1(-pe)
        b = a * math.exp(-pe)
    else:
        b = -pe / -math.expm1(pe)
        a = b * math.exp(pe)
    return max(0., ps * (a * blood_c - b * tissue_c))


def oncotic_pa(albumin_kg_m3):
    protein_g_dl = 1.6 * albumin_kg_m3 / 10
    return 133.322387415 * (2.1 * protein_g_dl + .16 * protein_g_dl**2 + .009 * protein_g_dl**3)


class RegionalSkinTransport:
    @classmethod
    def from_sources(cls, root, unit_index=0):
        root = Path(root).resolve()
        circuit = root / 'data/derived/coupling/native-skin-circuit.json'
        micro = root / 'data/derived/canonical/microvascular.json'
        native = json.loads(circuit.read_text())['model']
        units = json.loads(micro.read_text())['units']
        unit = units[unit_index]
        if 'forearm skin' not in unit['territory']: raise ValueError('forearm skin material attachment required')
        nodes = {n['name']: n['properties'] for n in native['native_nodes']}
        paths = {p['name']: p['properties'] for p in native['native_paths']}
        expected_units={'Pressure':'Pa','Volume':'m^3','Resistance':'Pa s/m^3','Compliance':'m^3/Pa','PressureSource':'Pa'}
        def quantity(props,k):
            item=props[k]
            if item['si_unit']!=expected_units[k] or not math.isfinite(item['si_value']): raise ValueError('invalid native SI quantity')
            return item['si_value']
        def nv(n, k): return quantity(nodes[n],k)
        def pv(n, k): return quantity(paths[n],k)
        if len(unit['edges'])!=len(unit['radius_m']) or any(r<=0 for r in unit['radius_m']): raise ValueError('invalid attached lumen geometry')
        lumen = sum(math.pi * r*r * math.dist(unit['positions_m'][a], unit['positions_m'][b])
                    for (a,b),r in zip(unit['edges'], unit['radius_m']))
        scale = lumen / nv('Skin1', 'Volume')
        state_path = Path(native['source']['path'])
        if hashlib.sha256(state_path.read_bytes()).hexdigest() != native['source']['sha256']:
            raise ValueError('native source hash mismatch')
        tree = ET.parse(state_path).getroot()
        def local(e): return e.tag.split('}')[-1]
        def child(e, name): return next(c for c in e if local(c) == name)
        compartments = {}
        for e in tree.iter():
            if local(e) in ('LiquidCompartment', 'TissueCompartment'):
                compartments[child(e, 'Name').text] = e
        def concentration(name):
            for e in compartments[name]:
                if local(e) == 'SubstanceQuantity' and child(e, 'Substance').text == 'Albumin':
                    c = child(e, 'Concentration')
                    if c.attrib['unit'] != 'ug/mL': raise ValueError('unexpected albumin unit')
                    return float(c.attrib['value']) * .001
            raise ValueError('native albumin missing')
        tissue_mass = child(compartments['SkinTissue'], 'TotalMass')
        if tissue_mass.attrib['unit'] != 'kg': raise ValueError('unexpected tissue mass unit')
        p = dict(scale=scale, external_reservoir_volume_multiplier=1000., arterial_pa=nv('Aorta1','Pressure'), venous_pa=nv('VenaCava','Pressure'),
                 pressure_reference_pa={k:nv(n,'Pressure') for k,n in [('blood','Skin1'),('interstitium','SkinE3'),('lymph','Lymph')]},
                 compliance_m3_pa={k:pv(n,'Compliance')*scale for k,n in [('blood','Skin1ToGround'),('interstitium','SkinE3ToGround'),('lymph','LymphToGround')]},
                 resistance_pa_s_m3={k:pv(n,'Resistance')/scale for k,n in [('arterial_supply','Aorta1ToSkin1'),('venous_return','Skin1ToSkin2'),('filtration','SkinE1ToSkinE2'),('lymph_uptake','SkinL1ToSkinL2'),('lymph_return','LymphToVenaCava')]},
                 lymph_drive_pa=pv('SkinE3ToSkinL1','PressureSource'),
                 albumin_ps_m3_s=.03306*1e-6/60*float(tissue_mass.attrib['value'])*scale,
                 albumin_reflection=.954, fluid_oncotic_reflection=1., max_internal_dt_s=1.)
        volumes = {k:nv(n,'Volume')*scale for k,n in [('blood','Skin1'),('interstitium','SkinE3'),('lymph','Lymph')]}
        concentrations = dict(blood=concentration('SkinVasculature'), interstitium=concentration('SkinTissueExtracellular'), lymph=concentration('Lymph'))
        # Finite external reservoirs close the ledger; their pressures are prescribed.
        volumes.update(arterial_reservoir=sum(volumes.values())*1000, venous_reservoir=sum(volumes.values())*1000)
        concentrations.update(arterial_reservoir=concentrations['blood'], venous_reservoir=concentration('VenaCava'))
        self = cls()
        self.parameters=p; self.volumes=volumes; self.mass={k:v*concentrations[k] for k,v in volumes.items()}
        self.initial_volumes=dict(volumes); self.initial_mass=dict(self.mass)
        self.time_s=0.; self.integrated={k:0. for k in p['resistance_pa_s_m3']}; self.integrated_mass={k:0. for k in self.integrated}
        self.attachment=dict(unit_id=unit['id'], territory=unit['territory'], material_attachment=unit['material_attachment'],
                             boundary_association=unit['boundary_association'], geometric_lumen_volume_m3=lumen)
        files=[circuit,micro,state_path,root/'ihm/assembly/skin_transport.py',root/'scripts/verify_body_transport.py',root/'scripts/build_body_transport.py',root/'docs/research/REGIONAL_SKIN_TRANSPORT.md']
        source_dir=root/'data/raw/physiology/biogears/projects/biogears/libBiogears/src/engine'
        files += [source_dir/'Systems/Diffusion.cpp',source_dir/'Systems/Tissue.cpp',source_dir/'Controller/BioGears.cpp']
        self.provenance=[dict(path=str(f.relative_to(root)),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]
        return self

    def rates(self, venous_delta_pa=0., inflow_fraction=1., lymph_obstruction=0.):
        p=self.parameters; r=p['resistance_pa_s_m3']; c={k:self.mass[k]/v for k,v in self.volumes.items()}
        pressures={k:p['pressure_reference_pa'][k]+(self.volumes[k]-self.initial_volumes[k])/p['compliance_m3_pa'][k]
                   for k in p['pressure_reference_pa']}
        blood,tissue,lymph=(pressures[k] for k in ('blood','interstitium','lymph'))
        venous=p['venous_pa']+venous_delta_pa
        q=dict(arterial_supply=inflow_fraction*(p['arterial_pa']-blood)/r['arterial_supply'],
               venous_return=(blood-venous)/r['venous_return'],
               filtration=(blood-tissue-oncotic_pa(c['blood'])+oncotic_pa(c['interstitium']))/r['filtration'],
               lymph_uptake=(1-lymph_obstruction)*max(0.,tissue+p['lymph_drive_pa']-lymph)/r['lymph_uptake'],
               lymph_return=(1-lymph_obstruction)*max(0.,lymph-venous)/r['lymph_return'])
        ends=dict(arterial_supply=('arterial_reservoir','blood'),venous_return=('blood','venous_reservoir'),
                  filtration=('blood','interstitium'),lymph_uptake=('interstitium','lymph'),lymph_return=('lymph','venous_reservoir'))
        j={k:v*c[ends[k][0] if v>=0 else ends[k][1]] for k,v in q.items()}
        j['filtration']=albumin_flux_kg_s(q['filtration'],c['blood'],c['interstitium'],p['albumin_ps_m3_s'],p['albumin_reflection'])
        return q,j,ends,pressures

    def step(self, dt_s, *, venous_delta_pa=0., inflow_fraction=1., lymph_obstruction=0.):
        if not all(math.isfinite(x) for x in (dt_s,venous_delta_pa,inflow_fraction,lymph_obstruction)) or dt_s<=0 or not 0<=inflow_fraction<=1 or not 0<=lymph_obstruction<=1:
            raise ValueError('finite positive dt and intervention fractions in [0,1] required')
        remaining=dt_s
        while remaining > max(1e-12,dt_s*1e-14):
            q,j,ends,_=self.rates(venous_delta_pa,inflow_fraction,lymph_obstruction)
            h=min(remaining,self.parameters['max_internal_dt_s'])
            # Bound aggregate donor losses before simultaneous conservative transfers.
            for values,rates in ((self.volumes,q),(self.mass,j)):
                outgoing={k:0. for k in values}
                for k,rate in rates.items(): outgoing[ends[k][0] if rate>=0 else ends[k][1]] += abs(rate)
                for k,rate in outgoing.items():
                    if rate>0: h=min(h,.1*values[k]/rate)
            if not math.isfinite(h) or h<=1e-14: raise FloatingPointError('transport exhausted a reservoir or became stiff')
            for values,rates,integrated in ((self.volumes,q,self.integrated),(self.mass,j,self.integrated_mass)):
                for k,rate in rates.items():
                    a,b=ends[k]; moved=h*rate
                    values[a]-=moved; values[b]+=moved; integrated[k]+=moved
            remaining-=h
        self.time_s+=dt_s
        return self.to_dict()

    def to_dict(self):
        # Regional change versus integrated reservoir exchanges avoids cancellation
        # in the much larger external reservoir totals.
        local=('blood','interstitium','lymph')
        dv=sum(self.volumes[k]-self.initial_volumes[k] for k in local)
        dm=sum(self.mass[k]-self.initial_mass[k] for k in local)
        def net(z): return z['arterial_supply']-z['venous_return']-z['lymph_return']
        _,_,ends,_=self.rates()
        def node_residual(values,initial,integrated):
            expected={k:0. for k in values}
            for k,(a,b) in ends.items(): expected[a]-=integrated[k]; expected[b]+=integrated[k]
            return {k:values[k]-initial[k]-expected[k] for k in values}
        return deepcopy(dict(schema='ihm.regional_skin_transport.v1',native_blood_storage_connected=False,
            coupling_mode='isolated_regional_native_law_transfer',time_s=self.time_s,parameters=self.parameters,
            volumes_m3=self.volumes,albumin_kg=self.mass,initial_volumes_m3=self.initial_volumes,initial_albumin_kg=self.initial_mass,
            integrated_flows_m3=self.integrated,integrated_albumin_kg=self.integrated_mass,attachment=self.attachment,provenance=self.provenance,
            balance=dict(volume_residual_m3=dv-net(self.integrated),albumin_residual_kg=dm-net(self.integrated_mass),
                         node_volume_residual_m3=node_residual(self.volumes,self.initial_volumes,self.integrated),
                         node_albumin_residual_kg=node_residual(self.mass,self.initial_mass,self.integrated_mass),
                         full_reservoir_volume_residual_m3=math.fsum(self.volumes[k]-self.initial_volumes[k] for k in self.volumes),
                         full_reservoir_albumin_residual_kg=math.fsum(self.mass[k]-self.initial_mass[k] for k in self.mass))))

    @classmethod
    def from_dict(cls, data):
        if data.get('schema')!='ihm.regional_skin_transport.v1' or data.get('native_blood_storage_connected') is not False:
            raise ValueError('unsupported transport checkpoint/ownership')
        self=cls(); d=deepcopy(data)
        for attr,key in [('parameters','parameters'),('volumes','volumes_m3'),('mass','albumin_kg'),('initial_volumes','initial_volumes_m3'),
                         ('initial_mass','initial_albumin_kg'),('integrated','integrated_flows_m3'),('integrated_mass','integrated_albumin_kg'),
                         ('attachment','attachment'),('provenance','provenance'),('time_s','time_s')]: setattr(self,attr,d[key])
        if not all(math.isfinite(v) and v>0 for v in self.volumes.values()) or not all(math.isfinite(v) and v>=0 for v in self.mass.values()):
            raise ValueError('nonpositive/nonfinite checkpoint storage')
        expected={'blood','interstitium','lymph','arterial_reservoir','venous_reservoir'}
        if any(set(z)!=expected for z in (self.volumes,self.mass,self.initial_volumes,self.initial_mass)):
            raise ValueError('checkpoint compartment mismatch')
        if not all(math.isfinite(v) and v>0 for v in self.initial_volumes.values()) or not all(math.isfinite(v) and v>=0 for v in self.initial_mass.values()):
            raise ValueError('invalid initial checkpoint storage')
        def finite_tree(value):
            if isinstance(value,dict):return all(finite_tree(v) for v in value.values())
            return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)
        if not finite_tree(self.parameters):raise ValueError('nonfinite checkpoint parameter')
        edges={'arterial_supply','venous_return','filtration','lymph_uptake','lymph_return'}
        if any(set(z)!=edges or not finite_tree(z) for z in (self.integrated,self.integrated_mass)):
            raise ValueError('invalid integrated checkpoint ledger')
        if set(self.parameters['resistance_pa_s_m3'])!=edges or any(set(self.parameters[k])!={'blood','interstitium','lymph'} for k in ('pressure_reference_pa','compliance_m3_pa')):
            raise ValueError('checkpoint constitutive topology differs')
        if not math.isfinite(self.time_s) or self.time_s<0: raise ValueError('invalid checkpoint time')
        for key in ('compliance_m3_pa','resistance_pa_s_m3'):
            if not all(math.isfinite(v) and v>0 for v in self.parameters[key].values()): raise ValueError('invalid hydraulic coefficient')
        for key in ('albumin_ps_m3_s','max_internal_dt_s'):
            if not math.isfinite(self.parameters[key]) or self.parameters[key]<=0: raise ValueError('invalid transport coefficient')
        if not 0<=self.parameters['albumin_reflection']<=1: raise ValueError('invalid reflection coefficient')
        return self


def build_skin_transport(root):
    """Materialize one canonical forearm region and record actual experiments."""
    from scripts.verify_body_transport import verify
    root=Path(root)
    initial=RegionalSkinTransport.from_sources(root).to_dict()
    artifact=dict(schema='ihm.body_skin_transport.v1',seed=417,stochastic_steps=False,model=initial,
                  experiments=verify(root),parameter_status='native source law and snapshot coefficients; regional scale and lymph allocation are uncalibrated priors',
                  native_blood_storage_connected=False,ownership='native full-skin remains native; isolated regional replica has private boundary reservoirs',
                  limitations=['Synthetic microvascular lumen volume determines extensive scale; no measured forearm tissue volume or vessel registration.',
                               'Whole-body lymph storage scaled by skin lumen fraction is an explicit allocation prior.',
                               'Intracellular exchange, sweating, active lymph pumping changes, inflammation and blood rheology are omitted.',
                               'Native albumin flux is positive-only; oncotic fluid reflection 1 differs from albumin reflection 0.954 by native design.',
                               'Pressure boundaries are prescribed; external reservoirs track fluid and albumin but do not feed back to native engine.'])
    output=root/'data/derived/canonical/skin-transport.json'
    output.write_text(json.dumps(artifact,indent=2,allow_nan=False)+'\n')
    return artifact
