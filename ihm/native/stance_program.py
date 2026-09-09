"""Digest-bound reference scheduling for the accepted engineering weight transfer.

This module returns reference states, equilibrium excitations and feedback gains.
It never writes body coordinates or applies forces. Only the single accepted
left-transfer segment is currently loadable; segment storage allows later schemas.
"""
from copy import deepcopy
from pathlib import Path
from types import MappingProxyType
import hashlib,io,json,math
import numpy as np
SOURCE_BYTES=Path(__file__).read_bytes()
def _finite(value):
    if isinstance(value,(bool,np.bool_)) or not np.isscalar(value) or not math.isfinite(value):raise ValueError('Finite nonboolean program time required')
    return float(value)
def _digest(raw):return hashlib.sha256(raw).hexdigest()
def _arrays(raw):
    with np.load(io.BytesIO(raw),allow_pickle=False) as data:return {k:data[k].copy() for k in data.files}

class StanceReferenceProgram:
    @classmethod
    def from_bundle(cls,root,bundle_path='data/models/engineering_weight_transfer_v1'):
        root=Path(root).resolve();self=cls();self.source_bytes={}
        def read(relative,digest=None,retain=True):
            path=Path(relative)
            if path.is_absolute() or '..' in path.parts:raise ValueError('Workspace-relative program source required')
            path=root/path
            if any(p.is_symlink() for p in (path,*path.parents)):raise ValueError('Program source may not redirect through symlinks')
            raw=path.read_bytes()
            if digest is not None and _digest(raw)!=digest:raise ValueError('Stance program source digest differs: '+str(relative))
            if retain:self.source_bytes[str(relative)]=raw
            return raw
        bundle=Path(bundle_path);manifest_raw=read(str(bundle/'manifest.json'));manifest=json.loads(manifest_raw)
        if manifest.get('schema')!='ihm.engineering-weight-transfer-bundle.v1':raise ValueError('Unsupported stance program bundle')
        for path,digest in manifest['files'].items():read(path,digest,retain=False)
        program_path=str(bundle/'program.json');program_raw=read(program_path,manifest['files'][program_path]);program=json.loads(program_raw)
        if program.get('schema')!='ihm.engineering-weight-transfer-program.v1' or program['baseline_binding']!=manifest['baseline_binding']:raise ValueError('Program and manifest baseline bindings disagree')
        if program['runtime_model_substitution'] or program['external_root_forces'] or program['prescribed_motion']:raise ValueError('Only reference-based programs are supported')
        binding=program['baseline_binding'];base=Path(binding['bundle_path'])
        read(str(base/'manifest.json'),binding['manifest_sha256'])
        reg=json.loads(read(str(base/'registration.json'),binding['registration_sha256']))
        read(reg['model_path'],binding['native_model_sha256'])
        policy=_arrays(read(str(base/'linearization.npz'),binding['policy_sha256']))
        target=_arrays(read(program['target_path'],manifest['files'][program['target_path']]))
        goal=_arrays(read(program['target_policy_path'],manifest['files'][program['target_policy_path']]))
        if str(policy['model_sha256'].item())!=binding['native_model_sha256'] or reg['model_sha256']!=binding['native_model_sha256'] or str(target['source_model_sha256'].item())!=binding['native_model_sha256']:raise ValueError('Program baseline native model binding differs')
        if str(goal['model_sha256'].item())!=str(target['target_model_sha256'].item()):raise ValueError('Endpoint gain model binding differs')
        proof_path=str(bundle/'provenance/model_equivalence.json');proof=json.loads(read(proof_path,manifest['files'][proof_path]))
        if proof.get('same_except_default_activation') is not True or proof.get('runtime_model_substitution') is not False or proof['baseline_model_sha256']!=binding['native_model_sha256'] or proof['target_model_sha256']!=str(target['target_model_sha256'].item()):raise ValueError('Endpoint model equivalence proof does not match program')
        for key in ('state_names','muscle_names'):
            if not np.array_equal(policy[key],target[key]) or not np.array_equal(policy[key],goal[key]):raise ValueError('Program native state chart or muscle order differs')
        self.target_mass_kg=_finite(binding['target_mass_kg'])
        if self.target_mass_kg<=0 or any(float(d['target_mass_kg'])!=self.target_mass_kg for d in (policy,target,goal)) or reg['target_mass_kg']!=self.target_mass_kg:raise ValueError('Program patient mass differs')
        self.dt_s=_finite(program['sampling_interval_s'])
        if self.dt_s!=.01 or float(policy['dt_s'])!=.01 or float(goal['dt_s'])!=.01:raise ValueError('Accepted program requires10ms sampling')
        self.state_names=tuple(policy['state_names'].tolist());self.muscle_names=tuple(policy['muscle_names'].tolist())
        if len(set(self.state_names))!=len(self.state_names) or len(set(self.muscle_names))!=len(self.muscle_names):raise ValueError('Unique state/muscle names required')
        self.x0=policy['x0'];self.u0=policy['u0'];self.K0=policy['K'];self.x_goal=target['x_target'];self.u_goal=target['u_target'];self.K_goal=goal['K']
        n,m=len(self.state_names),len(self.muscle_names)
        if any(a.shape!=shape or not np.isfinite(a).all() for a,shape in ((self.x0,(n,)),(self.x_goal,(n,)),(self.u0,(m,)),(self.u_goal,(m,)),(self.K0,(m,n)),(self.K_goal,(m,n)))):raise ValueError('Invalid program array shape or finite values')
        if any(np.any((u<0)|(u>1)) for u in (self.u0,self.u_goal)):raise ValueError('Program equilibrium excitation outside[0,1]')
        if not np.array_equal(goal['x0'],self.x_goal) or not np.array_equal(goal['u0'],self.u_goal):raise ValueError('Endpoint gain equilibrium differs')
        if not np.array_equal(target['x0'],self.x_goal) or not np.array_equal(target['u0'],self.u_goal):raise ValueError('Target x0/u0 alias interpretation changed')
        self.start_s=_finite(program['start_s']);self.transition_s=_finite(program['transition_s'])
        if self.start_s!=1. or self.transition_s!=5. or isinstance(program['fraction'],bool) or program['fraction']!=1.:raise ValueError('Only accepted1s hold/5s full transfer is supported')
        values={path.rsplit('/',1)[0]:i for i,path in enumerate(self.state_names) if path.endswith('/value')}
        self.speed_pairs=tuple((i,values[path.rsplit('/',1)[0]]) for i,path in enumerate(self.state_names) if path.endswith('/speed'))
        self.segments=(MappingProxyType({'start_s':self.start_s,'duration_s':self.transition_s,'name':'left_weight_transfer_80'}),)
        self.identity={'schema':'ihm.stance-reference-program.v1','bundle_manifest_sha256':_digest(manifest_raw),'program_sha256':_digest(program_raw),
            'baseline_binding':deepcopy(binding),'implementation_sha256':_digest(SOURCE_BYTES)}
        self.model_sha256=_digest(json.dumps(self.identity,sort_keys=True).encode());self.time_s=0.
        for a in (self.x0,self.u0,self.K0,self.x_goal,self.u_goal,self.K_goal):a.flags.writeable=False
        return self
    def reference_at(self,elapsed_s):
        t=_finite(elapsed_s)
        if not 0<=t<=1e9:raise ValueError('Program elapsed time outside[0,1e9]')
        segment=self.segments[0];phase=float(np.clip((t-segment['start_s'])/segment['duration_s'],0,1))
        blend=10*phase**3-15*phase**4+6*phase**5
        rate=(30*phase**2-60*phase**3+30*phase**4)/segment['duration_s']
        delta=self.x_goal-self.x0;reference=self.x0+blend*delta
        for speed,position in self.speed_pairs:reference[speed]+=rate*delta[position]
        equilibrium=self.u0+blend*(self.u_goal-self.u0);gain=self.K0+blend*(self.K_goal-self.K0)
        telemetry={'elapsed_s':t,'segment_index':0,'segment':segment['name'],'phase':'initial_hold' if t<self.start_s else 'transition' if t<self.start_s+self.transition_s else 'endpoint_hold',
            'phase_fraction':phase,'blend':blend,'blend_rate_per_s':rate,'finished_transition':t>=self.start_s+self.transition_s,
            'walking_demonstrated':False,'motor_owner':'engineering_reference_program'}
        return {'state_reference':reference,'equilibrium_excitations':equilibrium,'gain':gain,'telemetry':telemetry}
    def advance(self,dt_s):
        dt=_finite(dt_s)
        if abs(dt-self.dt_s)>1e-12:raise ValueError('Stance program must advance by10ms')
        result=self.reference_at(self.time_s)
        updated=self.time_s+dt
        if updated>1e9:raise ValueError('Program clock exhausted')
        self.time_s=updated;return result
    def checkpoint(self):return {'schema':'ihm.stance-reference-program-state.v1','model_sha256':self.model_sha256,'time_s':self.time_s}
    def restore(self,checkpoint):
        if not isinstance(checkpoint,dict) or checkpoint.get('schema')!='ihm.stance-reference-program-state.v1' or checkpoint.get('model_sha256')!=self.model_sha256:raise ValueError('Stance program checkpoint identity mismatch')
        value=_finite(checkpoint['time_s'])
        if not 0<=value<=1e9:raise ValueError('Program checkpoint time outside range')
        self.time_s=value
    def retain_sources(self,output):
        out=Path(output)/'stance-reference-program';out.mkdir(parents=True,exist_ok=True)
        retained={}
        for origin,raw in self.source_bytes.items():
            name=_digest(raw)[:16]+'-'+Path(origin).name;(out/name).write_bytes(raw);retained[origin]={'retained_name':name,'sha256':_digest(raw)}
        (out/'stance_program.py').write_bytes(SOURCE_BYTES)
        (out/'manifest.json').write_text(json.dumps({'identity':self.identity,'sources':retained},indent=2)+'\n')
        return str(out)
