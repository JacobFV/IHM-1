"""Experimental multi-muscle stance imitation through persistent IBM E/I dynamics.

Privileged complete native state is projected to disjoint sensory sites. A
bias-free motor decoder reads cortical state; the explicit equilibrium baseline
is retained when association edges are severed. This supplies no balance claim.
"""
from pathlib import Path
import numpy as np
import torch
from torch import nn
from .stance_lqr import NativeStanceLQR


class CorticalStancePolicy(nn.Module):
    def __init__(self,dyn,sensory_sites,motor_sites,*,state_names,muscle_names,x0,u0,state_scale,seed=47,reference_normalization=False,encoder_kind="random_projection"):
        super().__init__();self.dyn=dyn;self.reference_normalization=bool(reference_normalization)
        self.state_names=list(state_names);self.muscle_names=list(muscle_names)
        if not self.state_names or len(set(self.state_names))!=len(self.state_names):raise ValueError('Unique state order required')
        if not self.muscle_names or len(set(self.muscle_names))!=len(self.muscle_names):raise ValueError('Unique muscle order required')
        n,m=len(self.state_names),len(self.muscle_names)
        for name,value,shape in [('x0',x0,(n,)),('u0',u0,(m,)),('state_scale',state_scale,(n,))]:
            value=torch.as_tensor(value,dtype=dyn.embed.dtype)
            if value.shape!=shape or not torch.isfinite(value).all():raise ValueError('Invalid '+name+' dimensions/values')
            self.register_buffer(name,value.clone())
        if (self.state_scale<=0).any() or (self.u0<.01).any() or (self.u0>1).any():raise ValueError('Invalid scaling or native baseline excitation')
        for ports in (sensory_sites,motor_sites):
            if ports.ndim!=1 or ports.dtype!=torch.long or len(torch.unique(ports))!=len(ports) or (ports<0).any() or (ports>=dyn.n).any():
                raise ValueError('Unique in-range integral cortical site indices required')
        self.register_buffer('sensory_sites',sensory_sites.clone());self.register_buffer('motor_sites',motor_sites.clone())
        if not len(sensory_sites) or not len(motor_sites) or set(sensory_sites.tolist())&set(motor_sites.tolist()):
            raise ValueError('Nonempty disjoint cortical ports required')
        generator=torch.Generator().manual_seed(seed)
        if encoder_kind == 'signed_identity':
            if len(sensory_sites)!=2*n or len(motor_sites)<m+1:
                raise ValueError('Full-state signed ports need 2N sensory and at least M+1 motor sites')
            projection=torch.cat((torch.eye(n),-torch.eye(n)),dim=0)
        elif encoder_kind == 'random_projection':
            projection=torch.randn(len(sensory_sites),n,generator=generator)/(n**.5)
        else:raise ValueError('Unknown cortical encoder kind')
        self.encoder_kind=encoder_kind
        self.register_buffer('encoder',projection.to(dyn.embed.dtype))
        self.decoder=nn.Linear(len(motor_sites),m,bias=False,dtype=dyn.embed.dtype)
        nn.init.normal_(self.decoder.weight,std=.005)
        for parameter in dyn.parameters():parameter.requires_grad_(False)
        dyn.embed.requires_grad_(True)

    def state(self,batch=1):
        state=tuple(value.to(self.dyn.embed.dtype) for value in self.dyn.init_state(batch,self.dyn.embed.device))
        return (*state,*(x.clone() for x in state)) if self.reference_normalization else state

    def advance(self,state_vector,state,*,ticks=10,sever=False,weights=None,availability=1.,temperature_factor=1.):
        if state_vector.ndim!=2 or state_vector.shape[1]!=len(self.state_names) or not torch.isfinite(state_vector).all():
            raise ValueError('Exact finite native state vector width required')
        if type(ticks) is not int or ticks<1 or ticks>1000:raise ValueError('Bounded integral1ms ticks required')
        if not np.isfinite(availability) or not 0<=availability<=1 or not np.isfinite(temperature_factor) or not .01<=temperature_factor<=10:raise ValueError('Bounded physiology coupling required')
        normalized=((state_vector-self.x0)/self.state_scale).clamp(-5,5)
        drive=torch.zeros(len(state_vector),self.dyn.n,device=state_vector.device,dtype=state_vector.dtype)
        drive[:,self.sensory_sites]=availability*12*torch.tanh(normalized@self.encoder.T)
        weights=self.dyn.edge_weights() if weights is None else weights
        if sever:weights=weights*0
        expected=8 if self.reference_normalization else 4
        if len(state)!=expected:raise ValueError('Cortical reference-state schema mismatch')
        response=state[:4];reference=state[4:] if self.reference_normalization else None
        zero_drive=torch.zeros_like(drive)
        for _ in range(ticks):
            response=self.dyn.step(response,drive,.001*temperature_factor,weights)
            if reference is not None:reference=self.dyn.step(reference,zero_drive,.001*temperature_factor,weights)
        state=(*response,*reference) if reference is not None else response
        motor=state[1][:,self.motor_sites]
        if reference is not None: motor=motor-reference[1][:,self.motor_sites]
        # Algebraically the same mean-free readout, evaluated as differences
        # from one motor site so identical severed motor rates give exact zero
        # rather than a floating-point reduction residual.
        centered=motor-motor[:,:1]
        readout=self.decoder.weight-self.decoder.weight.mean(-1,keepdim=True)
        correction=.2*torch.tanh(torch.nn.functional.linear(10*centered,readout))
        output=(self.u0+correction).clamp(.01,1.)
        if output.shape!=(len(state_vector),len(self.muscle_names)):raise ValueError('Motor width mismatch')
        return output,state

    def snapshot_vector(self,snapshot):
        if set(snapshot.get('muscles',{}))!=set(self.muscle_names):raise ValueError('Exact native muscle catalog required')
        values=[]
        coordinates=snapshot.get('coordinates',snapshot.get('joints',{}))
        for path in self.state_names:
            parts=path.strip('/').split('/');name,variable=parts[-2:]
            if parts[0]=='jointset' and len(parts)==4 and variable in ('value','speed'):
                value=coordinates[name][variable]
            elif parts[0]=='forceset' and len(parts)==3 and variable in ('activation','fiber_length'):
                value=snapshot['muscles'][name]['fiber_length_m' if variable=='fiber_length' else variable]
            else:raise ValueError('Unsupported native state path')
            if isinstance(value,bool) or not np.isfinite(value):raise ValueError('Invalid native feedback')
            values.append(float(value))
        return torch.tensor([values],dtype=self.x0.dtype)


def load_cortical_stance(path,*,model_sha256=None,dt_s=None):
    """Restore exact dimensions and trained weights against retained IBM source."""
    import hashlib,io
    path=Path(path);raw=path.read_bytes()
    artifact=torch.load(io.BytesIO(raw),map_location='cpu',weights_only=True)
    if artifact.get('schema')!='ihm.ibm-cortical-stance.v1':raise ValueError('Wrong cortical stance artifact schema')
    provenance=artifact['provenance']
    interval=provenance.get('dt_s')
    if isinstance(interval,bool) or not isinstance(interval,(int,float)) or not np.isfinite(interval) or not .001<=interval<=.02 or abs(interval/.001-round(interval/.001))>1e-9:raise ValueError('Invalid artifact neural exchange interval')
    mass=provenance.get('target_mass_kg')
    if mass is not None and (isinstance(mass,bool) or not isinstance(mass,(int,float)) or not np.isfinite(mass) or not 0<mass<=1000):raise ValueError('Invalid artifact patient mass')
    if type(provenance.get('reference_normalization',False)) is not bool:raise ValueError('Invalid neural reference declaration')
    if model_sha256 is not None and provenance['model_sha256']!=model_sha256:raise ValueError('Stance model identity mismatch')
    if dt_s is not None and (isinstance(dt_s,bool) or not isinstance(dt_s,(int,float)) or not np.isfinite(dt_s) or abs(dt_s-provenance['dt_s'])>1e-12):raise ValueError('Stance sample interval mismatch')
    raw_source=(path.parent/'pretrain_video_loop.py').read_bytes()
    if hashlib.sha256(raw_source).hexdigest()!=provenance['cortical_source_sha256']:raise ValueError('Cortical equation source mismatch')
    namespace={'__name__':'_ihm_cortical_stance_equations'}
    exec(compile(raw_source,str(path.parent/'pretrain_video_loop.py'),'exec'),namespace)
    state=artifact['state_dict'];n,d=state['dyn.embed'].shape;k=state['dyn.idx'].shape[1]
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(0)
        dyn=namespace['CorticalDynamics'](n,d,k,'cpu').to(state['dyn.embed'].dtype)
        policy=CorticalStancePolicy(dyn,state['sensory_sites'],state['motor_sites'],
            state_names=provenance['state_names'],muscle_names=provenance['muscle_names'],
            x0=state['x0'],u0=state['u0'],state_scale=state['state_scale'],
            reference_normalization=provenance.get('reference_normalization',False),
            encoder_kind=provenance.get('encoder_kind','random_projection'))
    policy.load_state_dict(state,strict=True)
    if policy.encoder_kind=='signed_identity':
        width=len(policy.state_names)
        expected=torch.cat((torch.eye(width,dtype=policy.encoder.dtype),-torch.eye(width,dtype=policy.encoder.dtype)))
        if not torch.equal(policy.encoder,expected):raise ValueError('Signed identity encoder matrix differs from declared contract')
    if any(not torch.isfinite(v).all() for v in policy.state_dict().values()):raise ValueError('Nonfinite trained cortical stance model')
    artifact['artifact_sha256']=hashlib.sha256(raw).hexdigest()
    return policy.eval(),artifact


class PersistentStanceCommands:
    """Evaluation-only stateful command adapter; full runtime adoption is separate."""
    def __init__(self,policy,*,dt_s=.01):
        if not np.isfinite(dt_s) or abs(dt_s/.001-round(dt_s/.001))>1e-10 or not .001<=dt_s<=.02:
            raise ValueError('Integral1ms exchange required')
        self.policy=policy;self.ticks=round(dt_s/.001);self.cortical_state=policy.state()
        with torch.no_grad():self.weights=policy.dyn.edge_weights().detach()
    def commands(self,snapshot,*,sever=False):
        vector=self.policy.snapshot_vector(snapshot)
        with torch.no_grad():
            output,state=self.policy.advance(vector,self.cortical_state,ticks=self.ticks,sever=sever,weights=self.weights)
        if not torch.isfinite(output).all() or not all(torch.isfinite(x).all() for x in state):raise ValueError('Nonfinite cortical stance rollout')
        self.cortical_state=state
        # float32(.01) is slightly below the native double activation floor.
        # Enforce the physical output domain again after conversion to SI-port
        # Python scalars; this is rounding protection, not an extra controller.
        return dict(zip(self.policy.muscle_names,(max(.01,min(1.,float(v))) for v in output[0])))
