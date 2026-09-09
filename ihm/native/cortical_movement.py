"""Experimental cortical execution of a disclosed phase reference program.

Both feedback error and reference feedforward are sensory current inputs. The
fixed tonic native baseline is never replaced by a scheduled muscle command.
No live LQR gain occurs in advance(); phase decoders are fitted offline.
"""
from copy import deepcopy
from pathlib import Path
import hashlib,io
import numpy as np
import torch
from torch import nn

class CorticalMovementPolicy(nn.Module):
    def __init__(self,dyn,*,state_names,muscle_names,x0,u0,state_scale,feedforward_scale=10.):
        super().__init__();self.dyn=dyn;self.state_names=list(state_names);self.muscle_names=list(muscle_names)
        n,m=len(state_names),len(muscle_names);width=n+m
        if len(set(state_names))!=n or len(set(muscle_names))!=m or not n or not m or dyn.n-width<m+1:raise ValueError('Exact catalogs and disjoint cortical capacity required')
        for name,value,shape in [('x0',x0,(n,)),('u0',u0,(m,)),('state_scale',state_scale,(n,))]:
            tensor=torch.as_tensor(value,dtype=torch.float64)
            if tensor.shape!=shape or not torch.isfinite(tensor).all():raise ValueError('Invalid '+name)
            self.register_buffer(name,tensor.clone())
        if (self.state_scale<=0).any() or (self.u0<.01).any() or (self.u0>1).any() or not np.isfinite(feedforward_scale) or feedforward_scale<=0:raise ValueError('Invalid operating scales/baseline')
        self.register_buffer('feedforward_scale',torch.full((m,),float(feedforward_scale),dtype=torch.float64))
        self.register_buffer('sensory_sites',torch.arange(width));self.register_buffer('motor_sites',torch.arange(width,dyn.n))
        self.register_buffer('decoder_start',torch.zeros(m,dyn.n-width,dtype=torch.float64))
        self.register_buffer('decoder_end',torch.zeros(m,dyn.n-width,dtype=torch.float64))
        self.dyn.double()
        for parameter in self.dyn.parameters():parameter.requires_grad_(False)

    def state(self,batch=1):
        state=tuple(x.double() for x in self.dyn.init_state(batch,self.dyn.embed.device))
        return (*state,*(x.clone() for x in state))

    def advance(self,vector,reference,feedforward,state,*,blend,ticks=10,sever=False,weights=None):
        n,m=len(self.state_names),len(self.muscle_names)
        if vector.ndim!=2 or vector.shape[1]!=n or reference.shape!=vector.shape or feedforward.shape!=(len(vector),m):raise ValueError('Exact cortical movement input widths required')
        if any(not torch.isfinite(x).all() for x in (vector,reference,feedforward)) or not np.isfinite(blend) or not 0<=blend<=1:raise ValueError('Finite bounded movement inputs required')
        if type(ticks) is not int or not 1<=ticks<=1000 or len(state)!=8:raise ValueError('Integral1ms ticks and eight neural tensors required')
        error=(vector-reference)/self.state_scale
        requested=(feedforward-self.u0)/self.feedforward_scale
        features=torch.cat((error,requested),dim=1).clamp(-5,5)
        drive=torch.zeros(len(vector),self.dyn.n,dtype=torch.float64)
        drive[:,self.sensory_sites]=12*torch.tanh(features)
        weights=self.dyn.edge_weights() if weights is None else weights
        if sever:weights=weights*0
        response,baseline=state[:4],state[4:];zero=torch.zeros_like(drive)
        for _ in range(ticks):
            response=self.dyn.step(response,drive,.001,weights)
            baseline=self.dyn.step(baseline,zero,.001,weights)
        motor=response[1][:,self.motor_sites]-baseline[1][:,self.motor_sites]
        anchored=motor-motor[:,:1]
        decoder=self.decoder_start+float(blend)*(self.decoder_end-self.decoder_start)
        decoder=decoder-decoder.mean(-1,keepdim=True)
        output=(self.u0+.2*torch.tanh(torch.nn.functional.linear(10*anchored,decoder))).clamp(.01,1.)
        return output,(*response,*baseline)

    def snapshot_vector(self,snapshot):
        if set(snapshot.get('muscles',{}))!=set(self.muscle_names):raise ValueError('Movement native muscle catalog differs')
        coordinates=snapshot.get('coordinates',snapshot.get('joints',{}));values=[]
        for path in self.state_names:
            parts=path.strip('/').split('/');name,variable=parts[-2:]
            if parts[0]=='jointset' and len(parts)==4 and variable in ('value','speed'):value=coordinates[name][variable]
            elif parts[0]=='forceset' and len(parts)==3 and variable in ('activation','fiber_length'):value=snapshot['muscles'][name]['fiber_length_m' if variable=='fiber_length' else variable]
            else:raise ValueError('Unsupported native state path')
            if isinstance(value,bool) or not np.isfinite(value):raise ValueError('Invalid native movement feedback')
            values.append(float(value))
        return torch.tensor([values],dtype=torch.float64)


def load_cortical_movement(path):
    path=Path(path);raw=path.read_bytes();artifact=torch.load(io.BytesIO(raw),map_location='cpu',weights_only=True)
    if artifact.get('schema')!='ihm.cortical-movement.v1':raise ValueError('Wrong movement artifact schema')
    provenance=artifact['provenance'];equations=(path.parent/'pretrain_video_loop.py').read_bytes()
    if hashlib.sha256(equations).hexdigest()!=provenance['cortical_source_sha256']:raise ValueError('Cortical source differs')
    namespace={'__name__':'_movement_IBM_equations'};exec(compile(equations,str(path.parent/'pretrain_video_loop.py'),'exec'),namespace)
    state=artifact['state_dict'];n,d=state['dyn.embed'].shape;k=state['dyn.idx'].shape[1]
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(0);dyn=namespace['CorticalDynamics'](n,d,k,'cpu').double()
        policy=CorticalMovementPolicy(dyn,state_names=provenance['state_names'],muscle_names=provenance['muscle_names'],x0=state['x0'],u0=state['u0'],state_scale=state['state_scale'])
    policy.load_state_dict(state,strict=True)
    if any(not torch.isfinite(value).all() for value in policy.state_dict().values()):raise ValueError('Nonfinite movement artifact')
    if not torch.equal(policy.sensory_sites,torch.arange(len(policy.state_names)+len(policy.muscle_names))) or not torch.equal(policy.motor_sites,torch.arange(len(policy.sensory_sites),n)):raise ValueError('Unexpected movement port allocation')
    artifact['artifact_sha256']=hashlib.sha256(raw).hexdigest()
    return policy.eval(),artifact
