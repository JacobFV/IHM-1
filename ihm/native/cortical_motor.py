"""Experimental ankle readout through actual IBM E/I cortical dynamics.

Disjoint sensory/motor sites require the association edges to carry feedback.
This is a synthetic PD imitation experiment, not a biological motor decoder.
"""
import torch
from torch import nn


class CorticalAnklePolicy(nn.Module):
    def __init__(self, dyn, sensory_sites, motor_sites):
        super().__init__()
        self.dyn = dyn
        self.register_buffer('sensory_sites', sensory_sites.clone())
        self.register_buffer('motor_sites', motor_sites.clone())
        # Fixed signed encoding of position error and angular velocity.
        basis = torch.zeros(len(sensory_sites), 2)
        for i in range(len(sensory_sites)):
            basis[i, (i // 2) % 2] = 1 if i % 2 == 0 else -1
        self.register_buffer('encoder', basis)
        self.decoder = nn.Linear(len(motor_sites), 1, bias=False)
        nn.init.normal_(self.decoder.weight, std=.05)
        for parameter in dyn.parameters(): parameter.requires_grad_(False)
        dyn.embed.requires_grad_(True)

    def state(self, batch=1):
        return self.dyn.init_state(batch, self.dyn.embed.device)

    def advance(self, observations, state, ticks=20, sever=False, availability=1., temperature_factor=1., supplemental_drive=0., weights=None):
        x = observations / observations.new_tensor([.2, 1.])
        x = x.clamp(-2.,2.)
        drive = torch.zeros(observations.shape[0], self.dyn.n, device=observations.device)
        drive[:, self.sensory_sites] = 12*(1 + x @ self.encoder.T) + supplemental_drive
        weights = self.dyn.edge_weights() if weights is None else weights
        if sever: weights = weights*0
        for _ in range(ticks):
            h=.001*temperature_factor
            updated=self.dyn.step(state, drive*availability, h, weights*availability)
            correction=h*20*self.dyn.w_ee*state[1]/self.dyn.r_max*(1-availability)/self.dyn.tau_m
            state=(updated[0]-correction,*updated[1:])
        motor = state[1][:, self.motor_sites]
        centered = motor - motor.mean(-1,keepdim=True)
        differential = .5*torch.tanh(self.decoder(10*centered)[:,0])
        output = torch.stack((differential.clamp_min(0),(-differential).clamp_min(0)),-1)
        return output, state


def load_cortical_policy(path):
    """Load trained copied E/I parameters against retained exact IBM equations."""
    from pathlib import Path
    import hashlib
    import io
    path = Path(path)
    raw = path.read_bytes()
    artifact = torch.load(io.BytesIO(raw), map_location='cpu', weights_only=True)
    if artifact.get('schema') != 'ihm.ibm-cortical-ankle.v1':
        raise ValueError('Unknown cortical motor artifact')
    source = (path.parent/'pretrain_video_loop.py').read_bytes()
    expected = artifact['source_identity']['implementation_sha256']
    if hashlib.sha256(source).hexdigest() != expected:
        raise ValueError('Cortical equation source differs from training')
    namespace = {'__name__':'_ihm_trained_cortical_motor'}
    exec(compile(source, str(path.parent/'pretrain_video_loop.py'), 'exec'), namespace)
    state = artifact['state_dict']
    n,d = state['dyn.embed'].shape
    k = state['dyn.idx'].shape[1]
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(0)
        dyn = namespace['CorticalDynamics'](n,d,k,'cpu')
        policy = CorticalAnklePolicy(dyn,state['sensory_sites'],state['motor_sites'])
    policy.load_state_dict(state,strict=True)
    if any(not torch.isfinite(t).all() for t in policy.state_dict().values()):
        raise ValueError('Nonfinite cortical policy')
    policy.eval()
    artifact['artifact_sha256'] = hashlib.sha256(raw).hexdigest()
    return policy, artifact


class PersistentCorticalCommands:
    """Small evaluation adapter: actual E/I state persists between native ticks."""
    def __init__(self, policy, dt_s=.02):
        self.policy=policy
        self.cortical_state=policy.state()
        self.dt_s=dt_s
    def commands(self, state, target, sever=False):
        q=state.get('joints',state.get('coordinates',{}))['ankle_angle_r']
        x=torch.tensor([[target-q['value'],q['speed']]],dtype=torch.float32)
        with torch.no_grad():
            u,self.cortical_state=self.policy.advance(x,self.cortical_state,
                ticks=round(self.dt_s/.001),sever=sever)
        commands={name:0. for name in state['muscles']}
        commands['tibant_r']=float(u[0,0]);commands['soleus_r']=float(u[0,1])
        return commands
