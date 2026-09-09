"""Reduced IBM kernel motor materialization, with fixed ports and no readout bias.

This uses IBM's normalized-embedding tanh association rule, not its E/I time
integrator. Only a copied subset of the fused embedding is trainable. The target
is an externally supplied instruction; position/velocity are native feedback.
"""
import hashlib,io
from pathlib import Path
import torch
from torch import nn
from torch.nn import functional as F
SOURCE_BYTES=Path(__file__).read_bytes()

class KernelAnklePolicy(nn.Module):
    def __init__(self, embedding):
        super().__init__()
        if embedding.ndim!=2 or embedding.shape[0]!=8:raise ValueError('Eight kernel sites required')
        self.embed=nn.Parameter(embedding.clone().float())
    def forward(self, observations, sever=False):
        # Four distinct sensory sites (positive/negative error and velocity),
        # four motor sites. Fixed signed readout; no trainable adapter or bias.
        z=F.normalize(self.embed,dim=-1)
        kernel=torch.tanh(2*(z[4:]@z[:4].T))
        if sever:kernel=kernel*0
        sensory=torch.stack((observations[...,0],-observations[...,0],
                             observations[...,1],-observations[...,1]),-1)
        motor=sensory@kernel.T
        differential=(motor[...,0]-motor[...,1]+motor[...,2]-motor[...,3])*.5
        return torch.stack((differential.clamp(0,1),(-differential).clamp(0,1)),-1)
    def commands(self,state,target,sever=False):
        q=state['coordinates']['ankle_angle_r']
        x=torch.tensor([[target-q['value'],q['speed']]],dtype=torch.float32)
        with torch.no_grad():u=self(x,sever=sever)[0].tolist()
        commands={name:0. for name in state['muscles']}
        commands['tibant_r']=u[0];commands['soleus_r']=u[1]
        return commands

def source_embedding(checkpoint):
    path=Path(checkpoint);raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    checkpoint=torch.load(io.BytesIO(raw),map_location='cpu',weights_only=False)
    embedding=checkpoint['dyn.embed']
    # Uniform deterministic site selection, explicitly not anatomical mapping.
    indices=torch.linspace(0,len(embedding)-1,8).long()
    return embedding[indices],{'checkpoint_path':str(path.resolve()),'checkpoint_sha256':digest,
        'site_indices':indices.tolist(),'selection':'uniform index subset; not anatomical ports',
        'original_shape':list(embedding.shape),'sources':checkpoint.get('sources',[])}

def load_motor_policy(path):
    """Load only the separate reduced motor artifact; never rewrite IBM source."""
    return load_motor_policy_bytes(Path(path).read_bytes())

def load_motor_policy_bytes(raw):
    data=torch.load(io.BytesIO(raw),map_location='cpu',weights_only=False)
    if data.get('schema')!='ihm.ibm-native-ankle-kernel.v1':raise ValueError('Unknown motor kernel schema')
    embedding=data['embed']
    if not torch.isfinite(embedding).all():raise ValueError('Nonfinite motor embedding')
    return KernelAnklePolicy(embedding).eval(),data
