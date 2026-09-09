"""Experimental muscle-only feedback from an exact native local linearization.

This controller makes no nonlinear balance guarantee. Its output is muscle
excitation, subject to physical clipping; it does not prescribe coordinates.
"""
from pathlib import Path
import numpy as np


class NativeStanceLQR:
    def __init__(self, artifact, *, model_sha256=None, dt_s=None, target_mass_kg=None):
        """Wrappers should pass model_sha256 and dt_s to bind the exact plant.

        Native snapshots do not contain a source model digest. Catalog equality
        is always checked; it cannot distinguish different poses/parameters of
        the same muscle model. Passing the digest supplies that binding.
        """
        with np.load(Path(artifact), allow_pickle=False) as data:
            self.K=data['K'].copy();self.x0=data['x0'].copy();self.u0=data['u0'].copy()
            self.state_names=data['state_names'].tolist();self.muscle_names=data['muscle_names'].tolist()
            self.model_sha256=str(data['model_sha256'].item()) if 'model_sha256' in data else None
            self.dt_s=float(data['dt_s'].item()) if 'dt_s' in data else None
            self.minimum_activation=data['minimum_activation'].copy() if 'minimum_activation' in data else None
            self.target_mass_kg=float(data['target_mass_kg'].item()) if 'target_mass_kg' in data else None
        if not isinstance(self.state_names,list) or not isinstance(self.muscle_names,list):
            raise ValueError('LQR names must be one-dimensional arrays')
        for names in (self.state_names,self.muscle_names):
            if not names or any(not isinstance(name,str) or not name for name in names) or len(set(names))!=len(names):
                raise ValueError('LQR names must be unique nonempty strings')
        n,m=len(self.state_names),len(self.muscle_names)
        if self.K.shape != (m,n) or self.x0.shape!=(n,) or self.u0.shape!=(m,):
            raise ValueError('LQR artifact dimensions disagree')
        if not all(np.isfinite(v).all() for v in (self.K,self.x0,self.u0)) or not np.any(self.K):
            raise ValueError('LQR artifact is nonfinite or has no solved feedback')
        if np.any((self.u0<0)|(self.u0>1)):raise ValueError('LQR baseline excitation outside [0,1]')
        if self.minimum_activation is not None and (self.minimum_activation.shape!=(m,) or not np.isfinite(self.minimum_activation).all() or np.any((self.minimum_activation<0)|(self.minimum_activation>1))):
            raise ValueError('Invalid native minimum activation vector')
        if self.dt_s is not None and (not np.isfinite(self.dt_s) or self.dt_s<=0):
            raise ValueError('Invalid LQR sampling interval')
        if self.target_mass_kg is not None and (not np.isfinite(self.target_mass_kg) or self.target_mass_kg<=0):
            raise ValueError('Invalid LQR target mass')
        if model_sha256 is not None and model_sha256!=self.model_sha256:
            raise ValueError('LQR source model identity mismatch')
        if dt_s is not None and (isinstance(dt_s,(bool,np.bool_)) or self.dt_s is None or not np.isclose(dt_s,self.dt_s,rtol=0,atol=1e-12)):
            raise ValueError('LQR sampling interval mismatch')
        if target_mass_kg is not None and (isinstance(target_mass_kg,(bool,np.bool_)) or self.target_mass_kg is None or not np.isclose(target_mass_kg,self.target_mass_kg,rtol=0,atol=1e-9)):
            raise ValueError('LQR target mass mismatch or unknown artifact mass')
        for path in self.state_names:
            parts=path.strip('/').split('/')
            if len(parts)!=4 and len(parts)!=3:raise ValueError('Unsupported native state path: '+path)
            if not ((parts[0]=='jointset' and len(parts)==4 and parts[-1] in ('value','speed')) or
                    (parts[0]=='forceset' and len(parts)==3 and parts[-1] in ('activation','fiber_length') and parts[-2] in self.muscle_names)):
                raise ValueError('Unsupported native state path: '+path)

    def state_vector(self, snapshot):
        """Require every state by name; never replace missing afference by zero."""
        if set(snapshot['muscles'])!=set(self.muscle_names):
            raise ValueError('Native snapshot muscle catalog differs from LQR artifact')
        values=[]
        for path in self.state_names:
            parts=path.strip('/').split('/');name,variable=parts[-2:]
            if parts[0]=='jointset':
                value=snapshot['coordinates'][name][variable]
            elif parts[0]=='forceset' and variable in ('activation','fiber_length'):
                value=snapshot['muscles'][name]['fiber_length_m' if variable=='fiber_length' else variable]
            else:
                raise ValueError('Unsupported native state path: '+path)
            values.append(float(value))
        result=np.asarray(values)
        if not np.isfinite(result).all():raise ValueError('Nonfinite native feedback')
        return result

    def commands(self, snapshot, *, gain=1.0):
        if isinstance(gain,(bool,np.bool_)):raise ValueError('Boolean feedback gain is invalid')
        gain=float(gain)
        if not np.isfinite(gain) or gain<0:raise ValueError('Finite nonnegative gain required')
        error=self.state_vector(snapshot)-self.x0
        raw=self.u0-gain*self.K@error
        clipped=np.clip(raw,0.0,1.0)
        floor=self.minimum_activation
        return dict(zip(self.muscle_names,map(float,clipped))), {
            'clipped_muscles':int(np.count_nonzero(raw!=clipped)),
            'raw_excitation_min':float(raw.min()),'raw_excitation_max':float(raw.max()),
            'state_error_norm_mixed_units':float(np.linalg.norm(error)),
            'minimum_activation_known':floor is not None,
            'requested_below_minimum_activation_muscles':None if floor is None else int(np.count_nonzero(raw<floor)),
            'wire_below_minimum_activation_muscles':None if floor is None else int(np.count_nonzero(clipped<floor)),
            'minimum_requested_activation_margin':None if floor is None else float(np.min(raw-floor)),
            'controller':'local native LQR; no nonlinear acceptance implied',
        }
