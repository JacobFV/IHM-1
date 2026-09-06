"""Work-conjugate projection of external chest loads into the native circuit.

Native lung/chest constitutive elements retain storage and recoil ownership.
The spatial port returns only external load, never a second passive recoil.
The modal displacement distribution is an explicit generic anatomical prior.
"""
import numpy as np


class RespiratoryLoadPort:
    def __init__(self,payload):
        a=np.asarray(payload['volume_jacobian_m2'],float)
        k=np.asarray(payload['stiffness_n_m'],float)
        if a.shape!=(3,) or k.shape!=(3,3) or not np.isfinite(a).all() or not np.isfinite(k).all():
            raise ValueError('Finite three-mode respiratory operators required')
        if not np.allclose(k,k.T) or np.linalg.eigvalsh(k).min()<=0 or (a<=0).any():
            raise ValueError('Positive respiratory operators required')
        direction=np.linalg.solve(k,a)
        self.b=direction/float(a@direction)
        self.bindings={}
        for binding in payload['bindings']:
            ident=binding['entity_id'];basis=np.asarray(binding['translation_basis'],float)
            if ident in self.bindings or basis.shape!=(3,3) or not np.isfinite(basis).all():
                raise ValueError('Unique finite respiratory binding required')
            self.bindings[ident]=basis@self.b

    def evaluate(self,forces):
        generalized=0.
        for force in forces:
            ident=force['id'];vector=np.asarray(force['force_n'],float)
            if ident not in self.bindings or vector.shape!=(3,) or not np.isfinite(vector).all():
                raise ValueError('Force requires a known respiratory owner and finite vector')
            generalized+=float(vector@self.bindings[ident])
        return {'external_pressure_pa':-generalized,
            'mode_displacement_per_volume_m_per_m3':self.b.tolist(),
            'entity_displacement_per_volume_m_per_m3':{k:v.tolist() for k,v in self.bindings.items()},
            'ownership':'External load only; native chest/lung recoil and gas storage unchanged',
            'geometry':'Static modal distribution from canonical anatomy and generic stiffness prior',
            'work_convention':'External work on body = -external_pressure_pa * native_lung_volume_change_m3',
            'limitations':['No independent thoracic tissue recoil returned through this port',
                'Entity centroid translations, not articulated rib rotation or a resolved pressure field']}
