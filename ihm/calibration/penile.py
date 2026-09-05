"""Published human ex vivo constitutive fits, with explicit anatomical scope.

Khorshidi et al. (2024), equations1–3/Table1. Returns energy per reference
volume (J/m³) and first Piola stress (Pa), conjugate to deformation gradient.
Geometry, fiber architecture, density and living tissue calibration are separate.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/measurements/biomechanics/khorshidi_2024.json'

class PenileMaterial:
    def __init__(self,tissue,evidence,fiber):
        self.tissue=tissue;self._evidence=deepcopy(evidence)
        self._p=deepcopy(evidence['materials'][tissue])
        self._fiber=None if fiber is None else np.array(fiber,float,copy=True)

    @property
    def parameters(self):return deepcopy(self._p)

    def describe(self):
        return dict(tissue=self.tissue,parameters=self.parameters,
                    fiber_direction=None if self._fiber is None else self._fiber.tolist(),
                    doi=self._evidence['doi'],cohort=deepcopy(self._evidence['cohort']),
                    limitations=deepcopy(self._evidence['limitations']),
                    estimate_kind=self._evidence['estimate_kind'],parameter_covariance=None,
                    source_sha256=self._evidence['source_sha256'],whole_body_mapping=None)

    def initial_moduli(self):
        p=self._p
        return dict(shear_pa=sum(p['mu_pa']) if p['law']=='ogden' else 2*p['C10_pa'],
                    bulk_pa=2/p['D_pa_inverse'],
                    scope='isotropic matrix tangent; anisotropic fiber tangent excluded')

    def energy_piola(self,deformation):
        F=np.asarray(deformation,float)
        if F.shape[-2:]!=(3,3) or not np.isfinite(F).all():raise ValueError('Finite 3x3 deformation gradients required')
        J=np.linalg.det(F)
        if np.any(J<=0):raise ValueError('Positive deformation determinant required')
        try:
            with np.errstate(over='raise',invalid='raise',divide='raise'):
                inverse_t=np.swapaxes(np.linalg.inv(F),-1,-2);p=self._p;D=p['D_pa_inverse']
                if p['law']=='ogden':
                    u,s,vt=np.linalg.svd(F);bar=s/J[...,None]**(1/3)
                    energy=np.zeros_like(J);principal=np.zeros_like(s)
                    for mu,alpha in zip(p['mu_pa'],p['alpha']):
                        power=bar**alpha
                        energy+=2*mu/alpha**2*(power.sum(axis=-1)-3)
                        principal+=2*mu/alpha*(power-power.mean(axis=-1,keepdims=True))/s
                    stress=(u*principal[...,None,:])@vt
                else:
                    factor=J**(-2/3);I1=factor*np.sum(F*F,axis=(-2,-1))
                    gradient_I1=2*factor[...,None,None]*F-(2/3)*I1[...,None,None]*inverse_t
                    energy=p['C10_pa']*(I1-3);stress=p['C10_pa']*gradient_I1
                    if p['law']=='hgo':
                        Fa=F@self._fiber;I4=factor*np.sum(Fa*Fa,axis=-1)
                        gradient_I4=2*factor[...,None,None]*Fa[..., :,None]*self._fiber-(2/3)*I4[...,None,None]*inverse_t
                        k=p['kappa'];strain=np.maximum(0,k*(I1-3)+(1-3*k)*(I4-1))
                        exponential=np.exp(p['k2']*strain*strain)
                        energy+=p['k1_pa']/(2*p['k2'])*np.expm1(p['k2']*strain*strain)
                        stress+=(p['k1_pa']*strain*exponential)[...,None,None]*(k*gradient_I1+(1-3*k)*gradient_I4)
                if p['law']=='hgo':
                    energy+=((J*J-1)/2-np.log(J))/D
                    stress+=((J*J-1)/D)[...,None,None]*inverse_t
                else:
                    energy+=(J-1)**2/D
                    stress+=(2*(J-1)*J/D)[...,None,None]*inverse_t
                if not np.isfinite(energy).all() or not np.isfinite(stress).all():raise ValueError('Nonfinite constitutive state')
                return energy,stress
        except (FloatingPointError,np.linalg.LinAlgError) as error:
            raise ValueError('Constitutive deformation is numerically unsupported') from error

def penile_material(tissue,*,fiber_direction=None,evidence_root=ROOT):
    evidence_root=Path(evidence_root).resolve()
    evidence=json.loads((evidence_root/'data/measurements/biomechanics/khorshidi_2024.json').read_text())
    if tissue not in evidence['materials']:raise ValueError('No fitted tissue parameters for '+str(tissue))
    source=(evidence_root/evidence['source_pdf']).resolve()
    if not source.is_relative_to(evidence_root) or hashlib.sha256(source.read_bytes()).hexdigest()!=evidence['source_sha256']:
        raise ValueError('Primary full-text identity changed')
    fiber=None
    if tissue=='tunica_albuginea':
        if fiber_direction is None:raise ValueError('Explicit reference-space fiber direction required')
        fiber=np.asarray(fiber_direction,float)
        if fiber.shape!=(3,) or not np.isfinite(fiber).all() or not np.isclose(np.linalg.norm(fiber),1.,rtol=0,atol=1e-8):
            raise ValueError('Unit reference-space fiber direction required')
    elif fiber_direction is not None:raise ValueError('Isotropic material cannot accept a fiber mapping')
    return PenileMaterial(tissue,evidence,fiber)
