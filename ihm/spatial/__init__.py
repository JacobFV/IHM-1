"""Physical coordinate frames, rigid registration and conservative exchange maps."""
from dataclasses import dataclass, field
import numpy as np

LENGTH_UNITS={'m':1.,'cm':.01,'mm':.001,'um':1e-6}

@dataclass(frozen=True)
class Frame:
    id: str
    unit: str
    specimen: str
    reference: str=''
    def __post_init__(self):
        if not self.id or not self.specimen or self.unit not in LENGTH_UNITS:
            raise ValueError('frame needs identity, specimen and an explicit supported length unit')
    @property
    def meters_per_unit(self):return LENGTH_UNITS[self.unit]

@dataclass(frozen=True)
class Registration:
    source: Frame
    target: Frame
    matrix_m: np.ndarray
    rms_residual_m: float
    reference: str
    evidence: dict=field(default_factory=dict)
    def __post_init__(self):
        m=np.array(self.matrix_m,dtype=float,copy=True)
        if m.shape!=(4,4) or not np.isfinite(m).all() or not np.allclose(m[3],[0,0,0,1],atol=1e-12,rtol=0):
            raise ValueError('registration requires finite homogeneous 4x4 matrix')
        r=m[:3,:3]
        if not np.allclose(r.T@r,np.eye(3),rtol=1e-10,atol=1e-10) or not np.isclose(np.linalg.det(r),1,atol=1e-10):
            raise ValueError('physical registration requires a proper rigid rotation; scaling needs a separate deformation model')
        if not np.isfinite(self.rms_residual_m) or self.rms_residual_m<0 or not self.reference:
            raise ValueError('registration needs residual and source reference')
        m.flags.writeable=False;object.__setattr__(self,'matrix_m',m)
    def apply(self,points):
        p=_points(points)*self.source.meters_per_unit
        return (p@self.matrix_m[:3,:3].T+self.matrix_m[:3,3])/self.target.meters_per_unit
    def inverse(self):
        return Registration(self.target,self.source,np.linalg.inv(self.matrix_m),self.rms_residual_m,self.reference,dict(self.evidence))
    def to_dict(self):
        from dataclasses import asdict
        return {'source':asdict(self.source),'target':asdict(self.target),'matrix_m':self.matrix_m.tolist(),
                'rms_residual_m':self.rms_residual_m,'reference':self.reference,'evidence':self.evidence}
    @classmethod
    def from_dict(cls,d):return cls(Frame(**d['source']),Frame(**d['target']),d['matrix_m'],d['rms_residual_m'],d['reference'],d.get('evidence',{}))

def _points(p):
    p=np.asarray(p,dtype=float)
    if p.ndim!=2 or p.shape[1]!=3 or len(p)==0 or not np.isfinite(p).all():raise ValueError('expected nonempty finite Nx3 points')
    return p

def fit_landmarks(source_points,target_points,source,target,*,reference,weights=None,allow_cross_specimen=False):
    """Weighted proper-rotation Kabsch in SI; does not invent correspondences."""
    x=_points(source_points)*source.meters_per_unit;y=_points(target_points)*target.meters_per_unit
    if x.shape!=y.shape or len(x)<3:raise ValueError('matching sets need at least three landmarks')
    if source.specimen!=target.specimen and not allow_cross_specimen:raise ValueError('cross-specimen registration requires explicit authorization in the fit')
    w=np.ones(len(x)) if weights is None else np.asarray(weights,dtype=float)
    if w.shape!=(len(x),) or not np.isfinite(w).all() or (w<=0).any():raise ValueError('landmark weights must be finite positive')
    w=w/w.max();w=w/w.sum();cx=w@x;cy=w@y;xc=x-cx;yc=y-cy
    if min(np.linalg.matrix_rank(xc),np.linalg.matrix_rank(yc))<2:raise ValueError('collinear landmarks cannot identify a rotation')
    u,s,vt=np.linalg.svd((xc*w[:,None]).T@yc);d=np.eye(3);d[2,2]=np.linalg.det(vt.T@u.T)
    r=vt.T@d@u.T;t=cy-r@cx;residual=x@r.T+t-y
    m=np.eye(4);m[:3,:3]=r;m[:3,3]=t
    rms=float(np.sqrt(np.sum(w*np.sum(residual**2,axis=1))))
    return Registration(source,target,m,rms,reference,{'landmarks':len(x),'weights':w.tolist(),'residuals_m':residual.tolist(),
        'singular_values_m2':s.tolist(),'cross_specimen':source.specimen!=target.specimen,
        'validated_out_of_sample':False,'uncertainty':'fit residuals; not independent landmark measurement uncertainty'})

@dataclass(frozen=True)
class ConservativeMap:
    weights: np.ndarray
    def __post_init__(self):
        w=np.array(self.weights,float,copy=True)
        if w.ndim!=2 or not w.size or not np.isfinite(w).all() or (w<0).any() or not np.allclose(w.sum(axis=1),1,atol=1e-12,rtol=0):
            raise ValueError('exchange weights must be nonnegative with rows summing to one')
        w.flags.writeable=False;object.__setattr__(self,'weights',w)
    def intensive(self,target_values):
        """Interpolate target intensive values onto source exchange ports."""
        v=np.asarray(target_values,float)
        if v.shape!=(self.weights.shape[1],) or not np.isfinite(v).all():raise ValueError('invalid target values')
        return self.weights@v
    def extensive(self,source_amounts):
        """Distribute signed extensive source amounts; preserves their total."""
        v=np.asarray(source_amounts,float)
        if v.shape!=(self.weights.shape[0],) or not np.isfinite(v).all():raise ValueError('invalid source amounts')
        return self.weights.T@v

def conservative_map(weights):return ConservativeMap(weights)
