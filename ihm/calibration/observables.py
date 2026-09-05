"""Apical extracellular potential gradients, never cell membrane voltages."""
import numpy as np


def convert_field(values, from_unit, to_unit='V/m'):
    factors={'V/m':1.,'mV/mm':1.,'V/mm':1000.,'mV/m':.001}
    if from_unit not in factors or to_unit not in factors:
        raise ValueError('electric field requires voltage per distance units')
    values=np.asarray(values,float)
    if not np.isfinite(values).all(): raise ValueError('field must be finite')
    return values*factors[from_unit]/factors[to_unit]


def lateral_field(surface_potential_v, positions_m, edges):
    v=np.asarray(surface_potential_v,float); positions=np.asarray(positions_m,float)
    if v.ndim!=1 or positions.shape!=(len(v),3) or not np.isfinite(v).all() or not np.isfinite(positions).all():
        raise ValueError('matching finite surface potentials and 3D positions required')
    out=[]
    for a,b in edges:
        if not isinstance(a,(int,np.integer)) or not isinstance(b,(int,np.integer)) or not 0<=a<len(v) or not 0<=b<len(v):
            raise ValueError('edge indices invalid')
        length=np.linalg.norm(positions[b]-positions[a])
        if length<=0: raise ValueError('field requires nonzero distance')
        out.append((v[a]-v[b])/length)
    return np.asarray(out)


def skin_surface_field(model, patch):
    """Bind existing skin model's apical potential states to directional V/m."""
    ids=[f'skin.surface_{i}.potential' for i in range(len(patch.surface))]
    potentials=[model.mean[model.components.index(name)] for name in ids]
    return lateral_field(potentials,[site.position_m for site in patch.surface],[(a,b) for a,b,_ in patch.surface_edges])
