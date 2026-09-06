"""Bounded zero-thickness node/triangle tunneling guard for linear trajectories.

This computes coplanarity roots, not collision response or full cloth CCD.
Edge, persistent coplanarity, tangency and degenerate trajectories are unresolved.
"""
import numpy as np
from functools import reduce
from numpy.polynomial import polynomial as poly


def swept_vertex_face(point_start,point_end,triangle_start,triangle_end):
    a,b=np.asarray(point_start,float),np.asarray(point_end,float)
    t,u=np.asarray(triangle_start,float),np.asarray(triangle_end,float)
    if a.shape!=(3,) or b.shape!=(3,) or t.shape!=(3,3) or u.shape!=(3,3) or not all(np.isfinite(x).all() for x in (a,b,t,u)):
        raise ValueError('Finite point and triangle endpoints required')
    # Translation/scale normalization reduces dimensional polynomial conditioning.
    scale=max(float(np.linalg.norm(t[1:]-t[0],axis=1).max()),float(np.linalg.norm(u[1:]-u[0],axis=1).max()))
    if scale<=1e-15:return dict(status='unresolved',reason='degenerate_triangle')
    origin=t[0].copy();a=(a-origin)/scale;b=(b-origin)/scale;t=(t-origin)/scale;u=(u-origin)/scale
    e=t[1]-t[0];f=t[2]-t[0];de=u[1]-u[0]-e;df=u[2]-u[0]-f
    normal=np.array([np.cross(e,f),np.cross(de,f)+np.cross(e,df),np.cross(de,df)])
    norm2=reduce(poly.polyadd,(poly.polymul(normal[:,j],normal[:,j]) for j in range(3)))
    extrema=[0.,1.]+[float(r.real) for r in poly.polyroots(poly.polyder(norm2)) if abs(r.imag)<1e-9 and 0<r.real<1]
    if min(poly.polyval(r,norm2) for r in extrema)<=1e-20:
        return dict(status='unresolved',reason='triangle_degenerates_during_interval')
    relative=np.array([a-t[0],b-a-(u[0]-t[0])])
    coefficients=reduce(poly.polyadd,(poly.polymul(relative[:,j],normal[:,j]) for j in range(3)))
    tolerance=1e-10*max(1.,float(np.max(np.abs(coefficients))))
    if np.max(np.abs(coefficients))<=tolerance:return dict(status='unresolved',reason='coplanar_or_near_coplanar_interval')
    while len(coefficients)>1 and abs(coefficients[-1])<1e-14*np.max(np.abs(coefficients)):coefficients=coefficients[:-1]
    events=[]
    for root in poly.polyroots(coefficients):
        if abs(root.imag)>1e-7 or root.real< -1e-10 or root.real>1+1e-10:continue
        if abs(root.imag)>1e-10:return dict(status='unresolved',reason='near_multiple_root')
        s=float(np.clip(root.real,0,1));tri=t+s*(u-t);p=a+s*(b-a)
        edge=np.stack((tri[1]-tri[0],tri[2]-tri[0]),axis=1)
        weights=np.linalg.lstsq(edge,p-tri[0],rcond=None)[0];beta=np.r_[1-weights.sum(),weights]
        if np.min(beta)<-1e-9:continue
        if np.min(beta)<=1e-9:return dict(status='unresolved',reason='edge_or_vertex_event',fraction=s)
        slope=float(poly.polyval(s,poly.polyder(coefficients)))
        if abs(slope)<=tolerance:return dict(status='unresolved',reason='tangent_or_ill_conditioned_root',fraction=s)
        # An initially touching point moving away from outward-oriented skin
        # needs no new impact. Persistent/tangent contact was rejected above.
        if s<=1e-10 and slope>0:continue
        n=poly.polyval(s,normal);n/=np.linalg.norm(n)
        events.append(dict(status='crossing',fraction=s,barycentric=beta.tolist(),normal=n.tolist(),
                           point_m=(origin+scale*p).tolist(),direction='inward' if slope<0 else 'outward'))
    if len(events)>1:return dict(status='unresolved',reason='multiple_crossings',events=events)
    return events[0] if events else dict(status='clear',scope='no interior node-face crossing under linear zero-thickness motion')
