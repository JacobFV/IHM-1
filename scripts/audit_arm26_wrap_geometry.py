"""Source-only wrap work audit and exact cylindrical-branch reference geometry."""
from pathlib import Path
import hashlib,json,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
BASE='data/raw/mechanics/opensim-core/OpenSim/Simulation/'
FILES=[BASE+p for p in ['Wrap/WrapCylinder.cpp','Wrap/WrapEllipsoid.cpp','Wrap/PathWrapPoint.cpp','Model/GeometryPath.cpp','MomentArmSolver.cpp']]

def cylinder_branch(p1,p2,radius,sense=1):
    """Exact stationary infinite-cylinder geodesic on a selected winding branch.

Caller owns branch/quadrant/finite-cylinder admissibility. Points and result
are in the retained cylinder frame (axis +Z). No branch selection is inferred.
"""
    p1=np.asarray(p1,float);p2=np.asarray(p2,float);r=float(radius)
    if p1.shape!=(3,) or p2.shape!=(3,) or not np.isfinite([p1,p2]).all() or not np.isfinite(r) or r<=0 or sense not in (-1,1):raise ValueError('Invalid cylinder inputs')
    rho=np.array([np.linalg.norm(p1[:2]),np.linalg.norm(p2[:2])])
    if np.any(rho<=r):raise ValueError('Endpoints must lie strictly outside cylinder')
    phi=np.array([np.arctan2(p1[1],p1[0]),np.arctan2(p2[1],p2[0])])
    angles=phi+np.array([sense,-sense])*np.arccos(r/rho)
    theta=(sense*(angles[1]-angles[0]))%(2*np.pi)
    if min(theta,2*np.pi-theta)<1e-8:raise ValueError('Branch boundary needs explicit handling')
    straight=np.sqrt(rho*rho-r*r);arc=r*theta;horizontal=straight.sum()+arc;dz=p2[2]-p1[2]
    z1=p1[2]+dz*straight[0]/horizontal;z2=p1[2]+dz*(straight[0]+arc)/horizontal
    contacts=np.array([[r*np.cos(angles[0]),r*np.sin(angles[0]),z1],[r*np.cos(angles[1]),r*np.sin(angles[1]),z2]])
    gradients=np.array([(p1-contacts[0])/np.linalg.norm(p1-contacts[0]),(p2-contacts[1])/np.linalg.norm(p2-contacts[1])])
    return {'length':float(np.hypot(horizontal,dz)),'contacts':contacts,'endpoint_gradients':gradients,'wrap_length':float(np.hypot(arc,z2-z1)),'angle_rad':float(theta)}

def plane_section_geodesic_curvature(axes,basis,t):
    """Tangential curvature of a central plane section on a triaxial ellipsoid."""
    a=np.asarray(axes,float);e=np.asarray(basis,float)
    x=a*(e[0]*np.cos(t)+e[1]*np.sin(t));d=a*(-e[0]*np.sin(t)+e[1]*np.cos(t));dd=-x
    speed2=d@d;k=dd/speed2-d*(d@dd)/(speed2*speed2)
    normal=x/(a*a);normal/=np.linalg.norm(normal)
    return float(np.linalg.norm(k-normal*(k@normal)))

def audit():
    source={p:{'sha256':hashlib.sha256((ROOT/p).read_bytes()).hexdigest(),'bytes':(ROOT/p).stat().st_size} for p in FILES}
    path='data/derived/mechanics/whole_body_arm26_v2/subject_with_arms.osim';raw=(ROOT/path).read_bytes();root=ET.fromstring(raw);source[path]={'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
    wraps={w.get('name'):{'body':b.get('name'),'type':w.tag,'properties':{c.tag:c.text for c in w if len(c)==0}} for b in root.findall('.//BodySet/objects/Body') for w in b.findall('./WrapObjectSet/objects/*')}
    paths=[]
    for m in root.findall('.//ForceSet/objects/*'):
        if m.get('name','') in ['arm26_'+n+'_'+s for n in ('BRA','BIClong') for s in ('r','l')]:
            ws=m.findall('.//PathWrap');paths.append({'muscle':m.get('name'),'wrap_count':len(ws),'wraps':[{'method':w.findtext('method'),'range':w.findtext('range'),'object':wraps[w.findtext('wrap_object')]} for w in ws]})
    receipt='data/derived/native-wrap-work-ovst3fef/report.json';raw=(ROOT/receipt).read_bytes();r=json.loads(raw)
    return {'schema':'ihm.arm26-wrap-source-audit.v1','native_modified':False,'sources':source,'native_receipt':receipt,'native_receipt_sha256':hashlib.sha256(raw).hexdigest(),'native_evaluations':r['evaluations'],'paths':paths,
      'cylinder_correction':'Selected-branch analytic unrolled cylinder solution, exact axial tangency and length; preserve wrap transform/radius/quadrant/finite bounds, reject inadmissible branches.',
      'ellipsoid_correction':'Stationary constrained surface geodesic with free contact endpoints; integrated arc and endpoint tensions derived from same converged solution. Hybrid plane may seed but cannot constrain final path.',
      'not_claimed':['No actual corrected Arm26 native acceptance or quantified allocation of observed ellipsoid error among chord sum, plane selection and tangent residual.','No model or source-library patch; no muscle exclusion, force-parameter change, or moment-arm torque fudge.']}
if __name__=='__main__':
    out=ROOT/'data/research/arm26_wrap_geometry';out.mkdir(parents=True,exist_ok=True);r=audit();(out/'audit.json').write_text(json.dumps(r,indent=2)+'\n');print('Retained four-path source audit',out)
