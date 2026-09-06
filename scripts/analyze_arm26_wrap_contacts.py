"""Bounded offline contact/plane diagnostics from actual native observer JSONL."""
from pathlib import Path
import argparse,json
import numpy as np
from numpy.polynomial.legendre import leggauss
AXES=np.array([.027559229367297551,.02204738349383804,.02204738349383804])

def ellipsoid(row):
    p,q,a,b,c=[np.array(row[k],float) for k in ('p1','p2','r1','r2','c1_unscaled')]
    normal=np.cross(q-p,c-p);normal/=np.linalg.norm(normal);d=normal@p
    points=np.array(row['wrap_points']);plane_error=max(abs(points@normal-d))
    curvature=[]
    for x in points:
        g=2*x/(AXES*AXES);t=np.cross(normal,g);t/=np.linalg.norm(t)
        coefficients=np.linalg.solve([[g@g,g@normal],[g@normal,1.]],[-2*np.sum(t*t/(AXES*AXES)),0.])
        k=coefficients[0]*g+coefficients[1]*normal;ns=g/np.linalg.norm(g);curvature.append(np.linalg.norm(k-(k@ns)*ns))
    n=AXES*normal;y0=d*n/(n@n);radius=np.sqrt(1-d*d/(n@n))
    e1=a/AXES-y0;e1/=np.linalg.norm(e1);ns=n/np.linalg.norm(n);e2=np.cross(ns,e1)
    v=b/AXES-y0;angle=np.arctan2(v@e2,v@e1)
    if len(points)>2:
        first=points[1]/AXES-y0;direction=np.arctan2(first@e2,first@e1)
        if direction>0 and angle<0:angle+=2*np.pi
        if direction<0 and angle>0:angle-=2*np.pi
    projections=np.array([AXES*(y0+radius*e1),AXES*(y0+radius*(e1*np.cos(angle)+e2*np.sin(angle)))])
    lengths=[]
    for count in (32,64):
        nodes,weights=leggauss(count);theta=(nodes+1)*angle/2
        speed=np.linalg.norm(AXES[None,:]*radius*(-np.sin(theta)[:,None]*e1+np.cos(theta)[:,None]*e2),axis=1)
        lengths.append(float(abs(angle)/2*(weights@speed)))
    tangent=[]
    for x,end in [(a,p),(b,q)]:
        n=x/(AXES*AXES);n/=np.linalg.norm(n);v=end-x;v/=np.linalg.norm(v);tangent.append(abs(float(n@v)))
    return {'plane_residual_m':float(plane_error),'surface_implicit_residual_max':float(np.max(abs(np.sum((points/AXES)**2,axis=1)-1))),
       'endpoint_surface_tangency_dot':tangent,'planar_section_geodesic_curvature_max_per_m':float(max(curvature)),
       'projected_section_arc_m':lengths[-1],'quadrature_32_64_difference_m':abs(lengths[1]-lengths[0]),
       'contact_projection_displacement_max_m':float(np.max(np.linalg.norm(projections-np.array([a,b]),axis=1))),
       'projected_arc_minus_native_stored_m':lengths[-1]-row['wrap_length'],
       'scope':'Refined arc uses contacts projected to exact plane section; endpoint projection displacement is separate uncertainty. Nonzero section geodesic curvature indicates nonstationary surface route.'}

def cylinder(row):
    p,q,a,b=[np.array(row[k],float) for k in ('p1','p2','r1','r2')]
    incoming=(a-p)/np.linalg.norm(a-p);outgoing=(q-b)/np.linalg.norm(q-b);slope=(b[2]-a[2])/row['wrap_length']
    return {'incoming_axial_tangency_residual':float(incoming[2]-slope),'outgoing_axial_tangency_residual':float(outgoing[2]-slope),'minimum_cap_margin_m':float(.055118458734595102/2-max(abs(a[2]),abs(b[2])))}

def analyze(out):
    out=Path(out);result={}
    for mode in ('observed','candidate'):
        rows=[]
        for line in (out/f'{mode}.jsonl').read_text().splitlines():
            r=json.loads(line)
            if not r.get('wrapped'):continue
            item={'case':r['case'],'object':r['object'],'type':r['type']}
            try:item.update(ellipsoid(r) if 'Ellipsoid' in r['type'] else cylinder(r))
            except (ValueError,np.linalg.LinAlgError,FloatingPointError) as e:item['analysis_rejected']=str(e)
            rows.append(item)
        result[mode]=rows
    (out/'contact_analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();r=analyze(a.output);print({k:len(v) for k,v in r.items()})
