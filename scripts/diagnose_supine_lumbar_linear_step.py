"""Native finite-difference local reachability diagnostic at a frozen pose."""
from pathlib import Path
import sys,json,tempfile,xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import lsq_linear
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream

def main(path):
    source=Path(path).resolve();candidate=json.loads(source.read_text())
    q=candidate['coordinates'];act=candidate['activations']
    output=Path(tempfile.mkdtemp(prefix='supine-lumbar-linear-',dir=ROOT/'data/derived'))
    stream=NativeMechanicalStream(ROOT,output/'native',environment='supine',target_mass_kg=77.6122029,augmented_registration='data/derived/mechanics/whole_body_lumbar_current/registration.json')
    names=[n for n in q if n not in ('pelvis_ty','pelvis_tz','pelvis_list')];anames=list(act)
    def query(x):
        coordinates={**q,**dict(zip(names,x[:len(names)]))};activations=dict(zip(anames,x[len(names):]));args=['evaluate_static_pose',str(len(coordinates))]
        for n,v in coordinates.items():args.extend([n,str(float(v))])
        args.append(str(len(activations)))
        for n,v in activations.items():args.extend([n,str(float(v))])
        native=stream._request(' '.join(args));a=np.asarray(native['udot']);indices=[native['mobility_coordinate_names'].index(n) for n in ('pelvis_tilt','pelvis_list','pelvis_rotation')]
        return np.r_[a,100*np.asarray(native['com_acceleration_m_s2']),20*a[indices]],native
    try:
        x=np.array([q[n] for n in names]+[act[n] for n in anames]);r,native=query(x);j=np.zeros((len(r),len(x)))
        for i in range(len(x)):
            h=1e-5 if i<len(names) else 1e-4;y=x.copy();y[i]+=h
            j[:,i]=(query(y)[0]-r)/h
        step,*_=np.linalg.lstsq(j,-r,rcond=1e-8)
        singular=np.linalg.svd(j,compute_uv=False)
        radius=np.array([.005 if n=='pelvis_tx' else .03 for n in names]+[.03]*len(anames))
        ranges={c.get('name'):list(map(float,c.findtext('range').split())) for c in ET.parse(output/'native/inputs/subject_walk_scaled.osim').findall('.//Coordinate')}
        ranges.update(candidate.get('posture_bounds_rad',{}))
        bounds=np.array([ranges[n] for n in names]+[[.01,1.] for _ in anames]).T
        limited=lsq_linear(j,-r,bounds=(np.maximum(-radius,bounds[0]-x),np.minimum(radius,bounds[1]-x)),tol=1e-10,max_iter=200)
        report=dict(source=str(source),initial_residual_norm=float(np.linalg.norm(r)),singular_values=singular.tolist(),rank_at_relative1e_8=int(sum(singular>singular[0]*1e-8)),unbounded_predicted_residual_norm=float(np.linalg.norm(j@step+r)),unbounded_step=dict(zip(names+anames,map(float,step))),local_box_predicted_residual_norm=float(np.linalg.norm(j@limited.x+r)),local_box_step=dict(zip(names+anames,map(float,limited.x))),limitation='Local differential reachability only; finite steps may violate source bounds and this probe does not apply or accept a forward state.')
        (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if not k.endswith('step')},indent=2));print(output)
    finally:stream.close()
if __name__=='__main__':main(sys.argv[1])
