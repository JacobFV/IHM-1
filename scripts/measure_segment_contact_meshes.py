#!/usr/bin/env python3
"""What real segment surfaces cost, and whether the solver keeps their shape.

Three arms over the SAME stance pose, the same mass and the same joint stops, so
the only thing that differs is what the body touches the floor with:

  ``spheres``      what the plant does today: the 12 anatomically placed source
                   foot contacts, plus one COM sphere per non-foot segment whose
                   radius is inscribed in that segment's INERTIA ellipsoid.
  ``mesh_proxy``   the same 12 foot contacts, and the 16 inertia proxies replaced
                   by those segments' real surfaces.
  ``mesh_all``     every segment on its real surface, source foot spheres gone.

Gates, in the sense the log means -- against cases whose answer is known:

* **Momentum balance.** The engine reports
  ``m*(a_com - g) - (contact) - (external)``.  It is the arithmetic identity of
  the whole plant, so it is zero only if the new forces are summed correctly.  A
  contact term that is bookkept wrong shows up here and nowhere else.
* **Standing weight.** At rest on the floor the vertical contact force is the
  body's weight.  Reported next to ``m*g`` and next to the measured vertical
  acceleration, because the source pose is NOT an equilibrium and a number that
  misses weight while the body is accelerating is not a failure -- it is the
  body accelerating, and the two are only distinguishable together.
* **Interpenetration.** The lowest vertex of every contact mesh, in ground, at
  the start and at the end.  A surface that ends the run below the floor is a
  surface the solver let through.
* **Cost.** Median and worst wall seconds per ``advance``.  One advance is an
  error-controlled Simbody integration whose cost is not bounded by dt.
"""
from pathlib import Path
import argparse,json,math,statistics,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream

REGISTRATION='data/models/engineering_stance_v1/registration.json'
TARGET_MASS_KG=77.6122029
FLOOR_Y_M=0.

def read_obj(path):
    rows=[l.split()[1:4] for l in Path(path).read_text().splitlines() if l.startswith('v ')]
    return np.asarray(rows,dtype=float)

def lowest_points(state,bundle):
    """Lowest vertex of every contact mesh, in ground, at this state."""
    manifest=json.loads((bundle/'manifest.json').read_text())
    out={}
    for record in manifest['records']:
        vertices=read_obj(bundle/'meshes'/record['mesh_file'])
        transform=np.asarray(state['bodies'][record['body']]['transform_ground'])
        ground=vertices@transform[:3,:3].T+transform[:3,3]
        out[record['element']]=float(ground[:,1].min())
    return out

def summarise(state):
    contacts=[c for c in state['contacts'] if math.sqrt(sum(v*v for v in c['force_n']))>1e-9]
    total=np.zeros(3)
    for c in state['contacts']:total+=np.asarray(c['force_n'])
    return dict(contact_elements=len(state['contacts']),contact_elements_loaded=len(contacts),
                total_contact_force_n=[float(v) for v in total],
                reported_total_contact_force_n=state['contact_force_n'],
                vertical_contact_force_n=float(total[1]),
                momentum_balance_residual_n=state['momentum_balance_residual_n'],
                momentum_balance_residual_norm_n=float(np.linalg.norm(state['momentum_balance_residual_n'])),
                foot_contact_force_n=state['foot_contact_force_n'],
                contact_model=state['contact_model'],
                segment_contact_meshes=state.get('segment_contact_meshes',0),
                segment_contact_mesh_layer=state.get('segment_contact_mesh_layer'),
                pelvis_ty_m=state['coordinates']['pelvis_ty']['value'],
                loaded=[dict(name=c['name'],body=c['body_frame'],
                             force_n=float(math.sqrt(sum(v*v for v in c['force_n']))),
                             vertical_n=float(c['force_n'][1]),
                             geometry_type=c.get('geometry_type','sphere'),
                             mesh_faces=c.get('mesh_faces'))
                        for c in sorted(contacts,key=lambda c:-abs(c['force_n'][1]))])

def run(arm,bundle,replace_feet,steps,dt,work,pose,stops):
    out=work/('arm-'+arm)
    stream=NativeMechanicalStream(ROOT,out,environment='upright',target_mass_kg=TARGET_MASS_KG,
                                  initial_pose=pose,augmented_registration=REGISTRATION,
                                  coordinate_limits=stops,
                                  segment_contact_meshes=None if bundle is None else str(bundle.relative_to(ROOT)),
                                  segment_contact_replaces_source_feet=replace_feet)
    try:
        started=time.monotonic()
        initial=stream.snapshot()
        report=dict(arm=arm,bundle=None if bundle is None else str(bundle.relative_to(ROOT)),
                    replaces_source_feet=replace_feet,steps=steps,dt_s=dt,
                    weight_n=TARGET_MASS_KG*9.81,initial=summarise(initial))
        if bundle is not None:
            low=lowest_points(initial,bundle)
            report['initial']['lowest_mesh_vertex_above_floor_m']=min(low.values())-FLOOR_Y_M
            report['initial']['lowest_mesh_vertex_element']=min(low,key=low.get)
        costs=[]
        state=initial
        for _ in range(steps):
            mark=time.monotonic()
            state=stream.advance(dt)
            costs.append(time.monotonic()-mark)
        report['final']=summarise(state)
        if bundle is not None:
            low=lowest_points(state,bundle)
            report['final']['lowest_mesh_vertex_above_floor_m']=min(low.values())-FLOOR_Y_M
            report['final']['lowest_mesh_vertex_element']=min(low,key=low.get)
        report['cost']=dict(advances=len(costs),median_s=statistics.median(costs),
                            mean_s=statistics.fmean(costs),worst_s=max(costs),best_s=min(costs),
                            total_wall_s=time.monotonic()-started,
                            simulated_s=state['time_s'],
                            real_time_factor=state['time_s']/max(sum(costs),1e-12))
        return report
    finally:stream.close()

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--steps',type=int,default=50)
    parser.add_argument('--dt',type=float,default=.01)
    parser.add_argument('--work',default='data/derived/segment-contact-measure')
    parser.add_argument('--out',required=True)
    parser.add_argument('--arm',action='append',default=[])
    args=parser.parse_args()
    sys.path.insert(0,str(ROOT/'scripts'))
    import importlib.util
    spec=importlib.util.spec_from_file_location('crawl',ROOT/'scripts/crawl.py')
    crawl=importlib.util.module_from_spec(spec);spec.loader.exec_module(crawl)
    pose=json.loads((ROOT/'data/models/engineering_stance_v1/initial_pose.json').read_text())
    stops=crawl.joint_stops()
    work=ROOT/args.work;work.mkdir(parents=True,exist_ok=True)
    work=Path(work)/('run-'+str(int(time.time())));work.mkdir()
    arms={'spheres':(None,False),
          'mesh_proxy':(ROOT/'data/derived/segment-contact-meshes/stance-bone-proxy',False),
          'mesh_all':(ROOT/'data/derived/segment-contact-meshes/stance-bone-all',True)}
    selected=args.arm or list(arms)
    reports=[]
    for name in selected:
        bundle,replace=arms[name]
        print('== '+name,flush=True)
        try:reports.append(run(name,bundle,replace,args.steps,args.dt,work,pose,stops))
        except Exception as error:reports.append(dict(arm=name,failed=repr(error)))
        print(json.dumps(reports[-1].get('cost',reports[-1]),indent=2),flush=True)
    out=ROOT/args.out;out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(dict(schema='ihm.segment-contact-measurement.v1',
        pose='data/models/engineering_stance_v1/initial_pose.json',
        registration=REGISTRATION,target_mass_kg=TARGET_MASS_KG,steps=args.steps,dt_s=args.dt,
        work=str(work.relative_to(ROOT)),arms=reports),indent=2)+'\n')
    print(str(out))
