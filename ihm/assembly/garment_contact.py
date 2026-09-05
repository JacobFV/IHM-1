"""Physical generated-short-panel / source-shaped tissue experiment.

This is a constrained front-panel experiment, not whole-garment donning or
calibrated genital containment. The original generated garment is retained.
"""
from pathlib import Path
import hashlib,json,math,time
import numpy as np
from .clothing import Cloth
from .contact_dynamics import DynamicTetrahedra,NodeTriangleContact,step_coupled,_triangle_projection

BASE=Path(__file__).resolve().parents[2]

def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def prepare_panel(garment_path,*,prestrain=.04,stiffness_n_m=10.):
    source=json.loads(Path(garment_path).read_text());garment=next(g for g in source['garments'] if g['id']=='shorts')
    x=np.asarray(garment['positions'],float).reshape(-1,3);tri=np.asarray(garment['indices'],int).reshape(-1,3)
    region=(np.abs(x[:,0])<.075)&(x[:,1]>-.145)&(x[:,1]<.041)&(x[:,2]>.035)
    selected_faces=np.flatnonzero(np.all(region[tri],axis=1));faces=tri[selected_faces]
    source_nodes,reverse=np.unique(faces,return_inverse=True);x=x[source_nodes];tri=reverse.reshape(-1,3)
    # Target normals point toward the body, so material behind the panel has a
    # positive contact gap. This changes winding only, not garment geometry.
    normal=np.cross(x[tri[:,1]]-x[tri[:,0]],x[tri[:,2]]-x[tri[:,0]])
    if np.median(normal[:,2])>0:tri=tri[:,[0,2,1]]
    if np.any(np.cross(x[tri[:,1]]-x[tri[:,0]],x[tri[:,2]]-x[tri[:,0]])[:,2]>=0):raise ValueError('Panel is not an unambiguous anterior surface')
    edges=np.concatenate((tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]))
    unique,count=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True);fixed=np.unique(unique[count==1])
    center=x.mean(axis=0);reference=center+(1-prestrain)*(x-center)
    cloth=Cloth(reference,tri,areal_density_kg_m2=.18,edge_stiffness_n_m=stiffness_n_m)
    cloth.position_m=x.copy();cloth.fixed_nodes=fixed;cloth.material_owner_ids=('generated-shorts-front-panel',)
    return cloth,{'source_path':str(Path(garment_path).resolve()),'source_sha256':_sha(garment_path),
        'source_skin':source['source'],'garment_constructor_sha256':source['constructor_sha256'],
        'source_node_indices':source_nodes.tolist(),'source_triangle_indices':selected_faces.tolist(),
        'panel_node_count':len(x),'panel_face_count':len(tri),'fixed_boundary_node_count':len(fixed),
        'prestrain':prestrain,'areal_density_kg_m2':.18,'edge_stiffness_n_m':stiffness_n_m,
        'prestrain_basis':'uniform intrinsic edge-rest-length contraction from generated panel; explicit engineering prior',
        'support_basis':'extracted patch boundary fixed at generated-short positions; unresolved surrounding garment/pelvis represented by prescribed far field'}

def _penetration(tissue,panel,node_ids):
    maximum=0.;unresolved=0
    triangles=panel.position_m[panel.triangles]
    for i in node_ids:
        gap,_,_,inside,d2=_triangle_projection(tissue.position_m[i],triangles)
        k=int(np.argmin(d2))
        if not inside[k]:unresolved+=1;continue
        maximum=max(maximum,-float(gap[k]))
    return max(maximum,0.),unresolved

def run_garment_tissue_experiment(output_dir,*,seconds=.24,friction=True,contact=True,
                                dt_s=.00005,prestrain=.04,contact_iterations=2,
                                domain_dir=None,garment_path=None):
    if not isinstance(friction,bool) or not isinstance(contact,bool):raise ValueError('Boolean experimental switches required')
    for value in (seconds,dt_s,prestrain):
        if isinstance(value,bool) or not math.isfinite(value) or value<=0:raise ValueError('Positive finite experimental parameters required')
    if prestrain>.15 or seconds>.5 or dt_s>.00005 or contact_iterations not in (1,2,4,8):raise ValueError('Outside supported exploratory fixture range')
    if abs(seconds/dt_s-round(seconds/dt_s))>1e-8:raise ValueError('Final time must align to step clock')
    if abs(.005/dt_s-round(.005/dt_s))>1e-8 or abs(seconds/.005-round(seconds/.005))>1e-8:raise ValueError('Step and horizon must align to the 5 ms observation clock')
    out=Path(output_dir)
    if out.exists():raise ValueError('Choose a fresh garment/tissue experiment directory')
    domain=Path(domain_dir or BASE/'data/derived/material-domains/pelvis-0.004m');manifest=json.loads((domain/'manifest.json').read_text())
    for name,digest in manifest['artifacts'].items():
        if _sha(domain/name)!=digest:raise ValueError('Changed material domain artifact')
    data=np.load(domain/'pelvic-domain.npz',allow_pickle=False);x=data['vertices_m'];tetrahedra=data['tetrahedra'];labels=data['material_index']
    owners=tuple(m['source_id'] for m in manifest['material_regions'])
    tissue=DynamicTetrahedra(x,tetrahedra,mu_pa=data['mu_pa'],lambda_pa=data['lambda_pa'],density_kg_m3=data['density_kg_m3'],fixed_nodes=data['fixed_nodes'],material_owner_ids=owners)
    panel,panel_metadata=prepare_panel(garment_path or BASE/'data/derived/clothing/garments.json',prestrain=prestrain)
    if dt_s>min(tissue.max_explicit_dt_s,panel.max_explicit_dt_s):raise ValueError('Step exceeds fixture elastic bounds')
    glans=np.unique(tetrahedra[labels==2]);front=tissue.surface_nodes[x[tissue.surface_nodes,2]>.05]
    free_panel=np.setdiff1d(np.arange(len(panel.position_m)),panel.fixed_nodes)
    friction_static=.4 if friction else 0.;friction_kinetic=.3 if friction else 0.
    pairs=[NodeTriangleContact('tissue','panel',friction_static,friction_kinetic,.003,node_ids=front),
           NodeTriangleContact('panel','tissue',friction_static,friction_kinetic,.003,node_ids=free_panel)]*contact_iterations if contact else []
    out.mkdir(parents=True);(out/'source-snapshot').mkdir()
    code_sources=[Path(__file__),BASE/'ihm/assembly/contact_dynamics.py',BASE/'ihm/assembly/clothing.py',BASE/'ihm/assembly/mechanics_backend.py']
    for path in code_sources:(out/'source-snapshot'/path.name).write_bytes(path.read_bytes())
    configuration={'seconds':seconds,'dt_s':dt_s,'friction':friction,'contact':contact,'contact_iterations':contact_iterations,
        'friction_static':friction_static,'friction_kinetic':friction_kinetic,'friction_basis':'engineering coefficients; not calibrated penile skin/textile data',
        'load':{'start_s':.12,'end_s':.18,'total_force_n':[0,.02,0],'target':'glans nodes'},
        'gravity_m_s2':[0,0,0],'panel':panel_metadata,'tissue_manifest_path':str(domain/'manifest.json'),
        'tissue_manifest_sha256':_sha(domain/'manifest.json'),'material_owner_ids':owners,
        'source_sha256':{str(p.relative_to(BASE)):_sha(p) for p in code_sources}}
    (out/'configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    tissue_frames=[x.copy()];panel_frames=[panel.position_m.copy()];times=[0.];frames=[]
    tissue_velocities=[tissue.velocity_m_s.copy()];panel_velocities=[panel.velocity_m_s.copy()]
    tissue_contact_force=[np.zeros_like(x)];panel_contact_force=[np.zeros_like(panel.position_m)]
    interval_tissue_impulse=np.zeros_like(x);interval_panel_impulse=np.zeros_like(panel.position_m)
    initial_energy=tissue.elastic_forces(tissue.position_m)[1]+panel.elastic_forces(panel.position_m)[1]
    external_work=0.
    minimum_j=1.;momentum=0.;pair_residual=0.;pair_work_error=0.;count=0;edges=0;defect=0.;loss=0.;tissue_work=0.;panel_work=0.;tissue_impulse=np.zeros(3);absolute_tissue_impulse=0.
    elapsed=time.monotonic()
    with (out/'frames.jsonl').open('w') as frame_log:
        for step in range(round(seconds/dt_s)):
            at=step*dt_s;force=np.zeros_like(x)
            if .12<=at<.18:force[glans,1]=.02/len(glans)
            result=step_coupled({'tissue':tissue,'panel':panel},dt_s,contacts=pairs,
                                external_forces_n={'tissue':force},gravity_m_s2=(0,0,0))
            minimum_j=min(minimum_j,tissue.minimum_jacobian());momentum=max(momentum,float(np.linalg.norm(result['momentum_residual_ns'])))
            defect+=result['numerical_energy_defect_j'];loss+=result['contact_dissipation_j']
            external_work+=result['external_work_j']
            for pair,r in zip(pairs,result['contacts']):
                count+=r['contact_count'];edges+=r['unresolved_edge_contacts']
                pair_residual=max(pair_residual,float(np.linalg.norm(r['paired_impulse_residual_ns'])))
                pair_work_error=max(pair_work_error,abs(r['kinetic_transfer_a_j']+r['kinetic_transfer_b_j']+r['dissipation_j']))
                ia,wa,wb=(r['impulse_a_ns'],r['kinetic_transfer_a_j'],r['kinetic_transfer_b_j']) if pair.a_id=='tissue' else (r['impulse_b_ns'],r['kinetic_transfer_b_j'],r['kinetic_transfer_a_j'])
                impulse=ia.sum(axis=0);tissue_impulse+=impulse;absolute_tissue_impulse+=float(np.linalg.norm(impulse));tissue_work+=wa;panel_work+=wb
                interval_tissue_impulse+=ia
                interval_panel_impulse+=r['impulse_b_ns'] if pair.a_id=='tissue' else r['impulse_a_ns']
            if (step+1)%round(.005/dt_s)==0:
                penetration,unchecked=_penetration(tissue,panel,front)
                frame={'time_s':tissue.time_s,'glans_displacement_m':np.mean(tissue.position_m[glans]-x[glans],axis=0).tolist(),
                    'minimum_jacobian':tissue.minimum_jacobian(),'sampled_front_tissue_penetration_m':penetration,
                    'unresolved_nearest_edge_samples':unchecked,'numerical_energy_defect_j':defect,'contact_resolutions':count,
                    'contact_dissipation_j':loss,'tissue_interface_work_j':tissue_work,'panel_interface_work_j':panel_work,
                    'external_work_j':external_work,
                    'tissue_contact_impulse_ns':tissue_impulse.tolist(),'elastic_energy_j':result['elastic_energy_j'],'kinetic_energy_j':result['kinetic_energy_j']}
                frames.append(frame);frame_log.write(json.dumps(frame,allow_nan=False)+'\n');frame_log.flush()
                tissue_frames.append(tissue.position_m.copy());panel_frames.append(panel.position_m.copy());times.append(tissue.time_s)
                tissue_velocities.append(tissue.velocity_m_s.copy());panel_velocities.append(panel.velocity_m_s.copy())
                tissue_contact_force.append(interval_tissue_impulse/.005);panel_contact_force.append(interval_panel_impulse/.005)
                interval_tissue_impulse.fill(0.);interval_panel_impulse.fill(0.)
    final_penetration,unchecked=_penetration(tissue,panel,front)
    report={'status':'completed_exploratory_fixture','seconds':tissue.time_s,'wall_seconds':time.monotonic()-elapsed,
        'minimum_jacobian':minimum_j,'maximum_momentum_residual_ns':momentum,'maximum_pair_impulse_residual_ns':pair_residual,
        'maximum_pair_work_balance_error_j':pair_work_error,'contact_resolutions':count,'unresolved_edge_contact_candidates':edges,
        'tissue_contact_impulse_ns':tissue_impulse.tolist(),'tissue_contact_impulse_norm_ns':absolute_tissue_impulse,
        'tissue_interface_work_j':tissue_work,'panel_interface_work_j':panel_work,'contact_dissipation_j':loss,
        'numerical_energy_defect_j':defect,'final_glans_displacement_m':np.mean(tissue.position_m[glans]-x[glans],axis=0).tolist(),
        'initial_elastic_energy_j':initial_energy,'external_work_j':external_work,'final_total_energy_j':result['total_energy_j'],
        'energy_audit_closure_j':result['total_energy_j']-initial_energy-external_work+loss-defect,
        'final_sampled_front_penetration_m':final_penetration,'final_unresolved_edge_samples':unchecked,
        'whole_garment_containment_validated':False,'biological_calibration':False,
        'limitations':['Only an extracted generated-short front panel with fixed far field participates.','No tunica, skin/fascial shell or pressure-dependent porous erectile tissue law.','No full pelvis/thigh contact, cloth bending or continuous/self collision.','Nodal face samples do not certify absence of all triangle/edge intersections.','Effective tissue inertia follows the prior canonical mass allocation; density reconciliation remains outstanding.']}
    np.savez_compressed(out/'trajectory.npz',time_s=np.array(times),tissue_positions_m=np.array(tissue_frames),panel_positions_m=np.array(panel_frames),
        tissue_velocity_m_s=np.array(tissue_velocities),panel_velocity_m_s=np.array(panel_velocities),
        tissue_contact_force_n=np.array(tissue_contact_force),panel_contact_force_n=np.array(panel_contact_force),
        tissue_tetrahedra=tetrahedra,panel_triangles=panel.triangles,tissue_surface_triangles=tissue.triangles)
    report['artifacts_sha256']={p.name:_sha(p) for p in out.iterdir() if p.is_file()}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report
