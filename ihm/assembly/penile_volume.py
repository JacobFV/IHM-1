"""Separate source-volume CC/CS experiment with human ex vivo material fits.

No glans coefficient, tunica architecture or whole-body material handoff is
inferred. Old material/contact experiment artifacts remain unchanged.
"""
from pathlib import Path
import hashlib,json,math,shutil,time
import numpy as np
from .contact_dynamics import DynamicTetrahedra,step_coupled
from ..calibration.penile import penile_material

BASE=Path(__file__).resolve().parents[2]
TISSUES=('corpus_cavernosum','corpus_spongiosum')
SOURCE_IDS=('body-bp3d-FJ3132','body-bp3d-FJ3133')

def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

class PenileTetrahedra(DynamicTetrahedra):
    """Exact reference-volume energy/force assembly with explicit tissue labels.

    max_explicit_dt_s uses the initial isotropic tangent, not a guaranteed
    nonlinear stability bound. Each accepted step still requires positive J.
    """
    def __init__(self,vertices_m,tetrahedra,*,material_index,material_tissues,density_kg_m3,
                 fixed_nodes=(),material_owner_ids=()):
        labels=np.asarray(material_index);tissues=tuple(material_tissues)
        if not tissues or len(set(tissues))!=len(tissues) or any(t not in TISSUES for t in tissues):
            raise ValueError('Explicit supported CC/CS tissue laws required')
        if labels.shape!=(len(tetrahedra),) or labels.dtype.kind not in 'iu' or labels.min()<0 or labels.max()>=len(tissues):
            raise ValueError('Every tetrahedron must have a valid material index')
        if set(labels)!=set(range(len(tissues))):raise ValueError('Every declared tissue must own tetrahedra')
        self.material_index=labels.copy();self.material_tissues=tissues
        self.laws=tuple(penile_material(t) for t in tissues)
        moduli=[law.initial_moduli() for law in self.laws]
        shear=np.array([m['shear_pa'] for m in moduli])[labels]
        bulk=np.array([m['bulk_pa'] for m in moduli])[labels]
        super().__init__(vertices_m,tetrahedra,mu_pa=shear,lambda_pa=bulk-2*shear/3,
                         density_kg_m3=density_kg_m3,fixed_nodes=fixed_nodes,material_owner_ids=material_owner_ids)

    def elastic_forces(self,positions_m):
        F=self.region.deformation(positions_m);density=np.empty(len(F));piola=np.empty_like(F)
        for index,law in enumerate(self.laws):
            mask=self.material_index==index
            density[mask],piola[mask]=law.energy_piola(F[mask])
        h=self.region.volumes[:,None,None]*piola@np.swapaxes(self.region.inverse,1,2)
        local=np.concatenate((-h.sum(axis=2)[:,None,:],np.swapaxes(h,1,2)),axis=1)
        force=np.zeros_like(self.reference_position_m)
        for corner in range(4):np.add.at(force,self.region.tets[:,corner],-local[:,corner])
        return force,float(self.region.volumes@density)


def source_cc_cs_domain(parent=None):
    """Select exact retained CC/CS cells; do not repopulate removed glans cells."""
    directory=Path(parent or BASE/'data/derived/material-domains/pelvis-0.004m')
    manifest=json.loads((directory/'manifest.json').read_text())
    for name,digest in manifest['artifacts'].items():
        if _sha(directory/name)!=digest:raise ValueError('Parent material domain artifact changed')
    ids=tuple(m['source_id'] for m in manifest['material_regions'])
    if ids[:2]!=SOURCE_IDS:raise ValueError('Parent CC/CS material identities changed')
    for surface in manifest['source_surfaces']:
        if surface['id'] in SOURCE_IDS and _sha(BASE/surface['path'])!=surface['source_sha256']:
            raise ValueError('Source anatomy geometry changed')
    with np.load(directory/'pelvic-domain.npz',allow_pickle=False) as raw:
        selected=np.flatnonzero(np.isin(raw['material_index'],[0,1]))
        nodes,reverse=np.unique(raw['tetrahedra'][selected],return_inverse=True)
        fixed=np.flatnonzero(np.isin(nodes,raw['fixed_nodes']))
        arrays={'vertices_m':raw['vertices_m'][nodes], 'tetrahedra':reverse.reshape(-1,4),
                'material_index':raw['material_index'][selected], 'fixed_nodes':fixed,
                'density_kg_m3':raw['density_kg_m3'][selected],
                'generic_mu_pa':raw['mu_pa'][selected],'generic_lambda_pa':raw['lambda_pa'][selected],
                'parent_node_indices':nodes,'parent_tetrahedron_indices':selected}
    if len(fixed)<3:raise ValueError('Retained source domain lacks its prescribed support')
    return arrays,manifest,directory


def run_penile_volume(output_dir,*,law='published',dt_s=.00005,seconds=.03,parent=None):
    """Matched pulse/release fixture, actual finite-element states and work."""
    if law not in ('published','generic_prior'):raise ValueError('Unknown constitutive comparison')
    for value in (dt_s,seconds):
        if isinstance(value,bool) or not np.isfinite(value) or value<=0:raise ValueError('Positive finite times required')
    if dt_s>.00005 or seconds!=.03 or abs(.001/dt_s-round(.001/dt_s))>1e-8:
        raise ValueError('Use a 30 ms fixture and a step that divides its 1 ms sample interval')
    out=Path(output_dir)
    if out.exists():raise ValueError('Choose a new volume experiment directory')
    arrays,parent_manifest,directory=source_cc_cs_domain(parent)
    common=dict(density_kg_m3=arrays['density_kg_m3'],fixed_nodes=arrays['fixed_nodes'],material_owner_ids=SOURCE_IDS)
    x=arrays['vertices_m'];t=arrays['tetrahedra'];labels=arrays['material_index']
    body=(PenileTetrahedra(x,t,material_index=labels,material_tissues=TISSUES,**common) if law=='published' else
          DynamicTetrahedra(x,t,mu_pa=arrays['generic_mu_pa'],lambda_pa=arrays['generic_lambda_pa'],**common))
    if dt_s>body.max_explicit_dt_s:raise ValueError('Step exceeds initial-tangent restriction')
    load_nodes=np.setdiff1d(body.surface_nodes[x[body.surface_nodes,2]>=np.quantile(x[body.surface_nodes,2],.9)],body.fixed_nodes)
    if not len(load_nodes):raise ValueError('Distal source-surface loading region is empty')
    out.mkdir(parents=True)
    sources=[Path(__file__),BASE/'ihm/calibration/penile.py',BASE/'ihm/assembly/contact_dynamics.py',
             BASE/'ihm/assembly/mechanics_backend.py',BASE/'data/measurements/biomechanics/khorshidi_2024.json',
             BASE/'data/raw/biomechanics/human-penile-mechanics-2024/paper.pdf',directory/'manifest.json',directory/'pelvic-domain.npz']
    sources += [BASE/s['path'] for s in parent_manifest['source_surfaces'] if s['id'] in SOURCE_IDS]
    hashes={str(p.relative_to(BASE)):_sha(p) for p in sources}
    for path in sources:
        target=out/'inputs'/path.relative_to(BASE);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
        if _sha(target)!=hashes[str(path.relative_to(BASE))]:raise ValueError('Source changed during experiment capture')
    np.savez_compressed(out/'domain.npz',**arrays,load_nodes=load_nodes,boundary_triangles=body.triangles)
    configuration={'law':law,'dt_s':dt_s,'seconds':seconds,'sample_interval_s':.001,
                   'source_sha256':hashes,'material_tissues':TISSUES,'material_owner_ids':SOURCE_IDS,
                   'constitutive_evidence':[penile_material(tissue).describe() for tissue in TISSUES],
                   'frame':parent_manifest['frame'],'parent_domain':str(directory.relative_to(BASE)),
                   'spacing_m':parent_manifest['spacing_m'],'boundary_discretization_diagonal_m':parent_manifest['boundary_discretization_diagonal_m'],
                   'initial_tangent_dt_limit_s':body.max_explicit_dt_s,
                   'load':{'peak_force_n':[0,.02,0],'duration_s':.015,'profile':'sin(pi*t/duration) until duration, then zero',
                           'spatial_distribution':'equal nodal force over top anterior surface-position decile excluding fixed nodes'},
                   'support_basis':parent_manifest['support_basis'],'density_basis':'unchanged parent effective canonical mass allocation, not measured ex vivo density',
                   'glans_cells':'excluded; parent overlap priority preserved with no reassignment',
                   'whole_body_handoff_applied':False,'biological_validation':False,
                   'limitations':['No glans, tunica/fiber architecture, fascia/skin shell, resolved urethral lumen, pressure, perfusion or contact.',
                                  'CC/CS point estimates from three elderly fresh-frozen donors; no joint parameter covariance.',
                                  'D2=0 is interpreted as omitted higher-order volume term; original Abaqus runtime parity unresolved.',
                                  'Shared nodes bond interfaces; support and prescribed loading are engineering priors.',
                                  'Initial tangent timestep restriction does not bound nonlinear strain stiffening.',
                                  'Space discretization and physiological validity are not established by timestep refinement.']}
    (out/'configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    positions=[x.copy()];velocities=[body.velocity_m_s.copy()];times=[0.];frames=[]
    minimum_j=1.;work=0.;defect=0.;max_momentum=0.;initial=body.elastic_forces(x)[1];started=time.monotonic()
    interval_support=np.zeros(3);interval_load=np.zeros(3)
    with (out/'frames.jsonl').open('x') as stream:
        for step in range(round(seconds/dt_s)):
            at=step*dt_s;force=np.zeros_like(x)
            if at<.015:force[load_nodes,1]=.02*math.sin(math.pi*at/.015)/len(load_nodes)
            r=step_coupled({'tissue':body},dt_s,external_forces_n={'tissue':force},gravity_m_s2=(0,0,0))
            minimum_j=min(minimum_j,body.minimum_jacobian());work+=r['external_work_j'];defect+=r['numerical_energy_defect_j']
            interval_support+=r['support_impulse_ns'];interval_load+=force.sum(axis=0)*dt_s
            max_momentum=max(max_momentum,float(np.linalg.norm(r['momentum_residual_ns'])))
            if (step+1)%round(.001/dt_s)==0:
                record={key:r[key] for key in ('time_s','elastic_energy_j','kinetic_energy_j','total_energy_j')}
                record.update(minimum_jacobian=body.minimum_jacobian(),external_work_j=work,numerical_energy_defect_j=defect,
                              support_impulse_ns=interval_support.tolist(),mean_applied_force_n=(interval_load/.001).tolist(),
                              loaded_region_mean_displacement_m=(body.position_m[load_nodes]-x[load_nodes]).mean(axis=0).tolist())
                stream.write(json.dumps(record,allow_nan=False)+'\n');stream.flush();frames.append(record)
                positions.append(body.position_m.copy());velocities.append(body.velocity_m_s.copy());times.append((step+1)*dt_s)
                interval_support.fill(0.);interval_load.fill(0.)
    np.savez_compressed(out/'trajectory.npz',time_s=times,positions_m=positions,velocity_m_s=velocities)
    report={'status':'completed_separate_ex_vivo_law_fixture','law':law,'wall_seconds':time.monotonic()-started,
            'minimum_jacobian':minimum_j,'maximum_momentum_residual_ns':max_momentum,
            'external_work_j':work,'numerical_energy_defect_j':defect,'final_total_energy_j':r['total_energy_j'],
            'energy_audit_closure_j':r['total_energy_j']-initial-work-defect,
            'loaded_region_mean_displacement_m':frames[-1]['loaded_region_mean_displacement_m'],
            'vertices':len(x),'tetrahedra':len(t),'material_tetrahedron_counts':np.bincount(labels).tolist(),
            'mass_kg':float(body.mass_kg.sum()),'biological_validation':False}
    report['artifacts_sha256']={str(p.relative_to(out)):_sha(p) for p in out.rglob('*') if p.is_file()}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report
