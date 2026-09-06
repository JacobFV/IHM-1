"""Small native initialization acceptance, explicitly gated by --run-native."""
from pathlib import Path
import argparse,copy,json,os,subprocess,sys,tempfile,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.source_skin import file_sha256
from ihm.assembly.hair_dynamics import ElasticHairState
from materialize_hair_residual_native import materialize


def rigid_receipt(properties,state):
    transform=np.asarray(state['transform_ground']);rotation=transform[:3,:3];offset=rotation@np.asarray(properties['centroid_m']);center=transform[:3,3]+offset
    omega=np.asarray(state['angular_velocity_rad_s']);v=np.asarray(state['origin_velocity_m_s'])+np.cross(omega,offset);mass=properties['mass_kg'];inertia=rotation@np.asarray(properties['inertia_com_kg_m2'])@rotation.T
    momentum=mass*v
    return momentum,inertia@omega+np.cross(center,momentum),.5*mass*(v@v)+.5*omega@inertia@omega


def validate(snapshot,identity,base,partition):
    if snapshot.get('mass_transfer',{}).get('enabled',False):raise ValueError('Intake mass transfer port unexpectedly enabled')
    if snapshot['mass_scale']!=1. or set(snapshot['bodies'])!=set(identity['expected_native_bodies']):raise ValueError('Native scaled or changed the exact residual body set')
    if abs(snapshot['mass_kg']-identity['target_native_mass_kg'])>1e-12:raise ValueError('Native residual mass mismatch')
    originals={b['owner']:b for b in base['actual_frozen_native_inertia']['bodies']};before_p=np.zeros(3);before_l=np.zeros(3);before_e=0.;after_p=np.zeros(3);after_l=np.zeros(3);after_e=0.
    for owner,expected in identity['expected_native_bodies'].items():
        state=snapshot['bodies'][owner];moments=state['inertia_moments_kg_m2'];products=state['inertia_products_kg_m2'];tensor=np.diag(moments)
        for (i,j),v in zip(((0,1),(0,2),(1,2)),products):tensor[i,j]=tensor[j,i]=v
        actual={'mass_kg':state['mass_kg'],'centroid_m':state['mass_center_local_m'],'inertia_com_kg_m2':tensor.tolist()}
        for key in actual:
            if not np.allclose(actual[key],expected[key],rtol=0,atol=2e-14):raise ValueError('Loaded native residual tensor mismatch: '+owner+'/'+key)
        p,l,e=rigid_receipt(originals[owner],state);before_p+=p;before_l+=l;before_e+=e
        p,l,e=rigid_receipt(actual,state);after_p+=p;after_l+=l;after_e+=e
    native_e=after_e
    if abs(snapshot['kinetic_energy_j']-native_e)>2e-11:raise ValueError('Native kinetic energy disagrees with loaded body tensors')
    hair_states={};hair_p=np.zeros(3);hair_l=np.zeros(3);hair_e=0.
    for name,group in partition['groups'].items():
        prepared=group['prepared'];data=copy.deepcopy(prepared['strands']);x=np.asarray(data['centerlines_m']).reshape(-1,3);offsets=data['strand_offsets'];positions=np.zeros_like(x);velocities=np.zeros_like(x)
        for i,(a,b) in enumerate(zip(offsets[:-1],offsets[1:])):
            owner=prepared['guide_owner_names'][i];state=snapshot['bodies'][owner];embedding=np.asarray(originals[owner]['canonical_reference_to_body_local']);transform=np.asarray(state['transform_ground']);local=x[a:b]@embedding[:3,:3].T+embedding[:3,3]
            positions[a:b]=local@transform[:3,:3].T+transform[:3,3]
            velocities[a:b]=np.asarray(state['origin_velocity_m_s'])+np.cross(np.asarray(state['angular_velocity_rad_s']),positions[a:b]-transform[:3,3])
        data['centerlines_m']=positions.tolist();hair=ElasticHairState(data);hair.velocity_m_s[:]=velocities;hair.time_s=snapshot['time_s']
        p=(hair.mass_kg[:,None]*velocities).sum(0);l=np.cross(positions,hair.mass_kg[:,None]*velocities).sum(0);e=.5*float(np.sum(hair.mass_kg[:,None]*velocities**2))
        if abs(hair.energy_j()-e)>1e-18:raise ValueError('Undeclared initial beam deformation energy')
        hair_p+=p;hair_l+=l;hair_e+=e;hair_states[name]={'time_s':hair.time_s,'positions_native_ground_m':positions.tolist(),'velocities_native_ground_m_s':velocities.tolist(),'nodal_mass_kg':hair.mass_kg.tolist(),'sample_ids':prepared['sample_ids'],'owner_names':prepared['guide_owner_names'],'deformation':'Zero initially straight stress-free beam; rigid parent station velocity at every node'}
    dp=after_p+hair_p-before_p;dl=after_l+hair_l-before_l;de=after_e+hair_e-before_e
    if np.linalg.norm(dp)>2e-12 or np.linalg.norm(dl)>2e-12 or abs(de)>2e-12:raise ValueError('Native residual plus initial hair momentum/energy mismatch')
    return {'loaded_body_count':len(originals),'mass_scale':snapshot['mass_scale'],'represented_guides':94,'native_residual_mass_kg':snapshot['mass_kg'],'combined_mass_kg':snapshot['mass_kg']+identity['represented_hair_mass_kg'],
        'linear_momentum_residual_kg_m_s':dp.tolist(),'angular_momentum_residual_kg_m2_s':dl.tolist(),'kinetic_energy_residual_j':float(de),'native_kinetic_energy_j':native_e,'hair_kinetic_energy_j':hair_e,
        'initial_hair_states':hair_states,'native_interval_integrated':False,'common_interval_ready':False,'remaining':'Needs composed single-interval hair/body owner, adequate contact coverage and checkpoint/energy acceptance; initialization does not certify a coupled trajectory'}


def light():
    base=json.loads((ROOT/'data/research/hair_source_factory.json').read_bytes());partition=json.loads((ROOT/'data/research/hair_native_partition.json').read_bytes())
    with tempfile.TemporaryDirectory(prefix='hair-residual-light-',dir=ROOT/'data/derived') as tmp:
        identity=materialize(Path(tmp)/'inputs');bodies={};energy=0.
        for original in base['actual_frozen_native_inertia']['bodies']:
            owner=original['owner'];p=identity['expected_native_bodies'][owner];inertia=np.asarray(p['inertia_com_kg_m2']);state={'transform_ground':original['native_reference_transform'],'origin_velocity_m_s':[.013,.002,-.003],'angular_velocity_rad_s':[.004,-.006,.007]}
            energy+=rigid_receipt(p,state)[2];bodies[owner]={**state,'mass_kg':p['mass_kg'],'mass_center_local_m':p['centroid_m'],'inertia_moments_kg_m2':np.diag(inertia).tolist(),'inertia_products_kg_m2':[inertia[0,1],inertia[0,2],inertia[1,2]]}
        snapshot={'time_s':0.,'mass_scale':1.,'mass_kg':identity['target_native_mass_kg'],'kinetic_energy_j':energy,'bodies':bodies};receipt=validate(snapshot,identity,base,partition)
        assert receipt['hair_kinetic_energy_j']>0
        for mutate in ('scale','tensor'):
            bad=copy.deepcopy(snapshot)
            if mutate=='scale':bad['mass_scale']=1.0000001
            else:bad['bodies']['torso']['inertia_products_kg_m2'][0]+=1e-6
            try:validate(bad,identity,base,partition)
            except ValueError:pass
            else:raise AssertionError('Hidden rescale/tensor mismatch accepted')
    print('PASS isolated XML residual construction;22body mass/COM/full tensor checks; nonzero per-node rigid velocity momentum/angular momentum/energy closure; hidden rescale/tensor rejection; no native run')


def run_native():
    started=time.monotonic();output=Path(tempfile.mkdtemp(prefix='hair-residual-native-',dir=ROOT/'data/derived'));identity=materialize(output/'inputs');runtime=ROOT/'data/runtime/opensim'
    pointer=json.loads((ROOT/'data/runtime/mechanical-stream/latest.json').read_bytes());build=ROOT/pointer['build'];manifest=json.loads((build/'manifest.json').read_bytes());executable=build/'native_mechanical_stream'
    if not executable.is_file():raise ValueError('Expected existing native binary absent')
    expected=manifest['files'].get(str(executable.relative_to(ROOT)))
    if expected is None or file_sha256(executable)!=expected:raise ValueError('Existing native binary identity mismatch')
    libdirs=[runtime/'install/opensim/lib',runtime/'install/simbody/lib',*[runtime/'sysroot/usr/lib/aarch64-linux-gnu'/x for x in ('lapack','blas','')]]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',IHM_INSTANCE_MASS_MODE='0',LD_LIBRARY_PATH=':'.join(map(str,libdirs)))
    run=output/'native';run.mkdir();command=['prlimit','--as=1073741824','--cpu=30','--','nice','-n','10','taskset','-c','0',str(executable),str(output/'inputs'),str(run),'free',format(identity['target_native_mass_kg'],'.17g')]
    completed=subprocess.run(command,input='close\n',capture_output=True,text=True,env=env,cwd=run,timeout=45);(output/'stdout.log').write_text(completed.stdout);(output/'stderr.log').write_text(completed.stderr)
    if completed.returncode:raise RuntimeError('Native initial fixture failed; see '+str(output))
    rows=[json.loads(line[5:]) for line in completed.stdout.splitlines() if line.startswith('@IHM ')]
    if len(rows)!=1 or 'error' in rows[0]:raise ValueError('Expected exactly one initial native snapshot')
    snapshot=rows[0]
    if not any(np.linalg.norm(b['angular_velocity_rad_s'])>1e-6 for b in snapshot['bodies'].values()):raise ValueError('Nonzero rotational initialization missing')
    if snapshot['mass_transfer']['enabled'] or snapshot['mass_transfer']['last_sequence']!=0:raise ValueError('Intake mass port was used')
    base=json.loads((ROOT/'data/research/hair_source_factory.json').read_bytes());partition=json.loads((ROOT/'data/research/hair_native_partition.json').read_bytes());receipt=validate(snapshot,identity,base,partition)
    loaded=subprocess.run(['ldd',str(executable)],capture_output=True,text=True,env=env,timeout=5)
    libraries={p:file_sha256(Path(p)) for p in loaded.stdout.split() if p.startswith('/') and Path(p).is_file()}
    receipt.update(loaded_libraries=libraries,acceptance_script_sha256=file_sha256(Path(__file__)),schema='ihm.hair-residual-native-acceptance.v1',native_initialized=True,live_enabled=False,output=str(output.relative_to(ROOT)),wall_s=time.monotonic()-started,command=command,native_binary_sha256=expected,build_manifest_sha256=file_sha256(build/'manifest.json'),input_identity=identity)
    (output/'snapshot.json').write_text(json.dumps(snapshot,indent=2,allow_nan=False)+'\n');(output/'receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:receipt[k] for k in ('output','loaded_body_count','mass_scale','combined_mass_kg','kinetic_energy_residual_j','wall_s')}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run-native',action='store_true');args=parser.parse_args()
    run_native() if args.run_native else light()
