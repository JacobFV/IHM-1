"""Opt-in endpoint mass contract: light tests or coordinated native acceptance."""
from pathlib import Path
import argparse,copy,hashlib,json,os,resource,shutil,subprocess,sys,tempfile,time,unittest,threading
from types import SimpleNamespace
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.instance_mass import variant_identity,pointer_directory
from ihm.native.mechanical_stream import NativeMechanicalStream,NativeCommandRejected
ROOT=Path(__file__).resolve().parents[1]
VARIANT='data/runtime/opensim/variants/instance_mass_v1'
STATION=[.046177955611313035,.12240375804724452,-.037921654229177336]
ENERGY=('muscle_metabolic_energy_j','signed_active_fiber_work_j','muscle_heat_energy_j','positive_active_fiber_work_j','external_work_j')
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def normalized(s):s=copy.deepcopy(s);s.pop('kind',None);return s
class SelectionTests(unittest.TestCase):
    def test_requires_owned_complete_exact_library(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);v=root/'data/runtime/opensim/variants/example';v.mkdir(parents=True);lib=v/'libSimTKsimbody.so.3.9';lib.write_bytes(b'fixture library');manifest={'schema':'ihm.simbody-instance-mass-variant.v1','complete':True,'library_sha256':hashlib.sha256(lib.read_bytes()).hexdigest()};dump(v/'manifest.json',manifest)
            identity=variant_identity(root,v);self.assertEqual(pointer_directory(root,identity).parent.name,'instance-mass')
            lib.write_bytes(b'changed')
            with self.assertRaises(ValueError):variant_identity(root,v)
            with self.assertRaises(ValueError):variant_identity(root,root/'outside')
    def test_default_rejects_before_native_request(self):
        stream=NativeMechanicalStream.__new__(NativeMechanicalStream);stream.instance_mass_variant=None
        with self.assertRaisesRegex(ValueError,'explicitly selected'):stream.transfer_mass(sequence=1,owner='a',delta_mass_kg=0,station_m=STATION,velocity_source_m_s=[0,0,0])
    def test_python_payload_validation(self):
        stream=NativeMechanicalStream.__new__(NativeMechanicalStream);stream.instance_mass_variant={}
        base=dict(sequence=1,owner='meal',delta_mass_kg=.1,station_m=STATION,velocity_source_m_s=[0,0,0])
        for change in ({'sequence':True},{'sequence':0},{'sequence':2**64},{'owner':'bad name'},{'delta_mass_kg':.51},{'delta_mass_kg':float('nan')},{'body':'pelvis'}):
            with self.subTest(change=change),self.assertRaises(ValueError):stream.transfer_mass(**dict(base,**change))
    def test_uncertain_transport_poisoned_but_native_rejection_recoverable(self):
        for mode in ('short','parse','rejected'):
            with self.subTest(mode=mode):
                stream=NativeMechanicalStream.__new__(NativeMechanicalStream);stream.lock=threading.RLock();stream.closed=False
                stream.process=SimpleNamespace(stdin=SimpleNamespace(write=lambda x:len(x)-(mode=='short'),flush=lambda:None))
                def read():
                    if mode=='rejected':raise NativeCommandRejected('rolled back')
                    raise ValueError('malformed reply')
                stream._read=read;stream.close=lambda:setattr(stream,'closed',True)
                with self.assertRaises((ValueError,OSError)):stream._request('mass_transfer fixture')
                self.assertEqual(stream.closed,mode!='rejected')
    def test_wrong_success_kind_poisoned(self):
        stream=NativeMechanicalStream.__new__(NativeMechanicalStream);stream.lock=threading.RLock();stream.closed=False;stream.identity='test';stream.instance_mass_variant={}
        stream._request=lambda _:dict(kind='restored',mass_transfer=dict(reference_id='test',enabled=True,last_sequence=1,last_receipt=dict(sequence=1,zero=True)))
        stream.close=lambda:setattr(stream,'closed',True)
        with self.assertRaises(ValueError):stream.transfer_mass(sequence=1,owner='meal',delta_mass_kg=0,station_m=STATION,velocity_source_m_s=[0,0,0])
        self.assertTrue(stream.closed)
def native():
    from ihm.assembly.articulated import ArticulatedBodyPlant
    from verify_opensim_instance_mass_native import compare,Session
    started=time.monotonic();out=Path(tempfile.mkdtemp(prefix='native-local-mass-port-',dir=ROOT/'data/derived'))
    source_names=('scripts/native_mechanical_stream.cpp','scripts/native_local_mass_port.h','scripts/native_muscle_metabolism.h','scripts/native_static_pose.h','scripts/native_surface_foundation.h','scripts/native_bed_compression.h','scripts/build_native_mechanical_stream.py','ihm/native/mechanical_stream.py','ihm/native/instance_mass.py','ihm/assembly/articulated.py')
    frozen={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in source_names}
    base=NativeMechanicalStream(ROOT,out/'baseline',target_mass_kg=77.6122029,augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json')
    try:
        baseline=base.snapshot();assert baseline['mass_transfer']['enabled'] is False
        try:base.transfer_mass(sequence=1,owner='meal',delta_mass_kg=0,station_m=STATION,velocity_source_m_s=[0,0,0]);raise AssertionError('default enabled')
        except ValueError:pass
        base_point=base.body_point(body='torso',station_m=STATION)
    finally:base.close()
    plant=ArticulatedBodyPlant(ROOT,out/'plant',target_mass_kg=77.6122029,instance_mass_variant=VARIANT,augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json')
    try:
        n=plant.native;initial=n.snapshot();a=normalized(baseline);b=normalized(initial);a.pop('mass_transfer');b.pop('mass_transfer');parity=compare(a,b)
        assert len(initial['muscles'])==92;compare(base_point,n.body_point(body='torso',station_m=STATION))
        def request(dm,owner='meal',velocity=None):
            sequence=n.snapshot()['mass_transfer']['last_sequence']+1
            if velocity is None:velocity=plant.body_point(body='torso',station_m=STATION)['velocity_source_m_s']
            return plant.transfer_mass(sequence=sequence,owner=owner,delta_mass_kg=dm,station_m=STATION,velocity_source_m_s=velocity)
        def rejected(**change):
            before=normalized(n.snapshot());args=dict(sequence=before['mass_transfer']['last_sequence']+1,owner='meal',delta_mass_kg=.01,station_m=STATION,velocity_source_m_s=[0,0,0]);args.update(change)
            try:n.transfer_mass(**args);raise AssertionError('bad command accepted')
            except ValueError:pass
            compare(before,normalized(n._request('observe')))
        request(0);zero=n.snapshot();assert zero['mass_transfer']['owners']=={};compare(initial['bodies'],zero['bodies'])
        plant.advance(.001);pre=n.snapshot();pre_frame=plant.snapshot();saved=plant.checkpoint();request(.25,velocity=[.02,-.01,.03]);quarter=n.snapshot();request(.25,velocity=[.02,-.01,.03]);loaded=n.snapshot()
        assert loaded['mass_transfer']['owned_payload_mass_kg']==.5
        assert abs(loaded['bodies']['torso']['mass_kg']-initial['bodies']['torso']['mass_kg']-.5)<1e-12
        assert loaded['bodies']['torso']['model_baseline_mass_properties']==initial['bodies']['torso']['model_baseline_mass_properties']
        for key in ENERGY:assert loaded[key]==pre[key]
        assert loaded['metabolic_reference']==pre['metabolic_reference'];assert plant.snapshot()['mass_transfer']==loaded['mass_transfer'];assert plant.snapshot()['positive_muscle_work_j']==pre_frame['positive_muscle_work_j']
        for args in ({'delta_mass_kg':.001},{'sequence':loaded['mass_transfer']['last_sequence']},{'owner':'unknown','delta_mass_kg':-.01},{'station_m':[0,0,0],'delta_mass_kg':-.01},{'delta_mass_kg':-.1,'velocity_source_m_s':[100,0,0]}):rejected(**args)
        continuation=n.checkpoint();n.advance(.001);future=n.snapshot();n.restore(continuation);n.advance(.001);compare(normalized(future),normalized(n.snapshot()));n.restore(continuation);n.release(continuation)
        request(-.25);half_out=n.snapshot();request(-.25);empty=n.snapshot();assert empty['mass_transfer']['owned_payload_mass_kg']==0
        compare(initial['bodies']['torso']['mass_kg'],empty['bodies']['torso']['mass_kg']);rejected(delta_mass_kg=-.001)
        plant.restore(saved);request(.25,velocity=[.02,-.01,.03]);compare(normalized(quarter),normalized(n.snapshot()));plant.release(saved)
        execution=json.loads((n.output/'execution.json').read_text());inputs=n.output/'inputs';dump(out/'states.json',{'initial':initial,'loaded':loaded,'half_out':half_out,'empty':empty})
    finally:plant.close()
    # Compile an isolated deliberately faulting copy. No test hook is exposed in production.
    build=out/'fault-build';build.mkdir();manifest=execution['build'];original_build=ROOT/manifest['path']
    for name in ('native_mechanical_stream.cpp','native_local_mass_port.h','native_muscle_metabolism.h','native_static_pose.h','native_surface_foundation.h','native_bed_compression.h'):shutil.copyfile(original_build/name,build/name)
    source=build/'native_mechanical_stream.cpp';source.write_text('#define IHM_MASS_TEST_AFTER_SET throw std::runtime_error("injected post-setter failure")\n'+source.read_text())
    command=[str(build/'native_mechanical_stream') if x==str(original_build/'native_mechanical_stream') else str(source) if x==str(original_build/'native_mechanical_stream.cpp') else x for x in manifest['command']];dump(build/'command.json',command)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    with (build/'compile.log').open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    runtime=ROOT/'data/runtime/opensim';env.update(IHM_INSTANCE_MASS_MODE='1',IHM_MASS_REFERENCE_ID='fault-fixture',LD_LIBRARY_PATH=':'.join(str(p) for p in [ROOT/VARIANT,runtime/'install/opensim/lib',runtime/'install/simbody/lib',runtime/'sysroot/usr/lib/aarch64-linux-gnu/lapack',runtime/'sysroot/usr/lib/aarch64-linux-gnu/blas',runtime/'sysroot/usr/lib/aarch64-linux-gnu']))
    fault=Session(build/'native_mechanical_stream',inputs,out/'fault-session',env)
    try:
        before=fault.request('advance .001 0 0')
        try:fault.request('mass_transfer fault-fixture 1 meal torso .25 '+' '.join(map(str,STATION))+' 0 0 0');raise AssertionError('missing injected failure')
        except ValueError as e:assert 'injected post-setter failure' in str(e)
        after=fault.request('observe');compare(normalized(before),normalized(after));assert after['mass_transfer']['last_sequence']==0
        fault.request('checkpoint replay');forward=fault.request('advance .001 0 0');fault.request('restore replay');compare(normalized(forward),normalized(fault.request('advance .001 0 0')))
    finally:fault.close()
    assert all(hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==v for n,v in frozen.items())
    report={'passed':True,'schema':'ihm.native-local-mass-port-acceptance.v1','output_dir':str(out.relative_to(ROOT)),'source_sha256':frozen,'zero_mode_max_abs_difference':parity,'muscles':92,'checks':['baseline disabled','zero parity','body-point query','two cumulative additions','effective vs baseline properties','total cap','duplicate sequence','wrong owner','changed station','non-comoving removal rejection','co-moving staged outflow','inventory underflow','nonzero checkpoint continuation replay','owner sequence replay','separate preserved muscle ledger','post-setter failure rollback and continued replay'],'variant':execution['instance_mass_variant'],'build':manifest['path'],'wall_s':time.monotonic()-started,'maximum_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    dump(out/'report.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run-native',action='store_true');args=p.parse_args()
    if args.run_native:native()
    else:unittest.main(argv=[sys.argv[0]])
