"""Gated short coupled candidate fixture; no shared stream/runtime edits."""
from pathlib import Path
import argparse,copy,json,os,signal,selectors,subprocess,sys,tempfile,threading,time,uuid
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.assembly.source_skin import file_sha256
from ihm.assembly.hair_coupled_candidate import coupled_candidate
from scripts.materialize_hair_residual_native import materialize
from scripts.verify_hair_residual_native import validate,rigid_receipt


class CandidateStream(NativeMechanicalStream):
    def __init__(self,output,identity):
        self.output=output;output.mkdir();self.lock=threading.RLock();self.closed=False;self.tokens=set();self.identity=uuid.uuid4().hex;self.surface_sensor_identity={};self._buffer=b''
        runtime=ROOT/'data/runtime/opensim';pointer=json.loads((ROOT/'data/runtime/mechanical-stream/latest.json').read_bytes());build=ROOT/pointer['build'];manifest=json.loads((build/'manifest.json').read_bytes());exe=build/'native_mechanical_stream'
        if file_sha256(exe)!=manifest['files'][str(exe.relative_to(ROOT))]:raise ValueError('Native binary hash mismatch')
        self.binary_sha256=file_sha256(exe);self.log=(output/'engine.log').open('w')
        paths=[runtime/'install/opensim/lib',runtime/'install/simbody/lib',*[runtime/'sysroot/usr/lib/aarch64-linux-gnu'/x for x in ('lapack','blas','')]]
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',IHM_INSTANCE_MASS_MODE='0',LD_LIBRARY_PATH=':'.join(map(str,paths)))
        command=['prlimit','--as=1073741824','--cpu=30','--','nice','-n','10','taskset','-c','0',str(exe),str(ROOT/identity['destination']),str(output),'free',format(identity['target_native_mass_kg'],'.17g')]
        self.process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,env=env,cwd=output,bufsize=0);self._selector=selectors.DefaultSelector();self._selector.register(self.process.stdout,selectors.EVENT_READ)
        try:self.state=self._read()
        except BaseException:self.close();raise


def totals(snapshot,owner):
    p=np.zeros(3);l=np.zeros(3);energy=snapshot['kinetic_energy_j']+snapshot['potential_energy_j'];reg=owner.registration
    for b in snapshot['bodies'].values():
        tensor=np.diag(b['inertia_moments_kg_m2'])
        for (i,j),v in zip(((0,1),(0,2),(1,2)),b['inertia_products_kg_m2']):tensor[i,j]=tensor[j,i]=v
        bp,bl,_=rigid_receipt({'mass_kg':b['mass_kg'],'centroid_m':b['mass_center_local_m'],'inertia_com_kg_m2':tensor},b);p+=bp;l+=bl
    h=owner.hair;x=(h.position_m-reg.global_map[:3,3])@reg.basis;v=h.velocity_m_s@reg.basis;p+=(h.mass_kg[:,None]*v).sum(0);l+=np.cross(x,h.mass_kg[:,None]*v).sum(0)
    return {'momentum':p,'angular':l,'energy':energy+h.energy_j()}


def main(run=False):
    base=json.loads((ROOT/'data/research/hair_source_factory.json').read_bytes());partition=json.loads((ROOT/'data/research/hair_native_partition.json').read_bytes());registration=json.loads((ROOT/base['actual_frozen_native_inertia']['registration_path']).read_bytes())
    if not run:
        receipt=json.loads((ROOT/'data/research/hair_residual_native_acceptance.json').read_bytes());snapshot=json.loads((ROOT/receipt['output']/'snapshot.json').read_bytes());owner=coupled_candidate(partition,base,registration,snapshot)
        assert len(owner.hair.rods)==94 and len(owner.reference)==278 and owner.clamp_residual(snapshot)<1e-10
        saved=owner.checkpoint();owner.hair.position_m[-1,0]+=.001;owner.restore(saved);assert np.array_equal(owner.hair.position_m,saved['hair']['position_m'])
        print('PASS one94guide owner, exact278source nodes, initial moving material clamp, nonzero rigid velocities and checkpoint identity; no native run');return
    def timed_out(signum,frame):raise TimeoutError('Bounded45s coupled native acceptance timeout')
    signal.signal(signal.SIGALRM,timed_out);signal.alarm(45)
    output=Path(tempfile.mkdtemp(prefix='hair-coupled-native-',dir=ROOT/'data/derived'));identity=materialize(output/'inputs');native=CandidateStream(output/'native',identity);rows=[];started=time.monotonic()
    try:
        validate(native.snapshot(),identity,base,partition);owner=coupled_candidate(partition,base,registration,native.snapshot());initial=totals(native.snapshot(),owner)
        for interval in range(2):
            before=native.snapshot();result,report=owner.advance(native,1e-4,[],{},tolerance_m=1e-10,max_iterations=8);current=totals(result,owner)
            assert np.linalg.norm(report['body_impulse_ns'])>0
            assert np.linalg.norm(report['hair_body_impulse_residual_ns'])<1e-11
            assert abs((result['external_work_j']-before['external_work_j'])-report['native_hair_work_j'])<1e-12
            assert abs(result['time_s']-(interval+1)*1e-4)<1e-14 and abs(owner.hair.time_s-result['time_s'])<1e-14
            assert np.linalg.norm(current['momentum']-initial['momentum'])<1e-8 and np.linalg.norm(current['angular']-initial['angular'])<1e-8
            rows.append({'interval':interval,'native_time_s':result['time_s'],'hair_time_s':owner.hair.time_s,'linear_momentum_change':(current['momentum']-initial['momentum']).tolist(),'angular_momentum_change':(current['angular']-initial['angular']).tolist(),'total_mechanical_energy_change_j':current['energy']-initial['energy'],'native_external_work_j':result['external_work_j'],'native_signed_active_fiber_work_j':result['signed_active_fiber_work_j'],'native_muscle_heat_energy_j':result['muscle_heat_energy_j'],'hair_audit':owner.frame()['audit']})
        token=native.checkpoint();saved=owner.checkpoint();owner.advance(native,1e-4,[],{},tolerance_m=1e-10,max_iterations=8);expected=owner.hair.checkpoint();expected_native=native.snapshot();native.restore(token);owner.restore(saved)
        owner.advance(native,1e-4,[],{},tolerance_m=1e-10,max_iterations=8)
        assert np.array_equal(owner.hair.position_m,expected['position_m']) and np.array_equal(owner.hair.velocity_m_s,expected['velocity_m_s'])
        assert native.snapshot()['time_s']==expected_native['time_s'] and native.snapshot()['bodies']==expected_native['bodies']
        for key in ('external_work_j','signed_active_fiber_work_j','muscle_heat_energy_j','muscle_metabolic_energy_j'):assert native.snapshot()[key]==expected_native[key]
        native.restore(token);owner.restore(saved);before=native.snapshot()
        try:owner.advance(native,1e-4,[],{},tolerance_m=1e-30,max_iterations=1)
        except RuntimeError:pass
        else:raise AssertionError('Forced rejected interval accepted')
        assert native.snapshot()['time_s']==before['time_s'] and np.array_equal(owner.hair.position_m,saved['hair']['position_m'])
        assert native.snapshot()['bodies']==before['bodies']
        for key in ('external_work_j','signed_active_fiber_work_j','muscle_heat_energy_j','muscle_metabolic_energy_j'):assert native.snapshot()[key]==before[key]
        native.release(token);assert not native.tokens
        report={'schema':'ihm.hair-coupled-native-acceptance.v1','live_enabled':False,'output':str(output.relative_to(ROOT)),'binary_sha256':native.binary_sha256,'intervals':rows,'replay_exact_hair':True,'nonconvergence_atomic_rollback':True,'final_time_s':native.snapshot()['time_s'],'wall_s':time.monotonic()-started,'limitations':['Sparse94follicle faces only; no complete containment, edges,CCD,finite radius,self-contact or garment','Native muscles/passive tissue remain active; reported total mechanical change is not a closed whole-native energy certificate','Current22body reduction; scalp structural torso prior; no viewer synchronization']}
        (output/'receipt.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({'output':str(output),'final_time_s':report['final_time_s'],'wall_s':report['wall_s']}))
    finally:native.close();signal.alarm(0)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run-native',action='store_true');args=parser.parse_args();main(args.run_native)
