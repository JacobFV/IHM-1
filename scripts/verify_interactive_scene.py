"""Actual force response, ground friction, conservation and failed-step rollback."""
from pathlib import Path
import json
import argparse
import copy
import gzip
import hashlib
import tempfile
import types
from unittest.mock import patch
import numpy as np
from ihm.assembly.interactive_scene import Sphere,InteractiveScene,SceneSessions
import ihm.assembly.interactive_scene as scene_module


def light_verification():
    """Fault lifecycle with one fake mechanics owner; no whole-body trajectory."""
    class TinyMechanics:
        mutate_path=None
        def __init__(self,payload):
            self.payload=payload;self.ids=['test'];self.index={'test':0};self.mass=np.ones(1)
            self.x0=np.array([[.1,.2,.3]]);self.x=self.x0.copy();self.v=np.zeros_like(self.x);self.time=0.
            if self.mutate_path:self.mutate_path.write_text('{"changed_after_parse":true}')
        def step(self,dt,drivers):
            self.x[0,0]+=.00001;self.time+=dt
            return {'time_s':self.time,'entities':{'test':{'centroid_m':self.x[0].tolist()}},'audit':{}}
    with tempfile.TemporaryDirectory(prefix='scene-fault-tests-') as temporary:
        root=Path(temporary);model=root/'data/derived/canonical/mechanics.json';model.parent.mkdir(parents=True);raw=b'{"identity":"consumed-by-tiny-test"}';model.write_bytes(raw)
        receipts=[]
        for receipt in scene_module._IMPORT_SOURCES:
            destination=root/'ihm/assembly'/receipt['path'].name;destination.parent.mkdir(parents=True,exist_ok=True);destination.write_bytes(receipt['bytes'])
            receipts.append(dict(receipt,path=destination))
        with patch.object(scene_module,'BodyMechanics',TinyMechanics),patch.object(scene_module,'_IMPORT_SOURCES',tuple(receipts)):
            TinyMechanics.mutate_path=model
            manager=SceneSessions(root);initial=manager.create({});ident=initial['id'];scene=manager.sessions[ident]
            TinyMechanics.mutate_path=None
            assert (scene.output/'inputs/mechanics.json').read_bytes()==raw and scene.body.payload==json.loads(raw)
            assert initial['entities']['test']['centroid_m']==[.1,.2,.3]
            before=manager.current(ident);before['entities']['test']['centroid_m'][0]=999
            assert manager.current(ident)['entities']['test']['centroid_m'][0]==.1
            for sequence in [True,False,0.,-1,None]:
                try:scene.step({'sequence':sequence,'seconds':.002})
                except ValueError:pass
                else:raise AssertionError('Invalid sequence type accepted')
            ball=scene.objects['scene-ball'];badpoint=(ball.position+[1,0,0]).tolist()
            try:scene.step({'sequence':0,'seconds':.002,'forces':[{'id':'scene-ball','force_n':[1,0,0],'point_m':badpoint}]})
            except ValueError:pass
            else:raise AssertionError('Distant sphere force point accepted')
            files=list((scene.output/'events').glob('*.json.gz'));original={p.name:p.read_bytes() for p in files};position=scene.body.x.copy()
            for operation in ['fsync','link']:
                with patch.object(scene_module.os,operation,side_effect=OSError('injected '+operation+' failure')):
                    try:scene.step({'sequence':0,'seconds':.002})
                    except OSError:pass
                    else:raise AssertionError('Faulted event write accepted')
                assert scene.sequence==0 and scene.body.time==0 and np.array_equal(scene.body.x,position)
                assert {p.name:p.read_bytes() for p in (scene.output/'events').glob('*.json.gz')}==original
                assert not list((scene.output/'events').glob('.pending-*'))
            blocker=scene.output/'events/00000001-advance.json.gz';blocker.write_bytes(b'immutable-existing-record')
            try:scene.step({'sequence':0,'seconds':.002})
            except FileExistsError:pass
            else:raise AssertionError('Existing immutable event replaced')
            assert blocker.read_bytes()==b'immutable-existing-record' and scene.sequence==0;blocker.unlink()
            accepted=scene.step({'sequence':0,'seconds':.002,'forces':[{'id':'scene-ball','force_n':[1,0,0],'point_m':ball.position.tolist()}]})
            assert accepted['sequence']==1
            try:scene.step({'sequence':0,'seconds':.002})
            except ValueError:pass
            else:raise AssertionError('Stale scene sequence accepted')
            recovered=manager.current(ident);assert recovered['sequence']==1 and recovered['time_s']==.002
            # Surface ownership persists while a sphere rotates between ticks.
            ball=scene.objects['scene-ball'];ball.omega[:]=[0,0,3.];offset=np.array([ball.radius,0,0]);observed=[];original_step=Sphere.step
            def inspect_point(owner,dt,force,point,gravity,**options):
                observed.append(owner.rotation.T@(point-owner.position))
                return original_step(owner,dt,force,point,gravity,**options)
            with patch.object(Sphere,'step',inspect_point):
                moved=scene.step({'sequence':recovered['sequence'],'seconds':.02,'forces':[{'id':'scene-ball','force_n':[0,1,0],'point_m':(ball.position+offset).tolist()}]})
            assert moved['sequence']==2 and len(observed)==10 and np.allclose(observed,np.tile(offset,(10,1)),atol=1e-12)
            with patch.object(scene_module.os,'link',side_effect=OSError('injected close failure')):
                try:manager.command(ident,'close',{})
                except OSError:pass
                else:raise AssertionError('Faulted close accepted')
            assert ident in manager.sessions and not manager.current(ident)['closed']
            assert manager.command(ident,'close',{})['closed'] and ident not in manager.sessions
            assert manager.command(ident,'close',{})['closed']  # Retry after a lost close response.
            previous=None
            for index,p in enumerate(sorted((scene.output/'events').glob('*.json.gz'))):
                value=json.loads(gzip.decompress(p.read_bytes()));assert value['event_index']==index and value['previous_event_sha256']==previous
                assert p.stat().st_nlink==1
                previous=hashlib.sha256(p.read_bytes()).hexdigest()
            source=receipts[0]['path'];source.write_bytes(source.read_bytes()+b'\n')
            try:manager.create({})
            except ValueError:pass
            else:raise AssertionError('Hot-edited source accepted')
        # Compare actual loaded code to disk source without a server or body run.
        source=root/'probe.py';source.write_text('def f(): return 1\n');module=types.ModuleType('probe');module.__file__=str(source)
        exec(compile(source.read_bytes(),str(source),'exec'),module.__dict__);scene_module._loaded_source(module)
        source.write_text('def f(): return 2\n')
        try:scene_module._loaded_source(module)
        except ValueError:pass
        else:raise AssertionError('Loaded-code/source mismatch accepted')
    print(json.dumps({'passed':True,'scope':'Tiny fake-owner lifecycle and pure source identity tests; no whole-body trajectory',
                      'checks':['initial_centroids','detached_snapshot','strict_sequence','sphere_owner_point','rotating_material_point','fsync_rollback','atomic_publication_rollback','close_rollback','current_recovery','exact_consumed_bytes','event_chain','hot_source_rejection','loaded_code_source_match']},indent=2))


def main():
    ball=Sphere('test',[0.,1.,0.],radius=.05,mass=.4)
    initial=ball.kinetic()+ball.mass*9.81*ball.position[1]
    for _ in range(400):ball.step(.002,np.zeros(3),ball.position.copy(),np.array([0.,-9.81,0.]),axis=1,plane=0.)
    residual=ball.kinetic()+ball.mass*9.81*ball.position[1]-initial+ball.contact_loss-ball.projection_work
    assert abs(residual)<1e-10 and ball.position[1]>=ball.radius-1e-12
    assert ball.contact_loss>0
    sliding=Sphere('slide',[0.,.049,0.],radius=.05,mass=.4);sliding.velocity[:]=[1.,-1.,0.]
    before=sliding.kinetic();momentum=sliding.mass*sliding.velocity.copy()
    sliding.step(.002,np.zeros(3),sliding.position.copy(),np.zeros(3),axis=1,plane=0.)
    assert sliding.kinetic()<before and np.linalg.norm(sliding.omega)>0
    assert np.allclose(sliding.mass*sliding.velocity+sliding.support_impulse,momentum,atol=1e-12)
    spin=Sphere('spin',[0.,0.,0.]);force=np.array([1.,0.,0.]);point=np.array([0.,spin.radius,0.])
    spin.step(.002,force,point,np.zeros(3),axis=1,plane=-1,contact=False)
    assert np.isclose(spin.inertia*spin.omega[2],-spin.radius*.002)
    root=Path(__file__).resolve().parents[1]
    parent=Path(tempfile.mkdtemp(prefix='interactive-scene-',dir=root/'data/derived/audits'))
    scene=InteractiveScene(root,parent/'studio')
    target='body-bp3d-FJ3393';force=[1.,0.,0.];point=scene.body.x[scene.body.index[target]].tolist()
    assert scene.snapshot()['entities'][target]['centroid_m']==point
    output=scene.step(dict(sequence=0,seconds=.02,forces=[dict(id=target,force_n=force,point_m=point)]))
    impulse=np.sum(scene.body.mass[:,None]*scene.body.v,axis=0)
    assert np.allclose(impulse,[.02,0,0],atol=1e-9)
    assert output['entities'][target]['translation_m'][0]>0
    old=scene.body.x.copy();oldtime=scene.body.time
    for command in [dict(sequence=0,seconds=.02),dict(sequence=1,seconds=.021),
                    dict(sequence=1,forces=[dict(id='unknown',force_n=[1,0,0],point_m=[0,0,0])]),
                    dict(sequence=1,forces=[dict(id=target,force_n=[101,0,0],point_m=point)])]:
        try:scene.step(command)
        except ValueError:pass
        else:raise AssertionError('Invalid force command accepted')
        assert scene.sequence==1 and scene.body.time==oldtime and np.array_equal(scene.body.x,old)
    # Force a post-integration domain rejection and check transactional state.
    scene.body.x[:,0]+=1.;before=scene.body.x.copy()
    try:scene.step(dict(sequence=1,seconds=.02))
    except ValueError:pass
    else:raise AssertionError('Out-of-domain movement accepted')
    assert scene.sequence==1 and np.array_equal(scene.body.x,before) and scene.body.time==oldtime
    record=dict(passed=True,canonical_entities=len(scene.body.ids),external_impulse_ns=impulse.tolist(),
                sphere_energy_residual_j=float(residual),scene=str(parent.relative_to(root)),
                scope='Reduced body translation/affine mechanics and ideal supported environments; not full-body collision validation')
    (parent/'verification.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--light',action='store_true');args=parser.parse_args()
    light_verification() if args.light else main()
