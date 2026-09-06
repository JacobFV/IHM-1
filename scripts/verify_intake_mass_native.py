"""Actual two-owner intake endpoint fixture; not whole-body signed energy acceptance."""
from pathlib import Path
import argparse,hashlib,json,tempfile,time
import numpy as np
from ihm.native.session import SessionConfig,Meal
from ihm.native.coupled_session import SignedCoupledNativeSession
from ihm.assembly.articulated import ArticulatedBodyPlant
from ihm.assembly.intake_mass import IntakeMassBridge


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true');args=parser.parse_args()
    if not args.run:raise SystemExit('Explicit --run and coordinated native slot required')
    root=Path(__file__).resolve().parents[1]
    output=Path(tempfile.mkdtemp(prefix='intake-mass-bridge-',dir=root/'data/derived/audits'))
    paths=['ihm/assembly/intake_mass.py','ihm/assembly/articulated.py','ihm/native/mechanical_stream.py',
        'ihm/native/coupled_session.py','ihm/native/instance_mass.py','scripts/verify_intake_mass_native.py']
    sources={name:(root/name).read_bytes() for name in paths}
    for name,raw in sources.items():
        p=output/'inputs'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
    started=time.monotonic();native=plant=None
    try:
        reference=json.loads((root/'data/derived/systemic/exertion_v3/exercise/native/manifest.json').read_text())
        state=Path(reference['configuration']['state_path'])
        assert hashlib.sha256(state.read_bytes()).hexdigest()==reference['state_sha256']
        native=SignedCoupledNativeSession(SessionConfig(state_path=state,
            engine_variant='whole_body_integrity_gi_absorption',horizon_s=.04),output/'physiology')
        initial_native=native.snapshot();mass=initial_native['values']['patient_weight_kg']
        plant=ArticulatedBodyPlant(root,output/'mechanics',environment='free',target_mass_kg=mass,
            augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json',
            instance_mass_variant='data/runtime/opensim/variants/instance_mass_v1')
        initial=plant.snapshot();canonical=json.loads((output/'mechanics/canonical_mechanics.json').read_text())
        stomach=next(e for e in canonical['entities'] if e['name']=='stomach')
        registration=plant.registration;assert registration.rows[stomach['id']]['body']=='torso'
        transform=np.array(plant.native.snapshot()['bodies']['torso']['transform_ground'])
        station=(np.linalg.inv(transform)@np.linalg.inv(registration.global_map)@np.r_[stomach['centroid_m'],1.])[:3].tolist()
        registration_hash=hashlib.sha256((output/'mechanics/registration.json').read_bytes()).hexdigest()
        binding={'canonical_entity_id':stomach['id'],'station_source_m':station,
            'registration_sha256':registration_hash,'registration_basis':registration.rows[stomach['id']],
            'scope':'Canonical stomach centroid on inferred rigid torso support; no internal organ deformation'}
        (output/'site.json').write_text(json.dumps(binding,indent=2)+'\n')
        bridge=IntakeMassBridge(plant,initial_native['intake'],body='torso',station_m=station,
            registration_identity=registration_hash,incoming_velocity_basis='co_moving_at_ingestion_assumption')
        frames=[]
        for i,water in enumerate((10.,20.)):
            queued=native.meal(Meal(name=f'actual_bridge_water_{i}',water_ml=water))
            assert bridge.apply(queued['intake']) is None
            before=plant.advance(.02)
            consumed=native.step(.02)
            result=bridge.apply(consumed['intake']);after=plant.snapshot()
            assert result['mechanical_transfer_applied'] and abs(result['mass_kg']-water*.001)<1e-12
            assert before['time_s']==after['time_s']
            for name,coordinate in before['joints'].items():
                assert coordinate['value']==after['joints'][name]['value']
                assert abs(coordinate['speed']-after['joints'][name]['speed'])<1e-10
            for key in ('muscle_metabolic_energy_j','signed_active_fiber_work_j','muscle_heat_energy_j','metabolic_reference'):
                assert before[key]==after[key],key+' changed during endpoint mass transfer'
            assert abs(after['effective_native_body_mass_kg']-initial['effective_native_body_mass_kg']-bridge.applied_mass)<1e-10
            receipt=result['mechanical_receipt']
            for key in ('generalized_impulse_residual','linear_momentum_residual_kg_m_s','angular_momentum_residual_kg_m2_s','constraint_velocity_error'):
                assert abs(receipt[key])<1e-8,(key,receipt[key])
            assert bridge.apply(native.snapshot()['intake']) is None
            frames.append({'consumed':consumed['intake'],'boundary':result,'time_s':after['time_s'],
                'mechanical_mass_kg':after['effective_native_body_mass_kg'],'bridge':bridge.snapshot()})
        assert bridge.sequence==2 and abs(bridge.applied_mass-.03)<1e-12
        assert all((root/name).read_bytes()==raw for name,raw in sources.items()),'Source changed during acceptance'
        native_process=native.process;mechanical_process=plant.native.process
        plant.close();plant=None;native.close();native=None
        assert native_process.poll() is not None and mechanical_process.poll() is not None
        report={'passed':True,'source_sha256':{name:hashlib.sha256(raw).hexdigest() for name,raw in sources.items()},
            'frames':frames,'wall_s':time.monotonic()-started,'both_processes_reaped':True,
            'scope':'Actual native consumed intake drives localized constraint-consistent mechanical mass at common endpoints; free-body fixture only, no signed metabolic coupling, supported equilibrium or excretion acceptance'}
        (output/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'passed':True,'output':str(output),'wall_s':report['wall_s'],'applied_mass_kg':bridge.applied_mass},indent=2))
    except BaseException as error:
        (output/'failure.json').write_text(json.dumps({'error':str(error),'type':type(error).__name__,'wall_s':time.monotonic()-started},indent=2)+'\n')
        raise
    finally:
        if plant:plant.close()
        if native:native.close(graceful=False)

if __name__=='__main__':main()
