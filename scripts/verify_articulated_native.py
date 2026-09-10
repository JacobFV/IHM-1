"""Explicitly gated, short native acceptance; never runs on import/default CLI."""
from pathlib import Path
import argparse,json,sys,tempfile,time,resource,traceback
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.articulated import ArticulatedBodyPlant
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG
from verify_native_signed_energy_ledger import ENERGIES,check_interval
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-native',action='store_true');parser.add_argument('--garments',action='store_true');parser.add_argument('--augmentation-registration',default='data/derived/mechanics/whole_body_arm26_v2/registration.json');args=parser.parse_args()
    if not args.run_native:raise SystemExit('Native acceptance requires an explicitly coordinated --run-native resource slot')
    output=Path(tempfile.mkdtemp(prefix='articulated-acceptance-',dir=ROOT/'data/derived'));started=time.monotonic();plant=None
    try:
        plant=ArticulatedBodyPlant(ROOT,output/'plant',environment='supine',target_mass_kg=MECHANICAL_TARGET_MASS_KG,augmented_registration=args.augmentation_registration,enable_garments=args.garments)
        # The entity count is read from the canonical mechanics the plant was built on,
        # not written here: it was 2408, then 2403 after the duplicate-surface collapse,
        # then 4000 once the display promotion landed, and a literal only ever gets
        # rewritten to match. The MECHANICAL_TARGET_MASS_KG kg is the native 22-body patient mass, which
        # is scaled independently of the canonical proxy partition and is deliberately
        # not the profile.
        expected=json.loads((ROOT/'data/derived/canonical/mechanics.json').read_bytes())['counts']['entities']
        initial=plant.snapshot();assert len(initial['muscles'])==92 and len(initial['entities'])==expected
        assert initial['total_muscle_metabolic_w']>0 and all(initial[k]==0 for k in ENERGIES)
        checkpoint=plant.checkpoint();dt=.00005 if args.garments else .002
        baseline=plant.advance(dt);plant.restore(checkpoint);repeat=plant.advance(dt)
        assert repeat['native_bodies']==baseline['native_bodies'] and all(repeat[k]==baseline[k] for k in (*ENERGIES,'metabolic_reference'))
        ledger=check_interval(initial,baseline,endpoint_check=not args.garments)
        plant.restore(checkpoint);active=plant.advance(dt,actuation={'arm26_BRA_r':.8})
        assert active['muscles']['arm26_BRA_r']['activation']>baseline['muscles']['arm26_BRA_r']['activation']
        assert active['total_muscle_metabolic_w']>baseline['total_muscle_metabolic_w']
        plant.restore(checkpoint);ident=plant.registration.groups['hand_r']['canonical_bones'][0];point=initial['entities'][ident]['centroid_m']
        forced=plant.advance(dt,[{'id':ident,'point_m':point,'force_n':[20.,0,0]}]);plant.release(checkpoint)
        assert np.linalg.norm(np.array(forced['entities'][ident]['translation_m'])-baseline['entities'][ident]['translation_m'])>1e-12
        for frame in (baseline,active,forced):
            assert np.linalg.norm(frame['audit']['momentum_balance_residual_source_n'])<1e-5
            assert np.isfinite(frame['total_muscle_metabolic_w']) and frame['muscle_metabolic_energy_j']>0
        plant.close()
        report={'passed':True,'dt_s':dt,'wall_s':time.monotonic()-started,'maximum_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                'muscles':len(initial['muscles']),'canonical_entities':len(initial['entities']),'metabolic_baseline_w':baseline['total_muscle_metabolic_w'],
                'metabolic_active_w':active['total_muscle_metabolic_w'],'garments_enabled':args.garments,'signed_energy_ledger':ledger,
                'scope':'Short actual native branching/checkpoint/metabolic/force acceptance; not settled support, full-body motion or full garment containment'}
        for name,frame in [('initial',initial),('baseline',baseline),('active',active),('forced',forced)]:
            (output/(name+'.json')).write_text(json.dumps(frame,allow_nan=False)+'\n')
        (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({**report,'output':str(output)},indent=2))
    except BaseException as error:
        if plant is not None:plant.close()
        (output/'failure.json').write_text(json.dumps({'passed':False,'error':str(error),'traceback':traceback.format_exc(),'wall_s':time.monotonic()-started,'maximum_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss},indent=2)+'\n')
        raise
    finally:
        if plant is not None:plant.close()
if __name__=='__main__':main()
