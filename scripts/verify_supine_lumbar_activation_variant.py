"""Compare actual initialized native tone defaults with a static candidate receipt."""
from pathlib import Path
import sys,json,tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream

def main(variant):
    p=(ROOT/variant).resolve();pose=json.loads((p/'initial_pose.json').read_text());out=Path(tempfile.mkdtemp(prefix='supine-lumbar-default-receipt-',dir=ROOT/'data/derived'))
    support=pose.get('support_geometry') or {};kw={}
    if support:kw=dict(surface_contact_manifest=support['manifest_path'],bed_material=support['bed_material'])
    s=NativeMechanicalStream(ROOT,out/'native',environment='supine',target_mass_kg=77.6122029,augmented_registration=str((p/'registration.json').relative_to(ROOT)),initial_pose=dict(sorted(pose['coordinates'].items())),**kw)
    try:
        initial=s.snapshot();(out/'initialized.json').write_text(json.dumps(initial,indent=2)+'\n')
        q=pose['coordinates'];args=['evaluate_static_pose',str(len(q))]
        for n,v in q.items():args.extend([n,str(v)])
        native=s._request(' '.join(args));reference=pose['native_static_residual']
        diffs={k:float(np.max(np.abs(np.array(native[k])-reference[k]))) for k in ('udot','contact_force_n','constrained_zero_acceleration_residual_mobility_force')}
        actual={n:initial['muscles'][n]['activation'] for n in pose['activations']}
        ad=max(abs(actual[n]-v) for n,v in pose['activations'].items());qd=max(abs(native['coordinates'][n]['value']-v) for n,v in q.items())
        tolerances=dict(udot=1e-6,contact_force_n=1e-5,constrained_zero_acceleration_residual_mobility_force=1e-5,q=1e-8,activation=1e-12)
        report=dict(variant=str(p.relative_to(ROOT)),support_geometry=support,max_absolute_differences=diffs,max_q_difference=qd,max_activation_difference=ad,actual_initial_activations=actual,tolerances=tolerances,
            tolerance_basis='q tolerance matches existing native initializer assembly acceptance. Force/acceleration thresholds test numerical initialization equivalence, not equilibrium.',
            passed=all(diffs[k]<tolerances[k] for k in diffs) and ad<tolerances['activation'] and qd<tolerances['q'],
            comparison_basis='Same solved q, selected activation defaults in actual initialized native plant vs copied-state static activation overrides; native muscle fibers equilibrated independently.',native=native,accepted_equilibrium=False)
        (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='native'},indent=2));print(out)
    finally:s.close()
if __name__=='__main__':main(sys.argv[1])
