"""Retain identity-bound native composition ledger; never mutate a native model."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def moments(props,transform=None):
    mass=float(props['mass_kg']);center=np.asarray(props['center_m']);inertia=np.asarray(props['inertia_kg_m2'])
    if transform is not None:
        t=np.asarray(transform);r=t[:3,:3]
        if not np.allclose(r.T@r,np.eye(3),atol=1e-10,rtol=0) or np.linalg.det(r)<=0:raise ValueError('Proper rigid body transform required')
        center=r@center+t[:3,3];inertia=r@inertia@r.T
    second=np.trace(inertia)/2*np.eye(3)-inertia+mass*np.outer(center,center)
    return {'mass_kg':mass,'H_first_moment_kg_m':(mass*center).tolist(),'Q_second_moment_kg_m2':second.tolist()}


def summed(rows):
    return {key:np.sum([np.asarray(row[key]) for row in rows],axis=0).tolist() for key in rows[0]}


def build():
    paths={'cervical_inertia':ROOT/'data/research/cervical_inertia/v2/manifest.json',
           'cervical_registration':ROOT/'data/research/cervical_registration/v2/recipe.json',
           'thoracic_anatomy':ROOT/'data/research/thoracic_anatomy/v3/manifest.json',
           'thoracic_mechanism':ROOT/'data/research/thoracic_mechanism/v2/manifest.json'}
    raw={k:p.read_bytes() for k,p in paths.items()};data={k:json.loads(v) for k,v in raw.items()}
    receipts={k:{'path':str(paths[k].relative_to(ROOT)),'sha256':hashlib.sha256(v).hexdigest(),'bytes':len(v)} for k,v in raw.items()}
    prior=data['cervical_inertia'];neck=data['cervical_registration'];thorax=data['thoracic_anatomy'];mechanism=data['thoracic_mechanism']
    if neck['inertial_prior_manifest_identity']['sha256']!=receipts['cervical_inertia']['sha256']:raise ValueError('Cervical prior identity differs')
    if mechanism['anatomy_manifest']['sha256']!=receipts['thoracic_anatomy']['sha256']:raise ValueError('Thorax anatomy identity differs')
    if neck['target_model_identity']['sha256']!=thorax['target_model_identity']['sha256']:raise ValueError('Native targets differ')
    cervical={key:{'body_name':body['name'],'body_to_target_torso_reference':body['body_to_target_torso_reference'],
                   'mass_properties_in_new_body_frame':body['mass_properties_in_new_body_frame'],
                   'reference_moments_in_current_torso':moments(body['mass_properties_in_new_body_frame'],body['body_to_target_torso_reference'])}
              for key,body in neck['bodies'].items()}
    materials={key:moments(thorax['entities'][key]['torso_frame_mass_prior']) for key in mechanism['materials']}
    core=moments(thorax['replace_reduced_torso_with']);neck_total=summed([b['reference_moments_in_current_torso'] for b in cervical.values()])
    thorax_total=summed([core,*materials.values()]);target=moments(prior['current_torso']);composed=summed([neck_total,thorax_total])
    for key in target:
        if not np.allclose(composed[key],target[key],atol=2e-12,rtol=0):raise ValueError('Combined mass/first/second moment mismatch')
        if not np.allclose(thorax_total[key],moments(prior['residual_torso'])[key],atol=2e-12,rtol=0):raise ValueError('Cervical residual mismatch')
    return {'schema':'thoracic-native-composition-plan-v1','native_activation_allowed':False,'source_receipts':receipts,
        'target_native_model_identity':neck['target_model_identity'],'reference':'Current native torso frame; neck q=0; thorax q=0; not a moving-pose inertia constraint',
        'moment_convention':'H=integral x dm; Q=integral x x^T dm about current torso origin. Q here is a second moment, not generalized force.',
        'current_native_torso':target,'cervical_bodies':cervical,'cervical_aggregate':neck_total,
        'thoracic_residual_core':core,'thoracic_material_partitions':materials,'thoracic_subsystem_aggregate':thorax_total,'composed_aggregate':composed,
        'composition_order':['Verify exact current native identity and all source receipts; apply no additional mass normalization.',
            'Create the nine cervical body frames and debit their full reference mass/H/Q once from the current torso.',
            'Replace the cervical-reduced torso by the thoracic core and existing material mechanism; do not add its subsystem mass to the unreduced torso.',
            'Attach both subsystems to one original torso frame/pose/twist and original upstream parent joint; assemble all parent cross-inertia terms in one generalized system.',
            'Preserve original 92 effectors and audit 78 donor paths/18 couplers/24 neck coordinates before any donor activation.',
            'Replace neck/head/rib/diaphragm source material and contact owners exactly once; retain explicit mappings and common global registration.',
            'Validate physical locked-rib wrenches separately if native joint reactions are requested; formal eliminated-coordinate loads cannot supply them.',
            'Transfer overlapping native chest compliance and ideal pressure drive only after signed displacement/work mapping acceptance.'],
        'unresolved_native_acceptance':neck['required_before_native_acceptance']+[
            'Implement cervical and thoracic dynamics in one velocity map; a standalone 20.235757 kg RHS omits neck-parent cross inertia.',
            'Reference moment conservation does not make the combined inertia constant when neck/thorax coordinates move.',
            'Validate free trajectory, mixed-map momentum, source contact ownership, native signed work and supported supine pose before opt-in use.']}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise ValueError('Fresh output directory required')
    plan=build();args.output.mkdir(parents=True);(args.output/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    print(json.dumps({k:plan[k]['mass_kg'] for k in ['current_native_torso','cervical_aggregate','thoracic_subsystem_aggregate','thoracic_residual_core','composed_aggregate']}))

if __name__=='__main__':main()
