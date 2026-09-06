"""Verify a scheduled actual native snapshot against the source replacement metric."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_mechanism import ThoracicMechanism
from ihm.assembly.thoracic_cervical_metric import CoupledThoracicMetric
from ihm.assembly.thoracic_native_operator import native_frame_to_body_twist_map,compose_native_operator
ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--snapshot-dir',type=Path,required=True);args=parser.parse_args();directory=args.snapshot_dir.resolve()
    receipt=json.loads((directory/'manifest.json').read_bytes());snap_raw=(directory/'snapshot.json').read_bytes()
    if not receipt['native_snapshot_run'] or hashlib.sha256(snap_raw).hexdigest()!=receipt['snapshot']['sha256']:raise ValueError('Native snapshot receipt mismatch')
    model_raw=(directory/'assembled_model.osim').read_bytes();native_sha=hashlib.sha256(model_raw).hexdigest();snapshot=json.loads(snap_raw)
    if snapshot['native_integrated'] or not snapshot['zero_applied_force_inverse_dynamics']:raise ValueError('Unexpected native snapshot scope')
    metric=CoupledThoracicMetric(ThoracicMechanism(ROOT/'data/research/thoracic_mechanism/v2/manifest.json'),ROOT/'data/research/thoracic_mechanism/native_composition_v1/plan.json')
    u=np.asarray(snapshot['u']);native=np.asarray(snapshot['mass_matrix']);bias=np.asarray(snapshot['inertial_bias']);j,beta=native_frame_to_body_twist_map(snapshot['torso_rotation_body_to_world'],snapshot['torso_frame_jacobian_world_angular_first'],snapshot['torso_frame_bias_world_angular_first'],u)
    result=compose_native_operator(metric,native,bias,j,beta,u,np.zeros(32),np.zeros(32),native_model_sha256=native_sha)
    n=len(u);matrix=result['mass_matrix'];active=[i for i in range(n+32) if i not in result['eliminated_internal_indices']]
    native_error=float(np.max(np.abs(matrix[:n,:n]-native)));bias_error=float(np.max(np.abs(result['inertial_bias'][:n]-bias)));eigen=np.linalg.eigvalsh(matrix[np.ix_(active,active)])
    if native_error>2e-10 or bias_error>2e-10 or eigen[0]<=0:raise ValueError('Actual native reference replacement or rank failed')
    if snapshot['body_count']!=22 or snapshot['source_muscle_count']!=92 or abs(snapshot['torso_mass_kg']-27.654676965260336)>1e-10:raise ValueError('Wrong native mechanical inventory')
    report={'schema':'thoracic-actual-native-operator-proof-v1','snapshot_sha256':hashlib.sha256(snap_raw).hexdigest(),'native_model_sha256':native_sha,
        'native_speed_count':n,'composed_speed_count':n+32,'native_constraints_projected':False,'native_integrated':False,
        'reference_native_mass_block_max_error':native_error,'reference_native_bias_max_error':bias_error,'active_eigenvalues':eigen.tolist(),
        'native_to_new_internal_cross_mass_norm':float(np.linalg.norm(matrix[:n,n:])),'retained_source_subsystem_mass_kg':result['composite_state_evaluation']['mass_kg'],
        'scope':'Actual unconstrained tree operator extraction and reference replacement; no integrated native composite/contact/physiology acceptance'}
    target=directory/'composition_proof.json'
    if target.exists():raise ValueError('Proof already exists')
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
