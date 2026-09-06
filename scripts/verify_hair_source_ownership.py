"""Tiny inferred ownership and actual-moment partition tests; no native job."""
from pathlib import Path
import copy,sys,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.hair_source_ownership import infer_source_owners,partition_represented_guides,classify_root_eligibility,select_guide_geometry


def rejects(call):
    try:call()
    except ValueError:return
    raise AssertionError('Invalid ownership/partition accepted')


def main():
    eligibility=classify_root_eligibility({'sample_ids':[10,11,12],'source_face_indices':[0,1,2]},[0,1,45],[0],selected_component_id=0)
    assert eligibility['eligible_guide_indices']==[0] and [r['sample_id'] for r in eligibility['ineligible_guides']]==[11,12]
    patch={'source_id':'skin','source_sha256':'a'*64,'source_node_indices':[4,5,6],'source_face_indices':[7],
           'positions_m':[[-.001,0,0],[.001,0,0],[.001,.002,0]],'triangles':[[0,1,2]]}
    registration={'groups':{'left':{'bounds_min_m':[-1,-1,-1],'bounds_max_m':[0,1,1]},'right':{'bounds_min_m':[0,-1,-1],'bounds_max_m':[1,1,1]}}}
    owners=infer_source_owners({'sample':patch},registration,registration_sha256='b'*64,ambiguity_margin_m=.003)
    rows=owners['bindings']['sample']['nodes']
    assert {r['owner'] for r in rows}=={'right'}
    assert owners['summary']['overridden_nodes']==1 and owners['summary']['ambiguous_nodes']==3
    assert owners['summary']['coherent_components']==1 and owners['native_enabled'] is False
    assert owners['bindings']['sample']['source_sha256']==patch['source_sha256']
    assert patch['source_node_indices']==[4,5,6] and patch['source_face_indices']==[7]
    other=copy.deepcopy(patch);other['positions_m'][0][0]=1
    rejects(lambda:infer_source_owners({'sample':patch,'mismatch':other},registration,registration_sha256='b'*64))
    strands={'centerlines_m':[0,0,0,0,0,.005,0,0,.01], 'strand_offsets':[0,3],'radius_m':[4e-5], 'density_kg_m3':1312.,'tensile_modulus_pa':7.11e9,'bending_modulus_pa':5.7e9}
    prepared={'sample':{'strands':strands,'guide_owner_names':['right'],'sample_ids':[123]}}
    native=[{'owner':'left','frame':'native-body-local:left','mass_kg':1.,'centroid_m':[0,0,0],'inertia_com_kg_m2':(np.eye(3)*.1).tolist(),'canonical_reference_to_body_local':np.eye(4).tolist()},
            {'owner':'right','frame':'native-body-local:right','mass_kg':1.,'centroid_m':[0,0,0],'inertia_com_kg_m2':(np.eye(3)*.1).tolist(),'canonical_reference_to_body_local':np.eye(4).tolist()}]
    result=partition_represented_guides(prepared,native)
    assert result['represented_guide_count']==1 and result['modified_owner_count']==1
    assert result['native_enabled'] is False
    assert abs(result['conservation']['mass_residual_kg'])<1e-15
    assert np.linalg.norm(result['conservation']['first_moment_residual_kg_m'])<1e-18
    assert np.linalg.norm(result['conservation']['second_moment_residual_kg_m2'])<1e-18
    invalid=copy.deepcopy(native);invalid[1]['inertia_com_kg_m2']=(np.eye(3)*1e-20).tolist()
    rejects(lambda:partition_represented_guides(prepared,invalid))
    duplicate=copy.deepcopy(native);duplicate[1]['owner']='left';rejects(lambda:partition_represented_guides(prepared,duplicate))
    geometry={'strands':strands,'attachment':{'sample_ids':[123],'barycentric':[[1,0,0]],'reference_triangles_m':[0]*9}}
    selected=select_guide_geometry(geometry,[0])
    assert np.array_equal(np.asarray(selected['strands']['centerlines_m']).reshape(-1,3),np.asarray(strands['centerlines_m']).reshape(-1,3))
    assert selected['strands']['radius_m']==strands['radius_m'] and selected['attachment']['sample_ids']==[123]
    from ihm.assembly.source_skin import file_sha256
    root=Path(__file__).resolve().parents[1];artifact=root/'data/research/hair_native_partition.json'
    if artifact.exists():
        report=json.loads(artifact.read_bytes())
        for path,digest in report['source_dependencies'].items():assert file_sha256(root/path)==digest,path
        assert report['partition']['represented_guide_count']==94 and report['native_enabled'] is False and report['installed'] is False
        assert [r['sample_id'] for r in report['groups']['body']['guide_eligibility']['ineligible_guides']]==[185247,145098]
        assert report['groups']['body']['population_eligibility']['ineligible_count']==11827
    print('PASS fixed exact nodes/faces, explicit ambiguity/coherence overrides, incompatible-source rejection, represented-guide mass/full-moment closure, invalid-inertia rejection and native disabled')


if __name__=='__main__':main()
