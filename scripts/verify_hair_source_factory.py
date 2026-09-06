"""Small exact-source factory checks; no anatomical geometry or native engine."""
from pathlib import Path
import copy,gzip,hashlib,json,sys,tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.source_skin import extract_source_patch,SourceSkinPatch
from ihm.assembly.hair_source_factory import map_retained_guides,prepare_hair_factory,guide_mass_properties,partition_native_inertia,transform_mass_properties


def rejects(call):
    try:call()
    except (ValueError,KeyError):return
    raise AssertionError('Invalid source mapping accepted')


def main():
    with tempfile.TemporaryDirectory(prefix='hair-source-') as folder:
        source=Path(folder)/'skin.json.gz';source.write_bytes(gzip.compress(json.dumps({'positions':[0,0,0,1,0,0,0,1,0,2,0,0,2,1,0],'indices':[0,1,2,1,3,4]}).encode()))
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        population={'ids':np.array([8,3]),'face_index':np.array([1,0]),'barycentric':np.array([[1,0,0],[1,0,0]]),'roots_m':np.array([[1.,0,0],[0,0,0]]),'normals':np.array([[0.,0,1],[0,0,1]]),'radius_m':np.array([4e-5,4e-5]),'length_m':np.array([.01,.01])}
        attachment={'skin_entity_id':'skin','source_geometry_sha256':digest,'sample_ids':[3],'barycentric':[[1,0,0]],'reference_triangles_m':[0,0,0,1,0,0,0,1,0]}
        strands={'centerlines_m':[0,0,0,0,0,.005,0,0,.01],'strand_offsets':[0,3],'radius_m':[4e-5],'density_kg_m3':1312.,'tensile_modulus_pa':7.11e9,'bending_modulus_pa':5.7e9}
        mapping=map_retained_guides(attachment,population,source_id='skin',source_sha256=digest)
        assert mapping['source_face_indices']==[0] and mapping['population_indices']==[1]
        extracted=extract_source_patch(source,digest,mapping['source_face_indices'],source_id='skin')
        binding={'source_id':'skin','source_sha256':digest,'registration_sha256':'a'*64,'basis':'explicit synthetic fixture','nodes':[{'source_node_index':i,'owner':'body'} for i in [0,1,2]]}
        patch=SourceSkinPatch.bind(extracted,binding,expected_registration_sha256='a'*64,owner_names=['body'])
        prepared=prepare_hair_factory({'attachment':attachment,'strands':strands},population,patch)
        assert prepared['follicle_triangles']==[0] and prepared['native_enabled'] is False
        assert prepared['source_face_indices']==[0] and prepared['source_node_indices']==[0,1,2]
        assert not patch.positions_m.flags.writeable
        rejects(lambda:extract_source_patch(source,'0'*64,[0],source_id='skin'))
        corrupt=copy.deepcopy(population);corrupt['ids']=np.array([3,3]);rejects(lambda:map_retained_guides(attachment,corrupt,source_id='skin',source_sha256=digest))
        corrupt=copy.deepcopy(attachment);corrupt['sample_ids']=[99];rejects(lambda:map_retained_guides(corrupt,population,source_id='skin',source_sha256=digest))
        corrupt=copy.deepcopy(binding);corrupt['nodes'][0]['owner']='unknown';rejects(lambda:SourceSkinPatch.bind(extracted,corrupt,expected_registration_sha256='a'*64,owner_names=['body']))
        corrupt=copy.deepcopy(binding);corrupt['nodes'].pop();rejects(lambda:SourceSkinPatch.bind(extracted,corrupt,expected_registration_sha256='a'*64,owner_names=['body']))
        corrupt=copy.deepcopy(binding);corrupt['registration_sha256']='b'*64;rejects(lambda:SourceSkinPatch.bind(extracted,corrupt,expected_registration_sha256='a'*64,owner_names=['body']))
        corrupt=copy.deepcopy(binding);corrupt['nodes'][0]['owner']='other'
        mixed=SourceSkinPatch.bind(extracted,corrupt,expected_registration_sha256='a'*64,owner_names=['body','other'])
        rejects(lambda:prepare_hair_factory({'attachment':attachment,'strands':strands},population,mixed))
        corrupt=copy.deepcopy(attachment);corrupt['reference_triangles_m'][0]=.001
        rejects(lambda:prepare_hair_factory({'attachment':corrupt,'strands':strands},population,patch))
        retained=source.read_bytes();source.write_bytes(b'changed after capture')
        assert extract_source_patch(retained,digest,[0],source_id='skin')==extracted
        rejects(lambda:extract_source_patch(source,digest,[0],source_id='skin'))
        properties=guide_mass_properties(strands,owner='body')
        assert abs(properties['mass_kg']-1312*np.pi*(4e-5)**2*.01)<1e-20
        native={'mass_kg':1.,'centroid_m':[0,0,0],'inertia_com_kg_m2':(np.eye(3)*.1).tolist(),'owner':'body','frame':'canonical-reference-m'}
        debit=partition_native_inertia(native,properties,partition_basis='explicit fixture allocation only')
        assert debit['native_enabled'] is False
        assert abs(debit['remaining']['mass_kg']+properties['mass_kg']-1)<1e-15
        assert np.allclose(np.array(debit['remaining']['first_moment_kg_m'])+properties['first_moment_kg_m'],[0,0,0],atol=1e-20)
        assert np.allclose(np.array(debit['remaining']['second_moment_origin_kg_m2'])+properties['second_moment_origin_kg_m2'],np.eye(3)*.05,atol=1e-18)
        transform=np.eye(4);transform[:3,3]=[.1,.2,.3]
        moved=transform_mass_properties(properties,transform,frame='body-local',owner='body')
        restored=transform_mass_properties(moved,np.linalg.inv(transform),frame=properties['frame'],owner='body')
        assert np.allclose(restored['second_moment_origin_kg_m2'],properties['second_moment_origin_kg_m2'],atol=1e-20)
        rejects(lambda:partition_native_inertia(native,moved,partition_basis='wrong frame'))
        negative={**properties,'mass_kg':-1};rejects(lambda:partition_native_inertia(native,negative,partition_basis='invalid'))
        tiny={**native,'mass_kg':properties['mass_kg']/2};rejects(lambda:partition_native_inertia(tiny,properties,partition_basis='invalid'))
        invalid_tensor={**native,'inertia_com_kg_m2':(np.eye(3)*1e-20).tolist()};rejects(lambda:partition_native_inertia(invalid_tensor,properties,partition_basis='invalid tensor'))
    print('PASS exact ID/face join, bounded source extraction, source/hash/owner rejection, no nearest fallback, inertial partition closure and no native enablement')


if __name__=='__main__':main()
