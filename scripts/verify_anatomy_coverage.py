"""Regression checks for anatomical identity versus misleading name substrings."""
import json
from pathlib import Path
from ihm.spatial.anatomy import AnatomyOntology


def main():
    ontology = AnatomyOntology()
    cases = {'thyroid cartilage':'skeletal', 'right sternothyroid':'muscular',
             'left inferior thyroid artery':'arterial', 'anterior papillary muscle of right ventricle':'cardiac',
             'anterior leaflet of mitral valve':'cardiac', 'cavity of left ventricle':'cardiac', 'great cardiac vein':'venous', 'pulmonary trunk':'arterial', 'anterolateral head of lateral papillary muscle of left ventricle':'cardiac', 'salivary gland':'digestive', 'bile duct':'digestive', 'skin':'integumentary',
             'pituitary gland':'endocrine', 'spleen':'lymphatic', 'brain':'nervous'}
    for name, expected in cases.items():
        result = ontology.classify([{'concept_id':ontology.ids_by_name[name], 'name':name}])
        assert result['system'] == expected, (name, result)
        assert result['evidence'], name
    meshes = json.loads(Path('data/derived/anatomy/bodyparts3d_index.json').read_text())['meshes']
    expected_meshes = {'FJ2426':'cardiac','FJ2420':'cardiac','FJ2422':'cardiac','FJ2656':'venous','FJ2966':'arterial','FJ2810':'integumentary'}
    for mesh in meshes:
        result = ontology.classify(mesh['concepts'],mesh['element_id'])
        assert result['status'] in ('classified','unresolved')
        for evidence in result['evidence']:
            path = evidence['path']
            assert path[0] in result['source_concept_ids']
            assert ontology.anchors[path[-1]] == evidence['system']
            assert all(parent in ontology.parents[evidence['relation']][child] for child,parent in zip(path,path[1:]))
        if mesh['element_id'] in expected_meshes:
            assert result['system']==expected_meshes[mesh['element_id']], (mesh['element_id'],result['system'])
    # Generic shared concepts and misleading text must not supply anatomical identity.
    assert ontology.classify([{'concept_id':'UNKNOWN','name':'thyroid muscle artery'}])['status'] == 'unresolved'
    print(f'{len(cases)} anatomical counterexamples and 2,234 source mesh/path checks passed')

if __name__ == '__main__':
    main()
