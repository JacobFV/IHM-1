"""Regenerate a per-mesh BP3D classification crosswalk and coverage report."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from ihm.spatial.anatomy import AnatomyOntology

OUT = Path('data/derived/anatomy/anatomy_coverage.json')

# Frozen historical display grouping, solely for measuring the correction.
LEGACY = [('cardiac',['heart','myocardium','cardiac','atrioventricular','wall of atrium','wall of ventricle']),('skeletal',['bone organ','tooth','cartilage']),('muscular',['muscle organ','muscle tissue']),('arterial',['artery','arterial']),('venous',['vein','venous']),('lymphatic',['lymph','spleen','thymus']),('nervous',['nerve','brain','spinal cord','ganglion']),('respiratory',['lung','bronch','trachea','larynx','nasal']),('urinary',['kidney','ureter','urinary','renal pelvis']),('reproductive',['testis','penis','prostate','seminal','epididym','ductus deferens']),('endocrine',['thyroid','adrenal','pituitary','pineal']),('digestive',['liver','pancreas','stomach','intestin','colon','rectum','esophag','gallbladder','duoden','jejun','ileum']),('integumentary',['skin','nail','hair','mammary']),('connective',['ligament','tendon','fascia','aponeurosis'])]


def legacy_system(concepts):
    text = ' | '.join(c['name'].lower() for c in concepts)
    return next((system for system, words in LEGACY if any(w in text for w in words)), 'other')


def build():
    source = Path('data/derived/anatomy/bodyparts3d_index.json')
    meshes = json.loads(source.read_text())['meshes']
    ontology = AnatomyOntology()
    freq = Counter(c['concept_id'] for m in meshes for c in m['concepts'])
    records = []
    for mesh in meshes:
        classification = ontology.classify(mesh['concepts'], mesh['element_id'])
        name = classification['preferred_name'] or min(mesh['concepts'],key=lambda c:(freq[c['concept_id']],-len(c['name'])))['name']
        records.append({'element_id':mesh['element_id'], 'name':name,
                        'source_sha256':mesh['sha256'], 'source_path':mesh['source_path'],
                        'source_faces':mesh['faces'], 'source_vertices':mesh['vertices'],
                        'source_bounds_mm':mesh['bounds_in_source_coordinates'],
                        'legacy_system':legacy_system(mesh['concepts']), **classification})
    summary = {'mesh_count':len(records), 'system_counts':dict(sorted(Counter(r['system'] for r in records).items())),
               'legacy_system_counts':dict(sorted(Counter(r['legacy_system'] for r in records).items())),
               'unresolved_count':sum(r['status']=='unresolved' for r in records),
               'changed_count':sum(r['system']!=r['legacy_system'] for r in records),
               'classification_basis_counts':dict(Counter(r['classification_basis'] for r in records)),
               'ambiguous_primary_count':sum(r['ambiguous_primary'] for r in records),
               'multisystem_count':sum(len(r['systems'])>1 for r in records)}
    skin = [r for r in records if any(ontology.names.get(cid)=='skin' for cid in r['source_concept_ids'])]
    summary['skin_source_meshes'] = [{k:r[k] for k in ('element_id','name','source_faces','source_vertices','source_bounds_mm','source_sha256')} for r in skin]
    report = {'schema_version':1, 'scope':'BodyParts3D 4.0 source surface coverage; not complete human anatomy',
              'index_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
              'ontology_sha256':ontology.ontology_sha256, 'crosswalk_sha256':ontology.crosswalk_sha256,
              'summary':summary, 'meshes':records}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2)+'\n')
    lines = ['# BodyParts3D anatomy coverage', '',
             'Regenerate with `.venv/bin/python scripts/audit_anatomy_coverage.py`; verify with `.venv/bin/python scripts/verify_anatomy_coverage.py`.', '',
             f"All {len(records):,} acquired source meshes have an explicit classification status. {summary['unresolved_count']} are unresolved; {summary['changed_count']} display groupings changed from the historical substring classifier.", '',
             'The crosswalk uses exact recorded FMA concept anchors and directed child-to-parent paths in the downloaded is-a and part-of tables. Is-a identity takes precedence over part-of context. All supported system tags remain available; a primary display group does not imply an exclusive biological function. Each mesh record retains source concept IDs, paths with labels, six source-table SHA-256 hashes, a crosswalk hash, and its source OBJ hash.', '',
             '| Primary system | Historical | Source crosswalk |', '|---|---:|---:|']
    for system,count in summary['system_counts'].items():
        lines.append(f"| {system} | {summary['legacy_system_counts'].get(system,0)} | {count} |")
    lines += ['', '## Corrected examples', '', '| Source mesh | Recorded name | Historical | Source crosswalk |', '|---|---|---|---|']
    changes = [r for r in records if r['system']!=r['legacy_system']]
    for system in summary['system_counts']:
        for r in [r for r in changes if r['system']==system][:3]:
            lines.append(f"| {r['element_id']} | {r['name']} | {r['legacy_system']} | {r['system']} |")
    lines += ['', '## Ambiguity and contextual membership', '', f"{summary['classification_basis_counts'].get('isa_identity',0)} meshes have is-a identity evidence; {summary['classification_basis_counts'].get('partof_context',0)} rely on part-of context. {summary['ambiguous_primary_count']} primary assignments have equally near competing system anchors.", '', 'Pancreatic structures retain digestive and endocrine membership. Source part-of aggregates also associate some cardiac structures with the systemic arterial system and liver with venous structures; explicit is-a anchors for valve cusps, chamber cavities, papillary muscle regions and liver segments take precedence. Aggregate tags describe recorded source associations, not an assertion that a valve is an artery or a liver segment is a vein.', '', '## Skin and remaining anatomical limits', '']
    for r in skin:
        bounds=r['source_bounds_mm']; spans=[round(bounds[1][i]-bounds[0][i],2) for i in range(3)]
        lines.append(f"Source skin `{r['element_id']}` has {r['source_vertices']:,} vertices and {r['source_faces']:,} triangular faces; bounds in the official millimeter frame are `{bounds}`, with axis spans `{spans}` mm. This is a source surface shell, not layered epidermis/dermis or a skin physiology model. Display decimation must be reported separately from these source counts.")
    lines += ['', 'BP3D lymphatic coverage here consists of spleen and two thymus lobes. No lymph-node or lymph-vessel meshes are present in these downloaded element tables. Three classified lymphatic surfaces do not constitute a lymphatic network. Additional sources must remain separate coordinate frames until registration is validated.', '',
              'The sensory display group covers eye surfaces, chambers, lacrimal structures and external ear. Nerves retain nervous identity; associated muscles retain muscular identity. Thyroid cartilage is skeletal, sternothyroid is muscular, and inferior thyroid artery is arterial. Papillary cardiac muscles retain cardiac identity. Generic anatomical terms without a recorded path remain explicitly unresolved rather than being inferred from substrings.', '',
              'This classification audit establishes source inventory and traceability. It does not establish tissue completeness, watertightness, subject calibration, or physiological validity. The downloaded BP3D ontology is a selected representation, not the complete FMA release.', '']
    Path('docs/ANATOMY_COVERAGE.md').write_text('\n'.join(lines))
    print(json.dumps(summary,indent=2))
    return report


if __name__ == '__main__':
    build()
