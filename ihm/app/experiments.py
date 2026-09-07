"""Read body-attached experiment artifacts without silently accepting stale code."""
from pathlib import Path
import hashlib
import json
import gzip

EXPERIMENTS={'forearm-touch':'forearm-touch.json','skin-transport':'skin-transport.json','skin-electric':'skin-electric.json','spectra':'regional-spectra.json','details':'details.json'}


def read_experiment(root,kind):
    root=Path(root).resolve()
    if kind=='hair-strands':
        path=root/'data/derived/hair/elastic_v3/manifest_fragment.json'
        data=json.loads(path.read_bytes())
        for source,expected in data['source_hashes'].items():
            file=(root/source).resolve()
            if not file.is_relative_to(root) or hashlib.sha256(file.read_bytes()).hexdigest()!=expected:
                raise ValueError('Hair source changed; rebuild materialization')
        for structure in data['structures']:
            file=(root/structure['geometry_path']).resolve()
            if not file.is_relative_to(root) or hashlib.sha256(file.read_bytes()).hexdigest()!=structure['geometry_sha256']:
                raise ValueError('Hair strand geometry changed')
        return data
    if kind=='provenance' or kind.startswith('provenance-'):
        directory=root/'data/derived/structure-provenance-candidate-v1'
        manifest=json.loads((directory/'manifest.json').read_bytes())
        index_path=(directory/'index.json').resolve()
        if hashlib.sha256(index_path.read_bytes()).hexdigest()!=manifest['outputs_sha256']['index.json']:
            raise ValueError('Provenance index changed; rebuild scripts/build_structure_provenance.py')
        index=json.loads(index_path.read_bytes())
        if kind=='provenance':
            # The index answers for the whole scene, so its heavy inputs are verified here rather than per structure.
            for source,expected in manifest['inputs_sha256'].items():
                file=(root/source).resolve()
                if not file.is_relative_to(root) or hashlib.sha256(file.read_bytes()).hexdigest()!=expected:
                    raise ValueError('Provenance input changed; rebuild '+source)
            return {k:v for k,v in index.items() if k!='record_sha256'}
        identity=kind.removeprefix('provenance-')
        expected=index['record_sha256'].get(identity)
        if expected is None:raise ValueError('No provenance record for that structure')
        records=(directory/'records').resolve()
        path=(records/(identity+'.json')).resolve()
        if not path.is_relative_to(records) or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('Provenance record changed; rebuild export')
        return json.loads(path.read_bytes())
    if kind=='conforming-domains' or kind.startswith('conforming-domain-'):
        directory=root/'data/derived/conforming-domain-display-v1'
        index=json.loads((directory/'index.json').read_bytes())
        wanted=None if kind=='conforming-domains' else kind.removeprefix('conforming-domain-')
        entries=index['domains'] if wanted is None else [e for e in index['domains'] if e['id']==wanted]
        if not entries:raise ValueError('Unknown conforming domain')
        for entry in entries:
            # The display is a boundary projection of a volume; a changed volume invalidates it.
            for source,expected in entry['source_hashes'].items():
                file=(root/source).resolve()
                if not file.is_relative_to(root) or hashlib.sha256(file.read_bytes()).hexdigest()!=expected:
                    raise ValueError('Conforming domain source changed; rebuild '+entry['id'])
        if wanted is None:return index
        path=(directory/entries[0]['display_path']).resolve()
        if not path.is_relative_to(directory.resolve()) or hashlib.sha256(path.read_bytes()).hexdigest()!=entries[0]['display_sha256']:
            raise ValueError('Conforming domain display changed; rebuild export')
        return json.loads(gzip.decompress(path.read_bytes()))
    if kind=='systemic':
        index=root/'data/derived/canonical/systemic-index.json'
        return json.loads(index.read_text()) if index.exists() else {'runs':[]}
    if kind.startswith('systemic-'):
        protocol=kind.removeprefix('systemic-')
        from ihm.assembly.systemic import PROTOCOLS
        if protocol not in PROTOCOLS:raise ValueError('Unknown systemic protocol')
        index=read_experiment(root,'systemic')
        entry=next((r for r in index['runs'] if r['id']==protocol),None)
        if entry is None:raise ValueError('Systemic experiment unavailable')
        path=(root/entry['path']).resolve()
        if not path.is_relative_to(root) or hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256']:
            raise ValueError('Systemic display changed; rebuild index')
        data=json.loads(path.read_text())
        for source,expected in data['runtime_sources'].items():
            file=(root/source).resolve()
            if not file.is_relative_to(root) or hashlib.sha256(file.read_bytes()).hexdigest()!=expected:
                raise ValueError('Systemic source changed; rebuild '+protocol)
        return data
    if kind not in EXPERIMENTS and kind!='garment-contact':raise ValueError('Unknown body experiment')
    if kind=='spectra':
        # An unchanged record may itself have stale model/anatomical inputs.
        for dependency in ('forearm-touch','skin-electric'):
            read_experiment(root,dependency)
    if kind=='garment-contact':
        directory=root/'data/derived/garment-tissue-display-v2'
        manifest=json.loads((directory/'manifest.json').read_bytes())
        path=(directory/manifest['display_path']).resolve()
        if not path.is_relative_to(directory.resolve()) or hashlib.sha256(path.read_bytes()).hexdigest()!=manifest['display_sha256']:
            raise ValueError('Garment contact display changed; rebuild export')
        data=json.loads(gzip.decompress(path.read_bytes()))
        data['display_identity']={'path':str(path.relative_to(root)),'sha256':manifest['display_sha256'],
                                  'position_quantization_max_error_m':manifest['position_quantization_max_error_m']}
    else:
        data=json.loads((root/'data/derived/canonical'/EXPERIMENTS[kind]).read_bytes())
    hashes={**data.get('runtime_sources',{}),**data.get('source_hashes',{})}
    hashes.update({r['path']:r['sha256'] for r in data.get('model',{}).get('provenance',[])})
    for path,expected in hashes.items():
        file=(root/path).resolve()
        if not file.is_relative_to(root) or hashlib.sha256(file.read_bytes()).hexdigest()!=expected:
            raise ValueError('Experiment source changed; rebuild '+kind)
    if kind=='garment-contact' and data.get('configuration',{}).get('tissue_manifest_path'):
        domain=Path(data['configuration']['tissue_manifest_path']).resolve()
        if not domain.is_relative_to(root) or hashes.get(str(domain.relative_to(root)))!=data['configuration']['tissue_manifest_sha256']:
            raise ValueError('Garment tissue resolution source is not bound')
        metadata=json.loads(domain.read_bytes())
        data['geometry_precision']={key:metadata[key] for key in ['geometry_basis','spacing_m','boundary_discretization_diagonal_m']}
        data['geometry_precision']['source_manifest_sha256']=data['configuration']['tissue_manifest_sha256']
    anchor=data.get('anchor',{})
    for path,expected in anchor.get('source_hashes',{}).items():
        file=(root/path).resolve()
        if not file.is_relative_to(root) or hashlib.sha256(file.read_bytes()).hexdigest()!=expected:
            raise ValueError('Experiment anatomical source changed; rebuild '+kind)
    if anchor.get('geometry_path'):
        file=(root/anchor['geometry_path']).resolve()
        if not file.is_relative_to(root) or hashlib.sha256(file.read_bytes()).hexdigest()!=anchor['geometry_sha256']:
            raise ValueError('Experiment anchor geometry changed')
    return data
