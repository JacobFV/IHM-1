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


WARDROBE='data/derived/wardrobe-v1/wardrobe.json'


def _wardrobe(root):
    """The 33 shrinkwrapped garments, read from the wardrobe the cloth run produced.

    The two offset shells under data/derived/clothing are the earlier pattern
    prior and are not served: they are neither registered onto the outer
    envelope nor cloth-simulated, so they cannot stand for this wardrobe.
    """
    root=Path(root).resolve()
    data=json.loads((root/WARDROBE).read_bytes())
    if data.get('garment_count')!=len(data['garments']):
        raise ValueError('Wardrobe declares a garment count its records do not match')
    return root,data


def _wardrobe_asset(root,relative,expected=None):
    """Assets are addressed inside the workspace and are served only when they hash as declared."""
    file=(root/relative).resolve()
    if not file.is_relative_to(root) or not file.is_file():return None
    if expected is not None and hashlib.sha256(file.read_bytes()).hexdigest()!=expected:
        raise ValueError('Wardrobe asset changed since the cloth run: '+str(relative))
    return file


def clothing_catalog(root):
    """Garment tiles read this: identity, label, occupied slots and thumbnail.

    Exclusivity here is slot-set intersection, not a single slot. A dress
    occupies torso_base and legs at once, so it excludes every garment that
    holds either. The wardrobe's own slot model carries that relation and is
    passed through untouched; nothing about it is decided in this function or
    in the front end.
    """
    root,data=_wardrobe(root)
    model=data['slot_model']
    occupies={entry['garment']:entry['occupies'] for entry in model['garments']}
    excludes={entry['garment']:entry['excludes'] for entry in model['garments']}
    garments=[]
    for garment in data['garments']:
        identity=garment['id']
        slots=garment.get('slots') or occupies.get(identity) or []
        if not slots:raise ValueError('Wardrobe garment declares no slot: '+identity)
        if occupies.get(identity) not in (None,slots):
            raise ValueError('Wardrobe slot model disagrees with the garment record for '+identity)
        thumbnail=garment.get('thumbnail')
        if thumbnail is not None and _wardrobe_asset(root,thumbnail) is None:
            raise ValueError('Garment thumbnail is missing: '+str(thumbnail))
        geometry=garment.get('geometry')
        if geometry is not None and _wardrobe_asset(root,geometry) is None:
            raise ValueError('Garment geometry is missing: '+str(geometry))
        garments.append({'id':identity,'label':garment.get('name') or identity,
                         'slots':slots,'slot':slots[0],'excludes':excludes.get(identity,[]),
                         'colour':garment.get('colour'),
                         'thumbnail':thumbnail,
                         'thumbnail_url':None if thumbnail is None else '/api/clothing/thumbnail/'+identity,
                         'geometry_url':None if geometry is None else '/api/clothing/geometry/'+identity,
                         'fit_mode':garment.get('fit_mode'),'loose_fit':garment.get('loose_fit',False),
                         'mass_kg':garment.get('mass_kg'),'author':garment.get('author'),'license':garment.get('license'),
                         'inside_body_vertices':garment.get('inside_body_vertices'),
                         'max_penetration_mm':garment.get('max_penetration_mm'),
                         'evidence_kind':'registered_cc0_garment_cloth_solve',
                         'physical_contact_solved':bool(garment.get('contact_ready'))
                                                   and bool(garment.get('contacted_body_in_run'))})
    return {'schema':'ihm.clothing-catalog.v2','source':str(Path(WARDROBE)),
            'body':data.get('body'),'registration':data.get('registration'),
            'draw_order':model['draw_order'],'draw_order_basis':model.get('draw_order_basis'),
            'empty_slots':model.get('empty_slots',[]),'unfilled_slots':data.get('unfilled_slots',{}),
            'exclusivity_model':'Slot sets carry exclusivity: two garments are mutually exclusive when the slots they occupy intersect, and combine when they do not.',
            'loose_fit_note':data.get('loose_fit_note'),'verification':data.get('verification'),
            'garments':garments}


def clothing_asset(root,ident,field):
    """Serve one garment's thumbnail or geometry, hash-checked against the wardrobe."""
    root,data=_wardrobe(root)
    garment=next((g for g in data['garments'] if g['id']==ident),None)
    if garment is None or not garment.get(field):return None
    return _wardrobe_asset(root,garment[field],garment.get(field+'_sha256'))
