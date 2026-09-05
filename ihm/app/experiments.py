"""Read body-attached experiment artifacts without silently accepting stale code."""
from pathlib import Path
import hashlib
import json

EXPERIMENTS={'forearm-touch':'forearm-touch.json','skin-transport':'skin-transport.json','skin-electric':'skin-electric.json','spectra':'regional-spectra.json','details':'details.json'}


def read_experiment(root,kind):
    root=Path(root).resolve()
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
    if kind not in EXPERIMENTS:raise ValueError('Unknown body experiment')
    if kind=='spectra':
        # An unchanged record may itself have stale model/anatomical inputs.
        for dependency in ('forearm-touch','skin-electric'):
            read_experiment(root,dependency)
    data=json.loads((root/'data/derived/canonical'/EXPERIMENTS[kind]).read_bytes())
    hashes={**data.get('runtime_sources',{}),**data.get('source_hashes',{})}
    hashes.update({r['path']:r['sha256'] for r in data.get('model',{}).get('provenance',[])})
    for path,expected in hashes.items():
        file=(root/path).resolve()
        if not file.is_relative_to(root) or hashlib.sha256(file.read_bytes()).hexdigest()!=expected:
            raise ValueError('Experiment source changed; rebuild '+kind)
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
