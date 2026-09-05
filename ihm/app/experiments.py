"""Read body-attached experiment artifacts without silently accepting stale code."""
from pathlib import Path
import hashlib
import json

EXPERIMENTS={'forearm-touch':'forearm-touch.json','skin-transport':'skin-transport.json','skin-electric':'skin-electric.json','spectra':'regional-spectra.json','details':'details.json'}


def read_experiment(root,kind):
    if kind not in EXPERIMENTS:raise ValueError('Unknown body experiment')
    root=Path(root).resolve()
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
