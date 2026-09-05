"""Regional evidence must bind the selected workspace and transitive inputs."""
from pathlib import Path
import hashlib
import json
import tempfile
from ihm.app.experiments import read_experiment
from ihm.assembly.regional_touch import run_touch

with tempfile.TemporaryDirectory() as td:
    root=Path(td)
    try: run_touch(root,divisions=(2,2,2))
    except FileNotFoundError as error: assert 'IBM artifact missing' in str(error),str(error)
    else: raise AssertionError('A workspace without an IBM import used another workspace')
    directory=root/'data/derived/canonical';directory.mkdir(parents=True)
    source=root/'mechanics.py';source.write_text('original')
    touch=directory/'forearm-touch.json'
    touch.write_text(json.dumps({'runtime_sources':{'mechanics.py':hashlib.sha256(source.read_bytes()).hexdigest()}}))
    electric=directory/'skin-electric.json';electric.write_text('{}')
    spectra=directory/'regional-spectra.json'
    spectra.write_text(json.dumps({'source_hashes':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (touch,electric)}}))
    read_experiment(root,'spectra')
    source.write_text('changed')
    for kind in ('forearm-touch','spectra'):
        try: read_experiment(root,kind)
        except ValueError as error: assert 'source changed' in str(error)
        else: raise AssertionError('Stale transitive evidence accepted: '+kind)
print('PASS: workspace-specific IBM requirement and transitive spectrum dependency rejection')
