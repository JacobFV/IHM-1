"""Publish completed local research artifacts to the same canonical-body inspector."""
import argparse
import hashlib
import json
from pathlib import Path
from ihm.assembly.systemic_projection import project_systemic, accepted_systemic_sources
from ihm.assembly.anatomy import write_json
from ihm.native import _sha

p=argparse.ArgumentParser();p.add_argument('directories',nargs='+');args=p.parse_args()
root=Path(__file__).resolve().parents[1]
canonical=root/'data/derived/canonical'
index_path=canonical/'systemic-index.json'
index=json.loads(index_path.read_text()) if index_path.exists() else {'runs':[]}
by_id={run['id']:run for run in index['runs']}
for directory in args.directories:
    source=Path(directory)/'systemic.json'
    acceptance=accepted_systemic_sources(root,source)
    display=project_systemic(root,source)
    display['runtime_sources'].update(acceptance)
    display['acceptance']='Matched causal execution and local-store checks passed; not experimental human validation'
    id=display['configuration']['protocol']
    # Immutable records leave the old index usable if a later projection fails.
    digest=hashlib.sha256(json.dumps(display,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    target=canonical/'systemic'/f'{id}-{digest}.json'
    write_json(target,display)
    by_id[id]=dict(id=id,label=id.replace('_',' + ').title(),duration_s=display['clock']['end_s']-display['clock']['start_s'],
                  sample_interval_s=display['clock']['sample_interval_s'],has_body_projection=display['has_body_projection'],
                  path=str(target.relative_to(root)),sha256=_sha(target))
    print(json.dumps(by_id[id]),flush=True)
write_json(index_path,{'runs':list(by_id.values())})
