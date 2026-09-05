"""Historical inputs remain reproducible; changed snapshot bytes fail closed."""
import json
from pathlib import Path
import tempfile
from ihm.native import _sha
from ihm.assembly.systemic_evidence import freeze_sources,resolve_sources

with tempfile.TemporaryDirectory() as directory:
    root=Path(directory);source=root/'solver.cpp';source.write_text('original numerical source\n')
    sources={'solver.cpp':_sha(source)};run=root/'experiment';run.mkdir()
    frozen=freeze_sources(root,run,sources)
    source.write_text('a future numerical source\n')
    assert resolve_sources(root,run,sources)==frozen
    assert freeze_sources(root,run,sources)==frozen
    try:freeze_sources(root,run,{'solver.cpp':_sha(source)})
    except ValueError:pass
    else:raise AssertionError('Existing experiment was repinned')
    copy=run/'inputs/solver.cpp';copy.write_text('tampered\n')
    try:resolve_sources(root,run,sources)
    except ValueError:pass
    else:raise AssertionError('Tampered historical source accepted')
    unknown=root/'unarchived';unknown.mkdir()
    try:resolve_sources(root,unknown,sources)
    except ValueError:pass
    else:raise AssertionError('Unarchived stale input accepted')
print('PASS immutable systemic inputs and version-independent historical replay')
