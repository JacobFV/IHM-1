"""Reviewed immutable IBM default; explicit None retains the legacy replay law."""
from pathlib import Path
import json,re
from .candidate import SourcePin,verify_pin

DEFAULT_SOURCE=object()
ACTIVE_COMMIT='398375cc993ac17907d4ddd9f6d53eed2fab31de'
ACTIVE_MANIFEST_SHA256='397b3fd3e20b7a15da0055bbba9c83fec6cfaee04135e34c99407da8280d300b'
ACTIVE_PACKAGE_SHA256='f969da521ef87ece047da189de17118538a650ebab65605a3254a6be38dffa60'

def resolve_ibm_candidate(root,selector):
    """Resolve only server-installed commit directories; never import candidate code."""
    if (not isinstance(selector,dict) or set(selector)!={'commit','manifest_sha256'}
        or not isinstance(selector['commit'],str) or not re.fullmatch('[0-9a-f]{40}',selector['commit'])
        or not isinstance(selector['manifest_sha256'],str) or not re.fullmatch('[0-9a-f]{64}',selector['manifest_sha256'])):
        raise ValueError('ibm_candidate requires exact commit and manifest_sha256 identities')
    root=Path(root).absolute()
    candidate=root/'data/derived/ibm-candidates'/selector['commit']
    try:
        # Check before resolve: SourcePin canonicalizes paths and would hide links.
        if any(path.is_symlink() for path in (candidate,*candidate.parents)):
            raise ValueError('Candidate directory ancestors may not be symlinks')
        entries=list(candidate.rglob('*'))
        if any(path.is_symlink() or not (path.is_file() or path.is_dir()) for path in entries):
            raise ValueError('Candidate entries must be ordinary files and directories')
        if {path.name for path in candidate.iterdir()}!={'source','manifest.json','source_pin.json'}:
            raise ValueError('Candidate root has missing or extra entries')
        raw=json.loads((candidate/'source_pin.json').read_bytes())
        if raw.get('artifact_dir')!=str(candidate) or raw.get('manifest_sha256')!=selector['manifest_sha256']:
            raise ValueError('Candidate pin does not match server location and requested identity')
        pin=SourcePin.load(candidate/'source_pin.json')
        if raw!=pin.to_dict():raise ValueError('Candidate pin has unexpected fields')
        manifest,_=verify_pin(pin)
        if manifest.get('donor_commit')!=selector['commit']:
            raise ValueError('Candidate commit does not match requested identity')
        expected={Path('source'),Path('manifest.json'),Path('source_pin.json')}
        for name in manifest['files']:
            relative=Path('source')/name
            expected.add(relative);expected.update(p for p in relative.parents if p!=Path('.'))
        if {path.relative_to(candidate) for path in entries}!=expected:
            raise ValueError('Candidate contains unexpected directories or files')
        return pin
    except (OSError,KeyError,TypeError,ValueError) as error:
        raise ValueError('Invalid server-owned IBM candidate: '+str(error)) from error


def resolve_source(root,selection=DEFAULT_SOURCE):
    if selection is None:return None
    if selection is DEFAULT_SOURCE:
        pin=resolve_ibm_candidate(root,{'commit':ACTIVE_COMMIT,'manifest_sha256':ACTIVE_MANIFEST_SHA256})
        if pin.package_sha256!=ACTIVE_PACKAGE_SHA256:raise ValueError('Active IBM package identity differs')
        return pin
    if not isinstance(selection,SourcePin):raise ValueError('Explicit SourcePin or None legacy selection required')
    verify_pin(selection)
    return selection
