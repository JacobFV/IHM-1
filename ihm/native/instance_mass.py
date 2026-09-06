"""Read-only selection of the opt-in Simbody State-mass library and adapter."""
from pathlib import Path
import hashlib,json

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def variant_identity(root,path):
    root=Path(root).resolve();variant=(root/path).resolve()
    if not variant.is_relative_to(root/'data/runtime/opensim/variants'):raise ValueError('Instance mass variant must be an owned opt-in variant')
    manifest_path=variant/'manifest.json';manifest=json.loads(manifest_path.read_text());library=variant/'libSimTKsimbody.so.3.9'
    if manifest.get('schema')!='ihm.simbody-instance-mass-variant.v1' or manifest.get('complete') is not True or sha(library)!=manifest.get('library_sha256'):
        raise ValueError('Incomplete or changed native instance mass variant')
    return {'path':str(variant.relative_to(root)),'manifest_sha256':sha(manifest_path),'library_sha256':manifest['library_sha256'],'library_path':str(library.relative_to(root))}
def pointer_directory(root,identity):return Path(root)/'data/runtime/mechanical-stream/instance-mass'/identity['library_sha256']
