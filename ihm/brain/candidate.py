"""Immutable source-only IBM commit capture and explicit shared package/law pins."""
from dataclasses import dataclass
import ast
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tarfile

SYMBOLS=('_sigmoid','wilson_cowan_excitatory','shunting_inhibition_rate')
NEURAL='ibm/processes/neural.py'
MAX_BYTES=32*1024*1024


def _sha(raw):return hashlib.sha256(raw).hexdigest()

def _encoded(value):return (json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()


def _laws(raw):
    tree=ast.parse(raw)
    nodes={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in SYMBOLS}
    if set(nodes)!=set(SYMBOLS):raise ValueError('Candidate lacks required regional neural symbols')
    return {name:_sha(ast.dump(nodes[name],include_attributes=False).encode()) for name in SYMBOLS}


@dataclass(frozen=True)
class SourcePin:
    artifact_dir: Path
    package_sha256: str
    neural_source_sha256: str
    manifest_sha256: str

    def __post_init__(self):
        object.__setattr__(self,'artifact_dir',Path(self.artifact_dir).resolve())
        for digest in (self.package_sha256,self.neural_source_sha256,self.manifest_sha256):
            if not isinstance(digest,str) or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
                raise ValueError('Explicit SHA256 source pin required')

    def to_dict(self):
        return {'schema':'ihm.ibm-source-pin.v1','artifact_dir':str(self.artifact_dir),
                'package_sha256':self.package_sha256,'neural_source_sha256':self.neural_source_sha256,
                'manifest_sha256':self.manifest_sha256}

    @classmethod
    def load(cls,path):
        value=json.loads(Path(path).read_bytes())
        if value.pop('schema',None)!='ihm.ibm-source-pin.v1':raise ValueError('Invalid source pin schema')
        return cls(**value)


def verify_pin(pin):
    """Return verified manifest and captured bytes without importing donor modules."""
    if not isinstance(pin,SourcePin):raise ValueError('SourcePin required for opt-in import')
    root=pin.artifact_dir;source=root/'source'
    raw=(root/'manifest.json').read_bytes()
    if _sha(raw)!=pin.manifest_sha256:raise ValueError('Candidate manifest changed')
    manifest=json.loads(raw)
    if manifest.get('schema')!='ihm.ibm-candidate.v1' or manifest.get('evidence_policy')!='compiled_source_priors_only':
        raise ValueError('Unrecognized candidate/evidence policy')
    # The donor's optional summary loaders walk upward for data/sources + ibm.
    # Refuse such ambient evidence rather than silently changing prior origin.
    if any((p/'data/sources').is_dir() and (p/'ibm').is_dir() for p in (source,*source.parents)):
        raise ValueError('Candidate may not inherit an ambient IBM evidence root')
    actual={p.relative_to(source).as_posix() for p in source.rglob('*') if p.is_file()}
    if actual!=set(manifest['files']):raise ValueError('Candidate contains missing/unmanifested files')
    snapshots={}
    for relative,digest in manifest['files'].items():
        path=source/relative
        if path.is_symlink() or not path.resolve().is_relative_to(source.resolve()):raise ValueError('Candidate source path escapes snapshot')
        snapshots[relative]=path.read_bytes()
        if _sha(snapshots[relative])!=digest:raise ValueError('Candidate source bytes changed: '+relative)
    package=_sha(json.dumps(manifest['files'],sort_keys=True).encode())
    if package!=pin.package_sha256 or package!=manifest['package_sha256'] or _sha(snapshots[NEURAL])!=pin.neural_source_sha256:
        raise ValueError('Shared package/regional neural pin mismatch')
    if _laws(snapshots[NEURAL])!=manifest['regional_law']['ast_sha256']:raise ValueError('Candidate regional law receipt differs')
    return manifest,snapshots


def capture_candidate(donor,destination,*,revision='HEAD',baseline_neural=None):
    """Capture committed ibm Python bytes; never modifies donor or existing output.

    Uncommitted edits are reported as metadata, never included. No datasets,
    learned parameters, training scripts or checkpoints are archived.
    """
    donor=Path(donor).resolve();destination=Path(destination).resolve()
    if destination.exists() or destination.is_relative_to(donor):raise ValueError('Fresh candidate outside donor checkout required')
    if not isinstance(revision,str) or not revision or revision.startswith('-'):raise ValueError('Invalid donor revision')
    def git(*args):return subprocess.check_output(['git','-C',str(donor),*args],stderr=subprocess.PIPE)
    commit=git('rev-parse','--verify',revision+'^{commit}').decode().strip()
    archive=git('archive','--format=tar',commit,'--','ibm')
    if len(archive)>MAX_BYTES:raise ValueError('Source archive exceeds bounded capture size')
    files={}
    with tarfile.open(fileobj=io.BytesIO(archive),mode='r:') as tar:
        for item in tar:
            if item.isdir():continue
            name=PurePosixPath(item.name)
            if not item.isfile() or name.is_absolute() or '..' in name.parts or name.parts[0]!='ibm' or name.suffix!='.py':
                raise ValueError('Candidate archive must contain ordinary IBM Python source only')
            if item.name in files:raise ValueError('Duplicate source archive member')
            files[item.name]=tar.extractfile(item).read()
    if 'ibm/__init__.py' not in files or NEURAL not in files:raise ValueError('Incomplete IBM source archive')
    hashes={name:_sha(raw) for name,raw in sorted(files.items())};laws=_laws(files[NEURAL])
    baseline=None
    if baseline_neural is not None:
        old=Path(baseline_neural).read_bytes();old_laws=_laws(old)
        baseline={'source_sha256':_sha(old),'ast_sha256':old_laws,
                  'changed_symbols':[name for name in SYMBOLS if laws[name]!=old_laws[name]],
                  'acceptance':'source delta only; changed law requires explicit opt-in and bounded behavior acceptance'}
    manifest={'schema':'ihm.ibm-candidate.v1','donor_commit':commit,'donor_path':str(donor),
        'donor_status':git('status','--porcelain','--untracked-files=no').decode().strip(),
        'capture':'git archive of exact commit; working-tree changes excluded','files':hashes,
        'package_sha256':_sha(json.dumps(hashes,sort_keys=True).encode()),
        'evidence_policy':'compiled_source_priors_only','training_or_datasets_included':False,
        'regional_law':{'source_path':NEURAL,'source_sha256':hashes[NEURAL],'ast_sha256':laws,'baseline_comparison':baseline},
        'promotion_status':'candidate_only; active source unchanged'}
    destination.mkdir(parents=True,exist_ok=False)
    try:
        for name,raw in files.items():
            path=destination/'source'/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
        raw=_encoded(manifest);(destination/'manifest.json').write_bytes(raw)
        pin=SourcePin(destination,manifest['package_sha256'],hashes[NEURAL],_sha(raw))
        (destination/'source_pin.json').write_bytes(_encoded(pin.to_dict()))
        verify_pin(pin)
        return pin
    except BaseException:
        shutil.rmtree(destination)
        raise
