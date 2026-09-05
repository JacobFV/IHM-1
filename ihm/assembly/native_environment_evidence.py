"""Detached native input archives; external paths are identities, never replay paths.

Resource capture occurs after native startup. Without a prior resource inventory,
matching startup bytes cannot be inferred retrospectively from these snapshots.
"""
from pathlib import Path
import datetime,hashlib,json,os,shutil,stat,struct,tempfile

RESOURCES=('patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml')
SCHEMA='ihm.native-environment.v1'

def _sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()

def _paths(root,native_directory):
    root=Path(root).resolve();native=Path(native_directory)
    native=(native if native.is_absolute() else root/native).resolve()
    if native==root or not native.is_relative_to(root) or not native.is_dir():
        raise ValueError('Native experiment directory must exist inside the workspace')
    return root,native

def _digest(value):
    if not isinstance(value,str) or len(value)!=64 or any(c not in '0123456789abcdef' for c in value):
        raise ValueError('Expected SHA-256 identity')
    return value

def _interpreter(executable):
    """Read PT_INTERP without executing an executable or trusting a shell parser."""
    with Path(executable).open('rb') as stream:
        header=stream.read(64)
        if len(header)<64 or header[:4]!=b'\x7fELF':raise ValueError('Native executable must be ELF')
        endian={1:'<',2:'>'}.get(header[5]);width=header[4]
        if endian is None or width not in (1,2):raise ValueError('Unsupported ELF encoding')
        phoff=struct.unpack_from(endian+('Q' if width==2 else 'I'),header,32 if width==2 else 28)[0]
        entsize,count=struct.unpack_from(endian+'HH',header,54 if width==2 else 42)
        if entsize<(56 if width==2 else 32) or count>4096:raise ValueError('Invalid ELF program headers')
        for i in range(count):
            stream.seek(phoff+i*entsize);entry=stream.read(entsize)
            if len(entry)!=entsize:raise ValueError('Truncated ELF program header')
            if struct.unpack_from(endian+'I',entry)[0]!=3:continue
            offset=struct.unpack_from(endian+('Q' if width==2 else 'I'),entry,8 if width==2 else 4)[0]
            size=struct.unpack_from(endian+('Q' if width==2 else 'I'),entry,32 if width==2 else 16)[0]
            if not 1<size<4096:raise ValueError('Invalid ELF interpreter length')
            stream.seek(offset);value=stream.read(size)
            if not value.endswith(b'\0'):raise ValueError('Invalid ELF interpreter string')
            path=Path(value[:-1].decode())
            if not path.is_absolute():raise ValueError('ELF interpreter must have an absolute identity')
            return path
    return None

def freeze_native_environment(root,native_directory):
    """Archive startup-hashed dependencies and current runtime resources once.

    Return only workspace-relative immutable input paths and SHA-256 digests.
    Existing successful archives are resolved, never repinned to current inputs.
    Failed attempts remain in distinct ``.environment-capture-*`` directories.
    """
    root,native=_paths(root,native_directory);final=native/'environment-inputs'
    if final.exists() or final.is_symlink():return resolve_native_environment(root,native)
    staging=Path(tempfile.mkdtemp(prefix='.environment-capture-',dir=native))
    began=datetime.datetime.now(datetime.timezone.utc).isoformat()
    try:
        native_manifest=native/'manifest.json';manifest_hash=_sha(native_manifest)
        manifest=json.loads(native_manifest.read_text())
        if manifest.get('schema')!='ihm.native-session.v1' or not manifest.get('dependency_sha256'):
            raise ValueError('Native startup dependency manifest required')
        blobs={}
        retained=None
        def input_version(path,expected):
            nonlocal retained
            path=Path(path).resolve();_digest(expected)
            if path.is_file() and _sha(path)==expected:return path
            if retained is None:
                retained={};prior=native.parent/'frozen-sources.json'
                if prior.exists():
                    from .systemic_evidence import resolve_sources
                    identities=json.loads(prior.read_text())['original_sources']
                    retained=resolve_sources(root,native.parent,identities)
            for name,digest in retained.items():
                if digest==expected:return root/name
            raise ValueError('Startup input changed and no validated historical copy exists: '+str(path))
        def snapshot(path,expected=None):
            path=Path(path).resolve()
            before=path.stat()
            if not stat.S_ISREG(before.st_mode):raise ValueError('Only regular input files can be archived')
            digest=_sha(path)
            if expected is not None and digest!=_digest(expected):raise ValueError('Input differs from native startup identity: '+str(path))
            relative='blobs/'+digest;target=staging/relative
            if not target.exists():
                target.parent.mkdir(exist_ok=True);shutil.copyfile(path,target)
            after=path.stat()
            identity=lambda s:(s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)
            if identity(before)!=identity(after) or _sha(path)!=digest or _sha(target)!=digest:
                raise ValueError('Input changed during native environment capture: '+str(path))
            blobs[relative]=digest
            return {'captured_from_absolute_path':str(path),'sha256':digest,'bytes':before.st_size,
                    'mode':stat.S_IMODE(before.st_mode),'blob':relative}
        dependencies=[]
        for name,digest in sorted(manifest['dependency_sha256'].items()):
            path=Path(name)
            if not path.is_absolute():raise ValueError('Native dependency identity must be absolute')
            dependencies.append({'startup_identity_path':name,**snapshot(input_version(path,digest),digest)})
        executable_path=Path(manifest['command'][0]).resolve()
        if not executable_path.is_relative_to(root):raise ValueError('Native executable must belong to the workspace')
        executable={'startup_identity_path':manifest['command'][0],**snapshot(input_version(executable_path,manifest['executable_sha256']),manifest['executable_sha256'])}
        interpreter_path=_interpreter(staging/executable['blob'])
        interpreter=snapshot(interpreter_path) if interpreter_path else None
        if interpreter:interpreter['identity_basis']='ELF PT_INTERP; bytes captured now, absent from legacy startup ldd identity table'
        resources=[]
        def inventory(logical,path,ancestors=()):
            resolved=path.resolve()
            if not resolved.is_relative_to(root):raise ValueError('Runtime resource escapes the workspace: '+logical)
            if resolved in ancestors:raise ValueError('Runtime resource symlink cycle: '+logical)
            mode=resolved.stat().st_mode
            if stat.S_ISDIR(mode):
                resources.append({'kind':'directory','logical_path':logical,'original_absolute_path':str(resolved),
                                  'mode':stat.S_IMODE(mode)})
                for child in sorted(resolved.iterdir()):inventory(logical+'/'+child.name,child,(*ancestors,resolved))
            elif stat.S_ISREG(mode):resources.append({'kind':'file','logical_path':logical,**snapshot(resolved)})
            else:raise ValueError('Special runtime resource is not a numerical input: '+logical)
        for name in RESOURCES:inventory(name,native/name)
        # Inventory identities include content hashes; a second walk catches additions,
        # deletion and mutation while copying, without walking logs, output or states.
        first=resources;resources=[]
        for name in RESOURCES:inventory(name,native/name)
        if first!=resources:raise ValueError('Runtime resource tree changed during capture')
        if _sha(native_manifest)!=manifest_hash:raise ValueError('Native startup manifest changed during capture')
        journal=native/'receipts.jsonl';prefix=journal.read_bytes() if journal.exists() else b''
        records=[json.loads(line) for line in prefix.splitlines() if line.strip()]
        commands=[r for r in records if 'command' in r]
        capture={'started_utc':began,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'stage':'post_start_after_commands' if commands else 'post_start_before_first_action',
            'command_count_at_capture':len(commands),'journal_prefix_bytes':len(prefix),
            'journal_prefix_sha256':hashlib.sha256(prefix).hexdigest(),
            'resources_verified_at_process_start':False,
            'interpretation':'Dependencies match startup hashes. Resource bytes are observed at capture time; their identity during earlier initialization/actions is not retroactively established.'}
        receipt={'schema':SCHEMA,'native_directory':str(native.relative_to(root)),
            'native_manifest_sha256':manifest_hash,'capture':capture,'dependencies':dependencies,
            'executable':executable,'elf_interpreter':interpreter,'resources':resources,'blobs':blobs,
            'resource_roots':list(RESOURCES),'capturer_sha256':_sha(__file__),
            'limitations':['Archive bytes do not establish checkpoint completeness, kernel/CPU equivalence, or bitwise replay.',
                'Dynamic dlopen inputs and environment variables are not reconstructed by the legacy ldd manifest.',
                'System libraries retain their original licensing; this local research archive grants no redistribution rights.']}
        (staging/'manifest.json').write_text(json.dumps(receipt,indent=2)+'\n')
        # Rename a newly created directory only; never merge or repair an existing archive.
        if final.exists() or final.is_symlink():raise ValueError('Another capture already owns this archive')
        staging.rename(final)
        return resolve_native_environment(root,native)
    except BaseException as error:
        if staging.exists():
            (staging/'failure.json').write_text(json.dumps({'started_utc':began,'error_type':type(error).__name__,'error':str(error)},indent=2)+'\n')
        raise

def resolve_native_environment(root,native_directory):
    """Validate archive closure without opening its original external identities."""
    root,native=_paths(root,native_directory);archive=native/'environment-inputs';receipt_path=archive/'manifest.json'
    if archive.is_symlink() or not archive.is_dir() or receipt_path.is_symlink():raise ValueError('Native environment archive is missing or indirect')
    receipt=json.loads(receipt_path.read_text())
    if receipt.get('schema')!=SCHEMA or receipt.get('native_directory')!=str(native.relative_to(root)):
        raise ValueError('Native environment identity mismatch')
    original=native/'manifest.json'
    if original.is_symlink() or _sha(original)!=receipt['native_manifest_sha256']:raise ValueError('Native startup manifest changed')
    startup=json.loads(original.read_text())
    dependencies={entry['startup_identity_path']:entry['sha256'] for entry in receipt['dependencies']}
    if len(dependencies)!=len(receipt['dependencies']) or dependencies!=startup['dependency_sha256']:
        raise ValueError('Archive must retain every startup dependency identity')
    if receipt['executable']['startup_identity_path']!=startup['command'][0] or receipt['executable']['sha256']!=startup['executable_sha256']:
        raise ValueError('Archive native executable identity differs from startup')
    result={str(original.relative_to(root)):_sha(original),str(receipt_path.relative_to(root)):_sha(receipt_path)}
    references=[*receipt['dependencies'],receipt['executable'],*[r for r in receipt['resources'] if r['kind']=='file']]
    if receipt['elf_interpreter'] is not None:references.append(receipt['elf_interpreter'])
    if any(r['blob']!='blobs/'+_digest(r['sha256']) for r in references):raise ValueError('Conflicting archive blob identities')
    expected={r['blob']:_digest(r['sha256']) for r in references}
    if expected!=receipt['blobs'] or set(receipt['resource_roots'])!=set(RESOURCES):raise ValueError('Incomplete environment archive identities')
    for name,digest in expected.items():
        if name!='blobs/'+digest:raise ValueError('Invalid content-addressed input path')
        path=archive/name
        if path.is_symlink() or not path.resolve().is_relative_to(archive.resolve()) or not path.is_file() or path.stat().st_nlink!=1 or _sha(path)!=digest:
            raise ValueError('Retained native environment input changed: '+name)
        result[str(path.relative_to(root))]=digest
    actual={str(p.relative_to(archive)) for p in archive.rglob('*') if p.is_file()}
    if actual!=set(expected)|{'manifest.json'}:raise ValueError('Unexpected files in immutable native environment archive')
    return result
