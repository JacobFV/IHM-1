"""Detached native input archives; external paths are identities, never replay paths.

Capture phase is explicit. Historical snapshots cannot retroactively establish
startup bytes; prestart snapshots can supply a detached resource tree to a future
process, whose actual launch/use must be recorded separately by its owner.
"""
from pathlib import Path,PurePosixPath
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

def freeze_native_environment(root,native_directory,*,before_start=False):
    """Archive startup-hashed dependencies and current runtime resources once.

    Return only workspace-relative immutable input paths and SHA-256 digests.
    Existing successful archives are resolved, never repinned to current inputs.
    Failed attempts remain in distinct ``.environment-capture-*`` directories.
    """
    if not isinstance(before_start,bool):raise ValueError('Capture phase must be an explicit boolean')
    root,native=_paths(root,native_directory);final=native/'environment-inputs'
    if final.exists() or final.is_symlink():return _resolve_native_archive(root,native)
    staging=Path(tempfile.mkdtemp(prefix='.environment-capture-',dir=native))
    began=datetime.datetime.now(datetime.timezone.utc).isoformat()
    try:
        def check_prestart():
            journal=native/'receipts.jsonl'
            if before_start and ((journal.exists() and journal.read_bytes().strip()) or (native/'runner_stdout.log').exists()):
                raise ValueError('Prestart capture requires no prior process output or action receipts')
        check_prestart()
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
        check_prestart()
        capture={'started_utc':began,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'stage':'pre_start' if before_start else 'post_start_after_commands' if commands else 'post_start_before_first_action',
            'command_count_at_capture':len(commands),'journal_prefix_bytes':len(prefix),
            'journal_prefix_sha256':hashlib.sha256(prefix).hexdigest(),
            'resources_verified_at_process_start':False,
            'interpretation':('Captured at the caller-declared boundary before Popen. A validated detached resource tree must be selected for startup; this capture does not attest subsequent process consumption.' if before_start else
                'Dependencies match startup hashes. Resource bytes are observed at capture time; their identity during earlier initialization/actions is not retroactively established.')}
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
        return _resolve_native_archive(root,native)
    except BaseException as error:
        if staging.exists():
            (staging/'failure.json').write_text(json.dumps({'started_utc':began,'error_type':type(error).__name__,'error':str(error)},indent=2)+'\n')
        raise

def _resolve_native_archive(root,native_directory):
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


def _resource_entries(receipt):
    """Validate a complete logical tree before any destination paths are written."""
    entries={}
    for entry in receipt['resources']:
        logical=entry.get('logical_path')
        if not isinstance(logical,str) or not logical or '\\' in logical:
            raise ValueError('Invalid resource logical path')
        path=PurePosixPath(logical)
        if path.is_absolute() or not path.parts or path.as_posix()!=logical or any(part in ('.','..') for part in path.parts) or path.parts[0] not in RESOURCES:
            raise ValueError('Unsafe resource logical path: '+logical)
        if logical in entries or entry.get('kind') not in ('file','directory'):
            raise ValueError('Duplicate or unsupported resource entry: '+logical)
        mode=entry.get('mode')
        if isinstance(mode,bool) or not isinstance(mode,int) or not 0<=mode<=0o7777:
            raise ValueError('Invalid resource mode: '+logical)
        entries[logical]=entry
    for name in RESOURCES:
        expected='file' if name in ('UCEDefs.conf','BioGearsConfiguration.xml') else 'directory'
        if entries.get(name,{}).get('kind')!=expected:raise ValueError('Missing resource root: '+name)
    for logical in entries:
        for parent in PurePosixPath(logical).parents:
            if str(parent)=='.':continue
            if entries.get(str(parent),{}).get('kind')!='directory':raise ValueError('Resource parent is not a retained directory: '+logical)
    return entries


def _owned_regular(path):
    if path.is_symlink() or not path.is_file():
        raise ValueError('Selected native input is missing or indirect: '+str(path))
    info=path.stat()
    if info.st_uid!=os.geteuid() or info.st_nlink!=1:
        raise ValueError('Selected native input is not independently owned: '+str(path))
    return _sha(path)


def _absolute_identity(value):
    """Parse recorded identities lexically; never inspect an old absolute path."""
    if not isinstance(value,str) or '\\' in value or '\0' in value:
        raise ValueError('Invalid absolute native input identity')
    path=PurePosixPath(value)
    if not path.is_absolute() or path.as_posix()!=value or any(p in ('.','..') for p in path.parts) or value.startswith('//'):
        raise ValueError('Unsafe absolute native input identity')
    return path


def resolve_native_execution_inputs(root,native_directory):
    """Validate pre-Popen input selection against detached files, not old paths.

    A selection receipt is evidence of selected inputs, not process tracing or
    proof that initialization completed. Relocation retains native_directory's
    workspace-relative layout and maps recorded native-local absolute identities.
    """
    root,native=_paths(root,native_directory)
    _resolve_native_archive(root,native)
    archive_manifest=native/'environment-inputs/manifest.json'
    archive=json.loads(archive_manifest.read_text())
    if archive.get('capture',{}).get('stage')!='pre_start':
        raise ValueError('Execution selection is required only for explicit prestart captures')
    launch_path=native/'execution-inputs.json';launch_hash=_owned_regular(launch_path)
    launch=json.loads(launch_path.read_text());startup=json.loads((native/'manifest.json').read_text())
    entries=_resource_entries(archive);tree=native/'runtime-resources'
    _validate_materialized_tree(tree,entries,_sha(archive_manifest))
    for key,path in (('native_manifest_sha256',native/'manifest.json'),
                     ('environment_manifest_sha256',archive_manifest),
                     ('resource_materialization_sha256',tree/'materialization.json')):
        if _digest(launch.get(key))!=_sha(path):raise ValueError('Execution selection receipt differs from '+key)
    recorded_tree=_absolute_identity(launch.get('selected_resource_tree'))
    if recorded_tree.name!='runtime-resources':raise ValueError('Selected resource tree is not native-local')
    recorded_native=recorded_tree.parent
    relative=PurePosixPath(archive['native_directory'])
    if tuple(recorded_native.parts[-len(relative.parts):])!=relative.parts:
        raise ValueError('Selected resource tree has a different native directory identity')
    # Existing absolute root links remain verifiable after moving an archive to
    # another workspace. Their text is checked; their old targets are never read.
    for name in RESOURCES:
        link=native/name
        if not link.is_symlink():raise ValueError('Selected native resource root is missing or not a selection link: '+name)
        target=os.readlink(link)
        if target not in (str(recorded_tree/name),str(tree/name),str(PurePosixPath('runtime-resources')/name)):
            raise ValueError('Native resource selection escapes the detached tree: '+name)
    command=startup.get('command',[]);configuration=startup.get('configuration',{})
    patient=configuration.get('patient')
    if not isinstance(patient,str) or not patient or any(c in patient for c in '/\\\0') or patient in ('.','..'):
        raise ValueError('Invalid selected patient identity')
    if len(command)!=4 or command[1]!=patient:
        raise ValueError('Native command does not match selected patient')
    selected=_absolute_identity(launch.get('selected_patient_input'))
    initial=startup.get('initialization_patient_used')
    if initial is True:
        if configuration.get('state_path') is not None or startup.get('state_sha256') is not None or command[2]!='-':
            raise ValueError('Patient initialization conflicts with state selection')
        local=PurePosixPath('runtime-resources')/'patients'/(patient+'.xml')
        expected_hash=startup.get('patient_identity_input_sha256')
    elif initial is False:
        local=PurePosixPath('input-state.xml')
        if not configuration.get('state_path') or _absolute_identity(command[2])!=recorded_native/local:
            raise ValueError('Native command does not select the retained state')
        expected_hash=startup.get('state_sha256')
        if expected_hash!=startup.get('patient_identity_input_sha256'):
            raise ValueError('State and patient input identities disagree')
    else:raise ValueError('Native patient initialization mode is missing')
    if selected!=recorded_native/local:raise ValueError('Selected patient input is not the declared native-local input')
    patient_path=native/str(local)
    if _digest(launch.get('selected_patient_input_sha256'))!=_digest(expected_hash) or _owned_regular(patient_path)!=expected_hash:
        raise ValueError('Selected patient input bytes differ from startup')
    result={str(launch_path.relative_to(root)):launch_hash,
            str((tree/'materialization.json').relative_to(root)):_sha(tree/'materialization.json'),
            str(patient_path.relative_to(root)):expected_hash}
    for name,entry in entries.items():
        if entry['kind']=='file':result[str((tree/name).relative_to(root))]=entry['sha256']
    return result


def resolve_native_environment(root,native_directory):
    """Resolve retained inputs; new prestart runs also require execution selection."""
    root,native=_paths(root,native_directory)
    result=_resolve_native_archive(root,native)
    receipt=json.loads((native/'environment-inputs/manifest.json').read_text())
    phase=receipt.get('capture',{}).get('stage')
    if phase not in ('pre_start','post_start_before_first_action','post_start_after_commands'):
        raise ValueError('Native environment capture phase is missing or unknown')
    if phase=='pre_start':
        result.update(resolve_native_execution_inputs(root,native))
    return result


def _validate_materialized_tree(tree,entries,archive_manifest_hash):
    """Inspect regular owned files without following any directory or file links."""
    if tree.is_symlink() or not tree.is_dir() or tree.stat().st_uid!=os.geteuid():
        raise ValueError('Detached native resource tree is missing, indirect, or not owned')
    receipt_path=tree/'materialization.json'
    if receipt_path.is_symlink() or not receipt_path.is_file() or receipt_path.stat().st_nlink!=1 or receipt_path.stat().st_uid!=os.geteuid():
        raise ValueError('Detached native resource receipt is indirect or shared')
    record=json.loads(receipt_path.read_text())
    expected_files={name:entry['sha256'] for name,entry in entries.items() if entry['kind']=='file'}
    if record!={'schema':'ihm.native-runtime-resources.v1','archive_manifest_sha256':archive_manifest_hash,
                'resource_file_sha256':expected_files,'resource_roots':list(RESOURCES)}:
        raise ValueError('Detached native resource identity changed')
    observed=set()
    def walk(directory):
        for path in directory.iterdir():
            name=str(path.relative_to(tree));mode=path.lstat()
            if stat.S_ISLNK(mode.st_mode) or mode.st_uid!=os.geteuid():raise ValueError('Indirect or nonowned native resource: '+name)
            if name=='materialization.json':continue
            entry=entries.get(name)
            if entry is None:raise ValueError('Unexpected native resource: '+name)
            observed.add(name)
            if stat.S_IMODE(mode.st_mode)!=(entry['mode']&0o777):raise ValueError('Native resource permissions changed: '+name)
            if entry['kind']=='directory':
                if not stat.S_ISDIR(mode.st_mode):raise ValueError('Native resource directory changed: '+name)
                walk(path)
            elif not stat.S_ISREG(mode.st_mode) or mode.st_nlink!=1 or _sha(path)!=entry['sha256']:
                raise ValueError('Detached native resource bytes or ownership changed: '+name)
    walk(tree)
    if observed!=set(entries):raise ValueError('Detached native resource tree is incomplete')
    return tree


def materialize_native_resources(root,native_directory):
    """Return an absolute Path containing validated, detached resource files.

    Reads only an already-resolved archive, never its original live paths. Every
    output file is a new regular inode owned by this user. Source resource links
    in native_directory are untouched: the process owner must explicitly select
    this tree before Popen. Existing materializations are validated, never repaired.
    """
    root,native=_paths(root,native_directory)
    _resolve_native_archive(root,native)
    archive=native/'environment-inputs';manifest_path=archive/'manifest.json'
    manifest_hash=_sha(manifest_path);receipt=json.loads(manifest_path.read_text())
    entries=_resource_entries(receipt);final=native/'runtime-resources'
    if final.exists() or final.is_symlink():return _validate_materialized_tree(final,entries,manifest_hash)
    staging=Path(tempfile.mkdtemp(prefix='.resource-materialization-',dir=native))
    try:
        # Parents remain writable during copying; recorded permission bits are
        # applied afterward. Setuid/setgid/sticky bits are never recreated.
        for name,entry in sorted(entries.items(),key=lambda item:(item[0].count('/'),item[0])):
            target=staging/name
            if entry['kind']=='directory':target.mkdir()
            else:
                source=archive/entry['blob']
                with target.open('xb') as output,source.open('rb') as source_stream:
                    shutil.copyfileobj(source_stream,output)
                if target.stat().st_nlink!=1 or _sha(target)!=entry['sha256']:
                    raise ValueError('Archived resource changed during materialization: '+name)
        record={'schema':'ihm.native-runtime-resources.v1','archive_manifest_sha256':manifest_hash,
                'resource_file_sha256':{name:entry['sha256'] for name,entry in entries.items() if entry['kind']=='file'},
                'resource_roots':list(RESOURCES)}
        (staging/'materialization.json').write_text(json.dumps(record,indent=2)+'\n')
        for name,entry in sorted(entries.items(),key=lambda item:item[0].count('/'),reverse=True):
            (staging/name).chmod(entry['mode']&0o777)
        if _sha(manifest_path)!=manifest_hash:raise ValueError('Archive receipt changed during materialization')
        _resolve_native_archive(root,native)
        _validate_materialized_tree(staging,entries,manifest_hash)
        if final.exists() or final.is_symlink():raise ValueError('Another materialization already owns this resource tree')
        staging.rename(final)
        return _validate_materialized_tree(final,entries,manifest_hash)
    except BaseException as error:
        if staging.exists():
            (staging/'failure.json').write_text(json.dumps({'error_type':type(error).__name__,'error':str(error)},indent=2)+'\n')
        raise
