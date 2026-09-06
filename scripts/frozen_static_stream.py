"""Pinned pre-mass static adapter with independently attested archived physics."""
from pathlib import Path
import argparse,hashlib,json,os,re,shutil,subprocess,tempfile,types

BUILD='data/runtime/mechanical-stream/build-mptn9l8d'
BUILD_MANIFEST_SHA='44b87ebe4f0defc5bc973b3e132dd60d631a8d02f02f55a18b893b2007aec1fa'
ADAPTER_REVISION='91c97b1'
ADAPTER_SHA='c53ec4b4c5ea6ca569b7e365e8ce044bfd8daffd185178c7977be9600a11c1cc'
GUARD='''        pointer=json.loads((self.root/'data/runtime/mechanical-stream/latest.json').read_text());build=self.root/pointer['build'];manifest=json.loads((build/'manifest.json').read_text())
        for path,digest in manifest['files'].items():
            if sha(self.root/path)!=digest:raise ValueError('Stale native mechanical build: '+path)'''
ENV="        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',LD_LIBRARY_PATH=':'.join(map(str,libdirs)))"
ENGINE="        engine_command=[str(build/'native_mechanical_stream'),str(source),str(self.output),environment,str(float(target_mass_kg))]"
RECEIPT="        (self.output/'execution.json').write_text(json.dumps(execution,indent=2)+'\\n')"


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(root):
    root=Path(root).resolve();origin=root/BUILD/'manifest.json'
    if sha(origin)!=BUILD_MANIFEST_SHA:raise ValueError('Pinned build manifest changed')
    build=json.loads(origin.read_text());adapter=subprocess.check_output(['git','show',ADAPTER_REVISION+':ihm/native/mechanical_stream.py'],cwd=root)
    if hashlib.sha256(adapter).hexdigest()!=ADAPTER_SHA:raise ValueError('Pinned adapter source changed')
    output=Path(tempfile.mkdtemp(prefix='frozen-static-stream-',dir=root/'data/derived'))
    (output/'engine').mkdir();(output/'libraries').mkdir();(output/'sources').mkdir()
    (output/'origin_build_manifest.json').write_bytes(origin.read_bytes());(output/'adapter_original.py').write_bytes(adapter)
    mapping={};copies={};aliases={}
    for original,expected in build['files'].items():
        source=root/original
        if original.startswith('scripts/'):
            source=root/BUILD/Path(original).name
            if build['files'].get(str(source.relative_to(root)))!=expected:raise ValueError('Missing matching compiled source snapshot')
        if sha(source)!=expected:raise ValueError('Pinned compiled dependency changed: '+original)
        if '.so' in source.name:
            target=output/'libraries'/source.name
            if target.exists() and sha(target)!=expected:raise ValueError('Conflicting library basename')
            shutil.copy2(source,target)
            dynamic=subprocess.check_output(['readelf','-d',str(target)],text=True)
            found=re.search(r'\(SONAME\).*\[(.*?)\]',dynamic)
            if found and found.group(1)!=target.name:
                alias=output/'libraries'/found.group(1)
                if alias.exists() and sha(alias)!=expected:raise ValueError('Conflicting library SONAME')
                if not alias.exists():alias.symlink_to(target.name)
                aliases[str(alias.relative_to(output))]=str(target.relative_to(output))
        elif source.name=='native_mechanical_stream':target=output/'engine'/source.name;shutil.copy2(source,target)
        else:target=output/'sources'/source.name;shutil.copy2(source,target)
        mapping[original]=str(target.relative_to(output));copies[str(target.relative_to(output))]=expected
    # Resolve the remaining OS runtime under the original search directories,
    # with copied declared physics libraries taking precedence. ldd traces loader
    # resolution only; it does not initialize a model or run the executable main.
    runtime=root/'data/runtime/opensim'
    search=[output/'libraries',runtime/'install/opensim/lib',runtime/'install/simbody/lib',
        runtime/'sysroot/usr/lib/aarch64-linux-gnu/lapack',runtime/'sysroot/usr/lib/aarch64-linux-gnu/blas',runtime/'sysroot/usr/lib/aarch64-linux-gnu']
    resolved=subprocess.check_output(['ldd',str(output/'engine/native_mechanical_stream')],env=dict(os.environ,LD_LIBRARY_PATH=':'.join(map(str,search))),text=True)
    (output/'loader_resolution.txt').write_text(resolved);runtime_files={}
    for path in re.findall(r'(?:=>\s*)?(/[^\s()]+)\s+\(',resolved):
        source=Path(path)
        if source.is_relative_to(output):continue
        target=output/'libraries'/source.name;expected=sha(source)
        if target.exists() and sha(target)!=expected:raise ValueError('Resolved runtime library collision')
        shutil.copy2(source,target);copies[str(target.relative_to(output))]=expected
        runtime_files[str(source)]=dict(copy=str(target.relative_to(output)),sha256=expected)
    loader=next((item['copy'] for path,item in runtime_files.items() if Path(path).name.startswith('ld-linux-')),None)
    if loader is None:raise ValueError('Native dynamic loader was not resolved')
    for binary in [output/'engine/native_mechanical_stream',*(output/'libraries').iterdir()]:
        dynamic=subprocess.check_output(['readelf','-d',str(binary)],text=True)
        if '(RPATH)' in dynamic:raise ValueError('Legacy RPATH could bypass frozen library search')
    original=adapter.decode()
    for pattern in (GUARD,ENV,ENGINE,RECEIPT):
        if original.count(pattern)!=1:raise ValueError('Pinned adapter patch boundary mismatch')
    patched=original.replace(GUARD,'        build,manifest=_frozen_attestation(self.root)')
    patched=patched.replace(ENV,ENV+"\n        env['LD_LIBRARY_PATH']=_frozen_library_path")
    patched=patched.replace(ENGINE,ENGINE+"\n        engine_command=[_frozen_loader,'--library-path',_frozen_library_path,*engine_command]").replace(RECEIPT,"        execution['frozen_archive_attestation']=copy.deepcopy(_frozen_receipt)\n"+RECEIPT)
    (output/'adapter_frozen.py').write_text(patched)
    copies.update({name:sha(output/name) for name in ('adapter_original.py','adapter_frozen.py','origin_build_manifest.json','loader_resolution.txt')})
    record=dict(schema='ihm.frozen-static-stream.v1',origin_build=BUILD,origin_build_manifest_sha256=BUILD_MANIFEST_SHA,
        adapter_revision=ADAPTER_REVISION,adapter_original_sha256=ADAPTER_SHA,files=copies,original_build_validation_map=mapping,
        library_aliases=aliases,runtime_files=runtime_files,dynamic_loader=loader,scope='Archived compiled source/header/binary/library validation; current worktree source is not the validation target',
        adapter_changes=['replace latest/worktree build check with full pinned archive attestation','point LD_LIBRARY_PATH and explicit copied loader to verified declared and OS runtime libraries','emit explicit archived-source attestation'],
        system_runtime_scope='Loader and resolved OS runtime libraries copied and verified in addition to all declared build libraries; kernel/vDSO and CPU remain host resources')
    (output/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    validate(output/'manifest.json');return output/'manifest.json'


def validate(path):
    path=Path(path).resolve();record=json.loads(path.read_text());folder=path.parent
    if record.get('schema')!='ihm.frozen-static-stream.v1' or record.get('origin_build')!=BUILD:raise ValueError('Unknown frozen stream identity')
    if record['origin_build_manifest_sha256']!=BUILD_MANIFEST_SHA or record['adapter_original_sha256']!=ADAPTER_SHA:raise ValueError('Frozen origin pin mismatch')
    for relative,expected in record['files'].items():
        target=(folder/relative).resolve()
        if not target.is_relative_to(folder) or sha(target)!=expected:raise ValueError('Frozen artifact hash mismatch: '+relative)
    origin=json.loads((folder/'origin_build_manifest.json').read_text())
    if sha(folder/'origin_build_manifest.json')!=BUILD_MANIFEST_SHA or sha(folder/'adapter_original.py')!=ADAPTER_SHA:raise ValueError('Frozen source attestation mismatch')
    if set(record['original_build_validation_map'])!=set(origin['files']):raise ValueError('Incomplete native dependency attestation')
    for original,expected in origin['files'].items():
        target=folder/record['original_build_validation_map'][original]
        if not target.resolve().is_relative_to(folder) or sha(target)!=expected:raise ValueError('Compiled dependency attestation mismatch: '+original)
    loader=folder/record['dynamic_loader']
    if not loader.resolve().is_relative_to(folder) or record['dynamic_loader'] not in record['files'] or not loader.name.startswith('ld-linux-'):
        raise ValueError('Unattested frozen dynamic loader')
    for original,item in record['runtime_files'].items():
        target=folder/item['copy']
        if not target.resolve().is_relative_to(folder) or record['files'].get(item['copy'])!=item['sha256'] or sha(target)!=item['sha256']:
            raise ValueError('Frozen OS runtime attestation mismatch')
    for alias,target in record['library_aliases'].items():
        if (folder/alias).resolve()!=(folder/target).resolve():raise ValueError('Frozen library alias mismatch')
    original=(folder/'adapter_original.py').read_text()
    expected=original.replace(GUARD,'        build,manifest=_frozen_attestation(self.root)').replace(ENV,ENV+"\n        env['LD_LIBRARY_PATH']=_frozen_library_path").replace(ENGINE,ENGINE+"\n        engine_command=[_frozen_loader,'--library-path',_frozen_library_path,*engine_command]").replace(RECEIPT,"        execution['frozen_archive_attestation']=copy.deepcopy(_frozen_receipt)\n"+RECEIPT)
    if (folder/'adapter_frozen.py').read_text()!=expected:raise ValueError('Unexpected frozen adapter transformation')
    return record,origin


def load(path):
    path=Path(path).resolve();record,origin=validate(path)
    module=types.ModuleType('ihm_frozen_static_stream')
    receipt=dict(manifest=str(path),manifest_sha256=sha(path),origin_build=BUILD,
        validation_scope=record['scope'],original_build_validation_map=record['original_build_validation_map'],
        adapter_original_sha256=ADAPTER_SHA,adapter_frozen_sha256=sha(path.parent/'adapter_frozen.py'),
        copied_libraries={name:expected for name,expected in record['files'].items() if name.startswith('libraries/')},
        dynamic_loader=record['dynamic_loader'],runtime_files=record['runtime_files'])
    def attestation(root):
        validate(path)
        return path.parent/'engine',origin
    module.__dict__.update(_frozen_attestation=attestation,_frozen_library_path=str(path.parent/'libraries'),_frozen_loader=str(path.parent/record['dynamic_loader']),_frozen_receipt=receipt)
    exec(compile((path.parent/'adapter_frozen.py').read_text(),str(path.parent/'adapter_frozen.py'),'exec'),module.__dict__)
    return module.NativeMechanicalStream


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--prepare',action='store_true');parser.add_argument('--verify',type=Path);args=parser.parse_args()
    if args.prepare:print(prepare(Path(__file__).resolve().parents[1]))
    elif args.verify:record,_=validate(args.verify);load(args.verify);print(json.dumps({'passed':True,'native_run':False,'files':len(record['files'])}))
    else:raise SystemExit('Choose source-only --prepare or --verify')
