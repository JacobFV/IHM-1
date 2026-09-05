"""Native dependencies and runtime resources survive mutation outside the archive."""
from pathlib import Path
import sys,tempfile,json,hashlib,shutil,argparse,os
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

RESOURCES=('patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def verify():
    try:from ihm.assembly.native_environment_evidence import freeze_native_environment,resolve_native_environment
    except ImportError:raise AssertionError('Immutable native environment API missing') from None
    with tempfile.TemporaryDirectory() as tmp,tempfile.TemporaryDirectory() as external:
        root=Path(tmp);native=root/'experiment/native';native.mkdir(parents=True)
        runtime=root/'runtime';runtime.mkdir();executable=root/'engine';shutil.copyfile('/bin/true',executable)
        library=Path(external)/'library.so';library.write_bytes(b'external linked bytes')
        for name in RESOURCES:
            target=runtime/name
            if '.' in name:target.write_text('runtime input '+name)
            else:
                target.mkdir();(target/'input.xml').write_text('runtime input '+name)
                (target/'empty').mkdir()
            (native/name).symlink_to(target,target_is_directory=target.is_dir())
        manifest={'schema':'ihm.native-session.v1','command':[str(executable)],'executable_sha256':sha(executable),'dependency_sha256':{str(library):sha(library)}}
        (native/'manifest.json').write_text(json.dumps(manifest));(native/'receipts.jsonl').write_text('')
        (native/'states').mkdir();(native/'states/private-output.xml').write_text('MUST NOT ARCHIVE')
        receipt=freeze_native_environment(root,native)
        archive=json.loads((native/'environment-inputs/manifest.json').read_text())
        assert archive['capture']['stage']=='post_start_before_first_action'
        assert archive['capture']['resources_verified_at_process_start'] is False
        assert len(archive['dependencies'])==1 and archive['elf_interpreter'] is not None
        assert all(str(root/p).startswith(str(native)) for p in receipt)
        assert not any('states/' in item['logical_path'] for item in archive['resources'])
        assert any(item['kind']=='directory' and item['logical_path'].endswith('/empty') for item in archive['resources'])
        # Detached bytes, not symlinks/hardlinks to mutable source inputs.
        library.write_bytes(b'updated external library');(runtime/'substances/input.xml').write_text('updated substance')
        executable.unlink()
        assert resolve_native_environment(root,native)==receipt
        assert freeze_native_environment(root,native)==receipt
        receipt_path=native/'environment-inputs/manifest.json';original_receipt=receipt_path.read_bytes()
        incomplete=json.loads(original_receipt);dependency=incomplete['dependencies'].pop()
        del incomplete['blobs'][dependency['blob']]
        missing_blob=native/'environment-inputs'/dependency['blob'];missing_bytes=missing_blob.read_bytes();missing_blob.unlink()
        # Even with a rewritten self-consistent receipt, startup declares this dependency.
        receipt_path.write_text(json.dumps(incomplete))
        try:resolve_native_environment(root,native)
        except ValueError:pass
        else:raise AssertionError('Startup dependency omitted from archive acceptance')
        receipt_path.write_bytes(original_receipt)
        missing_blob.write_bytes(missing_bytes)
        blob=native/'environment-inputs'/archive['dependencies'][0]['blob']
        original=blob.read_bytes();blob.write_bytes(b'corruption')
        try:resolve_native_environment(root,native)
        except ValueError:pass
        else:raise AssertionError('Corrupted native environment accepted')
        blob.write_bytes(original)
        # Matching bytes at an external path must never weaken archive containment.
        outside=Path(external)/'same-bytes';outside.write_bytes(original)
        blob.unlink();blob.symlink_to(outside)
        try:resolve_native_environment(root,native)
        except ValueError:pass
        else:raise AssertionError('External archive symlink accepted')
        blob.unlink();os.link(outside,blob)
        try:resolve_native_environment(root,native)
        except ValueError:pass
        else:raise AssertionError('Archive hardlink shares mutable external ownership')
        failed=root/'failed/native';failed.mkdir(parents=True)
        (failed/'manifest.json').write_text(json.dumps(manifest))
        try:freeze_native_environment(root,failed)
        except ValueError:pass
        else:raise AssertionError('Changed startup library was captured as original')
        assert list(failed.glob('.environment-capture-*/failure.json'))
        assert not (failed/'environment-inputs').exists()
        # An old executable can be recovered from the run's already-validated
        # systemic snapshot; a newer live executable must not be silently substituted.
        from ihm.assembly.systemic_evidence import freeze_sources
        historical=root/'historical/native';historical.mkdir(parents=True)
        shutil.copyfile(native/'environment-inputs'/archive['executable']['blob'],executable)
        frozen=freeze_sources(root,historical.parent,{'engine':sha(executable)})
        old_manifest={**manifest,'dependency_sha256':{str(library):sha(library)}}
        (historical/'manifest.json').write_text(json.dumps(old_manifest))
        (historical/'receipts.jsonl').write_text(json.dumps({'command':'STEP 1','before_ticks':0})+'\n')
        for name in RESOURCES:(historical/name).symlink_to(runtime/name,target_is_directory=(runtime/name).is_dir())
        executable.write_bytes(b'new executable version')
        try:freeze_native_environment(root,historical)
        except ValueError as error:raise AssertionError('Existing immutable executable could not support historical environment capture') from error
        captured=json.loads((historical/'environment-inputs/manifest.json').read_text())
        assert captured['capture']['stage']=='post_start_after_commands'
        assert captured['executable']['sha256']==old_manifest['executable_sha256']
    print('PASS detached native environment identity, timing, ownership and tamper checks')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--capture',type=Path);args=p.parse_args();verify()
    if args.capture:
        from ihm.assembly.native_environment_evidence import freeze_native_environment
        result=freeze_native_environment(Path(__file__).resolve().parents[1],args.capture)
        print(json.dumps({'capture':str(args.capture),'verified_files':len(result)},indent=2))
