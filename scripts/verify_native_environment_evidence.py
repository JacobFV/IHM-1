"""Native dependencies and runtime resources survive mutation outside the archive."""
from pathlib import Path
import sys,tempfile,json,hashlib,shutil,argparse,os
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

RESOURCES=('patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def verify():
    try:from ihm.assembly.native_environment_evidence import freeze_native_environment,resolve_native_environment,materialize_native_resources
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
    verify_prestart()

def verify_prestart():
    from ihm.assembly.native_environment_evidence import freeze_native_environment,resolve_native_environment,materialize_native_resources
    with tempfile.TemporaryDirectory() as temporary:
        root=Path(temporary);native=root/'new/native';native.mkdir(parents=True)
        runtime=root/'live';runtime.mkdir();engine=root/'engine';shutil.copyfile('/bin/true',engine)
        library=root/'library.so';library.write_bytes(b'library')
        for name in RESOURCES:
            source=runtime/name
            if '.' in name:source.write_text('known '+name)
            else:source.mkdir();(source/'input.xml').write_text('known '+name);(source/'empty').mkdir()
            (native/name).symlink_to(source,target_is_directory=source.is_dir())
        state=native/'input-state.xml';state.write_bytes(b'independent initial patient state')
        manifest={'schema':'ihm.native-session.v1','command':[str(engine),'input',str(state),'100'],
                  'configuration':{'patient':'input','state_path':str(root/'original-state.xml')},
                  'initialization_patient_used':False,'state_sha256':sha(state),
                  'patient_identity_input_sha256':sha(state),
                  'executable_sha256':sha(engine),'dependency_sha256':{str(library):sha(library)}}
        (native/'manifest.json').write_text(json.dumps(manifest));(native/'receipts.jsonl').write_text('')
        identity=freeze_native_environment(root,native,before_start=True)
        archive=native/'environment-inputs';receipt=json.loads((archive/'manifest.json').read_text())
        assert receipt['capture']['stage']=='pre_start'
        assert receipt['capture']['resources_verified_at_process_start'] is False,'Capture alone cannot attest process consumption'
        (runtime/'substances/input.xml').write_text('changed before materialization')
        tree=materialize_native_resources(root,native)
        assert isinstance(tree,Path) and tree.is_absolute() and tree.is_relative_to(native)
        assert (tree/'substances/input.xml').read_text()=='known substances'
        assert (tree/'substances/empty').is_dir()
        assert all(not p.is_symlink() and (not p.is_file() or p.stat().st_nlink==1) for p in tree.rglob('*'))
        assert materialize_native_resources(root,native)==tree
        try:resolve_native_environment(root,native)
        except ValueError:pass
        else:raise AssertionError('Prestart capture accepted without execution selection receipt')
        receipt_path=archive/'manifest.json';retained_receipt=receipt_path.read_bytes()
        for unsafe in ('../escape','/tmp/escape','patients/../escape','patients//input.xml','.'):
            changed=json.loads(retained_receipt)
            next(r for r in changed['resources'] if r['logical_path']=='substances/input.xml')['logical_path']=unsafe
            receipt_path.write_text(json.dumps(changed))
            try:materialize_native_resources(root,native)
            except ValueError:pass
            else:raise AssertionError('Unsafe archived resource path accepted: '+unsafe)
            receipt_path.write_bytes(retained_receipt)
        for mutation in ('duplicate','missing_root','file_parent'):
            changed=json.loads(retained_receipt)
            if mutation=='duplicate':changed['resources'].append(changed['resources'][0].copy())
            elif mutation=='missing_root':changed['resources']=[r for r in changed['resources'] if r['logical_path']!='patients']
            else:next(r for r in changed['resources'] if r['logical_path']=='substances/input.xml')['logical_path']='UCEDefs.conf/child'
            receipt_path.write_text(json.dumps(changed))
            try:materialize_native_resources(root,native)
            except ValueError:pass
            else:raise AssertionError('Malformed archived resource tree accepted: '+mutation)
            receipt_path.write_bytes(retained_receipt)
        # Actual consumer opens the detached path after the live resource changed.
        import subprocess
        consumed=subprocess.check_output([sys.executable,'-c','from pathlib import Path;print(Path("substances/input.xml").read_text())'],cwd=tree,text=True).strip()
        assert consumed=='known substances'
        original=(tree/'substances/input.xml').read_bytes()
        for tamper in ('bytes','symlink','hardlink'):
            target=tree/'substances/input.xml';target.unlink()
            if tamper=='bytes':target.write_bytes(b'changed owned tree')
            elif tamper=='symlink':target.symlink_to(runtime/'substances/input.xml')
            else:os.link(archive/next(r['blob'] for r in receipt['resources'] if r['logical_path']=='substances/input.xml'),target)
            try:materialize_native_resources(root,native)
            except ValueError:pass
            else:raise AssertionError('Materialized resource '+tamper+' corruption accepted')
            target.unlink();target.write_bytes(original)
        verify_execution_selection(root,native,tree)
        # A second call with different phase intent must never rewrite an archive.
        archive_bytes=(archive/'manifest.json').read_bytes()
        assert freeze_native_environment(root,native)==identity
        assert (archive/'manifest.json').read_bytes()==archive_bytes
        # Legacy archives retain their original observed-after-start phase even
        # when a later caller requests a prestart interpretation.
        previous=json.loads(archive_bytes);previous['capture']['stage']='post_start_before_first_action'
        previous['capture']['interpretation']='Legacy post-start capture; no retrospective startup claim.'
        receipt_path.write_text(json.dumps(previous));legacy_bytes=receipt_path.read_bytes()
        freeze_native_environment(root,native,before_start=True)
        assert receipt_path.read_bytes()==legacy_bytes
        receipt_path.write_bytes(archive_bytes)
        # Patient initialization has no saved state; it must bind the selected XML
        # inside the exact detached inventory and the command's '-' mode.
        initialized=root/'initialized/native';initialized.mkdir(parents=True)
        init_manifest={**manifest,'command':[str(engine),'input','-','100'],
                       'configuration':{'patient':'input','state_path':None},
                       'initialization_patient_used':True,'state_sha256':None,
                       'patient_identity_input_sha256':sha(tree/'patients/input.xml')}
        (initialized/'manifest.json').write_text(json.dumps(init_manifest))
        for name in RESOURCES:(initialized/name).symlink_to(tree/name,target_is_directory=(tree/name).is_dir())
        freeze_native_environment(root,initialized,before_start=True)
        initialized_tree=materialize_native_resources(root,initialized)
        verify_execution_selection(root,initialized,initialized_tree)
        late=root/'late/native';late.mkdir(parents=True);(late/'manifest.json').write_text(json.dumps(manifest))
        (late/'receipts.jsonl').write_text(json.dumps({'command':'STEP 1'})+'\n')
        try:freeze_native_environment(root,late,before_start=True)
        except ValueError:pass
        else:raise AssertionError('Post-action directory labeled as prestart')
        assert not (late/'environment-inputs').exists()
    print('PASS prestart archive, detached resource consumption and materialization tamper checks')

def verify_execution_selection(root,native,tree):
    from ihm.assembly.native_environment_evidence import resolve_native_environment,resolve_native_execution_inputs
    def rejected(label):
        try:resolve_native_environment(root,native)
        except (ValueError,FileNotFoundError):pass
        else:raise AssertionError('Accepted invalid execution selection: '+label)
    for name in RESOURCES:
        link=native/name;link.unlink();link.symlink_to(tree/name,target_is_directory=(tree/name).is_dir())
    manifest=json.loads((native/'manifest.json').read_text())
    state=tree/'patients'/f"{manifest['configuration']['patient']}.xml" if manifest['initialization_patient_used'] else native/'input-state.xml'
    launch={'native_manifest_sha256':sha(native/'manifest.json'),
            'environment_manifest_sha256':sha(native/'environment-inputs/manifest.json'),
            'resource_materialization_sha256':sha(tree/'materialization.json'),
            'selected_resource_tree':str(tree),'selected_patient_input':str(state),
            'selected_patient_input_sha256':sha(state)}
    path=native/'execution-inputs.json';path.write_text(json.dumps(launch));original=path.read_bytes()
    identity=resolve_native_environment(root,native)
    selected=resolve_native_execution_inputs(root,native)
    assert selected.items()<=identity.items()
    assert str(state.relative_to(root)) in selected and str(path.relative_to(root)) in selected
    assert str((tree/'substances/input.xml').relative_to(root)) in selected
    for key in ('native_manifest_sha256','environment_manifest_sha256','resource_materialization_sha256','selected_patient_input_sha256'):
        changed={**launch,key:'a'*64};path.write_text(json.dumps(changed));rejected(key);path.write_bytes(original)
    for key,value in (('selected_resource_tree',str(root/'outside/runtime-resources')),
                      ('selected_resource_tree',str(native/'other/../runtime-resources')),
                      ('selected_patient_input',str(root/'original-state.xml')),
                      ('selected_patient_input',str(native/'../input-state.xml'))):
        path.write_text(json.dumps({**launch,key:value}));rejected(key+value);path.write_bytes(original)
    # Files with correct bytes but shared or indirect ownership are not retained inputs.
    for target in (path,state,tree/'materialization.json',tree/'substances/input.xml'):
        saved=target.read_bytes();outside=root/'same-bytes';outside.write_bytes(saved)
        for mutation in ('missing','changed','symlink','hardlink'):
            target.unlink()
            if mutation=='changed':target.write_bytes(b'changed')
            elif mutation=='symlink':target.symlink_to(outside)
            elif mutation=='hardlink':os.link(outside,target)
            rejected(str(target)+mutation)
            if target.exists() or target.is_symlink():target.unlink()
            target.write_bytes(saved)
    link=native/'substances';link.unlink();link.symlink_to(root/'live/substances')
    rejected('live resource root');link.unlink();link.symlink_to(tree/'substances')
    tree.rename(native/'temporarily-absent-resources');rejected('missing selected tree')
    (native/'temporarily-absent-resources').rename(tree)
    assert resolve_native_environment(root,native)==identity
    # Relocate every selected file and preserve old absolute link text. Resolution
    # maps native-local identities and must not read any of those old targets.
    with tempfile.TemporaryDirectory() as relocated:
        moved_root=Path(relocated);moved=moved_root/native.relative_to(root)
        shutil.copytree(native,moved,symlinks=True)
        native.rename(native.with_name('hidden-original'))
        try:assert resolve_native_environment(moved_root,moved)==identity
        finally:native.with_name('hidden-original').rename(native)
    print('PASS execution selection receipt, copied state, detached resource binding and relocation')


def verify_held(native):
    """Consume a relocated copy of an actual retained archive; never rewrite it."""
    from ihm.assembly.native_environment_evidence import resolve_native_environment,materialize_native_resources
    import subprocess
    root=Path(__file__).resolve().parents[1];native=Path(native).resolve()
    original_identity=resolve_native_environment(root,native)
    receipt=json.loads((native/'environment-inputs/manifest.json').read_text())
    with tempfile.TemporaryDirectory(prefix='ihm-held-resource-consumer-') as temporary:
        detached_root=Path(temporary);copy=detached_root/receipt['native_directory'];copy.mkdir(parents=True)
        shutil.copyfile(native/'manifest.json',copy/'manifest.json')
        shutil.copytree(native/'environment-inputs',copy/'environment-inputs',symlinks=True)
        tree=materialize_native_resources(detached_root,copy)
        if receipt['capture']['stage']=='pre_start':
            shutil.copyfile(native/'execution-inputs.json',copy/'execution-inputs.json')
            if (native/'input-state.xml').exists():shutil.copyfile(native/'input-state.xml',copy/'input-state.xml')
            for name in RESOURCES:(copy/name).symlink_to(os.readlink(native/name))
            assert resolve_native_environment(detached_root,copy)==original_identity
        files=[r for r in receipt['resources'] if r['kind']=='file']
        for entry in files:
            resource=tree/entry['logical_path'];blob=copy/'environment-inputs'/entry['blob']
            assert sha(resource)==entry['sha256'] and resource.stat().st_ino!=blob.stat().st_ino
        observed=subprocess.check_output([sys.executable,'-c','from pathlib import Path;import hashlib;print(hashlib.sha256(Path("UCEDefs.conf").read_bytes()).hexdigest())'],cwd=tree,text=True).strip()
        expected=next(r['sha256'] for r in files if r['logical_path']=='UCEDefs.conf')
        assert observed==expected
    assert resolve_native_environment(root,native)==original_identity
    result={'native_directory':receipt['native_directory'],'capture_stage':receipt['capture']['stage'],
            'resource_files':len(files),'regular_detached_files_consumable':True,
            'original_archive_unchanged':True,'native_manifest_sha256':receipt['native_manifest_sha256'],
            'archive_manifest_sha256':sha(native/'environment-inputs/manifest.json'),'consumer_ucedefs_sha256':observed}
    print(json.dumps(result,indent=2));return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--capture',type=Path);p.add_argument('--held',type=Path,action='append',default=[]);args=p.parse_args();verify()
    if args.held:
        results=[verify_held(path) for path in args.held]
        out=Path(__file__).resolve().parents[1]/'artifacts/verification/native-environment-prestart';out.mkdir(parents=True,exist_ok=True)
        (out/'held-resources.json').write_text(json.dumps(results,indent=2)+'\n')
    if args.capture:
        from ihm.assembly.native_environment_evidence import freeze_native_environment
        result=freeze_native_environment(Path(__file__).resolve().parents[1],args.capture)
        print(json.dumps({'capture':str(args.capture),'verified_files':len(result)},indent=2))
