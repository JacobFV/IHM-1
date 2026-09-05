"""Retain exact native regression dependencies and original fixture bytes.

This is explicitly a post-execution archive. Library/fixture identities are
checked against the actual regression receipts before copies are retained.
"""
from pathlib import Path
import argparse,json,sys
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from scripts.audit_long_horizon_thermal import sha
from ihm.assembly.systemic_evidence import freeze_sources
from ihm.assembly.native_environment_evidence import freeze_native_environment

def archive(directory):
    folder=Path(directory).resolve()
    if not folder.is_relative_to(BASE/'data/derived/audits'):raise ValueError('In-workspace regression audit required')
    manifest=folder/'verification.json';master=json.loads(manifest.read_text())
    if master['status']!='passed':raise ValueError('Completed regression receipt required')
    archives={}
    for name,case in master['reports'].items():
        report_path=BASE/case['report_path'];report=json.loads(report_path.read_text())
        if sha(report_path)!=case['report_sha256']:raise ValueError('Regression receipt changed')
        receipt=report if 'dependency_sha256' in report else report['receipts'][master['variant']]
        compile_command=receipt.get('compile_command',report.get('compile_command',report.get('command')))
        exe=Path(compile_command[compile_command.index('-o')+1])
        expected_binary=receipt.get('binary_sha256',receipt.get('probe_sha256'))
        if sha(exe)!=expected_binary:raise ValueError('Regression executable changed')
        if 'state_path' in report:state=Path(report['state_path'])
        elif name=='depletion_integrity':state=BASE/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml'
        else:state=BASE/'data/derived/audits/thermal-evaporation-humidity-fresh-hour/states/native_stabilized.xml'
        if sha(state)!=receipt['state_sha256']:raise ValueError('Actual regression fixture changed')
        destination=folder/name/'native-inputs'
        if destination.exists():raise ValueError('Preserve existing regression input archive')
        destination.mkdir()
        for resource in ('patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml'):
            target=BASE/'data/runtime/physiology/biogears-build/runtime'/resource
            (destination/resource).symlink_to(target,target_is_directory=target.is_dir())
        sources={str(p.relative_to(BASE)):sha(p) for p in (manifest,report_path,state,exe,Path(__file__))}
        sources.update(case['retained_inputs'])
        freeze_sources(BASE,destination,sources)
        # command[0] identifies the historical executable for the environment
        # archiver. It is not claimed to be a newly executed invocation.
        descriptor=dict(schema='ihm.native-session.v1',adapter_kind='post_execution_regression_inputs',command=[str(exe)],
            executable_sha256=sha(exe),dependency_sha256=receipt['dependency_sha256'],library_sha256=master['library_sha256'],
            interpretation='Historical regression executable identity and actual fixture; resources captured after all test commands. No new native execution or startup-resource reconstruction is claimed.')
        (destination/'manifest.json').write_text(json.dumps(descriptor,indent=2)+'\n')
        (destination/'receipts.jsonl').write_text(json.dumps({'command':'ARCHIVE_AFTER_COMPLETED_REGRESSION','regression_receipt':str(report_path.relative_to(BASE))})+'\n')
        archives[name]=freeze_native_environment(BASE,destination)
    result=dict(status='passed',archives=archives,master_receipt_sha256=sha(manifest),source_sha256=sha(__file__))
    (folder/'native-input-archives.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':list(archives),'receipt':str((folder/'native-input-archives.json').relative_to(BASE))},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);archive(p.parse_args().directory)
