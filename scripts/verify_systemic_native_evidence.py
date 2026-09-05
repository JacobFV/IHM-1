"""Check native systemic binding against real held records and isolated corruption."""
from copy import deepcopy
import json
from pathlib import Path
import shutil
import tempfile
from ihm.assembly.systemic_native_evidence import resolve_systemic_execution

def main():
    root=Path(__file__).resolve().parents[1];checks=[]
    for name in ('rest','apnea'):
        source=root/'data/derived/systemic/respiratory_v3'/name
        data=json.loads((source/'systemic.json').read_text())
        paths=resolve_systemic_execution(root,source,data)
        assert len(paths)==5
        checks.append(name+'_real_native_actions_frames_and_checkpoints')
    source=root/'data/derived/systemic/respiratory_v3/apnea'
    with tempfile.TemporaryDirectory(prefix='systemic-native-evidence-',dir=root/'data/derived/audits') as temporary:
        target=Path(temporary)
        for relative in paths:
            destination=target/(root/relative).relative_to(source)
            destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/relative,destination)
        assert resolve_systemic_execution(root,target,data)
        def rejected(candidate,label):
            try:resolve_systemic_execution(root,target,candidate)
            except ValueError:checks.append(label)
            else:raise AssertionError(label)
        bad=deepcopy(data);bad['native_manifest']['library_sha256']='a'*64
        rejected(bad,'embedded_manifest_forgery_rejected')
        bad=deepcopy(data);bad['actions'][0]['kind']='exercise'
        rejected(bad,'unexecuted_action_relabel_rejected')
        bad=deepcopy(data);bad['frames'][-1]['values']['stomach_water_ml']=-1000
        rejected(bad,'changed_sample_rejected')
        # Even rewriting the derived sampled trace cannot change the native
        # observation that the acknowledged process actually emitted.
        frames=target/'frames.jsonl';original=frames.read_bytes()
        frames.write_text(''.join(json.dumps(frame)+'\n' for frame in bad['frames']))
        rejected(bad,'changed_sample_and_derived_trace_rejected');frames.write_bytes(original)
        journal=target/'native/receipts.jsonl';original=journal.read_bytes()
        records=[json.loads(line) for line in original.splitlines()]
        records[-1]['acknowledgment']['elapsed_s']=float('nan')
        journal.write_text(''.join(json.dumps(record)+'\n' for record in records))
        rejected(data,'nonfinite_native_clock_rejected');journal.write_bytes(original)
        state=next((target/'native/states').iterdir());state.write_bytes(b'changed state')
        rejected(data,'changed_checkpoint_rejected')
    report={'passed':True,'checks':checks}
    out=root/'artifacts/verification/systemic-native-evidence';out.mkdir(parents=True,exist_ok=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
