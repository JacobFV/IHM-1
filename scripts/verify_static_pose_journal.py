"""Source-only regression checks for bounded static-pose cache recovery."""
import json,tempfile
from pathlib import Path
import static_pose_journal as journal

def main():
    assert hasattr(journal,'PoseJournal'), 'Missing persistent candidate journal'
    with tempfile.TemporaryDirectory() as temporary:
        root=Path(temporary);identity={'protocol':'v1','native':'abc','gauges':{'y':1.}}
        first=journal.PoseJournal(root/'first',identity)
        entry={'coordinates':{'a':.25},'native':{'udot':[3.],'constrained_zero_acceleration_residual_mobility_force':[-6.]},'cost':2.}
        first.append([.25],entry)
        second=journal.PoseJournal(root/'second',identity,root/'first')
        assert second.cache[journal.pose_key([.25])]==entry
        assert journal.pose_key([.25000000000001]) not in second.cache
        try:journal.PoseJournal(root/'bad',{'protocol':'changed'},root/'first')
        except ValueError:pass
        else:raise AssertionError('Source-incompatible cache accepted')
        with (root/'first'/'candidates.jsonl').open('a') as stream:stream.write('{truncated')
        try:journal.PoseJournal(root/'broken',identity,root/'first')
        except ValueError:pass
        else:raise AssertionError('Incomplete cache silently accepted')
    print(json.dumps({'passed':True,'native_run':False,'checks':['full response recovery','exact pose keys','source mismatch rejection','truncated record rejection']}))
if __name__=='__main__':main()
