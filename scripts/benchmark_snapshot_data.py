"""Benchmark independent copies of two retained browser frames, not simulation speed."""
import argparse,json,time,statistics,hashlib,sys
from pathlib import Path
from copy import deepcopy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.snapshot_data import clone_snapshot_data

def main():
    parser=argparse.ArgumentParser();parser.add_argument('frames',type=Path);parser.add_argument('--output',type=Path);parser.add_argument('--repeats',type=int,default=5);args=parser.parse_args()
    if not 1<=args.repeats<=20:parser.error('repeats must be between 1 and 20')
    raw=args.frames.read_bytes();payload=json.loads(raw);frames=[payload['first'],payload['last']]
    reference=deepcopy(frames);result=clone_snapshot_data(frames)
    if result!=reference:raise AssertionError('Retained frame copy differs from deepcopy')
    result[0]['mechanics']['entities']['snapshot-copy-audit']={'mutation':True}
    if frames!=reference:raise AssertionError('Retained source frames were mutated')
    timings={}
    for name,clone in [('deepcopy',deepcopy),('clone_snapshot_data',clone_snapshot_data)]:
        samples=[]
        for _ in range(args.repeats):
            start=time.perf_counter();copied=clone(frames);samples.append(time.perf_counter()-start)
            if copied!=reference:raise AssertionError('Copy changed values')
        timings[name]={'samples_s':samples,'median_s':statistics.median(samples)}
    report={'scope':'Copy-only benchmark; no claim of native end-to-end speedup','input':str(args.frames),'input_sha256':hashlib.sha256(raw).hexdigest(),'frame_times_s':[f['time_s'] for f in frames],'exact_deepcopy_equality':True,'source_mutation_isolation':True,'timings':timings,'copy_speedup':timings['deepcopy']['median_s']/timings['clone_snapshot_data']['median_s']}
    if args.output:args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
