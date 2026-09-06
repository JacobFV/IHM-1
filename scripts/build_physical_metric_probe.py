"""Isolated copied-engine mass-operator diagnostic; no shared source changes."""
from pathlib import Path
import argparse,json,shutil
from build_effective_potential_probe import prepare,build,verify,ROOT,sha

def prepare_metric():
    path=prepare();m=json.loads(path.read_text());folder=path.parent
    header=ROOT/'scripts/native_physical_metric_probe.h';shutil.copy2(header,folder/header.name)
    source=folder/'native_probe.cpp';text=source.read_text()
    anchor='#include "native_effective_potential_probe.h"'
    if text.count(anchor)!=1:raise ValueError('Missing copied diagnostic include boundary')
    source.write_text(text.replace(anchor,'#include "native_physical_metric_probe.h"'))
    for p in (header,folder/header.name,source,Path(__file__),ROOT/'scripts/verify_native_physical_metric.py'):
        m['files'][str(p.relative_to(ROOT))]=sha(p)
    m['schema']='ihm.physical-metric-probe-build.v1'
    m['source_changes']='Copied effective probe plus native calcM/calcMInv observations only; all physics unchanged'
    path.write_text(json.dumps(m,indent=2)+'\n');verify(path);return path
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--build',type=Path);a=p.parse_args()
    if a.prepare:print(prepare_metric())
    elif a.build:build(a.build.resolve())
    else:raise SystemExit('Choose source preparation or coordinated build')
