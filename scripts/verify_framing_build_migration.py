"""Source-only exact-frame migration audit on retained native builds."""
import copy,json
from pathlib import Path
from framing_build_migration import audit
ROOT=Path(__file__).resolve().parents[1]

def main():
    previous=json.loads((ROOT/'data/derived/constrained-supine-vmmvo7v_/cache_identity.json').read_text())
    execution=json.loads((ROOT/'data/derived/pending-static-pose-rssrue1z/native/execution.json').read_text())
    current={**previous,'build_files':execution['build']['files']}
    result=audit(previous,current,ROOT);assert result['kind']=='exact-static-dispatch-framing-order-only'
    altered=copy.deepcopy(current);key=next(k for k in altered['build_files'] if k.endswith('native_static_pose.h'));altered['build_files'][key]='bad'
    try:audit(previous,altered,ROOT)
    except ValueError:pass
    else:raise AssertionError('Changed physical header accepted')
    print(json.dumps({'passed':True,'native_run':False,'audit':result,'changed_physics_rejected':True}))
if __name__=='__main__':main()
