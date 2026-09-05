"""Build a source-pinned forearm contact → IBM causal receptor experiment."""
from pathlib import Path
import sys,json,os,argparse
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--threads',type=int,default=1);args=parser.parse_args()
    if not 1<=args.threads<=32:parser.error('threads must be in1..32')
    os.environ['OPENBLAS_NUM_THREADS']=str(args.threads)
    from ihm.assembly.regional_touch import run_touch
    result=run_touch(ROOT)
    out=ROOT/'data/derived/canonical/forearm-touch.json'
    out.write_text(json.dumps(result,separators=(',',':'),allow_nan=False)+'\n')
    print(json.dumps({'output':str(out),'tetrahedra':len(result['geometry']['tetrahedra']),
        'reaction_n':result['states']['loaded']['indenter_reaction_n'],'minimum_jacobian':result['states']['loaded']['minimum_jacobian'],
        'force_residual_n':result['states']['loaded']['force_balance_residual_n'],'frames':len(result['frames'])}))
