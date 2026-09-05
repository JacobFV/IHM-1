"""Build reproducible native and optional measured temporal artifacts."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.temporal.atlas import build_atlas

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--native-csv',type=Path); p.add_argument('--output',type=Path)
    args=p.parse_args()
    atlas=build_atlas(Path(__file__).resolve().parents[1],args.native_csv,args.output)
    for run in atlas['runs']:
        d=run['predictor']['model']['diagnostics']
        print(f"{run['id']}: {run['samples']} samples, {len(run['variables'])} variables, {run['duration_s']:.3f}s, rank {d['data_rank']} reduced to {d['retained_rank']}, {d['holdout_samples']} chronological holdout samples")
