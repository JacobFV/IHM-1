"""Build the single-body contract and optionally materialize a native run."""
from pathlib import Path
import argparse,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.body import build,CanonicalBody
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--native-directory',type=Path);p.add_argument('--output',type=Path,default=ROOT/'data/derived/canonical/trajectory.json');p.add_argument('--output-hz',type=float,default=10);a=p.parse_args()
 payload=build(ROOT);print('Canonical body manifest:',payload['entity_count'],'entities')
 if a.native_directory:
  from ihm.assembly.temporal import build_body_spectra
  result=CanonicalBody.from_workspace(ROOT).simulate(a.native_directory,a.output,a.output_hz);print('Materialized',len(result['frames']),'synchronized body frames:',a.output)
  build_body_spectra(a.native_directory,a.output.with_name(a.output.stem+'-spectra.json'),payload)
