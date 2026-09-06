"""Capture an opt-in IBM source candidate; never replaces the active artifact."""
import argparse
import json
from pathlib import Path
from ihm.brain.candidate import capture_candidate
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT.parent/'IBM-1')
    parser.add_argument('--revision',default='HEAD')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    pin=capture_candidate(args.source,args.output,revision=args.revision,
        baseline_neural=ROOT/'data/derived/canonical/brain-sources/ibm-neural.py')
    print(json.dumps(pin.to_dict(),indent=2))
