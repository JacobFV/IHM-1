#!/usr/bin/env python3
"""Acquire official sources and index traceable model defaults / validation targets."""
from pathlib import Path
import argparse
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.forge.physiology_assets import collect, verify
parser=argparse.ArgumentParser()
parser.add_argument('--download',action='store_true',help='Clone missing official repositories')
parser.add_argument('--verify',action='store_true',help='Verify existing asset hashes and extraction provenance')
args=parser.parse_args()
base=Path(__file__).resolve().parents[1]
print(json.dumps(verify(base) if args.verify else collect(base,args.download),indent=2))
