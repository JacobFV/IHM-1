#!/usr/bin/env python3
"""Run validated native BioGears scenarios, preserving original trajectory and states."""
from pathlib import Path
import argparse,json,sys
base=Path(__file__).resolve().parents[1];sys.path.insert(0,str(base))
from ihm.native import NativeConfig,run_native,summarize
parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=60);parser.add_argument('--patient',default='StandardMale');parser.add_argument('--state',type=Path);parser.add_argument('--config',type=Path);parser.add_argument('--sample-hz',type=int,default=50);parser.add_argument('--summarize-only',action='store_true');parser.add_argument('--output-dir',type=Path,default=base/'data/derived/physiology/biogears_native_run');args=parser.parse_args()
config=NativeConfig.from_dict(json.loads(args.config.read_text())) if args.config else NativeConfig(seconds=args.seconds,patient=args.patient,state_path=args.state,sample_hz=args.sample_hz)
result=summarize(args.output_dir) if args.summarize_only else run_native(config,args.output_dir)
print(json.dumps({k:v for k,v in result.items() if k!='summary'},indent=2,default=str))
