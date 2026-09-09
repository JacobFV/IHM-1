"""Freeze the shared geometric skin embedding from retained registration."""
import argparse,gzip,json
from pathlib import Path
from types import SimpleNamespace
from ihm.assembly.continuous_surface_binding import materialize,DEFAULT_ASSET

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser();parser.add_argument('--registration',type=Path,required=True);parser.add_argument('--output',type=Path,default=root/DEFAULT_ASSET);args=parser.parse_args()
    manifest=json.loads(args.registration.read_text());payload=json.loads((args.registration.parent/'canonical_mechanics.json').read_text())
    registration=SimpleNamespace(groups=manifest['groups'],specs={e['id']:e for e in payload['entities']},manifest=lambda:manifest)
    result=materialize(root,registration);raw=json.dumps(result,separators=(',',':'),allow_nan=False).encode();args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_bytes(gzip.compress(raw,mtime=0))
    print(json.dumps(dict(path=str(args.output),binding_identity=result['binding_identity'],algorithm=result['algorithm'],bytes=len(raw),compressed_bytes=args.output.stat().st_size),indent=2))
