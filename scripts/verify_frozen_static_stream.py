"""Source-only frozen adapter identity and tamper-rejection fixtures."""
from pathlib import Path
import argparse,copy,json
from frozen_static_stream import validate,load


def main(path):
    path=Path(path);record,origin=validate(path);klass=load(path)
    assert callable(klass) and len(record['original_build_validation_map'])==len(origin['files'])
    temporary=path.parent/'tampered-fixture-manifest.json'
    try:
        for key in ['adapter_original.py','adapter_frozen.py',next(k for k in record['files'] if k.startswith('libraries/'))]:
            changed=copy.deepcopy(record);changed['files'][key]='0'*64;temporary.write_text(json.dumps(changed))
            try:validate(temporary)
            except ValueError:pass
            else:raise AssertionError('Tampered frozen record accepted: '+key)
        changed=copy.deepcopy(record);changed['original_build_validation_map'].pop(next(iter(changed['original_build_validation_map'])))
        temporary.write_text(json.dumps(changed))
        try:validate(temporary)
        except ValueError:pass
        else:raise AssertionError('Missing compiled dependency accepted')
    finally:temporary.unlink(missing_ok=True)
    print(json.dumps({'passed':True,'native_run':False,'compiled_dependencies_verified':len(origin['files']),'copied_files':len(record['files']),'tamper_checks':4}))
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('manifest',type=Path);args=parser.parse_args();main(args.manifest)
