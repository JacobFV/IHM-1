#!/usr/bin/env python3
"""Check actual contact export and fail-closed reader without changing held data."""
from pathlib import Path
import gzip
import hashlib
import json
import sys
import tempfile

BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from ihm.app.experiments import read_experiment


def main():
    data=read_experiment(BASE,'garment-contact')
    assert len(data['time_s'])==49 and data['report']['whole_garment_containment_validated'] is False
    assert data['geometry_precision']['spacing_m']==.004
    assert data['geometry_precision']['source_manifest_sha256']==data['configuration']['tissue_manifest_sha256']
    checks=['held_export_source_and_runtime_hashes']
    with tempfile.TemporaryDirectory(prefix='ihm-contact-reader-') as temporary:
        root=Path(temporary);directory=root/'data/derived/garment-tissue-display-v2';directory.mkdir(parents=True)
        source=root/'source.bin';source.write_bytes(b'retained source')
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        payload=gzip.compress(json.dumps({'source_hashes':{'source.bin':digest}}).encode())
        output=directory/'display.json.gz';output.write_bytes(payload)
        manifest={'display_path':'display.json.gz','display_sha256':hashlib.sha256(payload).hexdigest(),'position_quantization_max_error_m':0}
        (directory/'manifest.json').write_text(json.dumps(manifest))
        assert read_experiment(root,'garment-contact')['source_hashes']['source.bin']==digest
        for label,mutate in [
            ('changed_anatomical_source_rejected',lambda:source.write_bytes(b'changed source')),
            ('changed_display_bytes_rejected',lambda:output.write_bytes(b'changed display')),
            ('escaping_display_path_rejected',lambda:(directory/'manifest.json').write_text(json.dumps({**manifest,'display_path':'../outside.gz'}))),
        ]:
            source.write_bytes(b'retained source');output.write_bytes(payload);(directory/'manifest.json').write_text(json.dumps(manifest))
            mutate()
            try:read_experiment(root,'garment-contact')
            except ValueError:checks.append(label)
            else:raise AssertionError(label)
    report={'passed':True,'checks':checks,'display_identity':data['display_identity'],'frame_count':len(data['time_s'])}
    out=BASE/'artifacts/verification/garment-contact-viewer';out.mkdir(parents=True,exist_ok=True)
    (out/'reader.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
