"""Acquire one public BIDMC waveform record (bounded <8 MB), preserving raw bytes."""
from pathlib import Path
import hashlib
import json
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
BASE='https://physionet.org/files/bidmc/1.0.0/'

def main():
    out=ROOT/'data/raw/temporal/bidmc'; out.mkdir(parents=True,exist_ok=True)
    files=['bidmc_csv/bidmc_01_Signals.csv','bidmc01.hea','LICENSE','README','SHA256SUMS.txt']
    records=[]
    for name in files:
        path=out/Path(name).name
        if not path.exists():
            with urllib.request.urlopen(BASE+name,timeout=60) as r:
                data=r.read(8_000_001)
            if len(data)>8_000_000: raise ValueError('download exceeded size bound')
            path.write_bytes(data)
        records.append(dict(path=str(path.relative_to(ROOT)),url=BASE+name,bytes=path.stat().st_size,
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    sums=(out/'SHA256SUMS.txt').read_text()
    for name,item in zip(files,records):
        if name=='SHA256SUMS.txt': continue
        matches=[line.split()[0] for line in sums.splitlines() if line.split()[-1].lstrip('./')==name]
        if matches and matches[0]!=item['sha256']: raise ValueError('upstream checksum mismatch: '+name)
        item['upstream_checksum_verified']=bool(matches)
    (out/'provenance.json').write_text(json.dumps(dict(source='BIDMC PPG and Respiration Dataset 1.0.0',
        source_kind='human_measurement',doi='10.13026/C2208R',license='Open Data Commons Attribution License v1.0',
        citation='Pimentel et al. 2016, Towards a Robust Estimation of Respiratory Rate from Pulse Oximeters, DOI 10.1109/TBME.2016.2613124',
        selection='record 01 only; no subject-level generalization',files=records),indent=2))
    print(json.dumps(records,indent=2))

if __name__=='__main__': main()
