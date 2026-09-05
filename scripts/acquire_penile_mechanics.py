"""Retain institutional full texts and repository metadata without repinning bytes."""
from pathlib import Path
import datetime
import hashlib
import json
import subprocess
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
PAPERS=[
    (26563843,48384526,'human-penile-mechanics-2024','01bd56edefd308fffc01b743cec34bec'),
    (24271627,42604816,'penile-testing-design-2023','003d7fb3b4ffd12d5f13d337357c7e5c'),
]

def main():
    for article,file_id,label,expected_md5 in PAPERS:
        out=ROOT/'data/raw/biomechanics'/label
        if out.exists():
            receipt=json.loads((out/'receipt.json').read_text())
            for name,identity in receipt['artifacts'].items():
                if hashlib.sha256((out/name).read_bytes()).hexdigest()!=identity:
                    raise ValueError('Retained acquisition changed: '+str(out/name))
            print(json.dumps({'id':label,'status':'held_verified'}));continue
        url=f'https://api.figshare.com/v2/articles/{article}/versions/1'
        with urllib.request.urlopen(url,timeout=60) as response:metadata=response.read()
        detail=json.loads(metadata)
        file=next(f for f in detail['files'] if f['id']==file_id)
        if file['computed_md5']!=expected_md5:raise ValueError('Repository file identity changed')
        with urllib.request.urlopen(file['download_url'],timeout=60) as response:pdf=response.read()
        if not pdf.startswith(b'%PDF') or len(pdf)!=file['size'] or hashlib.md5(pdf).hexdigest()!=expected_md5:
            raise ValueError('Downloaded paper is not the expected repository file')
        out.mkdir(parents=True)
        (out/'metadata.json').write_bytes(metadata)
        (out/'paper.pdf').write_bytes(pdf)
        subprocess.run(['pdftotext','-layout',str(out/'paper.pdf'),str(out/'paper.txt')],check=True)
        receipt={'id':label,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 'repository_metadata_url':url,'download_url':file['download_url'],
                 'repository_file_id':file_id,'repository_version':1,'repository_md5':expected_md5,
                 'license':detail['license'],'extraction_command':'pdftotext -layout paper.pdf paper.txt',
                 'artifacts':{name:hashlib.sha256((out/name).read_bytes()).hexdigest()
                              for name in ('metadata.json','paper.pdf','paper.txt')}}
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps({'id':label,'bytes':len(pdf),'sha256':receipt['artifacts']['paper.pdf']}))

if __name__=='__main__':main()
