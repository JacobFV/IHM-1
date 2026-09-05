"""Acquire pinned Physiome CellML/generated math and execute its original 10-day example."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,sys,urllib.request
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
CELLML='https://models.physiomeproject.org/workspace/schlosser_selgrade_2000/rawfile/73464388213e0eeb5ac101b80c211544db027123/schlosser_selgrade_2000.cellml'
EXPOSURE='https://models.physiomeproject.org/exposure/38fcbccf8eacd4c93e5aa531625a8950/schlosser_selgrade_2000.cellml'
FILES={'schlosser_selgrade_2000.cellml':CELLML,'generated.py':EXPOSURE+'/@@cellml_codegen/Python/raw','codegen.html':EXPOSURE+'/@@cellml_codegen/Python','license_citation.html':EXPOSURE+'/license_citation'}


def collect(root):
    out=Path(root)/'data/raw/reproductive/schlosser_selgrade_2000';out.mkdir(parents=True,exist_ok=True);records=[]
    for name,url in FILES.items():
        path=out/name
        if not path.exists():
            with urllib.request.urlopen(url,timeout=30) as response:data=response.read(250001)
            if len(data)>250000:raise ValueError('source download exceeds bound')
            path.write_bytes(data)
        records.append(dict(path=str(path.relative_to(root)),url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size))
    metadata=dict(source='Schlosser and Selgrade 2000 Physiome CellML',workspace_commit='73464388213e0eeb5ac101b80c211544db027123',
        exposure=EXPOSURE,license='Creative Commons Attribution 3.0 Unported',license_url='https://creativecommons.org/licenses/by/3.0/',
        original_authors=['Paul M. Schlosser','James F. Selgrade'],cellml_author='Catherine Lloyd, University of Auckland',
        citation='A model of gonadotropin regulation during the menstrual cycle in women: qualitative features. Environmental Health Perspectives 108:873–881 (2000), PMID 11035997.',
        acquired_utc=datetime.now(timezone.utc).isoformat(),files=records)
    (out/'provenance.json').write_text(json.dumps(metadata,indent=2));return metadata

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];collect(root)
    from ihm.native.reproductive import build_reproductive
    result=build_reproductive(root);print(json.dumps(result['validation'],indent=2))
