"""Bounded primary CSF-source acquisition; preserve challenge evidence and honest gaps."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,sys,urllib.request,urllib.error
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.csf import classify_source_response,build_csf_status,csf_source_status,SOURCE_URL,REVIEW_SHA256
SOURCES=[
    ('original-full.html','https://journals.physiology.org/doi/full/10.1152/jappl.1997.82.4.1256','html','requested_original'),
    ('original-pdf.response','https://journals.physiology.org/doi/pdf/10.1152/jappl.1997.82.4.1256','pdf','requested_original'),
    ('companion-clinical.pdf','https://air.unimi.it/bitstream/2434/243920/2/1270.full.pdf','pdf','author_companion_not_complete_model'),
    ('pubmed.html','https://pubmed.ncbi.nlm.nih.gov/9104864/','html','bibliographic_metadata')]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--offline',action='store_true',help='Build using already acquired, hash-verified original source review');args=parser.parse_args()
    root=Path(__file__).resolve().parents[1];out=root/'data/raw/csf/ursino_lodi_1997';out.mkdir(parents=True,exist_ok=True);attempts=[]
    if args.offline:
        status=csf_source_status(root)
        if status['execution_available']:status=build_csf_status(root)
        else:(root/'data/derived/csf').mkdir(parents=True,exist_ok=True);(root/'data/derived/csf/index.json').write_text(json.dumps(dict(schema_version=1,models=[status])))
        print(json.dumps(dict(execution_available=status['execution_available'],models=len(status.get('models',[]))),indent=2));return
    for name,url,expected,role in SOURCES:
        path=out/name;http_status=None
        try:
            with urllib.request.urlopen(url,timeout=30) as response:
                data=response.read(2_000_001);http_status=response.status
        except urllib.error.HTTPError as error:
            http_status=error.code;data=error.read(2_000_001)
        except (OSError,TimeoutError) as error:
            attempts.append(dict(url=url,role=role,status='transport_error',error=str(error)));continue
        if len(data)>2_000_000:
            attempts.append(dict(url=url,role=role,status='size_bound_exceeded',http_status=http_status));continue
        # Preserve successive response bytes without replacing prior acquisition.
        sha=hashlib.sha256(data).hexdigest()
        if path.exists() and path.read_bytes()!=data:path=out/(path.stem+'-'+sha[:12]+path.suffix)
        if not path.exists():path.write_bytes(data)
        attempts.append(dict(url=url,role=role,status=classify_source_response(data,expected),http_status=http_status,
            path=str(path.relative_to(root)),sha256=sha,bytes=len(data)))
    manifest=dict(acquired_utc=datetime.now(timezone.utc).isoformat(),attempts=attempts,
        copyright='Original and companion articles copyright 1997 American Physiological Society; repository availability is not an open reuse license.',
        primary_companion_citation='Ursino, Lodi, Rossi and Stocchetti. Intracranial pressure dynamics in patients with acute brain damage. J Appl Physiol 82:1270–1282 (1997).',
        primary_companion_repository='https://air.unimi.it/handle/2434/243920',
        discovery_audit='Original publisher challenges preserved. Original author-upload indexed full text supplied Table 1 and Appendices A-C; implementation transcribes reviewed equations, not official executable code.',
        reviewed_original=dict(url=SOURCE_URL,path='data/raw/csf/ursino_lodi_1997/indexed-author-source-review.json',sha256=REVIEW_SHA256,hash_scope='indexed retrieval record, not original PDF',acquisition='Web indexed author-upload retrieval; retained separately from HTTP challenge bodies'))
    (out/'acquisition.json').write_text(json.dumps(manifest,indent=2));status=csf_source_status(root)
    if status['execution_available']:status=build_csf_status(root)
    else:
        target=root/'data/derived/csf';target.mkdir(parents=True,exist_ok=True)
        (target/'index.json').write_text(json.dumps(dict(schema_version=1,models=[status]),indent=2))
    print(json.dumps(dict(execution_available=status['execution_available'],attempts=attempts,missing=status.get('missing_evidence',[])),indent=2))

if __name__=='__main__':main()
