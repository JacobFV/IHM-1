"""reproducible acquisition of public scientific artifacts, not automatic calibration."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def fetch(url, path, *, max_bytes=512*1024**2, prefix=None):
    path = Path(path); receipt = path.with_name(path.name+'.receipt.json')
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and receipt.exists():
        record = json.loads(receipt.read_text())
        if record['url'] == url and path.stat().st_size == record['bytes'] and sha256(path) == record['sha256']:
            if record['bytes']>max_bytes: raise ValueError('cached asset exceeds current size budget')
            if prefix is not None:
                with path.open('rb') as stream:
                    if not stream.read(len(prefix)).startswith(prefix): raise ValueError('cached asset fails expected header')
            return record
    temporary = path.with_name(path.name+'.partial')
    request = urllib.request.Request(url, headers={'User-Agent': 'IHM-1 scientific-data-acquisition/0.2'})
    try:
        with urllib.request.urlopen(request, timeout=90) as response, temporary.open('wb') as stream:
            size = 0; digest = hashlib.sha256(); first = True
            while block := response.read(1024*1024):
                if first and prefix is not None and not block.startswith(prefix):
                    raise ValueError(f'unexpected content header for {url}')
                first = False; size += len(block)
                if size > max_bytes:
                    raise ValueError(f'asset exceeds acquisition budget: {url}')
                digest.update(block); stream.write(block)
            if size == 0:
                raise ValueError(f'empty asset: {url}')
            record = {'url': url, 'resolved_url': response.url, 'path': str(path), 'bytes': size,
                      'sha256': digest.hexdigest(), 'content_type': response.headers.get('Content-Type'),
                      'retrieved_at_utc': datetime.now(timezone.utc).isoformat(), 'status': 'acquired'}
        temporary.replace(path)
        receipt.write_text(json.dumps(record, indent=2)+'\n')
        return record
    finally:
        temporary.unlink(missing_ok=True)
