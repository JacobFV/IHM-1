"""Source-bound native pose evaluation journal; exact records, no interpolation."""
import hashlib,json,math
from pathlib import Path


def pose_key(values):
    values=list(map(float,values))
    if not all(map(math.isfinite,values)):raise ValueError('Nonfinite cached coordinate')
    return tuple(value.hex() for value in values)


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


class PoseJournal:
    def __init__(self,directory,identity,resume=None):
        self.directory=Path(directory);self.directory.mkdir(parents=True,exist_ok=True)
        self.identity=identity;self.cache={}
        (self.directory/'cache_identity.json').write_text(json.dumps(identity,indent=2,allow_nan=False)+'\n')
        if resume is not None:
            resume=Path(resume)
            if json.loads((resume/'cache_identity.json').read_text())!=identity:
                raise ValueError('Static cache source/protocol identity mismatch')
            for line in (resume/'candidates.jsonl').read_text().splitlines():
                record=json.loads(line);payload=record['payload']
                if record['sha256']!=digest(payload):raise ValueError('Static cache record hash mismatch')
                self.append(payload['values'],payload['entry'])

    def append(self,values,entry):
        key=pose_key(values)
        if key in self.cache:
            if self.cache[key]!=entry:raise ValueError('Conflicting native responses for exact cached pose')
            return
        payload={'values':list(map(float,values)),'entry':entry}
        record={'payload':payload,'sha256':digest(payload)}
        # One flushed line per actual response survives the bounded solver exception.
        with (self.directory/'candidates.jsonl').open('a') as stream:
            stream.write(json.dumps(record,allow_nan=False)+'\n');stream.flush()
        self.cache[key]=entry
