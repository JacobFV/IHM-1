"""Immutable numerical inputs let historical experiments survive code evolution."""
import json
from pathlib import Path
import shutil
from ihm.native import _sha


def freeze_sources(root, directory, sources):
    root, directory=Path(root).resolve(),Path(directory).resolve()
    if not directory.is_relative_to(root):
        raise ValueError('Retained systemic inputs must be inside the workspace')
    receipt_path=directory/'frozen-sources.json'
    if receipt_path.exists():
        receipt=json.loads(receipt_path.read_text())
        if receipt['original_sources']!=sources:
            raise ValueError('An input snapshot cannot be repinned')
        return resolve_sources(root,directory,sources)
    frozen={}
    for name,expected in sources.items():
        source=(root/name).resolve()
        if not source.is_relative_to(root) or _sha(source)!=expected:
            raise ValueError('Cannot archive changed systemic input: '+name)
        target=directory/'inputs'/name
        if not target.resolve().is_relative_to(directory/'inputs'):
            raise ValueError('Invalid input archive path')
        target.parent.mkdir(parents=True,exist_ok=True)
        temporary=target.with_name(target.name+'.copying')
        shutil.copyfile(source,temporary)
        if _sha(temporary)!=expected:
            temporary.unlink()
            raise ValueError('Systemic input changed during archiving: '+name)
        temporary.replace(target)
        frozen[name]=dict(path=str(target.relative_to(root)),sha256=expected)
    receipt=dict(schema='ihm.systemic-inputs.v1',original_sources=sources,copies=frozen,
                 interpretation='Exact input bytes; preserving historical execution does not establish biological validity')
    temporary=receipt_path.with_suffix('.writing')
    temporary.write_text(json.dumps(receipt,indent=2)+'\n')
    temporary.replace(receipt_path)
    return resolve_sources(root,directory,sources)


def resolve_sources(root,directory,sources):
    root, directory=Path(root).resolve(),Path(directory).resolve()
    receipt_path=directory/'frozen-sources.json'
    if not receipt_path.exists():
        resolved={}
        for name,expected in sources.items():
            path=(root/name).resolve()
            if not path.is_relative_to(root) or _sha(path)!=expected:
                raise ValueError('Systemic source changed without retained input bytes: '+name)
            resolved[name]=expected
        return resolved
    receipt=json.loads(receipt_path.read_text())
    if receipt.get('schema')!='ihm.systemic-inputs.v1' or receipt.get('original_sources')!=sources or set(receipt.get('copies',{}))!=set(sources):
        raise ValueError('Systemic input snapshot identities differ')
    resolved={str(receipt_path.relative_to(root)):_sha(receipt_path)}
    for name,expected in sources.items():
        copy=receipt['copies'][name]
        path=(root/copy['path']).resolve()
        if not path.is_relative_to(directory/'inputs') or copy['sha256']!=expected or _sha(path)!=expected:
            raise ValueError('Retained systemic input changed: '+name)
        resolved[str(path.relative_to(root))]=expected
    return resolved
