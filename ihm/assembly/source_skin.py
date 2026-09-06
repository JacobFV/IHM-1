"""Bounded exact-index source skin patches with explicitly supplied owners.

This does not infer anatomical owners or search for nearest faces. The reader
selects numbers from retained hash-bound flat JSON arrays without materializing
the complete source mesh. It supports the canonical geometry serialization.
"""
from dataclasses import dataclass
from pathlib import Path
import gzip,hashlib,io,json,re
import numpy as np


def file_sha256(path):
    if isinstance(path,bytes):return hashlib.sha256(path).hexdigest()
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(65536),b''):digest.update(chunk)
    return digest.hexdigest()


def _chunks(path):
    total=0
    if isinstance(path,bytes):
        binary=io.BytesIO(path);stream=gzip.GzipFile(fileobj=binary) if path.startswith(b'\x1f\x8b') else binary
        handle=io.TextIOWrapper(stream,encoding='utf8')
    else:
        opener=gzip.open if str(path).endswith('.gz') else open;handle=opener(path,'rt',encoding='utf8')
    with handle:
        for chunk in iter(lambda:handle.read(65536),''):
            total+=len(chunk)
            if total>256*1024*1024:raise ValueError('Source JSON exceeds bounded scan budget')
            yield chunk


def _member_start(path,key):
    pattern=re.compile(r'"'+re.escape(key)+r'"\s*:\s*');tail=''
    iterator=iter(_chunks(path))
    for chunk in iterator:
        text=tail+chunk;match=pattern.search(text)
        if match:return text[match.end():],iterator
        tail=text[-128:]
    raise ValueError('Missing retained JSON member: '+key)


def read_small_member(path,key,*,max_chars=2*1024*1024):
    """Read a small retained member while skipping large sibling render arrays."""
    text,iterator=_member_start(path,key);decoder=json.JSONDecoder()
    while True:
        try:
            value,end=decoder.raw_decode(text.lstrip())
            if end>max_chars:raise ValueError('Selected JSON member exceeds budget')
            return value
        except json.JSONDecodeError:
            if len(text)>max_chars:raise ValueError('Selected JSON member exceeds budget')
            try:text+=next(iterator)
            except StopIteration:raise ValueError('Truncated selected JSON member') from None


def _selected_numbers(path,key,indices,*,integer=False):
    selected=set(map(int,indices));output={};text,iterator=_member_start(path,key)
    text=text.lstrip()
    if not text.startswith('['):raise ValueError('Expected canonical flat source array')
    text=text[1:];tail='';index=0
    while True:
        complete=']' in text
        if complete:text=text.split(']',1)[0]
        pieces=(tail+text).split(',');tail='' if complete else pieces.pop()
        for token in pieces:
            token=token.strip()
            if not token and complete and index==0:continue
            if index in selected:
                if integer and not re.fullmatch(r'0|[1-9][0-9]*',token):raise ValueError('Invalid source node index')
                try:value=int(token) if integer else float(token)
                except ValueError:raise ValueError('Invalid retained source coordinate') from None
                if not np.isfinite(value):raise ValueError('Nonfinite source value')
                output[index]=value
            index+=1
        if complete:break
        try:text=next(iterator)
        except StopIteration:raise ValueError('Truncated source array') from None
    if set(output)!=selected:raise ValueError('Requested source index is absent')
    return output,index


def extract_source_patch(path,expected_sha256,source_face_indices,*,source_id):
    if not isinstance(source_id,str) or not source_id or not isinstance(expected_sha256,str) or not re.fullmatch('[0-9a-f]{64}',expected_sha256):raise ValueError('Explicit source ID and SHA-256 required')
    if file_sha256(path)!=expected_sha256:raise ValueError('Source geometry hash mismatch')
    faces=np.asarray(source_face_indices)
    if faces.ndim!=1 or faces.dtype.kind not in 'iu' or not 0<len(faces)<=4096 or np.any(faces<0):raise ValueError('Bounded nonnegative source face indices required')
    faces=np.unique(faces);flat=[3*int(f)+k for f in faces for k in range(3)]
    values,total_indices=_selected_numbers(path,'indices',flat,integer=True)
    if total_indices%3 or int(faces.max())>=total_indices//3:raise ValueError('Invalid source face array or index')
    source_triangles=np.array([[values[3*int(f)+k] for k in range(3)] for f in faces],int)
    nodes=np.unique(source_triangles);flat=[int(3*n+k) for n in nodes for k in range(3)]
    values,total_coordinates=_selected_numbers(path,'positions',flat)
    if total_coordinates%3:raise ValueError('Invalid source coordinate array length')
    positions=np.array([[values[int(3*n+k)] for k in range(3)] for n in nodes])
    if file_sha256(path)!=expected_sha256:raise ValueError('Source geometry changed during bounded extraction')
    return {'source_id':source_id,'source_sha256':expected_sha256,'source_node_indices':nodes.tolist(),'source_face_indices':faces.tolist(),
            'positions_m':positions.tolist(),'triangles':np.searchsorted(nodes,source_triangles).tolist(),
            'source_vertex_count':total_coordinates//3,'source_face_count':total_indices//3,
            'coverage':'Exact retained requested faces only; not full source contact coverage'}


@dataclass(frozen=True)
class SourceSkinPatch:
    source_id:str
    source_sha256:str
    registration_sha256:str
    owner_basis:str
    positions_m:np.ndarray
    triangles:np.ndarray
    source_node_indices:np.ndarray
    source_face_indices:np.ndarray
    owner_names:tuple
    node_owner_indices:np.ndarray

    @classmethod
    def bind(cls,patch,binding,*,expected_registration_sha256,owner_names):
        if not 0<len(patch['source_node_indices'])<=12288 or not 0<len(patch['source_face_indices'])<=4096:raise ValueError('Source patch correspondence exceeds bounded size')
        if binding.get('source_id')!=patch['source_id'] or binding.get('source_sha256')!=patch['source_sha256']:raise ValueError('Owner/source identity mismatch')
        if not re.fullmatch('[0-9a-f]{64}',expected_registration_sha256) or binding.get('registration_sha256')!=expected_registration_sha256:raise ValueError('Owner registration hash mismatch')
        if not isinstance(binding.get('basis'),str) or not binding['basis'].strip():raise ValueError('Explicit owner assignment basis required')
        names=tuple(owner_names)
        if not names or len(set(names))!=len(names) or any(not isinstance(n,str) or not n for n in names):raise ValueError('Unique registered owner names required')
        rows=binding.get('nodes',[]);ids=[r['source_node_index'] for r in rows]
        if any(isinstance(i,bool) or not isinstance(i,int) for i in ids) or len(set(ids))!=len(ids) or set(ids)!=set(patch['source_node_indices']):raise ValueError('Exact unambiguous owner for every retained source node required')
        by_id={r['source_node_index']:r['owner'] for r in rows}
        if any(owner not in names for owner in by_id.values()):raise ValueError('Unknown registered source owner')
        x=np.array(patch['positions_m'],float);tri=np.asarray(patch['triangles']);nodes=np.asarray(patch['source_node_indices']);faces=np.asarray(patch['source_face_indices'])
        if x.shape!=(len(nodes),3) or not np.isfinite(x).all() or tri.shape!=(len(faces),3) or tri.dtype.kind not in 'iu' or tri.min()<0 or tri.max()>=len(x):raise ValueError('Invalid exact source patch geometry')
        if nodes.dtype.kind not in 'iu' or faces.dtype.kind not in 'iu' or len(np.unique(nodes))!=len(nodes) or len(np.unique(faces))!=len(faces) or np.any(nodes<0) or np.any(faces<0):raise ValueError('Ambiguous source correspondence')
        owners=np.array([names.index(by_id[int(i)]) for i in nodes],int)
        arrays=[x,tri.copy(),nodes.copy(),faces.copy(),owners]
        for a in arrays:a.flags.writeable=False
        return cls(patch['source_id'],patch['source_sha256'],expected_registration_sha256,binding['basis'],*arrays[:4],names,arrays[4])
