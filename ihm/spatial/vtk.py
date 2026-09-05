"""Read original VTK XML numeric arrays with bounded per-array decompression.

The acquired VMR/OpenSim files use ASCII or inline/appended base64 with independent
VTK zlib blocks. Raw appended binary and other compressors fail explicitly.
"""
import base64
import mmap
from pathlib import Path
import xml.etree.ElementTree as ET
import zlib
import numpy as np

TYPES={'Float32':'f4','Float64':'f8','Int32':'i4','Int64':'i8','UInt8':'u1','UInt32':'u4','UInt64':'u8','Int8':'i1','Int16':'i2','UInt16':'u2'}

def decode(node,root,encoded=None):
    endian='<' if root.get('byte_order','LittleEndian')=='LittleEndian' else '>'
    dtype=np.dtype(endian+TYPES[node.get('type')]);fmt=node.get('format','ascii')
    if fmt=='ascii':a=np.fromstring(node.text or '',sep=' ',dtype=dtype)
    else:
        encoded=b''.join(encoded.split()) if isinstance(encoded,bytes) else ''.join((encoded or node.text or '').split())
        hdtype=np.dtype(endian+TYPES[root.get('header_type','UInt32')]);word=hdtype.itemsize
        if root.get('compressor'):
            if root.get('compressor')!='vtkZLibDataCompressor':raise ValueError('unsupported VTK compressor')
            first=base64.b64decode(encoded[:((word+2)//3)*4]);n=int(np.frombuffer(first[:word],dtype=hdtype)[0])
            if n>10_000_000:raise ValueError('invalid VTK compression block count')
            chars=((word*(3+n)+2)//3)*4
            header=np.frombuffer(base64.b64decode(encoded[:chars])[:word*(3+n)],dtype=hdtype)
            payload=base64.b64decode(encoded[chars:]);blocks=[];offset=0
            for size in header[3:]:
                size=int(size);blocks.append(zlib.decompress(payload[offset:offset+size]));offset+=size
            raw=b''.join(blocks)
            expected=(int(header[0])-1)*int(header[1])+int(header[2]) if n else 0
            if len(raw)!=expected:raise ValueError('VTK decompressed length mismatch')
        else:
            raw=base64.b64decode(encoded);length=int(np.frombuffer(raw[:word],dtype=hdtype)[0])
            raw=base64.b64decode(encoded[((word+2)//3)*4:])[:length] if len(raw)<length+word else raw[word:word+length]
            if len(raw)!=length:raise ValueError('VTK uncompressed length mismatch')
        a=np.frombuffer(raw,dtype=dtype).copy()
    components=int(node.get('NumberOfComponents','1'))
    if components<1 or len(a)%components:raise ValueError('VTK component count mismatch')
    a=a.reshape(-1,components) if components>1 else a
    if not np.isfinite(a).all():raise ValueError('nonfinite original VTK array')
    return a

def read_arrays(path,names=None):
    """Return association/name arrays; optional names avoids decoding unrelated time states."""
    with Path(path).open('rb') as file,mmap.mmap(file.fileno(),0,access=mmap.ACCESS_READ) as mm:
        start=mm.find(b'<AppendedData')
        root=ET.fromstring(mm[:start]+b'</VTKFile>') if start>=0 else ET.fromstring(mm[:])
        if start>=0:
            end=mm.find(b'>',start)
            if b'encoding="base64"' not in mm[start:end]:raise ValueError('raw appended VTK not supported')
            payload_start=mm.find(b'_',end)+1
            payload_end=mm.rfind(b'</AppendedData>')
        nodes=[]
        for piece in root.iter('Piece'):
            for section in piece:
                for node in section.findall('DataArray'):nodes.append((section.tag+'/'+node.get('Name','Points' if section.tag=='Points' else 'unnamed'),node))
        offsets=sorted(int(n.get('offset')) for _,n in nodes if n.get('offset') is not None)
        ends=dict(zip(offsets,offsets[1:]));result={}
        for key,node in nodes:
            if names is not None and key not in names and key.split('/',1)[0]+'/*' not in names:continue
            if key in result:raise ValueError('multiple VTK pieces need explicit merge before import')
            encoded=None
            if node.get('format')=='appended':
                if start<0:raise ValueError('missing appended payload')
                offset=int(node.get('offset'));stop=payload_start+ends[offset] if offset in ends else payload_end
                encoded=mm[payload_start+offset:stop]
            result[key]=decode(node,root,encoded)
        return result

def surface(path):
    d=read_arrays(path,names={'Points/*','Polys/connectivity','Polys/offsets'})
    point_arrays=[v for k,v in d.items() if k.startswith('Points/')]
    if len(point_arrays)!=1:raise ValueError('exactly one point coordinate array required')
    points=point_arrays[0];conn=d['Polys/connectivity'];offsets=d['Polys/offsets']
    if points.ndim!=2 or points.shape[1]!=3:raise ValueError('surface points must be 3D')
    if offsets.ndim!=1 or conn.ndim!=1 or not np.issubdtype(offsets.dtype,np.integer) or not np.issubdtype(conn.dtype,np.integer):raise ValueError('integer polygon connectivity required')
    if len(offsets) and (np.any(np.diff(np.r_[0,offsets])<3) or offsets[-1]!=len(conn)):raise ValueError('invalid polygon offsets')
    faces=[];start=0
    for stop in offsets:
        polygon=conn[start:int(stop)];start=int(stop)
        if len(polygon)<3:raise ValueError('invalid polygon')
        faces.extend((int(polygon[0]),int(polygon[j]),int(polygon[j+1])) for j in range(1,len(polygon)-1))
    f=np.asarray(faces,dtype=np.int64)
    if f.size and (f.min()<0 or f.max()>=len(points)):raise ValueError('invalid surface index')
    return points,f
