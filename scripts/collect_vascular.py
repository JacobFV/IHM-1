#!/usr/bin/env python3
"""Acquire and audit original VMR human geometry, boundary inputs and CFD results.

No solver runs, inferred velocity fields, unit conversion, or patient registration.
Original ZIPs are authoritative; index-only is offline and checks ZIP CRCs + SHA256.
"""
from __future__ import annotations
import argparse, base64, csv, hashlib, io, json, re, urllib.request, zipfile, zlib
from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

REVISION = '7a1c17dfc9ce4b34eec587fd12f709e8e7baea4e'
CASES = ('0001_H_AO_SVD', '0050_H_CERE_H', '0077_H_PULM_H')
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / 'data/raw/vascular/vmr'
DERIVED = PROJECT_ROOT / 'data/derived/vascular'
ASSETS = ['svprojects/' + c + '.zip' for c in CASES] + [
    'svresults/0050_H_CERE_H/0050_H_CERE_H_3D_RIGID_VTU.zip']
METADATA = ['dataset/dataset-svprojects.csv', 'dataset/dataset-svresults.csv',
            'dataset/file_sizes.csv', 'LICENSE', 'js/globalVar.js']
UNIT_NOTE = ('Native solver values retained. CGS is consistent with density 1.06 and '
             'viscosity 0.04, but file-specific units are not declared in VTK/flow files; '
             'do not silently treat this inference as confirmed metadata.')

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''): h.update(b)
    return h.hexdigest()

def fetch(url, destination, max_bytes=5_000_000_000):
    if destination.exists(): return
    partial = destination.with_suffix(destination.suffix + '.partial')
    with urllib.request.urlopen(url, timeout=90) as response, partial.open('wb') as out:
        total = 0
        while data := response.read(1024 * 1024):
            total += len(data)
            if total > max_bytes: raise ValueError('Download exceeds configured byte limit')
            out.write(data)
    partial.replace(destination)

def decode_array(node, root):
    """Decode inline VTK XML arrays without altering numbers or coordinates."""
    types = {'Float32':'f4','Float64':'f8','Int32':'i4','Int64':'i8',
             'UInt8':'u1','UInt32':'u4','UInt64':'u8','Int8':'i1','Int16':'i2','UInt16':'u2'}
    endian = '<' if root.get('byte_order', 'LittleEndian') == 'LittleEndian' else '>'
    dtype = np.dtype(endian + types[node.get('type')])
    fmt = node.get('format', 'ascii')
    if fmt == 'ascii': a = np.fromstring(node.text or '', sep=' ', dtype=dtype)
    elif fmt == 'binary':
        encoded = ''.join((node.text or '').split())
        hdtype = np.dtype(endian + types[root.get('header_type', 'UInt32')])
        word = hdtype.itemsize
        if root.get('compressor'):
            if root.get('compressor') != 'vtkZLibDataCompressor': raise ValueError('Unsupported compressor')
            first = base64.b64decode(encoded[:((word + 2)//3)*4])
            nblocks = int(np.frombuffer(first[:word], dtype=hdtype)[0])
            nchars = ((word * (3 + nblocks) + 2)//3)*4
            header = np.frombuffer(base64.b64decode(encoded[:nchars])[:word*(3+nblocks)], dtype=hdtype)
            compressed = base64.b64decode(encoded[nchars:])
            blocks, offset = [], 0
            for size in header[3:]:
                size = int(size); blocks.append(zlib.decompress(compressed[offset:offset+size])); offset += size
            raw = b''.join(blocks)
        else:
            raw = base64.b64decode(encoded)
            size = int(np.frombuffer(raw[:word], dtype=hdtype)[0])
            if len(raw) < size + word:
                nchars = ((word+2)//3)*4
                raw = base64.b64decode(encoded[nchars:])[:size]
            else: raw = raw[word:word+size]
        a = np.frombuffer(raw, dtype=dtype)
    else: raise ValueError('Array encoding requires a VTK reader: ' + fmt)
    components = int(node.get('NumberOfComponents', '1'))
    return a.reshape(-1, components) if components > 1 else a

def vtk_index(data, decode=False):
    # Appended raw payload need not be XML parseable: metadata lives before it.
    if b'<AppendedData' in data:
        prefix = data.split(b'<AppendedData', 1)[0]
        root = ET.fromstring(prefix + b'</VTKFile>')
    else: root = ET.fromstring(data)
    pieces, arrays, decoded = [], [], {}
    appended_start = None
    offsets = sorted(int(n.get('offset')) for n in root.iter('DataArray') if n.get('offset') is not None)
    offset_ends = dict(zip(offsets, offsets[1:]))
    if b'<AppendedData' in data:
        tag_start = data.index(b'<AppendedData')
        tag_end = data.index(b'>', tag_start)
        if b'encoding="base64"' in data[tag_start:tag_end]:
            appended_start = data.index(b'_', tag_end) + 1
    steps = sorted({int(m.group(1)) for n in root.iter('DataArray')
                    if (m := re.search(r'_(\d+)$', n.get('Name','')))})
    for piece in root.iter('Piece'):
        pieces.append({k: int(v) for k,v in piece.attrib.items() if k.startswith('NumberOf')})
        for section in piece:
            for node in section.findall('DataArray'):
                entry = {'association': section.tag, **node.attrib}
                arrays.append(entry)
                temporal_match = re.search(r'_(\d+)$', node.get('Name',''))
                first_time = not temporal_match or not steps or int(temporal_match.group(1)) == steps[0]
                if decode and first_time and section.tag in ('Points','PointData','CellData','Cells'):
                    try:
                        if node.get('format') == 'appended' and appended_start is not None:
                            offset = int(node.get('offset'))
                            end = offset_ends.get(offset)
                            payload = data[appended_start+offset:appended_start+end] if end is not None else data[appended_start+offset:data.rindex(b'</AppendedData>')]
                            inline = ET.Element('DataArray', {**node.attrib, 'format':'binary'})
                            inline.text = payload.decode('ascii').strip()
                            a = decode_array(inline, root)
                        else:
                            a = decode_array(node, root)
                        key = section.tag + '/' + node.get('Name', 'unnamed')
                        decoded[key] = a
                        entry['decoded_shape'] = list(a.shape)
                        entry['all_finite'] = bool(np.isfinite(a).all())
                        if a.size:
                            entry['min'] = np.min(a, axis=0).tolist()
                            entry['max'] = np.max(a, axis=0).tolist()
                    except (ValueError, KeyError, zlib.error) as exc: entry['decode_error'] = str(exc)
    return {'vtk_type': root.get('type'), 'pieces': pieces, 'arrays': arrays, 'temporal_step_suffixes':steps, 'decoded_temporal_step':steps[0] if decode and steps else None}, decoded

def props(root):
    return {p.get('key'): p.get('value') for p in root.findall('prop')}

def job_index(data):
    text = re.sub(r'<\?xml[^>]*\?>', '', data.decode())
    doc = ET.fromstring('<wrapper>' + text + '</wrapper>')
    job = doc.find('mitk_job')
    if job is None: return {}
    result = {'job_attributes': job.attrib, 'sections': {}, 'caps': []}
    for section in job.find('job'):
        if section.tag == 'cap_props':
            result['caps'] = [{'name': cap.get('name'), 'properties': props(cap)} for cap in section]
        else: result['sections'][section.tag] = props(section)
    return result

def inventory(require_all=True):
    DERIVED.mkdir(parents=True, exist_ok=True)
    catalog = {r['Name']:r for r in csv.DictReader((RAW/'dataset-svprojects.csv').open())}
    sizes = {r['Name']:int(r['Size']) for r in csv.DictReader((RAW/'file_sizes.csv').open())}
    report = {'schema_version':1, 'generated_utc':datetime.now(timezone.utc).isoformat(),
              'catalog_revision':REVISION, 'unit_policy':UNIT_NOTE, 'cases':[], 'archives':[],
              'caveats':['Separate donors and anatomy regions; no whole-body registration.',
                         'CFD outputs are simulated, not measured 4D flow MRI.',
                         'Boundary input waveforms are processed solver inputs; original scanner acquisition arrays are not identified.',
                         'Simulation job status is retained verbatim, including failed status.',
                         'Initial-condition and prescribed inlet arrays are not counted as simulated result frames.']}
    report['metadata_sources'] = [
        {'file': str((RAW/Path(path).name).relative_to(PROJECT_ROOT)),
         'url':'https://raw.githubusercontent.com/SimVascular/vascularmodel/'+REVISION+'/'+path,
         'revision':REVISION, 'sha256':sha256(RAW/Path(path).name),
         'bytes':(RAW/Path(path).name).stat().st_size} for path in METADATA]
    report['license'] = {'name':'VMR research and development data license',
                         'path':str((RAW/'LICENSE').relative_to(PROJECT_ROOT)),
                         'spdx':None, 'preserve_copyright':True}
    case_records = {}
    for case in CASES:
        r = {'case_id':case, 'source_metadata':catalog[case], 'geometry':[], 'boundary_conditions':[], 'flow_series':[], 'result_frames':[]}
        report['cases'].append(r); case_records[case]=r
    archive_members = []
    for source_path in ASSETS:
        archive = RAW/Path(source_path).name
        if not archive.exists():
            if require_all: raise FileNotFoundError(archive)
            continue
        results = source_path.startswith('svresults/')
        case = source_path.split('/')[1] if results else archive.stem
        target = case_records[case]
        if archive.stat().st_size != sizes[source_path]: raise ValueError('Catalog length mismatch: '+str(archive))
        print('Auditing', archive.name, flush=True)
        ar = {'file':str(archive.relative_to(PROJECT_ROOT)), 'url':'https://www.vascularmodel.com/'+source_path,
              'bytes':archive.stat().st_size, 'sha256':sha256(archive), 'catalog_bytes':sizes[source_path]}
        report['archives'].append(ar)
        with zipfile.ZipFile(archive) as z:
            # Every member is read: CRC verification happens during decompression.
            for info in z.infolist():
                if info.is_dir(): continue
                n=info.filename
                if '__MACOSX' in n or '/._' in n: continue
                data=z.read(info)
                entry={'archive':archive.name, 'member':n,'bytes':len(data), 'sha256':hashlib.sha256(data).hexdigest(), 'zip_crc32':format(info.CRC,'08x')}
                archive_members.append(entry)
                path=Path(n)
                if path.is_absolute() or '..' in path.parts: raise ValueError('Unsafe archive member')
                # Preserve all model assets; keep large result series in original archive.
                if not results or (path.suffix=='.vtu' and len(target['result_frames'])==0):
                    out=RAW/'extracted'/path
                    out.parent.mkdir(parents=True, exist_ok=True)
                    if not out.exists(): out.write_bytes(data)
                    entry['extracted_path']=str(out.relative_to(PROJECT_ROOT))
                if path.suffix in ('.vtu','.vtp'):
                    important = results and len(target['result_frames'])==0 or '/Meshes/' in n and path.suffix=='.vtu'
                    meta, decoded = vtk_index(data, decode=important)
                    meta.update(entry)
                    meta['evidence_kind']='simulated_cfd_result' if results else ('prescribed_boundary_profile' if path.name=='bct.vtp' else 'initial_condition' if 'initial.' in n else 'image_derived_geometry_or_mesh')
                    if results:
                        match=re.search(r'(\d+)\.(?:vtu|vtp)$',n)
                        meta['solver_step_from_filename']=int(match.group(1)) if match else None
                        meta['physical_time_seconds']=None
                        target['result_frames'].append(meta)
                    else: target['geometry'].append(meta)
                    if decoded and important:
                        name=case + ('_first_cfd_frame.npz' if results else '_mesh.npz')
                        np.savez_compressed(DERIVED/name, **decoded)
                        meta['numeric_arrays_path']=str((DERIVED/name).relative_to(PROJECT_ROOT))
                elif path.suffix=='.flow':
                    try:
                        a=np.loadtxt(io.BytesIO(data),ndmin=2)
                        if a.shape[1]!=2 or not np.isfinite(a).all(): raise ValueError('Expected finite two-column flow')
                        series={**entry,'evidence_kind':'processed_solver_boundary_input','measurement_status':'not_raw_measurement',
                                'time_unit':'native; expected seconds, unconfirmed in file', 'flow_unit':'native; expected cm^3/s under CGS, unconfirmed in file',
                                'samples':a.shape[0],'time_range':a[:,0][[0,-1]].tolist(),'flow_range':[float(a[:,1].min()),float(a[:,1].max())],
                                'strictly_increasing_time':bool(np.all(np.diff(a[:,0])>0)), 'values':a.tolist()}
                        target['flow_series'].append(series)
                    except ValueError as exc: target['flow_series'].append({**entry,'parse_error':str(exc)})
                elif path.suffix=='.sjb':
                    target['boundary_conditions'].append({**entry,'kind':'simvascular_job','data':job_index(data)})
                elif path.suffix=='.svpre' or path.name in ('solver.inp','resistance.dat','rcrt.dat','cort.dat'):
                    target['boundary_conditions'].append({**entry,'kind':path.name,'text':data.decode(errors='replace')})
        ar['zip_crc_verified']=True
    report['archive_member_count']=len(archive_members)
    report['archive_member_bytes']=sum(x['bytes'] for x in archive_members)
    (DERIVED/'vmr_archive_members.json').write_text(json.dumps(archive_members,indent=2)+'\n')
    report['local_raw_file_count']=sum(1 for p in RAW.rglob('*') if p.is_file() and not p.name.endswith('.partial'))
    report['local_raw_bytes']=sum(p.stat().st_size for p in RAW.rglob('*') if p.is_file() and not p.name.endswith('.partial'))
    (DERIVED/'vmr_index.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('local_raw_file_count','local_raw_bytes','archive_member_count','archive_member_bytes')},indent=2))

def verify_assets():
    checks = []
    for case in CASES:
        asset = DERIVED/(case+'_mesh.npz')
        if not asset.exists(): continue
        with np.load(asset) as z:
            pts, con, off, typ = (z[k] for k in ('Points/Points','Cells/connectivity','Cells/offsets','Cells/types'))
            record = {'case_id':case, 'points':len(pts), 'cells':len(off),
                      'coordinates_finite':bool(np.isfinite(pts).all()),
                      'tetrahedral_cells_only':bool(np.all(typ==10)),
                      'all_cells_have_four_vertices':bool(np.all(np.diff(np.r_[0,off])==4)),
                      'connectivity_indices_valid':bool(con.min()>=0 and con.max()<len(pts))}
            assert all(record[k] for k in ('coordinates_finite','tetrahedral_cells_only','all_cells_have_four_vertices','connectivity_indices_valid'))
            checks.append(record)
    result = {'mesh_checks':checks}
    field_asset = DERIVED/'0050_H_CERE_H_first_cfd_frame.npz'
    if field_asset.exists():
        with np.load(field_asset) as f, np.load(DERIVED/'0050_H_CERE_H_mesh.npz') as m:
            a, b = f['Cells/connectivity'].reshape(-1,4), m['Cells/connectivity'].reshape(-1,4)
            velocity, pressure = f['PointData/velocity_01010'], f['PointData/pressure_01010']
            check = {'coordinates_equal_project_mesh':bool(np.array_equal(f['Points/Points'],m['Points/Points'])),
                     'connectivity_exactly_equal_project_mesh':bool(np.array_equal(a,b)),
                     'same_tetrahedra_up_to_local_vertex_order':bool(np.array_equal(np.sort(a,axis=1),np.sort(b,axis=1))),
                     'element_ids_equal':bool(np.array_equal(f['CellData/GlobalElementID'],m['CellData/GlobalElementID'])),
                     'velocity_shape':list(velocity.shape), 'pressure_shape':list(pressure.shape),
                     'all_fields_finite':all(bool(np.isfinite(f[k]).all()) for k in f.files),
                     'velocity_l2_range':[float(np.linalg.norm(velocity,axis=1).min()),float(np.linalg.norm(velocity,axis=1).max())],
                     'pressure_range':[float(pressure.min()),float(pressure.max())]}
            assert check['coordinates_equal_project_mesh'] and check['same_tetrahedra_up_to_local_vertex_order'] and check['element_ids_equal'] and check['all_fields_finite']
            assert check['velocity_shape'] == [153082,3]
            result['cfd_first_step_checks'] = check
    files = [p for p in (PROJECT_ROOT/'data/raw/vascular').rglob('*') if p.is_file() and not p.name.endswith('.partial')]
    result.update({'verified_utc':datetime.now(timezone.utc).isoformat(), 'raw_file_count':len(files),
                   'raw_bytes':sum(p.stat().st_size for p in files),
                   'npz_assets':[{'path':str(p.relative_to(PROJECT_ROOT)), 'bytes':p.stat().st_size, 'sha256':sha256(p)} for p in sorted(DERIVED.glob('*.npz'))]})
    assert result['raw_bytes'] + sum(p.stat().st_size for p in DERIVED.rglob('*') if p.is_file()) < 5_000_000_000
    (DERIVED/'vmr_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index-only',action='store_true')
    parser.add_argument('--allow-incomplete',action='store_true')
    args=parser.parse_args()
    RAW.mkdir(parents=True,exist_ok=True)
    if not args.index_only:
        for path in METADATA: fetch('https://raw.githubusercontent.com/SimVascular/vascularmodel/'+REVISION+'/'+path,RAW/Path(path).name)
        for path in ASSETS: fetch('https://www.vascularmodel.com/'+path, RAW/Path(path).name)
    inventory(require_all=not args.allow_incomplete)
    verify_assets()
if __name__=='__main__':main()
