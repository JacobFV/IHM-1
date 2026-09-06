#!/usr/bin/env python3
"""Immutable header-only native regional transporter-graph ownership correction."""
import argparse,json,shlex,shutil
from pathlib import Path
from verify_native_regional_skin import ROOT,RUNTIME,sha
PARENT=RUNTIME/'variants/whole_body_integrity_regional_skin_sweat_v1'
LIBRARY_PIN='bc91cbab829c04bfa2df7490433af5bee752e7514332ac715e9f0bd11a17650a'
MANIFEST_PIN='9f10eeeb03dfd7b4930ac7bb0000b54729fd60cbe9dc65c8d83a2535f8bd9d0a'


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=RUNTIME/'variants/whole_body_integrity_regional_skin_graph_v1');args=parser.parse_args()
    if sha(PARENT/'manifest.json')!=MANIFEST_PIN or sha(PARENT/'libbiogears.so.8.0.0')!=LIBRARY_PIN:raise ValueError('Immutable sweat parent changed')
    parent=json.loads((PARENT/'manifest.json').read_text());objects=shlex.split((PARENT/'objects.rsp').read_text());cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    if len(objects)!=len(set(objects)) or set(objects)!=set(parent['object_sha256']):raise ValueError('Parent object inventory mismatch')
    for obj in objects:
        if sha(cwd/obj)!=parent['object_sha256'][obj]:raise ValueError('Parent object changed: '+obj)
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ['libbiogears.so.8.0.0','objects.rsp','native_signed_muscle_port.h','native_regional_species.h']:shutil.copyfile(PARENT/name,out/name)
    header=ROOT/'scripts/native_regional_skin.h';shutil.copyfile(header,out/header.name)
    manifest={**parent,'variant':out.name,'parent_variant':PARENT.name,'parent_manifest_sha256':MANIFEST_PIN,'parent_library_sha256':LIBRARY_PIN,
        'library_sha256':sha(out/'libbiogears.so.8.0.0'),'regional_skin_header_sha256':sha(header),
        'header_receipts':{**parent.get('header_receipts',{}),header.name:sha(header)},
        'builder_sha256':sha(Path(__file__)),'inherited_object_count':len(objects),'replaced_object':None,
        'compile_command':None,'link_command':None,'resources':{},
        'scope':'Header-only native Skin graph vertex ownership correction; all360native objects and regional sweat correction inherited',
        'production_ready':False}
    if sha(PARENT/'manifest.json')!=MANIFEST_PIN or sha(PARENT/'libbiogears.so.8.0.0')!=LIBRARY_PIN:raise ValueError('Parent changed while copying')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({'variant':out.name,'manifest_sha256':sha(out/'manifest.json'),'library_sha256':sha(out/'libbiogears.so.8.0.0'),'regional_skin_header_sha256':sha(header)}))


if __name__=='__main__':main()
