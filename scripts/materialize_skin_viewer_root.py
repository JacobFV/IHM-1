#!/usr/bin/env python3
"""Byte-budgeted independent geometry snapshot; never mutate shared sources."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from prepare_skin_layer_migration import encode,digest
from verify_skin_candidate_stage import verify
from ihm.assembly.body import CanonicalBody


def sha(path):
 with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def materialize(root,epoch,plan_path,output,budget=193000000):
 root=Path(root).resolve();epoch=Path(epoch).resolve();output=Path(output).resolve()
 if output.exists() or not output.is_relative_to(root/'data/derived') or output.is_relative_to(root/'data/derived/canonical'):raise ValueError('Fresh derived output required')
 verify(epoch);plan_raw=Path(plan_path).read_bytes();plan=json.loads(plan_raw)
 if not 0<plan['total_bytes']<=budget or sum(e['bytes'] for e in plan['files'].values())!=plan['total_bytes']:raise ValueError('Geometry byte budget exceeded')
 staged=epoch/'root';sources={}
 for relative,entry in plan['files'].items():
  path=(root/relative).resolve()
  if not path.is_relative_to(root) or not relative.startswith('data/derived/canonical/geometry/'):raise ValueError('Invalid geometry inventory path')
  stat=path.stat()
  if stat.st_size!=entry['bytes']:raise ValueError('Geometry size changed')
  sources[relative]=(path,stat.st_size,stat.st_mtime_ns,stat.st_ino)
 output.mkdir(parents=True);destination=output/'root';metadata={}
 for path in staged.rglob('*'):
  if not path.is_file():continue
  relative=str(path.relative_to(staged))
  if relative in plan['files']:continue
  target=destination/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
  if sha(path)!=sha(target):raise ValueError('Metadata copy differs')
  metadata[relative]=sha(target)
 copies={}
 for relative,entry in plan['files'].items():
  source,size,mtime,inode=sources[relative];before=sha(source)
  if before!=entry['sha256']:raise ValueError('Geometry source hash changed: '+relative)
  target=destination/relative;target.parent.mkdir(parents=True,exist_ok=True)
  if target.exists():raise ValueError('Existing geometry destination')
  subprocess.run(['cp','--reflink=auto','--preserve=mode,timestamps','--',str(source),str(target)],check=True)
  if source.stat().st_ino==target.stat().st_ino and source.stat().st_dev==target.stat().st_dev:raise ValueError('Shared writable inode forbidden')
  after=source.stat()
  if (after.st_size,after.st_mtime_ns,after.st_ino)!=(size,mtime,inode) or sha(source)!=before:raise ValueError('Source changed during copy')
  if sha(target)!=before:raise ValueError('Geometry destination differs')
  copies[relative]={'sha256':before,'bytes':size,'source_mtime_ns':mtime,'independent_inode':True}
 loaded=CanonicalBody.from_workspace(destination)
 report={'schema':'ihm.skin-viewer-root.v1','root':str(destination.relative_to(root)),
  'parent_epoch':str(epoch.relative_to(root)),'plan_sha256':digest(plan_raw),'geometry_files':copies,
  'geometry_bytes':plan['total_bytes'],'metadata':metadata,'entity_count':len(loaded.entities),
  'body_loader_passed':True,'method':'Sequential reflink-auto independent files; verified source before/after and destination hashes',
  'native_executed':False,'published':False,'source_mtime_preserved':True,
  'scope':'Canonical anatomy geometry and body loader complete; no app/browser/server deployment, native factory, or historical research trajectory rerun',
  'implementation_sha256':digest(Path(__file__).read_bytes())}
 (output/'acceptance.json').write_bytes(encode(report));return report

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--epoch',type=Path,required=True);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 r=materialize(a.root,a.epoch,a.plan,a.output);print(json.dumps({k:r[k] for k in ('root','geometry_bytes','entity_count','body_loader_passed','published')}))
