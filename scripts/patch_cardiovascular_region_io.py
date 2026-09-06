"""Isolated existing-field IO repair and exact legacy topology metadata migration."""
from pathlib import Path
import argparse,difflib,hashlib,json,re,tempfile,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
DONOR=ROOT/'data/raw/physiology/biogears'
IO=DONOR/'projects/biogears/libBiogears/src/io/cdm/Circuit.cpp'
CONTROLLER=DONOR/'projects/biogears/libBiogears/src/engine/Controller/BioGears.cpp'
IO_SHA='efbff735f5a2931e23c245b89234870acf2e9dd0511a00a9363664bf73da2eb3'
CONTROLLER_SHA='451ddda942a26e547a0ed4f1d5ba875fa7ca9e2c6275895b791f9c81e7a39166'
REGIONS=('Cerebral','Extrasplanchnic','Muscle','Splanchnic','Myocardium')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def patch_io(text):
 if hashlib.sha256(text.encode()).hexdigest()!=IO_SHA:raise ValueError('Changed pinned circuit IO source')
 reader='    io::Property::UnMarshall(in.ResistanceBaseline(), out.GetResistanceBaseline());'
 # The anchor occurs in electrical/fluid/thermal methods: isolate the fluid one.
 start=text.index('  void Circuit::UnMarshall(const CDM::FluidCircuitPathData&');end=text.index('  void Circuit::Marshall(const SEFluidCircuitPath&',start)
 block=text[start:end];assert block.count(reader)==1
 read='\n    out.InvalidateCardiovascularRegion();\n    if (in.CardiovascularRegion().present()) {\n      switch (in.CardiovascularRegion().get()) {\n'
 for region in REGIONS:read+=f'      case CDM::enumResistancePathType::{region}: out.SetCardiovascularRegion(SEResistancePathType::{region}); break;\n'
 read+='      default: throw std::runtime_error("Unknown serialized cardiovascular region");\n      }\n    }'
 text=text[:start]+block.replace(reader,reader+read)+text[end:]
 anchor='    Marshall(static_cast<const SECircuitPath<FLUID_CIRCUIT_PATH>&>(in), static_cast<CDM::CircuitPathData&>(out));'
 assert text.count(anchor)==1
 write='\n    out.CardiovascularRegion().reset();\n    if (in.HasCardiovascularRegion()) {\n      switch (in.GetCardiovascularRegion()) {\n'
 for region in REGIONS:write+=f'      case SEResistancePathType::{region}: out.CardiovascularRegion(CDM::enumResistancePathType(CDM::enumResistancePathType::{region})); break;\n'
 write+='      default: throw std::runtime_error("Unknown native cardiovascular region");\n      }\n    }'
 return '#include <stdexcept>\n'+text.replace(anchor,anchor+write)
def topology():
 text=CONTROLLER.read_text()
 if sha(CONTROLLER)!=CONTROLLER_SHA:raise ValueError('Changed pinned topology constructor')
 variables={};assigned={};deleted=[];assignments=0
 for number,line in enumerate(text.splitlines(),1):
  created=re.search(r'SEFluidCircuitPath& (\w+) = .*CreatePath\(.*BGE::\w+Path::(\w+)\);',line)
  if created:variables[created[1]]=created[2];assigned.pop(created[2],None)
  tagged=re.search(r'(\w+)\.SetCardiovascularRegion\(SEResistancePathType::(\w+)\);',line)
  if tagged:
   if tagged[1] not in variables:raise ValueError('Unresolved native path variable')
   assigned[variables[tagged[1]]]={'region':tagged[2],'assignment_line':number};assignments+=1
  removed=re.search(r'DeleteFluidPath\(BGE::\w+Path::(\w+)\);',line)
  if removed and removed[1] in assigned:deleted.append({'path':removed[1],'delete_line':number,'prior':assigned.pop(removed[1])})
 if assignments!=41 or len(assigned)!=35 or len(deleted)!=6:raise ValueError('Unexpected cardiovascular/cerebral/renal topology inventory')
 return assigned,deleted

def migrate(raw,expected_sha):
 if hashlib.sha256(raw).hexdigest()!=expected_sha:raise ValueError('Legacy snapshot differs from explicit retained identity')
 if b'<CardiovascularRegion>' in raw:raise ValueError('Migration accepts unlabeled legacy topology only')
 assigned,deleted=topology();text=raw.decode();tree=ET.fromstring(raw);paths={};circuits={}
 def local(e):return e.tag.split('}')[-1]
 for e in tree.iter():
  if local(e)=='FluidPath':
   values={local(x):x.text for x in e};name=values['Name']
   if name in paths:raise ValueError('Duplicate fluid path')
   paths[name]=values
  elif local(e)=='FluidCircuit':
   name=next(x.text for x in e if local(x)=='Name');circuits[name]=[x.text for x in e if local(x)=='Path']
 if len(paths)!=348:raise ValueError('Unsupported legacy fluid graph size')
 if not set(assigned)<=set(paths):raise ValueError('Source-assigned path missing from legacy graph')
 absent={'Aorta1ToBrain1','Brain1ToBrain2','LeftKidney1ToLeftKidney2','RightKidney1ToRightKidney2'}
 if absent&set(paths):raise ValueError('Expected cerebral/renal replacement topology absent')
 # Renal construction deletes then recreates these names without region labels.
 for name in ('Aorta1ToLeftKidney1','Aorta1ToRightKidney1'):
  if name not in paths or name in assigned:raise ValueError('Renal reconstructed path identity mismatch')
 graph={'paths':{name:{k:v for k,v in row.items() if k in ('SourceNode','TargetNode')} for name,row in paths.items()},'circuits':circuits}
 inserted=[];seen=set()
 def replace(match):
  block=match[0];name=re.search(r'<Name>([^<]+)</Name>',block)[1]
  if name not in assigned:return block
  if name in seen:raise ValueError('Duplicate migrated path')
  seen.add(name);anchor=re.search(r'<ResistanceBaseline\b[^>]*/>',block)
  if not anchor:raise ValueError('Assigned region has no native resistance baseline')
  tag='\n      <CardiovascularRegion>'+assigned[name]['region']+'</CardiovascularRegion>';inserted.append(tag)
  return block[:anchor.end()]+tag+block[anchor.end():]
 migrated=re.sub(r'<FluidPath>.*?</FluidPath>',replace,text,flags=re.S)
 if seen!=set(assigned):raise ValueError('Incomplete region migration')
 stripped=migrated
 for tag in inserted:stripped=stripped.replace(tag,'',1)
 if stripped!=text:raise ValueError('Migration changed nonmetadata bytes')
 ET.fromstring(migrated)
 return migrated.encode(),{'assigned_paths':assigned,'source_deleted_region_assignments':deleted,'graph_sha256':hashlib.sha256(json.dumps(graph,sort_keys=True).encode()).hexdigest(),'all_nonmetadata_bytes_preserved':True}

def prepare(output=None):
 out=Path(output) if output else Path(tempfile.mkdtemp(prefix='cardiovascular-region-io-',dir=ROOT/'data/derived/audits'));out.mkdir(parents=True,exist_ok=True)
 if any(out.iterdir()):raise ValueError('Fresh empty output required')
 parent=ROOT/'data/runtime/physiology/variants/whole_body_integrity_gi_absorption';pm=json.loads((parent/'manifest.json').read_text())
 if sha(parent/'libbiogears.so.8.0.0')!=pm['library_sha256']:raise ValueError('Changed retained parent library')
 source=IO.read_text();patched=patch_io(source);(out/'Circuit.cpp').write_text(patched)
 for name in ('Circuit.h','Property.h'):(out/name).write_bytes((IO.parent/name).read_bytes())
 (out/'Circuit.patch').write_text(''.join(difflib.unified_diff(source.splitlines(True),patched.splitlines(True),fromfile=str(IO),tofile=str(out/'Circuit.cpp'))))
 reference=json.loads((ROOT/'data/derived/systemic/exertion_v3/exercise/native/manifest.json').read_text());state=Path(reference['configuration']['state_path']);raw=state.read_bytes();migrated,receipt=migrate(raw,reference['state_sha256']);(out/'legacy.xml').write_bytes(raw);(out/'region_metadata.xml').write_bytes(migrated)
 receipt.update(schema='ihm.cardiovascular-region-io-overlay.v1',parent_library_sha256=pm['library_sha256'],parent_manifest_sha256=sha(parent/'manifest.json'),legacy_state_sha256=reference['state_sha256'],migrated_state_sha256=sha(out/'region_metadata.xml'),io_source_sha256=IO_SHA,constructor_source_sha256=CONTROLLER_SHA,patched_io_sha256=sha(out/'Circuit.cpp'),files={p.name:sha(p) for p in out.iterdir()},scope='Explicit structural metadata reconstruction for the exact retained 348-path graph and pinned construction order; no physiological history recovered, no dynamic scalar changed, no original state overwritten')
 (out/'manifest.json').write_text(json.dumps(receipt,indent=2)+'\n');return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output');a=p.parse_args();print(prepare(a.output))
