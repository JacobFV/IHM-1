"""Prepare only: compose exact accepted graph objects and matched regional probe."""
import argparse,json,shlex,shutil,tempfile
from pathlib import Path
from audit_regional_cardiovascular_membership import ROOT,PARENT,audit,sha
from build_regional_signed_native import step_source,ENGINE,ENGINE_PIN
RUNTIME=ROOT/'data/runtime/physiology'
DONORS={'baseline':'whole_body_integrity_cardio_region_io_v2','reader':'whole_body_integrity_cardio_reader_region_v1','prior':'whole_body_integrity_cardio_prior_region_v1'}
def prepare(out):
 membership=audit();assert membership['regional_parent_manifest_sha256']=='52de403e3dab2682f774755b7bf1370711bc5ab218244c581460b84a82aad130'
 pm=json.loads((PARENT/'manifest.json').read_text());assert sha(PARENT/'libbiogears.so.8.0.0')==pm['library_sha256'];assert sha(ENGINE)==ENGINE_PIN
 cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears';original=shlex.split((PARENT/'objects.rsp').read_text());assert len(original)==360
 assert all(sha(cwd/p)==pm['object_sha256'][p] for p in original)
 out=Path(out);out.mkdir(parents=True,exist_ok=False);recipes={}
 for label,name in DONORS.items():
  donor=RUNTIME/'variants'/name;dm=json.loads((donor/'manifest.json').read_text());assert sha(donor/'libbiogears.so.8.0.0')==dm['library_sha256'];assert sha(cwd/'../../../outputs/Release/lib/libbiogears_cdm.so.8.0.0')==dm['paired_cdm_library_sha256']
  objects=original.copy();replaced=[]
  for suffix in ('/io/cdm/Circuit.cpp.o','/Cardiovascular.cpp.o'):
   positions=[i for i,p in enumerate(original) if p.endswith(suffix)];assert len(positions)==1
   if suffix.startswith('/io'):candidate=donor/'Circuit.cpp.o'
   else:
    choices=[p for p in dm['object_sha256'] if p.endswith('/Cardiovascular.cpp.o')];assert len(choices)==1;candidate=cwd/choices[0]
   assert sha(candidate)==dm['object_sha256'][str(candidate)] if str(candidate) in dm['object_sha256'] else sha(candidate)==dm['io_object_sha256']
   objects[positions[0]]=str(candidate);replaced.append({'slot':positions[0],'original':original[positions[0]],'replacement':str(candidate),'sha256':sha(candidate)})
  (out/(label+'.objects.rsp')).write_text(' '.join(shlex.quote(p) for p in objects))
  recipes[label]={'donor_manifest_sha256':sha(donor/'manifest.json'),'replaced':replaced,'unchanged_objects':358,'paired_cdm_library_sha256':dm['paired_cdm_library_sha256']}
 names=('native_regional_coupled_engine.h','native_regional_skin.h','native_regional_species.h','native_body_ports.h','native_tissue_ports.h','native_signed_muscle_port.h','native_signed_vascular_prior.h')
 for name in names:shutil.copyfile(ROOT/'scripts'/name,out/name)
 (out/'native_regional_coupled_step.inc').write_text(step_source(ENGINE.read_text()))
 probe=(ROOT/'scripts/native_signed_cardiovascular_probe.cpp').read_text();probe=probe.replace('"native_coupled_engine.h"','"native_regional_coupled_engine.h"').replace('CoupledBioGearsEngine bg','RegionalCoupledBioGearsEngine bg').replace('bg.coupling_enabled=true;','bg.install_regions();')
 probe=probe.replace('auto values=body_ports(bg);','auto values=body_ports(bg);auto regional_values=bg.regional_ports();values.insert(regional_values.begin(),regional_values.end());')
 old='for(const auto& name:{"Aorta1ToMuscle1","Muscle1ToMuscle2","Muscle2ToVenaCava"})'
 new='std::vector<std::string> observed{"Aorta1ToMuscle1","Muscle1ToMuscle2","Muscle2ToVenaCava","Aorta1ToSkin1","Skin1ToSkin2"};for(size_t region=0;region<3;++region)for(const auto& source:ihm::NativeRegionalSkin::path_names)observed.push_back(ihm::NativeRegionalSkin::name(region,source));\n  for(const auto& name:observed)'
 assert old in probe;probe=probe.replace(old,new);(out/'probe.cpp').write_text(probe)
 record={'schema':'ihm.regional-cardiovascular-preparation.v1','membership':membership,'recipes':recipes,'files':{p.name:sha(p) for p in out.iterdir()},'source_sha256':{str(ROOT/'scripts'/n):sha(ROOT/'scripts'/n) for n in names},'builder_sha256':sha(__file__),'compiled':False,'scope':'Exact graph_v2 composition only; retain regional Diffusion and Energy objects, no sleep/GI/Tissue codec composition implied. Regional install follows native LoadState/SetUp; no second SetUp or region guesses.'}
 (out/'manifest.json').write_text(json.dumps(record,indent=2)+'\n');return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path);a=p.parse_args();target=a.output or Path(tempfile.mkdtemp(prefix='regional-cardiovascular-',dir=ROOT/'data/derived/audits'))/'prepared';print(prepare(target))
