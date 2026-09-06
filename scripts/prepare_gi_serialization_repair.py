#!/usr/bin/env python3
"""Extract source-pinned candidate methods and exact isolated upstream patch."""
from pathlib import Path
import difflib,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/raw/physiology/biogears/projects/biogears/libBiogears'
FILES=['src/io/biogears/BioGearsPhysiology.cpp','src/io/cdm/Physiology.cpp','src/cdm/system/physiology/SEGastrointestinalSystem.cpp']

def method(text,signature):
    start=text.index(signature);brace=text.index('{',start);depth=1;end=brace+1
    while depth:
        depth+=(text[end]=='{')-(text[end]=='}');end+=1
    return text[start:end]

def prepare(out, full_sources=True):
    original={name:(BASE/name).read_text() for name in FILES};fixed=dict(original);selected=[]
    name=FILES[0];sig='void BiogearsPhysiology::Marshall(const Gastrointestinal& in, CDM::BioGearsGastrointestinalSystemData& out)';old=method(fixed[name],sig)
    new=old[:-1]+'''  out.DrugTransitStates().clear();
    for (const auto& item : in.m_DrugTransitStates) {
      if (!item.second) continue; // A null lookup entry owns no recorded state.
      auto state = std::make_unique<CDM::DrugTransitStateData>();
      io::Physiology::Marshall(*item.second, *state);
      out.DrugTransitStates().push_back(std::move(state));
    }
  }''';fixed[name]=fixed[name].replace(old,new);selected.append(('io',new))
    name=FILES[1];sig='void Physiology::Marshall(const SEDrugTransitState& in, CDM::DrugTransitStateData& out)';old=method(fixed[name],sig)
    new=old.replace('std::unique_ptr<CDM::ScalarMassData>()','std::make_unique<CDM::ScalarMassData>()').replace('io::Property::Marshall(*in.m_TotalMassMetabolized, out.MassExcreted())','io::Property::Marshall(*in.m_TotalMassExcreted, out.MassExcreted())')
    new=new.replace('{\n','{\n    out.LumenDissolvedMasses().clear(); out.LumenSolidMasses().clear(); out.EnterocyteMasses().clear();\n',1)
    fixed[name]=fixed[name].replace(old,new);selected.append(('io',new))
    name=FILES[2];sig='void SEGastrointestinalSystem::Invalidate()';old=method(fixed[name],sig)
    new=old.replace('  SESystem::Invalidate();','  for (auto& item : m_DrugTransitStates) delete item.second;\n  m_DrugTransitStates.clear();\n  SESystem::Invalidate();')
    fixed[name]=fixed[name].replace(old,new);selected.append(('',new))
    sig='SEDrugTransitState* SEGastrointestinalSystem::NewDrugTransitState(const SESubstance* sub)';old=method(fixed[name],sig)
    new=sig+"\n{\n  if (!sub) throw CommonDataModelException(\"Missing drug substance\");\n  auto fresh = std::make_unique<SEDrugTransitState>(*sub);\n  auto& owned = m_DrugTransitStates[sub];\n  delete owned;\n  owned = fresh.release();\n  return owned;\n}"
    fixed[name]=fixed[name].replace('#include <biogears/cdm/system/physiology/SEGastrointestinalSystem.h>', '#include <memory>\n#include <biogears/cdm/system/physiology/SEGastrointestinalSystem.h>')
    fixed[name]=fixed[name].replace(old,new);selected.append(('',new))
    out.mkdir(parents=True,exist_ok=True)
    if full_sources:
        for name,body in fixed.items():
            path=out/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(body)
    patch=''.join(''.join(difflib.unified_diff(original[n].splitlines(True),fixed[n].splitlines(True),fromfile='a/'+n,tofile='b/'+n)) for n in FILES)
    patch=''.join('\n' if line==' \n' else line for line in patch.splitlines(True))
    (out/'repair.patch').write_text(patch)
    (out/'source_hashes.json').write_text(json.dumps({n:hashlib.sha256((BASE/n).read_bytes()).hexdigest() for n in FILES},indent=2)+'\n')
    return out/'repair.patch'
if __name__=='__main__':prepare(ROOT/'data/research/gi_serialization', full_sources=False)
