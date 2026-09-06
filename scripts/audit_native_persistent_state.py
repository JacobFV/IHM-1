#!/usr/bin/env python3
"""Source-only targeted save/load audit. Does not execute or patch native code.

Assertions pin reviewed control-flow facts; numeric witnesses replay only the
small documented branches, not BioGears trajectories. Findings are defects,
not an assertion that serialization is correct or that native replay ran.
"""
from pathlib import Path
import re,json,hashlib,argparse
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/raw/physiology/biogears/projects/biogears/libBiogears'


def read(name):return (BASE/name).read_text(errors='replace')

def body(text,signature):
    start=text.index(signature);begin=text.index('{',start);level=0
    for i in range(begin,len(text)):
        if text[i]=='{':level+=1
        elif text[i]=='}':
            level-=1
            if level==0:return text[begin:i+1]
    raise ValueError('Unclosed function '+signature)

def uncomment(text):
    return re.sub(r'//[^\n]*|/\*.*?\*/','',text,flags=re.S)

def anchor(path,token):
    text=read(path);index=text.index(token)
    return dict(path=str((BASE/path).relative_to(ROOT)),line=text[:index].count('\n')+1,token=token)


def audit():
    paths=dict(engine_io='src/io/biogears/BioGearsPhysiology.cpp',cdm_io='src/io/cdm/Physiology.cpp',
        gi_base='src/cdm/system/physiology/SEGastrointestinalSystem.cpp',
        tissue='src/engine/Systems/Tissue.cpp',gi='src/engine/Systems/Gastrointestinal.cpp',
        renal='src/engine/Systems/Renal.cpp',energy='src/engine/Systems/Energy.cpp',endocrine='src/engine/Systems/Endocrine.cpp',
        full_io='src/io/biogears/BioGears.cpp',system='include/biogears/engine/Controller/BioGearsSystem.h',engine='src/engine/Controller/BioGears.cpp')
    texts={k:read(p) for k,p in paths.items()}
    outer=uncomment(body(texts['engine_io'],'void BiogearsPhysiology::Marshall(const Gastrointestinal&'))
    base=uncomment(body(texts['cdm_io'],'void Physiology::Marshall(const SEGastrointestinalSystem&'))
    incoming=uncomment(body(texts['engine_io'],'void BiogearsPhysiology::UnMarshall(const CDM::BioGearsGastrointestinalSystemData&'))
    cat=uncomment(body(texts['gi'],'void Gastrointestinal::ProcessDrugCAT()'))
    invalid=uncomment(body(texts['gi_base'],'void SEGastrointestinalSystem::Invalidate()'))
    nested=uncomment(body(texts['cdm_io'],'void Physiology::Marshall(const SEDrugTransitState&'))
    assert 'DrugTransitStates' not in outer+base and 'in.DrugTransitStates()' in incoming
    assert all(s in cat for s in ['SetLumenSolidMasses','SetLumenDissolvedMasses','SetEnterocyteMasses','GetTotalMassExcreted().IncrementValue','totalEffluxToPortal_ug_Per_s'])
    assert 'm_DrugTransitStates' not in invalid
    assert 'out.BioGears::SetUp();' in texts['full_io']
    assert 'm_GastrointestinalSystem = Gastrointestinal::make_unique(*this);' in body(texts['engine'],'void BioGears::SetUp()')
    assert 'io::Property::Marshall(*in.m_TotalMassMetabolized, out.MassExcreted());' in nested
    setup=uncomment(body(texts['tissue'],'void Tissue::SetUp()'))
    burn=uncomment(body(texts['tissue'],'void Tissue::CalculateCompartmentalBurn()'))
    tissue_out=uncomment(body(texts['engine_io'],'void BiogearsPhysiology::Marshall(const Tissue&'))
    tissue_in=uncomment(body(texts['engine_io'],'void BiogearsPhysiology::UnMarshall(const CDM::BioGearsTissueSystemData&'))
    regions=['trunk','leftArm','rightArm','leftLeg','rightLeg']
    members=[f'm_{r}DeltaResistance_mmHg_s_Per_mL' for r in regions]+['m_compartmentSyndromeCount','m_baselineECFluidVolume_mL']+[f'm_{r}Escharotomy' for r in regions]
    for name in members:
        assert re.search(re.escape(name)+r'\s*=\s*(?:0\.0|false)',setup)
        assert name in burn and name not in tissue_in+tissue_out
    assert 'CalculateCompartmentalBurn();' in uncomment(body(texts['tissue'],'void Tissue::PreProcess()'))
    assert 'SetUp();' in body(texts['system'],'virtual void LoadState()')
    # A negative control: sodium averaging is saved, and the derived renal
    # permeability modifiers are recomputed before active-transport consumption.
    renalpre=body(texts['renal'],'void Renal::PreProcess()')
    assert 'CalculateReabsorptionFeedback();' in renalpre
    assert 'CalculateOsmoreceptorFeedback();' in body(texts['renal'],'void Renal::CalculateReabsorptionFeedback()')
    assert 'm_LeftReabsorptionPermeabilityModificationFactor =' in body(texts['renal'],'void Renal::CalculateOsmoreceptorFeedback()')
    assert 'in.SodiumConcentration_mg_Per_mL()' in texts['engine_io']
    for name in ['m_packOn','m_previousWeightPack_kg']:
        assert name in body(texts['engine_io'],'void BiogearsPhysiology::Marshall(const Energy&')
        assert name in body(texts['engine_io'],'void BiogearsPhysiology::UnMarshall(const CDM::BioGearsEnergySystemData&')
    # Stable source-only numerical witnesses; quantities chosen to exercise
    # source branches, not claimed as measured or native evolved states.
    old_baseline=10000.;current_ecf=18000.
    def creep(baseline):
        if current_ecf>1.75*baseline:return min(max(3.*current_ecf/baseline,3.5),10.)
        if current_ecf>1.5*baseline:return min(max(1.5*current_ecf/baseline,1.5),1.75)
        return 1.
    witnesses=dict(burn_baseline=dict(before_reload=creep(old_baseline),after_setup_rebaseline=creep(current_ecf)),
        burn_resistance=dict(accumulated_prior=2.,new_increment=.1,continuous_next_increment=2.1,reloaded_next_increment=.1),
        excreted_copy=dict(metabolized_ug=3.,excreted_ug=7.,emitted_excreted_ug=3.),
        fresh_gi_map=dict(source_entries=1,serialized_entries=0,fresh_loaded_entries=0),
        disclaimer='Branch-level source witnesses only; no native engine or serialization roundtrip executed')
    evidence={
        'gi_missing_transit_map':[anchor(paths['engine_io'],'void BiogearsPhysiology::Marshall(const Gastrointestinal&'),anchor(paths['cdm_io'],'void Physiology::Marshall(const SEGastrointestinalSystem&'),anchor(paths['gi'],'void Gastrointestinal::ProcessDrugCAT()')],
        'gi_direct_codec_stale_map_and_cleanup':[anchor(paths['gi_base'],'void SEGastrointestinalSystem::Invalidate()'),anchor(paths['gi_base'],'SEDrugTransitState* SEGastrointestinalSystem::NewDrugTransitState')],
        'gi_wrong_excreted_member':[anchor(paths['cdm_io'],'io::Property::Marshall(*in.m_TotalMassMetabolized, out.MassExcreted());')],
        'tissue_burn_history_reset':[anchor(paths['tissue'],'m_leftArmDeltaResistance_mmHg_s_Per_mL = 0.0;'),anchor(paths['tissue'],'void Tissue::CalculateCompartmentalBurn()'),anchor(paths['engine_io'],'void BiogearsPhysiology::UnMarshall(const CDM::BioGearsTissueSystemData&')]}
    return dict(schema_version=1,kind='reviewed_source_control_flow_audit',native_execution=False,source_defects_confirmed=list(evidence),
        evidence=evidence,burn_history_members=members,branch_witnesses=witnesses,
        negative_controls=['Full engine load recreates GI: stale-map issue applies to direct same-system decode/cleanup, while full engine suffers omitted-state loss.',
            'GI digestion/transit parameter arrays are SetUp constants, stomach contents serialized, decrement flag explicitly enabled after load.',
            'Renal running averages and afferent feedback state serialized; permeability multipliers recomputed in PreProcess before transport.',
            'Energy pack history and blood averages serialized; m_Test only debug probe.',
            'Endocrine private substance pointers/setup constants are caches; AverageBiologicalDebt has assignment but no active read.',
            'Tissue metabolic averages/stored nutrients serialized; hepatic O2/CO2 scratch transfer consumed/reset after Hepatic::Process within full tick; not a settled checkpoint gap.',
            'Tissue totalFatConsumed_g only disabled debug probe; lastFatigueTime unused; not demonstrated physiology history gaps.'],
        source_hashes={str((BASE/p).relative_to(ROOT)):hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in paths.values()})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    text=json.dumps(audit(),indent=2)+'\n'
    if args.output:args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(text)
    else:print(text,end='')
