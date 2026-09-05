"""Independent retained audit of the unmodified native counterion probes."""
from pathlib import Path
import argparse,csv,json,sys,tempfile
import numpy as np
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from scripts.audit_long_horizon_thermal import sha
from ihm.assembly.systemic_evidence import freeze_sources,resolve_sources
from ihm.assembly.native_environment_evidence import _resolve_native_archive,materialize_native_resources

def audit(directory):
    directory=Path(directory).resolve()
    out=Path(tempfile.mkdtemp(prefix='meal-counterion-audit-',dir=BASE/'data/derived/audits'))
    modes=('rest','water_only','na_only','nutrients_only','na_water')
    records={};sources={str(Path(__file__).relative_to(BASE)):sha(Path(__file__))};results={}
    expected_time=np.r_[0.,.02,np.arange(10.,601.,10.)]
    for mode in modes:
        work=directory/mode;manifest=json.loads((work/'manifest.json').read_text())
        execution=json.loads((work/'execution.json').read_text())
        assert execution['returncode']==0 and manifest['configuration']['mode']==mode
        sources.update(resolve_sources(BASE,work,manifest['source_hashes']))
        # This direct C++ diagnostic predates NativeSession's execution-selection
        # receipt. Verify its archive and detached tree without manufacturing one.
        sources.update(_resolve_native_archive(BASE,work))
        tree=materialize_native_resources(BASE,work)
        assert (work/'substances').resolve()==tree/'substances'
        paths=[work/name for name in ('electrolytes.csv','states/final_state.xml','manifest.json','execution.json','receipts.jsonl')]
        sources.update({str(p.relative_to(BASE)):sha(p) for p in paths})
        rows=list(csv.DictReader((work/'electrolytes.csv').open()))
        data={key:np.array([float(row[key]) for row in rows]) for key in rows[0]}
        assert np.array_equal(data['time_s'],expected_time)
        for key,values in data.items():
            invalid=np.flatnonzero(~np.isfinite(values))
            assert len(invalid)==0 or (key=='respiratory_cycle_fraction' and invalid.tolist()==[0])
            if key.endswith('.mass_g'):assert np.all(values>=0)
        assert np.max(np.abs(data['native_sid_mmol_l']-data['vc_sid_reconstructed_mmol_l']))<1e-10
        assert np.all(data['SmallIntestineChyme.Chloride.mass_g']==0)
        records[mode]=data
        add_na=1. if mode in ('na_only','na_water') else 0.
        gut_initial=data['stomach_sodium_g'][0]+data['SmallIntestineChyme.Sodium.mass_g'][0]
        net_gi_na= gut_initial+add_na-data['stomach_sodium_g'][-1]-data['SmallIntestineChyme.Sodium.mass_g'][-1]
        keys=['arterial_ph','arterial_co2_mmhg','native_sid_mmol_l','vc_sid_reconstructed_mmol_l',
              'VenaCava.Sodium.mmol_l','VenaCava.Chloride.mmol_l','VenaCava.Potassium.mmol_l','VenaCava.Lactate.mmol_l',
              'Aorta.Bicarbonate.mmol_l','stomach_water_ml','stomach_sodium_g',
              'SmallIntestineChyme.Sodium.mass_g','SmallIntestineChyme.Chloride.mass_g']
        results[mode]={'samples':len(rows),'seconds':600,'net_gi_to_vascular_na_g':float(net_gi_na),
            'net_gi_to_vascular_na_mmol':float(net_gi_na*1000/22.9898),
            'observations':{key:dict(initial=float(data[key][0]),final=float(data[key][-1]),minimum=float(data[key].min()),maximum=float(data[key].max())) for key in keys}}
    for mode,data in records.items():
        for key,values in data.items():
            assert np.allclose(values[0],records['rest'][key][0],rtol=0,atol=0,equal_nan=True)
        for key,entry in results[mode]['observations'].items():
            entry['final_minus_rest']=float(data[key][-1]-records['rest'][key][-1])
    # Diagnostic contrasts, not thresholds for biological acceptance.
    assert results['nutrients_only']['net_gi_to_vascular_na_g']>20*results['rest']['net_gi_to_vascular_na_g']
    assert results['nutrients_only']['observations']['arterial_ph']['final_minus_rest']>.03
    assert abs(results['na_water']['observations']['native_sid_mmol_l']['final_minus_rest'])<1e-5
    donor=BASE/'data/raw/physiology/biogears/projects/biogears/libBiogears'
    source_files=[BASE/'data/runtime/physiology/variants/whole_body_integrity_depletion/Gastrointestinal.cpp',
        donor/'src/engine/Systems/BloodChemistry.cpp',donor/'src/engine/Systems/Saturation.cpp',
        donor/'include/biogears/cdm/patient/SENutrition.h',
        BASE/'data/runtime/physiology/biogears-build/runtime/substances/Sodium.xml',
        BASE/'data/runtime/physiology/biogears-build/runtime/substances/Chloride.xml',
        BASE/'data/runtime/physiology/biogears-build/runtime/substances/Calcium.xml']
    sources.update({str(p.relative_to(BASE)):sha(p) for p in source_files})
    freeze_sources(BASE,out,sources)
    result=dict(passed=True,kind='counterion_omission_diagnostic_not_biological_validation',
        experiment=str(directory.relative_to(BASE)),source_hashes=sources,results=results,
        sodium_molar_mass_g_per_mol=22.9898,one_g_na_mmol=1000/22.9898,
        initial_state_has_unpaired_nutrition_sodium_g=1.,initial_state_has_nutrition_water_ml=500.,
        interpretation='Net GI sodium is the difference of stomach+chyme stores after declared addition, including secretion minus reabsorption; not gross cotransport or complete whole-body charge balance.',
        limits=['No chloride input supported by SENutrition; no counterion supplied or synthesized in this probe.',
                'Chyme chloride remains zero; vascular chloride changes through other transport/dilution.',
                'Charged sodium transport is paired in mass between GI and vascular owners; acid-base electroneutrality is not audited by that mass identity.',
                'Custom adapter has frozen pre-start archives and detached resource hashes, but not NativeSession execution-inputs.json; actual loader access was not traced.',
                'Ten-minute perturbations identify a pathway defect, not long-horizon electrolyte calibration.',
                'Source pH includes carbon dioxide, bicarbonate, buffers and other ion terms; SID change does not isolate every pH contributor.'])
    (out/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(out/'audit.json');return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);args=parser.parse_args();audit(args.directory)
