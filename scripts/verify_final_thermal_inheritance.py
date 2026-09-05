"""Run held GI/renal/dry/energy/depletion regressions on the final thermal layer.

Original verifiers and native probes remain unchanged. Two historical comparison
drivers are retained as explicit adapted copies because their old negative-control
variant names and parity requirements do not describe this new thermal layer.
"""
from pathlib import Path
import contextlib,difflib,importlib.util,io,json,sys,tempfile
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from scripts.audit_long_horizon_thermal import sha
from ihm.assembly.systemic_evidence import freeze_sources

VARIANT='whole_body_integrity_evaporation_humidity'
STATE=BASE/'data/derived/audits/thermal-evaporation-humidity-fresh-hour/states/native_stabilized.xml'

def run():
    out=Path(tempfile.mkdtemp(prefix='thermal-final-inheritance-',dir=BASE/'data/derived/audits'))
    library=BASE/'data/runtime/physiology/variants'/VARIANT/'libbiogears.so.8.0.0'
    assert sha(library)=='ab485812968b17e16915b04a6e5144bb1e1ec7f5bf78213936c3334d087c0a5a'
    reports={};originals={}
    print(out,flush=True)
    for name in ('gi_integrity','renal_integrity','dry_gi','energy_integrity','depletion_integrity'):
        original=BASE/'scripts'/('verify_native_'+name+'.py');before=original.read_text();originals[str(original.relative_to(BASE))]=sha(original)
        case=out/name;case.mkdir();copied=case/original.name
        if name in ('gi_integrity','renal_integrity','dry_gi'):
            copied.write_text(before)
            spec=importlib.util.spec_from_file_location('held_'+name,original);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
            module.AUDIT=case
            with (case/'driver_stdout.log').open('w') as log,contextlib.redirect_stdout(log):result=module.verify(VARIANT,STATE)
            passed=result['passed'] if isinstance(result,dict) else result
            report_path=case/VARIANT/'report.json'
            adaptations=['Fresh output root and explicit final library/state arguments through existing API. Probe and assertions unchanged.']
        else:
            after=before
            prefix='energy' if name=='energy_integrity' else 'depletion'
            old=f"out=Path(tempfile.mkdtemp(prefix='{prefix}-integrity-',dir=ROOT/'data/derived/audits'))"
            assert after.count(old)==1
            after=after.replace(old,'out=Path('+repr(str(case))+')')
            if prefix=='energy':
                after=after.replace("for variant in ['whole_body_integrity_gi_water','whole_body_integrity_energy']:","for variant in ["+repr(VARIANT)+"]:")
                after=after.replace("fixed=variant.endswith('_energy')","fixed=True")
                parity=" checks['rest_parity']=results['whole_body_integrity_gi_water']['rest']==results['whole_body_integrity_energy']['rest']"
                assert after.count(parity)==1;after=after.replace(parity,' # Historical cross-variant rest parity is not applicable after thermal law changes.')
                after=after.replace("state=ROOT/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml'",'state=Path('+repr(str(STATE))+')')
                adaptations=['Run final variant only; classify inherited energy correction as fixed.','Use final fresh state for matched rest/exercise/stop.','Omit historical cross-variant rest parity because thermal equations differ; retain all within-variant branch, source, stop, store and causal checks.']
            else:
                after=after.replace("['whole_body_integrity_depletion']",'['+repr(VARIANT)+']').replace("fixed=variant.endswith('_depletion')",'fixed=True')
                adaptations=['Run final variant fixed-only; classify inherited depletion correction as fixed.','Retain original known stomach-water/clock fixture for the 1205-second depletion regression.']
            copied.write_text(after)
            (case/'driver.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile=str(original.relative_to(BASE)),tofile=str(copied.relative_to(BASE)))))
            namespace={'__name__':'retained_adapted_test','__file__':str(original)}
            exec(compile(after,str(copied),'exec'),namespace)
            previous=sys.argv
            try:
                sys.argv=[str(copied)]+(['--fixed-only'] if prefix=='depletion' else [])
                with (case/'driver_stdout.log').open('w') as log,contextlib.redirect_stdout(log):namespace['main']()
            finally:sys.argv=previous
            report_path=case/'report.json';passed=json.loads(report_path.read_text())['passed']
        assert passed,name
        assert sha(original)==originals[str(original.relative_to(BASE))]
        report=json.loads(report_path.read_text())
        sources={str(p.relative_to(BASE)):sha(p) for p in (original,copied,Path(__file__),library,library.parent/'manifest.json',STATE)}
        frozen=freeze_sources(BASE,case,sources)
        reports[name]={'passed':True,'report_path':str(report_path.relative_to(BASE)),'report_sha256':sha(report_path),'adaptations':adaptations,'retained_inputs':frozen,'case_count':len(report.get('rows',report.get('checks',{})))}
        print(json.dumps({'case':name,'passed':True,'report_path':str(report_path.relative_to(BASE))}),flush=True)
    result={'status':'passed','variant':VARIANT,'library_sha256':sha(library),'reports':reports,'original_verifier_sha256':originals,
        'limitations':['Source branch and bounded integration regressions, not six-hour final-variant multisystem acceptance.','Depletion uses its original fixture; energy/renal/GI/dry checks use the final fresh state.']}
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','output_dir':str(out)},indent=2));return out

if __name__=='__main__':run()
