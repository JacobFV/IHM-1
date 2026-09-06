"""Extract held primary human TEP evidence and run small prior scenarios.

Uses retained original CC-BY 4.0 XML/XLSX bytes; no network or native engine.
"""
from pathlib import Path
import sys,json,hashlib,zipfile,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.epithelial_electrodiffusion import EpithelialPatch,condition_tep_magnitude,prior_ensemble,fit_tape_strip_response
DIRECTORY=ROOT/'data/research/skin_epithelial'
SOURCE='https://doi.org/10.1371/journal.pone.0219198'
PREPARATION='human_posterior_forearm_one_adult_female_Abe2019_S4'

def extract():
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(DIRECTORY/'abe-2019-s1.xlsx') as z:
        ss=[''.join(s.itertext()) for s in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        cells={}
        for c in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//m:c',ns):
            v=c.find('m:v',ns)
            if v is not None:cells[c.get('r')]=ss[int(v.text)] if c.get('t')=='s' else v.text
    assert cells['B47']=='S4 Fig' and cells['E47']=='TEP [mV]'
    rows=[]
    for i in range(48,52):
        value=float(cells[f'E{i}'])
        rows.append(dict(tape_strips=int(cells[f'C{i}']),reported_tep_mV=value,tewl_g_m2_h=float(cells[f'D{i}']),
            value_v=abs(value)*1e-3,quantity='tep_magnitude',preparation=PREPARATION,source_id=SOURCE,
            source_cells=f'Sheet1!C{i}:E{i}',sign_status='Magnitude only: negative S4 worksheet values; text gives basal-positive porcine convention, human wiring polarity unresolved'))
    return rows


def main():
    rows=extract()
    raw_paths=['abe-2019.xml','abe-2019-s1.xlsx']
    receipt=dict(reference=SOURCE,authors='Abe et al.',year=2019,license='CC BY 4.0',license_url='https://creativecommons.org/licenses/by/4.0/',accessed_date='2026-09-05',
        acquisition_urls=['https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6620005/fullTextXML','https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0219198.s006'],
        raw_files=[dict(path=str((DIRECTORY/p).relative_to(ROOT)),bytes=(DIRECTORY/p).stat().st_size,sha256=hashlib.sha256((DIRECTORY/p).read_bytes()).hexdigest()) for p in raw_paths],
        participant_count=1,preparation=PREPARATION,rows=rows,
        limitations=['Only four supplied human states (0,1,4,7 strips); protocol mentions 10 but no corresponding worksheet row; missing is not imputed.',
            'Remaining workbook experiments are porcine, excluded from human extraction.',
            'No human current, capacitance, ionic concentration or millisecond transient data in this worksheet.',
            'TEP measurement sign unresolved for human supplement: condition magnitude only, declare model polarity separately.'])
    (DIRECTORY/'source-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    # All numbers below are declared sensitivity priors except initial |TEP|.
    # Independent artificial baths; volumes do NOT allocate additional native Skin.
    scenarios=[]
    for c_scale,g_scale in [(0.5,0.5),(0.5,2),(2,0.5),(2,2)]:
        p=EpithelialPatch(volumes_m3=[1e-12,1e-15,1e-12],
            concentrations_mol_m3=[[140,4,110],[15,140,30],[140,4,110]],
            capacitance_f=np.array([1e-11,1e-11])*c_scale,
            conductance_s=np.array([[1e-11,1e-12,1e-11],[1e-12,1e-10,1e-11],[1e-12,1e-12,1e-12]])*g_scale,
            initial_potential_v=[0,-.03,0],temperature_k=310.15,
            active_flux_mol_s=[[0,0,0],[3e-18,-2e-18,0],[0,0,0]],
            provenance=dict(preparation=PREPARATION,parameters='Unidentified sensitivity priors, including cellular Vm, bath ions, volumes, 37C temperature, conductances, capacitance and powered 3Na/2K transfer',
                capacitance_scale_prior=c_scale,conductance_scale_prior=g_scale,cell_vm_prior_v=-.03,
                transfer_scope='Source-conditioned alternative preparation, not anatomical whole-skin materialization'))
        p=condition_tep_magnitude(p,rows[0],basal_positive=True)
        scenarios.append(p)
    run=prior_ensemble(scenarios,.5,21)
    # Retain explicit coefficient arrays per scenario; provenance alone is not enough.
    for p,r in zip(scenarios,run['members']):
        r['parameters']={name:getattr(p,name).tolist() for name in ['volumes_m3','concentrations_mol_m3','capacitance_f','conductance_s','initial_potential_v','active_flux_mol_s']}
        r['parameters']['temperature_k']=p.temperature_k
    out=dict(schema_version=1,id='source_conditioned_finite_epithelial_demo',source_receipt_sha256=hashlib.sha256((DIRECTORY/'source-receipt.json').read_bytes()).hexdigest(),
        descriptive_voltage_fit=fit_tape_strip_response(rows),ensemble=run,
        validation_scope='Physics audits and within-participant held strip state only; no calibrated RC or human trajectory prediction',
        runtime_sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['ihm/assembly/epithelial_electrodiffusion.py','scripts/build_epithelial_electrodiffusion.py']})
    (DIRECTORY/'source-conditioned-demo.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(rows=len(rows),participants=1,fit=out['descriptive_voltage_fit'],scenario_audits=[r['audit'] for r in run['members']]),indent=2))
if __name__=='__main__':main()
