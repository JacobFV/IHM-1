"""Native compliance port contract and optional actual circuit perturbation."""
from pathlib import Path
import argparse
import json
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native import NativeConfig,run_native,load_trajectory

class RespiratoryPortTests(unittest.TestCase):
    def test_explicit_compliance_configuration(self):
        c=NativeConfig.from_dict({'chest_compliance_l_cmH2O':.2})
        self.assertEqual(c.chest_compliance_l_cmH2O,.2)
        for value in [0,-1,True,float('nan'),10]:
            with self.assertRaises(ValueError):NativeConfig.from_dict({'chest_compliance_l_cmH2O':value})

def experiment():
    from ihm.assembly.body import write_json
    import tempfile
    base=ROOT/'data/derived/audits/respiratory-port'
    base.mkdir(parents=True,exist_ok=True)
    runroot=Path(tempfile.mkdtemp(prefix='experiment-',dir=base))
    records={}
    for label,compliance in [('compliant',.2),('stiff',.1)]:
        directory=runroot/label
        summary=run_native(NativeConfig(seconds=12,patient='IHMGenericMale',engine_variant='saturation_bounds_heatflux',chest_compliance_l_cmH2O=compliance),directory)
        trajectory=load_trajectory(directory)
        port=load_trajectory(directory/'respiratory_port.csv')
        total=np.array(port['values']['LeftChestCompliance(L/cmH2O)'])+np.array(port['values']['RightChestCompliance(L/cmH2O)'])
        np.testing.assert_allclose(total[1:],compliance,rtol=1e-10)
        lung=np.array(trajectory['values']['TotalLungVolume(mL)'])
        records[label]={'directory':str(directory.relative_to(ROOT)),'compliance_l_cmH2O':compliance,'volume_excursion_ml':float(np.ptp(lung[200:])),'flow_excursion_l_s':float(np.ptp(port['values']['AirwayFlow(L/s)'][200:])), 'native_csv_sha256':summary['csv_sha256'],'port_receipt':summary['respiratory_port']}
    assert abs(records['compliant']['volume_excursion_ml']-records['stiff']['volume_excursion_ml'])>1,records
    result={'records':records,'verified':'Changing actual chest compliance changes native pressure-flow ventilation. Gas and blood remain natively integrated.','limitations':['Prescribed constitutive parameter, not full 3D pressure/work feedback.','No empirical calibration established.']}
    write_json(runroot/'verification.json',result);write_json(base/'latest.json',result);print(json.dumps(result,indent=2))

def resume_experiment():
    from ihm.native import summarize
    import tempfile
    state=ROOT/'data/derived/audits/respiratory-port/experiment-xl1lr93f/compliant/states/native_final.xml'
    directory=Path(tempfile.mkdtemp(prefix='resume-',dir=ROOT/'data/derived/audits/respiratory-port'))
    config=NativeConfig(seconds=.02,patient='IHMGenericMale',state_path=str(state),engine_variant='saturation_bounds_heatflux')
    summary=run_native(config,directory)
    port=load_trajectory(directory/'respiratory_port.csv')
    assert abs(port['time_s'][0]-12.)<1e-8 and abs(port['time_s'][-1]-12.02)<1e-8
    assert summarize(directory)['respiratory_port']['csv_sha256']==summary['respiratory_port']['csv_sha256']
    summarize(directory,config)  # Retain the full continuation configuration.
    print('Saved native state continuation clock PASS:',directory)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--experiment',action='store_true');p.add_argument('--resume',action='store_true');a=p.parse_args()
    result=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(RespiratoryPortTests))
    if not result.wasSuccessful():raise SystemExit(1)
    if a.experiment:experiment()
    if a.resume:resume_experiment()
