"""Actual native path query must preserve state and activation trajectories."""
from pathlib import Path
import json
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream


def main():
    output=Path(tempfile.mkdtemp(prefix='native-moment-arms-',dir=ROOT/'data/derived'))
    plant=NativeMechanicalStream(ROOT,output/'native',environment='free',target_mass_kg=70)
    try:
        initial=plant.snapshot()
        checkpoint=plant.checkpoint()
        baseline=plant.advance(.002)
        plant.restore(checkpoint)
        result=plant.moment_arms(muscles=['soleus_r','tibant_r'],coordinates=['ankle_angle_r'])
        arms=result['moment_arms_m']
        assert arms['soleus_r']['ankle_angle_r']*arms['tibant_r']['ankle_angle_r']<0, 'Antagonists need opposite ankle moments'
        assert result['time_s']==initial['time_s'] and plant.snapshot()['muscles']==initial['muscles']
        repeat=plant.advance(.002)
        assert repeat['bodies']==baseline['bodies'] and repeat['muscles']==baseline['muscles'], 'Query changed trajectory'
        try:plant.moment_arms(muscles=['soleus_r'],coordinates=['pelvis_ty'])
        except ValueError:pass
        else:raise AssertionError('Translational coordinate accepted as metre moment arm')
        plant.release(checkpoint)
        report={'passed':True,'moment_arms':result,'scope':'Read-only source path derivatives and exact native branch equivalence'}
        (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report|{'output':str(output)},indent=2))
    finally:plant.close()


if __name__=='__main__':main()
