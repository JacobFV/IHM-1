"""Fit the source-audited human wound-field phenotype; retain holdout and exclusions."""
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.calibration.skin import build_skin_calibration
if __name__=='__main__':
    result=build_skin_calibration(Path(__file__).resolve().parents[1])
    print(json.dumps(dict(split=result['split'],metrics=result['metrics'],parameters=dict(zip(result['feature_names'],result['fit']['parameters']))),indent=2))
