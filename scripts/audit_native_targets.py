"""Compare native final-minute windows to explicitly parsed upstream targets."""
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.calibration.native_targets import audit_native_targets
if __name__=='__main__':
    result=audit_native_targets(Path(__file__).resolve().parents[1])
    print(json.dumps({r['id']:r['coverage'] for r in result['runs']},indent=2))
