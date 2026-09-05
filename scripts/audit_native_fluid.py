from pathlib import Path
import json,argparse
from ihm.native.budgets import audit_fluid_budget
p=argparse.ArgumentParser();p.add_argument('run_dir',type=Path);a=p.parse_args()
r=audit_fluid_budget(a.run_dir);(a.run_dir/'fluid-budget.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:r[k] for k in ['circuit_volume_change_m3','stomach_water_change_m3','combined_volume_change_m3','relative_combined_change']},indent=2))
