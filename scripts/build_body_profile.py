from pathlib import Path
import json
from ihm.assembly.profile import build_profile
print(json.dumps(build_profile(Path(__file__).resolve().parents[1]),indent=2))
