"""Build the recorded canonical forearm electrical experiment."""
import json
from pathlib import Path
from ihm.assembly.skin_bioelectric import build_skin_electric


if __name__ == '__main__':
    root = Path(__file__).resolve().parents[1]
    artifact = build_skin_electric(root)
    output = root/'data/derived/canonical/skin-electric.json'
    output.write_text(json.dumps(artifact, separators=(',', ':'), allow_nan=False)+'\n')
    print(output)
