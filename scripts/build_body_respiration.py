#!/usr/bin/env python3
"""Build the source-bound reduced thoracic model separately from body assembly."""
from pathlib import Path
import json
import sys

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
from ihm.assembly.respiration import build
from ihm.assembly.body import write_json

if __name__ == '__main__':
    payload = build(BASE)
    write_json(BASE/'data/derived/canonical/respiration.json', payload)
    print(json.dumps({'model_id': payload['model_id'], 'bound_entities': len(payload['bindings']), 'calibrated': False}))
