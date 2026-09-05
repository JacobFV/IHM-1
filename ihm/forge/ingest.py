"""normalized physical measurements; raw modalities require explicit operators."""
from dataclasses import dataclass
import csv
import json
import math
import hashlib
from pathlib import Path
import numpy as np


@dataclass(frozen=True)
class SourceCard:
    id: str
    kind: str
    components: tuple[str, ...]
    reference: str
    species: str = 'human'
    status: str = 'local'

    def __post_init__(self):
        if not self.id or not self.reference or self.kind not in ('measured', 'derived', 'synthetic'):
            raise ValueError('source requires id, reference and explicit evidence kind')
        if not self.components:
            raise ValueError('source must declare component bindings')

    @classmethod
    def read(cls, path):
        data = json.loads(Path(path).read_text())
        data['components'] = tuple(data['components'])
        return cls(**data)


@dataclass(frozen=True)
class Evidence:
    source: str
    record: str
    subject: str
    time: float
    component: str
    value: float
    variance: float
    kind: str
    reference: str = 'unspecified: caller-supplied evidence'
    species: str = 'human'
    bindings: tuple[str, ...] = ()
    data_sha256: str = ''

    def __post_init__(self):
        if not all((self.source, self.record, self.subject, self.component)):
            raise ValueError('evidence needs source, record, subject and component')
        if not all(math.isfinite(v) for v in (self.time, self.value, self.variance)) or self.time < 0 or self.variance <= 0:
            raise ValueError('evidence requires finite value/time and positive variance')
        if self.species != 'human':
            raise ValueError('nonhuman measurements require an explicit transfer model')
        if self.kind not in ('measured', 'derived', 'synthetic'):
            raise ValueError('unknown evidence kind')

    @property
    def key(self):
        return (self.source, self.record, self.component)


def convert(value, variance, unit, target, component):
    scale, offset = 1., 0.
    if unit != target:
        conversions = {('kPa', 'mmHg'): (7.50061683, 0),
                       ('%', '1'): (.01, 0), ('mV', 'V'): (.001, 0),
                       ('uV', 'V'): (1e-6, 0), ('degF', 'degC'): (5/9, -32*5/9),
                       ('mL', 'L'): (.001, 0), ('mM', 'mmol/L'): (1, 0)}
        if component in ('metabolic.glucose', 'device.glucose') and (unit, target) == ('mg/dL', 'mmol/L'):
            scale = 1 / 18.01559
        elif (unit, target) in conversions:
            scale, offset = conversions[unit, target]
        else:
            raise ValueError(f'no declared conversion from {unit} to {target} for {component}')
    return float(value)*scale + offset, float(variance)*scale**2


def load(path, card, registry):
    """CSV/JSON long records or NPZ columns with the same explicit schema.

    NPZ stores one row per time and component; no implicit resampling or waveform
    reduction. Variance is supplied by the source, never inferred from amplitude.
    """
    if card.status != 'local' or card.species != 'human':
        raise ValueError('catalog-only or nonhuman data cannot condition a human subject')
    if any(c not in registry.components for c in card.components):
        raise ValueError('source card has unknown component binding')
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    if path.suffix == '.csv':
        with path.open(newline='') as f:
            rows = list(csv.DictReader(f))
    elif path.suffix == '.json':
        rows = json.loads(path.read_text())
        if not isinstance(rows, list):
            raise ValueError('JSON measurements must be a list')
    elif path.suffix == '.npz':
        with np.load(path, allow_pickle=False) as data:
            arrays = {k: data[k] for k in ('id', 'subject', 'time', 'component', 'value', 'unit', 'variance')}
            if any(a.ndim != 1 for a in arrays.values()) or len({len(a) for a in arrays.values()}) != 1:
                raise ValueError('NPZ columns must be equally sized vectors')
            rows = [dict(zip(arrays, values)) for values in zip(*arrays.values())]
    else:
        raise ValueError('supported data formats: .csv, .json, .npz')
    result, seen = [], set()
    for row in rows:
        if row.get('value') is None or str(row.get('value')).strip() == '':
            continue
        component = str(row['component'])
        if component not in card.components:
            raise ValueError(f'undeclared source binding: {component}')
        value, variance = convert(row['value'], row['variance'], str(row['unit']),
                                  registry.components[component].unit, component)
        e = Evidence(card.id, str(row['id']), str(row['subject']), float(row['time']),
                     component, value, variance, card.kind, card.reference, card.species, card.components, digest)
        if e.key in seen:
            raise ValueError(f'duplicate record: {e.key}')
        seen.add(e.key); result.append(e)
    return sorted(result, key=lambda e: e.time)
