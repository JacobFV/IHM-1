"""Measured joint targets and cross-model CMC seeds, never prescribed motion.

The Rajagopal coordinates and Gait2392 excitations are different subjects/trials.
Neither recording certifies native free-contact gait or periodic endpoints.
Only independent joint angles and anatomically mapped muscle excitations leave this module;
pelvis coordinates, residual actuators, and prescribed GRFs are excluded.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import math
import xml.etree.ElementTree as ET
import numpy as np

_COORDS = 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/coordinates.sto'
_CMC = 'data/raw/anatomy/opensim-models/source/Pipelines/Gait2392_Simbody/OutputReference/ResultsCMC/subject01_walk1_controls.sto'
_MODEL = 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled.osim'
_ALIASES = dict(glut_med='glmed', glut_min='glmin', glut_max='glmax', bifemlh='bflh',
    bifemsh='bfsh', sar='sart', rect_fem='recfem', vas_int='vasint', vas_lat='vaslat',
    vas_med='vasmed', med_gas='gasmed', lat_gas='gaslat', tib_ant='tibant',
    tib_post='tibpost', flex_dig='fdl', flex_hal='fhl', ext_dig='edl',
    ext_hal='ehl', per_brev='perbrev', per_long='perlong', add_long='addlong',
    add_brev='addbrev', peri='piri')


def _read(path):
    lines = path.read_text().splitlines()
    end = next(i for i, line in enumerate(lines) if line.strip() == 'endheader')
    metadata = dict(line.split('=', 1) for line in lines[:end] if '=' in line)
    if metadata.get('inDegrees') != 'no':
        raise ValueError(f'Expected explicit inDegrees=no: {path}')
    names = lines[end + 1].split()
    data = np.array([list(map(float, line.split())) for line in lines[end + 2:] if line.strip()])
    if names[0] != 'time' or len(set(names)) != len(names) or data.ndim != 2 or data.shape[1] != len(names):
        raise ValueError(f'Malformed storage columns: {path}')
    duplicate = np.r_[False, np.diff(data[:, 0]) == 0]
    if any(not np.array_equal(data[i], data[i - 1]) for i in np.flatnonzero(duplicate)):
        raise ValueError(f'Conflicting duplicate source times: {path}')
    duplicate_count = int(duplicate.sum())
    data = data[~duplicate]
    if len(data) < 2 or not np.isfinite(data).all() or not (np.diff(data[:, 0]) > 0).all():
        raise ValueError(f'Invalid storage time/data: {path}')
    # Source coordinates retain a stale nColumns=40; actual rows/header agree.
    return names, data, dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        source_metadata=metadata, actual_columns=len(names), rows=len(data), identical_duplicate_rows_removed=duplicate_count,
        time_range_s=[float(data[0, 0]), float(data[-1, 0])])


def _interpolate(data, t):
    t = float(t)
    if not math.isfinite(t) or not data[0, 0] <= t <= data[-1, 0]:
        raise ValueError(f'Source time {t} outside [{data[0, 0]}, {data[-1, 0]}]; no implicit looping')
    i = min(max(int(np.searchsorted(data[:, 0], t, side='right')) - 1, 0), len(data) - 2)
    fraction = (t - data[i, 0]) / (data[i + 1, 0] - data[i, 0])
    return data[i] + fraction * (data[i + 1] - data[i])


@dataclass
class GaitReference:
    coordinate_names: tuple
    coordinate_data: np.ndarray
    muscle_columns: dict
    muscle_data: np.ndarray
    native_muscles: tuple
    provenance: dict

    @classmethod
    def load(cls, root=None):
        root = Path(root) if root is not None else Path(__file__).resolve().parents[2]
        coordinates, q, qp = _read(root / _COORDS)
        controls, u, up = _read(root / _CMC)
        model = ET.parse(root / _MODEL)
        dependent = {x.findtext('dependent_coordinate_name').strip()
                     for x in model.findall('.//CoordinateCouplerConstraint')
                     if x.findtext('isEnforced', 'true').strip() == 'true'}
        excluded_dependents = {}
        for name in sorted(dependent):
            column = next((i for i, path in enumerate(coordinates) if path.endswith('/' + name + '/value')), None)
            if column is not None:
                excluded_dependents[name] = dict(reason='Enforced dependent coordinate; native constraint determines its value',
                    source_range=[float(q[:, column].min()), float(q[:, column].max())])
                if name in ('knee_angle_r_beta', 'knee_angle_l_beta'):
                    primary = coordinates.index('/jointset/walker_knee_' + name.removesuffix('_beta')[-1] + '/' + name.removesuffix('_beta') + '/value')
                    error = float(np.max(np.abs(q[:, column] - np.degrees(q[:, primary]))))
                    excluded_dependents[name]['maximum_error_vs_primary_degrees'] = error
                    excluded_dependents[name]['source_unit_anomaly'] = 'Dependent beta equals degrees(primary radians), despite global inDegrees=no' if error < 1e-8 else 'Relationship changed; dependent channel excluded without conversion'
        native = tuple(x.attrib['name'] for x in model.iter()
                       if 'Muscle' in x.tag and 'name' in x.attrib)
        selected = []
        for i, path in enumerate(coordinates[1:], 1):
            parts = path.split('/')
            if len(parts) != 5 or parts[1] != 'jointset' or parts[-1] != 'value':
                raise ValueError(f'Unexpected coordinate path {path}')
            if not parts[-2].startswith('pelvis_') and parts[-2] not in dependent:
                selected.append((i, parts[-2]))
        if len({name for _, name in selected}) != len(selected):
            raise ValueError('Duplicate joint names')
        mapping = {}
        source_names = {}
        for i, path in enumerate(controls[1:], 1):
            if not path.endswith('.excitation'):
                raise ValueError(f'Unexpected CMC control shape {path}')
            source = path.removesuffix('.excitation')
            target = source
            for old, new in _ALIASES.items():
                if source.startswith(old):
                    target = new + source[len(old):]
                    break
            if target in native:
                if target in mapping:
                    raise ValueError(f'Ambiguous muscle mapping {target}')
                if np.any(u[:, i] < 0) or np.any(u[:, i] > 1):
                    raise ValueError(f'Invalid excitation {source}')
                mapping[target] = i
                source_names[target] = source
        if len(native) != 80 or len(mapping) != 72:
            raise ValueError('Source muscle inventory changed; audit mapping before transfer')
        return cls(tuple(name for _, name in selected), q[:, [0] + [i for i, _ in selected]],
                   mapping, u, native, dict(coordinates=qp, muscle_seed=up,
                   muscle_mapping=source_names, excluded_dependent_coordinates=excluded_dependents, joint_units='radians', excitation_units='dimensionless [0,1]',
                   periodic=False, synchronized_trials=False,
                   limitations=['CMC source uses external GRFs and residual assistance; neither transferred.',
                                'Different plants and subjects: excitations are optimization seeds, not validated native controls.',
                                'Adductor magnus compartments are not one-to-one; missing values remain None.',
                                'No asserted stride period or endpoint continuity; phase maps finite recording ranges.']))

    def joint_targets(self, source_time_s):
        """Joint-name -> radian angle. Feed a controller; do not assign plant q."""
        values = _interpolate(self.coordinate_data, source_time_s)[1:]
        return dict(zip(self.coordinate_names, map(float, values)))

    def muscle_seed(self, source_time_s, muscle_names=None):
        """Muscle-name -> excitation or None; extra native muscles stay explicit."""
        values = _interpolate(self.muscle_data, source_time_s)
        names = self.native_muscles if muscle_names is None else tuple(muscle_names)
        if len(set(names)) != len(names):
            raise ValueError('Duplicate requested muscle names')
        if any(not isinstance(name, str) or not name for name in names):
            raise ValueError('Expected bare nonempty muscle names')
        if any('/' in name or ':' in name or '.' in name for name in names):
            raise ValueError('Expected bare muscle names, not paths or receptor IDs')
        return {name: float(values[self.muscle_columns[name]]) if name in self.muscle_columns else None for name in names}

    def sample_phase(self, phase):
        """Explicit normalized RECORDING phase, not an inferred gait cycle."""
        phase = float(phase)
        if not math.isfinite(phase) or not 0 <= phase <= 1:
            raise ValueError('Recording phase must be in [0,1]; no wrapping')
        qt = self.coordinate_data[0, 0] + phase * np.ptp(self.coordinate_data[:, 0])
        ut = self.muscle_data[0, 0] + phase * np.ptp(self.muscle_data[:, 0])
        return dict(joint_targets=self.joint_targets(qt), muscle_seed=self.muscle_seed(ut),
                    coordinate_source_time_s=float(qt), muscle_source_time_s=float(ut))
