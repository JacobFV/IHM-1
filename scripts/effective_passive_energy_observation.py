"""Source-derived diagnostic correction; never writes native energy or forces."""
from pathlib import Path
import hashlib
import math
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
THELEN_SOURCE = ROOT / 'data/raw/mechanics/opensim-core/OpenSim/Actuators/Thelen2003Muscle.cpp'


def fiber_primitive(x, lo, fiso, k, strain):
    if not all(math.isfinite(v) for v in (x, lo, fiso, k, strain)) or min(lo, fiso, k, strain) <= 0:
        raise ValueError('Invalid source passive primitive parameters')
    d = max(0., x - 1.)
    return fiso * lo * (strain / k * math.expm1(k * d / strain) - d) / math.expm1(k)


class PassiveEnergyObservation:
    """Bind the known getter defect to actual source and model parameters."""
    def __init__(self, model):
        model = Path(model)
        source = THELEN_SOURCE.read_text()
        if 'calcfpefisoPE(mli.fiberLength)' not in source:
            raise ValueError('Retained Thelen getter defect no longer matches correction')
        self.parameters = {
            m.get('name'): {k: float(m.findtext(k)) for k in
                ('optimal_fiber_length', 'max_isometric_force', 'KshapePassive', 'FmaxMuscleStrain')}
            for m in ET.parse(model).getroot().iter('Thelen2003Muscle')}
        self.provenance = dict(model_sha256=hashlib.sha256(model.read_bytes()).hexdigest(),
            thelen_source_sha256=hashlib.sha256(THELEN_SOURCE.read_bytes()).hexdigest(),
            correction_basis='Thelen source calcfpefisoPE: normalized fiber length replaces erroneous metre argument; tendon primitive unchanged',
            scope='Diagnostic passive potential observation only; native getter, forces and physiological ledger unchanged')

    def muscle(self, observation):
        raw = observation['passive_energy_j']
        correction = 0.
        if observation['type'] == 'Thelen2003Muscle':
            p = self.parameters[observation['name']]
            lo = p['optimal_fiber_length']; lf = observation['fiber_length_m']
            args = (lo, p['max_isometric_force'], p['KshapePassive'], p['FmaxMuscleStrain'])
            correction = fiber_primitive(lf / lo, *args) - fiber_primitive(lf, *args)
        elif observation['type'] != 'Millard2012EquilibriumMuscle':
            raise ValueError('Unsupported passive energy observation type')
        return dict(name=observation['name'], raw_native_passive_energy_j=raw,
            source_passive_correction_j=correction, corrected_passive_energy_j=raw + correction)

    def observe(self, response):
        rows = [self.muscle(m) for m in response['muscles']]
        return dict(provenance=self.provenance, muscles=rows,
            raw_native_passive_energy_j=sum(m['raw_native_passive_energy_j'] for m in rows),
            source_passive_correction_j=sum(m['source_passive_correction_j'] for m in rows),
            corrected_passive_energy_j=sum(m['corrected_passive_energy_j'] for m in rows))
