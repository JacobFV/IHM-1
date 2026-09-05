"""closed-loop four-chamber circulation coupled to single-compartment breathing.

Structural equations follow time-varying elastance/Windkessel and respiratory
resistance-compliance models. Defaults are illustrative, not a fitted subject.
Fixed-rate phase oscillators supply pacemaking and respiratory drive; neither
ECG electrophysiology nor chemoreflex regulation is claimed.
"""
from dataclasses import dataclass
import math
import numpy as np

NAMES = ('right_atrium', 'right_ventricle', 'pulmonary_artery', 'pulmonary_vein',
         'left_atrium', 'left_ventricle', 'systemic_artery', 'systemic_vein')
COMPONENTS = tuple('cardiopulmonary.'+name+'.volume' for name in NAMES) + (
    'cardiopulmonary.lung.volume', 'cardiopulmonary.cardiac.phase', 'cardiopulmonary.respiratory.phase')


@dataclass(frozen=True)
class Parameters:
    heart_rate_bpm: float = 70.
    respiratory_rate_per_min: float = 12.
    systemic_resistance_mmHg_s_mL: float = 1.1
    pulmonary_resistance_mmHg_s_mL: float = .12
    systemic_venous_resistance_mmHg_s_mL: float = .04
    pulmonary_venous_resistance_mmHg_s_mL: float = .03
    valve_resistance_mmHg_s_mL: float = .01
    left_ventricular_emax_mmHg_mL: float = 2.5
    right_ventricular_emax_mmHg_mL: float = .8
    airway_resistance_cmH2O_s_L: float = 2.
    lung_compliance_L_cmH2O: float = .2
    pleural_baseline_cmH2O: float = -5.
    respiratory_effort_cmH2O: float = 3.
    lung_relaxed_volume_L: float = 1.5
    posture: str = 'supine'

    def __post_init__(self):
        for name, value in vars(self).items():
            if name == 'posture':
                if value != 'supine':
                    raise ValueError('only supine posture is implemented')
            elif not math.isfinite(value) or (name != 'pleural_baseline_cmH2O' and value <= 0):
                raise ValueError(f'invalid cardiopulmonary parameter: {name}')
        if self.left_ventricular_emax_mmHg_mL <= .06 or self.right_ventricular_emax_mmHg_mL <= .03:
            raise ValueError('ventricular maximum elastance must exceed diastolic minimum')
        if not 20 <= self.heart_rate_bpm <= 200 or not 2 <= self.respiratory_rate_per_min <= 60:
            raise ValueError('pacing outside the declared numerical regime')


def initial_state(p):
    # Sum of the eight blood volumes is exactly 5000 mL.
    lung = p.lung_relaxed_volume_L-p.pleural_baseline_cmH2O*p.lung_compliance_L_cmH2O
    return np.array([60., 120., 150., 400., 60., 120., 700., 3390., lung, 0., 0.])


def activation(phase, start, duration):
    local = (phase-start) % 1.
    return .5*(1-math.cos(2*math.pi*local/duration)) if local < duration else 0.


def signals(x, p):
    v = x[:8]; cardiac_phase, respiratory_phase = x[9], x[10]
    # Inspiratory effort is smooth; exhalation follows passive elastic recoil.
    pleural_cm = p.pleural_baseline_cmH2O-p.respiratory_effort_cmH2O*(1-math.cos(2*math.pi*respiratory_phase))/2
    pleural_mm = pleural_cm*.735559
    ventricular_activation = activation(cardiac_phase, 0., .4)
    atrial_activation = activation(cardiac_phase, .78, .2)
    elastance = np.array([.1+.15*atrial_activation,
                         .03+(p.right_ventricular_emax_mmHg_mL-.03)*ventricular_activation,
                         .2, .05, .1+.15*atrial_activation,
                         .06+(p.left_ventricular_emax_mmHg_mL-.06)*ventricular_activation,
                         1/1.5, 1/100.])
    unstressed = np.array([10., 10., 100., 300., 10., 10., 600., 2800.])
    pressure = elastance*(v-unstressed)
    pressure[:6] += pleural_mm
    # Every outgoing flow is the next compartment's incoming flow in the loop.
    resistance = np.array([p.valve_resistance_mmHg_s_mL,
                           p.valve_resistance_mmHg_s_mL,
                           p.pulmonary_resistance_mmHg_s_mL,
                           p.pulmonary_venous_resistance_mmHg_s_mL,
                           p.valve_resistance_mmHg_s_mL,
                           p.valve_resistance_mmHg_s_mL,
                           p.systemic_resistance_mmHg_s_mL,
                           p.systemic_venous_resistance_mmHg_s_mL])
    flow = (pressure-np.roll(pressure, -1))/resistance
    # Four pressure-operated ideal valves. Vascular resistors allow reversal.
    flow[[0, 1, 4, 5]] = np.maximum(flow[[0, 1, 4, 5]], 0.)
    alveolar_pressure = pleural_cm+(x[8]-p.lung_relaxed_volume_L)/p.lung_compliance_L_cmH2O
    airflow = -alveolar_pressure/p.airway_resistance_cmH2O_s_L
    return {'pressure_mmHg': pressure, 'flow_mL_s': flow, 'airflow_L_s': airflow,
            'lung_volume_L': x[8], 'pleural_pressure_cmH2O': pleural_cm,
            'alveolar_pressure_cmH2O': alveolar_pressure,
            'systemic_arterial_pressure_mmHg': pressure[6], 'aortic_flow_mL_s': flow[5],
            'net_circulatory_flow_mL_s': float(np.sum(np.roll(flow, 1)-flow))}


def rhs(time, x, p):
    s = signals(x, p); flow = s['flow_mL_s']
    return np.r_[np.roll(flow, 1)-flow, s['airflow_L_s'], p.heart_rate_bpm/60., p.respiratory_rate_per_min/60.]


def register_cardiopulmonary(r, p):
    from ihm.fields import Component, Support
    from ihm.anatomy import Partition
    from ihm.topologies import Topology
    from ihm.processes import Process
    regions = tuple('cp_'+name for name in NAMES) + ('cp_lung', 'cp_sinoatrial', 'cp_brainstem')
    r.add(Support('cardiopulmonary_blood_lumens', 'supine_compartments', regions[:8]))
    r.add(Support('cardiopulmonary_airspace', 'supine_compartments', (regions[8],)))
    r.add(Support('cardiopulmonary_pacemakers', 'supine_compartments', regions[9:]))
    r.add(Partition('cardiopulmonary_anatomy', {region: {'cardiac' if i in (0,1,4,5) else 'vascular' if i < 8 else 'respiratory' if i in (8,10) else 'cardiac': 1.} for i, region in enumerate(regions)}))
    initial = initial_state(p)
    for i, c in enumerate(COMPONENTS):
        support = 'cardiopulmonary_blood_lumens' if i < 8 else 'cardiopulmonary_airspace' if i == 8 else 'cardiopulmonary_pacemakers'
        unit = 'mL' if i < 8 else 'L' if i == 8 else 'cycle'
        r.add(Component(c, support, regions[i], unit, initial[i], max(.01, abs(initial[i])*.1)))
    inputs = []
    for i in range(11):
        indices = sorted({(i-1)%8, i, (i+1)%8, 9, 10}) if i < 8 else [8,10] if i == 8 else [i]
        inputs.append(indices)
    pairs = {(regions[j], regions[i]) for i, indices in enumerate(inputs) for j in indices}
    r.add(Topology('cardiopulmonary_cycle', tuple(sorted(pairs))))
    for i, c in enumerate(COMPONENTS):
        r.add(Process('cycle.'+c, tuple(COMPONENTS[j] for j in inputs[i]), c, 'cardiopulmonary_cycle',
                      tuple(0. for _ in inputs[i]), 0., 0.,
                      'mechanistic_form;illustrative_parameters;deterministic;supine', 'cardiopulmonary_lumped'))
