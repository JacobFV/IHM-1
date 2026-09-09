"""Explicit coupling of the reduced SI mechanics and somatic peripheral model.

This is a numerical plant, not validated joint biomechanics. Commands enter the
peripheral delay/activation model; resulting activation drives the NEXT interval.
Local spindle/tendon rates are transduced signals, not arrived cortical input.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

from .mechanics import BodyMechanics
from .peripheral import BodyPeripheral


class MechanicalPeripheral:
    def __init__(self, peripheral, mechanics):
        self.peripheral = peripheral
        self.mechanics = mechanics
        self.ids = set(peripheral.bindings)
        if not self.ids <= {m['id'] for m, _ in mechanics.muscles}:
            raise ValueError('Peripheral and mechanics muscle IDs differ')
        self.activations = {}

    @classmethod
    def from_directory(cls, directory):
        directory = Path(directory)
        peripheral = json.loads((directory / 'peripheral.json').read_text())
        mechanics = json.loads((directory / 'mechanics.json').read_text())
        return cls(BodyPeripheral.from_dict(peripheral), BodyMechanics.from_dict(mechanics))

    def step(self, dt_s, stimuli=None, brain_state=None, drivers=None, blocked_nerves=()):
        drivers = dict(drivers or {})
        if 'activation' in drivers:
            raise ValueError('Activation is owned by peripheral motor output')
        drivers['activation'] = self.activations
        # Dynamic mechanics arrays mutate in place; all other runtime fields are
        # replaced on advance. Keep immutable assembly/operators shared. Solver
        # cache insertion is harmless and does not carry physical state.
        saved_mechanics = dict(vars(self.mechanics))
        for name in ('x', 'v', 'omega', 'r', 'deformation'):
            saved_mechanics[name] = saved_mechanics[name].copy()
        saved_peripheral = deepcopy(vars(self.peripheral))
        try:
            mechanical = self.mechanics.step(dt_s, drivers)
            out = self.peripheral.step(dt_s, stimuli, mechanical, brain_state, blocked_nerves)
            for key in ('spindle_rates_hz', 'tendon_rates_hz', 'proprioceptor_rates_hz'):
                if set(out[key]) != self.ids:
                    raise ValueError(f'{key} must contain every bare muscle ID')
        except Exception:
            vars(self.mechanics).clear()
            vars(self.mechanics).update(saved_mechanics)
            vars(self.peripheral).clear()
            vars(self.peripheral).update(saved_peripheral)
            raise
        self.activations = out['motor_activations']
        out['mechanical_state'] = mechanical
        return out
