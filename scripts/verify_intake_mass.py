"""Consumed native intake is an external boundary, not a scheduled mass guess."""
from copy import deepcopy
import unittest
from ihm.assembly.intake_mass import consumed_intake_delta


def initial():
    return {'schema':'ihm.native-consumed-intake.v1','owner_epoch':'fixture-owner',
        'consumed_count':0,'cumulative_mass_kg':0.,
        'cumulative_native':dict.fromkeys(('carbohydrate_g','protein_g','fat_g','sodium_g','calcium_mg','water_ml'),0.),
        'pending':None,'last_consumed':None}


def consumed():
    value=initial();parts=dict(value['cumulative_native'],carbohydrate_g=1.,protein_g=1.,fat_g=1.,water_ml=10.)
    value.update(consumed_count=1,cumulative_mass_kg=.013,cumulative_native=parts,
        last_consumed={'consumed_count':1,'meal_sequence':2,'advance_sequence':3,
            'interval_start_tick':0,'interval_end_tick':1,'native_start_s':3600.,'native_end_s':3600.02,
            'native_payload':{'name':'fixture','mass_kg':.013,**parts}})
    return value


class Checks(unittest.TestCase):
    def test_actual_consumption_delta_and_no_mutation(self):
        a,b=initial(),consumed();saved=deepcopy((a,b));delta=consumed_intake_delta(a,b)
        self.assertEqual(delta['mass_kg'],.013)
        self.assertEqual(delta['boundary_id'],['fixture-owner',1,2])
        self.assertEqual((a,b),saved)
        delta['native_payload']['water_ml']=0
        self.assertEqual(b['last_consumed']['native_payload']['water_ml'],10.)

    def test_queued_meal_adds_no_mass_and_reobservation_is_noop(self):
        a=initial();queued=deepcopy(a)
        queued['pending']={'meal_sequence':2,'native_payload':consumed()['last_consumed']['native_payload']}
        self.assertIsNone(consumed_intake_delta(a,queued))
        b=consumed();self.assertIsNone(consumed_intake_delta(b,deepcopy(b)))

    def test_epoch_count_and_component_tampering_rejected(self):
        for key,value in [('owner_epoch','different'),('consumed_count',2),('consumed_count',True),('cumulative_mass_kg',.014)]:
            b=consumed();b[key]=value
            with self.subTest(key=key,value=value),self.assertRaises(ValueError):consumed_intake_delta(initial(),b)
        b=consumed();b['cumulative_native']['water_ml']=11
        with self.assertRaises(ValueError):consumed_intake_delta(initial(),b)
        b=initial();b['cumulative_mass_kg']=.001
        with self.assertRaises(ValueError):consumed_intake_delta(initial(),b)

    def test_consumption_requires_later_advance_and_valid_clock(self):
        for key,value in [('advance_sequence',2),('interval_end_tick',0),('native_end_s',3600.03),('meal_sequence',False)]:
            b=consumed();b['last_consumed'][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):consumed_intake_delta(initial(),b)

    def test_cumulative_next_meal_and_no_internal_absorption_credit(self):
        a=consumed();b=deepcopy(a)
        b['consumed_count']=2;b['cumulative_mass_kg']=.026
        b['cumulative_native']={k:2*v for k,v in a['cumulative_native'].items()}
        b['last_consumed'].update(consumed_count=2,meal_sequence=4,advance_sequence=5,
            interval_start_tick=1,interval_end_tick=2,native_start_s=3600.02,native_end_s=3600.04)
        self.assertEqual(consumed_intake_delta(a,b)['mass_kg'],.013)
        reset=deepcopy(b);reset['last_consumed'].update(native_start_s=0.,native_end_s=.02)
        with self.assertRaises(ValueError):consumed_intake_delta(a,reset)
        b['last_consumed']['meal_sequence']=2
        with self.assertRaises(ValueError):consumed_intake_delta(a,b)

if __name__=='__main__':unittest.main()
