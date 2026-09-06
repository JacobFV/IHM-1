"""Endpoint mass feedback cannot replay a consumed physiology boundary."""
from copy import deepcopy
import unittest
from ihm.assembly.intake_mass import IntakeMassBridge
from scripts.verify_intake_mass import initial,consumed

class Plant:
    def __init__(self):
        self.state={'time_s':.02,'mass_transfer':{'enabled':True,'reference_id':'mechanical-epoch',
            'last_sequence':0,'owned_payload_mass_kg':0.,'owners':{}}};self.calls=[];self.fail=False
    def snapshot(self):return deepcopy(self.state)
    def body_point(self,**kwargs):
        return {'kind':'body_point','time_s':self.state['time_s'],'velocity_source_m_s':[.1,0.,0.],**kwargs}
    def transfer_mass(self,**kwargs):
        self.calls.append(kwargs)
        if self.fail:raise RuntimeError('Uncertain command result')
        m=self.state['mass_transfer'];m.update(last_sequence=kwargs['sequence'],owned_payload_mass_kg=m['owned_payload_mass_kg']+kwargs['delta_mass_kg'])
        m['last_receipt']={k:kwargs[k] for k in ('sequence','owner','body','delta_mass_kg')}
        m['owners'][kwargs['owner']]={'body':kwargs['body'],'station_m':kwargs['station_m'],'mass_kg':m['owned_payload_mass_kg']}
        return self.snapshot()

class Checks(unittest.TestCase):
    def bridge(self,plant):
        return IntakeMassBridge(plant,initial(),body='torso',station_m=[.04,.12,-.03],
            registration_identity='explicit-fixture-registration',incoming_velocity_basis='co_moving_at_ingestion_assumption')
    def test_once_only_endpoint_transfer_and_identity(self):
        plant=Plant();bridge=self.bridge(plant);r=bridge.apply(consumed())
        self.assertTrue(r['mechanical_transfer_applied']);self.assertEqual(plant.calls[0]['delta_mass_kg'],.013)
        self.assertEqual(plant.calls[0]['velocity_source_m_s'],[.1,0.,0.])
        self.assertIsNone(bridge.apply(consumed()));self.assertEqual(len(plant.calls),1)
        self.assertEqual(bridge.snapshot()['last_boundary']['boundary_id'],['fixture-owner',1,2])
    def test_false_point_or_missing_owner_receipt_rejected(self):
        plant=Plant();bridge=self.bridge(plant)
        plant.body_point=lambda **kw:{'kind':'body_point','time_s':.02,'body':'humerus_l','station_m':[100.,0.,0.],'velocity_source_m_s':[50.,0.,0.]}
        with self.assertRaises(ValueError):bridge.apply(consumed())
        self.assertEqual(plant.calls,[])
        for enabled in (False,True):
            plant=Plant();bridge=self.bridge(plant);real=plant.transfer_mass
            def corrupt(**kw):
                result=real(**kw);result['mass_transfer'].update(enabled=enabled,owners={});return result
            plant.transfer_mass=corrupt
            with self.subTest(enabled=enabled),self.assertRaises(ValueError):bridge.apply(consumed())
            self.assertTrue(bridge.snapshot()['failed'])
            self.assertIsNotNone(bridge.snapshot()['pending_boundary'])

    def test_zero_mass_boundary_needs_no_native_impulse(self):
        plant=Plant();bridge=self.bridge(plant);zero=consumed()
        zero['cumulative_mass_kg']=0.;zero['cumulative_native']=dict.fromkeys(zero['cumulative_native'],0.)
        zero['last_consumed']['native_payload'].update(zero['cumulative_native'],mass_kg=0.)
        result=bridge.apply(zero)
        self.assertFalse(result['mechanical_transfer_applied'])
        self.assertTrue(result['boundary_accounted'])
        self.assertEqual(plant.calls,[])

    def test_uncertainty_terminal_no_retry(self):
        plant=Plant();bridge=self.bridge(plant);plant.fail=True
        with self.assertRaises(RuntimeError):bridge.apply(consumed())
        plant.fail=False
        with self.assertRaisesRegex(RuntimeError,'failed'):bridge.apply(consumed())
        self.assertEqual(len(plant.calls),1)
        self.assertEqual(bridge.snapshot()['pending_boundary']['boundary_id'],['fixture-owner',1,2])
    def test_epoch_sequence_time_or_budget_mismatch_before_mutation(self):
        for mutate in (lambda p:p.state['mass_transfer'].update(reference_id='changed'),
                       lambda p:p.state['mass_transfer'].update(last_sequence=1),
                       lambda p:p.state.update(time_s=.04)):
            plant=Plant();bridge=self.bridge(plant);mutate(plant)
            with self.assertRaises(ValueError):bridge.apply(consumed())
            self.assertEqual(plant.calls,[])
        plant=Plant();bridge=self.bridge(plant);big=consumed()
        big['cumulative_mass_kg']=.6;big['cumulative_native']=dict.fromkeys(big['cumulative_native'],0.)
        big['cumulative_native']['water_ml']=600.;big['last_consumed']['native_payload'].update(big['cumulative_native'],mass_kg=.6)
        with self.assertRaises(ValueError):bridge.apply(big)
        self.assertEqual(plant.calls,[])
    def test_historical_consumption_or_missing_explicit_assumption_rejected(self):
        with self.assertRaises(ValueError):IntakeMassBridge(Plant(),consumed(),body='torso',station_m=[0,0,0],registration_identity='fixture',incoming_velocity_basis='co_moving_at_ingestion_assumption')
        with self.assertRaises(ValueError):IntakeMassBridge(Plant(),initial(),body='torso',station_m=[0,0,0],registration_identity='fixture',incoming_velocity_basis='measured')

if __name__=='__main__':unittest.main()
