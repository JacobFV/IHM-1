"""Causal environment checks; synthetic fixtures are explicitly not native validation."""
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from ihm.assembly.environment_dynamics import (EnvironmentDynamics, SpringMesh, RigidProp,
                                               box_contact, resolve_selection)

ROOT=Path(__file__).resolve().parents[1]


def fixture(points, *, gravity=(0,0,0), soft=(), rigid=(), static=()):
    owner=EnvironmentDynamics.__new__(EnvironmentDynamics)
    owner.ids=['bone']*len(points);owner.offsets=np.array(points,float);owner.previous=None
    owner.soft=list(soft);owner.rigid=list(rigid);owner.static=list(static);owner.gravity=np.array(gravity,float)
    owner.axis=2;owner.base_plane=-1.;owner.time_s=0.;owner.contacts=[];owner.last_impulse=np.zeros(3)
    owner.selection={};owner.placements=[];owner.sources={}
    return owner

ENTITY={'bone':{'centroid_m':[0,0,0],'rotation_matrix':np.eye(3).tolist()}}


class Tests(unittest.TestCase):
    def test_selection_rejects_wrong_base_and_unknown_objects(self):
        for selection in ({'scene':'bedroom'},{'objects':['missing']},{'objects':['ball-small']*17}):
            with self.assertRaises(ValueError):resolve_selection(ROOT,'upright',selection)
        selection,options=resolve_selection(ROOT,'supine',{'scene':'bedroom','bed_support_model':'bed-support-skin-quadrature','mattress_material':'mattress-soft'})
        self.assertEqual(options['bed_material'],'SM');self.assertTrue(options['surface_contact_manifest'])
        with self.assertRaises(ValueError):resolve_selection(ROOT,'supine',{'mattress_material':'mattress-soft'})

    def test_box_interior_and_exterior_contact_normals(self):
        d,n=box_contact(np.array([[1.01,0,0],[.9,0,0],[2,0,0]]),np.array([-1,-1,-1]),np.array([1,1,1]),.02)
        np.testing.assert_allclose(d,[.01,.12,0]);np.testing.assert_allclose(n[:,0],1)

    def test_cloth_body_reciprocal_impulse_and_checkpoint(self):
        mesh=SpringMesh('blanket-1','cloth',[-.1,-.1,.02],[.1,.1,.03],1.5)
        owner=fixture([[0,0,0]],soft=[mesh]);checkpoint=owner.checkpoint()
        ports=owner.advance(.02,ENTITY)
        self.assertTrue(ports);self.assertLess(sum(p['force_n'][2] for p in ports),0)
        # Internal spring forces cancel; damping is the only external loss.
        momentum=np.sum(mesh.v,axis=0)*mesh.mass
        np.testing.assert_allclose(momentum,-owner.last_impulse,rtol=.01,atol=1e-8)
        self.assertGreater(np.max(mesh.x[:,2]),.03)
        accepted=owner.frame();owner.restore(checkpoint);owner.advance(.02,ENTITY)
        self.assertEqual(owner.frame(),accepted)

    def test_pillow_lattice_deforms_and_recovers(self):
        mesh=SpringMesh('pillow-1','soft',[-.2,-.15,-.12],[.2,.15,0],.8)
        owner=fixture([[0,0,.015]],soft=[mesh]);before=mesh.x.copy()
        force=owner.advance(.02,ENTITY)
        self.assertTrue(force);self.assertLess(np.min(mesh.x[:,2]-before[:,2]),0)
        np.testing.assert_array_equal(mesh.x[mesh.fixed],mesh.rest[mesh.fixed])
        error=np.linalg.norm(mesh.x-mesh.rest)
        owner.offsets[:]=[0,0,10]
        for _ in range(100):owner.advance(.02,ENTITY)
        self.assertLess(np.linalg.norm(mesh.x-mesh.rest),error)

    def test_rigid_contact_transmits_force_and_torque(self):
        prop=RigidProp('block-1',{'bounds_m':{'min':[-.1]*3,'max':[.1]*3},'mass_kg':1.},np.zeros(3))
        owner=fixture([[.105,.06,0]],rigid=[prop]);ports=owner.advance(.02,ENTITY)
        self.assertGreater(ports[0]['force_n'][0],0);self.assertLess(prop.x[0],0)
        self.assertGreater(np.linalg.norm(prop.omega),0)
        np.testing.assert_allclose(prop.v*prop.mass,-owner.last_impulse,atol=1e-10)
        np.testing.assert_allclose(prop.rotation.T@prop.rotation,np.eye(3),atol=1e-10)

    def test_free_props_exchange_equal_opposite_impulses(self):
        record={'bounds_m':{'min':[-.1]*3,'max':[.1]*3},'mass_kg':1.,'radius_m':.1}
        a=RigidProp('a',record,np.array([-.09,0,0]));b=RigidProp('b',record,np.array([.09,0,0]))
        owner=fixture([[0,10,0]],rigid=[a,b]);owner.advance(.02,ENTITY)
        self.assertLess(a.v[0],0);self.assertGreater(b.v[0],0)
        np.testing.assert_allclose(a.v*a.mass+b.v*b.mass,0,atol=1e-10)

    def test_cloth_and_prop_exchange_impulse(self):
        mesh=SpringMesh('cloth','cloth',[-.1,-.1,.085],[.1,.1,.09],1.5)
        prop=RigidProp('ball',{'bounds_m':{'min':[-.1]*3,'max':[.1]*3},'mass_kg':1.,'radius_m':.1},np.zeros(3))
        owner=fixture([[0,10,0]],soft=[mesh],rigid=[prop]);owner.advance(.02,ENTITY)
        self.assertLess(prop.v[2],0);self.assertGreater(mesh.v[:,2].sum(),0)
        np.testing.assert_allclose(prop.v*prop.mass+mesh.v.sum(axis=0)*mesh.mass,0,atol=.003)

    def test_fixed_furniture_opposes_body_penetration(self):
        owner=fixture([[.095,0,0]],static=[(np.full(3,-.1),np.full(3,.1),'chair')])
        ports=owner.advance(.02,ENTITY);self.assertGreater(ports[0]['force_n'][0],0)
        owner.offsets[:]=[1,0,0];self.assertFalse(owner.advance(.02,ENTITY))

    def test_environment_force_reaches_embodied_plant_and_rolls_back(self):
        from scripts.verify_embodied_runtime import Plant,Neural,Native,Exchange,Load
        from ihm.assembly.embodied import EmbodiedRuntime
        class RecordingPlant(Plant):
            def advance(self,dt_s,forces=(),actuation=None):
                self.received=forces
                return super().advance(dt_s,forces,actuation)
        env=fixture([[.095,0,0]],static=[(np.full(3,-.1),np.full(3,.1),'chair')])
        plant=RecordingPlant();body=EmbodiedRuntime(plant,Neural(),Native(),Exchange(),Load(),environment_dynamics=env)
        body.mechanical_state['entities']=ENTITY
        result=body.step({'seconds':.02})
        self.assertTrue(plant.received);self.assertTrue(result['environment_state']['contact_count'])
        self.assertGreater(result['respiratory_load']['external_pressure_pa'],0)
        # A pre-native failure restores the environment along with the plant.
        plant=RecordingPlant();env=fixture([[.095,0,0]],static=[(np.full(3,-.1),np.full(3,.1),'chair')])
        body=EmbodiedRuntime(plant,Neural(),Native(),Exchange(),Load(),environment_dynamics=env);body.mechanical_state['entities']=ENTITY
        def fail(*a,**k):raise RuntimeError('fixture neural failure')
        body.neural.step=fail
        with self.assertRaises(RuntimeError):body.step({'seconds':.02})
        self.assertEqual(env.time_s,0);self.assertEqual(env.contacts,[]);self.assertIsNone(env.previous)


if __name__=='__main__':unittest.main()
