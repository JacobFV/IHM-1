"""Causal environment checks; synthetic fixtures are explicitly not native validation."""
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from ihm.assembly.environment_dynamics import (EnvironmentDynamics, SpringMesh, RigidProp,
                                               box_contact, resolve_selection, top_support_heights)

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
    def test_world_frame_body_force_ports_are_canonical(self):
        from ihm.assembly.world_frame import GravityAlignedWorldFrame
        transform=np.eye(4);transform[:3,:3]=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        frame=GravityAlignedWorldFrame(transform,source_to_world_rotation=np.eye(3),source_gravity_m_s2=[0,0,0])
        point=frame.points_to_canonical([.095,0,0]);owner=fixture([point],static=[(np.full(3,-.1),np.full(3,.1),'wall')]);owner.world_frame=frame
        ports=owner.advance(.02,ENTITY);self.assertTrue(ports)
        np.testing.assert_allclose(ports[0]['point_m'],point,atol=1e-12)
        self.assertGreater(ports[0]['force_n'][1],0);self.assertAlmostEqual(ports[0]['force_n'][0],0)
        np.testing.assert_allclose(np.sum([p['force_n'] for p in ports],axis=0)*.02,owner.last_impulse,atol=1e-12)

    def test_canonical_object_force_binds_once_in_internal_world(self):
        from ihm.assembly.world_frame import GravityAlignedWorldFrame
        transform=np.eye(4);transform[:3,:3]=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        frame=GravityAlignedWorldFrame(transform,source_to_world_rotation=np.eye(3),source_gravity_m_s2=[0,0,0])
        prop=RigidProp('ball',{'bounds_m':{'min':[-.1]*3,'max':[.1]*3},'mass_kg':1.,'radius_m':.1},np.zeros(3))
        owner=fixture([frame.points_to_canonical([0,10,0])],rigid=[prop]);owner.world_frame=frame
        raw=[{'id':'ball','force_n':frame.vectors_to_canonical([1,0,0]).tolist(),'point_m':frame.points_to_canonical([0,.1,0]).tolist()}]
        bound=owner.validate_object_forces(raw);np.testing.assert_allclose(bound[0]['force_n'],[1,0,0]);np.testing.assert_allclose(bound[0]['local_point_m'],[0,.1,0])
        owner.advance(.02,ENTITY,bound);np.testing.assert_allclose(prop.v,[.02,0,0],atol=1e-12);self.assertLess(prop.omega[2],0)

    def test_top_support_is_bounded_and_does_not_lift_cloth_onto_ceiling(self):
        mattress=(np.array([-.6,-.9,-.44]),np.array([.6,.9,-.24]),'mattress')
        ceiling=(np.array([-2,-2,1.59]),np.array([2,2,1.67]),'ceiling')
        points=np.array([[0,0,.2],[1,0,.2],[0,0,1.7]])
        np.testing.assert_allclose(top_support_heights(points,-.86,[mattress,ceiling]),[-.24,-.86,1.67])

    def test_object_mattress_does_not_duplicate_native_body_support(self):
        owner=fixture([[0,0,-.3]])
        owner.object_static=[(np.array([-.6,-.9,-.44]),np.array([.6,.9,-.24]),'mattress')]
        owner.object_plane=-.86
        self.assertEqual(owner.advance(.02,ENTITY),[])

    def test_ball_outside_bed_reaches_room_floor(self):
        prop=RigidProp('ball',{'bounds_m':{'min':[-.05]*3,'max':[.05]*3},'mass_kg':.4,'radius_m':.05},np.array([1.,0,0]))
        owner=fixture([[0,10,0]],gravity=(0,0,-9.81),rigid=[prop]);owner.object_plane=-.86
        owner.object_static=[(np.array([-.6,-.9,-.44]),np.array([.6,.9,-.24]),'mattress')]
        lowest=0.;bounce=False
        for _ in range(40):
            previous=prop.v[2];owner.advance(.02,ENTITY);lowest=min(lowest,prop.x[2])
            bounce |= previous<0 and prop.v[2]>0
        self.assertLess(lowest,-.7);self.assertTrue(bounce)
        self.assertGreaterEqual(prop.x[2],-.81-1e-12)

    def test_object_force_material_binding_and_rollback(self):
        record={'bounds_m':{'min':[-.1]*3,'max':[.1]*3},'mass_kg':1.,'radius_m':.1}
        prop=RigidProp('ball',record,np.zeros(3));owner=fixture([[0,10,0]],rigid=[prop])
        raw=[{'id':'ball','force_n':[1,0,0],'point_m':[0,.1,0]}]
        ports=owner.validate_object_forces(raw);checkpoint=owner.checkpoint()
        owner.advance(.02,ENTITY,ports)
        np.testing.assert_allclose(prop.v,[.02,0,0],atol=1e-12)
        self.assertLess(prop.omega[2],0)
        accepted=owner.frame();owner.restore(checkpoint);owner.advance(.02,ENTITY,ports)
        self.assertEqual(owner.frame(),accepted)
        for bad in ([{'id':'missing','force_n':[1,0,0],'point_m':[0,0,0]}],
                    [{'id':'ball','force_n':[101,0,0],'point_m':[0,0,0]}],
                    [{'id':'ball','force_n':[1,0,0],'point_m':[2,0,0]}],
                    [{'id':'ball','force_n':[float('nan'),0,0],'point_m':[0,0,0]}]):
            with self.assertRaises(ValueError):owner.validate_object_forces(bad)
        self.assertEqual(owner.frame(),accepted)

    def test_cloth_applied_force_uses_mobile_material_node(self):
        mesh=SpringMesh('blanket','cloth',[-.2,-.2,0],[.2,.2,.001],1.5)
        mesh.fixed[0]=True;owner=fixture([[0,10,0]],soft=[mesh]);node=len(mesh.x)//2
        raw=[{'id':mesh.id,'force_n':[0,0,.1],'point_m':mesh.x[node].tolist()}]
        ports=owner.validate_object_forces(raw);self.assertEqual(ports[0]['node_index'],node)
        owner.advance(.02,ENTITY,ports)
        self.assertGreater(mesh.x[node,2],.001)
        self.assertGreater(mesh.v[:,2].sum()*mesh.mass,0)
        self.assertEqual(mesh.frame()['fixed_nodes'],[0])
        with self.assertRaises(ValueError):owner.validate_object_forces([{'id':mesh.id,'force_n':[0,0,.1],'point_m':mesh.x[0].tolist()}])
        before=owner.frame()
        with self.assertRaises(ValueError):owner.advance(.02,ENTITY,[{'id':mesh.id,'force_n':[0,0,.1],'node_index':-1}])
        self.assertEqual(owner.frame(),before)

    def test_ball_drop_rebound_and_sustained_settling(self):
        record={'bounds_m':{'min':[-.1]*3,'max':[.1]*3},'mass_kg':.4,'radius_m':.1}
        prop=RigidProp('ball',record,np.array([0,0,.1]))
        owner=fixture([[0,10,0]],gravity=(0,0,-9.81),rigid=[prop])
        initial_height=prop.x[2]-(owner.base_plane+prop.radius)
        initial_energy=.5*prop.mass*(prop.v@prop.v)+prop.mass*9.81*prop.x[2]
        rebound=[];seen_impact=False;max_residual=0.
        for _ in range(500):
            previous_v=prop.v[2];owner.advance(.02,ENTITY)
            if prop.v[2]>0 and previous_v<0:seen_impact=True
            if seen_impact and prop.v[2]>0:rebound.append(prop.x[2]-(owner.base_plane+prop.radius))
            energy=.5*prop.mass*(prop.v@prop.v)+prop.mass*9.81*prop.x[2]+prop.contact_dissipation_j
            max_residual=max(max_residual,abs(energy-initial_energy))
        self.assertTrue(seen_impact)
        self.assertAlmostEqual(max(rebound)/initial_height,.75**2,delta=.003)
        self.assertLess(max_residual,1e-8)
        self.assertAlmostEqual(prop.x[2],owner.base_plane+prop.radius,places=9)
        np.testing.assert_allclose(prop.v,0,atol=1e-10)
        self.assertEqual(prop.frame()['contact_material']['restitution'],.75)

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

    def test_partial_world_exchange_failure_restores_all_owners(self):
        from copy import deepcopy
        from scripts.verify_embodied_runtime import Plant,Neural,Native,Exchange,Load
        from ihm.assembly.embodied import EmbodiedRuntime
        class FailingPlant(Plant):
            def advance(self,dt_s,forces=(),actuation=None):
                frame=super().advance(dt_s,forces,actuation)
                if len(self.commands)==3:raise RuntimeError('third mechanical substep failed')
                frame['entities']=deepcopy(ENTITY)
                return frame
        plant=FailingPlant();neural=Neural();native=Native()
        env=fixture([[.095,0,0]],static=[(np.full(3,-.1),np.full(3,.1),'chair')])
        body=EmbodiedRuntime(plant,neural,native,Exchange(),Load(),environment_dynamics=env)
        body.mechanical_state['entities']=deepcopy(ENTITY)
        before=deepcopy(body.mechanical_state);commands=deepcopy(body.next_excitation)
        with self.assertRaisesRegex(RuntimeError,'third mechanical substep'):
            body.step({'seconds':.02})
        self.assertEqual(plant.t,0);self.assertEqual(plant.commands,[])
        self.assertEqual(neural.t,0);self.assertEqual(native.t,0)
        self.assertEqual(env.time_s,0);self.assertEqual(env.contacts,[])
        self.assertIsNone(env.previous)
        self.assertEqual(body.mechanical_state,before)
        self.assertEqual(body.next_excitation,commands)
        self.assertEqual(body.time_s,0);self.assertFalse(body.failed)

    def test_environment_force_reaches_embodied_plant_and_rolls_back(self):
        from scripts.verify_embodied_runtime import Plant,Neural,Native,Exchange,Load
        from ihm.assembly.embodied import EmbodiedRuntime
        class RecordingPlant(Plant):
            def advance(self,dt_s,forces=(),actuation=None):
                self.received=forces
                frame=super().advance(dt_s,forces,actuation)
                frame['entities']=ENTITY
                return frame
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
