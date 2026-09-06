"""Source-only executable mechanism checks; no native process or stiffness."""
import unittest
import json,hashlib
from unittest.mock import patch
from scripts.build_thoracic_mechanism import read_json_receipt,verified_gzip_json
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_mechanism import closed_volume,ThoracicMechanism

ROOT=Path(__file__).resolve().parents[1]
class ThoracicMechanismTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=ThoracicMechanism(ROOT/'data/research/thoracic_mechanism/v2/manifest.json')
    def test_closed_tetra_volume_gradient_and_reaction(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1.]])
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        volume,g=closed_volume(v,f)
        self.assertAlmostEqual(volume,1/6)
        np.testing.assert_allclose(g.sum(0),0,atol=1e-14)
        np.testing.assert_allclose(np.cross(v,g).sum(0),0,atol=1e-14)
        h=1e-7;vp=v.copy();vp[1,0]+=h;vm=v.copy();vm[1,0]-=h
        self.assertAlmostEqual(g[1,0],(closed_volume(vp,f)[0]-closed_volume(vm,f)[0])/(2*h),places=9)
        with self.assertRaises(ValueError):closed_volume(v,f[:3])


    def test_actual_full_mass_and_kinetic_energy(self):
        m=self.model;q=np.zeros(26);velocity=np.linspace(-.03,.04,32)
        for k in m.locked:velocity[6+k]=0
        zero=m.kinetic(q,velocity)
        prior=json.loads((ROOT/'data/research/cervical_inertia/v2/manifest.json').read_bytes())['residual_torso']
        for key in ['mass_kg','center_m','inertia_kg_m2']:np.testing.assert_allclose(zero[key],prior[key],rtol=0,atol=1e-10)
        np.testing.assert_allclose(zero['mass_matrix'],zero['mass_matrix'].T,atol=1e-12)
        self.assertEqual(len(zero['active_indices']),30)
        self.assertGreater(zero['active_eigenvalues'][0],0)
        self.assertAlmostEqual(zero['kinetic_energy_J'],zero['direct_material_energy_J'],places=12)
        q[0]=.003;q[24]=.0005;q[25]=.002
        moved=m.kinetic(q,velocity)
        self.assertAlmostEqual(moved['mass_kg'],prior['mass_kg'],places=11)
        self.assertAlmostEqual(moved['kinetic_energy_J'],moved['direct_material_energy_J'],places=12)
        self.assertGreater(moved['active_eigenvalues'][0],0)
    def test_actual_cavity_finite_difference_and_pressure_work(self):
        m=self.model;q=np.zeros(26);q[0]=.001;q[25]=.001
        result=m.cavity(q,pressure_pa=75,gas_reference_volume_m3=.003)
        finite=[];h=1e-7
        for k in range(26):
            if k in m.locked:finite.append(0.);continue
            step=np.zeros(26);step[k]=h
            finite.append((m.cavity(q+step)['geometric_volume_m3']-m.cavity(q-step)['geometric_volume_m3'])/(2*h))
        np.testing.assert_allclose(result['volume_jacobian_m3_per_coordinate'],finite,rtol=1e-5,atol=2e-10)
        speed=np.linspace(-.01,.02,26)
        for k in m.locked:speed[k]=0
        nodal_velocity=np.einsum('ijk,k->ij',result['material_jacobian'],speed)
        nodal_power=float(np.sum(result['pressure_nodal_forces']*nodal_velocity))
        self.assertAlmostEqual(nodal_power,float(result['pressure_generalized_force']@speed),places=12)
        self.assertAlmostEqual(nodal_power,75*float(np.array(finite)@speed),places=9)
        np.testing.assert_allclose(result['pressure_resultant_force_N'],0,atol=1e-11)
        np.testing.assert_allclose(result['pressure_resultant_moment_Nm'],0,atol=1e-11)
    def test_moving_material_anchors_and_external_reaction(self):
        m=self.model;q=np.zeros(26);q[:24]=.002;q[24]=.001;q[25]=.002
        for k in m.locked:q[k]=0
        for ident,material in m.material.items():
            if material['entry']['map_kind']!='moving_anchors':continue
            b=material['bindings'];x,j=m.material_state(ident,q,b['anchor_source_vertices'])
            expected,_=m.driver(b['anchor_reference'],b['anchor_driver'],q)
            np.testing.assert_allclose(x,expected,atol=1e-12)
        ident=next(iter(m.material));x,j=m.material_state(ident,q,[0]);force=np.array([2,-3,4.]);load=m.point_load(ident,0,q,force)
        speed=np.linspace(-.02,.02,32)
        for k in m.locked:speed[6+k]=0
        point_speed=speed[:3]+np.cross(speed[3:6],x[0])+j[0]@speed[6:]
        self.assertAlmostEqual(float(force@point_speed),float(load['generalized_force']@speed),places=12)
        np.testing.assert_allclose(load['fixed_parent_reaction_force_N']+force,0)
        np.testing.assert_allclose(load['fixed_parent_reaction_moment_Nm']+np.cross(x[0],force),0)
    def test_retained_lung_and_license_receipts(self):
        m=self.model;c=m.recipe['cavity']
        for row in c['source_lobes']+[c['upstream_license']]:
            raw=(m.path.parent/row['retained_copy']).read_bytes()
            self.assertEqual(hashlib.sha256(raw).hexdigest(),row['sha256'])
            self.assertEqual(len(raw),row['bytes'])

    def test_exact_source_node_ids(self):
        m=self.model;ident=next(iter(m.material));q=np.zeros(26)
        for bad in [-1,.9,True,np.float64(1.),len(m.material[ident]['reference'])]:
            with self.assertRaises(ValueError):m.point_load(ident,bad,q,[1,2,3])
            with self.assertRaises(ValueError):m.material_state(ident,q,[bad])
        with self.assertRaises(ValueError):m.material_state(ident,q,[True,1])

    def test_read_once_and_verify_before_parse(self):
        raw=b'{"value":1}'
        with patch.object(Path,'read_bytes',side_effect=[raw,b'{"value":2}']) as read:
            parsed,identity=read_json_receipt(ROOT/'unused.json')
        self.assertEqual(read.call_count,1)
        self.assertEqual(parsed,{'value':1})
        self.assertEqual(identity['sha256'],hashlib.sha256(raw).hexdigest())
        with patch.object(Path,'read_bytes',return_value=b'bad gzip'),patch('scripts.build_thoracic_mechanism.gzip.decompress') as decode:
            with self.assertRaises(ValueError):verified_gzip_json(ROOT/'unused.gz','0'*64)
            decode.assert_not_called()

    def test_locked_geometry_and_domain_rejection(self):
        m=self.model;q=np.zeros(26);q[11]=.001
        with self.assertRaises(ValueError):m.coordinates(q)
        q[11]=0;q[0]=.06
        with self.assertRaises(ValueError):m.coordinates(q)

if __name__=='__main__':unittest.main()
