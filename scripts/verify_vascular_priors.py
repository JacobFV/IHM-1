"""Bounded source-integrity, density, topology and passive-closure checks."""
from pathlib import Path
import json, tempfile, math
import numpy as np
import sys, unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.vascular_priors import load_registry, materialize_muscle, solve_passive
class Checks(unittest.TestCase):
 def test_primary_registry(self):
  r=load_registry(ROOT);p=r['organs']['Muscle']['cohorts']['young_men']
  self.assertEqual(p['density']['mean'],331);self.assertEqual(p['density']['sd'],94)
  self.assertEqual(p['domain_area']['unit'],'um2');self.assertFalse(r['organs']['Liver']['materialization_supported'])
 def test_source_tampering_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);registry=load_registry(ROOT)
   path=root/'data/sources/organ_microvascular_priors.json';path.parent.mkdir(parents=True);path.write_text(json.dumps(registry))
   source=root/registry['sources'][0]['path'];source.parent.mkdir(parents=True);source.write_text('changed')
   with self.assertRaises(ValueError):load_registry(root)
 def test_overflow_rejected(self):
  with self.assertRaises(ValueError):materialize_muscle(ROOT,section_mm=(1e308,1e308))
  with self.assertRaises(ValueError):solve_passive([[0,1]],[1e-320],{0:1.,1:0.})
 def test_deterministic_density_not_measured_anatomy(self):
  a=materialize_muscle(ROOT,cohort='young_men',section_mm=(.2,.2),extent_mm=.5,seed=19)
  b=materialize_muscle(ROOT,cohort='young_men',section_mm=(.2,.2),extent_mm=.5,seed=19)
  self.assertEqual(a,b);self.assertEqual(a['density_constraint']['plane_intersections'],13)
  self.assertLessEqual(abs(a['density_constraint']['residual_per_mm2']),.5/.04+1e-8)
  self.assertIsNone(a['radius_m']);self.assertIsNone(a['flow_solution']);self.assertFalse(a['independent_store'])
  self.assertEqual(a['native_owner'],'Muscle');self.assertIn('domain_area_distribution',a['unmet_constraints'])
 def test_bad_or_oversized_requests_fail(self):
  for k in [{'section_mm':(0,.2)},{'section_mm':(20,20)},{'seed':True},{'extent_mm':float('nan')},{'density_per_mm2':-1}]:
   with self.assertRaises(ValueError):materialize_muscle(ROOT,**k)
 def test_passive_analytic_parallel_and_conservation(self):
  r=solve_passive([[0,1],[1,2],[0,2]],[2.,2.,4.],{0:8.,2:0.})
  self.assertAlmostEqual(r['pressure_pa'][1],4.);self.assertEqual(r['flow_m3_per_s'],[2.,2.,2.])
  self.assertLess(r['max_internal_residual_m3_per_s'],1e-12);self.assertAlmostEqual(sum(r['boundary_outflow_m3_per_s'].values()),0.)
 def test_materialized_graph_passive_closure(self):
  g=materialize_muscle(ROOT,section_mm=(.1,.1),seed=7)
  r=solve_passive(g['edges'],[1e15]*len(g['edges']),{0:4000.,1:1000.})
  self.assertLess(r['max_internal_residual_m3_per_s'],1e-20)
  self.assertGreater(r['boundary_outflow_m3_per_s']['0'],0.)
 def scenario(self):
  return {'radius_mode':'human_quadriceps_1988_equal_area','diameter_state':'shrinkage_corrected','capillary_radius_cv':.2,'supply_radius_multiplier':2.,'return_radius_multiplier':2.4,'viscosity_pa_s':.003,'pressure_boundaries_pa':{0:4000.,1:1000.}}
 def test_source_conditioned_radius_and_flow_realization(self):
  scenario=self.scenario();g=materialize_muscle(ROOT,seed=13,hydraulic_scenario=scenario)
  self.assertEqual(g,materialize_muscle(ROOT,seed=13,hydraulic_scenario=scenario))
  r=np.array(g['radius_m']);k=np.array(g['edge_kind']);self.assertTrue(np.all(r>0))
  expected=math.sqrt(5.31*3.62)*.5e-6
  self.assertAlmostEqual(float(r[k==2].mean())/expected,1.)
  self.assertEqual(g['radius_conditioning']['target_sd_is_measured'],False)
  f=g['flow_solution'];self.assertLess(f['max_internal_residual_m3_per_s'],1e-20)
  pressure=np.array(f['pressure_pa']);self.assertGreaterEqual(pressure.min(),1000-1e-8);self.assertLessEqual(pressure.max(),4000+1e-8)
  edges=np.array(g['edges']);self.assertTrue(np.all(pressure[edges[:,0]]>=pressure[edges[:,1]]-1e-8))
  self.assertAlmostEqual(sum(f['boundary_outflow_m3_per_s'].values()),0.,delta=1e-20)
  self.assertGreater(g['geometric_lumen_volume_m3'],0.);self.assertFalse(g['independent_store'])
 def test_engineering_option_explicit_and_radius_scaling(self):
  scenario=self.scenario();scenario.update(radius_mode='engineering_lognormal',capillary_radius_mean_m=3e-6,capillary_radius_cv=0.)
  g=materialize_muscle(ROOT,hydraulic_scenario=scenario)
  scenario['capillary_radius_mean_m']=6e-6;h=materialize_muscle(ROOT,hydraulic_scenario=scenario)
  self.assertAlmostEqual(h['flow_solution']['boundary_outflow_m3_per_s']['0']/g['flow_solution']['boundary_outflow_m3_per_s']['0'],16.)
  self.assertAlmostEqual(h['geometric_lumen_volume_m3']/g['geometric_lumen_volume_m3'],4.)
 def test_incomplete_hydraulics_rejected(self):
  for value in [{},{'radius_mode':'engineering_lognormal'},dict(self.scenario(),capillary_radius_cv=-1),dict(self.scenario(),radius_mode='unknown')]:
   with self.assertRaises(ValueError):materialize_muscle(ROOT,hydraulic_scenario=value)
 def test_disconnected_unanchored_and_invalid_resistance_fail(self):
  for edges,res,bound in [([[0,1],[2,3]],[1.,1.],{0:1.,1:0.}),([[0,1]],[0.],{0:1.,1:0.})]:
   with self.assertRaises(ValueError):solve_passive(edges,res,bound)
if __name__=='__main__':unittest.main()
