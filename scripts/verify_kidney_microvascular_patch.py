"""Bounded kidney-specific topology and multi-boundary hydraulic checks."""
from pathlib import Path
import sys,unittest,math,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.kidney_microvascular_patch import materialize_kidney_patch
class Checks(unittest.TestCase):
 def scenario(self):return {'viscosity_pa_s':.003,'peritubular_radius_m':5e-6,'arterial_pressure_pa':12000.,'venous_pressure_pa':1000.,'bowman_pressure_pa':2500.,'interstitial_pressure_pa':500.,'glomerular_hydraulic_conductance_m3_s_pa':1e-17,'peritubular_hydraulic_conductance_m3_s_pa':1e-18}
 def test_source_conditioning_and_reproducibility(self):
  g=materialize_kidney_patch(ROOT,owner='RightKidney',seed=11,scenario=self.scenario())
  self.assertEqual(g,materialize_kidney_patch(ROOT,owner='RightKidney',seed=11,scenario=self.scenario()))
  self.assertEqual(g['conditioning']['lobules'],7);self.assertEqual(g['conditioning']['afferent_diameter_um'],21.5)
  self.assertAlmostEqual(g['realized']['filtration_capillary_length_m'],.0095,places=11)
  self.assertGreater(g['realized']['cycle_rank'],0);self.assertFalse(g['independent_store']);self.assertIsNone(g['solute_fluxes'])
 def test_multiboundary_conservation_and_serial_beds(self):
  g=materialize_kidney_patch(ROOT,scenario=self.scenario());f=g['hydraulics']
  self.assertLess(f['max_internal_residual_m3_per_s'],1e-20)
  self.assertAlmostEqual(sum(f['boundary_outflow_m3_per_s'].values()),0.,delta=1e-20)
  self.assertEqual(set(g['boundary_nodes']),{'arterial','venous','bowman','interstitial'})
  self.assertGreater(g['exchange_flows_m3_s']['bowman'],0.)
  pressure=np.array(f['pressure_pa']);self.assertGreaterEqual(pressure.min(),500.-1e-7);self.assertLessEqual(pressure.max(),12000.+1e-7)
  self.assertIn('efferent_arteriole',g['edge_kind']);self.assertIn('peritubular_capillary',g['edge_kind'])
 def test_physical_lengths_and_radii(self):
  g=materialize_kidney_patch(ROOT,scenario=self.scenario())
  for points,length,radius in zip(g['edge_polylines_m'],g['edge_length_m'],g['radius_m']):
   self.assertAlmostEqual(float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum()),length,places=13);self.assertGreater(radius,0.)
 def test_human_ptc_area_conditioning(self):
  scenario=self.scenario();del scenario['peritubular_radius_m'];scenario['peritubular_radius_mode']='human_control_area_equivalent'
  g=materialize_kidney_patch(ROOT,scenario=scenario)
  self.assertAlmostEqual(g['peritubular_radius_conditioning']['realized_radius_m'],math.sqrt(71/math.pi)*1e-6,places=14)
  self.assertFalse(g['peritubular_radius_conditioning']['measured_radius_distribution'])
 def test_zero_exchange_ports_do_not_create_filtration(self):
  scenario=self.scenario();scenario.update(glomerular_hydraulic_conductance_m3_s_pa=0.,peritubular_hydraulic_conductance_m3_s_pa=0.)
  g=materialize_kidney_patch(ROOT,scenario=scenario)
  self.assertEqual(g['exchange_flows_m3_s'],{'bowman':0.,'interstitial':0.})
  self.assertIsNone(g['boundary_nodes']['bowman']);self.assertIsNone(g['boundary_nodes']['interstitial'])
  self.assertAlmostEqual(sum(g['hydraulics']['boundary_outflow_m3_per_s'].values()),0.,delta=1e-20)
  # Removing the efferent arteriole must disconnect the two serial beds.
  adjacency={i:[] for i in range(len(g['positions_m']))}
  for (a,b),kind in zip(g['edges'],g['edge_kind']):
   if kind!='efferent_arteriole':adjacency[a].append(b);adjacency[b].append(a)
  seen={0};pending=[0]
  while pending:
   for node in adjacency[pending.pop()]:
    if node not in seen:seen.add(node);pending.append(node)
  self.assertNotIn(g['boundary_nodes']['venous'],seen)
 def test_invalid_requests_fail(self):
  for kw in [{'owner':'Brain'},{'seed':True},{'scenario':{}},{'scenario':dict(self.scenario(),peritubular_radius_m=-1)}]:
   with self.assertRaises(ValueError):materialize_kidney_patch(ROOT,**kw)
if __name__=='__main__':unittest.main()
