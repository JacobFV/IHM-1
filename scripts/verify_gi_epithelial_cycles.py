#!/usr/bin/env python3
"""Bounded finite-preparation regression; no native runtime."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.assembly.gi_epithelial_cycles import GIPreparation


def preparation(changes=(), **kwargs):
    n=np.full((3,5),1e-6)
    for c,s,v in changes:n[c,s]=v
    return GIPreparation(n, [1e-6]*3, [1.0,1.0], [0,0,0],310,
        {'preparation':'synthetic finite GI pools','parameters':'test priors, not measured human kinetics'},
        atp_cycles_mol=kwargs.get('atp',1e-6), atp_free_energy_j_mol=kwargs.get('work',50000))


def main():
    p=preparation([(1,0,1e-8),(1,3,1e-8)]); start=p.moles.copy()
    r=p.apply('sglt1',1e-8)
    assert r['extent_mol']==1e-8
    np.testing.assert_allclose(p.moles-start,[[-2e-8,0,0,-1e-8,0],[2e-8,0,0,1e-8,0],[0,0,0,0,0]],atol=1e-22)
    # A cotransporter cannot proceed uphill without a declared work source.
    p=preparation();before=p.moles.copy()
    try:p.apply('sglt1',1e-8)
    except ValueError:pass
    else:raise AssertionError('uphill cotransport accepted')
    np.testing.assert_array_equal(before,p.moles)
    p=preparation([(0,3,1e-10),(1,3,1e-14),(1,0,1e-14)])
    r=p.apply('sglt1',1);assert r['extent_mol']==1e-10
    assert p.moles[0,3]==0
    assert p.apply('sglt1',1)['extent_mol']==0
    p=preparation([(1,0,1e-10),(2,1,1e-8)],work=200000)
    r=p.apply('nak_pump',1);assert r['extent_mol']<=1e-10/3
    assert p.moles.min()>=0
    assert r['atp_cycles_mol']==r['extent_mol']
    assert abs(r['pump_net_outward_charge_c']-r['extent_mol']*96485.33212)<1e-16
    # Finite reservoir budget and pump work are both enforced transactionally.
    p=preparation(atp=0)
    assert p.apply('nak_pump',1)['extent_mol']==0
    p=preparation(work=1e-12)
    before=p.moles.copy()
    try:p.apply('nak_pump',1e-8)
    except ValueError:pass
    else:raise AssertionError('unfunded pump accepted')
    np.testing.assert_array_equal(before,p.moles)
    p=preparation([(0,0,2e-6)])
    r=p.apply('paracellular_na',1e-9)
    assert r['free_energy_change_j']<0
    p=preparation([(2,0,2e-6)])
    r=p.apply('paracellular_na',-1e-9)
    assert r['extent_mol']<0 and r['free_energy_change_j']<0
    p=preparation([(1,4,1e-14),(1,0,1e-7)])
    initial=p.moles.sum(axis=0)
    p.apply('b0at1',1e-10);p.apply('nak_pump',1e-11)
    np.testing.assert_allclose(p.moles.sum(axis=0),initial,atol=1e-21)
    a=p.audit();assert abs(a['energy_residual_j'])<1e-15
    assert abs(a['global_charge_change_c'])<1e-15
    assert a['maximum_species_residual_mol']<1e-20
    print('PASS: stoichiometry, tail/zero donors, pump ATP/work bounds, reversible passive flux, conservation and energy ledger')

if __name__=='__main__':main()
