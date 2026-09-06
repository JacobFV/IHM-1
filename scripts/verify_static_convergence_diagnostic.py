"""Offline analytic fixtures; no native model execution."""
import numpy as np
from diagnose_static_convergence import linear_diagnostic

def main():
    A=np.array([[1.,0.]])
    d=linear_diagnostic(np.eye(2),np.array([2.,3.]),A,np.zeros(1))
    assert abs(d['linearized_equality_residual_floor_norm']-2)<1e-12
    assert abs(d['projected_half_cost_gradient_norm']-3)<1e-12
    assert d['tangent_rank']==1
    d=linear_diagnostic(np.eye(2),np.array([2.,0.]),A,np.zeros(1))
    assert d['projected_half_cost_gradient_norm']==0
    assert d['linearized_equality_residual_floor_norm']==2
    print('PASS: exact constrained residual floor and projected gradient')
if __name__=='__main__':main()
