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
    from bounded_static_root import constrained_local_step
    step,_=constrained_local_step(np.eye(4),-np.ones(4),np.zeros(4),np.array([[-1.,1.]]*4),np.eye(3,4),np.zeros(3),radius=.0075)
    assert np.max(abs(step))<=.0075+1e-12 and np.max(abs(step[:3]))<1e-12
    print('PASS: exact constrained residual floor, projected gradient and smaller-box source-preserving solve')
if __name__=='__main__':main()
