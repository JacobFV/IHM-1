"""No-native bound-adjacent initialization and local Newton-step fixtures."""
import numpy as np
import bounded_static_root as root

def main():
    assert hasattr(root,'interior_origin'),'Missing explicit interior origin'
    bounds=np.array([[0.,.1],[-.1,.1]])
    q=root.interior_origin([5e-13,0.],bounds)
    assert q[0]>=1e-12 and q[1]==0.
    d=root.local_step(np.eye(2),[.01,-.02],q,bounds)
    assert abs(d[1]-.02)<1e-8 and np.max(abs(d))<=.03
    assert np.all(q+d>=bounds[:,0]) and np.all(q+d<=bounds[:,1])
    d=root.local_step(np.eye(2),[-10.,-10.],q,bounds)
    assert np.max(abs(d))<=.03+1e-14
    assert hasattr(root,'constrained_local_step'),'Missing force/moment-constrained local step'
    step,diagnostic=root.constrained_local_step(np.eye(4),[-1.,-1.,-1.,-1.],np.zeros(4),np.array([[-1.,1.]]*4),
        np.array([[1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.]]),np.array([.001,-.002,.003]))
    assert np.allclose(step[:3],[-.001,.002,-.003],atol=1e-9) and abs(step[3]-.03)<1e-8
    assert diagnostic['linear_support_max']<1e-9
    assert hasattr(root,'balanced_backtracked_trial'),'Missing nonlinear support correction/filter'
    def observed(q):return {'cost':1-q[3],'support_constraints':[q[0]+q[3]**2,q[1],q[2]],'gauge_residual':[0.,0.,0.]}
    candidate,entry=root.balanced_backtracked_trial(np.zeros(4),np.array([0.,0.,0.,.03]),np.array([[-1.,1.]]*4),1.,observed,
        np.array([[1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.]]),[0,1,2])
    assert candidate is not None and max(map(abs,entry['support_constraints']))<=1e-4
    assert candidate[0]<0 and candidate[3]==.03
    candidate,entry=root.balanced_backtracked_trial(np.zeros(4),np.array([0.,0.,0.,.03]),np.array([[-1.,1.]]*4),1.,
        lambda q:{'cost':.5,'support_constraints':[0.,0.,0.],'gauge_residual':[.001,0.,0.]},np.eye(3,4),[0,1,2])
    assert candidate is None, 'Nonzero held-gauge residual was accepted'
    print('PASS: explicit interior origin, near-bound progress, strict source and local bounds')
if __name__=='__main__':main()
