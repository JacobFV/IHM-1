"""Bounded local re-solves: nonlinear rejection shrinks box, not stale direction."""
import numpy as np
import bounded_static_root as root

def main():
    A=np.eye(3,4);bounds=np.array([[-1.,1.]]*4);q=np.zeros(4);seen=[]
    def observe(x):
        seen.append(x.copy())
        return dict(cost=2. if x[3]>.004 else .5,support_constraints=[0.,0.,0.],gauge_residual=[0.,0.,0.])
    x,e,radius,trace=root.recomputed_box_trial(np.eye(4),np.array([0.,0.,0.,-1.]),q,bounds,A,np.zeros(3),1.,observe,[0,1,2],.0075)
    assert radius==.00375 and abs(x[3]-.00375)<1e-10
    assert len(trace)==2 and trace[0]['accepted'] is False and trace[1]['accepted'] is True
    assert len(seen)==2 and max(abs(seen[0]))<=.0075+1e-12
    # Support repair must stay inside the active smaller box too.
    def impossible(x):return dict(cost=.5,support_constraints=[.02,0.,0.],gauge_residual=[0.,0.,0.])
    x,e=root.balanced_backtracked_trial(q,np.array([0.,0.,0.,.0075]),bounds,1.,impossible,A,[0,1,2],radius=.0075,maximum_backtracks=1)
    assert x is None
    print('PASS: actual-objective rejection recomputes smaller bounded QP; support repair respects box')
if __name__=='__main__':main()
