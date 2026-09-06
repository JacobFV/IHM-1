import numpy as np
from diagnose_coupled_static_residual import tangent_descent,minimax_step

def main():
    # First coordinate is at its upper bound; second still permits descent.
    q=np.zeros(3);bounds=np.array([[-1.,0.],[-1.,1.],[-1.,1.]])
    result=tangent_descent(np.array([-1.,1.,0.]),np.array([[0.,0.,1.]]),q,bounds,.01)
    assert result['half_cost_directional_derivative']<-.009
    assert result['step'][0]<=0 and result['step'][2]==0
    step=minimax_step(np.array([[1.],[1.]]),np.array([2.,0.]),np.zeros((0,1)),np.zeros(0),np.zeros(1),np.array([[-2.,2.]]),1.)
    assert abs(step[0]+1.)<1e-12
    print('PASS: source-bounded support-tangent descent and analytic minimax solution')
if __name__=='__main__':main()
