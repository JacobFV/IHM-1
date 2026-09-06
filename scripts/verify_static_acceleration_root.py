"""Source-only all-mobility static-root objective and acceptance fixtures."""
import json
import numpy as np
import solve_constrained_supine_pose as solver

def main():
    assert hasattr(solver,'acceleration_residual'),'Missing all-mobility acceleration objective'
    n={'udot':[0.,0.,1500.],'mobility_rotational':[False,True,True],
       'constraint_position_error':0.,'constraint_velocity_error':0.,'constraint_acceleration_error':0.}
    e={'native':n,'support_constraints':[0.,0.,0.],'gauge_residual':[0.,0.,0.]}
    assert np.array_equal(solver.acceleration_residual(n),[0.,0.,1500.])
    assert not solver.static_converged(e), 'Balanced weight hides toe acceleration'
    e['native']={**n,'udot':[0.,0.,0.]};assert solver.static_converged(e)
    e['gauge_residual']=[0.,.001,0.];assert not solver.static_converged(e)
    e['gauge_residual']=[0.,0.,0.];e['native']['constraint_position_error']=.001
    assert not solver.static_converged(e)
    try:solver.acceleration_residual({**n,'udot':[0.,float('nan'),0.]})
    except ValueError:pass
    else:raise AssertionError('Nonfinite root residual accepted')
    print(json.dumps({'passed':True,'native_run':False,'checks':['tiny-inertia acceleration retained','all mobility threshold','gauge rejection','constraint error rejection','nonfinite rejection']}))
if __name__=='__main__':main()
