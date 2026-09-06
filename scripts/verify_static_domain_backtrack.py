"""Rejected copied-state material trials backtrack; other failures remain fatal."""
import numpy as np
import bounded_static_root as root

def main():
    assert hasattr(root,'backtracked_trial'),'Missing domain-aware local trial handling'
    seen=[]
    def evaluate(q):
        seen.append(float(q[0]))
        if q[0]>.02:raise root.StaticDomainRejection('equal-pressure skin/bed solution exceeds retained domains')
        return {'cost':1.-q[0]}
    q,entry=root.backtracked_trial(np.array([0.]),np.array([.03]),np.array([[-.1,.1]]),1.,evaluate)
    assert seen==[.03,.015] and q[0]==.015 and entry['cost']<1
    def other(q):raise ValueError('unrelated mechanical failure')
    try:root.backtracked_trial(np.array([0.]),np.array([.03]),np.array([[-.1,.1]]),1.,other)
    except ValueError:pass
    else:raise AssertionError('Unrelated failure swallowed')
    print('PASS: exact-domain trial rejection, bounded half-step retry, unrelated failure propagation')
if __name__=='__main__':main()
