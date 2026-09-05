"""Independent analytic Gaussian conditioning tests across units and deficient ranks."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.human import PopulationBelief


def belief(cov,mean=None):
    n=len(cov)
    return PopulationBelief([str(i) for i in range(n)],['u']*n,np.zeros(n) if mean is None else mean,cov,{})


def main():
    c=np.array([[1,.7,.2],[.7,2,.3],[.2,.3,3.]])
    observations={'0':(1,.1,'u'),'1':(2,.01,'u')}
    actual=belief(c).condition(observations)
    h=np.eye(3)[:2];r=np.diag([.1,.01]);k=np.linalg.solve(h@c@h.T+r,(c@h.T).T).T
    np.testing.assert_allclose(actual.mean,k@[1,2],rtol=1e-12,atol=1e-13)
    np.testing.assert_allclose(actual.cov,c-k@h@c,rtol=1e-12,atol=1e-13)
    for scale in [np.array([1e-12,1e12,1e-3]),np.array([1e12,1e-12,1e3])]:
        scaled=belief(c*scale[:,None]*scale[None,:]).condition({str(i):(v*scale[i],noise*scale[i]**2,'u') for i,(v,noise) in enumerate([(1,.1),(2,.01)])})
        np.testing.assert_allclose(scaled.mean/scale,actual.mean,rtol=1e-12,atol=1e-13)
        np.testing.assert_allclose(scaled.cov/scale[:,None]/scale[None,:],actual.cov,rtol=1e-12,atol=1e-13)
    zero=belief(np.zeros((3,3)),[2,3,4]).condition({'0':(100,1e-20,'u')})
    np.testing.assert_array_equal(zero.mean,[2,3,4]);np.testing.assert_array_equal(zero.cov,0)
    unseen=belief(np.diag([1.,2.,3.])).condition({'0':(1,1e-20,'u')})
    np.testing.assert_allclose(unseen.mean,[1,0,0],atol=1e-15)
    np.testing.assert_allclose(unseen.cov,np.diag([1e-20,2,3]),rtol=1e-12,atol=0)
    # Five exact aliases expose positive eigensolver roundoff mistaken for rank.
    for scale in [np.ones(5),np.array([1e-9,1e9,1,1e-3,1e3])]:
        posterior=belief(np.outer(scale,scale)).condition({str(i):(scale[i],1e-20*scale[i]**2,'u') for i in range(5)})
        np.testing.assert_allclose(posterior.mean/scale,1,atol=1e-14)
        np.testing.assert_allclose(posterior.cov/scale[:,None]/scale[None,:],np.full((5,5),2e-21),rtol=1e-10,atol=0)
    print('Population conditioning: analytic correlated update, mixed unit invariance, rank zero, unobserved latent variance and tiny-noise exact aliases PASS')

if __name__=='__main__':main()
