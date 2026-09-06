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
    print('PASS: explicit interior origin, near-bound progress, strict source and local bounds')
if __name__=='__main__':main()
