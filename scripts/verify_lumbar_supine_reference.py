import numpy as np
from prepare_lumbar_supine_reference import retained_indices

def main():
    face=np.array([3,0,4,2]);components=[0,1,45,0,0]
    keep=retained_indices(face,components,0)
    assert keep.tolist()==[0,1,2] and face[keep].tolist()==[3,0,4]
    assert np.array_equal(retained_indices(face,[0]*5,0),np.arange(4))
    print('PASS: excluded seam omitted with original order and face identities preserved')
if __name__=='__main__':main()
