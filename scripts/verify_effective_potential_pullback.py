import numpy as np
from verify_native_effective_potential import pullback

def main():
    names=['root','knee','beta'];independent=['root','knee'];coupling={'beta':('knee',1.)}
    physical=pullback(names,independent,coupling,[2.,3.,4.]);assert np.array_equal(physical,[2.,7.])
    reaction=pullback(names,independent,coupling,[0.,-5.,5.]);assert np.array_equal(reaction,[0.,0.])
    print('PASS: dependent-coordinate work retained and ideal coupler reaction cancels')
if __name__=='__main__':main()
