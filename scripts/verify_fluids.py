"""volume conservation is independent of physiological parameter calibration."""
import numpy as np
from ihm.materialize.fluids import FluidNetwork, fluid_model, flows
network = FluidNetwork.example()
m = fluid_model(network)
assert np.allclose(np.sum(m.a, axis=0), 0)
assert abs(m.b.sum()) < 1e-12
initial_volume = m.mean.sum()
m.advance(10)
assert np.isclose(m.mean.sum(), initial_volume, atol=1e-10)
assert np.linalg.eigvalsh(m.cov).min() >= -1e-8
assert len(flows(m, network)['mL_per_s']) == len(network.edges)
print('verified: vascular/interstitial/lymphatic network conserves fluid volume')
