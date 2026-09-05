"""Unit conversions at the body/IBM boundary; mechanics is never an implicit gain."""
import numpy as np


def pressure_to_indentation_um(pressure_pa, *, stiffness_pa_per_m):
    """Linear foundation assumption p=k*d, supplied by the caller, not IBM.

    k has units Pa/m (not a Young modulus). This local normal indentation model
    does not establish skin thickness, contact area, or constitutive calibration.
    """
    stiffness = np.asarray(stiffness_pa_per_m, dtype=float)
    pressure = np.asarray(pressure_pa, dtype=float)
    if not np.all(np.isfinite(stiffness)) or np.any(stiffness <= 0):
        raise ValueError('stiffness_pa_per_m must be finite and positive')
    if not np.all(np.isfinite(pressure)):
        raise ValueError('pressure must be finite')
    return pressure / stiffness * 1e6


def displacement_mm_to_indentation_um(displacement_mm):
    """Donor mechanical.displacement is mm; receptor transfer accepts um."""
    return np.asarray(displacement_mm, dtype=float) * 1000
