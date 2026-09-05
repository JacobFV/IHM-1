"""Conditioned human skin/cotton measurements, not a universal friction law.

Temel, Johnson & Lloyd (2022), doi:10.1007/s11249-021-01560-5.
Means and SD are separately retained; SD is never used as a coefficient,
confidence interval, or independent static/dynamic sampling distribution.
"""
from copy import deepcopy
import json
from pathlib import Path

_DATA = json.loads((Path(__file__).resolve().parents[2] /
                    'data/measurements/textile/temel_2022.json').read_text())


def measurement_conditions():
    """Return a fresh description of the measured protocol, in SI where possible.

    This identifies an evidence cohort. It does not assert that arbitrary
    garment contacts satisfy its load, speed, hydration or temperature regime.
    """
    return deepcopy(_DATA['conditions'])


def textile_friction(region, *, fabric, conditions):
    """Return reported dimensionless mean/SD for an exactly named evidence case.

    Only chest and dorsal_forearm have numerically transcribed coefficients.
    The required complete conditions prevent implicit ambient/wet defaults.
    No interpolation, anatomical reassignment or fabric substitution is made.
    """
    if not isinstance(region,str) or region not in _DATA['regions']:
        raise ValueError('unsupported region: numeric evidence covers chest and dorsal_forearm')
    if fabric != _DATA['fabric_id']:
        raise ValueError('unsupported fabric: use temel_2022_cotton_single_jersey')
    if not isinstance(conditions,dict) or conditions != _DATA['conditions']:
        raise ValueError('unsupported conditions: require the complete measurement_conditions() protocol')
    result=deepcopy(_DATA['regions'][region])
    result.update(region=region,fabric=fabric,conditions=measurement_conditions(),
                  participants_n=_DATA['cohort']['participants_n'],cohort=deepcopy(_DATA['cohort']),
                  fabric_properties=deepcopy(_DATA['fabric']),doi=_DATA['doi'],
                  estimate_kind='reported_cohort_mean_and_sd',unit='dimensionless',
                  limitations=deepcopy(_DATA['limitations']))
    return result
