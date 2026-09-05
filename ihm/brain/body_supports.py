"""Explicit discrete sampling of canonical body sites, without relabeling head geometry."""
import numpy as np


def body_surface_geometry(sites_m):
    from ibm.materialize import DiscreteGeometry, GeometrySet
    sites = np.asarray(sites_m, dtype=float)
    if sites.ndim != 2 or sites.shape[1] != 3 or not len(sites) or not np.isfinite(sites).all():
        raise ValueError('sites_m must be a nonempty finite (sites, 3) coordinate array')
    return GeometrySet({'body_surface': DiscreteGeometry(
        'body_surface', 'ihm_canonical_body', sites * 1000,
        ids=tuple(f'ihm-body-site-{i}' for i in range(len(sites))),
        source='IHM canonical body sites, metres converted to millimetres; discrete sampling')},
        subject='ihm_canonical_body')


SUPPORT_GAPS = (
    'head_volume mechanical/thermal carriers are not whole-body supports',
    'body_surface uses explicit canonical point samples, not a reconstructed receptor census',
    'afferent/efferent tract geometry and target circuits are not supplied by this adapter',
    'thermal, nociceptive, chemical and visceral execution is not validated here',
)
