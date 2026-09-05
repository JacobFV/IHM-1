"""IHM-1: one implicit substrate, multiple explicit predictor views."""
from ihm.body import body
from ihm.materialize import Request, materialize

from ihm.human import ImplicitHuman

__all__ = ['body', 'Request', 'materialize', 'ImplicitHuman']
