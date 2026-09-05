"""physical components on named supports; priors are epistemic beliefs."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Support:
    id: str
    frame: str
    regions: tuple[str, ...]


@dataclass(frozen=True)
class Component:
    id: str
    support: str
    region: str
    unit: str
    mean: float
    std: float
    provenance: str = 'illustrative_weak_prior'
