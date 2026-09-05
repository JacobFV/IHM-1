"""affine pressure on physical state: dx_out/dt += bias + weights @ x_in."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Process:
    id: str
    inputs: tuple[str, ...]
    output: str
    topology: str
    weights: tuple[float, ...]
    bias: float
    noise: float
    provenance: str = 'illustrative_weak_prior'
    form: str = 'affine'
