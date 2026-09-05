"""allowed compartment interactions; strengths belong to processes."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Topology:
    id: str
    edges: tuple[tuple[str, str], ...]
