"""independent partitions may overlap; membership within a system sums to <= 1."""
from dataclasses import dataclass
from types import MappingProxyType


@dataclass(frozen=True)
class Partition:
    id: str
    memberships: dict[str, dict[str, float]]

    def __post_init__(self):
        object.__setattr__(self, 'memberships', MappingProxyType({
            region: MappingProxyType(dict(weights)) for region, weights in self.memberships.items()}))
