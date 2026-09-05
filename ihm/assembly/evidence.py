"""Evidence ancestry with explicit unknown uncertainty and parameter conditions."""
from dataclasses import dataclass
from .contracts import identity, finite, dimension


@dataclass(frozen=True)
class Uncertainty:
    status: str
    variance: float | None = None
    unit: str | None = None
    reason: str | None = None

    def __post_init__(self):
        if self.status == 'unknown':
            if self.variance is not None: raise ValueError('Unknown uncertainty cannot have a variance')
            identity(self.reason)
        elif self.status == 'quantified':
            finite(self.variance)
            if self.variance < 0: raise ValueError('Negative variance')
            dimension(self.unit)
        else: raise ValueError('Uncertainty must be unknown or quantified')


@dataclass(frozen=True)
class EvidenceNode:
    id: str
    kind: str
    dependencies: tuple[str, ...]
    uncertainty: Uncertainty

    def __post_init__(self):
        identity(self.id, *self.dependencies)
        object.__setattr__(self, 'dependencies', tuple(self.dependencies))
        if self.kind not in {'source','repair','synthesis','fit','reduction'}: raise ValueError('Unknown evidence kind')
        if not isinstance(self.uncertainty, Uncertainty): raise ValueError('Explicit uncertainty required')


@dataclass(frozen=True)
class ParameterCard:
    id: str
    value: float
    unit: str
    conditions: str
    evidence_ids: tuple[str, ...]
    identifiability: str
    fit_data_ids: tuple[str, ...]
    holdout_data_ids: tuple[str, ...]
    uncertainty: Uncertainty

    def __post_init__(self):
        identity(self.id, self.conditions, self.identifiability)
        finite(self.value); dimension(self.unit)
        if not self.evidence_ids: raise ValueError('Parameter evidence required')
        for name in ('evidence_ids', 'fit_data_ids', 'holdout_data_ids'):
            object.__setattr__(self, name, tuple(getattr(self, name)))
            identity(*getattr(self, name))
        if set(self.fit_data_ids) & set(self.holdout_data_ids): raise ValueError('Fit and holdout data overlap')
        if not isinstance(self.uncertainty, Uncertainty): raise ValueError('Explicit uncertainty required')


class EvidenceGraph:
    def __init__(self): self.nodes = {}

    def add(self, node: EvidenceNode):
        if node.id in self.nodes: raise ValueError('Duplicate evidence identity')
        # Topological insertion rejects cycles and unresolved ancestry atomically.
        if node.id in node.dependencies or set(node.dependencies)-self.nodes.keys(): raise ValueError('Cyclic or missing evidence dependency')
        self.nodes[node.id] = node

    def ancestors(self, node_id):
        if node_id not in self.nodes: raise ValueError('Unknown evidence identity')
        result, pending = set(), list(self.nodes[node_id].dependencies)
        while pending:
            key = pending.pop()
            if key not in result: result.add(key); pending.extend(self.nodes[key].dependencies)
        return result

    def uncertainty(self, node_id):
        ancestry = self.ancestors(node_id)
        if any(self.nodes[k].uncertainty.status == 'unknown' for k in ancestry | {node_id}):
            return Uncertainty('unknown', reason='Unquantified uncertainty in evidence ancestry')
        if ancestry:
            return Uncertainty('unknown', reason='No joint covariance and propagation model supplied for dependency graph')
        return self.nodes[node_id].uncertainty

    def validate_parameter(self, card: ParameterCard):
        if set(card.evidence_ids + card.fit_data_ids + card.holdout_data_ids)-self.nodes.keys(): raise ValueError('Unknown parameter evidence')
        if card.uncertainty.status == 'quantified':
            if card.uncertainty.unit != card.unit: raise ValueError('Parameter uncertainty unit mismatch')
            if any(self.uncertainty(k).status == 'unknown' for k in card.evidence_ids+card.fit_data_ids): raise ValueError('Quantified parameter depends on unknown uncertainty')
