"""Fail-closed registration before a coupled scenario can run."""
from dataclasses import dataclass
from typing import Protocol, Any
from .contracts import StateSpec, Interaction, MaterialPoint, Exchange, ObservationRequest, Frame, identity, dimension


@dataclass(frozen=True)
class Port:
    id: str
    interface_id: str
    owner: str
    state_id: str
    quantity: str
    unit: str

    def __post_init__(self):
        identity(self.id, self.interface_id, self.owner, self.state_id, self.quantity)
        if dimension(self.unit) != self.quantity: raise ValueError('Port dimension mismatch')


class EngineAdapter(Protocol):
    def checkpoint(self) -> Any: ...
    def restore(self, checkpoint: Any) -> None: ...
    def advance(self, start_s: float, end_s: float, incoming: tuple[Exchange, ...], event_queue: list) -> Any: ...


class ExplicitBodyModel(Protocol):
    """Execution boundary for future materializers; no implicit source replay."""
    def initialize(self, reference_conditions: Any) -> Any: ...
    def advance(self, checkpoint: Any, end_time_s: float, interventions: Any) -> tuple[Any, Any, Any]: ...
    def observe(self, checkpoint: Any, material_points: tuple[MaterialPoint, ...], quantities: tuple[str, ...]) -> Any: ...
    def project_for_display(self, checkpoint: Any, view_request: Any) -> Frame: ...


class Materializer(Protocol):
    def materialize(self, request: ObservationRequest) -> ExplicitBodyModel: ...


class ContractRegistry:
    def __init__(self, body_version, entity_ids):
        identity(body_version, *entity_ids)
        self.body_version, self.entity_ids = body_version, frozenset(entity_ids)
        self.states, self.interactions, self.ports = {}, {}, {}
        self._validated_signature = None

    def validate_point(self, point: MaterialPoint):
        if point.body_version != self.body_version or point.entity_id not in self.entity_ids:
            raise ValueError('Unknown body version or material entity')

    def register_state(self, state: StateSpec):
        if state.id in self.states: raise ValueError('Duplicate state or dual state ownership: ' + state.id)
        if state.support not in self.entity_ids: raise ValueError('Unknown state support')
        self.states[state.id] = state

    def register_interaction(self, interaction: Interaction):
        if interaction.id in self.interactions: raise ValueError('Duplicate interface identity')
        if {interaction.endpoint_a, interaction.endpoint_b} - self.entity_ids: raise ValueError('Unknown interaction endpoint')
        self.interactions[interaction.id] = interaction

    def register_port(self, port: Port):
        if port.id in self.ports: raise ValueError('Duplicate port identity')
        state = self.states.get(port.state_id)
        interface = self.interactions.get(port.interface_id)
        if state is None or interface is None: raise ValueError('Unregistered state or interface')
        if state.owner != port.owner or dimension(state.unit) != dimension(port.unit): raise ValueError('Port ownership or dimension mismatch')
        if state.support not in {interface.endpoint_a, interface.endpoint_b}: raise ValueError('Port state outside interface')
        self.ports[port.id] = port

    def _signature(self):
        return (tuple(self.states.items()), tuple(self.interactions.items()), tuple(self.ports.items()))

    def validate_scenario(self, owners, law_ids, parameter_ids, evidence_ids):
        self._validated_signature = None
        if not self.states or set(owners) != {s.owner for s in self.states.values()}: raise ValueError('Scenario must register exactly its state owners')
        for i in self.interactions.values():
            if i.unresolved_reason or i.law_id not in law_ids or i.parameters_id not in parameter_ids or not i.evidence_ids or set(i.evidence_ids)-set(evidence_ids):
                raise ValueError('Unresolved interface dependencies: ' + i.id)
            ports = [p for p in self.ports.values() if p.interface_id == i.id]
            if len({p.owner for p in ports}) < 2: raise ValueError('Interface needs ports for both owners')
            for p in ports:
                if not any(q.owner != p.owner and q.quantity == p.quantity and self.states[q.state_id].support != self.states[p.state_id].support for q in ports): raise ValueError('Unpaired port dimension or anatomical endpoint')
        self._validated_signature = self._signature()

    def require_ready(self, owners):
        if self._validated_signature != self._signature() or self._validated_signature is None:
            raise ValueError('Validate scenario dependencies before enabling integration')
        if set(owners) != {s.owner for s in self.states.values()}: raise ValueError('Engine/state ownership mismatch')

    def validate_exchange(self, exchange: Exchange):
        if not isinstance(exchange, Exchange): raise ValueError('Integrated Exchange required; rate samples are not amounts')
        for owner in (exchange.source_owner, exchange.target_owner):
            if not any(p.interface_id == exchange.interface_id and p.owner == owner and p.quantity == exchange.quantity for p in self.ports.values()):
                raise ValueError('Exchange has unregistered owner/interface/dimension')
        source_ports=[p for p in self.ports.values() if p.interface_id==exchange.interface_id and p.owner==exchange.source_owner and p.quantity==exchange.quantity]
        target_ports=[p for p in self.ports.values() if p.interface_id==exchange.interface_id and p.owner==exchange.target_owner and p.quantity==exchange.quantity]
        if not any(self.states[p.state_id].support!=self.states[q.state_id].support for p in source_ports for q in target_ports):
            raise ValueError('Exchange must connect opposite anatomical endpoints')
