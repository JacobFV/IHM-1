"""Versioned material identities and SI exchange contracts.

These records specify software obligations, not scientifically validated laws.
Units are an intentionally small explicit vocabulary; conversion is external.
"""
from dataclasses import dataclass
from typing import Any, Mapping
import math

UNITS = {'1':'dimensionless', 'm':'length', 'm2':'area', 'm3':'volume',
         's':'time', 'Hz':'frequency', 'kg':'mass', 'mol':'species',
         'C':'charge', 'kg*m/s':'momentum', 'J':'energy', 'N':'force',
         'Pa':'pressure', 'K':'temperature', 'V':'voltage', 'A':'current',
         'kg/s':'mass_rate', 'mol/s':'species_rate', 'm3/s':'volume_rate', 'W':'power'}
QUANTITIES = {'mass':'kg', 'species':'mol', 'charge':'C', 'momentum':'kg*m/s', 'energy':'J', 'volume':'m3'}


def identity(*values):
    if any(not isinstance(v, str) or not v.strip() for v in values):
        raise ValueError('Nonempty explicit identities required')


def finite(value, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or (positive and value <= 0):
        raise ValueError('Finite numeric value required' + (' and must be positive' if positive else ''))


def dimension(unit):
    if unit not in UNITS: raise ValueError('Unregistered unit: ' + str(unit))
    return UNITS[unit]


@dataclass(frozen=True)
class MaterialPoint:
    body_version: str
    entity_id: str
    element_id: str
    local_coordinates: tuple[float, ...]

    def __post_init__(self):
        identity(self.body_version, self.entity_id, self.element_id)
        object.__setattr__(self, 'local_coordinates', tuple(self.local_coordinates))
        if not 1 <= len(self.local_coordinates) <= 4: raise ValueError('Expected 1..4 reference coordinates')
        for x in self.local_coordinates: finite(x)


@dataclass(frozen=True)
class StateSpec:
    id: str
    owner: str
    support: str
    unit: str
    conserved_quantity_or_none: str | None
    spatial_discretization: str
    temporal_representation: str
    validity_domain: str

    def __post_init__(self):
        identity(self.id, self.owner, self.support, self.spatial_discretization, self.temporal_representation, self.validity_domain)
        dim = dimension(self.unit)
        if self.conserved_quantity_or_none is not None and (self.conserved_quantity_or_none not in QUANTITIES or dim != self.conserved_quantity_or_none):
            raise ValueError('Conserved state dimension mismatch')


@dataclass(frozen=True)
class Interaction:
    id: str
    endpoint_a: str
    endpoint_b: str
    kind: str
    orientation: str
    law_id: str | None
    parameters_id: str | None
    active_domain: str
    evidence_ids: tuple[str, ...] = ()
    unresolved_reason: str | None = None

    def __post_init__(self):
        identity(self.id, self.endpoint_a, self.endpoint_b, self.orientation, self.active_domain)
        if self.kind not in {'joint','attachment','contact','sliding','flow','diffusion','neural_afferent','neural_efferent','receptor','gap_junction','thermal'}:
            raise ValueError('Unknown interaction kind')
        object.__setattr__(self, 'evidence_ids', tuple(self.evidence_ids))
        identity(*self.evidence_ids)
        if not self.unresolved_reason: identity(self.law_id, self.parameters_id)


@dataclass(frozen=True)
class Exchange:
    interface_id: str
    interval_start_s: float
    interval_end_s: float
    quantity: str
    amount_SI: float
    source_owner: str
    target_owner: str
    work_J_or_none: float | None = None

    def __post_init__(self):
        identity(self.interface_id, self.source_owner, self.target_owner)
        for v in (self.interval_start_s, self.interval_end_s, self.amount_SI): finite(v)
        if self.interval_end_s <= self.interval_start_s: raise ValueError('Positive exchange interval required')
        if self.source_owner == self.target_owner: raise ValueError('Exchange requires distinct owners')
        if self.quantity not in QUANTITIES: raise ValueError('Unknown conserved exchange dimension')
        if self.work_J_or_none is not None: finite(self.work_J_or_none)

    def signed_entries(self):
        """Bookkeeping debits/credits; adapters must independently audit storage."""
        work = self.work_J_or_none or 0.
        return ((self.source_owner, -self.amount_SI, -work), (self.target_owner, self.amount_SI, work))


@dataclass(frozen=True)
class RateSample:
    interface_id: str
    time_s: float
    quantity: str
    rate_SI_per_s: float

    def __post_init__(self):
        identity(self.interface_id)
        finite(self.time_s); finite(self.rate_SI_per_s)
        if self.quantity not in QUANTITIES: raise ValueError('Unknown rate dimension')


@dataclass(frozen=True)
class ObservationRequest:
    body_version: str
    region: str
    observables: tuple[str, ...]
    physical_resolution: float
    temporal_bandwidth_Hz: float
    error_target: float
    compute_budget: Mapping[str, float]
    intervention: Mapping[str, Any] | None = None

    def __post_init__(self):
        identity(self.body_version, self.region)
        object.__setattr__(self, 'observables', tuple(self.observables))
        if not self.observables: raise ValueError('Observables required')
        identity(*self.observables)
        for v in (self.physical_resolution, self.temporal_bandwidth_Hz, self.error_target): finite(v, positive=True)
        if not self.compute_budget: raise ValueError('Explicit compute budget required')
        for k,v in self.compute_budget.items(): identity(k); finite(v, positive=True)


@dataclass(frozen=True)
class MaterializationReceipt:
    id: str
    data_hashes: Mapping[str, str]
    code_hashes: Mapping[str, str]
    solver_versions: Mapping[str, str]
    configuration: Mapping[str, Any]
    seed: int | None
    state_owners: Mapping[str, str]
    transfers: tuple[str, ...]
    error_estimates: Mapping[str, Any]
    uncertainty: Mapping[str, Any]
    validity_domain: str
    playback_kind: str = 'cosimulation'

    def __post_init__(self):
        identity(self.id, self.validity_domain)
        if self.playback_kind not in {'cosimulation','recorded_source_replay'}: raise ValueError('Explicit execution semantics required')
        for mapping in (self.data_hashes, self.code_hashes, self.solver_versions, self.state_owners):
            for k,v in mapping.items(): identity(k,v)
        if self.seed is not None and (isinstance(self.seed, bool) or not isinstance(self.seed, int)): raise ValueError('Integer seed or None required')


@dataclass(frozen=True)
class Frame:
    model_time_s: float
    receipt_id: str
    transforms: Mapping[str, Any]
    deformation_fields: Mapping[str, Any]
    physiological_fields: Mapping[str, Any]
    neural_fields: Mapping[str, Any]
    observations: Mapping[str, Any]
    audit: Mapping[str, Any]
    playback_kind: str

    def __post_init__(self):
        finite(self.model_time_s); identity(self.receipt_id)
        if self.playback_kind not in {'cosimulation','recorded_source_replay'}: raise ValueError('Explicit frame execution semantics required')
