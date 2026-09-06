"""Admission contract for a future native muscle supply preview.

This module neither estimates substrate availability nor limits force. The held
native adapter does not yet produce SupplyPreview. A missing preview fails
closed; snapshots and previous-step unmet demand are not supply certificates.
All energies are joules under the native chemical-demand convention, not ATP.
"""
from dataclasses import dataclass
import math
import sys

@dataclass(frozen=True)
class TrialBinding:
    physiology_state_id: str
    mechanics_state_id: str
    control_id: str
    reference_id: str
    native_build_id: str
    sequence: int
    start_s: float
    end_s: float

@dataclass(frozen=True)
class MuscleDemandTrial:
    binding: TrialBinding
    delta_metabolic_j: float
    delta_heat_j: float
    delta_signed_work_j: float

@dataclass(frozen=True)
class SupplyPreview:
    binding: TrialBinding
    delta_metabolic_j: float
    native_base_muscle_j: float
    requested_muscle_j: float
    provided_muscle_j: float
    unmet_muscle_j: float
    # Must describe a candidate native replay from the bound pre-step state.
    # The producer must restore the complete native transaction after preview.
    origin: str = 'native-candidate-replay'

@dataclass(frozen=True)
class SupplyDecision:
    admissible: bool
    reason: str
    unmet_muscle_j: float | None = None

def _finite(*values):
    if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in values):
        raise ValueError('Finite numeric energy/time required')

def _same(a,b):
    return abs(a-b)<=128*sys.float_info.epsilon*max(1.,abs(a),abs(b))

def assess_supply(trial: MuscleDemandTrial, preview: SupplyPreview | None) -> SupplyDecision:
    """Validate an uncommitted trial; no state, force, chemistry or heat writes.

    Insufficient supply requires a new mechanically integrated candidate from
    the original joint checkpoint. It never authorizes editing performed work.
    """
    b=trial.binding
    if any(not isinstance(v,str) or not v for v in (b.physiology_state_id,b.mechanics_state_id,b.control_id,b.reference_id,b.native_build_id)):
        raise ValueError('Full native state/control/reference/build binding required')
    if isinstance(b.sequence,bool) or not isinstance(b.sequence,int) or b.sequence<=0:
        raise ValueError('Positive exchange sequence required')
    _finite(b.start_s,b.end_s,trial.delta_metabolic_j,trial.delta_heat_j,trial.delta_signed_work_j)
    if b.end_s<=b.start_s:raise ValueError('Positive trial interval required')
    if not _same(trial.delta_metabolic_j,trial.delta_heat_j+trial.delta_signed_work_j):
        raise ValueError('Signed chemical/heat/work incidence mismatch')
    if preview is None:return SupplyDecision(False,'native_supply_preview_unavailable')
    if preview.binding!=b:raise ValueError('Stale or mismatched native supply preview')
    if preview.origin!='native-candidate-replay':raise ValueError('Native candidate replay required')
    _finite(preview.delta_metabolic_j,preview.native_base_muscle_j,preview.requested_muscle_j,preview.provided_muscle_j,preview.unmet_muscle_j)
    if any(v<0 for v in (preview.native_base_muscle_j,preview.requested_muscle_j,preview.provided_muscle_j,preview.unmet_muscle_j)):
        raise ValueError('Absolute muscle supply budgets must be nonnegative')
    if not _same(preview.delta_metabolic_j,trial.delta_metabolic_j):raise ValueError('Preview demand differs from candidate')
    if not _same(preview.requested_muscle_j,preview.native_base_muscle_j+trial.delta_metabolic_j):
        raise ValueError('Native base plus signed increment differs from request')
    if not _same(preview.requested_muscle_j,preview.provided_muscle_j+preview.unmet_muscle_j):
        raise ValueError('Native provided/unmet chemical ledger mismatch')
    if not _same(preview.unmet_muscle_j,0.):
        return SupplyDecision(False,'reintegrate_candidate_from_joint_checkpoint',preview.unmet_muscle_j)
    return SupplyDecision(True,'chemical_supply_admissible_only',preview.unmet_muscle_j)
