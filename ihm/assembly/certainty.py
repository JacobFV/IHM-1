"""Materialization uncertainty and dependence-aware scalar evidence fusion.

Provenance class is not a calibrated probability. Unknown covariance stays unknown.
"""
from dataclasses import dataclass,asdict
from collections import defaultdict,Counter
import math
import numpy as np

@dataclass(frozen=True)
class Estimate:
    source_id:str
    dependency_group:str
    mean:float
    variance:float
    unit:str
    basis:str
    def __post_init__(self):
        if not all(isinstance(s,str) and s for s in (self.source_id,self.dependency_group,self.unit,self.basis)):
            raise ValueError('Explicit source, dependency, unit and evidence basis required')
        if any(isinstance(v,bool) or not isinstance(v,(float,int)) or not math.isfinite(v) for v in (self.mean,self.variance)) or self.variance<=0:
            raise ValueError('Finite estimate and positive quantified variance required')

def fuse(estimates):
    """Equal-weight covariance intersection within ancestry; independent groups fuse.

Group independence is an explicit caller assertion. This operation cannot turn
unquantified registration discrepancy into measurement uncertainty.
"""
    estimates=list(estimates)
    if not estimates or any(not isinstance(e,Estimate) for e in estimates):raise ValueError('Expected estimates')
    if len({e.unit for e in estimates})!=1:raise ValueError('Convert units explicitly before fusion')
    if len({e.source_id for e in estimates})!=len(estimates):raise ValueError('Duplicate evidence identity')
    groups=defaultdict(list)
    for e in estimates:groups[e.dependency_group].append(e)
    group_results=[]
    for name,items in groups.items():
        smallest=min(e.variance for e in items)
        weights=np.array([smallest/e.variance for e in items]);total=float(weights.sum())
        mean=float(sum((w/total)*e.mean for w,e in zip(weights,items)))
        variance=smallest*len(items)/total
        group_results.append(dict(dependency_group=name,mean=mean,variance=variance,source_ids=[e.source_id for e in items]))
    smallest=min(g['variance'] for g in group_results)
    weights=np.array([smallest/g['variance'] for g in group_results]);total=float(weights.sum())
    mean=float(sum((w/total)*g['mean'] for w,g in zip(weights,group_results)));variance=smallest/total
    if not math.isfinite(mean) or not math.isfinite(variance) or variance<=0:raise FloatingPointError('Fusion not representable')
    return dict(mean=mean,variance=variance,unit=estimates[0].unit,groups=group_results,evidence=[asdict(e) for e in estimates],
                method='equal-weight covariance intersection within dependency groups; independent Gaussian product between groups',
                assumptions=['Distinct dependency groups have independent errors.','Input variances are justified externally; precision does not validate model form.'])

def entity_certainty(entity):
    result=dict(entity.get('uncertainty',{}))
    result['evidence_kind']=entity.get('evidence_kind','unclassified')
    result['assumptions']=list(entity.get('assumptions',[]))
    result['dependency_group']=entity.get('provenance',{}).get('dependency_group')
    bounds=entity.get('bounds_m',{});values=np.array([bounds.get('min',[0,0,0]),bounds.get('max',[0,0,0])],float)
    bound=float(np.max(np.abs(values)))
    if not np.isfinite(bound):raise ValueError('Invalid physical bounds')
    rounded=np.nextafter(np.float32(bound),np.float32(np.inf))
    result['display_numerics']={**result.get('display_numerics',{}),
        'gpu_coordinate_type':'float32','gpu_coordinate_rounding_bound_m':float(np.spacing(rounded)),
        'rounding_bound_interpretation':'Conservative componentwise numeric bound, not biological accuracy.'}
    return result

def summarize(anatomy):
    return dict(entities=len(anatomy['entities']),by_evidence_kind=dict(Counter(e.get('evidence_kind','unclassified') for e in anatomy['entities'])),
        assumption_count=len(anatomy.get('assumption_ledger',[])),calibrated_confidence_score=None,
        uncertainty_axes=['biological','registration_or_synthesis','display_numerics'],
        rule='Source, registration and synthesis status are not confidence percentages; only justified variances enter quantitative fusion.')
