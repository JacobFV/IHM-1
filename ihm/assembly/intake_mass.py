"""Validate native consumed-intake boundaries without owning body mass.

Queued actions and internal GI redistribution never earn mass credit. A caller
must retain the returned boundary identity with its mechanical transaction;
this pure observer neither applies nor retries a mechanical impulse.
"""
from copy import deepcopy
import math

COMPONENTS = ('carbohydrate_g','protein_g','fat_g','sodium_g','calcium_mg','water_ml')
# Retained SENutrition::GetWeight convention: water mL contributes grams.
KG_FACTORS = (1e-3,1e-3,1e-3,1e-3,1e-6,1e-3)


def _number(value, label):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
        raise ValueError('Invalid consumed intake '+label)
    return value


def _integer(value,label,minimum=0):
    if type(value) is not int or value<minimum:
        raise ValueError('Invalid consumed intake '+label)
    return value


def _equal(a,b,label):
    if not math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12):
        raise ValueError('Inconsistent consumed intake '+label)


def _components(value,label):
    if not isinstance(value,dict) or set(value)!=set(COMPONENTS):
        raise ValueError('Invalid consumed intake '+label)
    return {k:_number(value[k],label+'.'+k) for k in COMPONENTS}


def _mass(parts):
    return math.fsum(parts[k]*factor for k,factor in zip(COMPONENTS,KG_FACTORS))


def _envelope(value):
    if not isinstance(value,dict) or value.get('schema')!='ihm.native-consumed-intake.v1':
        raise ValueError('Native consumed intake receipt required')
    epoch=value.get('owner_epoch')
    if not isinstance(epoch,str) or not epoch or len(epoch)>256:
        raise ValueError('Native intake owner epoch required')
    count=_integer(value.get('consumed_count'),'count')
    total=_number(value.get('cumulative_mass_kg'),'cumulative mass')
    parts=_components(value.get('cumulative_native'),'cumulative components')
    _equal(total,_mass(parts),'native mass convention')
    last=value.get('last_consumed')
    if count==0:
        if last is not None or total!=0 or any(parts.values()):raise ValueError('Nonempty initial intake receipt')
    else:
        if not isinstance(last,dict) or _integer(last.get('consumed_count'),'last count',1)!=count:
            raise ValueError('Missing latest native consumption')
        meal=_integer(last.get('meal_sequence'),'meal sequence',1)
        advance=_integer(last.get('advance_sequence'),'advance sequence',1)
        if advance<=meal:raise ValueError('Consumption must follow enqueue')
        start=_integer(last.get('interval_start_tick'),'start tick')
        end=_integer(last.get('interval_end_tick'),'end tick')
        if end<=start:raise ValueError('Consumption requires an advancing interval')
        t0=_number(last.get('native_start_s'),'native start')
        t1=_number(last.get('native_end_s'),'native end')
        if t1<=t0 or not math.isclose(t1-t0,(end-start)*.02,rel_tol=1e-8,abs_tol=1e-7):
            raise ValueError('Native intake clock disagreement')
        payload=last.get('native_payload')
        if not isinstance(payload,dict) or not isinstance(payload.get('name'),str):raise ValueError('Native payload required')
        composition=_components({k:payload.get(k) for k in COMPONENTS},'payload')
        _equal(_number(payload.get('mass_kg'),'payload mass'),_mass(composition),'payload mass convention')
    return epoch,count,total,parts,last


def consumed_intake_delta(previous,current):
    """Return one verified external mass boundary, or None for no consumption.

    Missing observations, skipped consumption counts, owner changes and altered
    counters reject. A fresh consumer must initialize from the native zero-count
    envelope rather than assume historical meals have not already been applied.
    """
    old_epoch,old_count,old_mass,old_parts,old_last=_envelope(previous)
    epoch,count,mass,parts,last=_envelope(current)
    if epoch!=old_epoch:raise ValueError('Native intake owner changed')
    if count==old_count:
        if mass!=old_mass or parts!=old_parts or last!=old_last:
            raise ValueError('Native intake changed without consumption')
        return None
    if count!=old_count+1:raise ValueError('Skipped or replayed native consumption')
    if old_last is not None:
        if last['meal_sequence']<=old_last['advance_sequence'] or last['interval_start_tick']<old_last['interval_end_tick']:
            raise ValueError('Native intake sequence or clock went backwards')
        gap_ticks=last['interval_start_tick']-old_last['interval_end_tick']
        if not math.isclose(last['native_start_s']-old_last['native_end_s'],gap_ticks*.02,rel_tol=1e-8,abs_tol=1e-7):
            raise ValueError('Native intake clock epoch changed')
    payload=last['native_payload']
    _equal(mass-old_mass,payload['mass_kg'],'mass boundary increment')
    for key in COMPONENTS:
        _equal(parts[key]-old_parts[key],payload[key],'component boundary increment '+key)
    return {'schema':'ihm.consumed-intake-mass-boundary.v1',
        'boundary_id':[epoch,count,last['meal_sequence']],
        'mass_kg':payload['mass_kg'],'native_payload':deepcopy(payload),
        'native_advance_sequence':last['advance_sequence'],
        'interval_start_tick':last['interval_start_tick'],'interval_end_tick':last['interval_end_tick'],
        'mechanical_transfer_applied':False,
        'basis':'Observed native consumption; mass uses retained SENutrition convention'}
