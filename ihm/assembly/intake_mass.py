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


def planned_nutrition_mass_kg(meal):
    """Capacity screening only; configured nutrition never earns mass credit."""
    return _mass(_components({k:getattr(meal,k,None) for k in COMPONENTS},'planned nutrition'))


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


class IntakeMassBridge:
    """Opt-in endpoint transfer from a fresh consumed-intake owner to one plant.

    The declared ingestion velocity is a co-moving assumption at the registered
    station. This is not measured swallowing momentum. Capacity is the native
    port's validated 0.5 kg total payload; no implicit splitting or excretion.
    """
    def __init__(self,plant,initial_intake,*,body,station_m,registration_identity,incoming_velocity_basis):
        import hashlib,json
        epoch,count,_,_,_=_envelope(initial_intake)
        if count!=0:raise ValueError('Intake bridge requires fresh zero-count native owner')
        if body!='torso':raise ValueError('Only validated torso intake target supported')
        if not isinstance(registration_identity,str) or not registration_identity:
            raise ValueError('Explicit intake station registration identity required')
        if incoming_velocity_basis!='co_moving_at_ingestion_assumption':
            raise ValueError('Explicit co-moving ingestion assumption required')
        self.station=self._vector(station_m)
        self.plant=plant;self.body=body;self.registration_identity=registration_identity
        self.incoming_velocity_basis=incoming_velocity_basis
        self.previous=deepcopy(initial_intake);self.failed=False;self.pending=None;self.last=None
        state=plant.snapshot();mass=state.get('mass_transfer',{})
        if mass.get('enabled') is not True or mass.get('last_sequence')!=0 or mass.get('owned_payload_mass_kg')!=0:
            raise ValueError('Fresh enabled mechanical mass owner required')
        self.reference=mass.get('reference_id')
        if not isinstance(self.reference,str) or not self.reference:raise ValueError('Mechanical mass reference required')
        self.owner='intake-'+hashlib.sha256(json.dumps([epoch,self.reference,body,self.station,registration_identity]).encode()).hexdigest()
        self.sequence=0;self.applied_mass=0.

    @staticmethod
    def _vector(value):
        if not isinstance(value,(list,tuple)) or len(value)!=3 or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in value):
            raise ValueError('Finite source-frame vector required')
        return [float(v) for v in value]

    def apply(self,current):
        if self.failed:raise RuntimeError('Intake mass bridge failed; no replay permitted')
        boundary=consumed_intake_delta(self.previous,current)
        if boundary is None:
            self.previous=deepcopy(current)
            return None
        self.pending=deepcopy(boundary)
        try:return self._apply_boundary(current,boundary)
        except BaseException:
            self.failed=True
            raise

    def _apply_boundary(self,current,boundary):
        state=self.plant.snapshot();mass=state.get('mass_transfer',{})
        if mass.get('enabled') is not True or mass.get('reference_id')!=self.reference or mass.get('last_sequence')!=self.sequence:
            raise ValueError('Mechanical intake owner or sequence changed')
        _equal(_number(mass.get('owned_payload_mass_kg'),'mechanical payload'),self.applied_mass,'mechanical payload owner')
        if self.applied_mass+boundary['mass_kg']>.5:
            raise ValueError('Consumed intake exceeds validated 0.5 kg mechanical payload capacity')
        if not math.isclose(_number(state.get('time_s'),'mechanical time'),boundary['interval_end_tick']*.02,rel_tol=0,abs_tol=1e-8):
            raise ValueError('Mechanical intake endpoint clock mismatch')
        point=self.plant.body_point(body=self.body,station_m=self.station)
        if point.get('kind')!='body_point' or point.get('body')!=self.body or self._vector(point.get('station_m'))!=self.station:
            raise ValueError('Mechanical point query site changed')
        if point.get('time_s')!=state['time_s']:raise ValueError('Mechanical point query clock changed')
        velocity=self._vector(point.get('velocity_source_m_s'))
        if boundary['mass_kg']>0:
            result=self.plant.transfer_mass(sequence=self.sequence+1,owner=self.owner,body=self.body,
                delta_mass_kg=boundary['mass_kg'],station_m=self.station,velocity_source_m_s=velocity)
            observed=result.get('mass_transfer',{});receipt=observed.get('last_receipt',{})
            if observed.get('enabled') is not True or observed.get('reference_id')!=self.reference or observed.get('last_sequence')!=self.sequence+1:
                raise ValueError('Mechanical intake acknowledgment identity mismatch')
            if receipt.get('sequence')!=self.sequence+1 or receipt.get('owner')!=self.owner or receipt.get('body')!=self.body:
                raise ValueError('Mechanical intake acknowledgment owner mismatch')
            owners=observed.get('owners',{});owner=owners.get(self.owner,{}) if isinstance(owners,dict) else {}
            if owner.get('body')!=self.body or self._vector(owner.get('station_m'))!=self.station:
                raise ValueError('Mechanical intake acknowledgment site mismatch')
            _equal(_number(owner.get('mass_kg'),'owner inventory'),self.applied_mass+boundary['mass_kg'],'owner inventory')
            _equal(_number(receipt.get('delta_mass_kg'),'acknowledged mass'),boundary['mass_kg'],'acknowledged boundary')
            _equal(_number(observed.get('owned_payload_mass_kg'),'acknowledged payload'),self.applied_mass+boundary['mass_kg'],'acknowledged payload')
            if result.get('time_s')!=state['time_s']:raise ValueError('Mass transfer unexpectedly advanced time')
            self.sequence+=1;self.applied_mass+=boundary['mass_kg']
        else:receipt=None
        report={**boundary,'mechanical_transfer_applied':boundary['mass_kg']>0,'boundary_accounted':True,'mechanical_reference_id':self.reference,
            'mechanical_owner':self.owner,'mechanical_receipt':deepcopy(receipt),
            'registration_identity':self.registration_identity,'station_source_m':list(self.station),
            'incoming_velocity_basis':self.incoming_velocity_basis,'incoming_velocity_source_m_s':velocity}
        self.previous=deepcopy(current);self.last=deepcopy(report);self.pending=None
        return report

    def snapshot(self):
        return {'schema':'ihm.intake-mass-bridge.v1','failed':self.failed,
            'mechanical_reference_id':self.reference,'mechanical_sequence':self.sequence,
            'applied_mass_kg':self.applied_mass,'last_boundary':deepcopy(self.last),
            'pending_boundary':deepcopy(self.pending),'capacity_kg':.5,
            'checkpoint_scope':'Observation only; cannot resume or rewind native consumption'}
