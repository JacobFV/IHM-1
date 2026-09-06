"""Native-owned receptor endpoint receipts and explicit inferred IBM allocation.

No receptor equations, baseline subtraction, extra brain step or efferent control.
"""
from copy import deepcopy
import hashlib,json,math

CHANNELS=('chemoreceptor','baroreceptor_carotid','baroreceptor_aortic','pulmonary_stretch')
IDENTITY_KEYS=('library_sha256','executable_sha256','state_sha256','manifest_sha256')

def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def finite(value,label,low=0.,high=1e9):
    if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value) or not low<=value<=high:raise ValueError('Invalid '+label)
    return float(value)

def identity(value):
    if not isinstance(value,dict) or set(value)!=set(IDENTITY_KEYS):raise ValueError('Exact native identity required')
    for x in value.values():
        if not isinstance(x,str) or len(x)!=64 or any(c not in '0123456789abcdef' for c in x):raise ValueError('SHA256 native identity required')
    return deepcopy(value)

def endpoint_receipt(snapshot,native_identity,source_sha256):
    """Bind actual observer fields to the accepted adapter endpoint, without defaults."""
    owner=identity(native_identity)
    if not isinstance(source_sha256,str) or len(source_sha256)!=64 or any(c not in '0123456789abcdef' for c in source_sha256):raise ValueError('Observer source SHA256 required')
    elapsed=finite(snapshot.get('elapsed_s'),'native elapsed time');origin=finite(snapshot.get('origin_s'),'native origin')
    native_time=finite(snapshot.get('time_s'),'native absolute time');tick=round(elapsed*50)
    sequence=snapshot.get('sequence')
    if type(sequence) is not int or sequence<0 or abs(elapsed-tick*.02)>1e-9 or abs(native_time-origin-elapsed)>1e-8:raise ValueError('Native afferent endpoint clock differs')
    values=snapshot.get('values',{})
    rates={name:finite(values.get('nervous.afferent.'+name+'_hz'),'native '+name+' Hz',0,1000) for name in CHANNELS}
    result={'schema':'ihm.native-afferent-endpoint.v1','identity':owner,'observer_source_sha256':source_sha256,
        'native_sequence':sequence,'endpoint_tick':tick,'elapsed_s':elapsed,'native_time_s':native_time,'origin_s':origin,'rates_hz':rates}
    result['sha256']=digest(result);return result

class NativeAfferentBridge:
    def __init__(self,brain,*,native_identity,source_sha256,origin_s,allocation):
        self.brain=brain
        self.native_identity=identity(native_identity);self.source_sha256=source_sha256;self.origin_s=finite(origin_s,'native origin')
        if not isinstance(allocation,dict) or set(allocation)!=set(CHANNELS):raise ValueError('Explicit allocation for every native afferent required')
        medulla={n['id'] for n in brain.data['nodes'] if n.get('autonomic_role')=='medulla'}
        self.allocation=deepcopy(allocation)
        for name,row in self.allocation.items():
            if not isinstance(row,dict) or not row or set(row)-medulla:raise ValueError('Allocation requires exact existing medulla IDs')
            if abs(math.fsum(finite(v,'allocation fraction',0,1) for v in row.values())-1)>1e-12:raise ValueError('Afferent allocation must sum to one')
        self.model_sha256=digest({'identity':self.native_identity,'observer_source_sha256':source_sha256,'origin_s':self.origin_s,
            'allocation':allocation,'brain_source_identity':brain.source_identity,'brain_data_sha256':digest(brain.data),
            'basis':'Absolute native firing allocated as additive IBM sensory Hz; inferred regional recruitment, no measured axon map'})
        self.time_s=0.;self.last_sequence=-1;self.last_receipt_sha256=None

    def step(self,dt_s,receipt,*,sensory_blocks=()):
        if finite(dt_s,'afferent exchange',.02,.02)!=.02:raise ValueError('Expected20ms exchange')
        if not isinstance(sensory_blocks,(list,tuple,set)) or any(not isinstance(k,str) for k in sensory_blocks) or set(sensory_blocks)-set(CHANNELS):raise ValueError('Unknown native sensory block')
        raw=deepcopy(receipt)
        if not isinstance(raw,dict):raise ValueError('Native afferent receipt required')
        claimed=raw.pop('sha256',None)
        if claimed!=digest(raw) or raw.get('schema')!='ihm.native-afferent-endpoint.v1':raise ValueError('Native afferent receipt hash differs')
        # Reuse the strict endpoint constructor to reject unexpected clock/rate shapes.
        if set(raw.get('rates_hz',{}))!=set(CHANNELS):raise ValueError('Missing native receptor channels')
        rebuilt=endpoint_receipt({'sequence':raw.get('native_sequence'),'elapsed_s':raw.get('elapsed_s'),'origin_s':raw.get('origin_s'),
            'time_s':raw.get('native_time_s'),'values':{'nervous.afferent.'+k+'_hz':v for k,v in raw['rates_hz'].items()}},raw.get('identity'),raw.get('observer_source_sha256'))
        if rebuilt!=receipt or raw['identity']!=self.native_identity or raw['observer_source_sha256']!=self.source_sha256 or raw['origin_s']!=self.origin_s:raise ValueError('Native afferent owner/source differs')
        if abs(raw['elapsed_s']-self.time_s)>1e-9 or raw['native_sequence']<=self.last_sequence:raise ValueError('Missing, stale or repeated native afferent endpoint')
        inputs={};ports=[]
        for name,rate in raw['rates_hz'].items():
            blocked=name in sensory_blocks
            for region,weight in self.allocation[name].items():
                contribution=0. if blocked else rate*weight
                inputs[region]=inputs.get(region,0.)+contribution
                ports.append({'channel':name,'region':region,'native_rate_hz':rate,'allocation_fraction':weight,'sensory_input_hz':contribution,'blocked':blocked})
        if any(v>1000 for v in inputs.values()):raise ValueError('Allocated native firing exceeds IBM sensory rate domain')
        start=self.time_s;self.time_s+=.02;self.last_sequence=raw['native_sequence'];self.last_receipt_sha256=claimed
        return {'schema':'ihm.native-afferent-input.v1','source_endpoint_s':start,'brain_interval_s':[start,self.time_s],
            'sensory_inputs_hz':inputs,'ports':ports,'receipt_sha256':claimed,'model_sha256':self.model_sha256,
            'basis':'Native receptor dynamics; inferred absolute-rate regional allocation; no efferent control'}

    def checkpoint(self):return {'model_sha256':self.model_sha256,'time_s':self.time_s,'last_sequence':self.last_sequence,'last_receipt_sha256':self.last_receipt_sha256}
    def restore(self,saved):
        if not isinstance(saved,dict) or set(saved)!=set(self.checkpoint()) or saved['model_sha256']!=self.model_sha256:raise ValueError('Afferent checkpoint identity differs')
        time=finite(saved['time_s'],'afferent checkpoint time');seq=saved['last_sequence'];sha=saved['last_receipt_sha256']
        if abs(time-round(time*50)*.02)>1e-9 or type(seq) is not int or seq< -1 or (time==0)!=(seq==-1) or (time==0)!=(sha is None):raise ValueError('Invalid afferent checkpoint clock')
        if sha is not None and (not isinstance(sha,str) or len(sha)!=64 or any(c not in '0123456789abcdef' for c in sha)):raise ValueError('Invalid afferent checkpoint receipt')
        self.time_s=time;self.last_sequence=seq;self.last_receipt_sha256=sha
