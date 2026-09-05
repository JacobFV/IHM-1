"""Held body interventions with explicit times, identities and modalities."""
import copy
import math
from .peripheral import BodyPeripheral


def validate_interventions(events, peripheral_data, end_s):
    if not isinstance(events,(list,tuple)) or len(events)>100:
        raise ValueError('Expected at most 100 body interventions')
    result=copy.deepcopy(list(events))
    probe=BodyPeripheral.from_dict(peripheral_data)
    previous=-1.
    for event in result:
        if not isinstance(event,dict) or set(event)-{'start_s','end_s','stimuli','motor_commands','blocked_nerves'}:
            raise ValueError('Unknown body intervention field')
        for key in ['start_s','end_s']:
            value=event.get(key)
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
                raise ValueError('Body intervention requires finite start/end seconds')
        if not 0<=event['start_s']<event['end_s']<=end_s or event['start_s']<previous:
            raise ValueError('Body interventions must be ordered and within run duration')
        previous=event['start_s']
        probe.step(.001,event.get('stimuli',{}),brain_state={'motor_commands':event.get('motor_commands',{})},blocked_nerves=event.get('blocked_nerves',()))
    # Reject conflicting simultaneous values; time order cannot silently win.
    for i,a in enumerate(result):
        for b in result[i+1:]:
            if max(a['start_s'],b['start_s'])<min(a['end_s'],b['end_s']):
                for key in ['stimuli','motor_commands']:
                    if set(a.get(key,{}))&set(b.get(key,{})):
                        raise ValueError('Overlapping body intervention targets')
    return result


def inputs_at(events,time_s):
    out={'stimuli':{},'motor_commands':{},'blocked_nerves':set()}
    for event in events:
        if event['start_s']<=time_s<event['end_s']:
            out['stimuli'].update(event.get('stimuli',{}))
            out['motor_commands'].update(event.get('motor_commands',{}))
            out['blocked_nerves'].update(event.get('blocked_nerves',()))
    out['blocked_nerves']=sorted(out['blocked_nerves'])
    return out
