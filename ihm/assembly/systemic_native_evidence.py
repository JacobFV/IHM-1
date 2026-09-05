"""Bind systemic observations and actions to retained native execution records."""
import json
import math
from pathlib import Path
from ihm.native import _sha

def resolve_systemic_execution(root,directory,data):
    root=Path(root).resolve();directory=Path(directory).resolve();native=directory/'native'
    manifest=native/'manifest.json';frames_path=directory/'frames.jsonl';journal=native/'receipts.jsonl'
    if json.loads(manifest.read_text())!=data['native_manifest']:
        raise ValueError('Embedded native manifest differs from retained execution')
    frames=[json.loads(line) for line in frames_path.read_text().splitlines()]
    if frames!=data['frames']:raise ValueError('Systemic observations differ from retained sampled frames')
    samples={round(frame['time_s']*50):frame['values'] for frame in frames}
    sources={str(p.relative_to(root)):_sha(p) for p in (manifest,frames_path,journal)}
    seen=set();actions=[];saved=[];pending=None;sequence=0;ticks=0;closed=False;last_op=None
    for line in journal.open():
        record=json.loads(line)
        if 'command' in record:
            if pending is not None or closed:raise ValueError('Incomplete native command sequence')
            words=record['command'].split();sequence+=1
            if len(words)<2 or int(words[0])!=sequence or record['before_ticks']!=ticks:
                raise ValueError('Invalid native command identity or clock')
            pending=words
        elif 'acknowledgment' in record:
            ack=record['acknowledgment']
            if pending is None or ack['sequence']!=sequence:raise ValueError('Unpaired native acknowledgment')
            if any(isinstance(ack.get(k),bool) or not isinstance(ack.get(k),(float,int)) or not math.isfinite(ack[k]) for k in ('elapsed_s','time_s','origin_s')):
                raise ValueError('Invalid native acknowledgment clock')
            if abs(ack['time_s']-ack['origin_s']-ack['elapsed_s'])>1e-4:
                raise ValueError('Native absolute/elapsed clocks differ')
            op=pending[1];arguments=pending[2:]
            if op=='step':
                if len(arguments)!=1 or int(arguments[0])<=0:raise ValueError('Invalid native step')
                ticks+=int(arguments[0])
            elif op in ('meal','apnea','exercise'):
                if op=='meal':
                    if len(arguments)!=7:raise ValueError('Invalid native meal receipt')
                    value=dict(zip(('name','carbohydrate_g','protein_g','fat_g','sodium_g','calcium_mg','water_ml'),
                                   [arguments[0],*[float(x) for x in arguments[1:]]]))
                else:
                    if len(arguments)!=1:raise ValueError('Invalid native action receipt')
                    value=float(arguments[0])
                actions.append(dict(time_s=ticks*.02,kind=op,value=value))
            elif op not in ('save','snapshot','quit'):raise ValueError('Unknown native command in systemic execution')
            if abs(ack['elapsed_s']-ticks*.02)>1e-8 or ack['status']!=('closed' if op=='quit' else 'ok'):
                raise ValueError('Native acknowledgment did not complete the expected command')
            if ticks in samples and ticks not in seen:
                if ack.get('values')!=samples[ticks]:raise ValueError('Sample differs from first native observation at its time')
                seen.add(ticks)
            closed=op=='quit';pending=None;last_op=op
        elif 'saved_state_sha256' in record:
            if pending is not None or last_op!='save' or record['sequence']!=sequence:raise ValueError('Unpaired native saved state')
            # Recorded absolute execution roots may have been relocated. Only
            # the prescribed unique state filename within this run is selected.
            expected=f'native_session_{sequence:08d}.xml'
            if Path(record['path']).name!=expected:raise ValueError('Unexpected saved state identity')
            path=native/'states'/expected
            if _sha(path)!=record['saved_state_sha256']:raise ValueError('Native saved state changed')
            sources[str(path.relative_to(root))]=record['saved_state_sha256'];saved.append(record['saved_state_sha256'])
        else:raise ValueError('Unknown native execution record')
    if pending is not None or not closed or ticks!=round(data['configuration']['seconds']*50):
        raise ValueError('Incomplete native systemic execution')
    if seen!=samples.keys() or actions!=data['actions']:raise ValueError('Systemic samples/actions differ from native execution')
    if len(saved)!=2 or saved[0]!=data['initial_state_sha256'] or saved[-1]!=data['final_state_sha256']:
        raise ValueError('Systemic checkpoint identity mismatch')
    return sources
