"""Pure scheduler checks; no native process or physiological approximation."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ihm.assembly.intake_schedule import IntakeEvent, IntakeSchedule
from ihm.native.session import Meal


def rejects(call):
    try:
        call()
    except (ValueError, RuntimeError):
        return
    raise AssertionError('Expected rejection')


def main():
    food = IntakeEvent('food', .04, Meal(name='breakfast', carbohydrate_g=20, protein_g=5))
    water = IntakeEvent('water', .08, Meal(name='drink', water_ml=250))
    schedule = IntakeSchedule([water, food], horizon_s=1)
    copy = IntakeSchedule.from_checkpoint(json.loads(json.dumps(schedule.checkpoint())))
    assert copy.due(1) == ()
    assert copy.due(2) == (food,)
    rejects(copy.checkpoint)
    rejects(lambda: copy.due(2))
    rejects(lambda: copy.record_accepted('water', {}))
    rejects(lambda: copy.record_accepted('food', {'status': 'ok', 'sequence': True, 'elapsed_s': .04, 'pending_meal': True}))
    rejects(lambda: copy.record_accepted('food', {'status': 'ok', 'sequence': 3, 'elapsed_s': .06, 'pending_meal': True}))
    receipt = copy.record_accepted('food', {'status': 'ok', 'sequence': 3, 'elapsed_s': .04, 'pending_meal': True})
    assert receipt.outcome == 'accepted' and receipt.issued_tick == 2
    assert copy.due(2) == ()
    assert copy.due(4) == (water,)
    copy.record_uncertain('water', 'native acknowledgment timed out')
    rejects(lambda: copy.due(5))
    rejects(copy.checkpoint)
    rejects(lambda: copy.record_accepted('water', {}))
    assert [r.outcome for r in copy.receipts] == ['accepted', 'uncertain']
    rejects(lambda: IntakeSchedule([food, food]))
    for seconds in (.021, float('nan'), float('inf'), -1, True, 86400.02):
        rejects(lambda seconds=seconds: IntakeEvent('bad', seconds, Meal(water_ml=1)))
    for amount in (float('nan'), float('inf'), -1, 10001, True):
        rejects(lambda amount=amount: IntakeEvent('bad', 0, Meal(water_ml=amount)))
    for tick in (.5, True, -1, 51):
        rejects(lambda tick=tick: schedule.due(tick))
    rejects(lambda: IntakeSchedule([water], horizon_s=.04))
    schedule.due(1)
    rejects(lambda: schedule.due(0))
    snapshot = schedule.checkpoint()
    snapshot['issued'] = ['food']
    rejects(lambda: IntakeSchedule.from_checkpoint(snapshot))
    simultaneous = IntakeSchedule([IntakeEvent('z', 0, water.meal), IntakeEvent('a', 0, food.meal)])
    assert simultaneous.due(0)[0].event_id == 'a'
    rejects(lambda: simultaneous.record_uncertain('a', ''))
    simultaneous.record_uncertain('a', 'caller failed after reservation')
    accepted = IntakeSchedule([food, water], horizon_s=1)
    assert accepted.due(3) == (food,)  # Late delivery records actual issuing tick.
    accepted.record_accepted('food', {'status': 'ok', 'sequence': 3, 'elapsed_s': .06, 'pending_meal': True})
    assert accepted.due(4) == (water,)
    rejects(lambda: accepted.record_accepted('water', {'status': 'ok', 'sequence': 3, 'elapsed_s': .08, 'pending_meal': True}))
    accepted.record_accepted('water', {'status': 'ok', 'sequence': 5, 'elapsed_s': .08, 'pending_meal': True})
    assert accepted.due(50) == ()
    assert accepted.receipts[0].scheduled_tick == 2 and accepted.receipts[0].issued_tick == 3
    live = IntakeSchedule([], horizon_s=1)
    live.add_events([food], 2)
    assert IntakeSchedule.from_checkpoint(live.checkpoint()).snapshot() == live.snapshot()
    assert live.snapshot()['events'][0]['state'] == 'queued'
    live.due(2)
    assert live.snapshot()['events'][0]['state'] == 'issued'
    rejects(lambda: live.add_events([water], 2))
    live.record_accepted('food', {'status': 'ok', 'sequence': 1, 'elapsed_s': .04, 'pending_meal': True})
    before = live.snapshot()
    rejects(lambda: live.add_events([water, food], 2))
    assert live.snapshot() == before
    live.add_events([water], 3)
    assert [e['state'] for e in live.snapshot()['events']] == ['accepted', 'queued']
    rejects(lambda: live.add_events([IntakeEvent('old', 0, water.meal)], 3))
    live.due(4)
    live.record_uncertain('water', 'unknown delivery')
    assert live.snapshot()['events'][1]['state'] == 'uncertain'
    rejects(lambda: live.add_events([], 4))
    print('PASS: timed food/water, identity, validation, monotonic ticks, explicit receipts, uncertainty, checkpoint boundary')


if __name__ == '__main__':
    main()
