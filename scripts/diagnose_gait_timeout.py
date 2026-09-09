#!/usr/bin/env python3
"""Why does one 10 ms native advance stop returning at t=1.83 s?

The search-J gait run dies on ``TimeoutError: Native mechanical stream response
timed out``, reproducibly, at the same instant.  That message is a *Python* read
deadline (``ihm/native/mechanical_stream.py``), not an engine diagnosis: it only
says the engine did not answer within 120 s.  The engine answers one ``advance``
by building a fresh ``Manager`` and running an error-controlled Simbody
integration over the interval, so its cost is not bounded by ``dt``.

This script times every advance, keeps the engine log, and raises the read
deadline so the engine is allowed to finish or to fail on its own terms.  It
dumps the plant state entering each slow advance so the cause is measured rather
than guessed.
"""
from __future__ import annotations

import argparse, json, math, os, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--params', default='data/derived/gait-search/searchJ_best_params.json')
    ap.add_argument('--horizon', type=float, default=2.2)
    ap.add_argument('--timeout', type=float, default=1800.0)
    ap.add_argument('--slow-s', type=float, default=2.0,
                    help='an advance slower than this is dumped in full')
    ap.add_argument('--out', default='data/derived/gait-timeout-diagnosis')
    a = ap.parse_args()
    os.environ['IHM_NATIVE_RESPONSE_TIMEOUT_S'] = str(a.timeout)

    import walk_gait as wg  # noqa: E402  (after the env var is set)

    out = ROOT / a.out
    out.mkdir(parents=True, exist_ok=True)
    trace_path = out / 'advance_trace.jsonl'
    trace = trace_path.open('w')

    base = wg.NativeMechanicalStream
    slow = a.slow_s

    class Timed(base):
        def advance(self, dt_s, forces=(), actuation=None):
            before = self.state
            t0 = time.time()
            try:
                result = base.advance(self, dt_s, forces=forces, actuation=actuation)
            except BaseException as exc:
                wall = time.time() - t0
                trace.write(json.dumps({
                    'time_s': before['time_s'], 'wall_s': wall, 'failed': True,
                    'exception': type(exc).__name__ + ': ' + str(exc)[:300],
                    'entering_state': summarize(before, actuation),
                }) + '\n')
                trace.flush()
                raise
            wall = time.time() - t0
            row = {'time_s': before['time_s'], 'wall_s': wall, 'failed': False}
            if wall >= slow:
                row['entering_state'] = summarize(before, actuation)
            trace.write(json.dumps(row) + '\n')
            trace.flush()
            return result

    wg.NativeMechanicalStream = Timed
    params = json.loads((ROOT / a.params).read_text())
    if 'best_report' in params:
        params = params['best_report']['params']
    full = dict(wg.SEED_PARAMS)
    full.update(params)
    policy = wg.Policy()
    work = ROOT / 'data/derived/gait-work/diagnose'
    report, frames = wg.rollout(full, policy, horizon_s=a.horizon, record=True,
                                work_dir=work, wall_budget_s=1e9)
    trace.close()
    report.pop('params', None)
    report.pop('phase_events', None)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    engine_log = work / 'engine.log'
    if engine_log.exists():
        (out / 'engine.log').write_text(engine_log.read_text()[-20000:])
    print(json.dumps({k: report[k] for k in
                      ('simulated_s', 'wall_s', 'stop_reason', 'failure',
                       'pelvis_forward_travel_m', 'max_clipped_muscles')}, indent=2))
    print('trace:', trace_path)


def summarize(state, actuation):
    coords = {k: {'value': v['value'], 'speed': v['speed']}
              for k, v in state['coordinates'].items()}
    speeds = sorted(((abs(v['speed']), k) for k, v in coords.items()), reverse=True)[:8]
    muscles = state['muscles']

    def worst(field, n=8):
        rows = [(abs(m.get(field, 0.0)), name) for name, m in muscles.items()
                if isinstance(m.get(field), (int, float))]
        rows.sort(reverse=True)
        return [{'muscle': name, field: value} for value, name in rows[:n]]

    fields = sorted({k for m in muscles.values() for k in m})
    contacts = [{'name': c['name'],
                 'force_n': math.sqrt(sum(x * x for x in c['force_n'])),
                 'center_m': c['center_m']}
                for c in state['contacts']]
    contacts.sort(key=lambda c: -c['force_n'])
    return {
        'coordinates': coords,
        'fastest_coordinates': [{'coordinate': k, 'abs_speed': v} for v, k in speeds],
        'muscle_fields': fields,
        'shortest_fibers': sorted(
            ({'muscle': n, 'fiber_length_m': m['fiber_length_m']}
             for n, m in muscles.items() if 'fiber_length_m' in m),
            key=lambda r: r['fiber_length_m'])[:8],
        'fastest_fibers': worst('fiber_velocity_m_s'),
        'highest_activation': sorted(
            ({'muscle': n, 'activation': m.get('activation')} for n, m in muscles.items()
             if isinstance(m.get('activation'), (int, float))),
            key=lambda r: -r['activation'])[:8],
        'contacts': contacts[:10],
        'foot_contact_force_n': dict(state.get('foot_contact_force_n', {})),
        'commanded_max_excitation': None if not actuation else max(actuation.values()),
        'commanded_at_ceiling': None if not actuation else
            sum(1 for v in actuation.values() if v >= 0.999),
    }


if __name__ == '__main__':
    main()
