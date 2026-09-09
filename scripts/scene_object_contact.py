#!/usr/bin/env python3
"""Can the body feel the ball that is already in its scene?

The environment catalogue's own `not_selectable` list answers no, and says why:

    body-object contact -- "No engine solves it.  The interactive scene declares
    body_object_contact=false and objects contact only the environment plane."

That is true of `ihm/assembly/interactive_scene.py`, which integrates its
0.065 m / 0.4 kg ball in Python with `sphere_plane_step` against a half space and
never against the body.  It is NOT true of the machinery available: the installed
Simbody carries `CollisionDetectionAlgorithm::SphereSphere`, `OpenSim::HuntCrossleyForce`
exposes it over named `ContactGeometry`, and the native plant already carries 28
`ContactSphere` elements over 20 bodies in the `upright` environment.  So the gap
is a connection that was never made, not a solver that does not exist.

This script makes it and measures it.  Two experiments:

``drop``   the catalogue's own ball released above a prone body, which is
           unambiguous: a named body element either reports a contact force or it
           does not.
``push``   the ball resting on the floor in the path of a crawling body, which is
           the one that matters: contact that changes where the object ends up.

Everything reported is a measured force between two named contact elements.  The
body elements are the engine's inertia-inscribed proxies, not an anatomical skin
surface, and the report says so.
"""
from __future__ import annotations

import argparse, json, math, shutil, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402
import crawl  # noqa: E402

CATALOGUE = ROOT / 'data/derived/environment-catalogue-v1/catalogue.json'
WORK = ROOT / 'data/derived/crawl-work'
OUT = ROOT / 'data/derived/scene-object-contact'


def catalogue_ball(ident='ball-small'):
    """Radius and mass as the environment catalogue declares them."""
    entry = next(o for o in json.loads(CATALOGUE.read_text())['objects'] if o['id'] == ident)
    bounds = entry['bounds_m']
    radius = (bounds['max'][0] - bounds['min'][0]) / 2
    return {'catalogue_id': ident, 'label': entry['label'], 'radius_m': radius,
            'mass_kg': entry['mass_kg'],
            'scene_position_m': [(a + b) / 2 for a, b in zip(bounds['min'], bounds['max'])],
            'catalogue_geometry_sha256': entry['geometry_sha256']}


def run(mode, seconds, ball, position, params, trace_every=1):
    out = WORK / ('scene-' + mode)
    shutil.rmtree(out, ignore_errors=True)
    native = NativeMechanicalStream(
        ROOT, out, environment='upright', target_mass_kg=crawl.TARGET_MASS_KG,
        initial_pose=crawl.PRONE_POSE,
        augmented_registration='data/models/engineering_stance_v1/registration.json',
        coordinate_limits=crawl.joint_stops(),
        scene_objects=[{'id': 'ball_small', 'radius_m': ball['radius_m'],
                        'mass_kg': ball['mass_kg'], 'position_m': list(position)}])
    frames, touches = [], []
    started = time.time()
    failed = False
    try:
        state = native.snapshot()
        muscles = sorted(state['muscles'])
        pattern = crawl.CrawlPattern(params, muscles) if mode == 'push' else None
        rest = {k: (0.3 if k.startswith('pro_sup') else 0.0)
                for k in [a + '_' + s for a in crawl.ARM_PORT for s in 'rl'] + list(crawl.LUMBAR_PORT)}
        start = list(state['scene_objects'][0]['position_source_m'])
        for i in range(int(round(seconds / crawl.DT))):
            t_s = state['time_s']
            if pattern is None:
                excitation = {m: 0.02 for m in muscles}
                targets = rest
            else:
                excitation = dict(zip(muscles, map(float, pattern.excitation(t_s))))
                targets = pattern.limb_targets(t_s)
            obj = state['scene_objects'][0]
            if i % trace_every == 0:
                frames.append({'time_s': t_s,
                               'ball_position_m': obj['position_source_m'],
                               'ball_velocity_m_s': obj['velocity_source_m_s'],
                               'body_contact_elements_touching': obj['body_contact_elements_touching'],
                               'body_contact_magnitude_n': obj['body_contact_magnitude_n'],
                               'pelvis_tx': state['coordinates']['pelvis_tx']['value']})
            for c in obj['body_contacts']:
                touches.append({'time_s': t_s, 'body_element': c['body_element'],
                                'body_frame': c['body_frame'], 'magnitude_n': c['magnitude_n'],
                                'force_on_object_n': c['force_on_object_n']})
            state = native.advance(crawl.DT, actuation=excitation,
                                   coordinate_actuation=crawl.limb_pd(state, targets))
        final = list(state['scene_objects'][0]['position_source_m'])
        execution = json.loads((out / 'execution.json').read_text())
    except BaseException:
        failed = True
        raise
    finally:
        native.close()
        if failed:
            # Keep the engine log: a native termination is diagnosable only from it.
            kept = OUT / ('failed-' + mode)
            OUT.mkdir(parents=True, exist_ok=True)
            shutil.rmtree(kept, ignore_errors=True)
            shutil.copytree(out, kept, ignore=shutil.ignore_patterns('inputs'))
            print('kept engine output at', kept, file=sys.stderr)
        shutil.rmtree(out, ignore_errors=True)
    return {'frames': frames, 'touches': touches, 'start_m': start, 'final_m': final,
            'wall_s': time.time() - started,
            'scene_contact_material': execution['scene_contact_material'],
            'scene_object_contact_basis': execution['scene_object_contact_basis']}


def summarize(result, mass_kg):
    """Per-element contact, and the object's momentum change it has to explain.

    The impulse is the VECTOR integral of the contact force, not a running sum of
    its magnitude.  A sum of magnitudes over a rolling contact is 36x the actual
    momentum change here, because most of the force is normal to the motion and
    is carried straight into the floor; reporting it as an impulse would be the
    same shape of error as every entry in the ledger.
    """
    touches = result['touches']
    by_element = {}
    total_impulse = [0.0, 0.0, 0.0]
    for t in touches:
        row = by_element.setdefault(t['body_element'], {
            'body_element': t['body_element'], 'body_frame': t['body_frame'],
            'contact_frames': 0, 'peak_n': 0.0, 'impulse_on_object_n_s': [0.0, 0.0, 0.0],
            'first_contact_s': t['time_s'], 'last_contact_s': t['time_s']})
        row['contact_frames'] += 1
        row['peak_n'] = max(row['peak_n'], t['magnitude_n'])
        for k in range(3):
            row['impulse_on_object_n_s'][k] += t['force_on_object_n'][k] * crawl.DT
            total_impulse[k] += t['force_on_object_n'][k] * crawl.DT
        row['last_contact_s'] = t['time_s']
    for row in by_element.values():
        row['impulse_magnitude_n_s'] = math.sqrt(sum(v * v for v in row['impulse_on_object_n_s']))
    rows = sorted(by_element.values(), key=lambda r: -r['peak_n'])
    moved = [b - a for a, b in zip(result['start_m'], result['final_m'])]
    speeds = [math.sqrt(sum(v * v for v in f['ball_velocity_m_s'])) for f in result['frames']]
    return {
        'contacted_body_elements': rows,
        'total_contact_frames': len(touches),
        'peak_contact_force_n': max((t['magnitude_n'] for t in touches), default=0.0),
        'contact_impulse_on_object_n_s': total_impulse,
        'contact_impulse_magnitude_n_s': math.sqrt(sum(v * v for v in total_impulse)),
        'object_peak_speed_m_s': max(speeds, default=0.0),
        'object_peak_momentum_kg_m_s': mass_kg * max(speeds, default=0.0),
        'ball_start_m': result['start_m'], 'ball_final_m': result['final_m'],
        'ball_displacement_m': moved,
        'ball_horizontal_displacement_m': math.hypot(moved[0], moved[2]),
        'wall_s': result['wall_s'],
        'scene_contact_material': result['scene_contact_material'],
        'basis': result['scene_object_contact_basis'],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', choices=('drop', 'push', 'both'), default='both')
    ap.add_argument('--seconds', type=float, default=1.6)
    ap.add_argument('--push-seconds', type=float, default=8.0)
    ap.add_argument('--push-x', type=float, default=0.75)
    ap.add_argument('--control-x', type=float, default=3.0,
                    help='out-of-reach arm: the same run with the ball beyond the body')
    ap.add_argument('--params', default='data/derived/crawl-search/searchB.json')
    ap.add_argument('--out', default=str(OUT))
    a = ap.parse_args()
    ball = catalogue_ball()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    report = {'schema': 'ihm.scene-object-contact.v1', 'ball': ball,
              'body_contact_elements': 'the engine\'s 28 upright contact elements over 20 bodies; '
                                       'the non-foot ones are inertia-inscribed proxies at segment '
                                       'centres of mass, not an anatomical skin surface'}
    params = dict(crawl.SEED)
    loaded = json.loads((ROOT / a.params).read_text()) if Path(ROOT / a.params).exists() else {}
    params.update(loaded.get('best_report', {}).get('params', loaded.get('params', {})))
    def write():
        (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    if a.mode in ('drop', 'both'):
        # released above the prone torso proxy, whose centre settles near x = 0.367
        result = run('drop', a.seconds, ball, [0.367, 0.90, 0.0], params)
        report['drop'] = summarize(result, ball['mass_kg'])
        (out / 'drop_frames.json').write_text(json.dumps(result['frames']) + '\n')
        write()
    if a.mode in ('push', 'both'):
        result = run('push', a.push_seconds, ball, [a.push_x, ball['radius_m'], 0.0], params)
        report['push'] = summarize(result, ball['mass_kg'])
        report['push']['crawl_parameters'] = params
        (out / 'push_frames.json').write_text(json.dumps(result['frames']) + '\n')
        write()
        # The control the result needs: the identical run with the ball placed
        # beyond the body's reach.  If it moves here too, the push means nothing.
        control = run('control', a.push_seconds, ball,
                      [a.control_x, ball['radius_m'], 0.0], params)
        report['control_out_of_reach'] = summarize(control, ball['mass_kg'])
        report['control_out_of_reach']['placed_at_x_m'] = a.control_x
    write()
    printable = {k: v for k, v in report.items() if k not in ('drop', 'push')}
    for key in ('drop', 'push', 'control_out_of_reach'):
        if key in report:
            printable[key] = {k: v for k, v in report[key].items()
                              if k not in ('scene_contact_material', 'basis', 'crawl_parameters')}
    print(json.dumps(printable, indent=2))
    print('wrote', (out / 'report.json').relative_to(ROOT))


if __name__ == '__main__':
    main()
