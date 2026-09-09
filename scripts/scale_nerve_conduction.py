"""A taller person has longer nerves, and therefore slower conduction.

    .venv/bin/python scripts/scale_nerve_conduction.py --stature-m 2.03
    .venv/bin/python scripts/scale_nerve_conduction.py --scale 1.0    # identity

This is the gate that says the parametrization reached the body rather than a
scale factor in a config file.

Conduction velocity is set by axon diameter and internodal myelin length.  Those
are cellular dimensions; they do not know how tall their owner is, and
`peripheral.json` declares all 17 fibre-class velocities as fixed constants with
sources.  Route length is anatomy and scales with stature.  So the delay --
`route / velocity` -- takes the length exponent exactly, and **every conduction
delay in the body lengthens in proportion to stature**.  At 2.03 m against the
source subject's 1.7973 m that is +12.95% on all 1,743 recorded delays plus every route the
fibre-class join derives one from.

Why it matters beyond being correct.  `docs/MILESTONES.md` records that lumping
every peripheral delay to step 0 costs as much as deleting an entire fibre group
from the interoceptive pathway.  Delays are load-bearing for the brain, so a body
of a different size genuinely has a different periphery, and the two facts
together are what makes stature a parameter of the *organism* rather than a
rendering scale.

**The independent arm.**  Multiplying 146 recorded route lengths by `s` proves
nothing about whether the body scaled -- it proves the multiplication ran.  So
every route length is instead RECOMPUTED, from geometry, by the method the route
itself declares:

    muscle bindings (249)  max(0.03, 1.15 * |anatomy centroid - relay|)
    receptor patches (16)  max(0.03, 1.15 * |patch position - relay|)
    polyline nerves (94)   sum of authored polyline segment lengths
    the other 52 nerves    upper median of their own bindings' lengths

and the endpoint of a muscle binding is an entity centroid read out of
`anatomy.json` -- a 17.8 MB artifact built by a different script, which neither
knows nor cares that anything is being scaled.  All 249 + 16 + 94 + 52 reproduce
their declared length exactly at the base, so the recomputation is entitled to be
believed at the scaled one.  Scale the routes without scaling the anatomy and
this arm fails on 249 of them while the delay arithmetic stays perfectly
self-consistent.

**The floor.**  The binding method is `max(0.03, ...)`, and a max is not
homogeneous: below 30 mm the recomputed length would stop scaling while the
recorded one carried on.  Whether that ever bites is measured, not hoped --
the shortest route is 58.1 mm, so the floor becomes active only below s = 0.516,
against 0.779 at the bottom of the declared stature range.
"""
from pathlib import Path
import argparse, hashlib, json, math, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm import body_scaling as bs                                  # noqa: E402
from ihm.body_parameters import MECHANICAL_STATURE_M, resolve       # noqa: E402

PERIPHERAL = 'data/derived/canonical/peripheral.json'
ANATOMY = 'data/derived/canonical/anatomy.json'
BINDING_FLOOR_M = 0.03
BINDING_FACTOR = 1.15

#: The delay fields, and the quantity each takes.  `central_delay_s` is here to
#: be excluded on purpose: it is synaptic and central transport, not a peripheral
#: route, and this scaler has no basis for moving it.  Leaving it out of the
#: table would look identical to forgetting it.
DELAY_FIELDS = {
    'motor_delay_s': 'conduction_delay',
    'afferent_delay_s': 'conduction_delay',
    'delays_s': 'conduction_delay',
    'central_delay_s': None,
}
LENGTH_FIELDS = ('path_length_m', 'rest_path_length_m', 'radius_m')

#: Parameters of the peripheral reduction that are lengths or velocities.
#: Velocities are invariant -- that is the whole argument -- and are listed so
#: the reader can see they were considered and held rather than missed.
PARAMETER_QUANTITY = {
    'tactile_velocity_m_s': None, 'warm_velocity_m_s': None,
    'cold_velocity_m_s': None, 'central_afferent_delay_s': None,
    'receptor_tau_s': None, 'activation_tau_s': None,
    'pressure_gain_hz_pa': None, 'stretch_gain_hz': None,
    'temperature_gain_hz_C': None, 'baseline_skin_temperature_C': None,
    'max_receptor_rate_hz': None,
}

#: Headline consequences, named in advance so the report cannot be read as
#: whatever came out.  Each is a route and a fibre class this body already
#: declares, and the expected scaled value is the base times s -- nothing else.
HEADLINE = (
    ('peripheral-nerve-left-vagus', 'c',
     'the vagal C fibre: visceral sensation, the slowest thing in the body'),
    ('peripheral-nerve-left-optic', 'retinal_magno', 'optic, magnocellular'),
    ('peripheral-nerve-left-optic', 'retinal_parvo', 'optic, parvocellular'),
    ('peripheral-nerve-left-optic', 'retinal_konio', 'optic, koniocellular'),
    ('peripheral-nerve-left-greater_splanchnic', 'c', 'greater splanchnic C fibre'),
    ('peripheral-nerve-left-median', 'abeta', 'median A-beta: touch'),
    ('peripheral-nerve-left-sciatic_tibial', 'ia', 'sciatic Ia: the stretch reflex'),
)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def velocities(peripheral):
    """One typical velocity per fibre class.  Invariant, and that is the point."""
    return {name: float(entry['typical'])
            for name, entry in peripheral['fibre_velocity_m_s'].items()}


def delays_for(length_m, classes, velocity):
    return {c: length_m / velocity[c] for c in classes if c in velocity}


# ---------------------------------------------------------------------------
# Recomputation from geometry -- the independent arm
# ---------------------------------------------------------------------------

def recompute_lengths(peripheral, centroids, scale, origin=(0.0, 0.0, 0.0),
                      anatomy_scale=None):
    """Every route length, rebuilt from scaled geometry by its declared method.

    Nothing here reads ``path_length_m``.  Positions come from the relay set,
    the receptor patches and -- for the 249 muscle bindings -- from anatomy
    entity centroids, and each is scaled before the distance is taken, so the
    length that comes out is a distance between scaled points rather than a
    scaled distance.
    """
    # `anatomy_scale` exists only so the sabotage arm can scale the nerves and
    # leave the body at its old size; in every real run it is `scale`.
    anatomy_scale = scale if anatomy_scale is None else anatomy_scale
    relay = {r['id']: bs.scale_point(r['position_m'], scale, origin)
             for r in peripheral['relays']}
    out = {'muscle_bindings': {}, 'receptor_patches': {}, 'nerves': {},
           'floor_active': 0}

    def endpoint_to_relay(point, relay_id):
        raw = BINDING_FACTOR * math.dist(point, relay[relay_id])
        if raw < BINDING_FLOOR_M:
            out['floor_active'] += 1
        return max(BINDING_FLOOR_M, raw)

    per_nerve = {}
    for binding in peripheral['muscle_bindings']:
        centroid = centroids.get(binding['endpoint_id'])
        if centroid is None:
            continue
        length = endpoint_to_relay(bs.scale_point(centroid, anatomy_scale, origin),
                                   binding['relay_id'])
        out['muscle_bindings'][binding['muscle_id']] = length
        per_nerve.setdefault(binding['nerve_id'], []).append(length)
    for patch in peripheral['receptor_patches']:
        length = endpoint_to_relay(bs.scale_point(patch['position_m'], scale, origin),
                                   patch['relay_id'])
        out['receptor_patches'][patch['id']] = length
        per_nerve.setdefault(patch['nerve_id'], []).append(length)
    for nerve in peripheral['nerves']:
        points = nerve.get('points_m')
        if points:
            out['nerves'][nerve['id']] = bs.polyline_length(
                [bs.scale_point(p, scale, origin) for p in points])
        else:
            seen = sorted(per_nerve.get(nerve['id'], []))
            if seen:
                out['nerves'][nerve['id']] = seen[len(seen) // 2]
    return out


# ---------------------------------------------------------------------------
# Scaling
# ---------------------------------------------------------------------------

def scale_record(record, scale):
    out = dict(record)
    for field in LENGTH_FIELDS:
        if out.get(field) is not None:
            out[field] = bs.scale('route_length', out[field], scale)
    for field, quantity in DELAY_FIELDS.items():
        if quantity is None or out.get(field) is None:
            continue
        value = out[field]
        out[field] = ({k: bs.scale(quantity, v, scale) for k, v in value.items()}
                      if isinstance(value, dict) else bs.scale(quantity, value, scale))
    if 'anchors' in out and out['anchors']:
        anchors = []
        for anchor in out['anchors']:
            anchor = dict(anchor)
            for key in ('point_m', 'source_point_ground_m'):
                if anchor.get(key):
                    anchor[key] = bs.scale_vector('length', anchor[key], scale)
            if anchor.get('projection_distance_m') is not None:
                anchor['projection_distance_m'] = bs.scale(
                    'length', anchor['projection_distance_m'], scale)
            anchors.append(anchor)
        out['anchors'] = anchors
    if out.get('points_m'):
        out['points_m'] = [bs.scale_point(p, scale) for p in out['points_m']]
    if out.get('position_m'):
        out['position_m'] = bs.scale_point(out['position_m'], scale)
    if out.get('max_isometric_force_n') is not None:
        out['max_isometric_force_n'] = bs.scale(
            'max_isometric_force', out['max_isometric_force_n'], scale)
    return out


def build(root, scale, anatomy_scale=None):
    root = Path(root)
    peripheral = json.loads((root / PERIPHERAL).read_text())
    anatomy = json.loads((root / ANATOMY).read_text())
    centroids = {e['id']: e['centroid_m'] for e in anatomy['entities']}
    velocity = velocities(peripheral)

    scaled = dict(peripheral)
    for key in ('nerves', 'relays', 'receptor_patches', 'muscle_bindings'):
        scaled[key] = [scale_record(r, scale) for r in peripheral[key]]
    scaled['stature_scale'] = scale
    scaled['scaling'] = dict(
        schema=bs.SCHEMA,
        route_length_exponent=bs.exponent('route_length'),
        conduction_delay_exponent=bs.exponent('conduction_delay'),
        conduction_velocity_exponent=bs.exponent('conduction_velocity'),
        velocity_basis=bs.PRIMITIVES['conduction_velocity']['basis'],
        held_invariant=sorted([k for k, v in DELAY_FIELDS.items() if v is None]
                              + [k for k, v in PARAMETER_QUANTITY.items() if v is None]),
        held_invariant_reason=(
            'central_delay_s and central_afferent_delay_s are synaptic and central '
            'transport, not peripheral route lengths, and nothing here fixes their '
            'scaling. The receptor and activation time constants are membrane '
            'properties. The three transduction velocities '
            '(tactile/warm/cold_velocity_m_s) are conduction velocities and are '
            'invariant for the same reason as every other one. The receptor gains '
            'are Hz per pascal and per degree, which are transduction properties '
            'of the ending and not of the body carrying it.'))

    gates, headline = run_gates(root, peripheral, scaled, centroids, velocity, scale,
                                anatomy_scale=anatomy_scale)
    report = dict(
        schema='ihm.nerve-conduction-scaling.v1', stature_scale=scale,
        routes=len(peripheral['nerves']),
        muscle_bindings=len(peripheral['muscle_bindings']),
        receptor_patches=len(peripheral['receptor_patches']),
        fibre_classes=len(velocity),
        delays_scaled=sum(len(r.get('delays_s') or {}) for r in scaled['muscle_bindings'])
                      + sum(1 for r in scaled['muscle_bindings'] for f in
                            ('motor_delay_s', 'afferent_delay_s') if r.get(f) is not None),
        conduction_delay=bs.explain('conduction_delay'),
        headline=headline, gates=gates,
        scaling=scaled['scaling'],
        limitations=[
            'Route lengths are authored schematic geometry, not dissected nerves. '
            'peripheral.json says so on every record and scaling does not improve '
            'it: a 12.95% longer schematic route is still a schematic route.',
            'Central and synaptic delays do not scale here. A whole reflex latency '
            'is therefore NOT 12.95% longer -- only its peripheral part is, and the '
            'report gives both so the difference is visible.',
            'One typical velocity per fibre class. The declared low/high bounds '
            'span a factor of 2-4 for several classes, which is larger than the '
            'entire stature effect; the stature effect is a shift of the whole '
            'distribution, not a claim to resolve within it.',
        ])
    return scaled, report


def run_gates(root, peripheral, scaled, centroids, velocity, scale,
              anatomy_scale=None):
    gates, headline = [], []

    def record(name, error, tolerance, note, **extra):
        gates.append(dict(gate=name, relative_error=error, tolerance=tolerance,
                          passed=(error <= tolerance), note=note, **extra))

    base_by_id = {n['id']: n for n in peripheral['nerves']}
    scaled_by_id = {n['id']: n for n in scaled['nerves']}

    # 1. Identity.
    if scale == 1.0:
        differences = sum(
            1 for key in ('nerves', 'relays', 'receptor_patches', 'muscle_bindings')
            for a, b in zip(peripheral[key], scaled[key]) if a != b)
        record('identity at s=1.0 reproduces the base exactly', float(differences), 0.0,
               'whole-record equality across all four record sets')

    # 2. Every delay takes exactly the length exponent.
    worst, n = 0.0, 0
    for a, b in zip(peripheral['muscle_bindings'], scaled['muscle_bindings']):
        for key, value in (a.get('delays_s') or {}).items():
            worst = max(worst, abs(b['delays_s'][key] / (value * scale) - 1.0)); n += 1
        for field in ('motor_delay_s', 'afferent_delay_s'):
            if a.get(field):
                worst = max(worst, abs(b[field] / (a[field] * scale) - 1.0)); n += 1
    record('every conduction delay scales by exactly s', worst, 1e-12,
           'delay = route / velocity, and velocity is invariant, so s**%g'
           % bs.exponent('conduction_delay'), delays=n)

    # 3. Delay consistency: the scaled delay must equal the scaled route over the
    #    UNCHANGED velocity. This is the arithmetic identity that would break if
    #    velocity had been scaled too -- the single most tempting mistake here,
    #    because scaling everything is easier than deciding what not to scale.
    worst, n = 0.0, 0
    for binding in scaled['muscle_bindings']:
        for key, value in (binding.get('delays_s') or {}).items():
            if key in velocity:
                worst = max(worst, abs(value * velocity[key]
                                       / binding['path_length_m'] - 1.0)); n += 1
    record('scaled delay = scaled route / UNSCALED velocity', worst, 1e-12,
           'if the velocities had been scaled along with everything else the '
           'delays would come out invariant and the whole result would vanish; '
           'this is the gate that says they did not', delays=n)

    # 4. THE INDEPENDENT ARM. Route lengths rebuilt from scaled geometry, with
    #    the muscle-binding endpoints read out of anatomy.json.
    base_recomputed = recompute_lengths(peripheral, centroids, 1.0)
    scaled_recomputed = recompute_lengths(peripheral, centroids, scale,
                                          anatomy_scale=anatomy_scale)
    for label, key, records in (
            ('muscle bindings', 'muscle_id', ('muscle_bindings', 'muscle_bindings')),
            ('receptor patches', 'id', ('receptor_patches', 'receptor_patches')),
            ('nerve routes', 'id', ('nerves', 'nerves'))):
        base_error, scaled_error, n = 0.0, 0.0, 0
        for original, moved in zip(peripheral[records[0]], scaled[records[1]]):
            name = original[key]
            got_base = base_recomputed[records[0]].get(name)
            got_scaled = scaled_recomputed[records[0]].get(name)
            if got_base is None or not original.get('path_length_m'):
                continue
            base_error = max(base_error, abs(got_base / original['path_length_m'] - 1.0))
            scaled_error = max(scaled_error, abs(got_scaled / moved['path_length_m'] - 1.0))
            n += 1
        record('base: %s reproduce their declared length from geometry' % label,
               base_error, 1e-9,
               'the recomputation earns the right to be believed at the scaled '
               'body by reproducing the base one exactly', routes=n)
        record('INDEPENDENT: %s rebuilt from SCALED geometry match' % label,
               scaled_error, 1e-9,
               'distances between scaled points, not scaled distances. For the '
               'muscle bindings the endpoint is an anatomy.json entity centroid, '
               'so this fails if the nerves were scaled and the anatomy was not.',
               routes=n)

    # 5. The floor. max(0.03, ...) is not homogeneous; measure whether it bites.
    shortest = min(r['path_length_m'] for r in peripheral['muscle_bindings'])
    record('the 30 mm route floor is inactive at this scale',
           float(scaled_recomputed['floor_active']), 0.0,
           'max(0.03, 1.15*d) stops being homogeneous below 30 mm. The shortest '
           'route is %.1f mm, so the floor activates only below s = %.3f, against '
           '%.3f at the bottom of the declared stature range.'
           % (1000 * shortest, BINDING_FLOOR_M / shortest, 1.40 / MECHANICAL_STATURE_M),
           shortest_route_mm=1000 * shortest,
           floor_active_below_scale=BINDING_FLOOR_M / shortest)

    # 6. The named headline consequences.
    composition = {}
    for binding in peripheral['muscle_bindings']:
        composition.setdefault(binding['nerve_id'], set()).update(
            (binding.get('delays_s') or {}).keys())
    worst = 0.0
    for nerve_id, fibre, description in HEADLINE:
        base = base_by_id.get(nerve_id)
        if base is None or fibre not in velocity:
            continue
        base_delay = base['path_length_m'] / velocity[fibre]
        got = scaled_by_id[nerve_id]['path_length_m'] / velocity[fibre]
        worst = max(worst, abs(got / (base_delay * scale) - 1.0))
        headline.append(dict(
            nerve=nerve_id, fibre_class=fibre, description=description,
            route_mm=1000 * base['path_length_m'],
            scaled_route_mm=1000 * scaled_by_id[nerve_id]['path_length_m'],
            velocity_m_s=velocity[fibre],
            delay_ms=1000 * base_delay, scaled_delay_ms=1000 * got,
            change_ms=1000 * (got - base_delay),
            change_percent=100 * (got / base_delay - 1.0)))
    record('the named headline delays scale by exactly s', worst, 1e-12,
           'named in advance in HEADLINE so the report cannot be read as whatever '
           'came out', consequences=len(headline))

    # 7. Separation. Three retinal populations on ONE optic nerve must separate
    #    in proportion: it is the same route, so their DIFFERENCES scale too.
    optic = [h for h in headline if h['nerve'] == 'peripheral-nerve-left-optic']
    if len(optic) >= 2:
        base_spread = max(h['delay_ms'] for h in optic) - min(h['delay_ms'] for h in optic)
        got_spread = (max(h['scaled_delay_ms'] for h in optic)
                      - min(h['scaled_delay_ms'] for h in optic))
        record('the optic nerve populations separate in proportion',
               abs(got_spread / (base_spread * scale) - 1.0), 1e-12,
               'one route, three velocities: the SPREAD between the fastest and '
               'slowest retinal population is a delay too and takes the same '
               'exponent. %.3f ms -> %.3f ms.' % (base_spread, got_spread),
               base_spread_ms=base_spread, scaled_spread_ms=got_spread)

    # 8. The invariant itself, asserted against the table rather than against
    #    this script's own output. Everything above rests on conduction velocity
    #    taking exponent 0; if it ever stops doing so, every delay silently
    #    becomes invariant and the whole result evaporates while every length
    #    still scales perfectly.
    record('conduction velocity carries exponent 0 in the scaling table',
           abs(bs.exponent('conduction_velocity')), 0.0,
           'the load-bearing invariant, checked at its source. '
           + bs.PRIMITIVES['conduction_velocity']['basis'],
           classes=len(velocity),
           velocities_unchanged=(peripheral['fibre_velocity_m_s']
                                 == scaled['fibre_velocity_m_s']))
    return gates, headline


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--stature-m', type=float)
    group.add_argument('--scale', type=float)
    parser.add_argument('--output', default=None)
    parser.add_argument('--sabotage', choices=('scale-the-velocities',
                                               'routes-without-anatomy'), default=None)
    args = parser.parse_args()

    stature = args.stature_m if args.stature_m is not None else args.scale * MECHANICAL_STATURE_M
    scale = resolve({'stature_m': stature})['derived']['stature_scale']
    out = ROOT / (args.output or ('data/derived/body-variants/stature_%s'
                                  % ('%.6f' % scale).replace('.', 'p')))
    out.mkdir(parents=True, exist_ok=True)

    if args.sabotage == 'scale-the-velocities':
        # The tempting mistake: scale everything, including the velocities, and
        # every delay comes out invariant while every length still scales.
        bs.PRIMITIVES['conduction_velocity']['exponent'] = 1
    # The other tempting mistake: scale the nerves and leave the body alone. The
    # delay arithmetic stays perfectly self-consistent and only the geometry arm
    # can see it, which is the whole reason that arm exists.
    scaled, report = build(
        ROOT, scale,
        anatomy_scale=1.0 if args.sabotage == 'routes-without-anatomy' else None)

    report['stature_m'] = stature
    report['sources'] = {p: _sha(ROOT / p) for p in (PERIPHERAL, ANATOMY)}
    report['script_sha256'] = _sha(Path(__file__))
    (out / 'peripheral_scaled.json').write_text(json.dumps(scaled, indent=1) + '\n')
    (out / 'conduction_report.json').write_text(json.dumps(report, indent=2) + '\n')

    for gate in report['gates']:
        print('%-5s %-64s %.3g (tol %.3g)'
              % ('PASS' if gate['passed'] else 'FAIL', gate['gate'],
                 gate['relative_error'], gate['tolerance']))
    print()
    print('%-34s %8s %8s %8s %8s' % ('', 'route mm', 'v m/s', 'delay ms', 'scaled ms'))
    for row in report['headline']:
        print('%-34s %8.1f %8.1f %8.2f %8.2f  %+.2f%%'
              % (row['description'][:34], row['route_mm'], row['velocity_m_s'],
                 row['delay_ms'], row['scaled_delay_ms'], row['change_percent']))
    print('\nwrote', out)
    if not all(g['passed'] for g in report['gates']):
        raise SystemExit('scaled conduction failed its own gates')


if __name__ == '__main__':
    main()
