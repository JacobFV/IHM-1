"""The brain/body signal graph the app's ring draws, and the magnitudes on it.

Two endpoints are served from here.

``/api/brain/graph`` is the static structure: the systems that can be
highlighted, and every route or materialization that touches one. It is derived
from declarations, never authored here:

- ``data/derived/canonical/peripheral.json`` — this repository's own nerve
  routes, each with a measured ``path_length_m`` over the body mesh, a relay,
  and its evidence flags. Fibre-class *composition* is explicitly not owned by
  this file (``fibre_velocity_scope``: "Composition is owned by IBM"), so the
  classes come from IBM-1's ``ibm/topologies/nerve.py`` and the delay is
  ``measured length / typical velocity of that class`` — one route, one delay
  per class, which is the whole point of the join.
- IBM-1's ``ibm/embodiment.py`` and ``ibm/anatomy/muscles.py`` for the port
  count and the cord's muscle/segment mapping. Those are read with ``ast``
  rather than imported, the way ``ihm/assembly/ibm_controller.py`` already does,
  so nothing here depends on IBM-1 being installed. If IBM-1 is absent the
  affected counts come back ``null`` with a reason: an unavailable number is
  reported as unavailable, never substituted.

``/api/brain/state`` is the per-edge magnitude. **The honest part.** Most of
these paths are near zero and a few are not, and the interface sizes its arrows
by this number, so what it means has to be exact. Every edge carries a
``magnitude_kind``:

``measured``
    A number this project actually measured, with the receipt named in
    ``source``. ``0.0005`` on the cortical sheet's internal sensory-to-motor hop
    is a measurement, not an absence of one: 0.03-0.07% of the driven signal
    arrives at the disjoint readout (IBM-1 ``docs/LOG.md``, 2026-09-09). A dead
    path drawn from a measurement is a result; the arrow should look dead.
``live``
    Read out of the current embodied frame this tick.
``declared``
    A gain or latency the model declares — a cord arc's +0.40 at 30 ms — which
    is a parameter, not a flow. It is not a magnitude and is not drawn as one.
``unmeasured``
    No measurement exists for this edge. Distinct from a measured zero, and the
    interface must not let the two look alike.
"""
from __future__ import annotations

import ast
import json
from pathlib import Path

IBM_ROOT = Path.home() / 'Documents/IBM-1'

#: Systems a reader can highlight. A system is a thing signal is exchanged
#: *with*; the ring around it is everything currently touching it.
SYSTEMS = (
    {'id': 'cortical_sheet', 'label': 'Cortical sheet',
     'kind': 'materialization',
     'note': '1,024 materialized sites of the IBM E/I equations, integrating at 1 ms '
             'and exchanging commands with native mechanics at 10 ms.'},
    {'id': 'precentral', 'label': 'Precentral (motor)', 'kind': 'region',
     'note': 'The motor readout population. Disjoint from the sensory sites, '
             'which is the measurement below.'},
    {'id': 'postcentral', 'label': 'Postcentral (somatosensory)', 'kind': 'region',
     'note': 'Where somatic afference lands. Every muscle binding in this body '
             'targets it.'},
    {'id': 'cord', 'label': 'Spinal cord', 'kind': 'process',
     'note': 'Per-segment alpha/gamma pools closing three declared reflex arcs '
             'between descending drive and muscle.'},
    {'id': 'upper_limb', 'label': 'Upper limbs', 'kind': 'body',
     'note': 'Muscle bindings and cutaneous receptor patches of both arms.'},
    {'id': 'lower_limb', 'label': 'Lower limbs', 'kind': 'body',
     'note': 'Muscle bindings and cutaneous receptor patches of both legs.'},
    {'id': 'viscera', 'label': 'Viscera', 'kind': 'body',
     'note': 'Visceral routes — vagus, splanchnics, pelvic. Topology only: '
             'there is no autonomic controller.'},
    {'id': 'special_sense', 'label': 'Special sense', 'kind': 'body',
     'note': 'Optic, cochlear, vestibular and olfactory routes. Declared, with '
             'measured lengths; no transduction runs on them here.'},
)
SYSTEM_IDS = {s['id'] for s in SYSTEMS}

#: Region ids in the peripheral declarations, mapped to the systems above.
_REGION_SYSTEM = {'upper_limb': 'upper_limb', 'lower_limb': 'lower_limb'}
_BRAIN_SYSTEM = {'postcentral': 'postcentral', 'precentral': 'precentral'}


def _literal(tree, name):
    """A module-level literal, read without importing the module."""
    for node in tree.body:
        target = getattr(node, 'target', None)
        targets = [target] if target is not None else getattr(node, 'targets', [])
        for t in targets:
            if getattr(t, 'id', None) == name:
                return node.value
    raise KeyError(name)


def _resolve(node, env):
    """literal_eval, but a bare Name resolves against already-read literals.

    ``TRUNK_COMPOSITION`` is written as ``{"median": MIXED, ...}``; refusing to
    resolve MIXED would mean retyping the composition table here, which is
    exactly the drift this module exists to avoid.
    """
    if isinstance(node, ast.Name):
        return env[node.id]
    if isinstance(node, ast.Dict):
        return {_resolve(k, env): _resolve(v, env) for k, v in zip(node.keys, node.values)}
    if isinstance(node, (ast.Tuple, ast.List)):
        out = []
        for element in node.elts:
            if isinstance(element, ast.Starred):
                out.extend(_resolve(element.value, env))
            else:
                out.append(_resolve(element, env))
        return tuple(out) if isinstance(node, ast.Tuple) else out
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        return _resolve(node.left, env) + _resolve(node.right, env)
    return ast.literal_eval(node)


def read_ibm(ibm_root=None):
    """IBM-1's declarations: fibre composition, port count, cord mapping.

    Returns ``{'available': False, 'reason': ...}`` rather than substituting
    anything when IBM-1 is not on this machine.
    """
    root = Path(ibm_root or IBM_ROOT)
    nerve = root / 'ibm/topologies/nerve.py'
    muscles = root / 'ibm/anatomy/muscles.py'
    embodiment = root / 'ibm/embodiment.py'
    missing = [str(p) for p in (nerve, muscles, embodiment) if not p.is_file()]
    if missing:
        return {'available': False,
                'reason': 'IBM-1 declarations not readable here: ' + ', '.join(missing)}
    try:
        tree = ast.parse(nerve.read_bytes())
        env = {}
        for name in ('SENSORY', 'CUTANEOUS', 'MUSCULAR', 'MIXED', 'AUTONOMIC'):
            env[name] = _resolve(_literal(tree, name), env)
        composition = _resolve(_literal(tree, 'TRUNK_COMPOSITION'), env)
        trunk_mm = _resolve(_literal(tree, 'TRUNK_LENGTH_MM'), env)
        innervation = ast.literal_eval(_literal(ast.parse(muscles.read_bytes()), 'INNERVATION'))
        sensor_ports = _resolve(_literal(ast.parse(embodiment.read_bytes()), 'SENSOR_PORTS'), {})
    except (SyntaxError, ValueError, KeyError, TypeError) as error:
        return {'available': False, 'reason': f'IBM-1 declarations unreadable: {error}'}
    return {'available': True, 'root': str(root), 'composition': composition,
            'trunk_length_mm': trunk_mm, 'innervation': innervation,
            'sensor_ports': sensor_ports}


#: cord.py's own segment list, C1-Co1.
LEVELS = ('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8',
          't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 't10', 't11', 't12',
          'l1', 'l2', 'l3', 'l4', 'l5', 's1', 's2', 's3', 's4', 's5', 'co1')

#: The four arcs cord.py closes, read from its own table. Gains and latencies,
#: which are declared parameters and not flows: they are never drawn as width.
ARCS = (
    ('stretch', 'Monosynaptic stretch · Ia → alpha, same segment', 0.40, 0.030),
    ('reciprocal', 'Reciprocal inhibition · Ia → antagonist alpha', -0.25, 0.032),
    ('autogenic', 'Autogenic inhibition · Ib → alpha', -0.15, 0.034),
    ('renshaw', 'Renshaw recurrent · alpha → alpha', -0.20, 0.004),
)


def _bare(nerve_id):
    parts = nerve_id.split('-')
    if len(parts) >= 4 and parts[0] == 'peripheral' and parts[1] == 'nerve':
        return '_'.join(parts[3:]) or parts[-1], parts[2]
    return nerve_id, 'unknown'


def _system_for_route(record, region_hint):
    kind = record.get('route_kind')
    if kind == 'special_sense':
        return 'special_sense'
    if kind == 'visceral':
        return 'viscera'
    return _REGION_SYSTEM.get(region_hint, 'cord')


def _mean(points):
    points = [p for p in points if isinstance(p, (list, tuple)) and len(p) == 3]
    if not points:
        return None
    return [sum(p[i] for p in points) / len(points) for i in range(3)]


def _anchors(root, peripheral):
    """Where each system sits in the canonical body frame, for the leader lines.

    Read from the declarations that already carry positions — the registered
    cortical nodes in ``brain.json``, the peripheral relays, and the receptor
    patches — all in ``bodyparts3d-display-m``, the frame the scene draws in.
    These are drawing anchors: a leader has to land somewhere, and landing it on
    the declared centroid is better than landing it on a guess.
    """
    relays = peripheral.get('relays', [])
    patches = peripheral.get('receptor_patches', [])
    nerves = peripheral.get('nerves', [])
    out = {}
    try:
        brain = json.loads((root / 'data/derived/canonical/brain.json').read_bytes())
        nodes = brain.get('nodes', [])
    except (OSError, ValueError):
        nodes = []
    label = lambda node: ((node.get('ibm_partition') or {}).get('label'))
    cortical = [n['position_m'] for n in nodes
                if n.get('kind') == 'cortical_population' and n.get('position_m')]
    out['cortical_sheet'] = _mean(cortical)
    for region in ('precentral', 'postcentral'):
        out[region] = _mean([n['position_m'] for n in nodes
                             if label(n) == region and n.get('position_m')]) or out['cortical_sheet']
    out['cord'] = _mean([r['position_m'] for r in relays
                         if r.get('kind') == 'spinal_segment_group' and r.get('position_m')])
    # A bilateral system averaged across both sides anchors on the midline,
    # which is on the trunk and not on either limb. One side is the honest
    # place to land the line, and the label already says "limbs", plural.
    for region in ('upper_limb', 'lower_limb'):
        out[region] = _mean([p['position_m'] for p in patches
                             if p.get('sensorimotor_region') == region
                             and p.get('side') == 'left' and p.get('position_m')])
    endpoint = lambda n: (n.get('points_m') or [None])[0]
    out['viscera'] = _mean([endpoint(n) for n in nerves if n.get('route_kind') == 'visceral'])
    out['special_sense'] = _mean([endpoint(n) for n in nerves
                                  if n.get('route_kind') == 'special_sense'])
    return out


def brain_graph(root, ibm_root=None):
    """The static structure: systems, contributors, and what each one carries."""
    root = Path(root)
    peripheral = json.loads((root / 'data/derived/canonical/peripheral.json').read_bytes())
    ibm = read_ibm(ibm_root)
    velocity = peripheral.get('fibre_velocity_m_s') or {}
    composition = ibm.get('composition') or {}
    contract = peripheral.get('route_contract') or {}
    aliases = contract.get('aliases') or {}

    # Which body region each nerve serves, taken from the bindings that name it
    # rather than from the nerve's name.
    region = {}
    for record in peripheral.get('muscle_bindings', []) + peripheral.get('receptor_patches', []):
        hint = record.get('sensorimotor_region')
        if hint:
            region.setdefault(record.get('nerve_id'), hint)
    binding_count, patch_count = {}, {}
    for record in peripheral.get('muscle_bindings', []):
        binding_count[record.get('nerve_id')] = binding_count.get(record.get('nerve_id'), 0) + 1
    for record in peripheral.get('receptor_patches', []):
        patch_count[record.get('nerve_id')] = patch_count.get(record.get('nerve_id'), 0) + 1

    relay_position = {r['id']: r.get('position_m') for r in peripheral.get('relays', [])}
    edges, joined, unjoined = [], 0, []
    for record in peripheral.get('nerves', []):
        name, side = _bare(record['id'])
        # The alias table is the declaration's own (`route_contract.aliases`),
        # and a name it lists as unmatched stays unmatched: an unjoined route is
        # reported, never guessed into the nearest trunk.
        classes = composition.get(aliases.get(name, name))
        length = record.get('path_length_m')
        if not isinstance(length, (int, float)) or length <= 0:
            continue
        fibres = []
        for cls in classes or ():
            typical = (velocity.get(cls) or {}).get('typical')
            if typical:
                fibres.append({'fibre_class': cls, 'velocity_m_s': typical,
                               'delay_s': length / typical})
        if fibres:
            joined += 1
        else:
            unjoined.append(record['id'])
        system = _system_for_route(record, region.get(record['id']))
        edges.append({
            'id': record['id'],
            'label': record.get('name', record['id']),
            'kind': 'route',
            'route_kind': record.get('route_kind') or 'somatic',
            'side': side,
            'system': system,
            # A somatic route reaches postcentral through its relay; that is what
            # the bindings themselves declare as brain_target_id.
            'target': 'postcentral' if system in ('upper_limb', 'lower_limb', 'cord') else 'cortical_sheet',
            'relay_id': record.get('relay_id'),
            'anchor_m': (record.get('points_m') or [None])[0] or relay_position.get(record.get('relay_id')),
            'path_length_m': length,
            'path_length_scope': record.get('path_length_scope'),
            'length_method': record.get('length_method'),
            'fibres': sorted(fibres, key=lambda f: f['delay_s']),
            'fibre_source': 'IBM-1 ibm/topologies/nerve.py TRUNK_COMPOSITION'
                            if fibres else None,
            'muscle_ports': binding_count.get(record['id'], 0),
            'receptor_ports': patch_count.get(record['id'], 0),
            'downstream_brain_target_ids': record.get('downstream_brain_target_ids') or [],
            'measured_axon_geometry': bool(record.get('measured_axon_geometry')),
            'evidence_kind': record.get('evidence_kind'),
            'geometry_kind': record.get('geometry_kind'),
            'runtime_support': record.get('runtime_support'),
            'limitations': record.get('limitations') or [],
            # A declared route with no measurement of what it is carrying. This
            # is NOT a measured zero, and the interface must not draw it as one.
            'magnitude': None,
            'magnitude_kind': 'unmeasured',
            'measurement': 'No per-route magnitude is measured. Length and delay '
                           'are declared; traffic on this route is not observed.',
            'source': 'data/derived/canonical/peripheral.json',
        })

    # -- the internal edges: what the measurements are actually about --------
    innervation = ibm.get('innervation') or {}
    mapped = [m for m, entry in innervation.items()
              if any(r in LEVELS for r in entry[1])] if innervation else []
    segments = sorted({r for entry in innervation.values() for r in entry[1] if r in LEVELS},
                      key=LEVELS.index) if innervation else []
    ports = None
    if innervation and ibm.get('sensor_ports') is not None:
        ports = {'motor_alpha': len(innervation),
                 'motor_gamma': sum(1 for e in innervation.values() if e[2] > 0),
                 'plant_in': 3 * len(innervation),
                 'sensor_in': len(ibm['sensor_ports'])}
        ports['total'] = sum(ports.values())

    edges.append({
        'id': 'cortex.postcentral_to_precentral',
        'label': 'Sensory sites → motor sites, across the sheet',
        'kind': 'materialization', 'route_kind': 'cortical',
        'side': 'both', 'system': 'postcentral', 'target': 'precentral',
        'anchor_m': None,
        'fibres': [], 'muscle_ports': 0, 'receptor_ports': 0,
        'magnitude': 0.0005, 'magnitude_kind': 'measured',
        'magnitude_range': [0.0003, 0.0007],
        'measurement': 'With the drive entering postcentral and the command read '
                       'from precentral — two disjoint site populations — 0.03–0.07% '
                       'of the driven signal arrives at the readout. The ratio scales '
                       'linearly with association gain rather than compounding: the '
                       'signal crosses in one weak hop and never propagates.',
        'source': 'IBM-1 docs/LOG.md, 2026-09-09',
        'limitations': ['Severing the kernel changed the motor command by nothing, '
                        'because nothing was getting through to sever.'],
    })
    edges.append({
        'id': 'cortex.stance_correction',
        'label': 'Cortical sheet → all 98 muscle commands (stance policy)',
        'kind': 'materialization', 'route_kind': 'cortical',
        'side': 'both', 'system': 'cortical_sheet', 'target': 'cord',
        'fibres': [], 'muscle_ports': 98, 'receptor_ports': 0,
        'magnitude': 1.0, 'magnitude_kind': 'measured',
        'measurement': 'The one load-bearing path on this graph. Persistent E/I '
                       'dynamics own all 98 corrective commands: the intact arm '
                       'holds the pushed body at 0.17 mm final COM, the severed arm '
                       'falls at 2.90 s. Matched trials, same native checkpoint.',
        'source': 'docs/IBM_CORTICAL_STANCE.md',
        'limitations': ['Severing removes the dynamics, not the learning. A kernel '
                        'with its site rows permuted recovers the same push (1.7452 '
                        'against 1.7844 mm peak COM) and its severed arm falls at the '
                        'same 2.90 s, so this says the cortex is the motor owner and '
                        'says nothing about whether the trained content contributed.'],
    })
    edges.append({
        'id': 'cortex.trained_content',
        'label': 'Trained kernel content → stance recovery',
        'kind': 'materialization', 'route_kind': 'cortical',
        'side': 'both', 'system': 'cortical_sheet', 'target': 'cord',
        'fibres': [], 'muscle_ports': 0, 'receptor_ports': 0,
        'magnitude': 0.022, 'magnitude_kind': 'measured',
        'measurement': 'A site-permuted kernel — learning destroyed, every marginal '
                       'statistic preserved — recovers the same push at 1.7452 mm '
                       'peak COM against the trained kernel’s 1.7844 mm: a 2.2% '
                       'difference, and both severed arms fall on the same trajectory '
                       'at the same 2.90 s.',
        'source': 'docs/IBM_CURRICULUM16_KERNEL.md',
    })
    for name, label, gain, delay in ARCS:
        edges.append({
            'id': f'cord.arc.{name}',
            'label': label,
            'kind': 'arc', 'route_kind': 'spinal',
            'side': 'both', 'system': 'cord', 'target': 'cord',
            'fibres': [], 'muscle_ports': len(mapped), 'receptor_ports': 0,
            'declared_gain': gain, 'loop_delay_s': delay,
            'magnitude': None, 'magnitude_kind': 'declared',
            'measurement': f'Declared arc: gain {gain:+.2f} at a {delay * 1000:.0f} ms '
                           'loop delay. A parameter, not a flow — a live peak appears '
                           'here only while a body is running.',
            'source': 'IBM-1 ibm/processes/cord.py, read from ibm/processes/spinal.py',
        })

    anchors = _anchors(root, peripheral)
    systems = [dict(system, anchor_m=anchors.get(system['id']),
                    anchor_scope='Declared centroid in the canonical body frame; a '
                                 'place to land a leader line, not a claim about '
                                 'where a function lives.')
               for system in SYSTEMS]
    return {
        'schema': 'ihm.brain-graph.v1',
        'systems': systems,
        'edges': edges,
        'counts': {
            'routes': len(peripheral.get('nerves', [])),
            'routes_with_fibre_classes': joined,
            'muscle_bindings': len(peripheral.get('muscle_bindings', [])),
            'receptor_patches': len(peripheral.get('receptor_patches', [])),
            'relays': len(peripheral.get('relays', [])),
            'unsupported_muscles': len(peripheral.get('unsupported_muscles', [])),
            'ports': ports,
            'cord_muscles_declared': len(innervation) or None,
            'cord_muscles_mapped': len(mapped) or None,
            'cord_segments_declared': len(LEVELS),
            'cord_segments_used': len(segments) or None,
        },
        'routes_without_fibre_classes': unjoined,
        'ibm_declarations': {'available': ibm['available'],
                             'reason': ibm.get('reason'),
                             'root': ibm.get('root')},
        'fibre_velocity_m_s': velocity,
        'fibre_velocity_scope': peripheral.get('fibre_velocity_scope'),
        'route_contract': peripheral.get('route_contract'),
        'limitations': peripheral.get('limitations') or [],
        'magnitude_contract': {
            'measured': 'A measurement with its receipt named in source. Arrow width '
                        'is this number. A measured near-zero is drawn near-zero.',
            'live': 'Read out of the current embodied frame this tick.',
            'declared': 'A declared gain or latency. Not a flow; never drawn as width.',
            'unmeasured': 'No measurement exists for this edge. Not the same as a '
                          'measured zero and must not be drawn as one.',
        },
    }


def _finite(value):
    return value if isinstance(value, (int, float)) and value == value and abs(value) != float('inf') else None


def brain_state(frame=None):
    """Per-edge magnitude for this tick.

    With no live body there is nothing per-tick to report, so the measured
    static values stand and the payload says so — ``live: false`` and
    ``basis: 'measured-not-live'``. Nothing is animated to fill the gap.
    """
    edges = {}
    live = False
    neural = (frame or {}).get('neural') or {}
    controller = (frame or {}).get('controller') or neural.get('controller') or {}
    if frame and frame.get('schema') == 'ihm.embodied-frame.v1':
        live = True
        excitations = [abs(v) for v in (neural.get('motor_excitations') or {}).values()
                       if _finite(v) is not None]
        peak = max(excitations) if excitations else None
        if peak is not None:
            edges['cortex.stance_correction'] = {
                'magnitude': min(1.0, peak), 'magnitude_kind': 'live',
                'measurement': f'Live motor excitation peak {peak:.4g} over '
                               f'{len(excitations)} effectors this tick.',
            }
        correction = _finite((neural.get('cortical_stance') or {}).get('max_cortical_correction'))
        if correction is not None:
            edges['cortex.postcentral_to_precentral'] = {
                'magnitude': min(1.0, abs(correction)), 'magnitude_kind': 'live',
                'measurement': f'Live cortical correction peak {correction:.4g}. This is '
                               'the correction the sheet emits, not the fraction of '
                               'sensory drive that crossed it.',
            }
        for name, _label, _gain, _delay in ARCS:
            value = _finite((neural.get('arc_max') or {}).get(name))
            available = (controller.get('arc_availability') or {}).get(name)
            if value is None:
                edges[f'cord.arc.{name}'] = {
                    'magnitude': None, 'magnitude_kind': 'declared',
                    'measurement': 'Inactive for this controller.' if available is False
                                   else 'This controller does not report an arc peak.',
                }
            else:
                edges[f'cord.arc.{name}'] = {
                    'magnitude': min(1.0, abs(value)), 'magnitude_kind': 'live',
                    'measurement': f'Live arc peak {value:.4g} this tick.',
                }
    return {
        'schema': 'ihm.brain-state.v1',
        'live': live,
        'basis': 'live-frame' if live else 'measured-not-live',
        'note': 'Per-tick magnitudes read from the running body.' if live else
                'No body is running, so no edge has a live magnitude. The graph’s '
                'measured values stand as measured; nothing here is animated to '
                'stand in for a signal that is not flowing.',
        'time_s': (frame or {}).get('time_s') if live else None,
        'sequence': (frame or {}).get('sequence') if live else None,
        'controller': controller.get('kind') if live else None,
        'severed': bool(controller.get('sever')) if live else None,
        'edges': edges,
    }
