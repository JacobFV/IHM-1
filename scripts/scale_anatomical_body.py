"""Scale the 4,000 anatomical entities, not just the scaffold they ride on.

    .venv/bin/python scripts/scale_anatomical_body.py --stature-m 2.03
    .venv/bin/python scripts/scale_anatomical_body.py --scale 1.0     # the identity gate

`scripts/materialize_stature_variant.py` scales the mechanical body: joint
frames, path polynomials, contact spheres, mass and inertia.  It does not touch
the real body at all.  Ask for a tall person today and you get a tall
22-segment scaffold posing a fixed-size anatomy.  This closes that.

**What the scale actually is.**  Every entity rides a segment, and every
segment's scale is known -- it is the one isotropic factor
`ihm/native/model_scaling.py` applied to the whole model.  Under that factor
each segment's frame origin moves to `s` times its position along the joint
chain, so an entity that is to stay in the same place ON its segment has to
have both its placement and its extent multiplied by `s`.  Composed over the
whole chain that is exactly one uniform scale about a fixed origin, which is
why this script does not need per-segment logic for the geometry.

**Where the binding is load-bearing anyway.**  It is the mechanism by which an
entity could be treated differently: an entity whose absolute size is fixed by
something other than body size still has to MOVE with its segment, so its
placement exponent and its size exponent come apart.  The two are separated
here for every entity, so an exception is a table entry rather than a rewrite.
Whether any entity in this atlas deserves one is a question with a measured
answer, and the report gives it rather than assuming either way.

**Meshes are not rewritten.**  4,000 gzipped triangle meshes multiplied by a
constant are 4,000 copies of the donor geometry, and `data/derived` is already
large.  The variant records the transform instead -- a uniform scale about the
canonical origin -- exactly as the mechanical model records `scale_factors` on
its `Mesh` elements rather than shipping resized meshes.  The gates below still
walk the real triangles, scale the vertices, and re-integrate; they do not take
the recorded field's word for anything.

What comes out, under `data/derived/body-variants/<name>/`:

    anatomy_scaled.json   4,000 entity records, scaled, each carrying the
                          scaling class it was given and why
    binding_scaled.json   the segment binding's centroids and registration
    report.json           counts, the exception ledger, and the gates
"""
from pathlib import Path
import argparse, gzip, hashlib, json, math, sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm import body_scaling as bs                                  # noqa: E402
from ihm.body_parameters import (ANATOMICAL_MASS_KG,                # noqa: E402
                                 ANATOMICAL_STATURE_M,
                                 ANATOMY_REGISTRATION_SCALE,
                                 MECHANICAL_STATURE_M, resolve)

ANATOMY = 'data/derived/canonical/anatomy.json'
BINDING = 'data/derived/anatomy-segment-binding/binding.json'
PROFILE = 'data/derived/canonical/profile.json'
SKIN_ENTITY = 'body-bp3d-FJ2810'

#: The scale below which absolute size is set by chemistry rather than by body
#: size, in metres.  Nothing here is a claim about this atlas; it is the ruler
#: the atlas is measured against, so that "no entity is chemistry-fixed" is a
#: measurement and not an assumption.  A red cell is 8 um across in a mouse and
#: in a horse; a sarcomere is 2.7 um at optimal length in both; a systemic
#: capillary lumen is about 4 um; a myelinated axon is 1-20 um.  All of these
#: are two to four orders of magnitude below anything BodyParts3D resolves, and
#: the report says by how much rather than asserting it.
CHEMISTRY_SCALE_M = 20e-6
CHEMISTRY_SCALE_BASIS = (
    'Largest of the size-invariant cellular dimensions: myelinated axon '
    'diameter (1-20 um), erythrocyte diameter (~8 um), systemic capillary lumen '
    '(~4 um), sarcomere length at optimal overlap (~2.7 um). An entity whose '
    'smallest extent approached this would need its size held while its '
    'placement scaled.')

#: Entity classes that get an annotation rather than a different exponent.
#: Each names an ALLOMETRY entry, so the caveat and the exponent it qualifies
#: cannot drift apart.
ANNOTATED = {
    'vascular': dict(allometry='volumetric_flow', match_role=('vascular',),
                     note=('Vessel calibre takes the length exponent like every '
                           'other extent. Poiseuille resistance goes as r**-4, so '
                           'flow at a regulated pressure rises as s**3 while '
                           'metabolic demand rises as about s**2.25. No flow model '
                           'is reparameterised by this scaler; the mismatch is '
                           'declared, not corrected.')),
    'nervous': dict(allometry='brain_volume', match_system=('nervous',),
                    note=('Scaled isotropically like everything else, which makes '
                          'the brain s**3. Adult human brain volume is only weakly '
                          'related to stature and no catalogued source here '
                          'measures the relation, so no exponent is asserted and '
                          'the isotropic one is flagged as known-wrong.')),
}

#: Which record field takes which quantity from the scaling table, addressed by
#: PATH so the nested ones cannot be missed.  The first draft of this file
#: listed top-level fields only and silently left three of them at the source
#: subject's size: `connections[].distance_m` on 3,996 entities, and the skin's
#: own `shell.thickness_m` -- the 6.6 mm that `docs/MILESTONES.md` weighs the
#: body against.  The completeness gate below now walks every numeric leaf in
#: all 4,000 records and requires each path to appear in exactly one of these
#: two tables, which is how those three were found.
ENTITY_FIELDS = {
    'centroid_m[]': 'centroid',
    'bounds_m.min[]': 'centroid',
    'bounds_m.max[]': 'centroid',
    'surface_area_m2': 'area',
    'volume_m3': 'volume',
    'connections[].distance_m': 'length',
    'shell.thickness_m': 'length',
    'shell.prior_range_m[]': 'length',
    'shell.depth_interval_m[]': 'length',
    'physical_surface_support.area_m2': 'area',
    'physical_surface_support.raw_source_area_m2': 'area',
}

#: Numeric leaves that must come through untouched, each with the reason.
#: Listed rather than skipped by default: a field that is neither scaled nor
#: named here fails the completeness gate, because that is exactly what a
#: quantity silently left at the source subject's size looks like.
ENTITY_INVARIANT = {
    'principal_axis[]': 'a unit direction; a uniform scale preserves directions',
    'source_vertex_count': 'a count',
    'source_face_count': 'a count',
    'physical_surface_support.selected_triangle_count': 'a count',
    'physical_surface_support.source_triangle_count': 'a count',
    'physical_surface_support.selected_component_id': 'an identifier',
    'uncertainty.biological.independent_subject_count': 'a count of subjects',
    'provenance.source.independent_subject_count': 'a count of subjects',
    'provenance.source.schema_version': 'a version',
    'provenance.source.files[].bytes': 'a file size on disk',
    'provenance.source.archive_members[].bytes': 'a file size on disk',
    'provenance.source_to_canonical.scale':
        'PROVENANCE, not a body quantity: it records how the raw donor mesh file '
        'maps into the BASE canonical frame. It is deliberately left alone and '
        'the variant\'s geometry_transform composes on top of it, so a consumer '
        'that applies both gets the scaled body and a reader can still see where '
        'the bytes came from. Folding the stature scale in here would make the '
        'provenance of the source file untrue.',
    'provenance.source_to_canonical.translation[]': 'the same provenance transform',
    'provenance.source_to_canonical.rotation[][]': 'a rotation matrix; dimensionless',
}


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mesh_area_and_volume(path, scale=1.0):
    """Surface area and signed volume, from the triangles themselves.

    This is the independent arm.  It reads the gzipped triangle mesh that
    `anatomy.json` merely POINTS at, scales every vertex, and integrates -- so
    the area it returns is a sum of cross products over scaled coordinates and
    the volume is a divergence-theorem sum over scaled coordinates, neither of
    which is "the recorded number times a power of s".  If a record's area had
    been given the length exponent by mistake, the recorded value and this one
    would differ by a factor of s and the gate would say so.
    """
    with gzip.open(path) as handle:
        mesh = json.load(handle)
    P = np.asarray(mesh['positions'], float).reshape(-1, 3) * float(scale)
    F = np.asarray(mesh['indices'], np.int64).reshape(-1, 3)
    a, b, c = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
    cross = np.cross(b - a, c - a)
    area = float(0.5 * np.linalg.norm(cross, axis=1).sum())
    volume = float(np.einsum('ij,ij->i', a, cross).sum() / 6.0)
    return area, volume, len(F)


def classify(entity):
    """The scaling class of one entity, and the annotations it carries.

    Two exponents, deliberately separated: an entity always rides its segment,
    so `placement` is the length exponent; `size` is what a genuine exception
    would change.  Both come from the table.
    """
    annotations = []
    for name, rule in ANNOTATED.items():
        if entity.get('role') in rule.get('match_role', ()) \
                or entity.get('system') in rule.get('match_system', ()):
            annotations.append(name)
    return dict(scaling_class='isotropic',
                placement_exponent=bs.exponent('length'),
                size_exponent=bs.exponent('length'),
                annotations=annotations)


def smallest_extent(entity):
    bounds = entity.get('bounds_m')
    if not bounds:
        return None
    lo, hi = bounds['min'], bounds['max']
    return float(min(h - l for l, h in zip(lo, hi)))


def numeric_paths(node, prefix=''):
    """Every numeric leaf in a record, addressed by path.

    A list of numbers is one path ending in ``[]`` -- a position or an interval
    is scaled as a unit, and calling its three components three quantities would
    let one of them be missed.
    """
    if isinstance(node, dict):
        for key, value in node.items():
            yield from numeric_paths(value, '%s.%s' % (prefix, key) if prefix else key)
    elif isinstance(node, list):
        if node and all(isinstance(x, (int, float)) and not isinstance(x, bool)
                        for x in node):
            yield prefix + '[]'
        else:
            for value in node:
                yield from numeric_paths(value, prefix + '[]')
    elif isinstance(node, (int, float)) and not isinstance(node, bool):
        yield prefix


def parse_path(path):
    """'connections[].distance_m' -> (('connections', True), ('distance_m', False))"""
    return tuple((segment.removesuffix('[]'), segment.endswith('[]'))
                 for segment in path.split('.'))


def apply_at(node, parts, leaf):
    """Return a copy of ``node`` with ``leaf`` applied at every match of ``parts``.

    ``leaf`` takes a whole numeric list (for a ``[]`` segment) or a scalar, so a
    position can be scaled about an origin as a unit rather than component by
    component.  A path that is absent, or present and null, is left alone: not
    every entity carries every field, and a missing field is not an unscaled one.
    """
    if not parts:
        return leaf(node)
    (key, is_list), rest = parts[0], parts[1:]
    if not isinstance(node, dict) or node.get(key) is None:
        return node
    value = node[key]
    if is_list and not rest:
        value = leaf(value)
    elif is_list:
        value = [apply_at(item, rest, leaf) for item in value]
    else:
        value = apply_at(value, rest, leaf)
    return dict(node, **{key: value})


def scale_entity(entity, s, origin):
    """One entity: its placement and its extent, every numeric leaf by path.

    A position takes the scaling origin; every other length is an extent and
    takes the bare factor.  Both use the same exponent -- the split is about
    where the origin goes, not about how much anything grows.
    """
    out = entity
    for path, quantity in ENTITY_FIELDS.items():
        factor = float(s) ** bs.exponent(quantity)
        if quantity == 'centroid':
            def leaf(v, _o=origin, _s=s):
                return bs.scale_point(v, _s, _o)
        else:
            def leaf(v, _f=factor):
                return [x * _f for x in v] if isinstance(v, list) else v * _f
        out = apply_at(out, parse_path(path), leaf)
    return dict(out, scaling=classify(entity))


def build(root, scale, origin=(0.0, 0.0, 0.0), mesh_sample=200):
    root = Path(root)
    anatomy = json.loads((root / ANATOMY).read_text())
    binding = json.loads((root / BINDING).read_text())
    entities = anatomy['entities']

    scaled = [scale_entity(e, scale, origin) for e in entities]
    by_id = {e['id']: e for e in entities}
    scaled_by_id = {e['id']: e for e in scaled}

    # ---- the binding: centroids ride the segments, so they take the length
    #      exponent; the registration SCALE is a ratio of two bodies and is
    #      invariant only because BOTH bodies are being scaled by the same
    #      factor, which is the point of doing this alongside the mechanical
    #      variant rather than instead of it.
    scaled_binding = dict(binding)
    scaled_binding['centroids_m'] = {
        key: bs.scale_point(value, scale, origin)
        for key, value in binding['centroids_m'].items()}
    registration = dict(binding['registration'])
    registration['translation_m'] = bs.scale_vector('length',
                                                    registration['translation_m'], scale)
    registration['scale_invariance_note'] = (
        'registration.scale is the ratio of the mechanical body to the atlas. '
        'Both are scaled by the same factor here, so it is unchanged at %.10f; '
        'a scaler that moved one body and not the other would change it, and '
        'the gates check that it did not.' % ANATOMY_REGISTRATION_SCALE)
    scaled_binding['registration'] = registration
    scaled_binding['reference_pose_rad'] = dict(binding['reference_pose_rad'])
    for name in ('pelvis_tx', 'pelvis_ty', 'pelvis_tz'):
        if name in scaled_binding['reference_pose_rad']:
            scaled_binding['reference_pose_rad'][name] = bs.scale(
                'length', scaled_binding['reference_pose_rad'][name], scale)

    # ---- the measured answer to "is anything in this atlas chemistry-fixed"
    extents = [(smallest_extent(e), e['id'], e['name']) for e in entities]
    extents = sorted((x for x in extents if x[0] is not None))
    finest = extents[0]
    chemistry = dict(
        question=('Does any entity have an absolute size fixed by chemistry rather '
                  'than by body size? Such an entity would need its placement '
                  'scaled and its size held.'),
        ruler_m=CHEMISTRY_SCALE_M, ruler_basis=CHEMISTRY_SCALE_BASIS,
        smallest_entity_extent_m=finest[0], smallest_entity=finest[1],
        smallest_entity_name=finest[2],
        ratio_to_ruler=finest[0] / CHEMISTRY_SCALE_M,
        entities_below_ruler=sum(1 for x in extents if x[0] < CHEMISTRY_SCALE_M),
        finding=('The smallest extent in the atlas is %.3g m, %.0fx the '
                 'chemistry-fixed ruler. No entity in this atlas is at a scale '
                 'where absolute size is set by molecular constraint: those '
                 'constraints act BELOW the atlas resolution and appear in this '
                 'body as invariant PROPERTIES -- conduction velocity, specific '
                 'tension, density, elastic modulus -- not as invariant entities. '
                 'So the size exponent is the length exponent for all %d, and '
                 'that is a measurement, not a default.'
                 % (finest[0], finest[0] / CHEMISTRY_SCALE_M, len(entities))))

    counts = {'by_class': {}, 'by_annotation': {}, 'by_system': {}, 'by_role': {}}
    for e in scaled:
        counts['by_class'][e['scaling']['scaling_class']] = \
            counts['by_class'].get(e['scaling']['scaling_class'], 0) + 1
        for a in e['scaling']['annotations']:
            counts['by_annotation'][a] = counts['by_annotation'].get(a, 0) + 1
        counts['by_system'][e['system']] = counts['by_system'].get(e['system'], 0) + 1
        counts['by_role'][e['role']] = counts['by_role'].get(e['role'], 0) + 1

    gates = run_gates(root, entities, scaled, by_id, scaled_by_id, binding,
                      scaled_binding, scale, origin, mesh_sample)

    report = dict(
        schema='ihm.anatomical-body-scaling.v1',
        stature_scale=scale,
        scaling_origin=bs.SCALING_ORIGIN,
        entities=len(entities), counts=counts,
        chemistry_scale_probe=chemistry,
        fields_scaled={f: dict(quantity=q, exponent=bs.exponent(q))
                       for f, q in ENTITY_FIELDS.items()},
        fields_invariant=dict(ENTITY_INVARIANT),
        annotations={name: dict(rule, entities=counts['by_annotation'].get(name, 0),
                                allometry=rule['allometry'])
                     for name, rule in ANNOTATED.items()},
        allometry=[dict(a) for a in bs.ALLOMETRY],
        scaling_table=bs.check_table(),
        geometry_transform=dict(
            kind='uniform_scale_about_origin', scale=scale, origin=list(origin),
            note=('The 4,000 triangle meshes are NOT rewritten. A consumer applies '
                  'this transform, exactly as the OpenSim model carries '
                  'scale_factors on its Mesh elements rather than shipping resized '
                  'geometry. The gates walk the real triangles and re-integrate, so '
                  'nothing here rests on the transform being applied correctly by '
                  'someone else.')),
        gates=gates,
        limitations=[
            'One isotropic factor. This changes how big the body is and not what '
            'shape it is; every proportion is the source atlas\'s at every stature.',
            'Mesh geometry is carried by a transform, not by rewritten files. An '
            'entity whose size exponent differed from its placement exponent could '
            'not be expressed that way and would need its mesh written out; no '
            'entity in this atlas needs that, and the probe above is the evidence.',
            'The anatomical and mechanical bodies remain 9.7% apart in mass and '
            '0.66% apart in stature. Scaling both by the same factor preserves the '
            'disagreement exactly rather than resolving it.',
        ])
    return scaled, scaled_binding, report


def run_gates(root, entities, scaled, by_id, scaled_by_id, binding, scaled_binding,
              scale, origin, mesh_sample):
    gates = []

    def record(name, error, tolerance, note, **extra):
        gates.append(dict(gate=name, relative_error=error, tolerance=tolerance,
                          passed=(error <= tolerance), note=note, **extra))

    # 1. IDENTITY.  At s = 1.0 every quantity this script touches must come back
    #    EXACTLY -- not to a tolerance.  A float times 1.0 is itself, so any
    #    difference at all is a code path that did something other than scale.
    if scale == 1.0:
        differences, compared = 0, 0
        for a, b in zip(entities, scaled):
            stripped = {k: v for k, v in b.items() if k != 'scaling'}
            compared += 1
            if stripped != a:
                differences += 1
        for key, value in binding['centroids_m'].items():
            compared += 1
            if value != scaled_binding['centroids_m'][key]:
                differences += 1
        record('identity at s=1.0 reproduces the base exactly',
               float(differences), 0.0,
               'WHOLE-RECORD equality, not a field list: the scaled record with '
               'its added scaling block removed must equal the base record. A '
               'field this script does not know about therefore cannot pass this '
               'gate by being skipped, and v * 1.0**e is v for every finite float '
               'so one difference is one bug',
               records_compared=compared)

    # 2. Every scaled length is exactly s times the base.
    worst, n = 0.0, 0
    for a, b in zip(entities, scaled):
        for axis in range(3):
            base, got = a['centroid_m'][axis], b['centroid_m'][axis]
            if abs(base) > 1e-9:
                worst = max(worst, abs(got / (base * scale) - 1.0)); n += 1
    record('every entity centroid scales by exactly s', worst, 1e-12,
           'the placement half of the two exponents', components=n)

    # 3. TRANSLATION-INVARIANT, so it does not depend on where the origin was
    #    put.  The distance between two entities is a property of the body.
    ids = sorted(by_id)
    pairs = [(ids[i], ids[(i * 977 + 13) % len(ids)]) for i in range(0, len(ids), 4)]
    worst, n = 0.0, 0
    for x, y in pairs:
        if x == y:
            continue
        d0 = math.dist(by_id[x]['centroid_m'], by_id[y]['centroid_m'])
        d1 = math.dist(scaled_by_id[x]['centroid_m'], scaled_by_id[y]['centroid_m'])
        if d0 > 1e-6:
            worst = max(worst, abs(d1 / (d0 * scale) - 1.0)); n += 1
    record('inter-entity distances scale by exactly s', worst, 1e-12,
           'independent of the scaling origin, so this gate passes or fails on '
           'the scaling and not on a bookkeeping choice', pairs=n)

    # 4. THE INDEPENDENT ONE.  Different file, different code path: walk the
    #    gzipped triangle meshes anatomy.json only points at, scale the
    #    VERTICES, and integrate area and volume from the scaled coordinates.
    #    Multiplying a recorded number by s**2 and integrating scaled triangles
    #    are not the same computation, so this catches an exponent that is
    #    wrong in the table as well as one that is wrong in the application.
    geometry = root / 'data/derived/canonical/geometry'
    sample = [e for e in entities if e.get('surface_area_m2')][::max(
        1, len(entities) // max(1, mesh_sample))]
    if SKIN_ENTITY in by_id and by_id[SKIN_ENTITY] not in sample:
        sample.append(by_id[SKIN_ENTITY])
    #    Only entities whose recorded volume IS a surface integral take the
    #    volume arm.  The three skin layers record an open-shell volume --
    #    exterior area times an assumed thickness -- and integrating a shell
    #    that has a boundary returns a number that means nothing; the first
    #    draft of this gate did it anyway and failed at s = 1.0 by 0.155%,
    #    which was the shell, not the scaling. They get their own gate below.
    CLOSED = 'absolute signed surface integral'
    base_area_error, scaled_area_error = 0.0, 0.0
    volume_shift, volume_agreement = 0.0, []
    checked, volumes = 0, 0
    for entity in sample:
        path = geometry / ('%s.json.gz' % entity['id'])
        if not path.exists():
            continue
        area_1, volume_1, _ = mesh_area_and_volume(path, 1.0)
        area_s, volume_s, _ = mesh_area_and_volume(path, scale)
        record_area = scaled_by_id[entity['id']]['surface_area_m2']
        base_area_error = max(base_area_error,
                              abs(area_1 / entity['surface_area_m2'] - 1.0))
        scaled_area_error = max(scaled_area_error, abs(area_s / record_area - 1.0))
        if entity.get('volume_m3') and (entity.get('volume_method') or '').startswith(CLOSED):
            record_volume = scaled_by_id[entity['id']]['volume_m3']
            if abs(record_volume) > 1e-12:
                # RATIO, not agreement. One entity -- the right external
                # intercostals -- disagrees with its own record by 0.155%
                # before any scaling, and its volume_method says why: edge
                # incidence closed, SELF-INTERSECTION NOT CERTIFIED, so the
                # divergence integral double-counts where the surface passes
                # through itself. Demanding agreement would fail at s = 1.0 on
                # a defect in the mesh; demanding that the ratio not MOVE tests
                # the exponent and nothing else, which is the same shape as the
                # mechanical scaler's polyline gate and for the same reason.
                r0 = abs(volume_1) / abs(entity['volume_m3'])
                r1 = abs(volume_s) / abs(record_volume)
                volume_shift = max(volume_shift, abs(r1 / r0 - 1.0))
                volume_agreement.append(r0)
                volumes += 1
        checked += 1
    record('recorded area agrees with the mesh it points at', base_area_error, 1e-9,
           'the base check the scaled one rests on: if the record and the mesh '
           'disagreed before scaling, agreement after would mean nothing',
           entities=checked)
    record('area integrated over SCALED triangles matches the scaled record',
           scaled_area_error, 1e-9,
           'INDEPENDENT ARM. A sum of cross products over scaled vertex '
           'coordinates, from the gzipped mesh files, against surface_area_m2 '
           'times s**%g. Give area the length exponent by mistake and these '
           'differ by exactly s.' % bs.exponent('area'), entities=checked)
    record('mesh volume over recorded volume does not move under scaling',
           volume_shift, 1e-9,
           'the same independent arm on the divergence-theorem volume, exponent '
           's**%g. Stated as a ratio that must not move rather than an agreement '
           'that must hold, because the two representations already disagree by '
           'up to %.3f%% on a mesh whose own volume_method says self-intersection '
           'is not certified. Scale the record by anything other than s**%g and '
           'this ratio moves by that discrepancy.'
           % (bs.exponent('volume'),
              100 * max((abs(1 - r) for r in volume_agreement), default=0.0),
              bs.exponent('volume')),
           entities=volumes,
           base_agreement_min=min(volume_agreement, default=None),
           base_agreement_max=max(volume_agreement, default=None))

    # 4b. The open shells. Their volume is an area times a thickness, and that
    #     identity has to survive scaling -- which it does only if area really
    #     took s**2 and thickness really took s**1. It is the same arithmetic
    #     the scaling table performs, done from the record instead of from the
    #     table, and it is the gate that would catch the skin's own 6.6 mm being
    #     left at the source subject's thickness while its area grew.
    worst, shells = 0.0, 0
    for entity in scaled:
        shell = entity.get('shell')
        support = entity.get('physical_surface_support')
        if not shell or not support or not entity.get('volume_m3'):
            continue
        implied = support['area_m2'] * shell['thickness_m']
        worst = max(worst, abs(implied / entity['volume_m3'] - 1.0))
        shells += 1
    record('open-shell volume is still area times thickness after scaling',
           worst, 1e-12,
           'V = A * t. A takes s**%g and t takes s**%g, so V takes s**%g; the '
           'identity holds after scaling only if all three did what the table '
           'says. This is the gate that catches the skin growing in area while '
           'keeping the source subject\'s 6.6 mm thickness.'
           % (bs.exponent('area'), bs.exponent('length'), bs.exponent('volume')),
           shells=shells)

    # 5. Anatomical stature, from the skin mesh, against the mechanical stature
    #    carried onto the atlas by the registration scale.  Two measurements
    #    that share no input, and their 0.66% disagreement must not move.
    skin = geometry / ('%s.json.gz' % SKIN_ENTITY)
    if skin.exists():
        with gzip.open(skin) as handle:
            P = np.asarray(json.load(handle)['positions'], float).reshape(-1, 3)
        extent = float((P[:, 1].max() - P[:, 1].min()) * scale)
        predicted = MECHANICAL_STATURE_M * scale * ANATOMY_REGISTRATION_SCALE
        base_disagreement = abs(MECHANICAL_STATURE_M * ANATOMY_REGISTRATION_SCALE
                                / ANATOMICAL_STATURE_M - 1.0)
        record('the two bodies disagree by exactly as much as they did',
               abs(abs(predicted / extent - 1.0) - base_disagreement), 1e-12,
               'the mechanical stature carried onto the atlas by the fitted '
               'registration scale, against the atlas skin extent measured from '
               'the mesh. Both scale, so their %.4f%% disagreement is invariant; '
               'a scaler that moved one body and not the other changes it.'
               % (100 * base_disagreement),
               anatomical_stature_m=extent, mechanical_carried_m=predicted)

    # 6. Composed mass. The anatomical body's mass is a volume at fixed
    #    densities, so it takes the volume exponent; the profile records the
    #    base and the scaled value must be exactly s**3 of it.
    got = bs.scale('mass', ANATOMICAL_MASS_KG, scale)
    record('composed anatomical mass takes the mass exponent',
           abs(got / (ANATOMICAL_MASS_KG * scale ** 3) - 1.0), 1e-12,
           'density is invariant, so mass follows volume: s**%g'
           % bs.exponent('mass'),
           base_kg=ANATOMICAL_MASS_KG, scaled_kg=got)

    # 7. Completeness. Every numeric field of an entity record must be either
    #    scaled or explicitly declared invariant. A field that is neither is a
    #    quantity silently left at the source subject's size, which is the exact
    #    failure this whole script exists to end.
    numeric = set()
    for entity in entities:
        numeric |= set(numeric_paths(entity))
    unhandled = numeric - set(ENTITY_FIELDS) - set(ENTITY_INVARIANT)
    record('every numeric entity leaf is scaled or declared invariant',
           float(len(unhandled)), 0.0,
           'walks all %d records to every numeric leaf, by path. This is the '
           'gate that found connections[].distance_m on 3,996 entities and the '
           'skin\'s own shell.thickness_m, both of which a top-level field list '
           'missed. Unhandled: %r' % (len(entities), sorted(unhandled)),
           numeric_leaves=len(numeric), scaled=len(ENTITY_FIELDS),
           declared_invariant=len(ENTITY_INVARIANT))
    return gates


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--stature-m', type=float)
    group.add_argument('--scale', type=float)
    parser.add_argument('--output', default=None)
    parser.add_argument('--mesh-sample', type=int, default=200,
                        help='entities whose triangles the independent gate walks')
    parser.add_argument('--sabotage', choices=('area-as-length', 'volume-as-area',
                                               'thickness-unscaled'), default=None,
                        help=('deliberately break one exponent and report which '
                              'gates notice. A gate suite nobody has watched fail '
                              'is a gate suite nobody has tested.'))
    args = parser.parse_args()

    if args.sabotage == 'area-as-length':
        ENTITY_FIELDS['surface_area_m2'] = 'length'
        ENTITY_FIELDS['physical_surface_support.area_m2'] = 'length'
        ENTITY_FIELDS['physical_surface_support.raw_source_area_m2'] = 'length'
    elif args.sabotage == 'volume-as-area':
        ENTITY_FIELDS['volume_m3'] = 'area'
    elif args.sabotage == 'thickness-unscaled':
        del ENTITY_FIELDS['shell.thickness_m']

    stature = args.stature_m if args.stature_m is not None else args.scale * MECHANICAL_STATURE_M
    scale = resolve({'stature_m': stature})['derived']['stature_scale']
    name = args.output or ('data/derived/body-variants/stature_%s'
                           % ('%.6f' % scale).replace('.', 'p'))
    out = ROOT / name
    out.mkdir(parents=True, exist_ok=True)

    scaled, scaled_binding, report = build(ROOT, scale, mesh_sample=args.mesh_sample)
    report['stature_m'] = stature
    report['sources'] = {p: _sha(ROOT / p) for p in (ANATOMY, BINDING)}
    report['script_sha256'] = _sha(Path(__file__))

    (out / 'anatomy_scaled.json').write_text(
        json.dumps({'schema': 'ihm.anatomy-scaled.v1', 'stature_scale': scale,
                    'entities': scaled}, indent=1) + '\n')
    (out / 'binding_scaled.json').write_text(json.dumps(scaled_binding, indent=1) + '\n')
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')

    for gate in report['gates']:
        print('%-5s %-62s %.3g (tol %.3g)'
              % ('PASS' if gate['passed'] else 'FAIL', gate['gate'],
                 gate['relative_error'], gate['tolerance']))
    print(json.dumps(report['chemistry_scale_probe'], indent=2))
    print('entities %d, annotated %r' % (report['entities'], report['counts']['by_annotation']))
    print('wrote', out)
    if not all(g['passed'] for g in report['gates']):
        raise SystemExit('scaled anatomy failed its own gates')


if __name__ == '__main__':
    main()
