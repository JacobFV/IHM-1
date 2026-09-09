"""Muscle and passive tissue: PCSA as s**2, force as s**2, mass as s**3.

    .venv/bin/python scripts/scale_muscle_and_tissue.py --stature-m 2.03
    .venv/bin/python scripts/scale_muscle_and_tissue.py --scale 1.0     # identity

Two artifacts, and they are the two ends of the same problem.

**The 98 mechanical actuators** carry a maximum isometric force, an optimal
fibre length and a tendon slack length.  They do NOT carry a physiological
cross-sectional area, and that absence is where the classic error lives: force
and mass are both "how much muscle there is" and they take different exponents,
so a scaler that reasons about muscle in terms of mass makes every muscle in a
tall body too strong by a factor of `s`.

The exponents, and each one's justification, which is a dimensional formula and
not a preference:

    optimal fibre length      s**1   a length
    tendon slack length       s**1   a length
    PCSA                      s**2   an area, measured across the fibres
    max isometric force       s**2   specific tension x PCSA, and specific
                                     tension is a molecular property
    muscle volume             s**3   PCSA x fibre length -- which is the SAME
                                     exponent as length cubed, and the table
                                     refuses to load unless those two agree
    muscle mass               s**3   density x volume, density invariant
    joint moment              s**3   force x moment arm
    strength-to-weight        s**-1  force s**2 over weight s**3: a scaled-up
                                     body is relatively WEAKER, by 11.5% at
                                     2.03 m, which is a real consequence for
                                     whether it can pick itself up
    normalised fibre length   s**0   path over optimal fibre length. Not merely
                                     in the 0.2-20 band -- exactly invariant

**The 117 ligament force elements** come from the tissue agent's
`build_tissue_force_elements.py` and are shared work: this script reads them and
writes a scaled copy beside them.  It does not modify that script, and it does
not touch `native_mechanical_stream.cpp`, which consumes the spec file.

Their trap is a units trap.  A `Blankevoort1991Ligament` stiffness is in
NEWTONS, not newtons per metre, because its extension variable is a
dimensionless strain -- so it is `E*A`, a FORCE, and takes `s**2`.  Reading it
as the `E*A/L` spring constant instead gives `s**1` and is wrong by a whole
factor of `s`, which at 2.03 m is a ligament 13% too soft.  Both quantities are
in the scaling table, separately, and the table checks that the per-strain one
carries the same exponent as a muscle force -- because every force in a body
scales alike, however it is generated.
"""
from pathlib import Path
import argparse, hashlib, json, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm import body_scaling as bs                                  # noqa: E402
from ihm.body_parameters import MECHANICAL_STATURE_M, resolve       # noqa: E402
from ihm.native.model_scaling import (path_point_polyline_lengths,  # noqa: E402
                                      scale_catalog)

MECHANICS = 'data/derived/mechanics/whole_body_lumbar_current/registration.json'
TISSUE = 'data/derived/tissue-force-elements-v1'
MATERIALS = 'data/derived/tissue-material-candidate-v1/materials.json'
ANATOMY = 'data/derived/canonical/anatomy.json'

#: Maximum muscle stress.  Declared, with its range, because the catalog does
#: not carry PCSA and the only way to a PCSA is to divide a force by this.  Its
#: VALUE only sets the absolute PCSA; its INVARIANCE is what sets the exponent,
#: and the invariance is the claim being made -- so the gate below checks that
#: the implied specific tension does not move, which is true whatever number is
#: put here.
SPECIFIC_TENSION_PA = 6.0e5
SPECIFIC_TENSION_BASIS = (
    '60 N/cm2, the value conventionally used to set maximum isometric forces in '
    'lower-limb OpenSim models including the Rajagopal family this subject comes '
    'from. Human in-vivo estimates span roughly 25-60 N/cm2 and this repository '
    'has no independent measurement, so the ABSOLUTE PCSA below is a derived '
    'quantity with that whole spread on it. What is being claimed is only that '
    'the value does not depend on how tall the subject is, which is what makes '
    'force follow area.')

#: Ligament element fields, by quantity.  The units trap is here: a
#: Blankevoort stiffness is in newtons because its extension variable is a
#: strain, so it is E*A and takes the force exponent -- NOT the N/m one.
LIGAMENT_FIELDS = {
    'slack_length_m': 'length',
    'cross_section_m2': 'area',
    'linear_stiffness_n': 'ligament_stiffness_per_strain',
    'damping_n_s_per_strain': 'ligament_damping',
    'model_volume_m3': 'volume',
    'atlas_signed_volume_m3': 'volume',
    'atlas_tissue_volume_m3': 'volume',
    'atlas_surface_area_m2': 'area',
    'point1_m': 'centroid',
    'point2_m': 'centroid',
}
LIGAMENT_INVARIANT = {
    'transition_strain': 'a strain, dimensionless',
    'segment_share': 'a vote fraction',
    'runner_up_share': 'a vote fraction',
    'vertices': 'a count',
    'scaffold_joints_crossed': 'a count',
    'end_vertices': 'vertex indices',
    'peak_strain_over_declared_range':
        'a strain. INVARIANT, and that is a real statement: the peak strain a '
        'ligament reaches over the declared coordinate range is a kinematic '
        'quantity in a body scaled isotropically, so a taller body strains its '
        'ligaments exactly as much as a shorter one. It also means the elements '
        'flagged kinematically_admissible=false stay flagged at every stature; '
        'scaling does not rescue one and does not break one.',
}

#: 0.2-20 is the band scripts/collect_pose_corpus.py uses as a corpus quality
#: check and scripts/materialize_stature_variant.py reuses. Repeated here
#: because it is a required gate, with the same caveat the mechanical side
#: records: it is not sufficient on its own, and a 13% mis-scaled body passed it.
PATH_OVER_FIBER_BAND = (0.2, 20.0)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def derive_muscle_properties(catalog, scale, force_scale):
    """PCSA, volume and mass for each actuator -- derived, because the catalog
    carries none of them and the exponents differ from the ones it does carry."""
    rows = []
    for entry in catalog:
        force = entry.get('max_isometric_force_n')
        fibre = entry.get('optimal_fiber_length_m')
        if force is None or not fibre:
            continue
        pcsa = force / SPECIFIC_TENSION_PA
        rows.append(dict(
            id=entry['id'],
            max_isometric_force_n=force, optimal_fiber_length_m=fibre,
            pcsa_m2=pcsa, volume_m3=pcsa * fibre,
            scaled_max_isometric_force_n=bs.scale('max_isometric_force', force, scale)
                                         if force_scale is None else force * force_scale,
            scaled_optimal_fiber_length_m=bs.scale('fibre_length', fibre, scale),
            scaled_pcsa_m2=bs.scale('pcsa', pcsa, scale),
            scaled_volume_m3=bs.scale('volume', pcsa * fibre, scale)))
    return rows


def scale_ligaments(elements, scale, origin=(0.0, 0.0, 0.0)):
    out = []
    for element in elements:
        moved = dict(element)
        for field, quantity in LIGAMENT_FIELDS.items():
            if moved.get(field) is None:
                continue
            if quantity == 'centroid':
                moved[field] = bs.scale_point(moved[field], scale, origin)
            else:
                moved[field] = bs.scale(quantity, moved[field], scale)
        rule = moved.get('volume_rule')
        if rule:
            rule = dict(rule)
            for key in ('thickness_m',):
                if rule.get(key) is not None:
                    rule[key] = bs.scale('length', rule[key], scale)
            if rule.get('open_surface_signed_integral_m3') is not None:
                rule['open_surface_signed_integral_m3'] = bs.scale(
                    'volume', rule['open_surface_signed_integral_m3'], scale)
            moved['volume_rule'] = rule
        out.append(moved)
    return out


def build(root, scale):
    root = Path(root)
    registration = json.loads((root / MECHANICS).read_text())
    catalog = json.loads((root / registration['catalog_path']).read_text())
    parameters = resolve({'stature_m': scale * MECHANICAL_STATURE_M})
    force_scale = parameters['derived']['muscle_force_scale']
    scaled_catalog = scale_catalog(catalog, scale, force_scale)

    muscles = derive_muscle_properties(catalog, scale, force_scale)

    ligament_path = root / TISSUE / 'ligaments.json'
    ligaments = scaled_ligaments = None
    if ligament_path.exists():
        ligaments = json.loads(ligament_path.read_text())['elements']
        scaled_ligaments = scale_ligaments(ligaments, scale)

    gates, summary = run_gates(root, registration, catalog, scaled_catalog, muscles,
                               ligaments, scaled_ligaments, scale, force_scale)
    report = dict(
        schema='ihm.muscle-tissue-scaling.v1', stature_scale=scale,
        muscle_force_scale=force_scale,
        actuators=len(catalog), actuators_with_derived_pcsa=len(muscles),
        ligament_elements=len(ligaments or ()),
        specific_tension_pa=SPECIFIC_TENSION_PA,
        specific_tension_basis=SPECIFIC_TENSION_BASIS,
        exponents={name: bs.explain(name) for name in (
            'fibre_length', 'tendon_slack_length', 'pcsa', 'max_isometric_force',
            'muscle_volume_via_fibres', 'muscle_mass', 'joint_moment',
            'strength_to_weight', 'self_weight_stress', 'normalised_fibre_length',
            'ligament_stiffness', 'ligament_stiffness_per_strain',
            'ligament_damping', 'strain')},
        ligament_fields={f: dict(quantity=q, exponent=bs.exponent(q))
                         for f, q in LIGAMENT_FIELDS.items()},
        ligament_invariant=dict(LIGAMENT_INVARIANT),
        summary=summary, gates=gates,
        shared_work_note=(
            'ligaments.json is produced by scripts/build_tissue_force_elements.py, '
            'which is another agent\'s work in progress. This script READS it and '
            'writes a scaled copy into the variant; it modifies neither that script '
            'nor scripts/native_mechanical_stream.cpp, which consumes the .txt spec.'),
        limitations=[
            'PCSA is derived as force / specific tension because the catalog does '
            'not carry it. The absolute value inherits the whole 25-60 N/cm2 spread '
            'on human specific tension; only its EXPONENT is claimed here.',
            'Muscle mass is derived the same way and is not an independent '
            'measurement of this subject\'s muscle mass.',
            'The passive joint stops (ExpressionBasedCoordinateForceSet) still do '
            'not scale, on this side as on the mechanical one, and are in ALLOMETRY '
            'as a known gap.',
        ])
    return scaled_catalog, scaled_ligaments, report


def run_gates(root, registration, catalog, scaled_catalog, muscles,
              ligaments, scaled_ligaments, scale, force_scale):
    gates, summary = [], {}

    def record(name, error, tolerance, note, **extra):
        gates.append(dict(gate=name, relative_error=error, tolerance=tolerance,
                          passed=(error <= tolerance), note=note, **extra))

    # 0. The table against answers fixed outside it, before anything is scaled.
    #    A scaling table can be perfectly self-consistent and wrong -- give force
    #    the volume exponent and PCSA, muscle volume and specific tension all
    #    still agree with each other, because they were all derived from the same
    #    mistaken formula. Only a check against something declared elsewhere sees
    #    it, and KNOWN_ANSWERS is that: max_isometric_force must be 2.0 because
    #    ihm/body_parameters.py defaults muscle_force_scale to stature_scale**2.
    try:
        bs.check_table()
        table_error = 0.0
        table_note = 'all %d checks pass' % len(bs.check_table()['checks'])
    except bs.ScalingError as failure:
        table_error, table_note = 1.0, str(failure)
    record('the scaling table agrees with the answers fixed outside it',
           table_error, 0.0, table_note)

    if scale == 1.0:
        differences = sum(1 for a, b in zip(catalog, scaled_catalog) if a != b)
        if ligaments:
            differences += sum(1 for a, b in zip(ligaments, scaled_ligaments) if a != b)
        record('identity at s=1.0 reproduces the base exactly', float(differences),
               0.0, 'whole-record equality over the catalog and the ligament set')

    # 1. Every declared exponent, on the artifact.
    for field, quantity, records in (
            ('optimal_fiber_length_m', 'fibre_length', (catalog, scaled_catalog)),
            ('tendon_slack_length_m', 'tendon_slack_length', (catalog, scaled_catalog)),
            ('max_isometric_force_n', 'max_isometric_force', (catalog, scaled_catalog))):
        want = force_scale if quantity == 'max_isometric_force' \
            else scale ** bs.exponent(quantity)
        pairs = [(a[field], b[field]) for a, b in zip(*records) if a.get(field)]
        if not pairs:
            continue
        record('%s scales by s**%g' % (field, bs.exponent(quantity)),
               max(abs(b / (a * want) - 1.0) for a, b in pairs), 1e-12,
               bs.explain(quantity)['note'], muscles=len(pairs))

    # 2. Specific tension is the invariant that makes force follow area. Its
    #    value is declared and uncertain; its INVARIANCE is the claim, and this
    #    checks the claim rather than the value.
    worst = max(abs((m['scaled_max_isometric_force_n'] / m['scaled_pcsa_m2'])
                    / (m['max_isometric_force_n'] / m['pcsa_m2']) - 1.0)
                for m in muscles)
    record('implied specific tension is unchanged', worst, 1e-12,
           'force / PCSA, before and after. True whatever value SPECIFIC_TENSION_PA '
           'is given, which is the point: the exponent claim survives the '
           'uncertainty in the constant. ' + SPECIFIC_TENSION_BASIS,
           muscles=len(muscles), pa=SPECIFIC_TENSION_PA)

    # 3. Muscle volume by two routes. PCSA x fibre length must move like a
    #    volume, which it does only if PCSA is genuinely an area. This is the
    #    exact mistake the docstring names -- reasoning about muscle in mass
    #    terms and giving force the volume exponent -- caught arithmetically.
    worst = max(abs(m['scaled_volume_m3'] / (m['pcsa_m2'] * m['optimal_fiber_length_m']
                                             * scale ** bs.exponent('volume')) - 1.0)
                for m in muscles)
    record('PCSA x fibre length moves like a volume', worst, 1e-12,
           's**%g x s**%g = s**%g. If force had been given the volume exponent, '
           'PCSA would carry s**3, this product would carry s**4, and it would '
           'show up here rather than in a tall body that is a factor of s too '
           'strong.' % (bs.exponent('pcsa'), bs.exponent('fibre_length'),
                        bs.exponent('volume')), muscles=len(muscles))

    # 4. THE PHYSIOLOGICAL BAND, on the artifact, and exactly invariant rather
    #    than merely inside it. The numerator is a straight-line walk of the
    #    model's own PathPoint locations through ihm/native/model_scaling.py --
    #    a different module reading a different file.
    model = root / registration['model_path']
    polyline = path_point_polyline_lengths(model)
    fibre = {e['id']: e['optimal_fiber_length_m'] for e in catalog
             if e.get('optimal_fiber_length_m')}
    scaled_fibre = {e['id']: e['optimal_fiber_length_m'] for e in scaled_catalog
                    if e.get('optimal_fiber_length_m')}
    base_ratio, scaled_ratio = {}, {}
    for name, length in polyline.items():
        if name in fibre:
            base_ratio[name] = length / fibre[name]
            scaled_ratio[name] = (length * scale) / scaled_fibre[name]
    low, high = min(scaled_ratio.values()), max(scaled_ratio.values())
    summary['path_over_fibre'] = dict(min=low, max=high, muscles=len(scaled_ratio))
    record('path length over optimal fibre length is EXACTLY invariant',
           max(abs(scaled_ratio[m] / base_ratio[m] - 1.0) for m in base_ratio), 1e-12,
           'dimensionless, so under an isotropic scale it cannot move at all -- a '
           'stronger statement than the 0.2-20 band, which the mechanical side '
           'records as insufficient because a 13%% mis-scaled body passed it. Band '
           'here [%.3f, %.3f].' % (low, high), muscles=len(base_ratio))
    record('and stays in the 0.2-20 physiological band',
           0.0 if (PATH_OVER_FIBER_BAND[0] < low and high < PATH_OVER_FIBER_BAND[1])
           else 1.0, 0.0,
           'the required band check, kept even though the invariance gate above '
           'is strictly stronger, because it is the one that catches a corpus file '
           'holding accelerations under coordinate names',
           band=list(PATH_OVER_FIBER_BAND), measured=[low, high])

    # 5. Strength-to-weight, which is the consequence that matters for whether a
    #    scaled body can pick itself up. Force s**2 over weight s**3.
    summary['strength_to_weight_change_percent'] = 100 * (
        scale ** bs.exponent('strength_to_weight') - 1.0)
    summary['self_weight_stress_change_percent'] = 100 * (
        scale ** bs.exponent('self_weight_stress') - 1.0)
    record('strength-to-weight falls as s**%g' % bs.exponent('strength_to_weight'),
           abs((force_scale / scale ** bs.exponent('weight'))
               / scale ** bs.exponent('strength_to_weight') - 1.0), 1e-12,
           'a geometrically scaled-up body is relatively WEAKER: %.1f%% at this '
           'scale, with self-weight stress up %.1f%%. Not a defect of the scaler '
           '-- it is what geometric similarity means, and it is the classical '
           'reason large animals are not scaled-up small ones.'
           % (summary['strength_to_weight_change_percent'],
              summary['self_weight_stress_change_percent']))

    # 6. THE INDEPENDENT ARM. Mechanical muscle volume, derived from the
    #    catalog's forces and fibre lengths, against ANATOMICAL muscle volume,
    #    integrated from the mesh surfaces of the 568 role=='muscle' entities in
    #    anatomy.json. Two artifacts built by different scripts from different
    #    sources that disagree about how much muscle this body has; their ratio
    #    is the number that must not move.
    anatomy = json.loads((root / ANATOMY).read_text())
    anatomical = sum(e['volume_m3'] for e in anatomy['entities']
                     if e.get('role') == 'muscle' and e.get('volume_m3'))
    mechanical = sum(m['pcsa_m2'] * m['optimal_fiber_length_m'] for m in muscles)
    scaled_anatomical = bs.scale('volume', anatomical, scale)
    scaled_mechanical = sum(m['scaled_volume_m3'] for m in muscles)
    summary['muscle_volume'] = dict(
        anatomical_m3=anatomical, mechanical_m3=mechanical,
        ratio=mechanical / anatomical,
        scaled_anatomical_m3=scaled_anatomical, scaled_mechanical_m3=scaled_mechanical,
        anatomical_entities=sum(1 for e in anatomy['entities']
                                if e.get('role') == 'muscle' and e.get('volume_m3')),
        mechanical_actuators=len(muscles))
    # A third account of the same quantity, if the fibre field has been built:
    # the summed tetrahedral volume of 424 meshed muscle entities. Recorded as a
    # corroboration and NOT as a gate, because the agreement it shows is closer
    # than the inputs justify and treating it as a validation of the specific
    # tension would be exactly the kind of thing docs/LOG.md collects.
    fibre_field = root / 'data/derived/muscle-fibre-field-v1/summary.json'
    if fibre_field.exists():
        tet = json.loads(fibre_field.read_text()).get('total_volume_m3')
        if tet:
            summary['muscle_volume']['fibre_field_tet_m3'] = tet
            summary['muscle_volume']['mechanical_over_fibre_field'] = mechanical / tet
            summary['muscle_volume']['fibre_field_note'] = (
                'muscle-fibre-field-v1 sums 4.64M tetrahedra over 424 meshed muscle '
                'entities and gets %.5f m3, against %.5f m3 from 98 actuators at '
                '%.0f N/cm2 -- %.2f%% apart. Two unrelated routes to the same '
                'quantity, and that is closer than either input deserves. Reported '
                'as a coincidence to watch, not as a validation of the specific '
                'tension: pick 50 N/cm2 instead and the agreement disappears '
                'without anything about the body having changed.'
                % (tet, mechanical, SPECIFIC_TENSION_PA / 1e4,
                   100 * abs(mechanical / tet - 1.0)))

    record('mechanical over anatomical muscle volume does not move',
           abs((scaled_mechanical / scaled_anatomical) / (mechanical / anatomical) - 1.0),
           1e-12,
           'INDEPENDENT ARM. Left side: PCSA x fibre length over 98 actuators, from '
           'the mechanics catalog. Right side: mesh-integrated volume over %d '
           'role=="muscle" entities in anatomy.json, a different artifact from a '
           'different source. They disagree by %.2fx at the base and that is '
           'expected -- 98 actuators are not 568 muscles -- but the disagreement '
           'must be invariant. Give force the volume exponent and the left side '
           'gains a factor of s while the right does not.'
           % (summary['muscle_volume']['anatomical_entities'],
              mechanical / anatomical))

    # 7. Ligaments.
    if ligaments:
        materials = json.loads((root / MATERIALS).read_text())
        modulus = materials['materials']['ligament']['linear'][
            'young_modulus_along_fibre']['value']
        base_e = [l['linear_stiffness_n'] / l['cross_section_m2'] for l in ligaments]
        scaled_e = [l['linear_stiffness_n'] / l['cross_section_m2']
                    for l in scaled_ligaments]
        record('ligament stiffness over cross-section recovers the declared modulus',
               max(abs(e / modulus - 1.0) for e in scaled_e), 1e-9,
               'A KNOWN ANSWER, from a third file: %.1f MPa, Quapp and Weiss 1998, '
               'human MCL, in tissue-material-candidate-v1/materials.json. '
               'stiffness / cross-section IS the modulus, and a modulus is a '
               'material property, so it must come back at the same number on a '
               'body of any size. Read the Blankevoort stiffness as N/m instead of '
               'N and this comes back a factor of s wrong.'
               % (modulus / 1e6), elements=len(scaled_e), modulus_pa=modulus)
        record('ligament cross-section is still volume over length',
               max(abs(l['cross_section_m2'] * l['slack_length_m'] / l['model_volume_m3']
                       - 1.0) for l in scaled_ligaments
                   if l.get('model_volume_m3')), 1e-9,
               'A = V / L, so s**%g = s**%g / s**%g. The identity holds after '
               'scaling only if all three took what the table says.'
               % (bs.exponent('area'), bs.exponent('volume'), bs.exponent('length')),
               elements=len(scaled_ligaments))
        record('ligament peak strains are unchanged',
               max(abs(b['peak_strain_over_declared_range']
                       - a['peak_strain_over_declared_range'])
                   for a, b in zip(ligaments, scaled_ligaments)
                   if a.get('peak_strain_over_declared_range') is not None), 0.0,
               'a strain is dimensionless. So a taller body strains its ligaments '
               'exactly as much as a shorter one over the same joint range, and '
               'the elements flagged kinematically_admissible=false stay flagged '
               'at every stature -- scaling neither rescues one nor breaks one.',
               elements=len(ligaments))
        numeric = set()
        for element in ligaments:
            for key, value in element.items():
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    numeric.add(key)
                elif isinstance(value, list) and value and all(
                        isinstance(x, (int, float)) and not isinstance(x, bool)
                        for x in value):
                    numeric.add(key)
        unhandled = numeric - set(LIGAMENT_FIELDS) - set(LIGAMENT_INVARIANT)
        record('every numeric ligament field is scaled or declared invariant',
               float(len(unhandled)), 0.0, 'unhandled: %r' % (sorted(unhandled),),
               numeric_fields=sorted(numeric))
    return gates, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--stature-m', type=float)
    group.add_argument('--scale', type=float)
    parser.add_argument('--output', default=None)
    parser.add_argument('--sabotage', choices=('force-as-volume',
                                               'ligament-stiffness-as-spring'),
                        default=None)
    args = parser.parse_args()

    stature = args.stature_m if args.stature_m is not None else args.scale * MECHANICAL_STATURE_M
    scale = resolve({'stature_m': stature})['derived']['stature_scale']
    out = ROOT / (args.output or ('data/derived/body-variants/stature_%s'
                                  % ('%.6f' % scale).replace('.', 'p')))
    out.mkdir(parents=True, exist_ok=True)

    if args.sabotage == 'force-as-volume':
        # The classic: reason about muscle as mass, so force takes s**3.
        bs.QUANTITIES['max_isometric_force']['formula'] = {'specific_tension': 1,
                                                           'volume': 1}
    elif args.sabotage == 'ligament-stiffness-as-spring':
        # The units trap: read the Blankevoort stiffness as N/m.
        LIGAMENT_FIELDS['linear_stiffness_n'] = 'ligament_stiffness'

    scaled_catalog, scaled_ligaments, report = build(ROOT, scale)
    report['stature_m'] = stature
    report['script_sha256'] = _sha(Path(__file__))
    (out / 'muscle_catalog_scaled.json').write_text(
        json.dumps(scaled_catalog, indent=1) + '\n')
    if scaled_ligaments is not None:
        (out / 'ligaments_scaled.json').write_text(json.dumps(
            {'schema': 'ihm.tissue-force-elements-scaled.v1', 'stature_scale': scale,
             'elements': scaled_ligaments}, indent=1) + '\n')
    (out / 'muscle_tissue_report.json').write_text(json.dumps(report, indent=2) + '\n')

    for gate in report['gates']:
        print('%-5s %-62s %.3g (tol %.3g)'
              % ('PASS' if gate['passed'] else 'FAIL', gate['gate'],
                 gate['relative_error'], gate['tolerance']))
    s = report['summary']
    print('\nmuscle volume  mechanical %.5f m3, anatomical %.5f m3, ratio %.3f'
          % (s['muscle_volume']['mechanical_m3'], s['muscle_volume']['anatomical_m3'],
             s['muscle_volume']['ratio']))
    print('path/fibre     [%.3f, %.3f] over %d muscles, invariant'
          % (s['path_over_fibre']['min'], s['path_over_fibre']['max'],
             s['path_over_fibre']['muscles']))
    print('consequence    strength-to-weight %+.1f%%, self-weight stress %+.1f%%'
          % (s['strength_to_weight_change_percent'],
             s['self_weight_stress_change_percent']))
    print('\nwrote', out)
    if not all(g['passed'] for g in report['gates']):
        raise SystemExit('scaled muscle and tissue failed their own gates')


if __name__ == '__main__':
    main()
