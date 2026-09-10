"""What every quantity in the body does when the body changes size.

``ihm/native/model_scaling.py`` scales the 22-segment mechanical scaffold.  This
file scales **the body** -- the 4,000 anatomical entities, the 146 nerve routes
and their per-fibre-class conduction delays, the 1,326 skin patches, and the
muscle and tissue properties -- and it exists mostly to make one class of error
impossible.

The error is **an exponent chosen for convenience**.  Areas go as ``s**2`` and
volumes as ``s**3`` and everybody knows that; the ones that are actually got
wrong are the mixed quantities, where the temptation is to type a number that
looks about right.  Is a ligament's linear stiffness ``s``?  ``s**2``?  Invariant?
A conduction delay?  A moment arm?  A receptor density?  Typing any of those is
guessing, and a guess that happens to be right is indistinguishable from a guess
that happens to be wrong.

So **no exponent in this file is typed**.  Every one is *computed* from a
dimensional formula over a small set of primitives, and the only place a number
is written by hand is the primitive table -- where each entry is either 1 (a
length, because the scale is a linear map) or 0 (a material or chemical property,
because a taller person is made of the same stuff).  ``exponent('ligament_
stiffness')`` returns 1 because ``E*A/L`` expands to ``0 + 2 - 1``, not because
anyone decided it should.

The primitives are the whole argument, and they are stated as claims that can be
contradicted rather than as facts:

**geometric (exponent 1).**  ``length``.  One isotropic linear factor, the same
assumption ``stature_m`` already declares for the mechanical body and for exactly
the same reason: it is what makes every musculotendon path length scale by
``s`` at every pose.  It is a modelling assumption, and ``ALLOMETRY`` below is
the list of places where it is known to be false, each with the size of the error
it makes.

**material and chemical (exponent 0).**  Density, elastic modulus, muscle
specific tension, nerve conduction velocity, blood viscosity, regulated arterial
pressure.  These are properties of the tissue and not of the specimen.  The
load-bearing one is ``conduction_velocity``: axon diameter and internodal
myelin length are cellular dimensions and do not know how tall their owner is,
so a taller person's nerves are longer at the same velocity and **every
conduction delay in the body lengthens in proportion to stature**.  That is the
one physiological consequence of size that this repository can already measure
end to end, because ``docs/MILESTONES.md`` records that lumping conduction
delays costs as much as deleting an entire fibre group.

**counting (exponent 0).**  ``count``.  The number of discrete structures is a
property of the anatomy, not of its size: a taller person has the same 31 spinal
nerve pairs, the same 206 bones and -- see ``scripts/scale_skin_patches.py`` --
the same 1,326 innervated skin patches.  The consequence is that every *density*
per unit area or volume is NOT invariant, and comes out of the algebra with a
negative exponent rather than being overlooked.
"""
from __future__ import annotations

import math

SCHEMA = 'ihm.body-scaling.v1'


# ---------------------------------------------------------------------------
# The primitives.  These are the only hand-written exponents in this module.
# ---------------------------------------------------------------------------

PRIMITIVES: dict[str, dict] = {
    'length': dict(
        exponent=1, kind='geometric', unit='m',
        basis=('The isotropic linear factor itself. Under one uniform scale about '
               'a fixed origin every point maps to s times its position, so every '
               'distance between two points on the body maps to s times itself. '
               'This is a modelling assumption about how bodies differ in size, '
               'not a measurement; ALLOMETRY records where it is known to fail.')),
    'angle': dict(
        exponent=0, kind='geometric', unit='rad',
        basis=('A uniform scale is a conformal map and preserves angles exactly. '
               'This is mathematics, not physiology, and it is why the mechanical '
               'scaler leaves rotational coordinate ranges and '
               'pennation_angle_at_optimal alone.')),
    'fraction': dict(
        exponent=0, kind='geometric', unit='1',
        basis='Dimensionless by construction: a ratio of two like quantities.'),
    'count': dict(
        exponent=0, kind='counting', unit='1',
        basis=('The number of discrete anatomical structures is a property of the '
               'anatomy and not of its size. A taller person has the same number '
               'of bones, spinal roots, nerve trunks and innervated skin patches. '
               'Every density derived from a count therefore carries a NEGATIVE '
               'exponent, which is the point of declaring this rather than '
               'leaving it implicit.')),
    'density': dict(
        exponent=0, kind='material', unit='kg/m3',
        basis=('Tissue composition. The anatomical mass 70.7713 kg is composed by '
               'ihm/assembly/profile.py from a voxel partition at sourced '
               'per-constituent densities; those densities are properties of fat, '
               'muscle, bone and blood, not of the specimen carrying them.')),
    'elastic_modulus': dict(
        exponent=0, kind='material', unit='Pa',
        basis=("A material property of the tissue. This repository's own ligament "
               'along-fibre modulus -- 332.2 MPa, Quapp and Weiss 1998, human MCL, '
               'in data/derived/tissue-material-candidate-v1/materials.json -- is '
               'a property of collagen, and a taller person\'s collagen is the '
               'same collagen.')),
    'specific_tension': dict(
        exponent=0, kind='material', unit='Pa',
        basis=('Maximum muscle stress is set by the density of force-generating '
               'cross-bridges per unit area of myofilament lattice, which is a '
               'molecular arrangement. It is the reason maximum isometric force '
               'follows physiological cross-sectional area rather than muscle '
               'volume, and therefore the reason the mechanical scaler defaults '
               'muscle_force_scale to s**2.')),
    'conduction_velocity': dict(
        exponent=0, kind='material', unit='m/s',
        basis=('Set by axon diameter and internodal myelin length -- cellular '
               'dimensions, orders of magnitude below anything this atlas '
               'resolves. peripheral.json declares 17 fibre-class velocities as '
               'fixed constants and not one of them is a function of body size. '
               'THIS IS THE LOAD-BEARING INVARIANT: it is what makes a taller '
               "person's conduction delays longer rather than merely their "
               'nerves.')),
    'viscosity': dict(
        exponent=0, kind='material', unit='Pa s',
        basis='Blood composition; a property of the plasma and the cells in it.'),
    'regulated_pressure': dict(
        exponent=0, kind='regulated', unit='Pa',
        basis=('Mean arterial pressure is a regulated set-point, not a geometric '
               'consequence. Declared invariant because the baroreflex holds it '
               'so, which is a different kind of reason from the material '
               'invariants above and is labelled differently for that reason.')),
    'gravity': dict(
        exponent=0, kind='world', unit='m/s2',
        basis='A property of the world, not of the body.'),
    'authored_damping_time': dict(
        exponent=0, kind='authored', unit='s',
        basis=('The damping coefficient carried beside a ligament stiffness is an '
               'authored fraction with units of seconds, not a measured tissue '
               'property -- the same status ihm/native/model_scaling.py gives the '
               'Hunt-Crossley contact parameters, and it is held for the same '
               'reason: nothing in this repository fixes its scaling. Under '
               'gravitational dynamic similarity a characteristic time would go '
               'as s**0.5 (equal Froude number, which is why a taller person '
               'walks at a lower cadence); ALLOMETRY records that, and this '
               'scaler does not apply it, because rescaling time would change '
               'the meaning of every declared time constant in the body at '
               'once.')),
}


# ---------------------------------------------------------------------------
# The derived quantities.  Each is a FORMULA, never an exponent.
# ---------------------------------------------------------------------------
#: name -> (formula over primitives and other quantities, note)
QUANTITIES: dict[str, dict] = {
    # -- pure geometry ------------------------------------------------------
    'length':       dict(formula={'length': 1},  note='any distance, extent, radius, route or slack length'),
    'area':         dict(formula={'length': 2},  note='any surface area'),
    'volume':       dict(formula={'length': 3},  note='any enclosed volume'),
    'centroid':     dict(formula={'length': 1},  note='a position; scales about the declared scaling origin'),
    'moment_arm':   dict(formula={'length': 1},  note='OpenSim takes it as -dL/dq for rotational q, so it is a length'),
    'pennation_angle': dict(formula={'angle': 1},
                            note='pennation_angle_at_optimal. Invariant, and '
                                 'ihm/native/model_scaling.py already leaves it alone '
                                 'for exactly this reason'),
    'joint_range':  dict(formula={'angle': 1},
                         note='a rotational coordinate range; invariant, and the same '
                              'entry on the mechanical side'),
    'surface_normal': dict(formula={'angle': 1},
                           note='a unit direction. Invariant under a uniform scale, '
                                'which is why the 1,326 patch normals are copied '
                                'rather than transformed'),

    # -- mass and inertia ---------------------------------------------------
    'mass':         dict(formula={'density': 1, 'volume': 1},
                         note='constant composition, so mass follows volume'),
    'weight':       dict(formula={'mass': 1, 'gravity': 1}, note='what the floor has to hold'),
    'inertia':      dict(formula={'mass': 1, 'length': 2},
                         note='second moment of mass; cross-checked against the '
                              'mechanical scaler, which uses mass_scale * s**2'),

    # -- muscle -------------------------------------------------------------
    'fibre_length':        dict(formula={'length': 1}, note='optimal fibre length'),
    'tendon_slack_length': dict(formula={'length': 1}, note='series-elastic slack length'),
    'pcsa':                dict(formula={'area': 1},
                                note='physiological cross-sectional area, measured '
                                     'perpendicular to the fibres'),
    'max_isometric_force': dict(formula={'specific_tension': 1, 'pcsa': 1},
                                note='force = stress x area; this is WHY force goes '
                                     'as s**2 and not as s**3'),
    'muscle_volume_via_fibres': dict(formula={'pcsa': 1, 'fibre_length': 1},
                                     note='must equal volume; the table checks it'),
    'muscle_mass':         dict(formula={'density': 1, 'muscle_volume_via_fibres': 1},
                                note='independent route to the same s**3'),
    'joint_moment':        dict(formula={'max_isometric_force': 1, 'moment_arm': 1},
                                note='what a muscle can do at a joint'),
    'contraction_velocity_absolute': dict(formula={'fibre_length': 1},
                                          note='max_contraction_velocity is declared in '
                                               'OPTIMAL FIBRE LENGTHS per second and is '
                                               'therefore invariant; the absolute m/s '
                                               'velocity it implies is not'),
    'normalised_fibre_length': dict(formula={'length': 1, 'fibre_length': -1},
                                    note='path length over optimal fibre length: the '
                                         'gate band 0.2-20 is on THIS, and it must be '
                                         'exactly invariant, not merely in band'),

    # -- passive tissue -----------------------------------------------------
    'ligament_stiffness': dict(formula={'elastic_modulus': 1, 'area': 1, 'length': -1},
                               note='k = E A / L for a linear axial element, in N/m; '
                                    'the cross-section grows faster than the length, so '
                                    'a bigger ligament is stiffer in N per metre'),
    'ligament_stiffness_per_strain': dict(formula={'elastic_modulus': 1, 'area': 1},
                                          note='E*A, in NEWTONS -- what a '
                                               'Blankevoort1991Ligament calls its '
                                               'stiffness, because its extension '
                                               'variable is a dimensionless strain '
                                               'rather than a length. It is therefore '
                                               'a force and must scale like one, and '
                                               'the table checks that it does. Reading '
                                               'it as the N/m entry above would be off '
                                               'by a whole factor of s.'),
    'ligament_damping': dict(formula={'ligament_stiffness_per_strain': 1,
                                      'authored_damping_time': 1},
                             note='N s per unit strain. The damping is an authored '
                                  'fraction of the stiffness with units of seconds, '
                                  'not a measured tissue property'),
    'ligament_force_at_strain': dict(formula={'ligament_stiffness': 1, 'length': 1},
                                     note='k x (strain x L): strain is dimensionless, '
                                          'so force goes as s**2 like every other force'),
    'strain':             dict(formula={'fraction': 1}, note='dimensionless by construction'),
    'stress':             dict(formula={'max_isometric_force': 1, 'area': -1},
                               note='invariant, as it must be: it is the specific '
                                    'tension it came from'),
    'self_weight_stress': dict(formula={'weight': 1, 'area': -1},
                               note="stress in a load-bearing cross-section under the "
                                    "body's OWN weight. NOT invariant: this is the "
                                    'classical reason large animals are not scaled-up '
                                    'small ones, and it comes out of the algebra rather '
                                    'than being remembered'),
    'strength_to_weight': dict(formula={'max_isometric_force': 1, 'weight': -1},
                               note='a geometrically scaled-up body is relatively '
                                    'WEAKER; directly relevant to whether it can pick '
                                    'itself up'),

    # -- nerve --------------------------------------------------------------
    'route_length':     dict(formula={'length': 1},
                             note='receptor or muscle endpoint to relay, over the body'),
    'conduction_delay': dict(formula={'route_length': 1, 'conduction_velocity': -1},
                             note='THE GATE. Longer nerves at the same velocity means '
                                  'a longer delay, in exact proportion to stature'),
    'central_delay':    dict(formula={}, note='synaptic and central transport; not a '
                                              'peripheral route length, and this scaler '
                                              'does not claim to scale it'),
    'reflex_latency':   dict(formula={'conduction_delay': 1},
                             note='afferent + efferent route, both scaling alike; the '
                                  'central component above does not'),

    # -- skin and receptors -------------------------------------------------
    'patch_area':       dict(formula={'area': 1}, note='one dermatome patch'),
    'patch_count':      dict(formula={'count': 1},
                             note='HELD FIXED BY CHOICE; see scripts/scale_skin_patches.py '
                                  'for the choice and its justification'),
    'receptor_density': dict(formula={'count': 1, 'area': -1},
                             note='the consequence of the two above: a taller body has '
                                  'the same number of afferent channels spread over more '
                                  'skin, so density falls as s**-2'),

    # -- vasculature --------------------------------------------------------
    'vessel_radius':       dict(formula={'length': 1}, note='calibre under isotropy'),
    'vascular_resistance': dict(formula={'viscosity': 1, 'length': 1, 'vessel_radius': -4},
                                note='Poiseuille: 8 mu L / (pi r**4). The r**4 is why a '
                                     'naive calibre scale is a physiological change and '
                                     'not a geometric one'),
    'volumetric_flow':     dict(formula={'regulated_pressure': 1, 'vascular_resistance': -1},
                                note='at a regulated driving pressure. Comes out at s**3, '
                                     'which ALLOMETRY records as the wrong answer'),

    # -- whole-body ---------------------------------------------------------
    'bmi': dict(formula={'mass': 1, 'length': -2},
                note='NOT invariant under geometric similarity, and in the population '
                     'it very nearly is; ALLOMETRY carries the measurement'),
    'body_surface_area': dict(formula={'area': 1}, note='exterior skin'),
    'surface_to_volume':  dict(formula={'area': 1, 'volume': -1},
                               note='falls as 1/s; thermal and metabolic consequences '
                                    'are not modelled here and this is recorded so that '
                                    'is visible'),
}


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------

class ScalingError(ValueError):
    """A quantity the scaling table cannot resolve, or a table that is wrong."""


def _resolve(name, stack=()):
    """(exponent, expansion over primitives) for one quantity."""
    if name in stack:
        raise ScalingError('Cyclic scaling formula: ' + ' -> '.join(stack + (name,)))
    if name in PRIMITIVES:
        return float(PRIMITIVES[name]['exponent']), {name: 1.0}
    entry = QUANTITIES.get(name)
    if entry is None:
        raise ScalingError('Unknown scaling quantity %r' % (name,))
    total, expansion = 0.0, {}
    for term, power in entry['formula'].items():
        sub_exponent, sub_expansion = _resolve(term, stack + (name,))
        total += sub_exponent * power
        for primitive, sub_power in sub_expansion.items():
            expansion[primitive] = expansion.get(primitive, 0.0) + sub_power * power
    return total, {k: v for k, v in expansion.items() if v}


def exponent(name):
    """The power of the linear scale ``s`` that ``name`` takes.

    Computed from the formula table every time.  There is deliberately no cache
    and no literal: the only way to change an exponent is to change a formula or
    a primitive, both of which carry a written justification.
    """
    return _resolve(name)[0]


def explain(name):
    """Everything a reader needs to disagree with an exponent."""
    value, expansion = _resolve(name)
    entry = QUANTITIES.get(name) or PRIMITIVES[name]
    return dict(
        quantity=name, exponent=value,
        formula={k: v for k, v in entry.get('formula', {name: 1}).items()},
        primitive_expansion=expansion,
        note=entry.get('note') or entry.get('basis'),
        primitive_basis={p: PRIMITIVES[p]['basis'] for p in expansion})


def scale(name, value, s):
    """``value`` at linear scale ``s``.  ``None`` passes through unchanged."""
    if value is None:
        return None
    return float(value) * float(s) ** exponent(name)


def scale_vector(name, values, s):
    factor = float(s) ** exponent(name)
    return [float(v) * factor for v in values]


# ---------------------------------------------------------------------------
# What geometric similarity gets wrong, with the size of the error
# ---------------------------------------------------------------------------
#: Each entry is a place where the isotropic assumption is known to be false for
#: a real human.  None of them is silently corrected.  A correction needs an
#: exponent, an exponent needs evidence, and where this repository holds the
#: evidence the entry says so and names the script that measured it; where it
#: does not, the entry says that too and the isotropic value is what the
#: materializer writes, flagged rather than fixed.
ALLOMETRY = (
    dict(quantity='mass',
         isotropic_exponent=3.0,
         measured_exponent=2.034,
         evidence='measured_in_repo',
         measured_by='scripts/measure_stature_allometry.py',
         source='NHANES 2017-2018 BMX_J + DEMO_J, 4,822 adults aged 20-79, '
                'survey-weighted by WTMEC2YR',
         detail=('Weighted regression of log weight on log standing height gives '
                 '2.034 +/- 0.083 pooled, 2.372 for men and 1.721 for women -- not '
                 '3. Isometry is rejected at z = -11.6. The '
                 'pooled slope lies BETWEEN the two within-sex slopes, so it is '
                 'not an artefact of mixing two groups with different means. This '
                 'is a between-person association and is confounded by adiposity, '
                 'which is exactly the right comparison for the question "if I ask '
                 'for a 2.03 m person, what should they weigh".'),
         consequence=('At stature 2.03 m the geometric-similarity mass is 111.9 kg '
                      'and the population-exponent mass is 99.2 kg, 12.8% apart.'),
         disposition=('APPLIED as the DEFAULT for mass_kg in '
                      'ihm/body_parameters.py, on the WITHIN-SEX slope (2.372 '
                      'male, 1.721 female) rather than the pooled 2.034, since '
                      'the body being resolved has a declared sex. Previously '
                      'not applied at all, which meant a stature-only request '
                      'kept the 77.6 kg literal and the stature scaling was '
                      'divided back out of the mass. mass_kg REMAINS an '
                      'independent knob -- an explicit request still wins '
                      'outright and is not checked against this curve -- so a '
                      'caller who wants s**3, or any other pair, still gets it '
                      'by naming the mass. This scaler itself still does not '
                      'move mass; the parameter layer does.')),
    dict(quantity='bmi',
         isotropic_exponent=1.0,
         measured_exponent=0.034,
         evidence='measured_in_repo',
         measured_by='scripts/measure_stature_allometry.py',
         source='the same NHANES sample',
         detail=('Geometric similarity says a taller person has a proportionally '
                 'higher BMI. The population says BMI is essentially independent '
                 'of stature -- exponent 0.034 pooled, z = -11.7 against the '
                 'predicted 1.0 -- which agrees with the '
                 "repository's own already-measured male/female BMI difference of "
                 'd = -0.03 in data/derived/anthropometry/.'),
         consequence='At 2.03 m, isotropy predicts BMI 24.0 -> 27.1.',
         disposition='NOT applied; reported as the falsification of the assumption.'),
    dict(quantity='volumetric_flow',
         isotropic_exponent=3.0,
         measured_exponent=2.25,
         evidence='published_prior_not_measured_here',
         source='Kleiber: metabolic rate ~ M**0.75, so cardiac output ~ M**0.75 = '
                's**2.25 under mass ~ s**3',
         detail=('Isotropic calibre scaling makes vascular resistance fall as '
                 's**-3 (the r**4 in Poiseuille beats the length), so flow at a '
                 'regulated pressure rises as s**3. Metabolic demand rises more '
                 'slowly. The mismatch is s**0.75, which is 9.6% of over-perfusion '
                 'at 2.03 m.'),
         consequence=('A 2.03 m isotropically scaled body is perfused about 9.6% '
                      'above what its metabolic scaling asks for.'),
         disposition=('NOT applied. No vessel radii are rewritten by this scaler at '
                      'all -- vascular entities take the same placement and extent '
                      'scaling as every other entity and no flow model is '
                      'reparameterised. Correcting it would need a measured '
                      'stature-to-calibre relation, which no catalogued source '
                      'here provides.')),
    dict(quantity='brain_volume',
         isotropic_exponent=3.0,
         measured_exponent=None,
         evidence='no_source_in_repository',
         source=None,
         detail=('Adult human brain volume is only weakly related to stature; the '
                 'intraspecific brain-body allometric exponent is far below the '
                 'interspecific 0.75 and nothing close to isometry. No catalogued '
                 'source in this repository measures it, so no exponent is '
                 'asserted.'),
         consequence=('The 146 nervous-system entities are scaled isotropically '
                      'like everything else, which at 2.03 m makes the brain 1.44x '
                      'its volume. That is very probably too big.'),
         disposition=('NOT corrected, and flagged in every variant. This is a known '
                      'wrong number, not an unexamined one.')),
    dict(quantity='receptor_density',
         isotropic_exponent=-2.0,
         measured_exponent=-2.0,
         evidence='chosen_and_consistent',
         source='patch count held fixed; see scripts/scale_skin_patches.py',
         detail=('Not an error -- an entry here because it is the one place where '
                 'a NEGATIVE exponent is the intended answer and would otherwise '
                 'look like a bug. Holding the 1,326 patches fixed while their '
                 'area grows is a choice, and it is the choice that keeps the '
                 "brain's afferent channel count invariant under body size."),
         consequence='At 2.03 m, afferent channels per m2 of skin fall by 21.6%.',
         disposition='Applied, deliberately.'),
    dict(quantity='authored_damping_time',
         isotropic_exponent=0.0,
         measured_exponent=0.5,
         evidence='published_prior_not_measured_here',
         source='gravitational dynamic similarity: equal Froude number '
                'v**2/(g L) makes a characteristic time go as sqrt(L/g)',
         detail=('This scaler rescales lengths and holds time. Two bodies of '
                 'different size moving dynamically alike under the same gravity '
                 'do NOT share a clock: the natural period goes as s**0.5, which '
                 'is why a taller person walks at a lower cadence at the same '
                 'dimensionless speed. 6.3% at 2.03 m.'),
         consequence=('Every declared time constant -- damping coefficients, '
                      'receptor and activation tau, the integrator step -- is the '
                      "source subject's. A scaled body is dynamically similar in "
                      'geometry and not in time.'),
         disposition=('NOT applied, deliberately. Rescaling time would change the '
                      'meaning of every time constant in the body at once, '
                      'including the conduction delays, which are already correct '
                      'for a different reason -- they come from a length and a '
                      'velocity, not from a clock.')),
    dict(quantity='passive_joint_stop_moment',
         isotropic_exponent=3.0,
         measured_exponent=None,
         evidence='no_source_in_repository',
         source=None,
         detail=('Inherited unchanged from ihm/native/model_scaling.py, which '
                 'records the same gap: the ExpressionBasedCoordinateForceSet '
                 'shoulder/elbow/hip stops are authored range-of-motion limits '
                 'rather than measured tissue properties, so nothing fixes their '
                 'scaling. Repeated here so a reader of the body-side scaling does '
                 'not have to find it on the mechanical side.'),
         consequence="A scaled body has the source subject's passive joint stops.",
         disposition='NOT applied; a known gap on both sides.'),
)


# ---------------------------------------------------------------------------
# The table has to check itself
# ---------------------------------------------------------------------------

#: Exponents whose value is fixed by something outside this table, so that the
#: table can be checked against an answer that was known before it was written.
#: Two of them are the mechanical scaler's own literals, which is the point:
#: this file and ``ihm/native/model_scaling.py`` were written independently and
#: must agree, or one of them is wrong.
KNOWN_ANSWERS = {
    'length': 1.0,
    'area': 2.0,
    'volume': 3.0,
    'angle': 0.0,
    'strain': 0.0,
    # ihm/native/model_scaling.py: body mass takes mass_scale, default s**3.
    'mass': 3.0,
    # ihm/native/model_scaling.py: inertia takes mass_scale * scale**2 = s**5.
    'inertia': 5.0,
    # ihm/body_parameters.py: muscle_force_scale defaults to stature_scale**2.
    'max_isometric_force': 2.0,
    # A stress is the specific tension it came from, so it cannot move.
    'stress': 0.0,
    # scripts/materialize_stature_variant.py gates path/optimal fibre length as
    # exactly invariant, not merely in band.
    'normalised_fibre_length': 0.0,
    # The headline consequence, and the reason this module exists.
    'conduction_delay': 1.0,
}


def check_table():
    """Validate the table against things that were true before it was written.

    Returns a report.  Raises on anything that does not hold, because a scaling
    table that disagrees with itself will produce plausible numbers forever.
    """
    report = {'schema': SCHEMA, 'primitives': len(PRIMITIVES),
              'quantities': len(QUANTITIES), 'checks': []}

    def check(name, got, want, note):
        ok = (got is None and want is None) or abs(got - want) < 1e-12
        report['checks'].append(dict(check=name, got=got, expected=want,
                                     passed=ok, note=note))
        if not ok:
            raise ScalingError('%s: got %r, expected %r (%s)' % (name, got, want, note))

    # 1. Every quantity resolves, and no formula is cyclic.
    for name in QUANTITIES:
        exponent(name)
    check('every quantity resolves', float(len(QUANTITIES)), float(len(QUANTITIES)),
          'a formula naming an unknown term or forming a cycle raises here')

    # 2. Exponents whose value was fixed elsewhere.
    for name, want in KNOWN_ANSWERS.items():
        check('known answer: ' + name, exponent(name), want,
              'fixed outside this table; see KNOWN_ANSWERS')

    # 3. The table must agree with ITSELF where two formulas describe the same
    #    physical quantity by different routes.  This is the check that catches
    #    a plausible-looking formula: muscle volume as PCSA x fibre length has
    #    to come out at the same exponent as volume as length cubed, and it only
    #    does if PCSA really is an area and fibre length really is a length.
    for a, b, why in (
            ('muscle_volume_via_fibres', 'volume',
             'V = PCSA x optimal fibre length, and separately V = L**3. If PCSA '
             'had been given the volume exponent by mistake these would differ by '
             'one factor of s and the discrepancy would be visible here rather '
             'than in a muscle that is silently the wrong mass.'),
            ('muscle_mass', 'mass', 'the same identity carried through density'),
            ('ligament_force_at_strain', 'max_isometric_force',
             'every force in the body scales alike, however it is generated; a '
             'passive tissue force and an active muscle force must not diverge'),
            ('weight', 'mass', 'g is a property of the world'),
            ('reflex_latency', 'conduction_delay',
             'a reflex latency is two route lengths at fixed velocities'),
            ('ligament_stiffness_per_strain', 'max_isometric_force',
             'a Blankevoort ligament stiffness is in NEWTONS, because its '
             'extension variable is a dimensionless strain. It is a force and '
             'must scale like every other force in the body; reading it as the '
             'N/m spring constant instead would be wrong by a whole factor of s'),
            ('ligament_damping', 'max_isometric_force',
             'the damping coefficient is authored in seconds and held, so the '
             'damping follows its stiffness exactly')):
        check('self-consistency: %s == %s' % (a, b), exponent(a), exponent(b), why)

    # 4. Every primitive is used by at least one quantity, and every primitive
    #    exponent is 0 or 1.  A primitive with an exponent of its own that is
    #    neither would be a derived quantity wearing a primitive's coat.
    used = set()
    for name in QUANTITIES:
        used |= set(_resolve(name)[1])
    unused = sorted(set(PRIMITIVES) - used)
    check('every primitive is used', float(len(unused)), 0.0,
          'unused: %r -- a primitive nothing consumes is dead justification' % (unused,))
    bad = [p for p, e in PRIMITIVES.items() if e['exponent'] not in (0, 1)]
    check('primitives are 0 or 1', float(len(bad)), 0.0,
          'a primitive with any other exponent is a derived quantity in disguise: %r' % (bad,))

    # 5. Allometry entries must name a quantity the table knows, so that the
    #    isotropic exponent they claim to contradict is the one actually used.
    for entry in ALLOMETRY:
        name = entry['quantity']
        if name in PRIMITIVES:
            check('allometry baseline: ' + name, float(PRIMITIVES[name]['exponent']),
                  float(entry['isotropic_exponent']),
                  'an allometry entry on a PRIMITIVE contradicts the primitive '
                  'exponent itself')
            continue
        if name not in QUANTITIES and name not in ('brain_volume', 'passive_joint_stop_moment'):
            raise ScalingError('ALLOMETRY names unknown quantity %r' % (name,))
        if name in QUANTITIES:
            check('allometry baseline: ' + name, exponent(name),
                  float(entry['isotropic_exponent']),
                  'the isotropic exponent an allometry entry contradicts must be '
                  'the one this table actually computes')
    report['allometry_entries'] = len(ALLOMETRY)
    return report


def table():
    """The whole table, resolved, for writing into a variant's provenance."""
    return dict(
        schema=SCHEMA,
        primitives={name: dict(entry) for name, entry in PRIMITIVES.items()},
        quantities={name: explain(name) for name in sorted(QUANTITIES)},
        allometry=[dict(entry) for entry in ALLOMETRY],
        known_answers=dict(KNOWN_ANSWERS))


# ---------------------------------------------------------------------------
# The scaling origin
# ---------------------------------------------------------------------------

#: The anatomical body is scaled by one uniform factor about the ORIGIN of the
#: canonical frame.  The choice of origin is recorded rather than hidden, and
#: the reason it does not matter is worth stating: every quantity this module
#: scales is either translation-invariant (an extent, an area, a volume, a
#: route length, a distance between two entities) or is a position that the
#: playback path re-derives from the segment pose anyway.  A gate checks the
#: translation-invariant half directly, because that check passes or fails
#: independently of where the origin was put.
SCALING_ORIGIN = 'canonical frame origin (0, 0, 0); see docs/BODY_PARAMETERS.md'


def scale_point(point, s, origin=(0.0, 0.0, 0.0)):
    return [float(o) + (float(p) - float(o)) * float(s) for p, o in zip(point, origin)]


def polyline_length(points):
    return float(sum(math.dist(a, b) for a, b in zip(points, points[1:])))
