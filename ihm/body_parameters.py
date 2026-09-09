"""One declared schema for the parameters that describe *which body* is being run.

Before this file the answer to "how tall is the body, how heavy, and what sex"
was spread over three places that disagree with each other, and one of the three
was a literal repeated in fifty-eight scripts.  Nothing here invents a new knob:
every entry names the call site that already consumes it, or says plainly that
nothing does.

Three statuses, and the difference between them is the whole point of the file:

``surfaced``
    Already a runtime knob before this schema existed.  The schema only gives it
    a name, a unit, a range and a provenance.  ``mass_kg`` is the case.
``implemented``
    Made into a knob by this work.  ``stature_m`` is the case; it is realised by
    ``ihm.native.model_scaling`` and materialised by
    ``scripts/materialize_stature_variant.py``.
``declared``
    Written down, with a domain, and **not** reaching the running body.  ``sex``
    is the case, and its domain has exactly one member because the simulated
    body has zero sex-specific anatomical entities.  A declared parameter is a
    statement of scope, never a capability.

The two bodies.  This programme runs an anatomical body (BodyParts3D + BioGears,
described by ``ihm/assembly/profile.py``) and a mechanical body (the OpenSim
Rajagopal subject driven by ``ihm.native.mechanical_stream``).  They are
different objects with different masses and different statures, and the schema
holds both rather than picking one and hiding the other.  ``docs/DISCONNECTS.md``
item 1 is the same seam seen from the geometry side.
"""
from pathlib import Path
import hashlib, json, math

# ---------------------------------------------------------------------------
# Measured constants.  Every one is recomputed by ``assert_measured_defaults``
# from the artifact named beside it, so none of them can drift away from its
# own evidence.  That is the pattern ``ihm/assembly/profile.py`` uses for
# ``MASS_KG`` and the reason its 70.7713 kg is trustworthy.
# ---------------------------------------------------------------------------

MECHANICAL_MODEL = ('data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/'
                    'example3DWalking/subject_walk_scaled.osim')
MECHANICAL_PATHSET = ('data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/'
                      'example3DWalking/subject_walk_scaled_FunctionBasedPathSet.xml')
ANATOMICAL_PROFILE = 'data/derived/canonical/profile.json'

#: Sum of the 22 ``<Body><mass>`` entries in the source ``.osim``.  This is what
#: the native engine divides into ``target_mass_kg`` to get its uniform
#: ``mass_scale`` (``scripts/native_mechanical_stream.cpp``, line 51).
MECHANICAL_SOURCE_MASS_KG = 85.26984854173146

#: The mass every mechanical script actually asks for.  It is a literal, and no
#: derivation for it exists anywhere in this repository: it is not the source
#: model's mass (85.270 kg), not the anatomical body's composed mass
#: (70.7713 kg), not the superseded BioGears StandardMale constant
#: (77.1107029 kg), and not recoverable from the trial's ground reaction force
#: (mean vertical GRF over ``grf_walk.mot`` gives 60.83 kg, which is a duty-cycle
#: artefact, not a weight).  It is most likely the AddBiomechanics subject's
#: measured mass, carried in by hand.  Recorded as unattributed rather than
#: guessed.
MECHANICAL_TARGET_MASS_KG = 77.6122029

#: Vertical distance, at the model's default pose, from the plane of the
#: AddBiomechanics ``*Ground`` virtual markers to the ``Head`` marker.  Measured
#: by ``ihm.native.model_scaling.head_marker_height_m``; see that function for
#: why this and not a bounding box.
MECHANICAL_STATURE_M = 1.7972725296074763

#: Head-to-toe extent of the canonical skin mesh (BodyParts3D ``FJ2810``), as
#: written into ``data/derived/canonical/profile.json`` by
#: ``ihm/assembly/profile.py``.
ANATOMICAL_STATURE_M = 1.7194712

#: Composed over the anatomical body's own measured interior; supersedes the
#: inherited BioGears StandardMale 77.1107029 kg.  See ``profile.py``.
ANATOMICAL_MASS_KG = 70.7713

#: Fitted similarity scale taking OpenSim body geometry onto the BodyParts3D
#: atlas (``scripts/bind_anatomy_to_segments.py``; ``docs/ANATOMY_SEGMENT_BINDING.md``).
ANATOMY_REGISTRATION_SCALE = 0.963

#: Where the sex-stratified proportions come from.  Written by
#: ``scripts/index_anthropometry.py``; survey-weighted over 4,883 NHANES
#: 2017-2018 adults aged 20-79, with the weighting checked against the published
#: CDC/NCHS means (175.4 against 175.3 cm for men, 161.3 against 161.3 for
#: women).  ``assert_measured_defaults`` recomputes every ratio below from it.
ANTHROPOMETRY_SUMMARY = 'data/derived/anthropometry/nhanes-2017-2018/summary.json'

#: Female/male ratios.  Each is a ratio of survey-weighted means, so it says how
#: this SUBJECT would be re-proportioned, not that a population mean replaces
#: the subject.  ``d`` beside each is the weighted standardised difference: a
#: large absolute difference can be a small difference in shape and the ratios
#: are where that shows.
FEMALE_OVER_MALE = {
    'stature':          0.9199420560444753,   # d = 1.98, the dominant term
    'mass':             0.8515774154511445,   # d = 0.62
    'hip_over_stature': 1.1264611737929469,   # d = -0.91, what stature cannot carry
    'upper_leg_over_stature': 0.97447397046321,   # d = 0.44
    'upper_arm_over_stature': 0.9937527178769803,  # d = 0.14, under 1%
    'waist_over_hip':   0.9194355202181121,   # d = 1.07, soft tissue only
}

#: Bodies scaled by ``leg_length_ratio``.  Isotropic per segment, and that is a
#: decision with a reason: ``walker_knee_r`` takes its three translation
#: functions in a frame rotated by (-1.64, 1.45, 1.57) rad inside ``femur_r``,
#: so a LONG-AXIS-ONLY femur scale has no componentwise meaning there and
#: ``ihm.native.anisotropic_scaling`` refuses it.  A uniform per-segment factor
#: is exact regardless of frame orientation.  It shortens the leg and slims it
#: together, which is a modelling choice, not a measurement.
LEG_SEGMENT_BODIES = ('femur_r', 'femur_l', 'tibia_r', 'tibia_l', 'patella_r',
                      'patella_l', 'talus_r', 'talus_l', 'calcn_r', 'calcn_l',
                      'toes_r', 'toes_l')


def _rel(x, y):
    return abs(x - y) / y


#: How far apart the two bodies are, as a fraction.  Printed rather than
#: absorbed, because it is the honest error bar on any statement of the form
#: "the body is N metres tall".
STATURE_DISAGREEMENT = _rel(MECHANICAL_STATURE_M * ANATOMY_REGISTRATION_SCALE,
                            ANATOMICAL_STATURE_M)
MASS_DISAGREEMENT = _rel(MECHANICAL_TARGET_MASS_KG, ANATOMICAL_MASS_KG)


# ---------------------------------------------------------------------------
# The schema
# ---------------------------------------------------------------------------

PARAMETERS = (
    dict(
        name='stature_m', unit='m', kind='continuous',
        range=(1.40, 2.05), default=MECHANICAL_STATURE_M, status='implemented',
        body='mechanical',
        basis=('Operationally, the vertical distance at the model default pose '
               'from the AddBiomechanics floor-level virtual markers '
               '(R/L.HeelGround, MT5Ground, ToeGround) to the Head marker. This '
               'is a stature PROXY, not a measured standing height: the Head '
               'marker sits on the head, not at the vertex. Its size is bounded '
               'below -- carried onto the atlas by the fitted registration scale '
               '%.3f it gives %.4f m against the atlas skin extent %.4f m, a '
               '%.2f%% disagreement.'
               % (ANATOMY_REGISTRATION_SCALE,
                  MECHANICAL_STATURE_M * ANATOMY_REGISTRATION_SCALE,
                  ANATOMICAL_STATURE_M, 100 * STATURE_DISAGREEMENT)),
        range_basis=('Measured, from NHANES 2017-2018 body measures over 4,883 '
                     'adults aged 20-79 (scripts/index_anthropometry.py). The '
                     'survey-weighted 0.5th-99.5th percentile of standing height '
                     'is 1.437-1.929 m over both sexes, 1.560-1.955 m for men '
                     'and 1.437-1.801 m for women. The declared bound is widened '
                     'to 1.40-2.05 m and the widening is a modelling choice, not '
                     'a measurement. The bound that actually bites is the '
                     'scaling: beyond roughly 0.8x-1.15x the source subject the '
                     'fitted path polynomials are being evaluated far outside the '
                     'coordinate box they were fitted on.'),
        source=MECHANICAL_MODEL,
        consumers=('ihm.native.model_scaling.scale_model',
                   'scripts/materialize_stature_variant.py',
                   'NativeMechanicalStream(augmented_registration=...)'),
        derived=('stature_scale',),
    ),
    dict(
        name='mass_kg', unit='kg', kind='continuous',
        range=(35.0, 160.0), default=MECHANICAL_TARGET_MASS_KG, status='surfaced',
        body='mechanical',
        basis=('Unattributed literal, replicated across 67 call sites in 58 files. See '
               'MECHANICAL_TARGET_MASS_KG for the four candidate derivations '
               'that were checked and rejected.'),
        range_basis=('Measured. The same NHANES sample gives a survey-weighted '
                     '0.5th-99.5th percentile weight of 45.1-164.5 kg over both '
                     'sexes. The declared bound 35-160 kg is close to that and '
                     'is not identical to it; the native engine itself accepts '
                     'any positive mass, so what the bound buys is that a typo '
                     'fails here rather than integrating.'),
        source=MECHANICAL_MODEL,
        consumers=('NativeMechanicalStream(target_mass_kg=...)',
                   'scripts/native_mechanical_stream.cpp line 52'),
        derived=('mass_scale',),
        limitation=('The native engine multiplies every body mass AND every '
                    'inertia tensor by the same scalar. Mass distribution is '
                    'therefore invariant under this knob: a 120 kg body has the '
                    'segment mass fractions of a 77.6 kg one, and its radii of '
                    'gyration are unchanged. Adiposity is not represented.'),
    ),
    dict(
        name='muscle_force_scale', unit='dimensionless', kind='continuous',
        range=(0.25, 4.0), default=None, status='implemented',
        body='mechanical',
        basis=('Defaults to stature_scale**2 -- geometric similarity, in which '
               'maximum isometric force follows physiological cross-sectional '
               'area. That is a modelling ASSUMPTION, not a measurement on this '
               'subject, and it is a knob precisely so it can be contradicted. '
               'Set it to 1.0 to scale the skeleton while leaving strength '
               'alone; the resulting body is relatively weaker or stronger by '
               'exactly stature_scale**2.'),
        range_basis='Engineering bounds.',
        source=MECHANICAL_MODEL,
        consumers=('ihm.native.model_scaling.scale_model',),
        derived=(),
    ),
    dict(
        name='pelvis_breadth_ratio', unit='dimensionless', kind='continuous',
        range=(0.85, 1.20), default=1.0, status='implemented',
        body='mechanical',
        basis=('Hip-joint-centre separation RELATIVE TO STATURE, as a multiple of '
               "the source subject's. Defining it against stature rather than as "
               'a raw geometric factor is what makes it the quantity NHANES '
               'actually measures, and it is what lets it compose with '
               'stature_m: the materializer solves for the raw medio-lateral '
               'factor that lands the measured ratio on the request. '
               'Realised as a medio-lateral (model z) factor on the pelvis body '
               'alone. It moves '
               'the pelvis mass centre, its wrap objects, the muscle path points '
               'and markers attached to it, and the parent offset frames of the '
               'hip, back and ground_pelvis joints -- so the hip joint centres '
               'separate and the femurs are the same femurs at a wider stance. '
               'Measured effect on the body: hip joint separation goes as the '
               'factor exactly, segment mass fractions shift, and 50 of the 80 '
               'fitted muscle paths change length while 30 do not. That last '
               'fact is the whole reason this parameter needed the refitter: no '
               'single polynomial coefficient factor can represent it.'),
        range_basis=('Engineering bounds around the measured span. The female/male '
                     'ratio of hip circumference over stature is %.4f (d = -0.91), '
                     'the largest shape difference in the NHANES sample after '
                     'size itself.' % FEMALE_OVER_MALE['hip_over_stature']),
        source=MECHANICAL_MODEL,
        consumers=('ihm.native.anisotropic_scaling.scale_model_anisotropic',
                   'scripts/materialize_proportional_variant.py'),
        derived=(),
        limitation=(
            'CORRESPONDENCE. NHANES measures hip CIRCUMFERENCE at the maximum '
            'protrusion of the buttocks; this parameter is a SKELETAL '
            'medio-lateral factor. Using the circumference ratio as the bone '
            'factor assumes soft tissue and bone scale together in that '
            'direction, and they do not have to. Bi-iliac breadth -- the '
            'measurement that would settle it -- is carried by no catalogued '
            'source in this repository, in either sex. The number is therefore '
            'an ESTIMATE with a stated derivation, and the parameter is a '
            'geometric knob that accepts any value.'),
    ),
    dict(
        name='shoulder_breadth_ratio', unit='dimensionless', kind='continuous',
        range=(0.85, 1.20), default=1.0, status='implemented',
        body='mechanical',
        basis=('Acromial-joint-centre separation RELATIVE TO STATURE, as a '
               "multiple of the source subject's, realised as a medio-lateral "
               '(model z) factor on the torso body alone. Measured effect: '
               'acromial separation goes as the solved factor exactly and the '
               'shoulder-over-hip breadth ratio follows.'),
        range_basis='Engineering bounds.',
        source=MECHANICAL_MODEL,
        consumers=('ihm.native.anisotropic_scaling.scale_model_anisotropic',
                   'scripts/materialize_proportional_variant.py'),
        derived=(),
        limitation=(
            'NO MEASURED SEX VALUE EXISTS HERE. BMX_J carries standing height, '
            'weight, upper leg length, upper arm length, arm, waist and hip '
            'circumference -- and no biacromial or bideltoid breadth. Every sex '
            'preset in this schema therefore leaves this parameter at 1.0, and '
            'that is a gap rather than a finding: shoulder:hip is one of the '
            'most visible dimorphic proportions and this repository cannot '
            'quantify it. ANSUR II carries biacromial breadth by sex and is not '
            'catalogued here.'),
    ),
    dict(
        name='leg_length_ratio', unit='dimensionless', kind='continuous',
        range=(0.85, 1.20), default=1.0, status='implemented',
        body='mechanical',
        basis=('(femur + tibia) length RELATIVE TO STATURE, as a multiple of the '
               "source subject's, realised as a uniform factor on each "
               'lower-limb segment -- femur, tibia, patella, talus, calcaneus, '
               'toes, both sides. Self-referential, because the legs are 48% of '
               'stature: shortening them shortens stature too, so the '
               'materializer solves the fixed point on the MEASURED model rather '
               'than assuming the raw factor is the ratio. Measured effect: the '
               'leg-over-stature ratio moves and stature itself does not.'),
        range_basis=('The measured target is the female/male ratio of upper leg '
                     'length over stature, %.5f (d = 0.44): women\'s thighs are '
                     '2.6%% shorter relative to stature. That is a small effect '
                     'and it is reported as small.'
                     % FEMALE_OVER_MALE['upper_leg_over_stature']),
        source=MECHANICAL_MODEL,
        consumers=('ihm.native.anisotropic_scaling.scale_model_anisotropic',
                   'scripts/materialize_proportional_variant.py'),
        derived=(),
        limitation=(
            'ISOTROPIC PER SEGMENT, so it slims the leg as it shortens it. A '
            'long-axis-only femur scale is refused by '
            'ihm.native.anisotropic_scaling: walker_knee_r takes its '
            'roll-glide translations in a frame rotated by (-1.64, 1.45, 1.57) '
            'rad inside femur_r, where componentwise scaling has no meaning. '
            'Separately, NHANES BMXLEG is inguinal crease to proximal tibia and '
            'is NOT the OpenSim hip-centre-to-knee-centre length; the ratio is '
            'transferred between them without a measured correspondence.'),
    ),
    dict(
        name='sex', unit='category', kind='enum',
        domain=('male', 'female'), default='male',
        status='implemented (proportions only; anatomy unchanged)',
        body='both',
        basis=("Measured by scripts/audit_sex_specific_anatomy.py. The body is "
               "male, and specifically so: 32 of the 4,000 segment-bound "
               "entities are male genital anatomy -- the three penile bodies, "
               "both testes, epididymides, deferent ducts, seminal vesicles and "
               "ejaculatory ducts, the prostate, and their vessels -- all bound "
               "to the pelvis segment. Female-specific entities: 0 of 4,000 and "
               "0 of 8,979 in the display atlas. Mammary gland, nipple and "
               "areola: 0 in every set, in EITHER sex; what exists is Z-Anatomy's "
               "'mammary region', a named patch of chest skin. NONE OF THAT IS "
               "CHANGED BY THIS PARAMETER, in either direction."),
        range_basis=(
            "sex='female' now sets four measured PROPORTIONS -- stature, mass, "
            'pelvis_breadth_ratio and leg_length_ratio -- each the female/male '
            'ratio of survey-weighted NHANES means applied to THIS subject, and '
            'each of which changes the running mechanical body in a way that is '
            'measured and gated (scripts/verify_body_proportions.py). It stopped '
            'raising when it started changing something. What it does NOT change '
            'is anatomy: there is still no female-specific entity, no mammary '
            'gland in either sex, and the male genital tract is still present '
            'and still bound to the pelvis. A female body needs new source '
            'geometry, and docs/BODY_PARAMETERS.md carries the acquisition path.'),
        source=ANATOMICAL_PROFILE,
        consumers=('resolve() -> parameters/derived/sex_realisation',
                   'scripts/materialize_proportional_variant.py'),
        limitation=(
            "sex='female' produces MALE ANATOMY AT FEMALE-TYPICAL PROPORTIONS. "
            'It is not a female body and must never be reported as one. Every '
            'resolved record carries a sex_realisation block naming exactly '
            'which quantities moved and which did not, so a caller that logs '
            'the record cannot quote the label without the caveat.'),
    ),
    dict(
        name='age_years', unit='year', kind='continuous',
        range=(18.0, 90.0), default=44.0, status='declared',
        body='anatomical',
        basis=('Inherited from the BioGears StandardMale patient prior; not '
               'measured on this specimen. Carried in '
               'data/derived/canonical/profile.json.'),
        range_basis='Adult range of the inherited prior; not a validated domain.',
        source=ANATOMICAL_PROFILE,
        consumers=(),
        limitation=('Written into the profile record and read by nothing. '
                    "``ihm/assembly/profile.py`` writes only ``Height`` into the "
                    'BioGears patient XML; age does not reach the physiology '
                    'engine.'),
    ),
    dict(
        name='body_fat_fraction', unit='fraction', kind='continuous',
        range=(0.05, 0.50), default=0.21, status='declared',
        body='anatomical',
        basis='Inherited BioGears StandardMale prior; not measured.',
        range_basis='Prior range; not validated.',
        source=ANATOMICAL_PROFILE,
        consumers=(),
        limitation=('Declared only. The anatomical mass 70.7713 kg is composed '
                    'from a voxel partition at sourced per-constituent '
                    'densities, and this fraction is not one of its inputs.'),
    ),
    dict(
        name='environment', unit='category', kind='enum',
        domain=('free', 'supine', 'upright'), default='supine', status='surfaced',
        body='mechanical', scenario=True,
        basis=('Gravity direction only, plus whether the foot contact set is '
               'instantiated. Not a body property; carried here because it is '
               'the same constructor argument and because callers reach for it '
               'in the same breath.'),
        range_basis='The three values the native engine accepts.',
        source=MECHANICAL_MODEL,
        consumers=('NativeMechanicalStream(environment=...)',),
    ),
)

BY_NAME = {p['name']: p for p in PARAMETERS}


class BodyParameterError(ValueError):
    """A body parameter request that the schema refuses."""


def _check_continuous(spec, value):
    value = float(value)
    lo, hi = spec['range']
    if not math.isfinite(value) or not lo <= value <= hi:
        raise BodyParameterError('%s must be a finite value in [%g, %g], got %r'
                                 % (spec['name'], lo, hi, value))
    return value


def sex_presets(sex):
    """Default parameter values for a sex, as measured ratios on THIS subject.

    Not population means substituted for the subject.  The source subject is one
    measured man; the female preset is that man's own measurements multiplied by
    the female/male ratio of survey-weighted NHANES means, which keeps whatever
    is idiosyncratic about him and changes only what the population says differs
    by sex.  Substituting the female mean outright would have silently swapped
    subjects, and this repository's ledger has that error in it already.
    """
    if sex == 'male':
        return {}
    if sex != 'female':
        raise BodyParameterError('No preset for sex=%r' % (sex,))
    return {
        'stature_m': MECHANICAL_STATURE_M * FEMALE_OVER_MALE['stature'],
        'mass_kg': MECHANICAL_TARGET_MASS_KG * FEMALE_OVER_MALE['mass'],
        'pelvis_breadth_ratio': FEMALE_OVER_MALE['hip_over_stature'],
        'leg_length_ratio': FEMALE_OVER_MALE['upper_leg_over_stature'],
        # shoulder_breadth_ratio is deliberately absent. BMX_J carries no
        # biacromial breadth, so there is no measured value to put here, and a
        # guessed one would be the only unmeasured number in the preset.
    }


def anisotropic_factors(resolved):
    """{body: (fx, fy, fz)} for ``ihm.native.anisotropic_scaling``.

    Empty when every proportional parameter is 1.0, which is what lets a
    stature-only request keep taking the exact isotropic path.
    """
    factors = {}
    pelvis = resolved.get('pelvis_breadth_ratio', 1.0)
    if pelvis != 1.0:
        factors['pelvis'] = (1.0, 1.0, pelvis)
    shoulder = resolved.get('shoulder_breadth_ratio', 1.0)
    if shoulder != 1.0:
        factors['torso'] = (1.0, 1.0, shoulder)
    leg = resolved.get('leg_length_ratio', 1.0)
    if leg != 1.0:
        for body in LEG_SEGMENT_BODIES:
            factors[body] = (leg, leg, leg)
    return factors


#: Measured by ``scripts/audit_sex_specific_anatomy.py``.  Held here so the
#: realisation block below states the absence in numbers rather than in prose.
SEX_ANATOMY_COUNTS = {'simulated_entities': 4000, 'display_entities': 8979,
                      'male_specific_simulated': 32, 'male_specific_display': 78,
                      'female_specific_simulated': 0, 'female_specific_display': 0,
                      'mammary_gland_nipple_areola_either_sex': 0}


def sex_realisation(sex, resolved, presets):
    """Exactly what ``sex`` changed and exactly what it did not.

    This block exists so that the label cannot be quoted without the caveat.
    ``docs/DIRECTION.md``: *a scaled male mesh with surfaces bolted onto it is
    not a female body, and must not be labelled one.*
    """
    changed = {}
    for name, value in sorted(presets.items()):
        default = BY_NAME[name]['default']
        changed[name] = {'from': default, 'to': value,
                         'applied': resolved[name] == value,
                         'ratio': (value / default) if default else None}
    return {
        'sex': sex,
        'what_changed': changed if sex != 'male' else {},
        'what_did_not_change': {
            'anatomical_entities': dict(
                SEX_ANATOMY_COUNTS,
                note=('The entity set is identical for both values of this '
                      'parameter. There is no female-specific entity to add and '
                      'the male genital tract is not removed.')),
            'mesh_source': ('BodyParts3D 4.0 is male-only; Z-Anatomy is male and '
                            'not an independent subject; the OpenSim models are '
                            '50th-percentile male. No catalogued source ships a '
                            'female whole-body mesh.'),
            'shoulder_breadth_ratio': ('left at 1.0 in every preset: BMX_J carries '
                                     'no biacromial or bideltoid breadth, so '
                                     'there is no measured value.'),
            'hormones': ('ihm/native/reproductive.py runs a menstrual-cycle '
                         'gonadotropin model whose E2, P4 and inhibin are '
                         'prescribed functions of time, coupled to nothing, on a '
                         'body with testes. This parameter does not reach it.'),
            'hair': ('ihm/assembly/hair_fields.py labels 7 of 21 fields '
                     'androgen-dependent as metadata strings with no hormone '
                     'input and no measured female values.'),
        },
        'honest_summary': (
            'Male anatomy at %s-typical proportions.' % sex
            if sex != 'male' else
            'The source subject: male anatomy at its own measured proportions.'),
    }


def resolve(request=None):
    """Validate a body-parameter request and return the resolved record.

    The returned mapping carries the resolved parameters, the derived scale
    factors the downstream materializers consume, and -- always -- the list of
    limitations attached to whichever parameters were actually set.  A caller
    that logs the record logs its own caveats.
    """
    request = dict(request or {})
    unknown = set(request) - set(BY_NAME)
    if unknown:
        raise BodyParameterError('Unknown body parameters: ' + ', '.join(sorted(unknown)))

    resolved, explicit = {}, sorted(request)
    sex = request.get('sex', BY_NAME['sex']['default'])
    if sex not in BY_NAME['sex']['domain']:
        raise BodyParameterError('sex=%r is outside the declared domain %r. %s'
                                 % (sex, BY_NAME['sex']['domain'],
                                    BY_NAME['sex'].get('limitation', '')))
    presets = sex_presets(sex)
    for spec in PARAMETERS:
        name = spec['name']
        if spec['kind'] == 'enum':
            value = request.get(name, spec['default'])
            if value not in spec['domain']:
                raise BodyParameterError(
                    '%s=%r is outside the declared domain %r. %s'
                    % (name, value, spec['domain'], spec.get('limitation', '')))
            resolved[name] = value
        elif name == 'muscle_force_scale':
            continue                      # depends on stature_scale; filled below
        else:
            # An explicit request always wins; the sex preset only supplies a
            # default, so a caller can ask for a female-proportioned body at a
            # stature of their choosing without the preset silently overriding it.
            fallback = presets.get(name, spec['default'])
            resolved[name] = _check_continuous(spec, request.get(name, fallback))

    stature_scale = resolved['stature_m'] / MECHANICAL_STATURE_M
    resolved['muscle_force_scale'] = (
        _check_continuous(BY_NAME['muscle_force_scale'], request['muscle_force_scale'])
        if 'muscle_force_scale' in request else stature_scale ** 2)

    derived = {
        'stature_scale': stature_scale,
        'mass_scale': resolved['mass_kg'] / MECHANICAL_SOURCE_MASS_KG,
        'bmi_kg_m2': resolved['mass_kg'] / resolved['stature_m'] ** 2,
        'muscle_force_scale': resolved['muscle_force_scale'],
        'anisotropic_factors': anisotropic_factors(resolved),
        'sex_realisation': sex_realisation(sex, resolved, presets),
    }
    limitations = [spec['limitation'] for spec in PARAMETERS
                   if spec.get('limitation') and (spec['name'] in request
                                                  or spec['status'] != 'declared')]
    if derived['anisotropic_factors']:
        limitations.append(
            'This request is ANISOTROPIC: %d bodies carry per-axis factors. The '
            'fitted muscle path polynomials CANNOT be carried across it by a '
            'coefficient scale -- 50 of 80 paths move and 30 do not -- so the '
            'variant must be built by scripts/materialize_proportional_variant.py, '
            'which refits them, and a coefficient-scaled path set beside an '
            'anisotropic model is a body whose muscles belong to a different '
            'skeleton.' % len(derived['anisotropic_factors']))
    limitations.append(
        'stature_m and mass_kg are INDEPENDENT knobs. Nothing constrains the '
        'implied BMI (%.1f kg/m2 here) and no population joint distribution is '
        'consulted, so an inconsistent pair is accepted without complaint.'
        % derived['bmi_kg_m2'])
    limitations.append(
        'The mechanical body and the anatomical body are different objects. '
        'They disagree by %.2f%% in stature and %.1f%% in mass; these '
        'parameters address the mechanical body except where body="anatomical".'
        % (100 * STATURE_DISAGREEMENT, 100 * MASS_DISAGREEMENT))
    return {'schema': 'ihm.body-parameters.v1', 'requested': explicit,
            'parameters': resolved, 'derived': derived, 'limitations': limitations}


def assert_measured_defaults(root):
    """Recompute every measured constant from its artifact; raise on drift.

    This is the guard, not the documentation.  A constant that stops matching
    the file it came from is a silent change of subject.
    """
    from .native.model_scaling import head_marker_height_m
    import xml.etree.ElementTree as ET
    root = Path(root)

    model = ET.parse(root / MECHANICAL_MODEL).getroot()
    mass = sum(float(b.findtext('mass')) for b in model.iter('Body'))
    if abs(mass - MECHANICAL_SOURCE_MASS_KG) > 1e-9:
        raise BodyParameterError('Source model mass drifted: %r vs %r'
                                 % (mass, MECHANICAL_SOURCE_MASS_KG))

    stature = head_marker_height_m(root / MECHANICAL_MODEL)
    if abs(stature - MECHANICAL_STATURE_M) > 1e-9:
        raise BodyParameterError('Measured stature drifted: %r vs %r'
                                 % (stature, MECHANICAL_STATURE_M))

    profile_path = root / ANATOMICAL_PROFILE
    checked_profile = profile_path.exists()
    if checked_profile:
        profile = json.loads(profile_path.read_text())
        for key, constant in (('height_m', ANATOMICAL_STATURE_M),
                              ('mass_kg', ANATOMICAL_MASS_KG)):
            if abs(float(profile[key]) - constant) > 1e-7:
                raise BodyParameterError('Anatomical %s drifted: %r vs %r'
                                         % (key, profile[key], constant))
        if profile['sex'] != BY_NAME['sex']['default']:
            raise BodyParameterError('Anatomical profile sex drifted')
    return {'source_mass_kg': mass, 'stature_m': stature,
            'anatomical_profile_checked': checked_profile,
            'anatomical_profile_sha256': (
                hashlib.sha256(profile_path.read_bytes()).hexdigest()
                if checked_profile else None)}
