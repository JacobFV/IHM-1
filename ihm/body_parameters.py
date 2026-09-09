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
        range_basis=('Engineering bounds, not a population interval. 1.40-2.05 m '
                     'is roughly the 0.1st to 99.9th percentile of adult human '
                     'stature, but no anthropometric table in this repository '
                     'defines it and none is cited. What the bounds actually '
                     'guard is the scaling itself: outside roughly 0.8x-1.15x '
                     'the source subject the muscle path polynomials are being '
                     'extrapolated far past the poses they were fitted on.'),
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
        range_basis=('Engineering bounds. The native engine accepts any positive '
                     'mass; these bounds only keep a typo from producing a body '
                     'that silently integrates.'),
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
        name='sex', unit='category', kind='enum',
        domain=('male',), default='male', status='declared',
        body='both',
        basis=("The anatomical body is the BodyParts3D adult male atlas and the "
               "mechanical body is a male Rajagopal subject. A search of all "
               "~4,000 anatomical entities for penis, testis, prostate, uterus, "
               "ovary, vagina, breast, mammary, labia and clitoris returns "
               "nothing: the simulated body has ZERO sex-specific entities. The "
               "domain of this parameter has one member for that reason and no "
               "other."),
        range_basis=('Widening the domain requires new source geometry, not new '
                     'code. docs/BODY_PARAMETERS.md sets out what specifically.'),
        source=ANATOMICAL_PROFILE,
        consumers=(),
        limitation=("Writing sex='female' would change a string in a JSON file "
                    'and nothing else. The schema refuses the value rather than '
                    'accepting it and quietly meaning nothing.'),
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
            resolved[name] = _check_continuous(spec, request.get(name, spec['default']))

    stature_scale = resolved['stature_m'] / MECHANICAL_STATURE_M
    resolved['muscle_force_scale'] = (
        _check_continuous(BY_NAME['muscle_force_scale'], request['muscle_force_scale'])
        if 'muscle_force_scale' in request else stature_scale ** 2)

    derived = {
        'stature_scale': stature_scale,
        'mass_scale': resolved['mass_kg'] / MECHANICAL_SOURCE_MASS_KG,
        'bmi_kg_m2': resolved['mass_kg'] / resolved['stature_m'] ** 2,
        'muscle_force_scale': resolved['muscle_force_scale'],
    }
    limitations = [spec['limitation'] for spec in PARAMETERS
                   if spec.get('limitation') and (spec['name'] in request
                                                  or spec['status'] != 'declared')]
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
