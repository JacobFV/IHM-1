"""Sourced tissue densities, and the rule that says what volume a canonical
surface actually holds as tissue.

Two questions decide an entity's mass, and until this module existed only one of
them had an owner.

  density      `scripts/build_body_mechanics.py` gave every non-bone entity
               1000 kg/m3 and every bone 1900 kg/m3. Neither is a tissue
               density: 1000 is water, and 1900 is marrow-free cortical bone,
               which cannot describe whole-bone atlas geometry because that
               geometry encloses the marrow cavity. Both numbers are replaced
               here by the values already sourced in
               `data/derived/tissue-material-candidate-v1/`, so the two lanes
               cite the same publications for the same tissue.

  volume       a CLOSED surface's signed integral is a measurement of what the
               surface encloses. An OPEN surface's signed integral is not. For
               a membrane authored as one open sheet -- a fascia, an
               intermuscular septum, a joint capsule, a serosa -- the signed
               integral is the volume of the structures the membrane WRAPS, not
               the membrane. `left fascia lata` is a 0.169 m2 sheet about a
               millimetre thick, and the shipped builder was carrying it as
               4.086 L, the volume of the thigh inside it. The 84 open
               connective-tissue sheets carried 29.6 L between them, 32% of the
               whole unscaled proxy mass, and that volume is the muscle it
               wraps, counted a second time.

The repository already owns this argument twice, in
`ihm.assembly.anatomy.SURFACE_REGION_TOKENS` ("Their enclosed volume is the skin
they lie on, so counting it as tissue would double-count that skin") and in the
shipped vascular rule (wall = surface area x 0.3 mm, lumen excluded). This
module states it once, for every representation, so a rebuild reaches the same
answer and the rule can be audited in one place.

Tier vocabulary is the repository's: `measured` is reserved for a value measured
on THIS specimen and no density below reaches it; `transferred` is a published
value carried in from a donor cohort or reference table; `assumed` is a value no
retrieved source constrains. Every membrane and viscus-wall thickness here is
`assumed` and is swept in
`scripts/build_tissue_material_assignment_candidate.py`.
"""

# --------------------------------------------------------------------------- sources

SOURCES = {
    'icru44_nist': {
        'url': 'https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html',
        'finding': 'ICRU Report 44 tissue compositions and densities as tabulated by NIST: '
                   'adipose 950, skeletal muscle 1050, whole blood 1060, soft tissue 1060, '
                   'brain grey and white matter 1040 kg/m3. Reference-man table values, not '
                   'measurements on this specimen.'},
    'icrp89': {
        'url': 'https://www.icrp.org/publication.asp?id=ICRP%20Publication%2089',
        'finding': 'ICRP Publication 89 Table 2.20 and para 434/514: whole fresh adult skeleton '
                   '1300 kg/m3, hydrated cortical bone 1900, trabecular bone with marrow 1100, '
                   'skin about 1100 kg/m3. The whole-skeleton value is the one that matches '
                   'whole-bone atlas geometry because that geometry encloses the marrow cavity.'},
    'kanematsu2015': {
        'url': 'https://pubmed.ncbi.nlm.nih.gov/26305733/',
        'finding': 'ICRP 110 reference phantom mass densities, including lung 384 kg/m3 and a '
                   'pooled skin/cartilage/spongiosa category at 1090 kg/m3. The pooled category '
                   'is not a cartilage measurement.'},
}

# Density rows this module is willing to defend, keyed by tissue class. Each row
# is (value_kg_m3, source_id, tier, note). Every note says what the number is and
# what it is not; the `assumed` rows exist because
# `scripts/candidate_tissue_materials.py` searched for a human value for that
# tissue and recorded that it found none.
DENSITY = {
    'adipose': (950., 'icru44_nist', 'transferred',
                'ICRU-44 adipose tissue; ICRP 110 adipose/marrow is also 950. Pure human fat at '
                '37 C is 900; adipose TISSUE is about 80% fat by mass in adults, which is why the '
                'tissue density exceeds the fat density. The same 950 is used by the fill rows of '
                'data/derived/interstitial-composition-prior-v1, so both lanes agree.'),
    'skeletal_muscle': (1050., 'icru44_nist', 'transferred',
                        'ICRU-44 skeletal muscle. Ward and Lieber measured 1112 and 1055 kg/m3 in '
                        'FIXED human muscle; no fresh human value was retrieved.'),
    'skin': (1100., 'icrp89', 'transferred',
             'ICRP 89 para 514, skin about 1.1 g/cm3, the value that report itself uses to convert '
             'mass thickness to linear thickness.'),
    'skeleton_with_marrow': (1300., 'icrp89', 'transferred',
                             'ICRP 89 Table 2.20 whole fresh adult skeleton. This replaces the '
                             '1900 kg/m3 the builder previously used, which is marrow-free '
                             'cortical bone and cannot describe a surface that encloses the '
                             'marrow cavity. Cross-check: the atlas bone geometry integrates to '
                             '5.91 L, so at 1300 the skeleton weighs 7.7 kg against the ICRP 89 '
                             'reference 10.5 kg; the shortfall is atlas bone VOLUME, not density, '
                             'and is recorded rather than absorbed by raising the density back.'),
    'cartilage': (1090., 'kanematsu2015', 'assumed',
                  'no cartilage-specific density was retrieved. ICRP 110 pools skin, cartilage and '
                  'spongiosa at 1090 kg/m3; that is a pooled category, not cartilage, so this row '
                  'is assumed rather than transferred.'),
    'dense_connective': (1060., 'icru44_nist', 'assumed',
                         'ICRU-44 soft tissue as the class fallback for tendon, ligament, fascia '
                         'and joint capsule. data/derived/tissue-material-candidate-v1 records '
                         'that no human whole-tendon, ligament or fascia density was retrieved; '
                         'the fascia constituent of interstitial-composition-prior-v1 uses the '
                         'same 1060, so the two lanes agree on the fallback.'),
    'nervous_tissue': (1040., 'icru44_nist', 'transferred',
                       'ICRU-44 brain grey and white matter. No peripheral-nerve density was '
                       'retrieved; central and peripheral share this row.'),
    'blood_and_vessel_wall': (1060., 'icru44_nist', 'transferred',
                              'ICRU-44 whole blood / soft tissue. No vessel-wall-specific value '
                              'was found; the vascular entities are lumen-shaped and carry a '
                              '0.3 mm wall.'),
    'lung': (384., 'kanematsu2015', 'transferred',
             'ICRP 110 reference lung, inflated. Cross-check: healthy whole-lung mass at '
             'functional residual capacity is 997 +/- 133 g by CT.'),
    'soft_tissue': (1060., 'icru44_nist', 'transferred',
                    'ICRU-44 soft tissue, the highest density in that table for any soft tissue. '
                    'Used for every soft organ with no organ-specific row.'),
    'lymphoid': (1030., 'icru44_nist', 'assumed',
                 'ICRU-44 soft tissue band, lower bound. No lymph-node density was retrieved; '
                 'this matches the lymphatic row of the candidate table system densities.'),
}

# --------------------------------------------------------------------------- representation

# A membrane is an anatomical sheet: it wraps, lines or partitions, and its own
# tissue is thin. Matched on the name because the canonical role vocabulary has
# no membrane role, exactly as `physical_role` matches 'intervertebral disk'.
# `tensor fasciae latae` is a muscle whose name contains 'fascia', so the tokens
# are head nouns and the muscle exclusion below is explicit.
MEMBRANE_TOKENS = ('fascia', 'aponeuros', 'interosseous membrane', 'obturator membrane',
                   'thyrohyoid membrane', 'cricothyroid membrane', 'perineal membrane',
                   'tectorial membrane', 'atlanto-occipital membrane', 'atlantoaxial membrane',
                   'articular capsule', 'joint capsule', 'capsule of', 'retinacul',
                   'intermuscular septum', 'intermuscular septa', 'crural septum',
                   'mesentery', 'mesocolon', 'mesoappendix', 'mesosalpinx', 'mesovarium',
                   'peritoneum', 'pleura', 'pericardium', 'omentum', 'dura mater', 'arachnoid',
                   'pia mater', 'tentorium', 'falx', 'diaphragma sellae', 'tendon sheath',
                   'synovial sheath', 'galea', 'periosteum', 'perichondrium', 'linea alba',
                   'iliotibial tract', 'tunica vaginalis', 'conjunctiva')

# Names that contain a membrane token but are not membranes.
MEMBRANE_EXCLUSIONS = ('tensor fasciae latae', 'tensor of fascia lata', 'fat pad',
                       'membranous part of', 'membranous urethra', 'membranous septum',
                       'tympanic membrane', 'mucous membrane', 'membranous labyrinth',
                       'basilar membrane', 'membrane of tympanum')

# A hollow viscus is authored in BodyParts3D as ONE closed surface enclosing wall
# and lumen together, so its enclosed volume includes gas or chyme that is not
# tissue. Measured on this specimen in
# `data/derived/interstitial-matrix-v1/lumen.json`: 120 of 225 hollow organs
# admit an inscribed sphere over 3 mm inside their own surface, stomach 39.1 mm.
HOLLOW_VISCUS_TOKENS = ('trachea', 'bronch', 'larynx', 'pharynx', 'esophagus', 'oesophagus',
                        'stomach', 'duodenum', 'jejunum', 'ileum', 'caecum', 'cecum', 'colon',
                        'rectum', 'anal canal', 'vermiform appendix', 'gall bladder', 'gallbladder',
                        'bile duct', 'cystic duct', 'hepatic duct', 'pancreatic duct',
                        'urinary bladder', 'ureter', 'urethra', 'uterus', 'uterine tube',
                        'vagina', 'seminal vesicle', 'ductus deferens', 'epididymis')

# Thicknesses. Every one is `assumed`: no measurement in this repository
# constrains any of them and no external measurement was obtainable in the
# session that wrote this module. The sweep bounds are the sensitivity envelope
# reported beside the assignment, not a probability interval.
THICKNESS = {
    'membrane': (0.0010, 'assumed', (0.0005, 0.0020),
                 'one millimetre of sheet. Human deep fascia is conventionally described as under '
                 'two millimetres and joint capsule and serosa are thinner still, but NO retrieved '
                 'source supplies a number here, so this is an assumed prior swept over '
                 '[0.5, 2.0] mm. Applied only where the surface is OPEN, because a closed sheet '
                 'already measures its own thickness.'),
    'hollow_viscus_wall': (0.0030, 'assumed', (0.0020, 0.0050),
                           'three millimetres of gut, airway or urinary wall. No retrieved source '
                           'supplies a number here; swept over [2, 5] mm. Applied as an upper '
                           'bound, so an organ whose enclosed volume is already thinner than its '
                           'own wall keeps its enclosed volume.'),
    'vessel_wall': (0.0003, 'assumed', (0.0002, 0.0005),
                    'the 0.3 mm wall the shipped builder already used for every vascular entity, '
                    'restated here so every thickness prior has one owner.'),
    'skin_boundary_carrier': (0.000001, 'numerical', (0.000001, 0.000001),
                              'one micrometre numerical carrier for the skin boundary. The '
                              'physical layer masses belong to the three skin-layer entities.'),
}


def is_membrane(name, role):
    """True where the entity is an anatomical sheet rather than a solid."""
    name = name.lower()
    if any(x in name for x in MEMBRANE_EXCLUSIONS):
        return False
    if 'intermuscular septum' in name or 'intermuscular septa' in name:
        return True
    if role in ('vascular', 'nerve', 'lymph_node_group', 'rigid_bone', 'skin', 'skin_layer',
                'surface_region', 'lymphatic_network', 'fluid_cavity'):
        return False
    return any(x in name for x in MEMBRANE_TOKENS)


def is_hollow_viscus(name, role):
    """True where the surface encloses a lumen as well as a wall."""
    name = name.lower()
    if role not in ('soft_organ', 'muscle'):
        return False
    if is_membrane(name, role):
        return False
    return any(x in name for x in HOLLOW_VISCUS_TOKENS)


def tissue_class(name, role):
    """The density class this entity's material is drawn from."""
    name = name.lower()
    if role == 'skin_layer':
        # The hypodermis IS subcutaneous adipose. It is the ONLY entity in the
        # canonical model that can carry body fat at an adipose density: a search
        # of the 4000 entity names finds no other adipose geometry beyond two
        # 0.6 microlitre infrapatellar fat pads. Giving it 1000 kg/m3 left the
        # whole subcutaneous depot with no adipose material anywhere.
        return 'adipose' if 'hypodermis' in name else 'skin'
    if role == 'skin':
        return 'skin'
    if role == 'rigid_bone':
        return 'skeleton_with_marrow'
    if role == 'muscle':
        return 'skeletal_muscle'
    if role == 'cartilage':
        return 'cartilage'
    if role in ('tendon', 'ligament', 'connective_tissue'):
        return 'dense_connective'
    if role == 'vascular':
        return 'blood_and_vessel_wall'
    if role == 'nerve':
        return 'nervous_tissue'
    if role == 'lymph_node_group':
        return 'lymphoid'
    if 'lung' in name:
        return 'lung'
    return 'soft_tissue'


def density(name, role):
    """(value_kg_m3, record) for this entity's tissue class."""
    key = tissue_class(name, role)
    value, source, tier, note = DENSITY[key]
    return value, {'tissue_class': key, 'value_kg_m3': value, 'source': source, 'tier': tier,
                   'note': note, 'measured_on_this_specimen': False}


def material_volume(name, role, prior_volume_m3, surface_area_m2, watertight,
                    thickness=None):
    """The volume this surface holds as tissue, and the rule that decided it.

    `prior_volume_m3` is what the caller's own geometric fallback produced -- a
    closed signed integral where the surface is closed, an open signed integral
    or a bounds ellipsoid where it is not. `watertight` says which. The rules
    below only ever REDUCE that number, and only where the surface's own
    topology says the enclosed volume is not the entity's tissue.
    """
    t = dict(THICKNESS) if thickness is None else {**THICKNESS, **thickness}
    area = float(surface_area_m2 or 0.)
    prior = float(prior_volume_m3 or 0.)
    if role == 'skin':
        value = area * t['skin_boundary_carrier'][0] or prior
        return value, {'rule': 'skin_boundary_carrier', 'thickness_m': t['skin_boundary_carrier'][0],
                       'thickness_tier': t['skin_boundary_carrier'][1],
                       'basis': '1 micrometer numerical carrier for skin boundary; physical layer '
                                'masses belong to skin-layer entities'}
    if role == 'vascular':
        value = area * t['vessel_wall'][0] or prior
        return value, {'rule': 'vessel_wall', 'thickness_m': t['vessel_wall'][0],
                       'thickness_tier': t['vessel_wall'][1],
                       'basis': 'surface area times assumed 0.3 mm wall thickness; lumen excluded'}
    if is_membrane(name, role):
        if watertight:
            return prior, {'rule': 'closed_membrane_measured', 'thickness_m': None,
                           'thickness_tier': None,
                           'basis': 'the surface is closed, so its signed integral measures the '
                                    'sheet itself and no thickness prior is needed'}
        sheet = area * t['membrane'][0]
        value = min(prior, sheet) if prior > 0 else sheet
        return value, {'rule': 'open_membrane_sheet', 'thickness_m': t['membrane'][0],
                       'thickness_tier': t['membrane'][1],
                       'open_surface_signed_integral_m3': prior,
                       'basis': 'the surface is OPEN, so its signed integral is the volume of the '
                                'structures this sheet wraps, not the sheet. Tissue volume is '
                                'surface area times an assumed membrane thickness, taken as an '
                                'upper bound against the open-surface prior'}
    if is_hollow_viscus(name, role):
        wall = area * t['hollow_viscus_wall'][0]
        if wall < prior:
            return wall, {'rule': 'hollow_viscus_wall', 'thickness_m': t['hollow_viscus_wall'][0],
                          'thickness_tier': t['hollow_viscus_wall'][1],
                          'enclosed_volume_including_lumen_m3': prior,
                          'lumen_volume_m3': prior - wall,
                          'basis': 'a BodyParts3D viscus is one closed surface enclosing wall AND '
                                   'lumen, so the enclosed volume carries gas or chyme at tissue '
                                   'density. Tissue volume is surface area times an assumed wall '
                                   'thickness; the remainder is lumen'}
        return prior, {'rule': 'hollow_viscus_thinner_than_its_wall', 'thickness_m': None,
                       'thickness_tier': None,
                       'basis': 'the enclosed volume is already below one wall thickness, so no '
                                'lumen is subtracted'}
    return prior, {'rule': 'solid', 'thickness_m': None, 'thickness_tier': None,
                   'basis': 'the enclosed volume is taken as tissue'}
