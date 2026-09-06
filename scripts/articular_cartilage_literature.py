#!/usr/bin/env python3
"""Sourced per-joint articular cartilage thickness, one layer per bone surface.

Every value is TRANSFERRED literature measured on other subjects. Nothing here
was measured on the BodyParts3D reference specimen the canonical atlas carries,
so every row is `subject_calibrated: false`. Rows are tagged:

  measured_literature  a published mean/representative value for that surface
  interpolated         a published value for a neighbouring surface of the same
                       joint, or a published range midpoint, applied to a
                       surface the source did not itself report
  assumed              no usable published value found; the number is an
                       engineering placeholder and is flagged as such

Thicknesses are per articular surface in millimetres. The joint target gap used
downstream is the sum of the two surfaces' thicknesses.
"""

SOURCES = {
    'shepherd1999': {
        'citation': 'Shepherd DET, Seedhom BB. Thickness of human articular cartilage in joints of the '
                    'lower limb. Ann Rheum Dis 1999;58(1):27-34.',
        'doi': '10.1136/ard.58.1.27', 'pmid': '10343537',
        'method': 'needle-probe indentation, 11 cadaveric sets of matched ankle/knee/hip',
        'subjects': 11, 'population': 'cadaveric adult'},
    'athanasiou1995': {
        'citation': 'Athanasiou KA, Niederauer GG, Schenck RC. Biomechanical topography of human ankle '
                    'cartilage. Ann Biomed Eng 1995;23(5):697-704.',
        'doi': '10.1007/BF02584467', 'pmid': '7503470',
        'method': 'creep indentation at 14 sites, 14 cadaveric ankles (7 pairs)',
        'subjects': 14, 'population': 'cadaveric adult'},
    'wyler2007': {
        'citation': 'Wyler A, Bousson V, Bergot C, et al. Hyaline cartilage thickness in radiographically '
                    'normal cadaveric hips: comparison of spiral CT arthrographic and macroscopic '
                    'measurements. Radiology 2007;242(2):441-449.',
        'doi': '10.1148/radiol.2422051393', 'pmid': '17255415',
        'method': 'spiral multidetector CT arthrography vs anatomic slices, 12 cadaveric hips',
        'subjects': 12, 'population': 'cadaveric adult, mean age 76.5 y'},
    'carlson2023': {
        'citation': 'Carlson CG, Chen A, Patterson K, Ablove RH. Glenohumeral cartilage thickness: '
                    'implications in prosthetic design and osteochondral allograft transplantation. '
                    'Cartilage 2023.',
        'doi': '10.1177/19476035231154504', 'pmid': None,
        'method': 'sectioned fresh cadaveric shoulders, 5 points per 5 mm coronal section',
        'subjects': 16, 'population': 'cadaveric adult'},
    'jbjsoa2019elbow': {
        'citation': 'Regional distribution of articular cartilage thickness in the elbow joint: a '
                    '3-dimensional study in elderly humans. JBJS Open Access 2019;4(3):e0011.',
        'doi': '10.2106/JBJS.OA.19.00011', 'pmid': '31592501',
        'method': '3-D laser surface scanning before and after cartilage removal, 20 cadaveric elbows',
        'subjects': 20, 'population': 'cadaveric elderly'},
    'graichen1998elbow': {
        'citation': 'Noninvasive analysis of cartilage volume and cartilage thickness in the human elbow '
                    'joint using MRI. (German) 1998.',
        'doi': None, 'pmid': '9728274',
        'method': 'MRI segmentation of elbow cartilage volume and thickness',
        'subjects': None, 'population': 'in vivo adult'},
    'yamabe2013radius': {
        'citation': 'Articular cartilage thickness at the distal radius: a cadaveric study. '
                    'J Hand Surg Am 2013.',
        'doi': None, 'pmid': '23810572',
        'method': 'direct cadaveric measurement of distal radial articular cartilage',
        'subjects': None, 'population': 'cadaveric adult'},
    'radius3d2020': {
        'citation': 'Cartilage and subchondral bone distributions of the distal radius: a 3-dimensional '
                    'analysis using cadavers. 2020.',
        'doi': None, 'pmid': '32860992',
        'method': '3-D reconstruction of cadaveric distal radius cartilage',
        'subjects': None, 'population': 'cadaveric adult'},
    'spannow2010': {
        'citation': 'Spannow AH, et al. Ultrasound measurement of joint cartilage thickness in large and '
                    'small joints in healthy children.',
        'doi': None, 'pmid': None,
        'method': 'B-mode ultrasound, healthy children',
        'subjects': None, 'population': 'HEALTHY CHILDREN - not adult; used only as an order-of-magnitude '
                                        'anchor for hand joints'},
    'si2002': {
        'citation': 'Sacral and iliac articular cartilage thickness and cellularity: relationship to '
                    'subchondral bone end-plate thickness and cancellous bone density. '
                    'Rheumatology 2002;41(4):375-380.',
        'doi': '10.1093/rheumatology/41.4.375', 'pmid': None,
        'method': 'histomorphometry of cadaveric sacroiliac joints',
        'subjects': None, 'population': 'cadaveric adult'},
    'tmj1996': {
        'citation': 'Morphometric investigation of condylar cartilage and disc thickness in the human '
                    'temporomandibular joint. 1996.',
        'doi': None, 'pmid': '8835815',
        'method': 'morphometry of cadaveric TMJ condylar cartilage and disc',
        'subjects': None, 'population': 'cadaveric adult'},
    'womack2011facet': {
        'citation': 'Womack W, et al. Finite element lumbar spine facet contact parameter predictions are '
                    'affected by the cartilage thickness distribution and initial joint gap size. '
                    'J Biomech Eng 2011.',
        'doi': None, 'pmid': '21744929',
        'method': 'serial sectioning of cadaveric lumbar facets, 3-D cartilage thickness maps',
        'subjects': None, 'population': 'cadaveric adult'},
}

# joint_label -> {bone_kind or '*': (thickness_mm, tier, source_key, note)}
THICKNESS_MM = {
    'tibiofemoral': {'femur': (2.15, 'measured_literature', 'shepherd1999', 'femoral condylar mean'),
                     'tibia': (2.30, 'interpolated', 'shepherd1999',
                               'mean of meniscus-covered 2.01 and uncovered 2.59')},
    'patellofemoral': {'patella': (2.50, 'interpolated', 'shepherd1999',
                                   'patellar cartilage is the thickest in the body; midpoint of the '
                                   'commonly quoted 2-3 mm band'),
                       'femur': (2.15, 'measured_literature', 'shepherd1999', 'trochlear/condylar mean')},
    'coxofemoral': {'femur': (1.43, 'interpolated', 'wyler2007', 'midpoint of 0.32-2.53 mm femoral head'),
                    'hip_bone': (2.04, 'interpolated', 'wyler2007', 'midpoint of 0.95-3.13 mm acetabulum')},
    'talocrural': {'tibia': (1.20, 'interpolated', 'athanasiou1995',
                             'distal tibial plafond, within the 0.95-1.45 mm ankle band'),
                   'fibula': (0.95, 'measured_literature', 'athanasiou1995', 'distal fibula, thinnest site'),
                   'talus': (1.45, 'measured_literature', 'athanasiou1995', 'posterolateral talar ridge, '
                                                                            'thickest talar site')},
    'subtalar': {'*': (1.00, 'assumed', 'athanasiou1995',
                       'no subtalar-specific source found; transferred down from the tibiotalar band')},
    'talonavicular': {'*': (1.00, 'assumed', 'athanasiou1995', 'as subtalar')},
    'calcaneocuboid': {'*': (0.90, 'assumed', 'athanasiou1995', 'as subtalar, one step thinner')},
    'cubonavicular': {'*': (0.80, 'assumed', 'athanasiou1995', 'small tarsal facet placeholder')},
    'cuneonavicular': {'*': (0.80, 'assumed', 'athanasiou1995', 'small tarsal facet placeholder')},
    'intertarsal_distal': {'*': (0.70, 'assumed', 'athanasiou1995', 'small tarsal facet placeholder')},
    'tarsometatarsal': {'*': (0.70, 'assumed', 'athanasiou1995', 'small tarsal facet placeholder')},
    'intermetatarsal': {'*': (0.60, 'assumed', 'athanasiou1995', 'small tarsal facet placeholder')},
    'metatarsophalangeal': {'*': (0.70, 'assumed', 'spannow2010',
                                  'no adult MTP source found; taken one notch below the paediatric MCP '
                                  'ultrasound value')},
    'metatarsosesamoid': {'*': (0.60, 'assumed', 'spannow2010', 'placeholder')},
    'interphalangeal_hallux': {'*': (0.50, 'assumed', 'spannow2010', 'placeholder')},
    'interphalangeal_foot_proximal': {'*': (0.45, 'assumed', 'spannow2010', 'placeholder')},
    'interphalangeal_foot_distal': {'*': (0.40, 'assumed', 'spannow2010', 'placeholder')},
    'glenohumeral': {'humerus': (1.77, 'measured_literature', 'carlson2023', 'central humeral head'),
                     'scapula': (2.28, 'interpolated', 'carlson2023',
                                 'glenoid: mean of central 1.69, superior 2.61, inferior 2.53')},
    'acromioclavicular': {'*': (1.50, 'assumed', 'carlson2023',
                                'fibrocartilaginous AC surfaces; no direct source retrieved')},
    'sternoclavicular': {'*': (1.50, 'assumed', 'carlson2023',
                               'fibrocartilaginous SC surfaces plus an articular disc; placeholder')},
    'humeroulnar': {'humerus': (1.27, 'measured_literature', 'jbjsoa2019elbow',
                                'thickest distal-humeral region, capitellotrochlear groove'),
                    'ulna': (0.90, 'measured_literature', 'graichen1998elbow', 'proximal ulna MRI mean')},
    'humeroradial': {'humerus': (1.40, 'measured_literature', 'graichen1998elbow', 'capitulum humeri MRI mean'),
                     'radius': (1.10, 'interpolated', 'graichen1998elbow',
                                'radial head; between the reported 0.9 mm ulna and 1.4 mm capitulum')},
    'radioulnar': {'*': (0.80, 'assumed', 'graichen1998elbow',
                         'proximal and distal radioulnar facets; placeholder')},
    'radiocarpal': {'radius': (0.70, 'measured_literature', 'yamabe2013radius',
                               'scaphoid fossa 0.7, lunate fossa 0.6, ridge 0.8; overall mean 0.6'),
                    '*': (0.70, 'interpolated', 'radius3d2020',
                          'proximal carpal row assumed to mirror its radial fossa')},
    'ulnocarpal': {'*': (0.60, 'assumed', 'yamabe2013radius',
                         'the ulnocarpal space is bridged by the TFCC, not by bone-on-bone cartilage')},
    'midcarpal': {'*': (0.60, 'assumed', 'yamabe2013radius', 'carpal facet placeholder')},
    'intercarpal_proximal': {'*': (0.55, 'assumed', 'yamabe2013radius', 'carpal facet placeholder')},
    'intercarpal_distal': {'*': (0.55, 'assumed', 'yamabe2013radius', 'carpal facet placeholder')},
    'carpometacarpal': {'*': (0.55, 'assumed', 'yamabe2013radius', 'carpal facet placeholder')},
    'carpometacarpal_thumb': {'*': (0.70, 'assumed', 'yamabe2013radius',
                                    'saddle joint, thicker than the rigid CMC 2-5 placeholder')},
    'intermetacarpal': {'*': (0.50, 'assumed', 'yamabe2013radius', 'carpal facet placeholder')},
    'metacarpophalangeal': {'*': (0.65, 'assumed', 'spannow2010',
                                  'paediatric ultrasound gives 1.52 mm; adult MCP cartilage is reported '
                                  'well under 1 mm, so this is an adult-scaled placeholder')},
    'interphalangeal_hand_proximal': {'*': (0.50, 'assumed', 'spannow2010',
                                            'paediatric PIP 0.73 mm scaled to adult')},
    'interphalangeal_hand_distal': {'*': (0.40, 'assumed', 'spannow2010', 'placeholder')},
    'interphalangeal_thumb': {'*': (0.50, 'assumed', 'spannow2010', 'placeholder')},
    'sacroiliac': {'sacrum': (1.81, 'measured_literature', 'si2002', 'sacral hyaline cartilage'),
                   'hip_bone': (0.80, 'measured_literature', 'si2002', 'iliac fibrocartilage')},
    'zygapophyseal': {'*': (0.50, 'assumed', 'womack2011facet',
                            'lumbar facet cartilage maps exist but no single mean was retrieved; '
                            '0.5 mm per surface is a placeholder')},
    'atlanto_occipital': {'*': (0.70, 'assumed', 'womack2011facet', 'placeholder')},
    'atlantoaxial': {'*': (0.70, 'assumed', 'womack2011facet', 'placeholder')},
    'costovertebral': {'*': (0.40, 'assumed', 'womack2011facet',
                             'costovertebral and costotransverse facets; placeholder')},
    'temporomandibular': {'*': (0.30, 'measured_literature', 'tmj1996',
                                'TMJ surfaces carry thin fibrocartilage (temporal 0.05-0.4 mm, condyle '
                                'similar) plus a 1-3 mm disc that is a separate body, not an offset shell')},
}

NON_ARTICULAR_CLASSES = ('dental', 'non_articular', 'other')

# Interfaces that carry no synovial cavity and no offset cartilage shell, but that
# must stay in contact for the skeleton to hold together. They enter the solve as
# zero-gap constraints and are reported separately from the cartilage joints.
BONDED_CLASSES = ('suture', 'amphiarthrosis', 'symphysis')
BONDED_TARGET_MM = 0.0
BONDED_NOTE = ('Cranial sutures, discovertebral interfaces and symphyses are modelled as bonded at zero '
               'gap. This is a structural constraint on the solve, not a literature cartilage thickness; '
               'the true interposed layers (sutural ligament, cartilaginous endplate, interpubic disc) are '
               'sub-millimetre to a few millimetres and are not modelled here.')


def target_gap_mm(label, kind_a, kind_b):
    """(target_gap_mm, tier, [rows]) for a joint, or None when unmodelled."""
    table = THICKNESS_MM.get(label)
    if table is None:
        return None
    rows, total, tiers = [], 0., []
    for kind in (kind_a, kind_b):
        value = table.get(kind, table.get('*'))
        if value is None:
            return None
        thickness, tier, source, note = value
        rows.append({'bone_kind': kind, 'thickness_mm': thickness, 'tier': tier,
                     'source': source, 'note': note})
        total += thickness
        tiers.append(tier)
    worst = 'assumed' if 'assumed' in tiers else ('interpolated' if 'interpolated' in tiers
                                                  else 'measured_literature')
    return round(total, 4), worst, rows
