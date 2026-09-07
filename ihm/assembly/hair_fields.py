"""Hair as anatomy: named integumentary fields on the measured outer body envelope.

Regions are Terminologia-Anatomica body surface labels carried verbatim from the
Z-Anatomy atlas and registered into the canonical frame, not axis-aligned boxes.
Densities and morphology are literature values with an explicit evidence tier;
nothing here is a rendering prior.
"""
from __future__ import annotations
import numpy as np
from scipy.spatial import cKDTree

# Hair-bearing footprints the atlas authors as their own geometry. Preferred over
# region lexicon membership where they exist, because they are the hair itself.
SHELL_FIELDS = {
    'scalp': ('Hairs of head',),
    'eyebrow': ('Hairs of eyebrow.l', 'Hairs of eyebrow.r'),
    'eyelash': ('Eyelashes.l', 'Eyelashes.r'),
    'pubic': ('Pubic hairs',),
}
SHELL_STANDOFF_M = .006
APPENDAGE_SHELLS = ('Hairs of head', 'Hairs of eyebrow.l', 'Hairs of eyebrow.r', 'Eyelashes.l', 'Eyelashes.r',
                    'Pubic hairs', 'Eyebrow.l', 'Eyebrow.r', 'Nail plate.l', 'Nail plate.r',
                    'Nail plate (foot).l', 'Nail plate (foot).r')

# Terminologia Anatomica surface region -> hair field. Every exported region is
# claimed exactly once; unclaimed regions are a build error, not a silent gap.
LEXICON = {
    'facial_vellus': ['Frontal region', 'Orbital region', 'Infra-orbital region', 'Nasal region', 'Nasolabial sulcus',
                      'Philtrum', 'Mentolabial sulcus', 'Angle of mouth', 'Oral region', 'Zygomatic region'],
    'beard': ['Buccal region', 'Mental region', 'Parotideomasseteric region', 'Submental triangle',
              'Submandibular triangle'],
    'auricular': ['Auricular region', 'Helix', 'Antihelix', 'Crura of antihelix', 'Scapha', 'Triangular fossa',
                  'Tragus', 'Antitragus', 'Intertragic incisure', 'Lobule of auricle', 'Concha of auricle',
                  'Cavity of concha', 'Cymba conchae', 'Apex of auricle', 'Auricular tubercle',
                  'Anterior notch of auricle', 'Posterior auricular groove', 'Eminentia conchae',
                  'Eminentia scaphae', 'Eminentia fossae triangularis', 'Fossa antihelica', 'Mastoid region'],
    'scalp': ['Parietal region', 'Occipital region', 'Temporal region'],
    'neck': ['Carotid triangle', 'Muscular triangle', 'Sternocleidomastoid region', 'Lateral region of neck',
             'Posterior region of neck', 'Greater supraclavicular fossa', 'Lesser supraclavicular fossa'],
    'chest': ['Presternal region', 'Pectoral region', 'Mammary region', 'Inframammary region',
              'Infraclavicular fossa', 'Deltopectoral triangle', 'Lateral region of thorax'],
    'abdomen': ['Epigastric region', 'Umbilical region', 'Umbilicus', 'Hypochondriac region',
                'Lateral region of abdomen', 'Hypogastric region', 'Inguinal region'],
    'back': ['Interscapular region', 'Scapular region', 'Infrascapular region', 'Vertebral region', 'Lumbar region',
             'Sacral region', 'Triangle of auscultation'],
    'gluteal': ['Gluteal region', 'Gluteal fold', 'Hip region'],
    'arm': ['Anterior region of arm', 'Posterior region of arm', 'Deltoid region', 'Lateral bicipital groove',
            'Medial bicipital groove', 'Anterior region of elbow', 'Posterior region of elbow', 'Cubital fossa'],
    'forearm': ['Anterior region of forearm', 'Posterior region of forearm', 'Lateral border of forearm',
                'Medial border of forearm', 'Anterior region of wrist', 'Posterior region of wrist'],
    'hand_dorsum': ['Dorsum of hand', 'Dorsal surfaces of digits of hand', 'Radial foveola'],
    'thigh': ['Anterior region of thigh', 'Posterior region of thigh', 'Femoral triangle', 'Anterior region of knee',
              'Posterior region of knee', 'Popliteal fossa'],
    'leg': ['Anterior region of leg', 'Posterior region of leg', 'Anterior region of ankle', 'Medial malleolus',
            'Lateral malleolus', 'Medial retromalleolar region', 'Lateral retromalleolar region'],
    'foot_dorsum': ['Dorsum of foot', 'Dorsal surfaces of digits of foot', 'Medial border of foot',
                    'Lateral border of foot', 'Metatarsal region'],
    'pubic': [],
    'perineal': ['Urogenital region', 'Anal region'],
    'glabrous': ['Palm', 'Palmar surfaces of digits of hand', 'Sole', 'Plantar surfaces of digits of foot',
                 'Heel region', 'Hallucial eminence', 'Medial part of longitudinal arch of foot',
                 'Lateral part of longitudinal arch of foot', 'Proximal transverse arch of foot',
                 'Distal transverse arch of foot', 'Perionyx', 'Perionyx (foot)', 'Tubercle of upper lip',
                 'Labial commissure'],
}

FIELDS = {
    'scalp': ('scalp hair', 'terminal', 'bilateral'),
    'eyebrow': ('eyebrow hair', 'terminal', 'bilateral'),
    'eyelash': ('eyelash', 'terminal', 'bilateral'),
    'beard': ('beard and moustache hair', 'terminal_androgen_dependent', 'bilateral'),
    'facial_vellus': ('facial vellus hair', 'vellus', 'bilateral'),
    'auricular': ('auricular hair', 'mixed', 'bilateral'),
    'neck': ('neck hair', 'mixed', 'bilateral'),
    'axillary': ('axillary hair', 'terminal_androgen_dependent', 'bilateral'),
    'chest': ('chest hair', 'mixed_androgen_dependent', 'bilateral'),
    'abdomen': ('abdominal hair', 'mixed_androgen_dependent', 'bilateral'),
    'back': ('back hair', 'mixed_androgen_dependent', 'bilateral'),
    'gluteal': ('gluteal hair', 'mixed', 'bilateral'),
    'arm': ('arm hair', 'mixed', 'bilateral'),
    'forearm': ('forearm hair', 'mixed', 'bilateral'),
    'hand_dorsum': ('dorsal hand and digit hair', 'mixed', 'bilateral'),
    'thigh': ('thigh hair', 'mixed', 'bilateral'),
    'leg': ('leg hair', 'mixed', 'bilateral'),
    'foot_dorsum': ('dorsal foot and digit hair', 'mixed', 'bilateral'),
    'pubic': ('pubic hair', 'terminal_androgen_dependent', 'median'),
    'perineal': ('perineal and anal hair', 'terminal_androgen_dependent', 'median'),
    'glabrous': ('glabrous integument', 'none', 'bilateral'),
}

ROLE_VOCABULARY = {
    'exact': False, 'role': 'skin_layer', 'anatomical_role': 'hair-bearing integumentary field',
    'basis': 'nearest available canonical term; both are non-volumetric fields carried on the acquired skin surface, '
             'neither is a volumetric organ',
    'rejected': ['soft_organ: the canonical atlas assigns it to bp3d-FJ2813 and bp3d-FJ2815, but a follicle population '
                 'field is not an organ; that assignment is carried here as a recorded defect, not reproduced',
                 'connective_tissue: hair is a keratinised epithelial appendage, not connective tissue',
                 'skin: reserved for the single acquired whole-body skin surface bp3d-FJ2810'],
    'vocabulary_gap': True, 'proposed_term': 'hair_field',
    'precedent': 'data/derived/joint-promotion-candidate-v1 role_vocabulary block',
}


def region_base(name):
    return name[:-2] if name.endswith(('.l', '.r')) else name


def lexicon_index():
    index = {}
    for field, regions in LEXICON.items():
        for region in regions:
            if region in index:
                raise ValueError('Region claimed twice: ' + region)
            index[region] = field
    return index


def area_sample(vertices, faces, spacing, rng):
    area = np.linalg.norm(np.cross(vertices[faces[:, 1]] - vertices[faces[:, 0]],
                                   vertices[faces[:, 2]] - vertices[faces[:, 0]]), axis=1)/2
    counts = np.maximum(1, np.ceil(area/(spacing*spacing)).astype(np.int64))
    index = np.repeat(np.arange(len(faces)), counts)
    u = rng.random((len(index), 2))
    root = np.sqrt(u[:, 0:1])
    weights = np.c_[1 - root, root*(1 - u[:, 1:2]), root*u[:, 1:2]]
    return np.einsum('ni,nij->nj', weights, vertices[faces[index]]), float(area.sum())


def face_areas(vertices, faces):
    return np.linalg.norm(np.cross(vertices[faces[:, 1]] - vertices[faces[:, 0]],
                                   vertices[faces[:, 2]] - vertices[faces[:, 0]]), axis=1)/2


def submesh(vertices, faces, mask):
    used = np.unique(faces[mask])
    remap = np.full(len(vertices), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    return vertices[used], remap[faces[mask]]


def assign_regions(centroids, region_samples, region_labels):
    distance, which = cKDTree(region_samples).query(centroids, workers=-1)
    return region_labels[which], distance


def shell_footprint(centroids, shell_samples, standoff):
    distance, _ = cKDTree(shell_samples).query(centroids, workers=-1)
    return distance < standoff, distance


def strand_mass_kg(count, radius_m, length_m, density_kg_m3):
    """One follicle carries at most one shaft; a circular section is the stated idealisation."""
    return float(count)*np.pi*radius_m*radius_m*length_m*density_kg_m3
