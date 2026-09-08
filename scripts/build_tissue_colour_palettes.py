"""Named colour palettes for the materialization, emitted as data.

CANDIDATE. Reads data/derived/app/manifest.json and the structure-provenance candidate;
writes only under --output. Nothing canonical is modified. Every input is hashed and every
hash is recorded. The manifest does not list itself.

  .venv/bin/python scripts/build_tissue_colour_palettes.py --self-test
  .venv/bin/python scripts/build_tissue_colour_palettes.py \
      --output data/derived/tissue-colour-palette-candidate-v1

What this produces
  1 A tissue-class assignment for all 7390 display structures. One rule fires per structure and
    the rule text is recorded, so the assignment is auditable rather than hand-listed.
  2 Palettes, each a {class -> colour} table plus {structure_id -> colour} overrides, resolved
    into a flat {structure_id -> hex} map so a viewer applies one lookup.
      didactic     the palette the app ships today, read back out of the manifest, unaltered.
      realistic    fresh in-vivo appearance. Not a preserved cadaver and not an atlas plate.
      evidence     the tier of the realistic palette's colour for that structure, so a viewer can
                   see at a glance which structures are measured and which are asserted.
  3 A mucosal mapping: which canonical entities carry mucous membrane or a transitional zone,
    which mucosal targets have no entity in this body, and why.
  4 A source list. Colours converted from published CIE L*a*b* record the paper, the cohort, the
    instrument, the illuminant and standard observer, and the conversion inputs. Colours with no
    measurement behind them are tiered transferred or synthesized and say what they rest on.

Tiers use the vocabulary already in ihm.structure-provenance.v1, applied to colour rather than
geometry: measured (colorimetry of this tissue in humans), transferred (colorimetry of a
different tissue, species or state, mapped across with the mapping stated), derived (computed
from a recorded measured value, e.g. a lightness offset), synthesized (constructed from a
descriptive source with no colorimetry). "assumed" in the ask maps to synthesized.

Skin is a parameter, not a default. The realistic palette carries a declared skin-tone option
set keyed on the Individual Typology Angle, and the palette's own skin entry is the option the
caller selected (default: none applied, the neutral mid class), never a hidden constant.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

MANIFEST = ROOT/'data/derived/app/manifest.json'
PROVENANCE_INDEX = ROOT/'data/derived/structure-provenance-candidate-v1/index.json'
PROVENANCE_MANIFEST = ROOT/'data/derived/structure-provenance-candidate-v1/manifest.json'
DEFAULT_OUT = ROOT/'data/derived/tissue-colour-palette-candidate-v1'

SCHEMA_INDEX = 'ihm.tissue-colour-palette-index.v1'
SCHEMA_PALETTE = 'ihm.tissue-colour-palette.v1'
SCHEMA_MANIFEST = 'ihm.tissue-colour-palette-manifest.v1'
PROVENANCE_SCHEMA = 'ihm.structure-provenance.v1'

TIERS = {
    'measured': 'CIE L*a*b* colorimetry of this tissue, in humans, in the state the palette claims to represent.',
    'transferred': 'Colorimetry of a different tissue, species, site or state, mapped onto this tissue. The mapping is stated.',
    'derived': 'Computed from a recorded measured value in this table by a stated operation, with no measurement of its own.',
    'synthesized': 'Constructed from a descriptive (non-colorimetric) source. No measurement of this colour exists here.'}


# ---------------------------------------------------------------- colour space

# CIE 1931 2-degree and 1964 10-degree tristimulus values of the standard illuminants, x100.
WHITE = {'D65_2': (95.047, 100.0, 108.883), 'D50_2': (96.422, 100.0, 82.521),
         'C_2': (98.074, 100.0, 118.232), 'A_2': (109.850, 100.0, 35.585),
         'D65_10': (94.811, 100.0, 107.304), 'D50_10': (96.720, 100.0, 81.427),
         'C_10': (97.285, 100.0, 116.145)}
BRADFORD = np.array([[0.8951, 0.2664, -0.1614], [-0.7502, 1.7135, 0.0367], [0.0389, -0.0685, 1.0296]])
# IEC 61966-2-1 sRGB, D65 white, linear XYZ -> linear RGB.
M_XYZ2RGB = np.array([[3.2404542, -1.5371385, -0.4985314],
                      [-0.9692660, 1.8760108, 0.0415560],
                      [0.0556434, -0.2040259, 1.0572252]])


def lab_to_xyz(lab, white):
    L, a, b = lab
    fy = (L+16)/116
    fx, fz = fy+a/500, fy-b/200
    f = lambda t: t**3 if t**3 > 216/24389 else (116*t-16)/(24389/27)
    return np.array([f(fx)*white[0], f(fy)*white[1], f(fz)*white[2]])/100.0


def bradford_adapt(xyz, src, dst):
    s = BRADFORD@(np.array(src)/100.0)
    d = BRADFORD@(np.array(dst)/100.0)
    return np.linalg.inv(BRADFORD)@np.diag(d/s)@BRADFORD@xyz


def lab_to_srgb(lab, illuminant='D65_2'):
    """CIE L*a*b* under the stated illuminant/observer -> 8-bit sRGB hex.

    sRGB is defined on D65/2 degree. A value measured under any other white point is chromatically
    adapted with the Bradford transform before encoding, and the adaptation is reported so a
    reader can tell an unadapted conversion from an adapted one.
    """
    if illuminant not in WHITE:
        raise ValueError('Unknown illuminant/observer: '+illuminant)
    xyz = lab_to_xyz(lab, WHITE[illuminant])
    adapted = illuminant != 'D65_2'
    if adapted:
        xyz = bradford_adapt(xyz, WHITE[illuminant], WHITE['D65_2'])
    linear = M_XYZ2RGB@xyz
    out_of_gamut = bool(np.any(linear < -1e-6) or np.any(linear > 1+1e-6))
    excursion = float(max(0.0, float(np.max(np.abs(np.clip(linear, 0, 1)-linear)))))
    linear = np.clip(linear, 0, 1)
    encoded = np.where(linear <= 0.0031308, 12.92*linear, 1.055*linear**(1/2.4)-0.055)
    eight = np.clip(np.round(encoded*255), 0, 255).astype(int)
    return ('#%02x%02x%02x' % tuple(int(v) for v in eight),
            {'illuminant_observer': illuminant, 'chromatic_adaptation': 'bradford_to_D65_2' if adapted else 'none',
             'transfer': 'IEC 61966-2-1 sRGB', 'out_of_srgb_gamut': out_of_gamut,
             'gamut_clip_linear_excursion': round(excursion, 6)})


def hex_to_lab(hex_colour):
    """Round-trip helper used by the self-test only; sRGB D65/2 -> L*a*b*."""
    v = np.array([int(hex_colour[i:i+2], 16)/255 for i in (1, 3, 5)])
    linear = np.where(v <= 0.04045, v/12.92, ((v+0.055)/1.055)**2.4)
    xyz = np.linalg.inv(M_XYZ2RGB)@linear*100.0
    w = np.array(WHITE['D65_2'])
    r = xyz/w
    f = np.where(r > 216/24389, np.cbrt(r), (24389/27*r+16)/116)
    return (float(116*f[1]-16), float(500*(f[0]-f[1])), float(200*(f[1]-f[2])))


# ---------------------------------------------------------------- tissue classes

def normalized(name):
    n = name.lower().strip()
    n = re.sub(r'\.(l|r)$', '', n)
    n = re.sub(r'^\((.*)\)$', r'\1', n)
    return n


def _w(*words):
    return re.compile(r'(?<![a-z])(?:'+'|'.join(words)+r')(?![a-z])')


_WORD_CACHE = {}


def has(text, *words):
    key = words
    if key not in _WORD_CACHE:
        _WORD_CACHE[key] = _w(*words)
    return _WORD_CACHE[key].search(text) is not None


def tissue_class(name, system):
    """One class per structure, first rule wins; returns (class, rule text).

    The mucosal and transitional set is tested before anything else, because that is exactly the
    set that would otherwise inherit a skin or system colour, which is the defect being fixed.
    """
    n = normalized(name)
    if n in ('lip', 'tubercle of upper lip', 'labial commissure', 'angle of mouth'):
        return 'mucosa_lip_vermilion', 'named lip / vermilion border entity'
    if has(n, 'gingiva'):
        return 'mucosa_gingival', 'name contains gingiva'
    if n == 'tongue':
        return 'mucosa_lingual', 'tongue body, dorsal mucosa dominates the visible surface'
    if n in ('soft palate', 'uvula of palate'):
        return 'mucosa_palatal', 'named palatal mucosa entity'
    if n == 'mucosa of nasal cavity':
        return 'mucosa_nasal', 'explicitly named nasal mucosa'
    if n == 'mucosa of stomach':
        return 'mucosa_gastric', 'explicitly named gastric mucosa'
    if n == 'glans penis':
        return 'mucosa_glans', 'glans surface is a transitional mucosal zone, not keratinized skin'
    if n == 'anal region':
        return 'mucosa_anal_margin', 'topographic region containing the anal margin and perianal skin'
    if has(n, 'tonsil'):
        return 'lymphoid_tonsil', 'palatine tonsil, mucosa-covered lymphoid tissue'
    if n in ('nasopharynx', 'oropharynx', 'laryngopharynx', 'pharynx'):
        return 'mucosa_pharyngeal', 'named pharyngeal segment, mucosal luminal surface'
    if n == 'spleen':
        return 'spleen', 'spleen proper, distinguished from splenic vessels and splenic nodes'
    if has(n, 'thymus'):
        return 'thymus', 'name contains thymus'
    if system == 'arterial':
        return 'artery_wall', 'arterial system membership; the visible surface is vessel wall, not lumen blood'
    if system == 'venous':
        return 'vein_wall', 'venous system membership; the visible surface is vessel wall, not lumen blood'
    if system == 'lymphatic':
        return ('lymph_node', 'lymphatic system, name contains node') if has(n, 'node', 'nodes') \
            else ('lymph_vessel', 'lymphatic system, not a node')
    if has(n, 'cornea'):
        return 'cornea', 'name contains cornea'
    if has(n, 'sclera'):
        return 'sclera', 'name contains sclera'
    if has(n, 'iris'):
        return 'iris', 'name contains iris'
    if has(n, 'lens'):
        return 'lens', 'name contains lens'
    if has(n, 'retina') or n.startswith('optic part of'):
        return 'retina', 'named retinal entity'
    if has(n, 'vitreous') or 'anterior chamber' in n:
        return 'ocular_humour', 'ocular fluid compartment'
    if has(n, 'choroid') and 'plexus' not in n:
        return 'choroid', 'ocular choroid'
    if 'corona ciliaris' in n:
        return 'ciliary_body', 'ciliary body'
    if has(n, 'lacrimal') and system == 'sensory':
        return 'gland_serous', 'lacrimal gland and its drainage apparatus'
    if n == 'external ear':
        return 'skin', 'auricle, presented as its skin envelope over auricular cartilage'
    if system in ('integumentary', 'hair'):
        if 'nail plate' in n or has(n, 'perionyx'):
            return 'nail_plate', 'nail plate or perionychium'
        if has(n, 'hair', 'hairs', 'eyelashes', 'eyebrow'):
            return 'hair', 'hair-bearing entity'
        if n == 'skin' or has(n, 'philtrum') or n in ('urogenital region', 'external ear'):
            return 'skin', 'whole-body skin, or a skin region with no distinct pigmentation claim'
        if has(n, 'palm', 'sole') or 'palmar surface' in n or 'plantar surface' in n:
            return 'skin_palmoplantar', 'palmar or plantar surface; markedly lighter than the rest of the integument'
        return 'skin_region', 'topographic surface region of the integument'
    if system == 'nervous':
        if has(n, 'nerve', 'ganglion', 'trunk', 'ramus', 'rami', 'root', 'roots') or (
                has(n, 'plexus') and 'choroid plexus' not in n):
            return 'peripheral_nerve', 'peripheral nerve, ganglion or plexus'
        if 'choroid plexus' in n:
            return 'choroid_plexus', 'choroid plexus'
        if 'ventricle' in n or has(n, 'aqueduct', 'cistern') or 'central canal' in n or 'subarachnoid' in n:
            return 'csf_space', 'cerebrospinal fluid space'
        if has(n, 'capsule', 'callosum', 'commissure', 'fornix', 'tract', 'stria', 'peduncle',
               'radiation', 'fasciculus', 'lemniscus', 'chiasm', 'brachium'):
            return 'brain_white', 'named white-matter tract or commissure'
        if 'spinal cord' in n or has(n, 'medulla', 'cord'):
            return 'brain_white', 'cord or medulla, surface is predominantly white matter'
        return 'brain_grey', 'central nervous entity that is not a named tract, fluid space or nerve'
    if has(n, 'tendon', 'tendons', 'aponeurosis', 'aponeuroses') or 'tendinous' in n:
        return 'tendon', 'tendon or aponeurosis'
    if system == 'muscular':
        return 'skeletal_muscle', 'muscular system membership'
    if system == 'cardiac':
        if 'cavity of' in n:
            return 'blood_cavity', 'cardiac chamber cavity, filled with blood not tissue'
        if has(n, 'cusp', 'leaflet') or 'valve' in n:
            return 'valve_leaflet', 'cardiac valve cusp or leaflet'
        return 'cardiac_muscle', 'cardiac system membership'
    if has(n, 'cartilage', 'cartilages') or 'intervertebral disc' in n or has(n, 'epiglottis'):
        return 'cartilage_hyaline', 'hyaline or elastic cartilage'
    if has(n, 'meniscus', 'menisci', 'labrum') or 'articular disc' in n or 'anulus fibrosus' in n:
        return 'cartilage_fibro', 'fibrocartilage'
    if system == 'skeletal':
        if has(n, 'tooth', 'teeth', 'incisor', 'canine', 'premolar', 'molar', 'dentine', 'enamel'):
            return 'tooth', 'tooth; the visible surface is enamel'
        return 'cortical_bone', 'skeletal system membership; the visible surface is cortical bone'
    if has(n, 'ligament', 'ligaments', 'retinaculum', 'retinacula'):
        return 'ligament', 'ligament or retinaculum'
    if has(n, 'bursa', 'bursae'):
        return 'serous_membrane', 'synovial bursa'
    if has(n, 'liver') or re.search(r'segment of liver|hepatovenous segment', n):
        return 'liver', 'liver, a lobe or a hepatovenous segment of it'
    if 'biliary' in n or has(n, 'bile') or 'hepatic duct' in n or 'cystic duct' in n:
        return 'bile_duct', 'biliary duct or duct tree'
    if has(n, 'gallbladder'):
        return 'gallbladder', 'gallbladder'
    if has(n, 'pancreas') or 'pancreatic' in n:
        return 'pancreas', 'pancreas or pancreatic duct'
    if has(n, 'parotid', 'submandibular', 'sublingual'):
        return 'gland_salivary', 'major salivary gland or its duct'
    if has(n, 'stomach'):
        return 'stomach_serosa', 'stomach; the exposed surface is serosa, not mucosa'
    if has(n, 'oesophagus', 'esophagus'):
        return 'oesophagus', 'oesophagus'
    if has(n, 'duodenum', 'jejunum', 'ileum') or 'small intestine' in n or 'ileocecal' in n:
        return 'small_intestine', 'small bowel segment'
    if has(n, 'colon', 'caecum', 'cecum', 'appendix', 'rectum') or 'taenia' in n or 'meso-appendix' in n:
        return 'large_intestine', 'large bowel segment or taenia'
    if has(n, 'omentum', 'mesentery', 'mesocolon', 'mesoappendix'):
        return 'mesentery_omentum', 'peritoneal fold; visually dominated by its adipose content'
    if has(n, 'bronchus', 'bronchi', 'trachea') or 'bronchial tree' in n:
        return 'airway_wall', 'conducting airway; tested before lung so a named segmental bronchus is not read as parenchyma'
    if has(n, 'lung'):
        return 'lung', 'lung lobe'
    if has(n, 'pleura', 'pericardium', 'peritoneum'):
        return 'serous_membrane', 'serous membrane'
    if has(n, 'kidney') or 'renal pelvis' in n:
        return 'kidney', 'kidney or renal pelvis'
    if has(n, 'ureter', 'urethra', 'bladder'):
        return 'urothelium_wall', 'urothelium-lined conduit or reservoir'
    if has(n, 'testis', 'testes'):
        return 'testis', 'testis; the visible surface is tunica albuginea'
    if has(n, 'epididymis'):
        return 'epididymis', 'epididymis'
    if 'deferent duct' in n or 'ductus deferens' in n or 'ejaculatory duct' in n:
        return 'duct_muscular', 'thick-walled muscular genital duct'
    if 'seminal' in n:
        return 'seminal_vesicle', 'seminal vesicle'
    if has(n, 'prostate'):
        return 'prostate', 'prostate'
    if 'corpus cavernosum' in n or 'corpus spongiosum' in n:
        return 'erectile_tissue', 'erectile body; the visible surface is tunica albuginea over cavernous tissue'
    if has(n, 'thyroid'):
        return 'gland_thyroid', 'thyroid gland'
    if 'parathyroid' in n:
        return 'gland_parathyroid', 'parathyroid gland'
    if has(n, 'adrenal') or 'suprarenal' in n:
        return 'gland_adrenal', 'adrenal gland'
    if 'hypophysis' in n or 'pituitary' in n:
        return 'gland_pituitary', 'pituitary gland or one of its lobes'
    if has(n, 'pineal'):
        return 'gland_pineal', 'pineal gland'
    if 'nasolacrimal duct' in n:
        return 'mucosa_nasal', 'nasolacrimal duct, lined by respiratory mucosa continuous with the nose'
    if has(n, 'fascia', 'fasciae', 'septum', 'membrane', 'sheath', 'raphe', 'tract'):
        return 'fascia', 'fascia, intermuscular septum or fibrous sheath'
    if system == 'connective':
        if has(n, 'fat', 'fatty', 'pad') or 'adipose' in n:
            return 'adipose', 'named fat pad'
        if has(n, 'capsule'):
            return 'fascia', 'fibrous articular capsule'
        return 'fascia', 'connective system membership with no more specific rule'
    if system in ('bioelectric', 'microvascular'):
        return 'field_overlay', 'computed field overlay, not a tissue surface'
    return 'unclassified', 'no rule matched'


# ---------------------------------------------------------------- sources

# Every identifier below was re-resolved during this build. See sources.json for the resolution
# receipts: a source with resolution_verified false must not carry a measured tier.
SOURCES = {}
COLOURS = {}
SKIN_TONES = {}
ENTITY_OVERRIDES = {}
UNMAPPED = []


def load(module_path):
    """Colour tables live in a sibling data module so this file stays reviewable."""
    global SOURCES, COLOURS, SKIN_TONES, ENTITY_OVERRIDES, UNMAPPED
    import tissue_colour_tables as t
    SOURCES, COLOURS = t.SOURCES, t.COLOURS
    SKIN_TONES, ENTITY_OVERRIDES, UNMAPPED = t.SKIN_TONES, t.ENTITY_OVERRIDES, t.UNMAPPED
    return t


# ---------------------------------------------------------------- build

def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(*args):
    try:
        out = subprocess.run(('git',)+args, cwd=ROOT, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def builder_record(relative):
    commit = git('log', '-1', '--format=%H', '--', relative)
    dirty = git('status', '--porcelain', '--', relative)
    return {'script': relative, 'script_sha256': sha256_file(ROOT/relative), 'commit': commit,
            'commit_covers_working_tree': (dirty == '' or dirty is None) and commit is not None,
            'uncommitted_changes': bool(dirty)}


def resolve_colour(entry):
    """A colour table entry -> hex plus the conversion receipt."""
    if 'lab' in entry:
        hexed, receipt = lab_to_srgb(tuple(entry['lab']), entry.get('illuminant_observer', 'D65_2'))
        return hexed, {'basis': 'converted_from_cielab', 'lab': list(entry['lab']), **receipt}
    if 'srgb_hex' in entry:
        return entry['srgb_hex'], {'basis': 'stated_directly_in_srgb',
                                   'lab_under_d65_2': [round(v, 3) for v in hex_to_lab(entry['srgb_hex'])],
                                   'illuminant_observer': 'D65_2', 'chromatic_adaptation': 'none',
                                   'transfer': 'IEC 61966-2-1 sRGB', 'out_of_srgb_gamut': False,
                                   'gamut_clip_linear_excursion': 0.0}
    raise ValueError('Colour entry carries neither lab nor srgb_hex: '+repr(entry))


def build_classes(structures):
    rows, rules = [], {}
    for s in structures:
        cls, rule = tissue_class(s['name'], s['system'])
        rules.setdefault(cls, set()).add(rule)
        rows.append({'structure_id': s['id'], 'name': s['name'], 'model_id': s['model_id'],
                     'system': s['system'], 'tissue_class': cls, 'rule': rule})
    return rows, {k: sorted(v) for k, v in rules.items()}


def build_palettes(structures, rows, tables, skin_tone):
    by_id = {r['structure_id']: r for r in rows}
    didactic_system = {s['id']: s['color'] for s in structures}
    systems = {s['id']: s for s in json.loads(MANIFEST.read_bytes())['systems']}

    palettes = {}

    # 1 didactic: exactly what the manifest already carries. Read back, never recomputed.
    palettes['didactic'] = {
        'id': 'didactic', 'label': 'Didactic',
        'description': 'The palette the app ships today. One colour per anatomical system, '
                       'chosen so systems separate at a glance. It is textbook convention and '
                       'does not claim to resemble tissue.',
        'represents': 'system membership, not appearance',
        'is_app_default': True,
        'assignment_basis': 'per system, read verbatim out of data/derived/app/manifest.json',
        'roles': {sid: {'label': s['name'], 'srgb_hex': s['color'], 'tier': 'synthesized',
                        'tier_basis': 'didactic convention; no colorimetric claim is made',
                        'sources': []}
                  for sid, s in systems.items()},
        'role_key': 'system',
        'colours': dict(didactic_system)}

    # 2 realistic
    roles, tier_of = {}, {}
    for cls, entry in tables.COLOURS.items():
        chosen = entry
        if entry.get('takes_skin_colour'):
            # skin and its topographic regions carry the selected option, never a constant.
            chosen = {**entry, **tables.SKIN_TONES[skin_tone]['colour']}
        hexed, receipt = resolve_colour(chosen)
        roles[cls] = {'label': entry['label'], 'srgb_hex': hexed, 'tier': chosen['tier'],
                      'tier_basis': chosen['tier_basis'], 'note': entry.get('note'),
                      'sources': chosen.get('sources', []), 'conversion': receipt,
                      'measured_state': entry.get('state')}
        tier_of[cls] = chosen['tier']
    colours, tiers = {}, {}
    for s in structures:
        cls = by_id[s['id']]['tissue_class']
        override = tables.ENTITY_OVERRIDES.get(s['id']) or tables.ENTITY_OVERRIDES.get(normalized(s['name']))
        if override is not None:
            hexed, receipt = resolve_colour(override)
            colours[s['id']] = hexed
            tiers[s['id']] = override['tier']
        else:
            colours[s['id']] = roles[cls]['srgb_hex']
            tiers[s['id']] = tier_of[cls]
    palettes['realistic'] = {
        'id': 'realistic', 'label': 'Realistic (fresh, in vivo)',
        'description': 'Fresh in-vivo appearance under neutral white light: the colour of living, '
                       'perfused tissue at operation, not a formalin-fixed cadaver and not an atlas '
                       'plate. Values converted from published CIE L*a*b* where that colorimetry '
                       'exists; everything else is tiered transferred or synthesized and says so.',
        'represents': 'fresh, perfused, in-vivo tissue under neutral white illumination',
        'not_represented': ['embalmed or formalin-fixed cadaveric tissue (markedly desaturated and browner)',
                            'plastinated specimens', 'atlas plate convention', 'histological stain colour'],
        'is_app_default': False,
        'assignment_basis': 'per tissue class, with per-entity overrides where the class is too coarse',
        'roles': roles, 'role_key': 'tissue_class',
        'skin_tone': {'selected': skin_tone, **tables.SKIN_TONES[skin_tone]},
        'colours': colours, 'tiers': tiers}

    # 3 evidence: the realistic palette's tier, made visible
    tier_colour = {'measured': '#2f6f4f', 'transferred': '#4f7fa8', 'derived': '#b08637', 'synthesized': '#a45a5a'}
    palettes['evidence'] = {
        'id': 'evidence', 'label': 'Evidence tier',
        'description': 'Not an appearance palette. Each structure is drawn in the tier of its '
                       'realistic-palette colour, so a viewer sees at a glance which structures rest '
                       'on colorimetry of that tissue and which are asserted. Included because a '
                       'realistic palette that cannot be audited on screen invites its weakest '
                       'entries to be read as strongly as its strongest.',
        'represents': 'the evidence tier of the realistic palette, per structure',
        'is_app_default': False,
        'assignment_basis': 'per structure, from the realistic palette tier',
        'roles': {tier: {'label': tier, 'srgb_hex': colour, 'tier': 'synthesized',
                         'tier_basis': 'legend colour, chosen for separability; no colorimetric claim',
                         'sources': [], 'note': TIERS[tier]}
                  for tier, colour in tier_colour.items()},
        'role_key': 'tier',
        'colours': {sid: tier_colour[t] for sid, t in tiers.items()}}
    return palettes


def build(out_dir, skin_tone='not_applied'):
    tables = load(None)
    out = Path(out_dir)
    (out/'palettes').mkdir(parents=True, exist_ok=True)

    manifest_data = json.loads(MANIFEST.read_bytes())
    structures = manifest_data['structures']
    rows, rules = build_classes(structures)
    palettes = build_palettes(structures, rows, tables, skin_tone)

    # sources.json: every identifier with its re-resolution receipt.
    sources = {'schema': PROVENANCE_SCHEMA, 'kind': 'colour_source_list',
               'identifier_policy': 'Every PMID and DOI here was re-resolved against PubMed or '
                                    'Crossref during research for this build and the resolved title '
                                    'checked against the one claimed. resolution_verified false means '
                                    'the identifier did not resolve to the claimed work; such a source '
                                    'may not carry a measured tier.',
               'tiers': TIERS, 'sources': tables.SOURCES}

    classes = {'schema': 'ihm.tissue-class-assignment.v1',
               'basis': 'first matching rule over (normalized name, system) from '
                        'scripts/build_tissue_colour_palettes.py:tissue_class',
               'class_rules': rules,
               'counts': {c: sum(1 for r in rows if r['tissue_class'] == c) for c in sorted(rules)},
               'structures': rows}

    mucosal_classes = sorted(c for c in rules if c.startswith('mucosa_') or c in ('lymphoid_tonsil',))
    mucosa = {'schema': 'ihm.mucosal-structure-mapping.v1',
              'intent': 'Which canonical entities carry mucous membrane or a mucocutaneous '
                        'transitional zone, so they are coloured from the mucosal set rather than '
                        'inheriting a skin or system colour.',
              'mapped': {c: [{'structure_id': r['structure_id'], 'name': r['name'], 'model_id': r['model_id']}
                             for r in rows if r['tissue_class'] == c]
                         for c in mucosal_classes},
              'mapped_structure_count': sum(1 for r in rows if r['tissue_class'] in mucosal_classes),
              'unmapped_targets': tables.UNMAPPED}

    index = {'schema': SCHEMA_INDEX,
             'default_palette': 'didactic',
             'default_policy': 'This build does not change the palette the app uses. didactic stays '
                               'the default; a caller selects another palette explicitly.',
             'skin_tone_options': tables.SKIN_TONES,
             'skin_tone_policy': 'Constitutive skin colour is a declared parameter with a stated '
                                 'basis, never a hard-coded default. The shipped realistic palette '
                                 'is built with the neutral option; a caller may request any option '
                                 'in this set.',
             'tiers': TIERS,
             'palettes': [{'id': p['id'], 'label': p['label'], 'description': p['description'],
                           'represents': p['represents'], 'is_app_default': p['is_app_default'],
                           'role_key': p['role_key'], 'role_count': len(p['roles']),
                           'structure_count': len(p['colours']),
                           'tier_counts': ({} if 'tiers' not in p else
                                           {t: sum(1 for v in p['tiers'].values() if v == t)
                                            for t in sorted(set(p['tiers'].values()))}),
                           'path': f'palettes/{p["id"]}.json'}
                          for p in palettes.values()]}

    artifacts = {}
    def write(relative, payload):
        path = out/relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(json.dumps(payload, indent=1, sort_keys=False).encode()+b'\n')
        artifacts[relative] = sha256_file(path)

    for p in palettes.values():
        write(f'palettes/{p["id"]}.json', {'schema': SCHEMA_PALETTE, **p})
    write('index.json', index)
    write('sources.json', sources)
    write('classes.json', classes)
    write('mucosa.json', mucosa)

    # Only files this build actually reads are gated inputs. The structure-provenance candidate is
    # referenced for its tier vocabulary and hashed separately, so a reader does not fail on its absence.
    inputs = {str(MANIFEST.relative_to(ROOT)): sha256_file(MANIFEST)}
    referenced = {str(extra.relative_to(ROOT)): sha256_file(extra)
                  for extra in (PROVENANCE_INDEX, PROVENANCE_MANIFEST) if extra.is_file()}
    tables_path = Path(tables.__file__).resolve()

    manifest = {'schema': SCHEMA_MANIFEST,
                'python': sys.version,
                'builder': builder_record(str(Path(__file__).resolve().relative_to(ROOT))),
                'colour_tables': builder_record(str(tables_path.relative_to(ROOT))),
                'inputs_sha256': inputs,
                'referenced_sha256': referenced,
                'referenced_note': 'Hashed for the record, not gated: this build reads the tier vocabulary of '
                                   'ihm.structure-provenance.v1 but no bytes from these files.',
                'outputs_sha256': artifacts,
                'manifest_not_listed_in_outputs': 'manifest.json' not in artifacts,
                'structure_count': len(structures),
                'tissue_class_count': len(rules),
                'palette_ids': [p['id'] for p in palettes.values()],
                'default_palette_unchanged': index['default_palette'] == 'didactic',
                'skin_tone_selected': skin_tone,
                'source_count': len(tables.SOURCES),
                'verified_source_count': sum(1 for s in tables.SOURCES.values() if s.get('resolution_verified')),
                'tier_counts': {t: sum(1 for v in palettes['realistic']['tiers'].values() if v == t)
                                for t in sorted(TIERS)}}
    manifest['self_test'] = self_test(out, palettes, rows, tables, manifest)
    (out/'manifest.json').write_bytes(json.dumps(manifest, indent=1).encode()+b'\n')
    return palettes, manifest


# ---------------------------------------------------------------- self-test

def self_test(out, palettes, rows, tables, manifest):
    checks = []
    def check(name, passed, detail=''):
        checks.append({'check': name, 'passed': bool(passed), 'detail': str(detail)})

    original = json.loads(MANIFEST.read_bytes())
    by_system = {s['id']: s['color'] for s in original['systems']}
    check('didactic palette is the shipped one, unaltered',
          all(palettes['didactic']['colours'][s['id']] == s['color'] for s in original['structures'])
          and all(palettes['didactic']['roles'][k]['srgb_hex'] == v for k, v in by_system.items()),
          f'{len(by_system)} system colours and {len(original["structures"])} structure colours compared verbatim')

    check('the app default is not changed by this build',
          manifest['default_palette_unchanged'] and not palettes['realistic']['is_app_default']
          and not palettes['evidence']['is_app_default'],
          'index.default_palette=didactic; realistic and evidence are opt-in')

    check('every structure has a colour in every palette',
          all(len(p['colours']) == len(original['structures']) for p in palettes.values()),
          {p['id']: len(p['colours']) for p in palettes.values()})

    check('every colour is a 6-digit lowercase sRGB hex',
          all(re.fullmatch(r'#[0-9a-f]{6}', c) for p in palettes.values() for c in p['colours'].values()),
          'checked across all palettes')

    check('no structure is left unclassified',
          not any(r['tissue_class'] == 'unclassified' for r in rows),
          f"{sum(1 for r in rows if r['tissue_class'] == 'unclassified')} unclassified")

    check('every tissue class has a colour table entry',
          set(r['tissue_class'] for r in rows) <= set(tables.COLOURS),
          sorted(set(r['tissue_class'] for r in rows) - set(tables.COLOURS)) or 'all covered')

    check('no colour table entry is dead',
          set(tables.COLOURS) <= set(r['tissue_class'] for r in rows) | {'skin'},
          sorted(set(tables.COLOURS) - set(r['tissue_class'] for r in rows)) or 'all used')

    # Lab -> sRGB conversion against the sRGB primaries, whose L*a*b* under D65/2 is definitional.
    reference = {(53.2408, 80.0925, 67.2032): '#ff0000', (87.7347, -86.1827, 83.1793): '#00ff00',
                 (32.2970, 79.1875, -107.8602): '#0000ff', (100.0, 0.0, 0.0): '#ffffff',
                 (0.0, 0.0, 0.0): '#000000', (53.5850, 0.0, 0.0): '#808080'}
    got = {lab: lab_to_srgb(lab)[0] for lab in reference}
    check('CIE L*a*b* -> sRGB reproduces the sRGB primaries and neutral axis',
          got == reference, f'{sum(got[k] == v for k, v in reference.items())}/{len(reference)} exact')

    roundtrip = max(abs(np.array(hex_to_lab(v)) - np.array(k)).max() for k, v in reference.items())
    check('sRGB -> L*a*b* inverts the forward transform', roundtrip < 0.02,
          f'max |dL*,da*,db*| = {roundtrip:.4f} over the reference set (8-bit quantisation floor)')

    d50 = lab_to_srgb((100.0, 0.0, 0.0), 'D50_2')
    check('a non-D65 measurement is chromatically adapted before encoding',
          d50[0] == '#ffffff' and d50[1]['chromatic_adaptation'] == 'bradford_to_D65_2',
          f'D50/2 perfect diffuser -> {d50[0]} via {d50[1]["chromatic_adaptation"]}')

    measured = [c for c in tables.COLOURS.values() if c['tier'] == 'measured']
    check('every measured colour cites at least one source that resolved',
          all(c.get('sources') and all(tables.SOURCES[s]['resolution_verified'] for s in c['sources'])
              for c in measured),
          f'{len(measured)} measured entries')

    check('no source that failed re-resolution backs a measured colour',
          not any(not s.get('resolution_verified') and
                  any(sid in c.get('sources', []) for c in measured)
                  for sid, s in tables.SOURCES.items()),
          f"{sum(1 for s in tables.SOURCES.values() if not s.get('resolution_verified'))} unresolved sources present")

    check('every cited source id exists in the source list',
          all(s in tables.SOURCES for c in list(tables.COLOURS.values())+list(tables.ENTITY_OVERRIDES.values())
              for s in c.get('sources', [])),
          'checked colour table and entity overrides')

    check('every colour carries a tier from the provenance vocabulary',
          all(c['tier'] in TIERS for c in tables.COLOURS.values())
          and all(c['tier'] in TIERS for c in tables.ENTITY_OVERRIDES.values()),
          sorted(TIERS))

    check('every colour states the basis for its tier',
          all(c.get('tier_basis') for c in tables.COLOURS.values())
          and all(c.get('tier_basis') for c in tables.ENTITY_OVERRIDES.values()))

    mucosal = [r for r in rows if r['tissue_class'].startswith('mucosa_')]
    check('the mucosal set is separated from skin and from its system colour',
          mucosal and all(palettes['realistic']['colours'][r['structure_id']]
                          != palettes['realistic']['roles']['skin']['srgb_hex'] for r in mucosal)
          and all(palettes['realistic']['colours'][r['structure_id']] != by_system[r['system']]
                  for r in mucosal),
          f'{len(mucosal)} mucosal structures, none equal to skin or to their didactic system colour')

    # Mucosal colours must actually read as mucosa: pink-to-red-brown, i.e. positive a*.
    mucosal_lab = {c: hex_to_lab(v['srgb_hex']) for c, v in palettes['realistic']['roles'].items()
                   if c.startswith('mucosa_')}
    check('every mucosal colour is on the red side of neutral',
          all(lab[1] > 6 for lab in mucosal_lab.values()),
          {c: round(lab[1], 1) for c, lab in mucosal_lab.items()})

    check('skin is a parameter with a declared option set, not a constant',
          len(tables.SKIN_TONES) >= 4 and all(o.get('basis') and o.get('colour')
                                              for o in tables.SKIN_TONES.values()),
          f'{len(tables.SKIN_TONES)} declared options')

    variants = {name: build_palettes(original['structures'], rows, tables, name)['realistic']
                ['roles']['skin']['srgb_hex'] for name in tables.SKIN_TONES}
    check('selecting a different skin tone actually changes the skin colour',
          len(set(variants.values())) == len(variants), variants)

    check('non-skin colours do not move when the skin tone changes',
          all(build_palettes(original['structures'], rows, tables, name)['realistic']['roles']['liver']['srgb_hex']
              == palettes['realistic']['roles']['liver']['srgb_hex'] for name in tables.SKIN_TONES),
          'liver held fixed across every skin-tone option')

    check('the realistic palette is visibly different from the didactic one',
          sum(palettes['realistic']['colours'][s['id']] != s['color'] for s in original['structures'])
          > 0.9*len(original['structures']),
          f"{sum(palettes['realistic']['colours'][s['id']] != s['color'] for s in original['structures'])}"
          f"/{len(original['structures'])} structures recoloured")

    check('the evidence palette partitions the realistic palette by tier',
          set(palettes['evidence']['colours'].values()) ==
          {palettes['evidence']['roles'][t]['srgb_hex'] for t in set(palettes['realistic']['tiers'].values())},
          manifest['tier_counts'])

    # The three surfaces the owner named by hand, each either coloured from the mucosal set or
    # named as unmappable with a reason. This is the check that answers the original complaint.
    named = {'lips': ('mucosa_lip_vermilion', None),
             'nipples': (None, 'areola and nipple'),
             'tip of penis': ('mucosa_glans', None)}
    resolved = {}
    for what, (cls, target) in named.items():
        if cls is not None:
            structures = [r for r in rows if r['tissue_class'] == cls]
            resolved[what] = f'{len(structures)} structures coloured {palettes["realistic"]["roles"][cls]["srgb_hex"]} ({palettes["realistic"]["roles"][cls]["tier"]})'
        else:
            entry = next((u for u in tables.UNMAPPED if u['target'] == target), None)
            resolved[what] = 'no entity: '+(entry['reason'][:60]+'...' if entry else 'NOT DECLARED')
    check('lips, nipples and the tip of the penis are each answered',
          all(structures for what, (cls, _) in named.items() if cls
              for structures in [[r for r in rows if r['tissue_class'] == cls]])
          and any(u['target'] == 'areola and nipple' for u in tables.UNMAPPED),
          resolved)

    check('no mucosal colour equals the endothelium-like system colour it replaced',
          all(palettes['realistic']['colours'][r['structure_id']] != by_system[r['system']]
              for r in rows if r['tissue_class'] in ('mucosa_glans', 'mucosa_lip_vermilion')),
          'glans no longer takes the flat reproductive colour '+by_system['reproductive']
          +'; lip no longer takes '+by_system['digestive'])

    check('the shipped display manifest was not written by this build',
          sha256_file(MANIFEST) == manifest['inputs_sha256'][str(MANIFEST.relative_to(ROOT))],
          'data/derived/app/manifest.json rehashed after the build and matches the hash taken before it')

    check('unmapped mucosal targets are named with a reason',
          tables.UNMAPPED and all(u.get('target') and u.get('reason') for u in tables.UNMAPPED),
          f'{len(tables.UNMAPPED)} named')

    written = sorted(p.relative_to(out).as_posix() for p in out.rglob('*.json') if p.name != 'manifest.json')
    check('manifest lists every artifact and does not list itself',
          sorted(manifest['outputs_sha256']) == written and 'manifest.json' not in manifest['outputs_sha256'],
          f'{len(written)} artifacts')

    check('every declared output hash matches the file on disk',
          all(sha256_file(out/rel) == h for rel, h in manifest['outputs_sha256'].items()),
          f"{len(manifest['outputs_sha256'])} files rehashed")

    return {'passed': all(c['passed'] for c in checks), 'checks': checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default=str(DEFAULT_OUT))
    parser.add_argument('--skin-tone', default='not_applied')
    parser.add_argument('--self-test', action='store_true', help='fail the process on any failed check')
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    palettes, manifest = build(args.output, args.skin_tone)
    report = manifest['self_test']
    for entry in report['checks']:
        print(('PASS ' if entry['passed'] else 'FAIL ')+entry['check']+' :: '+entry['detail'])
    print(f"palettes={len(palettes)} classes={manifest['tissue_class_count']} "
          f"structures={manifest['structure_count']} sources={manifest['source_count']} "
          f"verified={manifest['verified_source_count']} tiers={manifest['tier_counts']} "
          f"passed={report['passed']}")
    if args.self_test and not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
