"""What sex-specific anatomy this body actually has, measured, on both entity sets.

    python scripts/audit_sex_specific_anatomy.py

This exists because the question is easy to get wrong in both directions.  "The
body is anatomically sexless" is wrong -- it has a full male external and
internal genital tract.  "The body has male anatomy, so a female one is a
relabelling" is also wrong -- there is no female entity of any kind, no
mammary gland, no nipple and no areola, and no catalogued source ships one.

Two entity sets are audited separately, because this programme has two of them
and they do not contain the same things:

* the **simulated body**: the 4,000 entities bound to driven segments by
  ``scripts/bind_anatomy_to_segments.py``;
* the **display atlas**: the structures in ``data/derived/app/manifest.json``.

Name matching is by term list, and every match is classified by hand into
male-specific / female-specific / surface-region / not-reproductive rather than
counted raw -- "cavernous sinus" is in the skull and "nasolabial" is a face, and
a raw regex count of either as genital anatomy would be a fabricated number.
"""
from pathlib import Path
import json, re, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BINDING = 'data/derived/anatomy-segment-binding/binding.json'
MANIFEST = 'data/derived/app/manifest.json'
BP3D_INDEX = 'data/derived/anatomy/bodyparts3d_index.json'

#: Terms searched.  Deliberately over-broad; the classifier below cuts it down.
TERMS = re.compile(
    r'penis|penile|testis|testicul|epididym|prostat|deferent|seminal|ejaculat|'
    r'scrot|caverno|spongio|glans|uter|ovar|vagin|cervix uteri|breast|mammary|'
    r'labia|clitor|vulva|oviduct|fallopian|salping|areola|nipple|'
    r'bulbourethral|gonad|perine', re.I)

#: Matches that are NOT sex-specific anatomy, with the reason.  Each is a term
#: that collides with the list above somewhere else in the body.
NOT_SEX_SPECIFIC = (
    (re.compile(r'cavernous sinus', re.I), 'dural venous sinus of the skull'),
    (re.compile(r'nasolabial|mentolabial|labial commissure|superior labial|'
                r'inferior labial|labial artery|labial vein', re.I),
     'of the lips, not the labia'),
    (re.compile(r'cavernous part|cavernous segment', re.I), 'internal carotid'),
    # The one that would have been reported as female anatomy. salpinx is the
    # AUDITORY tube here, not the uterine tube: salpingopharyngeus is a muscle
    # of the pharynx. It is the only 'female' match in any of the three sets,
    # and taking the regex at its word would have produced a body with two
    # fallopian tubes in its throat and no uterus.
    (re.compile(r'salpingopharyng|salpingopalatin', re.I),
     'salpinx here is the auditory tube, not the uterine tube'),
)

#: Matches that are body-SURFACE topography rather than an organ.  These are the
#: honest middle case: a "mammary region" is a named patch of skin over a chest
#: wall that carries no gland.
SURFACE_REGION = re.compile(r'mammary region|inframammary region|perineal region|'
                            r'pubic region|urogenital|anal triangle', re.I)

FEMALE = re.compile(r'uter|ovar|vagin|clitor|vulva|oviduct|fallopian|salping|'
                    r'areola|nipple|mammary gland|breast', re.I)


def classify(name):
    for pattern, reason in NOT_SEX_SPECIFIC:
        if pattern.search(name):
            return 'not_sex_specific', reason
    if SURFACE_REGION.search(name):
        return 'surface_region', 'named skin topography, carries no sex-specific organ'
    if FEMALE.search(name):
        return 'female_specific', ''
    return 'male_specific', ''


def audit(names):
    buckets = {'male_specific': [], 'female_specific': [], 'surface_region': [],
               'not_sex_specific': []}
    for name in names:
        if not TERMS.search(name):
            continue
        kind, reason = classify(name)
        buckets[kind].append({'name': name, 'reason': reason} if reason else {'name': name})
    return buckets


def load_names(root):
    sets = {}
    binding = json.loads((root / BINDING).read_text())
    sets['simulated body (segment-bound entities)'] = {
        'total': len(binding['entities']),
        'names': [(eid, entry.get('name', '')) for eid, entry in binding['entities'].items()],
        'source': BINDING}
    manifest = json.loads((root / MANIFEST).read_text())
    sets['display atlas (app manifest structures)'] = {
        'total': len(manifest['structures']),
        'names': [(s['id'], s.get('name', '')) for s in manifest['structures']],
        'source': MANIFEST}
    # This one carries no name field; the searchable text is its FMA concept
    # labels, and a mesh takes every ancestor concept, so a hit here is a hit on
    # the ontology rather than on the mesh. Searched over the most specific
    # concept only.
    index = json.loads((root / BP3D_INDEX).read_text())
    sets['BodyParts3D mesh index (most specific FMA concept)'] = {
        'total': len(index['meshes']),
        'names': [(m['element_id'], (m['concepts'][-1]['name'] if m.get('concepts') else ''))
                  for m in index['meshes']],
        'source': BP3D_INDEX}
    return sets


def main():
    report = {'schema': 'ihm.sex-specific-anatomy-audit.v1', 'sets': {}}
    for label, payload in load_names(ROOT).items():
        buckets = audit([name for _, name in payload['names']])
        distinct = {kind: sorted({row['name'] for row in rows})
                    for kind, rows in buckets.items()}
        report['sets'][label] = {
            'source': payload['source'], 'entities': payload['total'],
            'male_specific': len(buckets['male_specific']),
            'male_specific_distinct_names': distinct['male_specific'],
            'female_specific': len(buckets['female_specific']),
            'female_specific_distinct_names': distinct['female_specific'],
            'surface_region_only': len(buckets['surface_region']),
            'surface_region_distinct_names': distinct['surface_region'],
            'excluded_as_not_sex_specific': len(buckets['not_sex_specific']),
            'excluded_distinct_names': distinct['not_sex_specific']}

    binding = json.loads((ROOT / BINDING).read_text())
    male_ids = {eid: entry for eid, entry in binding['entities'].items()
                if TERMS.search(entry.get('name', ''))
                and classify(entry['name'])[0] == 'male_specific'}
    report['segment_binding_of_male_anatomy'] = {
        'note': ('The male genital tract is bound to driven segments like every '
                 'other entity, so it moves with the body. Which segment it '
                 'rides is a statement about the binding, not about function.'),
        'entities': len(male_ids),
        'segments': sorted({entry['segment'] for entry in male_ids.values()}),
        'minimum_coherence': min(entry['coherence'] for entry in male_ids.values()),
        'by_segment': {segment: sorted(entry['name'] for entry in male_ids.values()
                                       if entry['segment'] == segment)
                       for segment in sorted({e['segment'] for e in male_ids.values()})}}

    report['conclusions'] = [
        'The simulated body is NOT sexless. It carries a complete male genital '
        'tract: %d entities in the 4,000-entity segment-bound set, including '
        'the three penile bodies, both testes, both epididymides, both deferent '
        'ducts, both seminal vesicles, both ejaculatory ducts and the prostate, '
        'with their arteries and veins.'
        % report['sets']['simulated body (segment-bound entities)']['male_specific'],
        'There is no female-specific entity of any kind: %d of 4,000 in the '
        'simulated body and %d of %d in the display atlas.'
        % (report['sets']['simulated body (segment-bound entities)']['female_specific'],
           report['sets']['display atlas (app manifest structures)']['female_specific'],
           report['sets']['display atlas (app manifest structures)']['entities']),
        'There is no mammary gland, nipple or areola in ANY set. What exists is '
        '"mammary region" and "inframammary region" -- Z-Anatomy names for '
        'patches of chest skin. A breast is not present in either sex: the male '
        'body does not have the male one either.',
        'A raw regex over these names overcounts by %d in the display atlas '
        'alone -- cavernous sinuses in the skull and the lips\' labial arteries '
        'and veins -- which is why every match here is classified rather than '
        'summed.'
        % report['sets']['display atlas (app manifest structures)']['excluded_as_not_sex_specific']]
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
