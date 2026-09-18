"""Authored endpoint-to-relay route hypotheses, never anatomical measurements.

Coordinates use the canonical body frame. Each entry is one representative
branch, not the extent of an entire branching nerve/plexus. Numbers are authored
landmarks; references support topology only, not these coordinates or lengths.

ANCHORED overrides the authored point for every endpoint whose label names an
organ the atlas carries and which the 18 Sep 2026 endpoint check found off that
organ (docs/BODY_PERIPHERAL.md): the authored trunk-and-viscera points sat
10-17 cm low, like the relays.  An anchored endpoint is a vertex of the named
BodyParts3D mesh chosen by the stated rule, so it is ON the structure by
construction; the authored point is kept below as the record of what it replaced
and, for `nearest_surface`, as the point whose nearest surface vertex is taken.
VISCERAL anchors are declared on the LEFT and the right route is the mirror image
of the left (`sides='mirror_left'`).  That is not anatomy -- the right vagus
reaches the posterior gastric wall and the right kidney sits lower than the
left -- it is the contract IBM-1 reads visceral routes under:
`ihm_bridge.visceral_routes` collapses each visceral trunk to one delay and
RAISES if its two sides differ.  With per-side organ anchors the vagus would be
448 mm left and 423 mm right and the brain's interoceptive loop would refuse
to load.  So the right-side record carries `endpoint_source.rule =
'mirror_of_left'`, its endpoint is NOT on the organ, and lifting that needs the
brain-side join to accept per-side lengths.  Phrenic is somatic and has no such
join; it is anchored per side.
"""
# name, relay group, representative distal endpoint, authored position (left)
SOMATIC = [
 ('sciatic_fibular','sacral','short-head biceps femoris region',[.11,-.32,-.055]),
 ('brachial_plexus','cervical','proximal arm',[.19,.43,-.02]),
 ('lumbar_plexus','lumbar','proximal anterior thigh',[.09,-.12,-.02]),
 ('sacral_plexus','sacral','proximal posterior thigh',[.09,-.20,-.06]),
 ('cervical_plexus','cervical','lateral neck',[.06,.57,.01]),
 ('phrenic','cervical','diaphragm',[.08,.23,.01]),
 ('pudendal','sacral','perineum',[.04,-.22,.01]),
 ('saphenous','lumbar','medial ankle',[.08,-.72,.02]),
 ('common_fibular','sacral','lateral knee',[.13,-.40,-.02]),
 ('superficial_radial','cervical','dorsoradial hand',[.32,-.01,-.01]),
 ('anterior_interosseous','cervical','deep anterior forearm',[.28,.13,.02]),
 ('posterior_interosseous','cervical','deep posterior forearm',[.28,.14,-.03]),
 ('palmar_digital','cervical','palmar finger',[.33,-.08,.03]),
 ('dorsal_digital','cervical','dorsal finger',[.33,-.08,-.01]),
 ('medial_cutaneous_arm','cervical','medial arm',[.18,.38,.02]),
 ('medial_cutaneous_forearm','cervical','medial forearm',[.25,.16,.02]),
 ('supraclavicular','cervical','skin over clavicle',[.13,.51,.04]),
 ('great_auricular','cervical','skin near ear',[.075,.68,.00]),
 ('lesser_occipital','cervical','lateral posterior scalp',[.06,.72,-.07]),
 ('transverse_cervical','cervical','anterior neck skin',[.045,.57,.045]),
 ('ansa_cervicalis','cervical','infrahyoid muscle region',[.025,.57,.035]),
 ('long_thoracic','cervical','serratus anterior region',[.15,.32,.02]),
 ('thoracodorsal','cervical','latissimus region',[.15,.28,-.08]),
 ('dorsal_scapular','cervical','rhomboid region',[.08,.43,-.09]),
 ('upper_subscapular','cervical','upper subscapularis region',[.15,.46,-.03]),
 ('lower_subscapular','cervical','lower subscapularis region',[.16,.39,-.03]),
 ('medial_pectoral','cervical','medial pectoral region',[.10,.40,.08]),
 ('lateral_pectoral','cervical','lateral pectoral region',[.16,.44,.07]),
 ('subcostal','thoracic','abdominal wall below rib twelve',[.15,.12,.07]),
 ('thoracoabdominal','thoracic','anterior abdominal wall',[.12,.16,.10]),
 ('genitofemoral','lumbar','groin',[.07,-.14,.045]),
 ('ilioinguinal','lumbar','inguinal region',[.08,-.12,.06]),
 ('iliohypogastric','lumbar','suprapubic wall',[.08,-.04,.08]),
 ('posterior_femoral_cutaneous','sacral','posterior thigh skin',[.10,-.35,-.07]),
 # Every route above is a ventral ramus. Without a dorsal ramus the entire
 # paravertebral skin strip and the deep back had no declared trunk, and a
 # dermatomal partition had to either leave it uninnervated or assert that the
 # back of the trunk is supplied by the intercostal nerve, which supplies the
 # front. One representative mid-thoracic posterior branch, not all 31 pairs.
 ('dorsal_ramus','thoracic','paravertebral skin and deep back',[.035,.36,-.105]),
]
VISCERAL = [
 ('vagus','solitary','representative gastric wall',[.06,.18,.03]),
 ('glossopharyngeal','solitary','posterior tongue',[.025,.64,.045]),
 ('greater_splanchnic','thoracic','foregut visceral endpoint',[.06,.18,.03]),
 ('lesser_splanchnic','thoracic','upper abdominal visceral endpoint',[.07,.12,.01]),
 ('least_splanchnic','thoracic','renal visceral endpoint',[.08,.10,-.03]),
 ('lumbar_splanchnic','lumbar','lower abdominal visceral endpoint',[.06,-.03,.03]),
 ('pelvic_splanchnic','sacral','pelvic visceral endpoint',[.04,-.16,.02]),
 ('sympathetic_chain','thoracic','representative thoracic visceral endpoint',[.06,.36,.00]),
]
# Explicit first relay (optic ends at LGN), not downstream cortical target.
SPECIAL = [
 ('optic','lateral_geniculate','retinal ganglion cell region',[.030,.715,.075],[.022,.710,.010],'pericalcarine'),
 ('cochlear','cochlear_nuclei','spiral ganglion region',[.060,.675,.000],[.040,.665,-.012],'transversetemporal'),
 ('vestibular','vestibular_nuclei','vestibular ganglion region',[.060,.680,.000],[.038,.670,-.012],'insula'),
 ('olfactory','olfactory_bulb','olfactory epithelium',[.012,.710,.065],[.012,.720,.060],'entorhinal'),
]

# name -> anchor rule on a named atlas entity.  Resolved in enrich_peripheral_routes.py.
ANCHORED = {
 'vagus': dict(sides='mirror_left', rule='stomach_wall', entity={'left': 'body-bp3d-FJ2564', 'right': 'body-bp3d-FJ2564'},
               label='gastric wall: the left vagus (anterior vagal trunk) on the anterior wall, '
                     'directly in front of the stomach centroid; right = mirror of left',
               was_off_mm='87 mm from the stomach surface (132 mm below its centroid)'),
 'greater_splanchnic': dict(sides='mirror_left', rule='stomach_wall', entity={'left': 'body-bp3d-FJ2564', 'right': 'body-bp3d-FJ2564'},
               label='foregut (gastric) wall via the coeliac plexus; the vagus\'s anterior wall point; '
                     'right = mirror of left',
               was_off_mm='shared the vagus point: 87 mm from the stomach surface'),
 'least_splanchnic': dict(sides='mirror_left', rule='kidney_hilum', entity={'left': 'body-bp3d-FJ3145', 'right': 'body-bp3d-FJ3147'},
               label='renal: the left kidney, surface vertex nearest the midline at its centroid '
                     'height (the hilum); right = mirror of left',
               was_off_mm='101 mm from the left kidney surface'),
 'pelvic_splanchnic': dict(sides='mirror_left', rule='nearest_surface', entity={'left': 'body-bp3d-FJ3149', 'right': 'body-bp3d-FJ3149'},
               label='pelvic viscera: urinary bladder wall, vertex nearest the authored point; '
                     'right = mirror of left',
               was_off_mm='171 mm below the bladder, in the upper thigh'),
 'phrenic': dict(sides='per_side', rule='diaphragm_dome', entity={'left': 'body-bp3d-FJ3131', 'right': 'body-bp3d-FJ3131'},
               label='diaphragm: highest diaphragm vertex within 10 mm of the authored lateral '
                     'offset (the dome the phrenic pierces)',
               was_off_mm='29 mm from the diaphragm surface, 72 mm below its centroid'),
}
