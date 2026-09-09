"""Authored endpoint-to-relay route hypotheses, never anatomical measurements.

Coordinates use the canonical body frame. Each entry is one representative
branch, not the extent of an entire branching nerve/plexus. Numbers are authored
landmarks; references support topology only, not these coordinates or lengths.
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
