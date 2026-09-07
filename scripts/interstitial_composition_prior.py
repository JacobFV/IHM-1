#!/usr/bin/env python3
"""Sourced composition prior for the 32.117 L of body interior that belongs to no entity.

    .venv/bin/python scripts/interstitial_composition_prior.py --self-test
    .venv/bin/python scripts/interstitial_composition_prior.py --summary
    .venv/bin/python scripts/interstitial_composition_prior.py --output data/derived/interstitial-composition-prior-v1
    .venv/bin/python scripts/interstitial_composition_prior.py --verify-identifiers   # network

The atlas accounts for 37.585 L of a measured 69.720 L interior. The remaining
32.117 L (46.08%), of which one connected component holds 30.907 L, is real
tissue defined by absence: dermis the canonical skin slab never gave a volume,
subcutaneous / visceral / intermuscular adipose, fascia, loose areolar connective
tissue and interstitial fluid. Its GEOMETRY is exact (it is the complement). This
script supplies the COMPOSITION and the mass that follows from it.

Measured here (arithmetic on this specimen's shipped geometry):
  ledger    envelope, exclusive occupancy, void, and the body mass a composed
            fill implies, at both shipped voxel partitions
  regional  the void's depth histogram and per-nearest-entity-system attribution
            filled by a stated greedy allocation

Transferred here (published values; nothing was measured on the BodyParts3D
specimen):
  sources      one record per publication with DOI/PMID, cohort, method and how
               the numbers were read; every identifier was resolved
  composition  per constituent: volume, density, mechanics, each cell tiered
               measured / transferred / derived / assumed / absent

Not established here: no fit to this specimen, no probability interpretation of
any +/-, no joint depth x system distribution (the input carries two marginals
only), and no measured density for loose connective tissue or interstitial fluid.
"""
from pathlib import Path
import argparse, hashlib, json, sys

ROOT = Path(__file__).resolve().parents[1]
UNMODELLED = ROOT / 'data/derived/unmodelled-volume/summary.json'
CANDIDATE_LEDGER = ROOT / 'data/derived/tissue-material-candidate-v1/mass_ledger.json'
CANDIDATE_MATERIALS = ROOT / 'data/derived/tissue-material-candidate-v1/materials.json'
CANDIDATE_SOURCES = ROOT / 'data/derived/tissue-material-candidate-v1/sources.json'
PROFILE = ROOT / 'data/derived/canonical/profile.json'
INPUTS = (UNMODELLED, CANDIDATE_LEDGER, CANDIDATE_MATERIALS, CANDIDATE_SOURCES, PROFILE)

# Same vocabulary as data/derived/tissue-material-candidate-v1 so the two tables merge.
TIERS = {
    'measured': 'computed from the BodyParts3D specimen geometry shipped in this repo.',
    'transferred': 'a published value measured on other subjects, cohorts or species.',
    'derived': 'arithmetic on a transferred or measured value; the operation is in the note.',
    'assumed': 'an engineering number with no source. Every one is listed in assumed_cells.',
    'absent': 'no defensible value was retrieved. The cell is deliberately empty.',
}

# --------------------------------------------------------------- sources
# `verification` says how the numbers were read. `identifier_check` records the
# live resolution performed on 2026-09-06; --verify-identifiers repeats it.
SOURCES = {
    'icrp89': {
        'citation': 'ICRP. Basic anatomical and physiological data for use in radiological protection: '
                    'reference values. ICRP Publication 89. Ann ICRP 2002;32(3-4):1-277.',
        'doi': '10.1016/S0146-6453(03)00002-2',
        'pmid': None,
        'species': 'human',
        'state': 'reference values, not a specimen',
        'method': 'Task Group synthesis over Western European and North American autopsy and imaging series',
        'n': 'reference adult male, 73 kg, 176 cm, 1.90 m2',
        'ages': 'adult 20-50 y',
        'verification': 'fulltext',
        'fulltext_read': {
            'url': 'https://radon-and-life.narod.ru/pub/ICRP_89.pdf',
            'sha256': '3a4df292e29ff8f79be3544a840ba34f159c3237814484fe7ad55764cdf899e1',
            'passages': ['Table 2.8 (organ and tissue masses)', 'Table 2.9 (height, mass, surface area)',
                         'Table 2.27 (epidermis, dermis, total skin)', 'Table 4.5 (distribution of body water)',
                         'Table 11.2 (separable dense connective tissue)',
                         'paras 100-105 (body fluids), 109 (body fat), 542-546 (adipose), 573-578 (connective tissue)'],
        },
        'identifier_check': {'resolver': 'crossref', 'status': 'resolved',
                             'returned_title': 'Basic anatomical and physiological data for use in radiological '
                                               'protection: reference values'},
        'correction': 'data/derived/tissue-material-candidate-v1/sources.json records pmid 14527029 for this '
                      'source. That PMID is Boecker BB, Radiat Prot Dosimetry 2003;105(1-4):571-4, a four-page '
                      'conference summary ABOUT ICRP 89, not ICRP 89. ICRP 89 itself is not in PubMed. This '
                      'record therefore carries no PMID.',
    },
    'icru44_nist': {
        'citation': 'Hubbell JH, Seltzer SM. Tables of X-Ray Mass Attenuation Coefficients and Mass '
                    'Energy-Absorption Coefficients, Table 2. NIST Standard Reference Database 126. Densities '
                    'are those of ICRU Report 44 (1989).',
        'doi': '10.18434/T4D01F',
        'pmid': None,
        'species': 'human',
        'state': 'reference composition, not a specimen',
        'method': 'standardised reference tissue compositions and mass densities',
        'n': 'reference table',
        'ages': 'adult reference',
        'verification': 'table (physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html, re-read 2026-09-06)',
        'identifier_check': {'resolver': 'datacite', 'status': 'resolved',
                             'returned_title': 'Tables of X-Ray Mass Attenuation Coefficients and Mass '
                                               'Energy-Absorption Coefficients, NIST Standard Reference '
                                               'Database 126',
                             'note': 'a DataCite DOI; api.crossref.org returns 404 for it, which is expected '
                                     'and is not a bad identifier'},
        'read_here': 'adipose tissue 0.950, whole blood 0.950 -> 1.060, skeletal muscle 1.050, soft tissue '
                     '1.060, liquid water 1.000 g/cm3. The table has NO skin and NO connective tissue entry.',
    },
    'shen2003': {
        'citation': 'Shen W, Wang Z, Punyanita M, Lei J, Sinav A, Kral JG, Imielinska C, Ross R, Heymsfield SB. '
                    'Adipose tissue quantification by imaging methods: a proposed classification. '
                    'Obes Res 2003;11(1):5-16.',
        'doi': '10.1038/oby.2003.3', 'pmid': '12529479',
        'species': 'human', 'state': 'review',
        'method': 'systematic review of >100 imaging publications; taxonomy of adipose tissue depots',
        'n': 'not applicable', 'ages': 'not applicable', 'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Adipose Tissue Quantification by Imaging Methods: A Proposed '
                                               'Classification'},
        'used_for': 'the depot taxonomy only (total adipose tissue = subcutaneous + visceral + intermuscular + '
                    'other internal). No number is carried from it.',
    },
    'shen2009': {
        'citation': 'Shen W, Punyanitya M, Silva AM, Chen J, Gallagher D, Sardinha LB, Allison DB, Heymsfield SB. '
                    'Sexual dimorphism of adipose tissue distribution across the lifespan: a cross-sectional '
                    'whole-body magnetic resonance imaging study. Nutr Metab (Lond) 2009;6:17.',
        'doi': '10.1186/1743-7075-6-17', 'pmid': '19371437',
        'species': 'human', 'state': 'in vivo',
        'method': '1.5 T whole-body MRI, 10 mm slices at 40 mm intervals, ~40 slices, manual segmentation of '
                  'total-body SAT and VAT; compartment volumes from slice areas',
        'n': '164 adult males (and 188 adult females), BMI 25.6 +/- 3.7 kg/m2',
        'ages': 'adults >= 18 y, mean 37.8 +/- 14.8 y', 'verification': 'fulltext (PMC2678136, Table 1)',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Sexual dimorphism of adipose tissue distribution across the '
                                               'lifespan: a cross-sectional whole-body magnetic resonance '
                                               'imaging study'},
    },
    'gallagher2005': {
        'citation': 'Gallagher D, Kuznia P, Heshka S, Albu J, Heymsfield SB, Goodpaster B, Visser M, Harris TB. '
                    'Adipose tissue in muscle: a novel depot similar in size to visceral adipose tissue. '
                    'Am J Clin Nutr 2005;81(4):903-910.',
        'doi': '10.1093/ajcn/81.4.903', 'pmid': '15817870',
        'species': 'human', 'state': 'in vivo',
        'method': '1.5 T whole-body MRI, 10 mm slices at 40 mm intervals; TAT, SAT, VAT and IMAT. IMAT is the '
                  'adipose visible between muscle groups and beneath the muscle fascia',
        'n': '111 men (39 African American, 13 Asian, 59 white) and 227 women',
        'ages': 'white men 44.5 +/- 16.3 y, BMI 25.8 +/- 3.8', 'verification': 'fulltext (PMC1482784, Table 1)',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Adipose tissue in muscle: a novel depot similar in size to '
                                               'visceral adipose tissue'},
    },
    'ruan2007': {
        'citation': 'Ruan XY, Gallagher D, Harris T, Albu J, Heymsfield S, Kuznia P, Heshka S. Estimating whole '
                    'body intermuscular adipose tissue from single cross-sectional magnetic resonance images. '
                    'J Appl Physiol 2007;102(2):748-754.',
        'doi': '10.1152/japplphysiol.00304.2006', 'pmid': '17053107',
        'species': 'human', 'state': 'in vivo',
        'method': '1.5 T whole-body MRI, 10 mm slices at 40 mm intervals; whole-body IMAT',
        'n': '39 Caucasian men, BMI 27.1 +/- 3.8', 'ages': '45.2 +/- 14.6 y',
        'verification': 'fulltext (PMC2758818, cohort table)',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Estimating whole body intermuscular adipose tissue from single '
                                               'cross-sectional magnetic resonance images'},
    },
    'kelley2000': {
        'citation': 'Kelley DE, Thaete FL, Troost F, Huwe T, Goodpaster BH. Subdivisions of subcutaneous '
                    'abdominal adipose tissue and insulin resistance. Am J Physiol Endocrinol Metab '
                    '2000;278(5):E941-E948.',
        'doi': '10.1152/ajpendo.2000.278.5.E941', 'pmid': '10780952',
        'species': 'human', 'state': 'in vivo',
        'method': 'cross-sectional abdominal CT; SAT partitioned at the fascial plane within subcutaneous fat '
                  'into superficial SAT (above) and deep SAT (below)',
        'n': '47 lean and obese glucose-tolerant men and women', 'ages': 'not read',
        'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Subdivisions of subcutaneous abdominal adipose tissue and '
                                               'insulin resistance'},
        'used_for': 'the existence and metabolic distinctness of the superficial/deep SAT split at the '
                    'membranous (Scarpa) fascial plane. The abstract gives NO volume ratio, so no depth split '
                    'fraction is carried from it.',
    },
    'storchle2018': {
        'citation': 'Stoerchle P, Mueller W, Sengeis M, Lackner S, Holasek S, Fuerhapter-Rieger A. Measurement '
                    'of mean subcutaneous fat thickness: eight standardised ultrasound sites compared to 216 '
                    'randomly selected sites. Sci Rep 2018;8(1):16268.',
        'doi': '10.1038/s41598-018-34213-0', 'pmid': '30389952',
        'species': 'human', 'state': 'in vivo',
        'method': 'standardised ultrasound of SAT thickness at 216 randomly distributed whole-body sites, '
                  '2160 measurements; whole-body SAT volume as body surface area x calibrated mean thickness',
        'n': '10 participants, BMI < 28.5 kg/m2', 'ages': 'not read', 'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Measurement of mean subcutaneous fat thickness: eight '
                                               'standardised ultrasound sites compared to 216 randomly '
                                               'selected sites'},
    },
    'machann2005': {
        'citation': 'Machann J, Thamer C, Schnoedt B, Haap M, Haring HU, Claussen CD, Stumvoll M, Fritsche A, '
                    'Schick F. Standardized assessment of whole body adipose tissue topography by MRI. '
                    'J Magn Reson Imaging 2005;21(4):455-462.',
        'doi': '10.1002/jmri.20292', 'pmid': '15778954',
        'species': 'human', 'state': 'in vivo',
        'method': 'T1-weighted whole-body MRI; standardised adipose tissue profiles along the body axis',
        'n': '80 volunteers (40 male, 40 female) at increased risk of type 2 diabetes', 'ages': 'not read',
        'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Standardized assessment of whole body adipose tissue topography '
                                               'by MRI'},
        'used_for': 'corroboration that whole-body adipose topography is region-dependent and that males carry '
                    'less subcutaneous and more visceral adipose than BMI-matched females. The abstract carries '
                    'no per-region litre values, so none is taken.',
    },
    'wiig2012': {
        'citation': 'Wiig H, Swartz MA. Interstitial fluid and lymph formation and transport: physiological '
                    'regulation and roles in inflammation and cancer. Physiol Rev 2012;92(3):1005-1060.',
        'doi': '10.1152/physrev.00037.2011', 'pmid': '22811424',
        'species': 'human and animal', 'state': 'review',
        'method': 'review of interstitial fluid composition, transport and biomechanics',
        'n': 'not applicable', 'ages': 'not applicable', 'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Interstitial Fluid and Lymph Formation and Transport: '
                                               'Physiological Regulation and Roles in Inflammation and Cancer'},
        'used_for': 'the definition of the interstitium as fluid, protein, solutes and extracellular matrix. '
                    'No density and no compartment volume is taken from the abstract.',
    },
    'aukland1993': {
        'citation': 'Aukland K, Reed RK. Interstitial-lymphatic mechanisms in the control of extracellular '
                    'fluid volume. Physiol Rev 1993;73(1):1-78.',
        'doi': '10.1152/physrev.1993.73.1.1', 'pmid': '8419962',
        'species': 'human and animal', 'state': 'review',
        'method': 'review of interstitial fluid pressure measurement and interstitial pressure-volume curves',
        'n': 'not applicable', 'ages': 'not applicable', 'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Interstitial-lymphatic mechanisms in the control of '
                                               'extracellular fluid volume'},
    },
    'alkhouli2013': {
        'citation': 'Alkhouli N, Mansfield J, Green E, Bell J, Knight B, Liversedge N, Tham JC, Welbourn R, '
                    'Shore AC, Kos K, Winlove CP. The mechanical properties of human adipose tissues and their '
                    'relationships to the structure and composition of the extracellular matrix. Am J Physiol '
                    'Endocrinol Metab 2013;305(12):E1427-E1435.',
        'doi': '10.1152/ajpendo.00111.2013', 'pmid': '24105412',
        'species': 'human', 'state': 'fresh surgical samples',
        'method': 'uniaxial tension to 30% strain with stress relaxation; initial and final tangent moduli',
        'n': '44 subjects; 19 paired subcutaneous/omental comparisons', 'ages': 'not read',
        'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'The mechanical properties of human adipose tissues and their '
                                               'relationships to the structure and composition of the '
                                               'extracellular matrix'},
    },
    'chakouch2015': {
        'citation': 'Chakouch MK, Charleux F, Bensamoun SF. Quantifying the elastic property of nine thigh '
                    'muscles using magnetic resonance elastography. PLoS One 2015;10(9):e0138873.',
        'doi': '10.1371/journal.pone.0138873', 'pmid': '26397730',
        'species': 'human', 'state': 'in vivo, at rest (passive)',
        'method': '1.5 T MR elastography; shear modulus of nine thigh muscles and of subcutaneous adipose',
        'n': '29 healthy volunteers', 'ages': 'mean 26 +/- 3.41 y', 'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Quantifying the Elastic Property of Nine Thigh Muscles Using '
                                               'Magnetic Resonance Elastography'},
        'note': 'MRE shear modulus is dynamic at the driver frequency, not a quasi-static neo-Hookean mu',
    },
    'sommer2013': {
        'citation': 'Sommer G, Eder M, Kovacs L, Pathak H, Bonitz L, Mueller C, Regitnig P, Holzapfel GA. '
                    'Multiaxial mechanical properties and constitutive modeling of human adipose tissue: a '
                    'basis for preoperative simulations in plastic and reconstructive surgery. Acta Biomater '
                    '2013;9(11):9036-9048.',
        'doi': '10.1016/j.actbio.2013.06.011', 'pmid': '23811521',
        'species': 'human', 'state': 'fresh abdominal adipose tissue',
        'method': 'biaxial tension and triaxial shear; an anisotropic (Gasser-Ogden-Holzapfel type) '
                  'hyperelastic model was fitted',
        'n': 'not read', 'ages': 'not read',
        'verification': 'abstract only; the fitted coefficients are NOT in the abstract and were not retrieved',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Multiaxial mechanical properties and constitutive modeling of '
                                               'human adipose tissue: A basis for preoperative simulations in '
                                               'plastic and reconstructive surgery'},
    },
    'calvo2018': {
        'citation': 'Calvo-Gallego JL, Dominguez J, Gomez Cia T, Gomez Ciriza G, Martinez-Reina J. Comparison '
                    'of different constitutive models to characterize the viscoelastic properties of human '
                    'abdominal adipose tissue. A pilot study. J Mech Behav Biomed Mater 2018;80:293-302.',
        'doi': '10.1016/j.jmbbm.2018.02.013', 'pmid': '29455039',
        'species': 'human', 'state': 'ex vivo, human abdominal adipose from plastic surgery',
        'method': 'uniaxial compression stress relaxation; quasi-linear and internal-variable viscoelastic '
                  'models each with four strain-energy functions (5-term polynomial, first-order Ogden, '
                  'isotropic Gasser-Ogden-Holzapfel, neo-Hookean plus exponential)',
        'n': 'pilot study; n not in the abstract', 'ages': 'not read',
        'verification': 'abstract only; the abstract states the internal-variables model with the Ogden '
                        'function fits best but gives NO coefficient values, and they were not retrieved',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Comparison of different constitutive models to characterize the '
                                               'viscoelastic properties of human abdominal adipose tissue. A '
                                               'pilot study'},
    },
    'geerligs2008': {
        'citation': 'Geerligs M, Peters GW, Ackermans PA, Oomens CW, Baaijens FP. Linear viscoelastic behavior '
                    'of subcutaneous adipose tissue. Biorheology 2008;45(6):677-688.',
        'doi': '10.3233/BIR-2008-0517', 'pmid': '19065014',
        'species': 'PORCINE', 'state': 'ex vivo, fresh and snap-frozen',
        'method': 'rotational rheometer, parallel plate, shear in the linear regime up to 0.1% strain',
        'n': 'not read', 'ages': 'not read', 'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Linear viscoelastic behavior of subcutaneous adipose tissue'},
    },
    'otsuka2018': {
        'citation': 'Otsuka S, Yakura T, Ohmichi Y, Ohmichi M, Naito M, Nakano T, Kawakami Y. Site specificity '
                    'of mechanical and structural properties of human fascia lata and their gender '
                    'differences: A cadaveric study. J Biomech 2018;77:69-75.',
        'doi': '10.1016/j.jbiomech.2018.06.018', 'pmid': '29970229',
        'species': 'human', 'state': 'cadaveric (embalmed)',
        'method': 'caliper thickness at four thigh sites; uniaxial tensile tests longitudinal and transverse',
        'n': '17 legs of 12 cadavers (6 male, 6 female)', 'ages': '75-92 y', 'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Site specificity of mechanical and structural properties of '
                                               'human fascia lata and their gender differences: A cadaveric '
                                               'study'},
    },
    'bonaldi2023': {
        'citation': 'Bonaldi L, Berardo A, Pirri C, Stecco C, Carniel EL, Fontanella CG. Mechanical '
                    'characterization of human fascia lata: uniaxial tensile tests from fresh-frozen cadaver '
                    'samples and constitutive modelling. Bioengineering (Basel) 2023;10(2):226.',
        'doi': '10.3390/bioengineering10020226', 'pmid': '36829719',
        'species': 'human', 'state': 'fresh-frozen at -80 C for under one year, first thaw, tested within 12 h',
        'method': 'ten preconditioning cycles at 1%/s, then uniaxial tension to failure at 0.5%/s; a '
                  'Holzapfel-Gasser-Ogden two-fibre-family model fitted per subject',
        'n': '4 cadavers (2 male, 2 female); anterior thigh, 10 cm distal to the ASIS',
        'ages': '54-89 y', 'verification': 'fulltext (PMC9952725, Tables 1 and 3)',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Mechanical Characterization of Human Fascia Lata: Uniaxial '
                                               'Tensile Tests from Fresh-Frozen Cadaver Samples and '
                                               'Constitutive Modelling'},
    },
    'berardo2024': {
        'citation': 'Berardo A, Bonaldi L, Stecco C, Fontanella CG. Biomechanical properties of the human '
                    'superficial fascia: site-specific variability and anisotropy of abdominal and thoracic '
                    'regions. J Mech Behav Biomed Mater 2024;157:106637.',
        'doi': '10.1016/j.jmbbm.2024.106637', 'pmid': '38914036',
        'species': 'human', 'state': 'cadaveric',
        'method': 'uniaxial tension to failure and stress relaxation, cranio-caudal and latero-medial, on '
                  'abdominal and thoracic (back) superficial fascia',
        'n': '4 subjects', 'ages': 'not read', 'verification': 'abstract',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'Biomechanical properties of the human superficial fascia: '
                                               'Site-specific variability and anisotropy of abdominal and '
                                               'thoracic regions'},
    },
    'woodard1986': {
        'citation': 'Woodard HQ, White DR. The composition of body tissues. Br J Radiol 1986;59(708):1209-1218.',
        'doi': '10.1259/0007-1285-59-708-1209', 'pmid': '3801800',
        'species': 'human', 'state': 'reference composition, not a specimen',
        'method': 'reassessment of ICRP (1975) composition data; water, lipid, protein, carbohydrate, ash, '
                  'elemental composition, mass and electron densities for 56 body tissues',
        'n': 'reference table', 'ages': 'healthy adults', 'verification': 'abstract only',
        'identifier_check': {'resolver': 'crossref+pubmed', 'status': 'resolved',
                             'returned_title': 'The composition of body tissues'},
        'used_for': 'NOTHING NUMERIC. This is the table that would supply a connective-tissue mass density (it '
                    'and ICRU Report 46 are the only reference tabulations that carry one), and it could not be '
                    'retrieved. It is recorded so the gap has an address rather than being silently filled.',
    },
}

# --------------------------------------------------------------- densities
# Every density used anywhere in this file, in one place, with its tier.
RHO = {
    'adipose': {'value': 950.0, 'tier': 'transferred', 'source': 'icru44_nist',
                'note': 'ICRU-44 adipose tissue, 0.950 g/cm3, re-read on the NIST table 2026-09-06. Adipose '
                        'TISSUE is about 80% fat in adults (ICRP 89 para 545), which is why it exceeds the '
                        'density of fat itself'},
    'skin': {'value': 1100.0, 'tier': 'transferred', 'source': 'icrp89',
             'note': 'ICRP 89 para 514, skin about 1.1 g/cm3; the value ICRP itself uses to convert mass '
                     'thickness to linear thickness. The NIST ICRU-44 table has no skin row'},
    'soft_tissue_max': {'value': 1060.0, 'tier': 'transferred', 'source': 'icru44_nist',
                        'note': 'ICRU-44 soft tissue and whole blood, 1.060 g/cm3. No ICRU-44 soft tissue '
                                'exceeds this. Used as the UPPER bound wherever a real density is absent'},
    'water': {'value': 1000.0, 'tier': 'transferred', 'source': 'icru44_nist',
              'note': 'ICRU-44 liquid water, 1.000 g/cm3. Used as the LOWER bound for the fluid-rich residual'},
    'residual_central': {'value': 1030.0, 'tier': 'assumed', 'source': None,
                         'note': 'ASSUMED. The midpoint of [water 1000, ICRU-44 soft tissue 1060] for the '
                                 'loose-connective-tissue-plus-interstitial-fluid residual. No measured '
                                 'density for loose areolar connective tissue or for interstitial fluid was '
                                 'retrieved. Every mass this file reports is also reported at both bounds'},
    'fascia': {'value': 1060.0, 'tier': 'transferred', 'source': 'icru44_nist',
               'note': 'PROXY. No density for human fascia was retrieved. ICRU-44 soft tissue 1.060 is used as '
                       'the upper bound of the soft-tissue band; a collagen-dense sheet plausibly exceeds it, '
                       'so this understates the fascia mass by an unknown amount. Fascia is 2.7% of the fill '
                       'volume so the effect on the total is under 0.03 kg'},
}

# --------------------------------------------------------------- literature depot volumes
# Adult male whole-body adipose tissue depots. All in m3 of TISSUE, not of fat.
DEPOTS = {
    'subcutaneous_adipose': {
        'volume_m3': 0.0165, 'tier': 'transferred', 'source': 'shen2009',
        'reported': 'SAT 16.5 +/- 6.9 L, 164 adult men, BMI 25.6 +/- 3.7, whole-body MRI',
        'corroboration': 'gallagher2005 white men (n=59, BMI 25.8): SAT 17.5 +/- 6.5 kg = 18.4 L at 950 kg/m3. '
                         'storchle2018 ultrasound over 216 whole-body sites: SAT 3.2-12.4 kg in ten subjects '
                         'with mean SAT thickness 3-10 mm, i.e. 4.9-13.3% of body mass, which brackets the '
                         '15.7 kg / 22% implied here from above only for the leanest subjects',
    },
    'visceral_adipose': {
        'volume_m3': 0.0021, 'tier': 'transferred', 'source': 'shen2009',
        'reported': 'VAT 2.1 +/- 1.8 L, same 164 adult men',
        'corroboration': 'gallagher2005 white men: VAT 2.6 +/- 1.9 kg = 2.7 L at 950 kg/m3',
    },
    'intermuscular_adipose': {
        'volume_m3': 0.0009, 'tier': 'transferred', 'source': 'ruan2007',
        'reported': 'whole-body IMAT 0.9 +/- 0.5 L, 39 Caucasian men, BMI 27.1 +/- 3.8, whole-body MRI. IMAT '
                    'is adipose between muscle groups and beneath the muscle fascia',
        'corroboration': 'gallagher2005 white men: IMAT 0.74 +/- 0.46 kg = 0.78 L at 950 kg/m3',
    },
}

# ICRP 89 reference adult male, read from the full text.
ICRP89 = {
    'total_body_mass_kg': 73.0, 'height_m': 1.76, 'body_surface_m2': 1.90,
    'adipose_tissue_kg': 18.2, 'separable_adipose_kg': 14.5,
    'interstitial_adipose_fraction_of_total_adipose': 0.08,
    'fat_fraction_of_adipose_tissue_adult': 0.80, 'fat_fraction_of_adipose_tissue_range': [0.60, 0.90],
    'non_essential_body_fat_kg': 14.6, 'essential_fat_fraction_of_lbm': 0.02,
    'skin_total_kg': 3.3, 'epidermis_kg': 0.120, 'dermis_kg': 3.180,
    'separable_dense_connective_tissue_kg': 2.6,
    'periarticular_fraction_of_body_mass': 0.02, 'tendon_and_fascia_fraction_of_body_mass': [0.010, 0.015],
    'water_fraction_of_lean_body_mass': 0.73,
    'body_water_distribution_fraction': {'plasma': 0.075, 'interstitial_fluid_and_lymph': 0.20,
                                         'dense_connective_tissue_and_cartilage': 0.075, 'bone': 0.075,
                                         'transcellular': 0.025, 'total_extracellular': 0.45,
                                         'total_intracellular': 0.55},
    'unaccounted_body_mass_fraction': 0.04,
    'citation_note': 'Table 2.8, Table 2.9 note c, Table 2.27, Table 4.5, Table 11.2, paras 105, 109, 545, '
                     '546, 578. Table 4.5 is Edelman (1961) as reported by Forbes (1987)',
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


# --------------------------------------------------------------- composition

def composition(void_m3, dermis_volume_m3):
    """Constituent volumes for the void. The residual is the closure term."""
    lbm_kg = ICRP89['total_body_mass_kg'] - ICRP89['non_essential_body_fat_kg']
    tbw_l = ICRP89['water_fraction_of_lean_body_mass'] * lbm_kg          # 73% of 58.4 kg
    isf_l = ICRP89['body_water_distribution_fraction']['interstitial_fluid_and_lymph'] * tbw_l

    lo, hi = ICRP89['tendon_and_fascia_fraction_of_body_mass']
    fascia_kg = 0.5 * (lo + hi) * ICRP89['total_body_mass_kg']
    fascia_m3 = fascia_kg / RHO['fascia']['value']

    rows = []
    rows.append({
        'constituent': 'dermis_and_epidermis',
        'volume_m3': dermis_volume_m3, 'tier': 'derived', 'source': 'icrp89',
        'basis': 'the canonical skin is a zero-thickness double-sided slab, so the real dermis occupies void. '
                 'Volume = envelope area 1.7812538727 m2 x 1.6 mm epidermis-plus-dermis thickness, taken from '
                 'data/derived/tissue-material-candidate-v1/mass_ledger.json skin_slab_defect',
        'density_kg_m3': RHO['skin']['value'], 'density_key': 'skin',
        'cross_check': 'ICRP 89 Table 2.27 gives adult male total skin 3300 g; this volume at 1100 kg/m3 '
                       'weighs %.3f kg' % (dermis_volume_m3 * RHO['skin']['value']),
    })
    for key in ('subcutaneous_adipose', 'visceral_adipose', 'intermuscular_adipose'):
        d = DEPOTS[key]
        rows.append({
            'constituent': key, 'volume_m3': d['volume_m3'], 'tier': d['tier'], 'source': d['source'],
            'basis': d['reported'], 'corroboration': d['corroboration'],
            'density_kg_m3': RHO['adipose']['value'], 'density_key': 'adipose',
        })
    rows.append({
        'constituent': 'fascia_and_separable_dense_connective_tissue',
        'volume_m3': fascia_m3, 'tier': 'derived', 'source': 'icrp89',
        'basis': 'ICRP 89 para 578: tendons and fascia are about 1-1.5%% of total body mass. Midpoint 1.25%% '
                 'of 73 kg = %.4f kg, divided by the proxy density %.0f kg/m3' % (fascia_kg,
                                                                                  RHO['fascia']['value']),
        'density_kg_m3': RHO['fascia']['value'], 'density_key': 'fascia',
        'unresolved': 'the atlas already names 301 ligaments, 38 fascia, 36 capsules and some tendons, so an '
                      'unknown part of this ICRP compartment already has geometry and is double counted here. '
                      'It is 2.7% of the fill volume, so the double count is at most 0.9 kg and probably much '
                      'less',
    })
    claimed = sum(r['volume_m3'] for r in rows)
    residual = void_m3 - claimed
    rows.append({
        'constituent': 'loose_connective_tissue_and_interstitial_fluid',
        'volume_m3': residual, 'tier': 'derived', 'source': None,
        'basis': 'CLOSURE TERM: void minus every constituent above. It carries loose areolar connective '
                 'tissue, extra-organ interstitial fluid and gel, unmodelled microvasculature, small nerves '
                 'and lymphatics, and the 0.549 L that voxelisation loses on thin structures',
        'density_kg_m3': RHO['residual_central']['value'], 'density_key': 'residual_central',
        'independent_cross_check': {
            'quantity': 'ICRP 89 interstitial fluid and lymph compartment',
            'arithmetic': 'lean body mass = 73.0 - 14.6 = %.1f kg; total body water = 73%% of that = %.3f L '
                          '(ICRP 89 para 105); interstitial fluid and lymph = 20%% of total body water = '
                          '%.3f L (ICRP 89 Table 4.5)' % (lbm_kg, tbw_l, isf_l),
            'value_m3': isf_l / 1000.0, 'closure_residual_m3': residual,
            'ratio': residual / (isf_l / 1000.0),
            'reading': 'the closure residual and the reference-man interstitial fluid compartment agree to '
                       'within 5%. That is SUGGESTIVE, NOT A PROOF: much interstitial fluid is intramuscular '
                       'and intra-organ and therefore already inside atlas geometry (which biases the ICRP '
                       'number up as a comparator), while the residual also holds connective-tissue solids '
                       'and unmodelled vessels (which biases it down). The two errors have opposite sign and '
                       'neither was measured',
        },
    })
    return rows, {'lean_body_mass_kg': lbm_kg, 'total_body_water_l': tbw_l,
                  'interstitial_fluid_and_lymph_l': isf_l, 'fascia_kg': fascia_kg}


# --------------------------------------------------------------- allocation

def waterfill(bins, order, demands):
    """Greedy allocation of constituent demands into ordered capacity bins.

    `bins` maps bin key -> capacity m3. `order` is the bin traversal order for
    each constituent, as a list of (constituent, [bin keys]). Returns
    alloc[bin][constituent] and any unmet demand.
    """
    cap = dict(bins)
    alloc = {b: {} for b in bins}
    unmet = {}
    for name in order:
        need = demands[name]
        for b in order[name]:
            if need <= 1e-15:
                break
            take = min(need, cap[b])
            if take > 0:
                alloc[b][name] = alloc[b].get(name, 0.0) + take
                cap[b] -= take
                need -= take
        if need > 1e-12:
            unmet[name] = need
    return alloc, cap, unmet


def regionalise_by_depth(summary, rows):
    """Fill the measured depth histogram with the composed constituents."""
    hist = summary['primary']['depth_histogram']
    bins, labels = {}, []
    for h in hist:
        if h['volume_m3'] <= 0:
            continue
        k = '%g-%s_mm' % (h['lo_mm'], ('inf' if h['hi_mm'] is None else '%g' % h['hi_mm']))
        bins[k] = h['volume_m3']
        labels.append(k)
    demands = {r['constituent']: r['volume_m3'] for r in rows}
    shallow, deep = labels, list(reversed(labels))
    order = {
        # dermis is the surface layer; subcutaneous adipose sits directly under it
        'dermis_and_epidermis': shallow,
        'subcutaneous_adipose': shallow,
        # visceral adipose is the deepest depot; intermuscular and fascia sit above it
        'visceral_adipose': deep,
        'intermuscular_adipose': deep,
        'fascia_and_separable_dense_connective_tissue': deep,
        # the loose-connective-tissue residual takes whatever is left
        'loose_connective_tissue_and_interstitial_fluid': shallow,
    }
    seq = ['dermis_and_epidermis', 'subcutaneous_adipose', 'visceral_adipose',
           'intermuscular_adipose', 'fascia_and_separable_dense_connective_tissue',
           'loose_connective_tissue_and_interstitial_fluid']
    alloc, cap, unmet = waterfill(bins, {k: order[k] for k in seq}, demands)
    rho = {r['constituent']: r['density_kg_m3'] for r in rows}
    out = []
    for k in labels:
        v = bins[k]
        parts = alloc[k]
        mass = sum(x * rho[c] for c, x in parts.items())
        out.append({
            'depth_bin': k, 'volume_m3': v,
            'fractions': {c: x / v for c, x in sorted(parts.items(), key=lambda kv: -kv[1])},
            'adipose_volume_fraction': sum(x for c, x in parts.items() if 'adipose' in c) / v,
            'density_kg_m3': mass / v if v else None, 'mass_kg': mass,
            'unallocated_m3': cap[k],
        })
    return {'rule': 'greedy water-filling. Shallow-to-deep: dermis, then subcutaneous adipose. '
                    'Deep-to-shallow: visceral adipose, then intermuscular adipose, then fascia. The loose '
                    'connective tissue and interstitial fluid residual then fills whatever capacity remains.',
            'literature_anchor': 'the 0-16 mm void is 15.962 L against a whole-body subcutaneous adipose '
                                 'volume of 16.5 +/- 6.9 L in 164 adult men (shen2009), which is the '
                                 'strongest independent check in this file: the shallow void and the '
                                 'measured SAT compartment are the same size to within 3%. kelley2000 '
                                 'establishes that subcutaneous adipose is itself split at a membranous '
                                 'fascial plane into a superficial and a deep compartment that behave '
                                 'differently, so the 4-8 mm and 8-16 mm bins are not one tissue; no volume '
                                 'ratio for that split was retrieved, so both bins carry the same material',
            'assumption': 'the input artifact carries the depth histogram and the per-system attribution as '
                          'two MARGINALS; their joint distribution is not shipped. The depth ORDER of each '
                          'constituent is an anatomical assumption, not a measurement. Only the bin volumes '
                          'and the constituent totals are constrained',
            'known_artifacts': [
                'the strict ordering leaves the 24-48 mm band 100% loose connective tissue and interstitial '
                'fluid with no adipose at all. That is an artifact of the rule, not anatomy: deep '
                'subcutaneous adipose over the abdomen, flank and buttock reaches well past 32 mm, and '
                'intermuscular adipose lives throughout that band. A joint depth x region field would smear '
                'adipose across 16-64 mm instead of terminating it at 24 mm',
                'visceral adipose is placed deepest-first, which puts it in the 64-128 mm bins. The '
                'per-system view places the same depot on digestive-attributed void instead and is the '
                'better-anchored of the two; where they disagree, prefer the system view for visceral fat '
                'and the depth view for subcutaneous fat',
            ],
            'bins': out, 'unmet_demand_m3': unmet}


def regionalise_by_system(summary, rows):
    """Fill the measured per-nearest-entity-system void with the same constituents."""
    sysvol = dict(summary['primary']['void_volume_by_nearest_entity_system_m3'])
    demands = {r['constituent']: r['volume_m3'] for r in rows}
    order = {
        'dermis_and_epidermis': ['integumentary'],
        'subcutaneous_adipose': ['integumentary', 'muscular', 'skeletal'],
        'visceral_adipose': ['digestive', 'urinary', 'reproductive', 'cardiac', 'venous', 'arterial'],
        'intermuscular_adipose': ['muscular', 'skeletal'],
        'fascia_and_separable_dense_connective_tissue': ['connective', 'muscular', 'skeletal'],
        'loose_connective_tissue_and_interstitial_fluid':
            sorted(sysvol, key=lambda s: -sysvol[s]),
    }
    seq = ['dermis_and_epidermis', 'subcutaneous_adipose', 'visceral_adipose',
           'intermuscular_adipose', 'fascia_and_separable_dense_connective_tissue',
           'loose_connective_tissue_and_interstitial_fluid']
    alloc, cap, unmet = waterfill(sysvol, {k: order[k] for k in seq}, demands)
    rho = {r['constituent']: r['density_kg_m3'] for r in rows}
    out = []
    for s in sorted(sysvol, key=lambda x: -sysvol[x]):
        v = sysvol[s]
        parts = alloc[s]
        mass = sum(x * rho[c] for c, x in parts.items())
        out.append({
            'nearest_entity_system': s, 'void_volume_m3': v,
            'fractions': {c: x / v for c, x in sorted(parts.items(), key=lambda kv: -kv[1])},
            'adipose_volume_fraction': sum(x for c, x in parts.items() if 'adipose' in c) / v if v else 0.0,
            'density_kg_m3': mass / v if v else None, 'mass_kg': mass, 'unallocated_m3': cap[s],
        })
    return {'rule': 'greedy allocation over the measured per-nearest-entity-system void, each constituent '
                    'given the system list anatomy licenses: dermis and subcutaneous adipose to integumentary '
                    'first, visceral adipose to digestive first, intermuscular adipose to muscular, fascia to '
                    'connective then muscular, and the residual to the remainder largest-first.',
            'assumption': 'nearest-entity attribution is a Voronoi label, not an anatomical claim of '
                          'ownership. A voxel nearest to skin is not necessarily subcutaneous fat',
            'literature_anchor': 'digestive-attributed void is 1.833 L against a measured visceral adipose '
                                 'volume of 2.1 +/- 1.8 L in the same 164 men (shen2009), so the mesenteric '
                                 'and omental void very nearly holds the whole visceral depot on its own. '
                                 'machann2005 corroborates that whole-body adipose topography is strongly '
                                 'region-dependent and that men carry less subcutaneous and more visceral '
                                 'adipose than BMI-matched women, but its abstract carries no per-region '
                                 'litre values so no number is taken from it',
            'systems': out, 'unmet_demand_m3': unmet}


# --------------------------------------------------------------- ledger

def ledger(summary, cand, profile, rows, aux):
    p = summary['primary']
    env_voxel = p['interior_volume_m3']
    env_div = summary['envelope']['divergence_volume_m3']
    occ8 = p['occupied_volume_m3']
    void8 = p['void_volume_m3']
    g = cand['geometry_mass']
    mean_rho_owned = g['candidate_mean_density_kg_m3']

    def fill_mass(void, residual_rho):
        m = 0.0
        for r in rows:
            rho = residual_rho if r['constituent'].startswith('loose_') else r['density_kg_m3']
            v = r['volume_m3']
            if r['constituent'].startswith('loose_'):
                v = void - sum(x['volume_m3'] for x in rows if not x['constituent'].startswith('loose_'))
            m += v * rho
        return m

    cases = {}
    for label, void, ent_vol, ent_mass, ent_basis in (
        ('voxel_8mm', void8, occ8, occ8 * mean_rho_owned,
         'exclusive 8 mm occupancy %.6f m3 x the candidate table mean owned density %.4f kg/m3'
         % (occ8, mean_rho_owned)),
        ('voxel_10mm', env_div - g['total_owned_volume_m3'], g['total_owned_volume_m3'], g['candidate_mass_kg'],
         'the per-system rows of data/derived/tissue-material-candidate-v1/mass_ledger.json, summed'),
    ):
        band = {}
        for bk, brho in (('lower_water_1000', RHO['water']['value']),
                         ('central_assumed_1030', RHO['residual_central']['value']),
                         ('upper_icru44_soft_1060', RHO['soft_tissue_max']['value'])):
            fm = fill_mass(void, brho)
            band[bk] = {'residual_density_kg_m3': brho, 'fill_mass_kg': fm,
                        'fill_mean_density_kg_m3': fm / void, 'total_body_mass_kg': ent_mass + fm}
        central = band['central_assumed_1030']
        cases[label] = {
            'entity_volume_m3': ent_vol, 'entity_mass_kg': ent_mass, 'entity_mass_basis': ent_basis,
            'void_volume_m3': void, 'residual_density_band': band,
            'composed_fill_mass_kg': central['fill_mass_kg'],
            'composed_fill_density_kg_m3': central['fill_mean_density_kg_m3'],
            'implied_total_body_mass_kg': central['total_body_mass_kg'],
            'fill_density_required_to_reach_profile_kg_m3': (profile['mass_kg'] - ent_mass) / void,
        }

    c8 = cases['voxel_8mm']
    fat_in_fill = ICRP89['fat_fraction_of_adipose_tissue_adult'] * sum(
        r['volume_m3'] * r['density_kg_m3'] for r in rows if 'adipose' in r['constituent'])
    total = c8['implied_total_body_mass_kg']

    # Siri needs a gas-free volume. Take the lung gas from the candidate table's own respiratory row.
    resp = next(r for r in g['rows'] if r['system'] == 'respiratory')
    lung_tissue_rho = 1050.0
    lung_gas_m3 = resp['owned_volume_m3'] - resp['candidate_mass_kg'] / lung_tissue_rho
    gasfree = env_voxel - lung_gas_m3
    d_gasfree = total / gasfree
    siri_gasfree = 495.0 / (d_gasfree / 1000.0) - 450.0
    d_profile = profile['mass_kg'] / env_voxel
    siri_profile_nogas = 495.0 / (d_profile / 1000.0) - 450.0
    d_profile_gasfree = profile['mass_kg'] / gasfree
    siri_profile_gasfree = 495.0 / (d_profile_gasfree / 1000.0) - 450.0

    return {
        'geometry': {
            'envelope_divergence_m3': env_div, 'envelope_voxel_interior_m3': env_voxel,
            'envelope_area_m2': summary['envelope']['surface_area_m2'],
            'exclusive_occupied_8mm_m3': occ8, 'void_8mm_m3': void8,
            'void_fraction': p['void_fraction_of_interior'],
            'largest_connected_void_m3': p['void_components']['largest_volume_m3'],
            'tier': 'measured', 'note': 'measured on the BodyParts3D specimen; divergence and 8 mm voxel '
                                        'interior agree to 0.026%',
        },
        'cases': cases,
        'arithmetic_8mm': [
            'envelope interior            %12.6f m3   measured' % env_voxel,
            'entity-claimed (8 mm excl.)  %12.6f m3   measured' % occ8,
            'void                         %12.6f m3   measured  (%.3f%%)' % (void8, 100 * void8 / env_voxel),
        ] + ['  %-46s %10.6f m3 x %6.1f = %8.4f kg' % (
            r['constituent'], r['volume_m3'], r['density_kg_m3'], r['volume_m3'] * r['density_kg_m3'])
             for r in rows] + [
            '  %-46s %10.6f m3            = %8.4f kg' % ('FILL TOTAL', void8, c8['composed_fill_mass_kg']),
            'fill mean density            %12.4f kg/m3' % c8['composed_fill_density_kg_m3'],
            'entity mass                  %12.4f kg' % c8['entity_mass_kg'],
            'TOTAL BODY MASS              %12.4f kg' % total,
            'profile mass                 %12.4f kg' % profile['mass_kg'],
            'deficit                      %12.4f kg   (%.2f%%)' % (
                profile['mass_kg'] - total, 100 * (profile['mass_kg'] - total) / profile['mass_kg']),
        ],
        'body_fat': {
            'adipose_tissue_in_fill_kg': sum(r['volume_m3'] * r['density_kg_m3'] for r in rows
                                             if 'adipose' in r['constituent']),
            'fat_in_that_adipose_kg': fat_in_fill,
            'fat_fraction_lower_bound': fat_in_fill / total,
            'lower_bound_note': 'a LOWER BOUND. It counts only fat inside the separable adipose of the fill, '
                                'at ICRP 89 para 545 80% fat by mass. It omits yellow marrow (inside atlas '
                                'bone geometry), the inseparable interstitial adipose ICRP puts at 8% of '
                                'total adipose tissue, and essential fat at 2% of lean body mass',
            'completed_estimate_kg': fat_in_fill
                                     + ICRP89['fat_fraction_of_adipose_tissue_adult']
                                       * ICRP89['interstitial_adipose_fraction_of_total_adipose']
                                       * ICRP89['adipose_tissue_kg']
                                     + ICRP89['fat_fraction_of_adipose_tissue_adult']
                                       * (ICRP89['adipose_tissue_kg'] - ICRP89['separable_adipose_kg']
                                          - ICRP89['interstitial_adipose_fraction_of_total_adipose']
                                            * ICRP89['adipose_tissue_kg'])
                                     + ICRP89['essential_fat_fraction_of_lbm'] * aux['lean_body_mass_kg'],
        },
        'densitometry': {
            'method': 'Siri two-compartment, %fat = 495/D - 450, D in g/cm3. Siri requires a GAS-FREE body '
                      'volume, so lung gas is removed first',
            'lung_gas_m3': lung_gas_m3,
            'lung_gas_basis': 'the candidate table respiratory row, %.6f m3 owned at %.0f kg/m3 = %.4f kg; at '
                              'a gas-free lung tissue-plus-blood density of %.0f kg/m3 that mass occupies '
                              '%.6f m3, so %.6f m3 is gas'
                              % (resp['owned_volume_m3'], resp['density_kg_m3'], resp['candidate_mass_kg'],
                                 lung_tissue_rho, resp['candidate_mass_kg'] / lung_tissue_rho, lung_gas_m3),
            'gas_free_volume_m3': gasfree,
            'composed_body': {'mass_kg': total, 'density_kg_m3': d_gasfree, 'siri_fat_percent': siri_gasfree},
            'profile_body': {'mass_kg': profile['mass_kg'],
                             'density_with_gas_kg_m3': d_profile,
                             'siri_fat_percent_with_gas': siri_profile_nogas,
                             'density_gas_free_kg_m3': d_profile_gasfree,
                             'siri_fat_percent_gas_free': siri_profile_gasfree},
            'reading': 'the composed body is internally consistent: %.1f%% fat by densitometry against a '
                       '%.1f%% compositional lower bound and a %.1f%% completed compositional estimate. The '
                       'profile is not: forcing 77.1107 kg into the same gas-free volume gives %.1f%% fat, '
                       'and into the raw envelope gives %.1f%%.'
                       % (siri_gasfree, 100 * fat_in_fill / total,
                          100 * (fat_in_fill
                                 + ICRP89['fat_fraction_of_adipose_tissue_adult']
                                   * ICRP89['interstitial_adipose_fraction_of_total_adipose']
                                   * ICRP89['adipose_tissue_kg']
                                 + ICRP89['fat_fraction_of_adipose_tissue_adult']
                                   * (ICRP89['adipose_tissue_kg'] - ICRP89['separable_adipose_kg']
                                      - ICRP89['interstitial_adipose_fraction_of_total_adipose']
                                        * ICRP89['adipose_tissue_kg'])
                                 + ICRP89['essential_fat_fraction_of_lbm'] * aux['lean_body_mass_kg']) / total,
                          siri_profile_gasfree, siri_profile_nogas),
        },
        'verdict': {
            'answer': 'NO. The profile 77.1107029 kg cannot be reconciled with this geometry, and it is the '
                      'PROFILE that must yield.',
            'composed_total_body_mass_kg': total,
            'partition_sensitivity_kg': [cases['voxel_10mm']['implied_total_body_mass_kg'], total],
            'residual_density_sensitivity_kg': [
                c8['residual_density_band']['lower_water_1000']['total_body_mass_kg'],
                c8['residual_density_band']['upper_icru44_soft_1060']['total_body_mass_kg']],
            'implied_bmi': total / profile['height_m'] ** 2,
            'profile_bmi': profile['mass_kg'] / profile['height_m'] ** 2,
            'grounds': [
                'the composed fill weighs %.4f kg over %.6f m3, a mean density of %.1f kg/m3. Reaching the '
                'profile would need the void to weigh %.1f kg/m3, above ICRU-44 soft tissue 1060 and every '
                'other soft tissue in that table; only bone is denser'
                % (c8['composed_fill_mass_kg'], void8, c8['composed_fill_density_kg_m3'],
                   c8['fill_density_required_to_reach_profile_kg_m3']),
                'the answer is insensitive to the two free choices: swapping the 8 mm partition for the 10 mm '
                'one moves the total by %.3f kg, and sweeping the residual density across the full '
                '[1000, 1060] band moves it by %.3f kg'
                % (abs(cases['voxel_10mm']['implied_total_body_mass_kg'] - total),
                   c8['residual_density_band']['upper_icru44_soft_1060']['total_body_mass_kg']
                   - c8['residual_density_band']['lower_water_1000']['total_body_mass_kg']),
                'the geometry is measured on this specimen and two independent integrations of it agree to '
                '0.026%; the profile mass is an inherited BioGears StandardMale constant whose own file '
                'records it as a "generic synthesized reference profile, not a measured patient". Only the '
                'height in that profile comes from the atlas',
                'the composed body reproduces a coherent body composition at %.1f kg, BMI %.1f, %.1f%% fat by '
                'densitometry. The profile at 77.1107 kg in this envelope is %.1f%% fat by the same relation, '
                'which is not a body'
                % (total, total / profile['height_m'] ** 2, siri_gasfree, siri_profile_nogas),
            ],
            'what_to_change': 'hold the geometry and the 21%% fat intent; replace mass_kg 77.1107029 with '
                              '%.4f kg (BMI %.2f). Alternatively keep 77.1107 kg and declare the profile a '
                              'different subject from the BodyParts3D specimen, in which case no mass ledger '
                              'over this geometry can ever close.'
                              % (total, total / profile['height_m'] ** 2),
        },
    }


# --------------------------------------------------------------- mechanics

def mechanics():
    """Per-constituent mechanics in the hyperelastic forms this project records."""
    return {
        'adipose_tissue': {
            'applies_to': ['subcutaneous_adipose', 'visceral_adipose', 'intermuscular_adipose'],
            'share_of_void': 'the three adipose depots are 19.500 L, 60.7% of the 32.117 L void',
            'density': dict(RHO['adipose'], unit='kg/m3'),
            'linear': {
                'young_modulus_initial_subcutaneous': {
                    'value': 1600.0, 'unit': 'Pa', 'tier': 'transferred', 'source': 'alkhouli2013',
                    'note': 'human subcutaneous, initial tangent 1.6 +/- 0.8 kPa'},
                'young_modulus_final_subcutaneous': {
                    'value': 11700.0, 'unit': 'Pa', 'tier': 'transferred', 'source': 'alkhouli2013',
                    'note': 'human subcutaneous at 30% strain, 11.7 +/- 6.4 kPa. The 7x rise from initial to '
                            'final is why a linear modulus is the wrong description'},
                'young_modulus_initial_omental': {
                    'value': 2900.0, 'unit': 'Pa', 'tier': 'transferred', 'source': 'alkhouli2013',
                    'note': 'human omental (visceral) 2.9 +/- 1.5 kPa initial, 32 +/- 15.6 kPa final. Visceral '
                            'adipose is roughly 1.8x stiffer than subcutaneous in the same study, which is the '
                            'only depot-resolved mechanical contrast retrieved'},
                'shear_modulus_mre_subcutaneous': {
                    'value': 3040.0, 'unit': 'Pa', 'tier': 'transferred', 'source': 'chakouch2015',
                    'note': 'human subcutaneous in vivo 3.04 +/- 0.12 kPa, the softest tissue in that series. '
                            'DYNAMIC at the driver frequency, not a quasi-static shear modulus'},
                'shear_modulus_rheometer': {
                    'value': 7500.0, 'unit': 'Pa', 'tier': 'transferred', 'source': 'geerligs2008',
                    'note': 'PORCINE subcutaneous, 7.5 kPa at 10 rad/s and 37 C in the linear regime up to '
                            '0.1% strain. Species transfer AND a small-strain dynamic modulus'},
                'poisson_ratio': {'value': None, 'unit': '1', 'tier': 'absent', 'source': None,
                                  'note': 'no measurement retrieved for adipose tissue of any species'},
            },
            'hyperelastic': {
                'neo_hookean_mu': {
                    'value': 533.3, 'unit': 'Pa', 'tier': 'derived', 'source': 'alkhouli2013',
                    'note': 'mu = E/3 under incompressibility, from the 1.6 kPa INITIAL tangent modulus of '
                            'human subcutaneous adipose. Valid only in the toe region; the same tissue is 11.7 '
                            'kPa at 30% strain, so a neo-Hookean fit to the full range would land far higher. '
                            'The source reports E, not mu'},
                'ogden': {'value': None, 'unit': None, 'tier': 'absent', 'source': 'calvo2018',
                          'note': 'calvo2018 tested five model families on human abdominal adipose and reports '
                                  'the internal-variables viscoelastic model with a first-order Ogden function '
                                  'as the best fit, but the abstract carries no coefficients and the full text '
                                  'was not retrieved. NOT filled with a plausible number'},
                'mooney_rivlin': {'value': None, 'unit': None, 'tier': 'absent', 'source': None,
                                  'note': 'no Mooney-Rivlin coefficient set for adipose of any species was '
                                          'retrieved'},
                'hgo': {'value': None, 'unit': None, 'tier': 'absent', 'source': 'sommer2013',
                        'note': 'sommer2013 fitted a Gasser-Ogden-Holzapfel type model to multiaxial human '
                                'abdominal adipose but the coefficients are not in the abstract and were not '
                                'retrieved'},
            },
            'viscoelastic': {
                'stress_relaxation': {'value': None, 'unit': None, 'tier': 'absent', 'source': 'alkhouli2013',
                                      'note': 'stress relaxation was performed in alkhouli2013 but no decay '
                                              'constant is in the abstract'},
            },
        },
        'loose_areolar_connective_tissue_and_interstitial_fluid': {
            'applies_to': ['loose_connective_tissue_and_interstitial_fluid'],
            'share_of_void': 'the closure residual is 8.906 L, 27.7% of the void, the second largest '
                             'constituent and the one with the emptiest row',
            'density': {'value': None, 'unit': 'kg/m3', 'tier': 'absent', 'source': None,
                        'note': 'NO measured density for loose areolar connective tissue or for interstitial '
                                'fluid was retrieved. The ledger uses the assumed midpoint 1030 of the '
                                '[ICRU-44 liquid water 1000, ICRU-44 soft tissue 1060] band and reports every '
                                'total at both bounds. woodard1986 and ICRU Report 46 are the two reference '
                                'tabulations that would carry a connective-tissue density and neither was '
                                'retrieved'},
            'linear': {
                'young_modulus': {'value': None, 'unit': 'Pa', 'tier': 'absent', 'source': None,
                                  'note': 'no modulus for loose areolar connective tissue was retrieved from '
                                          'any species. This is the single largest evidential hole in the '
                                          'fill: 27.7% of the void has no mechanical row at all'},
                'poisson_ratio': {'value': None, 'unit': '1', 'tier': 'absent', 'source': None, 'note': 'none'},
            },
            'hyperelastic': {'value': None, 'unit': None, 'tier': 'absent', 'source': None,
                             'note': 'none retrieved'},
            'poroelastic': {
                'interstitial_fluid_pressure': {
                    'value': [-533.0, 0.0], 'unit': 'Pa', 'tier': 'transferred', 'source': 'aukland1993',
                    'note': 'control interstitial fluid pressure in SOFT CONNECTIVE TISSUES is 0 to -4 mmHg '
                            '(0 to -533 Pa); zero or slightly positive in other tissues. -4 mmHg = -533.3 Pa'},
                'structure': {
                    'value': None, 'unit': None, 'tier': 'transferred', 'source': 'wiig2012',
                    'note': 'the interstitium is fluid, protein, solutes and extracellular matrix; the '
                            'gel-like consistency comes from glycosaminoglycans, mainly hyaluronan in soft '
                            'connective tissue (aukland1993). This licenses a poroelastic or biphasic '
                            'treatment but supplies no coefficient'},
                'hydraulic_conductivity': {'value': None, 'unit': 'm4/(N.s)', 'tier': 'absent',
                                           'source': 'wiig2012',
                                           'note': 'no numeric value was retrieved from the abstract'},
            },
        },
        'superficial_fascia': {
            'applies_to': ['fascia_and_separable_dense_connective_tissue'],
            'anatomy': 'the membranous layer within subcutaneous fat, Camper above and Scarpa below; it is the '
                       'plane kelley2000 uses to split superficial from deep subcutaneous adipose',
            'density': {'value': None, 'unit': 'kg/m3', 'tier': 'absent', 'source': None,
                        'note': 'no measurement retrieved; the ledger uses ICRU-44 soft tissue 1060 as a proxy'},
            'linear': {
                'young_modulus_abdomen_latero_medial': {
                    'value': 3190000.0, 'unit': 'Pa', 'tier': 'transferred', 'source': 'berardo2024',
                    'note': 'human abdominal superficial fascia, 3.19 +/- 1.62 MPa, n = 4 subjects'},
                'young_modulus_thorax_latero_medial': {
                    'value': 24870000.0, 'unit': 'Pa', 'tier': 'transferred', 'source': 'berardo2024',
                    'note': 'human thoracic (back) superficial fascia, 24.87 +/- 15.23 MPa. The thorax is '
                            'about 8x the abdomen, so one superficial-fascia modulus cannot serve the body'},
                'anisotropy_ratio': {
                    'value': 2.0, 'unit': '1', 'tier': 'transferred', 'source': 'berardo2024',
                    'note': 'latero-medial over cranio-caudal Young modulus close to 2 in both districts'},
                'ultimate_tensile_strength_abdomen': {
                    'value': 850000.0, 'unit': 'Pa', 'tier': 'transferred', 'source': 'berardo2024',
                    'note': '0.85 +/- 0.39 MPa abdominal against 6.13 +/- 3.11 MPa thoracic, latero-medial'},
                'strain_at_break': {
                    'value': [0.38, 0.47], 'unit': '1', 'tier': 'transferred', 'source': 'berardo2024',
                    'note': '38-47% in both regions with no clear direction dependence'},
                'residual_stress_fraction_300s': {
                    'value': [0.35, 0.38], 'unit': '1', 'tier': 'transferred', 'source': 'berardo2024',
                    'note': '35-38% of peak stress remains after 300 s, no significant direction or district '
                            'difference'},
            },
            'hyperelastic': {'value': None, 'unit': None, 'tier': 'absent', 'source': 'berardo2024',
                             'note': 'no constitutive fit is reported for superficial fascia'},
        },
        'deep_fascia': {
            'applies_to': ['fascia_and_separable_dense_connective_tissue'],
            'anatomy': 'the aponeurotic and epimysial sheets, fascia lata being the measured exemplar. The '
                       'extended atlas marks 38 fascia, 40 tendon sheaths and 78 bursae as the interfaces '
                       'where sliding is licensed',
            'density': {'value': None, 'unit': 'kg/m3', 'tier': 'absent', 'source': None,
                        'note': 'no measurement retrieved'},
            'geometry': {
                'thickness': {'value': [0.0002, 0.0011], 'unit': 'm', 'tier': 'transferred',
                              'source': 'otsuka2018',
                              'note': 'fascia lata 0.2-1.1 mm over four thigh sites; lateral 0.8 +/- 0.2 mm is '
                                      'significantly thicker than the other three at 0.2-0.3 mm'},
                'thickness_corroboration': {'value': 0.00068, 'unit': 'm', 'tier': 'transferred',
                                            'source': 'bonaldi2023',
                                            'note': 'anterior thigh 0.68 +/- 0.14 mm, fresh-frozen, a '
                                                    'different cohort and preservation from otsuka2018'},
            },
            'linear': {
                'young_modulus_longitudinal': {
                    'value': [71600000.0, 275900000.0], 'unit': 'Pa', 'tier': 'transferred',
                    'source': 'otsuka2018',
                    'note': '71.6-275.9 MPa across four thigh sites, highest laterally; stiffness 20-283 N/mm'},
                'young_modulus_transverse': {
                    'value': [3200000.0, 41900000.0], 'unit': 'Pa', 'tier': 'transferred', 'source': 'otsuka2018',
                    'note': '3.2-41.9 MPa, lowest laterally; stiffness 3-16 N/mm. Along-fibre over transverse '
                            'spans roughly 2x to 60x depending on site'},
                'young_modulus_longitudinal_freshfrozen': {
                    'value': [10000000.0, 191900000.0], 'unit': 'Pa', 'tier': 'transferred',
                    'source': 'bonaldi2023',
                    'note': 'per subject and per layer, 10.0 +/- 9.7 to 191.9 +/- 54.5 MPa longitudinal; '
                            '4.5 +/- 1.1 to 119.9 +/- 65.8 MPa transverse. The two layers of one specimen '
                            'differ by up to 10x, so a single-layer fascia model is a choice, not a default'},
                'ultimate_tensile_strength': {
                    'value': [500000.0, 12000000.0], 'unit': 'Pa', 'tier': 'transferred', 'source': 'bonaldi2023',
                    'note': '0.5 +/- 0.0 to 12.0 +/- 1.0 MPa across subjects, layers and directions'},
                'strain_at_break': {
                    'value': 0.164, 'unit': '1', 'tier': 'transferred', 'source': 'bonaldi2023',
                    'note': 'mean 16.4%, range 8.4-23.2%. Roughly a third of the superficial fascia value'},
            },
            'hyperelastic': {
                'law': 'Holzapfel-Gasser-Ogden with TWO fibre families, one per fibrous layer',
                'tier': 'transferred', 'source': 'bonaldi2023',
                'note': 'fitted per subject; the paper reports no pooled set, so all four are carried. Units '
                        'MPa for K, C1, k1; k2 and kappa dimensionless. kappa is the fibre dispersion, '
                        '0 = perfectly aligned',
                'subjects': [
                    {'subject': 'S1', 'K_MPa': 0.0754, 'C1_MPa': 0.6807, 'k1_1_MPa': 2.7894, 'k2_1': 55.9368,
                     'kappa_1': 0.0004, 'k1_2_MPa': 1.9971, 'k2_2': 578.9087, 'kappa_2': 0.1879},
                    {'subject': 'S2', 'K_MPa': 0.0163, 'C1_MPa': 0.3597, 'k1_1_MPa': 4.3575, 'k2_1': 173.3235,
                     'kappa_1': 0.0767, 'k1_2_MPa': 0.7668, 'k2_2': 197.6555, 'kappa_2': 0.0392},
                    {'subject': 'S3', 'K_MPa': 0.0560, 'C1_MPa': 0.1071, 'k1_1_MPa': 5.7467, 'k2_1': 79.3530,
                     'kappa_1': 0.1006, 'k1_2_MPa': 1.4018, 'k2_2': 141.6217, 'kappa_2': 0.0746},
                    {'subject': 'S4', 'K_MPa': 0.0431, 'C1_MPa': 0.4693, 'k1_1_MPa': 5.3137, 'k2_1': 42.1639,
                     'kappa_1': 0.0004, 'k1_2_MPa': 0.9011, 'k2_2': 115.2489, 'kappa_2': 0.0393},
                ],
                'spread': 'C1 spans 6.4x, k2_2 spans 5.0x and kappa_1 spans 250x across four cadavers aged '
                          '54-89. Any single parameter set is one of four, not a population value',
            },
        },
        'dermis': {
            'applies_to': ['dermis_and_epidermis'],
            'density': dict(RHO['skin'], unit='kg/m3'),
            'linear': {'note': 'the moduli are already in data/derived/tissue-material-candidate-v1/'
                               'materials.json under skin_layer (epidermis about 4 MPa, dermis about 40 kPa, '
                               'hypodermis about 15 kPa in vivo; excised whole skin 83.3 +/- 34.9 MPa). Not '
                               'duplicated here'},
            'why_it_is_in_the_void': 'the canonical skin is a zero-thickness double-sided slab of 3.5026 m2 '
                                     'against a 1.7813 m2 envelope, and the three skin-layer entities occupy '
                                     'no voxel in either shipped partition. The dermis therefore has mass with '
                                     'no place to be, and the place it belongs is the shallowest void',
        },
    }


# --------------------------------------------------------------- build / emit

def build():
    summary = json.loads(UNMODELLED.read_text())
    cand = json.loads(CANDIDATE_LEDGER.read_text())
    profile = json.loads(PROFILE.read_text())
    dermis_v = cand['skin_slab_defect']['skin_epidermis_plus_dermis_volume_m3']
    void = summary['primary']['void_volume_m3']
    rows, aux = composition(void, dermis_v)
    doc = {
        'schema': 'ihm.interstitial-composition-prior.v1',
        'status': 'CANDIDATE. Not canonical, not promoted, not consumed by any runtime. No canonical asset '
                  'is read for writing and none is modified.',
        'question': 'the 32.117 L of body interior that belongs to no entity has exact geometry and no '
                    'composition. This assigns it one.',
        'tier_vocabulary': TIERS,
        'sources': SOURCES,
        'reference_male': ICRP89,
        'densities': RHO,
        'composition': {
            'taxonomy': 'the depot names follow shen2003, the accepted imaging classification of adipose '
                        'tissue topography: total adipose tissue = subcutaneous + visceral + intermuscular + '
                        'other internal. woodard1986 is recorded in the densities block as the reference '
                        'tabulation that would supply the missing connective-tissue density',
            'void_volume_m3': void,
            'constituents': rows,
            'derived_quantities': aux,
            'volume_fractions': {r['constituent']: r['volume_m3'] / void for r in rows},
            'adipose_volume_fraction_of_void': sum(r['volume_m3'] for r in rows
                                                   if 'adipose' in r['constituent']) / void,
        },
        'mechanics': mechanics(),
        'mass_ledger': ledger(summary, cand, profile, rows, aux),
        'regional': {
            'by_depth_below_skin': regionalise_by_depth(summary, rows),
            'by_nearest_entity_system': regionalise_by_system(summary, rows),
        },
    }
    doc['tier_census'] = tier_census(doc)
    doc['not_established'] = [
        'nothing here was measured on the BodyParts3D specimen. The honest split is '
        'transferred-versus-absent, not verified-versus-assumed, and every transferred row inherits its '
        'donor cohort: 164 men at BMI 25.6 for the adipose depots, 4 cadavers aged 54-89 for the fascia '
        'constitutive fit, 12 cadavers aged 75-92 for the fascia moduli, a 73 kg 176 cm reference male for '
        'every ICRP quantity',
        'no +/- carried here is a probability interval. SAT is 16.5 +/- 6.9 L, a 42% coefficient of '
        'variation; that dispersion is real population spread and is NOT propagated into the ledger',
        'the depth histogram and the per-nearest-entity-system attribution are two MARGINALS of one 3-D '
        'field. Their joint distribution is not in the input, so both regionalisations are allocations under '
        'a stated ordering rule, not measurements',
        'no measured density exists in this table for loose areolar connective tissue, for interstitial '
        'fluid or for fascia. Three of six constituents rest on a proxy or an assumed density',
        'no hyperelastic coefficient set for human adipose tissue was retrieved. Two papers fitted one and '
        'neither put the coefficients in an accessible abstract',
        'the fascia constituent double counts an unknown part of the 301 ligaments, 38 fascia and 36 capsules '
        'the atlas already names',
        'nearest-entity attribution is a Voronoi label. A void voxel nearest to skin need not be '
        'subcutaneous adipose, and 0.549 L of the "void" is voxelisation loss on thin structures, not anatomy',
        'the ICRP 89 interstitial-fluid cross-check compares a whole-body compartment against an '
        'extra-organ residual. They are not the same quantity and their agreement to 5% is a coincidence of '
        'two errors with opposite sign until someone measures the split',
    ]
    return doc


def tier_census(doc):
    counts, absent_cells, assumed_cells = {}, [], []

    def walk(node, path):
        if isinstance(node, dict):
            if 'tier' in node and isinstance(node.get('tier'), str):
                t = node['tier']
                counts[t] = counts.get(t, 0) + 1
                if t == 'absent':
                    absent_cells.append(path)
                if t == 'assumed':
                    assumed_cells.append(path)
            for k, v in node.items():
                if k in ('sources', 'tier_vocabulary'):
                    continue
                walk(v, path + '/' + str(k))
        elif isinstance(node, list):
            for i, v in enumerate(node):
                walk(v, '%s[%d]' % (path, i))

    for key in ('densities', 'composition', 'mechanics', 'mass_ledger'):
        walk(doc[key], key)
    return {'cells_by_tier': dict(sorted(counts.items())),
            'absent_cells': sorted(absent_cells), 'assumed_cells': sorted(assumed_cells),
            'verified_versus_assumed': {
                'transferred_or_derived_or_measured': sum(v for k, v in counts.items()
                                                          if k in ('transferred', 'derived', 'measured')),
                'assumed': counts.get('assumed', 0), 'absent': counts.get('absent', 0)},
            'sources_cited': len(SOURCES),
            'identifiers_resolved': sum(1 for s in SOURCES.values()
                                        if s.get('identifier_check', {}).get('status') == 'resolved'),
            'note': 'no cell in this table is tier `measured` except the geometry block, which is computed '
                    'from this specimen. Everything compositional is transferred, derived or absent'}


def emit(outdir):
    d = Path(outdir)
    d.mkdir(parents=True, exist_ok=True)
    doc = build()
    artifacts = {
        'sources.json': {'sources': doc['sources'], 'tier_vocabulary': doc['tier_vocabulary']},
        'composition.json': {'schema': doc['schema'], 'reference_male': doc['reference_male'],
                             'densities': doc['densities'], 'composition': doc['composition'],
                             'mechanics': doc['mechanics']},
        'ledger.json': doc['mass_ledger'],
        'regional.json': doc['regional'],
    }
    for name, payload in artifacts.items():
        (d / name).write_text(json.dumps(payload, indent=2, sort_keys=False) + '\n')
    man = {
        'schema': doc['schema'], 'status': doc['status'], 'question': doc['question'],
        'inputs_sha256': {str(p.relative_to(ROOT)): sha(p) for p in INPUTS},
        'tier_vocabulary': doc['tier_vocabulary'],
        'tier_census': doc['tier_census'],
        'headline': {
            'void_m3': doc['composition']['void_volume_m3'],
            'adipose_volume_fraction_of_void': doc['composition']['adipose_volume_fraction_of_void'],
            'composed_fill_density_kg_m3': doc['mass_ledger']['cases']['voxel_8mm'][
                'composed_fill_density_kg_m3'],
            'implied_total_body_mass_kg': doc['mass_ledger']['verdict']['composed_total_body_mass_kg'],
            'profile_mass_kg': 77.1107029,
            'verdict': doc['mass_ledger']['verdict']['answer'],
        },
        'not_established': doc['not_established'],
        'canonical_assets_modified': False,
        'artifacts_sha256': {name: hashlib.sha256((d / name).read_bytes()).hexdigest()
                             for name in artifacts},
    }
    (d / 'manifest.json').write_text(json.dumps(man, indent=2) + '\n')
    return d, doc, man


# --------------------------------------------------------------- identifiers

def verify_identifiers():
    """Live re-resolution of every DOI and PMID. Network. Not part of --self-test."""
    import urllib.request, urllib.parse, time
    out = []
    for key, s in sorted(SOURCES.items()):
        rec = {'source': key, 'doi': s['doi'], 'pmid': s['pmid']}
        if s['doi']:
            ok = False
            for api, host in (('crossref', 'https://api.crossref.org/works/'),
                              ('datacite', 'https://api.datacite.org/dois/')):
                try:
                    with urllib.request.urlopen(host + urllib.parse.quote(s['doi']), timeout=30) as r:
                        j = json.load(r)
                    t = (j['message'].get('title') or [''])[0] if api == 'crossref' \
                        else (j['data']['attributes'].get('titles') or [{}])[0].get('title', '')
                    rec['doi_%s' % api] = {'status': 'resolved', 'title': t}
                    ok = True
                    break
                except Exception as e:
                    rec['doi_%s' % api] = {'status': 'error', 'error': str(e)}
            rec['doi_resolved'] = ok
        if s['pmid']:
            try:
                u = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id='
                     + s['pmid'])
                with urllib.request.urlopen(u, timeout=30) as r:
                    j = json.load(r)['result'][s['pmid']]
                rec['pmid_resolved'] = True
                rec['pmid_title'] = j['title']
                rec['pmid_journal'] = j['source']
            except Exception as e:
                rec['pmid_resolved'] = False
                rec['pmid_error'] = str(e)
        out.append(rec)
        time.sleep(0.25)
    bad = [r for r in out
           if (r.get('doi') and not r.get('doi_resolved')) or (r.get('pmid') and not r.get('pmid_resolved'))]
    return {'checked': len(out), 'failures': len(bad), 'records': out}


# --------------------------------------------------------------- self test

def self_test():
    doc = build()
    comp, led, reg = doc['composition'], doc['mass_ledger'], doc['regional']

    # every source is cited by at least one non-source block, and every cited source exists
    blob = json.dumps({k: v for k, v in doc.items() if k != 'sources'})
    used = {k for k in SOURCES if k in blob}
    unused = sorted(set(SOURCES) - used)
    assert unused == [], 'source recorded but never cited: %s' % unused
    for s in SOURCES.values():
        assert s['doi'] or s['pmid'], s['citation']
        assert s['identifier_check']['status'] == 'resolved', s['citation']
        assert s['verification'], s['citation']
    # the one correction this lane found
    assert 'correction' in SOURCES['icrp89'] and '14527029' in SOURCES['icrp89']['correction']
    assert SOURCES['icrp89']['pmid'] is None
    # species transfers are declared, not hidden
    assert SOURCES['geerligs2008']['species'] == 'PORCINE'

    # composition closes on the measured void exactly
    void = comp['void_volume_m3']
    assert abs(void - 0.032117248) < 1e-12, void
    assert abs(sum(r['volume_m3'] for r in comp['constituents']) - void) < 1e-12
    assert abs(sum(comp['volume_fractions'].values()) - 1.0) < 1e-12
    assert 0.55 < comp['adipose_volume_fraction_of_void'] < 0.65
    # the residual is a closure term and must stay positive, or the depots overfill the void
    resid = next(r for r in comp['constituents'] if r['constituent'].startswith('loose_'))
    assert resid['volume_m3'] > 0, 'literature depots overfill the measured void'
    x = resid['independent_cross_check']
    assert 0.90 < x['ratio'] < 1.10, x['ratio']

    # every constituent carries a density and a tier
    for r in comp['constituents']:
        assert r['density_kg_m3'] > 0 and r['tier'] in TIERS, r
        assert r['basis'], r

    # ledger arithmetic closes at both partitions and both density bounds
    for name, c in led['cases'].items():
        for bk, b in c['residual_density_band'].items():
            assert abs(b['fill_mass_kg'] / c['void_volume_m3'] - b['fill_mean_density_kg_m3']) < 1e-9
            assert abs(c['entity_mass_kg'] + b['fill_mass_kg'] - b['total_body_mass_kg']) < 1e-9
        lo = c['residual_density_band']['lower_water_1000']['total_body_mass_kg']
        hi = c['residual_density_band']['upper_icru44_soft_1060']['total_body_mass_kg']
        assert lo < c['implied_total_body_mass_kg'] < hi, (name, lo, hi)
        assert 68.0 < lo and hi < 74.0, (name, lo, hi)
        assert c['fill_density_required_to_reach_profile_kg_m3'] > 1100.0, name
    tot = led['verdict']['composed_total_body_mass_kg']
    assert 69.0 < tot < 73.0, tot
    assert led['verdict']['answer'].startswith('NO.')
    assert abs(led['cases']['voxel_10mm']['implied_total_body_mass_kg'] - tot) < 0.5, 'partition sensitive'
    assert led['verdict']['profile_bmi'] > led['verdict']['implied_bmi'] + 1.5

    # the composed body is densitometrically coherent and the profile is not
    dm = led['densitometry']
    assert 0.0 < dm['composed_body']['siri_fat_percent'] < 35.0, dm['composed_body']
    assert dm['profile_body']['siri_fat_percent_with_gas'] < 0.0, dm['profile_body']
    assert dm['lung_gas_m3'] > 0.0
    assert led['body_fat']['fat_fraction_lower_bound'] < dm['composed_body']['siri_fat_percent'] / 100.0

    # regionalisations consume the measured histograms with nothing left over and nothing unmet
    for key, block, vk, ck in (('depth', reg['by_depth_below_skin'], 'volume_m3', 'bins'),
                               ('system', reg['by_nearest_entity_system'], 'void_volume_m3', 'systems')):
        assert block['unmet_demand_m3'] == {}, (key, block['unmet_demand_m3'])
        tv = sum(b[vk] for b in block[ck])
        assert abs(tv - void) < 1e-9, (key, tv, void)
        for b in block[ck]:
            assert abs(sum(b['fractions'].values()) + b['unallocated_m3'] / b[vk] - 1.0) < 1e-9, (key, b)
        assert abs(sum(b['mass_kg'] for b in block[ck])
                   - led['cases']['voxel_8mm']['composed_fill_mass_kg']) < 1e-6, key
    # regionalisation earns its keep: the adipose fraction must actually vary
    fr = [b['adipose_volume_fraction'] for b in reg['by_depth_below_skin']['bins']]
    assert max(fr) - min(fr) > 0.9, fr
    dens = [b['density_kg_m3'] for b in reg['by_depth_below_skin']['bins']]
    assert max(dens) - min(dens) > 100.0, dens
    # the shallowest bin is dermis, not fat, and the deepest is visceral fat
    assert reg['by_depth_below_skin']['bins'][0]['fractions'].get('dermis_and_epidermis', 0) > 0.99
    assert reg['by_depth_below_skin']['bins'][-1]['fractions'].get('visceral_adipose', 0) > 0.99

    # the honest empty cells are present and named
    t = doc['tier_census']
    assert t['cells_by_tier'].get('absent', 0) >= 10, t['cells_by_tier']
    assert t['cells_by_tier'].get('assumed', 0) >= 1, t['cells_by_tier']
    assert any('loose_areolar' in c for c in t['absent_cells']), t['absent_cells']
    assert t['identifiers_resolved'] == t['sources_cited'] == len(SOURCES)
    assert len(doc['not_established']) >= 6

    # writing goes to a temporary root, NEVER the repo
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        d, _, man = emit(Path(tmp) / 'interstitial-composition-prior-v1')
        for name, digest in man['artifacts_sha256'].items():
            assert hashlib.sha256((d / name).read_bytes()).hexdigest() == digest, name
        assert man['canonical_assets_modified'] is False
        assert set(man['inputs_sha256']) == {str(p.relative_to(ROOT)) for p in INPUTS}
    for p in INPUTS:
        assert p.exists(), p

    report = {
        'status': 'passed',
        'assertions': 'source completeness and identifier resolution, composition closure on the measured '
                      'void, positive closure residual, ledger closure at two partitions x three densities, '
                      'densitometric coherence, both regionalisations consuming their measured histogram, '
                      'tier census, temporary-root write with digest match',
        'sources': len(SOURCES), 'identifiers_resolved': t['identifiers_resolved'],
        'cells_by_tier': t['cells_by_tier'],
        'void_m3': void,
        'adipose_volume_fraction_of_void': comp['adipose_volume_fraction_of_void'],
        'composed_fill_density_kg_m3': led['cases']['voxel_8mm']['composed_fill_density_kg_m3'],
        'implied_total_body_mass_kg': tot,
        'profile_mass_kg': 77.1107029,
        'canonical_assets_modified': False,
    }
    print(json.dumps(report, indent=2))
    return report


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--output', type=Path)
    ap.add_argument('--summary', action='store_true', help='print the composition table and the ledger')
    ap.add_argument('--verify-identifiers', action='store_true', help='re-resolve every DOI and PMID (network)')
    a = ap.parse_args()
    if not (a.self_test or a.output or a.summary or a.verify_identifiers):
        ap.error('Select --self-test, --summary, --verify-identifiers or --output')
    if a.self_test:
        self_test()
    if a.verify_identifiers:
        v = verify_identifiers()
        print(json.dumps(v, indent=2))
        if v['failures']:
            sys.exit(1)
    if a.summary:
        doc = build()
        print(json.dumps({'composition': doc['composition'], 'mass_ledger': doc['mass_ledger'],
                          'regional': doc['regional'], 'tier_census': doc['tier_census'],
                          'not_established': doc['not_established']}, indent=2))
    if a.output:
        d, doc, man = emit(a.output)
        print(json.dumps({'output_dir': str(d), 'artifacts': man['artifacts_sha256'],
                          'headline': man['headline'], 'cells_by_tier': man['tier_census']['cells_by_tier']},
                         indent=2))
