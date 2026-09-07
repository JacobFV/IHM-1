#!/usr/bin/env python3
"""Candidate literature tissue material table, an audit of the canonical one, and a mass ledger.

    .venv/bin/python scripts/candidate_tissue_materials.py --self-test
    .venv/bin/python scripts/candidate_tissue_materials.py --output data/derived/tissue-material-candidate-v1

Nothing here is canonical. `data/derived/canonical/mechanics.json` is read only and
is not modified; this script writes a candidate table beside it so the two can be
diffed before anything is promoted.

Measured here (arithmetic over shipped assets):
  audit      per-role value/basis/source census of the canonical material block,
             with the volume and mass fraction resting on each basis tier
  ledger     envelope, exclusive-occupancy and void volumes against the profile
             mass, with the density each allocation implies

Transferred here (published values, none measured on this specimen):
  sources    one record per publication, with method, cohort and identifier
  materials  per canonical `role`, plus per-entity overrides where a source is
             specific, each field tagged measured / transferred / assumed

Not established here: no fit to this specimen, no in-silico validation, no claim
that a transferred donor cohort represents the atlas subject, and no probability
interpretation of any range. `prior_range` stays a sensitivity envelope.
"""
from pathlib import Path
import argparse, hashlib, json, sys

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / 'data/derived/canonical/mechanics.json'
PROFILE = ROOT / 'data/derived/canonical/profile.json'
UNMODELLED = ROOT / 'data/derived/unmodelled-volume/summary.json'
PARTITION = ROOT / 'data/derived/material-domains/whole-body-0.01m/manifest.json'
KHORSHIDI = ROOT / 'data/measurements/biomechanics/khorshidi_2024.json'
INPUTS = (CANONICAL, PROFILE, UNMODELLED, PARTITION, KHORSHIDI)

# Tier vocabulary, strictest first. `measured` is reserved for a value measured on
# this specimen and is therefore unreachable for every literature row here.
TIERS = ('measured', 'transferred', 'assumed')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- audit

def audit(mech):
    """Per-role census of the canonical material block."""
    fields = ('young_modulus', 'poisson_ratio', 'shear_modulus', 'lame_lambda', 'density')
    roles, totals = {}, {'volume_m3': 0., 'mass_kg': 0.}
    basis_volume, basis_mass, source_volume = {}, {}, {}
    for e in mech['entities']:
        role = e['role']
        r = roles.setdefault(role, {'entities': 0, 'volume_m3': 0., 'mass_kg': 0.,
                                    'constitutive': set(), 'fields': {}})
        r['entities'] += 1
        v = e.get('material_volume_m3') or 0.
        m = e.get('mass_kg') or 0.
        r['volume_m3'] += v
        r['mass_kg'] += m
        totals['volume_m3'] += v
        totals['mass_kg'] += m
        r['constitutive'].add(e.get('constitutive'))
        material = e.get('material') or {}
        for f in fields:
            b = material.get(f)
            slot = r['fields'].setdefault(f, {'values': set(), 'basis': set(), 'sources': set(),
                                              'prior_range': set(), 'absent': 0})
            if b is None:
                slot['absent'] += 1
                if f == 'young_modulus':
                    basis_volume['absent'] = basis_volume.get('absent', 0.) + v
                    basis_mass['absent'] = basis_mass.get('absent', 0.) + m
                    source_volume['absent'] = source_volume.get('absent', 0.) + v
                continue
            slot['values'].add(b.get('value'))
            slot['basis'].add(b.get('basis'))
            slot['sources'].add(tuple(b.get('sources') or ()))
            pr = b.get('prior_range')
            slot['prior_range'].add(tuple(pr) if pr else None)
            if f == 'young_modulus':
                k = b.get('basis')
                basis_volume[k] = basis_volume.get(k, 0.) + v
                basis_mass[k] = basis_mass.get(k, 0.) + m
                sk = ','.join(b.get('sources') or ()) or 'no_source'
                source_volume[sk] = source_volume.get(sk, 0.) + v
    out = {}
    for role, r in sorted(roles.items()):
        fo = {}
        for f, slot in r['fields'].items():
            vals = sorted(v for v in slot['values'] if v is not None)
            fo[f] = {'distinct_values': vals, 'basis': sorted(slot['basis']),
                     'sources': sorted(','.join(s) or 'no_source' for s in slot['sources']),
                     'prior_ranges': sorted((list(p) if p else None for p in slot['prior_range']), key=str),
                     'entities_without_field': slot['absent']}
        out[role] = {'entities': r['entities'], 'volume_m3': r['volume_m3'], 'mass_kg': r['mass_kg'],
                     'constitutive': sorted(x for x in r['constitutive'] if x), 'fields': fo}
    tv, tm = totals['volume_m3'], totals['mass_kg']
    return {'roles': out, 'totals': totals,
            'young_modulus_basis_share': {k: {'volume_m3': v, 'volume_fraction': v / tv,
                                              'mass_kg': basis_mass.get(k, 0.),
                                              'mass_fraction': basis_mass.get(k, 0.) / tm}
                                          for k, v in sorted(basis_volume.items())},
            'young_modulus_source_share': {k: {'volume_m3': v, 'volume_fraction': v / tv}
                                           for k, v in sorted(source_volume.items())},
            'distinct_literature_sources_cited': sorted(mech['sources']),
            'note': 'basis strings are the canonical vocabulary (source_informed_prior, assumed_prior, '
                    'derived_from_E_nu); `absent` marks a field the entity has no entry for at all'}


# ---------------------------------------------------------------- sources
# One record per publication. `verification` says how the numbers below were
# read: `abstract` = NCBI efetch abstract text, `fulltext` = the article PDF or
# an open full text, `secondary` = a number quoted by another paper that was
# read, `table` = a published reference table. Nothing here was measured on the
# BodyParts3D reference specimen.
SOURCES = {
    'morrow2010': {
        'citation': 'Morrow DA, Haut Donahue TL, Odegard GM, Kaufman KR. Transversely isotropic tensile '
                    'material properties of skeletal muscle tissue. J Mech Behav Biomed Mater '
                    '2010;3(1):124-129.',
        'doi': '10.1016/j.jmbbm.2009.03.004', 'pmid': '19878911',
        'species': 'rabbit (New Zealand White)', 'state': 'fresh-frozen, thawed; aponeurosis dissected off',
        'method': 'quasi-static uniaxial tension, longitudinal extension / transverse extension / '
                  'longitudinal shear',
        'n': '18 extensor digitorum longus muscles from 9 rabbits',
        'ages': 'not reported', 'verification': 'fulltext (author PDF, Table 1)'},
    'quapp1998': {
        'citation': 'Quapp KM, Weiss JA. Material characterization of human medial collateral ligament. '
                    'J Biomech Eng 1998;120(6):757-763.',
        'doi': '10.1115/1.2834890', 'pmid': '10412460',
        'species': 'human', 'state': 'cadaveric',
        'method': 'uniaxial tension of punched specimens along and transverse to the collagen fibre '
                  'direction, optical strain analysis; three constitutive models fitted by nonlinear '
                  'regression (coefficients not printed in the abstract)',
        'n': '10 medial collateral ligaments', 'ages': 'not reported in abstract',
        'verification': 'abstract'},
    'maganaris2002': {
        'citation': 'Maganaris CN, Paul JP. Tensile properties of the in vivo human gastrocnemius tendon. '
                    'J Biomech 2002;35(12):1639-1646.',
        'doi': '10.1016/s0021-9290(02)00240-3', 'pmid': '12445617',
        'species': 'human', 'state': 'in vivo',
        'method': 'real-time ultrasonography of tendon elongation under isometric plantarflexion; force '
                  'from joint moment and in vivo moment arm; 87.5 N to 875 N loading ramp',
        'n': '6 men', 'ages': 'not reported in abstract', 'verification': 'abstract'},
    'maganaris1999': {
        'citation': 'Maganaris CN, Paul JP. In vivo human tendon mechanical properties. '
                    'J Physiol 1999;521 Pt 1:307-313.',
        'doi': '10.1111/j.1469-7793.1999.00307.x', 'pmid': '10562354',
        'species': 'human', 'state': 'in vivo',
        'method': 'ultrasonography of tibialis anterior tendon under isometric dorsiflexion',
        'n': 'not re-read here', 'ages': 'not re-read here',
        'verification': 'inherited from the canonical `tendon` source record, not re-read',
        'note': 'this is the single tendon source the canonical table already cites'},
    'hansen2006': {
        'citation': 'Hansen P, Bojsen-Moller J, Aagaard P, Kjaer M, Magnusson SP. Mechanical properties of '
                    'the human patellar tendon, in vivo. Clin Biomech 2006;21(1):54-58.',
        'doi': '10.1016/j.clinbiomech.2005.07.008', 'pmid': '16183183',
        'species': 'human', 'state': 'in vivo',
        'method': 'simultaneous ultrasonography of tibial and patellar movement during maximal 10 s ramp '
                  'isometric knee extension; reproducibility study',
        'n': 'not stated in abstract', 'ages': 'not stated in abstract', 'verification': 'abstract'},
    'chakouch2015': {
        'citation': 'Chakouch MK, Charleux F, Bensamoun SF. Quantifying the elastic property of nine thigh '
                    'muscles using magnetic resonance elastography. PLoS One 2015;10(9):e0138873.',
        'doi': '10.1371/journal.pone.0138873', 'pmid': '26397730',
        'species': 'human', 'state': 'in vivo, at rest (passive)',
        'method': '1.5 T MR elastography, five protocols; shear modulus of nine thigh muscles and of '
                  'subcutaneous adipose tissue',
        'n': '29 healthy volunteers', 'ages': 'mean 26 +/- 3.41 y', 'verification': 'abstract',
        'note': 'MRE shear modulus is a dynamic modulus at the driver frequency, not a quasi-static one; '
                'it is not directly the mu of a quasi-static neo-Hookean fit'},
    'chakouch2016': {
        'citation': 'Chakouch MK, Pouletaut P, Charleux F, Bensamoun SF. Viscoelastic shear properties of '
                    'in vivo thigh muscles measured by MR elastography. J Magn Reson Imaging '
                    '2016;43(6):1423-1433.',
        'doi': '10.1002/jmri.25105', 'pmid': '26605873',
        'species': 'human', 'state': 'in vivo, passive',
        'method': 'multifrequency MRE at 70/90/110 Hz, complex modulus fitted with four rheological '
                  'models; Zener best fit',
        'n': '5 volunteers, 4 muscles', 'ages': 'not stated in abstract', 'verification': 'abstract'},
    'gennisson2010': {
        'citation': 'Gennisson JL, Deffieux T, Mace E, Montaldo G, Fink M, Tanter M. Viscoelastic and '
                    'anisotropic mechanical properties of in vivo muscle tissue assessed by supersonic '
                    'shear imaging. Ultrasound Med Biol 2010;36(5):789-801.',
        'doi': '10.1016/j.ultrasmedbio.2010.02.013', 'pmid': '20420970',
        'species': 'human', 'state': 'in vivo',
        'method': 'supersonic shear imaging with the probe tilted relative to the fibre direction; Voigt '
                  'viscoelastic model',
        'n': 'not read', 'ages': 'not read',
        'verification': 'abstract read for method only; no numeric parallel-vs-perpendicular split was '
                        'retrieved, so no value is carried from this source'},
    'vanloocke2006': {
        'citation': 'Van Loocke M, Lyons CG, Simms CK. A validated model of passive muscle in '
                    'compression. J Biomech 2006;39(16):2999-3009.',
        'doi': '10.1016/j.jbiomech.2005.10.016', 'pmid': '16313914',
        'species': 'porcine', 'state': 'fresh and aged post-mortem',
        'method': 'quasi-static unconstrained uniaxial compression to 30% strain at 0/30/45/60/90 deg to '
                  'the fibre direction; transversely isotropic hyperelastic and strain-dependent-Young '
                  'modulus (SYM) models fitted',
        'n': 'not stated in abstract', 'ages': 'not stated in abstract',
        'verification': 'abstract for method; the two numeric moduli below are quoted from the Morrow '
                        '2010 discussion (secondary), not read in Van Loocke itself'},
    'takaza2013': {
        'citation': 'Takaza M, Moerman KM, Gindre J, Lyons G, Simms CK. The anisotropic mechanical '
                    'behaviour of passive skeletal muscle tissue subjected to large tensile strain. '
                    'J Mech Behav Biomed Mater 2013;17:209-220.',
        'doi': '10.1016/j.jmbbm.2012.09.001', 'pmid': '23127635',
        'species': 'porcine (longissimus dorsi)', 'state': 'ex vivo',
        'method': 'uniaxial tension along, perpendicular and at 30/45/60 deg to the fibre direction',
        'n': 'not stated in abstract', 'ages': 'not stated in abstract', 'verification': 'abstract'},
    'ward2005': {
        'citation': 'Ward SR, Lieber RL. Density and hydration of fresh and fixed human skeletal muscle. '
                    'J Biomech 2005;38(11):2317-2320.',
        'doi': '10.1016/j.jbiomech.2004.10.001', 'pmid': '16154420',
        'species': 'human', 'state': 'cadaveric, formaldehyde fixed then PBS rehydrated',
        'method': 'volume, water content and mass of fixed muscle samples rehydrated for 0-30 h',
        'n': '54 samples per fixation method; water content in living muscle n=4',
        'ages': 'not stated in abstract', 'verification': 'abstract'},
    'guenard1992': {
        'citation': 'Guenard H, Diallo MH, Laurent F, Vergeret J. Lung density and lung mass in '
                    'emphysema. Chest 1992;102(1):198-203.',
        'doi': '10.1378/chest.102.1.198', 'pmid': '1623752',
        'species': 'human', 'state': 'in vivo, whole-lung CT at functional residual capacity',
        'method': 'mean lung density times radiologic lung volume from CT',
        'n': '16 healthy subjects (plus 24 emphysema patients, not used here)',
        'ages': 'not stated in abstract', 'verification': 'abstract'},
    'icru44_nist': {
        'citation': 'Hubbell JH, Seltzer SM. Tables of X-Ray Mass Attenuation Coefficients and Mass '
                    'Energy-Absorption Coefficients, Table 2. NIST Standard Reference Database 126. '
                    'Densities are those of ICRU Report 44 (1989), Tissue Substitutes in Radiation '
                    'Dosimetry and Measurement.',
        'doi': '10.18434/T4D01F', 'pmid': None,
        'species': 'human', 'state': 'reference composition, not a specimen',
        'method': 'standardised reference tissue compositions and mass densities',
        'n': 'reference table', 'ages': 'adult reference',
        'verification': 'table (physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html)'},
    'icrp89': {
        'citation': 'ICRP. Basic anatomical and physiological data for use in radiological protection: '
                    'reference values. ICRP Publication 89. Ann ICRP 2002;32(3-4):1-277.',
        'doi': '10.1016/S0146-6453(03)00002-2', 'pmid': None,
        'pmid_note': 'ICRP 89 is not indexed in PubMed. PMID 14527029, previously recorded here, '
                     'is Boecker, Radiat Prot Dosimetry 2003;105(1-4):571-4 — a four-page conference '
                     'summary about ICRP 89, not the publication itself. Verified against esummary.',
        'species': 'human', 'state': 'reference values, not a specimen',
        'method': 'Task Group synthesis over Western European and North American autopsy and imaging '
                  'series; Table 2.8 organ masses, Table 2.9 height/mass/surface area, Table 2.20 '
                  'skeletal densities, paragraphs 434 and 514',
        'n': 'reference adult male, 73 kg, 176 cm, 1.90 m2 body surface',
        'ages': 'adult 20-50 y', 'verification': 'fulltext (Tables 2.8, 2.9, 2.20)'},
    'khorshidi2024': {
        'citation': 'Characterisation of human penile tissue properties using experimental testing '
                    'combined with multi-target inverse finite element modelling. Acta Biomater 2024.',
        'doi': '10.1016/j.actbio.2024.06.035', 'pmid': None,
        'species': 'human', 'state': 'fresh-frozen ex vivo, tested in PBS at 37 C',
        'method': 'multi-target inverse finite-element fit to whole-organ and specimen tests',
        'n': '3 donors, 10 whole-organ segments, 28 tunica and 32 cavernosum specimens',
        'ages': '75-89 y',
        'verification': 'already in this repo as data/measurements/biomechanics/khorshidi_2024.json'},
    'organs_pmid33176223': {
        'citation': 'Micromechanical poroelastic and viscoelastic properties of ex-vivo soft tissues. '
                    'J Biomech 2020;113:110090.',
        'doi': '10.1016/j.jbiomech.2020.110090', 'pmid': '33176223',
        'species': 'not human-confirmed here', 'state': 'ex vivo',
        'method': 'microindentation', 'n': 'not read here', 'ages': 'not read here',
        'verification': 'identifier confirmed against PubMed; this is the single source the canonical '
                        'table leans on for 29.4 L of soft tissue'},
}

SOURCES.update({
    'budday2017': {
        'citation': 'Budday S, Sommer G, Birkl C, Langkammer C, Haybaeck J, Kohnert J, Bauer M, Paulsen F, '
                    'Steinmann P, Kuhl E, Holzapfel GA. Mechanical characterization of human brain tissue. '
                    'Acta Biomater 2017;48:319-340.',
        'doi': '10.1016/j.actbio.2016.10.036', 'pmid': '27989920',
        'species': 'human', 'state': 'cadaveric, post-mortem interval under 24 h',
        'method': 'simple shear to gamma=0.2, unconfined compression to 10% strain and uniaxial tension to '
                  '10% strain with relaxation holds; one-term modified Ogden fitted simultaneously across '
                  'all three modes',
        'n': '10 brains, four regions each (corpus callosum, corona radiata, basal ganglia, cortex)',
        'ages': '54-81 y, mean 66 y', 'verification': 'fulltext (Tables 4 and 7)'},
    'gao2015': {
        'citation': 'Gao H, Li WG, Cai L, Berry C, Luo XY. Parameter estimation in a Holzapfel-Ogden law '
                    'for healthy myocardium. J Eng Math 2015;95(1):231-248.',
        'doi': '10.1007/s10665-014-9740-3', 'pmid': '26663931',
        'species': 'human', 'state': 'in vivo',
        'method': 'inverse estimation of the eight Holzapfel-Ogden orthotropic constants from cine MRI '
                  'and measured cavity pressure',
        'n': '3 healthy volunteers', 'ages': '22, 28, 31 y', 'verification': 'fulltext (Table 5)'},
    'sommer2015': {
        'citation': 'Sommer G, Schriefl AJ, Andrae M, Sacherer M, Viertler C, Wolinski H, Holzapfel GA. '
                    'Biomechanical properties and microstructure of human ventricular myocardium. '
                    'Acta Biomater 2015;24:172-192.',
        'doi': '10.1016/j.actbio.2015.06.031', 'pmid': '26141152',
        'species': 'human', 'state': 'ex vivo passive left ventricle',
        'method': 'planar biaxial extension and triaxial simple shear with second-harmonic-generation '
                  'microscopy of fibre and sheet orientation',
        'n': 'not read', 'ages': 'not read',
        'verification': 'abstract only; the orthotropic constants are in paywalled tables and were NOT '
                        'retrieved, so no number is carried from this source'},
    'rouviere2006': {
        'citation': 'Rouviere O, Yin M, Dresner MA, Rossman PJ, Burgart LJ, Fidler JL, Ehman RL. '
                    'MR elastography of the liver: preliminary results. Radiology 2006;240(2):440-448.',
        'doi': '10.1148/radiol.2402050606', 'pmid': '16864671',
        'species': 'human', 'state': 'in vivo',
        'method': 'MR elastography of the liver, shear stiffness',
        'n': 'healthy volunteers plus fibrosis patients', 'ages': 'not read',
        'verification': 'abstract'},
    'lee2013': {
        'citation': 'Lee DH, Lee JM, Han JK, Choi BI. MR elastography of healthy liver parenchyma: normal '
                    'value and reliability of the liver stiffness value measurement. J Magn Reson Imaging '
                    '2013;38(5):1215-1223.',
        'doi': '10.1002/jmri.23958', 'pmid': '23281116',
        'species': 'human', 'state': 'in vivo',
        'method': 'MR elastography, three region-of-interest methods',
        'n': '49 living liver donors with normal laboratory results', 'ages': 'not read',
        'verification': 'abstract'},
    'arda2011': {
        'citation': 'Arda K, Ciledag N, Aktas E, Aribas BK, Kose K. Quantitative assessment of normal '
                    'soft-tissue elasticity using shear-wave ultrasound elastography. AJR Am J Roentgenol '
                    '2011;197(3):532-536.',
        'doi': '10.2214/AJR.10.5449', 'pmid': '21862792',
        'species': 'human', 'state': 'in vivo',
        'method': 'shear-wave ultrasound elastography of normal organs',
        'n': '127 healthy volunteers', 'ages': 'mean 37.72 +/- 9.11 y', 'verification': 'abstract'},
    'pawlus2016': {
        'citation': 'Pawlus A, Sokolowska-Dabek D, Szymanska K, Inglot MS, Zaleska-Dorobisz U. Shear wave '
                    'elastography of the spleen: evaluation of spleen stiffness in healthy volunteers. '
                    'Abdom Radiol 2016;41(11):2169-2174.',
        'doi': '10.1007/s00261-016-0834-4', 'pmid': '27389244',
        'species': 'human', 'state': 'in vivo',
        'method': 'shear-wave elastography of the spleen',
        'n': '59 healthy volunteers', 'ages': 'not read', 'verification': 'fulltext'},
    'nelson2026': {
        'citation': 'Nelson TM, Eskandari M. Human lung parenchyma: tensile mechanics and the effects of '
                    'smoking. J R Soc Interface 2026;23(238):20250721.',
        'doi': '10.1098/rsif.2025.0721', 'pmid': '42120048',
        'species': 'human', 'state': 'ex vivo donor lungs',
        'method': 'uniaxial tension of parenchymal strips; final (large-strain) stiffness modulus',
        'n': '8 donor lungs', 'ages': 'not read', 'verification': 'abstract'},
    'tang2021': {
        'citation': 'Tang X, Wang L, Guo R, Huang S, Tang Y, Qiu L. Preliminary study on the influencing '
                    'factors of shear wave elastography for peripheral nerves in healthy population. '
                    'Sci Rep 2021;11:5582.',
        'doi': '10.1038/s41598-021-84900-8', 'pmid': '33692411',
        'species': 'human', 'state': 'in vivo',
        'method': 'shear-wave elastography; the authors report shear wave VELOCITY and explicitly decline '
                  'to report a Young modulus for these anisotropic nerves',
        'n': '105 healthy volunteers', 'ages': 'not read', 'verification': 'fulltext'},
    'teng2015': {
        'citation': 'Teng Z, Feng J, Zhang Y, Sutcliffe MPF, Huang Y, Brown AJ, Jing Z, Lu Q, Gillard JH. '
                    'Layer- and direction-specific material properties, extreme extensibility and ultimate '
                    'material strength of human abdominal aorta and aneurysm. Ann Biomed Eng '
                    '2015;43(11):2745-2759.',
        'doi': '10.1007/s10439-015-1323-6', 'pmid': '25905688',
        'species': 'human', 'state': 'fresh surgical samples, cryopreserved',
        'method': 'uniaxial tension at 0.01 mm/s, circumferential and axial, layer separated; modified '
                  'Mooney-Rivlin fitted (NOT Holzapfel-Gasser-Ogden)',
        'n': '8 normal aorta donors (plus 11 aneurysm patients, not used here)',
        'ages': '34.1 +/- 7.8 y for the normal group', 'verification': 'fulltext'},
    'holzapfel2005': {
        'citation': 'Holzapfel GA, Sommer G, Gasser CT, Regitnig P. Determination of layer-specific '
                    'mechanical properties of human coronary arteries with nonatherosclerotic intimal '
                    'thickening and related constitutive modeling. Am J Physiol Heart Circ Physiol '
                    '2005;289(5):H2048-H2058.',
        'doi': '10.1152/ajpheart.00934.2004', 'pmid': None,
        'species': 'human', 'state': 'autopsy, non-stenotic left anterior descending coronary arteries',
        'method': 'cyclic quasi-static uniaxial tension of separated intima, media and adventitia strips',
        'n': 'reported as 13 arteries / 78 strips, not verified here', 'ages': 'reported as 71.5 +/- 7.3 y',
        'verification': 'NOT retrieved; full text was inaccessible, so no number is carried from this '
                        'source even though it is the standard layer-specific human artery reference'},
    'sommer2013': {
        'citation': 'Sommer G, Eder M, Kovacs L, Pathak H, Bonitz L, Mueller C, Regitnig P, Holzapfel GA. '
                    'Multiaxial mechanical properties and constitutive modeling of human adipose tissue: '
                    'a basis for preoperative simulations in plastic and reconstructive surgery. '
                    'Acta Biomater 2013;9(11):9036-9048.',
        'doi': '10.1016/j.actbio.2013.06.011', 'pmid': '23811521',
        'species': 'human', 'state': 'fresh abdominal adipose tissue',
        'method': 'biaxial tension and triaxial shear, quasi-static and dynamic; a hyperelastic model was '
                  'fitted',
        'n': 'not read', 'ages': 'not read',
        'verification': 'abstract only; the fitted coefficients are not in the abstract and were NOT '
                        'retrieved, so no number is carried from this source'},
    'alkhouli2013': {
        'citation': 'Alkhouli N, Mansfield J, Green E, Bell J, Knight B, Liversedge N, Tham JC, Welbourn R, '
                    'Shore AC, Kos K, Winlove CP. The mechanical properties of human adipose tissues and '
                    'their relationships to the structure and composition of the extracellular matrix. '
                    'Am J Physiol Endocrinol Metab 2013;305(12):E1427-E1435.',
        'doi': '10.1152/ajpendo.00111.2013', 'pmid': '24105412',
        'species': 'human', 'state': 'fresh surgical samples',
        'method': 'uniaxial tension to 30% strain with stress relaxation; initial and final tangent moduli',
        'n': '44 subjects; 19 paired subcutaneous/omental comparisons', 'ages': 'not read',
        'verification': 'abstract'},
    'niannaidh2012': {
        'citation': 'Ni Annaidh A, Bruyere K, Destrade M, Gilchrist MD, Ottenio M. Characterization of the '
                    'anisotropic mechanical properties of excised human skin. J Mech Behav Biomed Mater '
                    '2012;5(1):139-148.',
        'doi': '10.1016/j.jmbbm.2011.08.016', 'pmid': '22100088',
        'species': 'human', 'state': 'excised back skin, refrigerated at 4 C, tested 2-9 days post mortem, '
                                     'no preconditioning',
        'method': 'uniaxial tension, ASTM D412 dogbone, strain rate 0.012 /s; mean specimen thickness '
                  '2.56 +/- 0.39 mm; no strain-energy function was fitted despite the title',
        'n': '56 samples from 7 donors (3 male, 4 female)', 'ages': '81-97 y, mean 89 +/- 6 y',
        'verification': 'fulltext'},
    'feng2022': {
        'citation': 'Feng X, Li GY, Ramier A, Eltony AM, Yun SH. In vivo stiffness measurement of '
                    'epidermis, dermis, and hypodermis using broadband Rayleigh-wave optical coherence '
                    'elastography. Acta Biomater 2022;146:295-305.',
        'doi': '10.1016/j.actbio.2022.04.030', 'pmid': '35470076',
        'species': 'human', 'state': 'in vivo, volar forearm',
        'method': 'broadband leaky Rayleigh-wave optical coherence elastography, 0.1-10 kHz; layer moduli '
                  'from surface-wave dispersion. The authors ASSUME Poisson ratio 0.4999 and density '
                  '1 g/cm3 per layer; neither is measured by this method',
        'n': 'healthy volunteers, count not read', 'ages': 'not read', 'verification': 'fulltext'},
    'pailler2008': {
        'citation': 'Pailler-Mattei C, Bec S, Zahouani H. In vivo measurements of the elastic mechanical '
                    'properties of human skin by indentation tests. Med Eng Phys 2008;30(5):599-606.',
        'doi': '10.1016/j.medengphy.2007.06.011', 'pmid': '17869160',
        'species': 'human', 'state': 'in vivo',
        'method': 'indentation with a two-layer elastic model accounting for subcutaneous layers',
        'n': 'not verified', 'ages': 'not verified',
        'verification': 'abstract read for method; the frequently quoted 9.5 +/- 2 kPa was NOT confirmed '
                        'against the abstract, so no number is carried from this source'},
    'reilly1975': {
        'citation': 'Reilly DT, Burstein AH. The elastic and ultimate properties of compact bone tissue. '
                    'J Biomech 1975;8(6):393-405.',
        'doi': '10.1016/0021-9290(75)90075-5', 'pmid': '1206042',
        'species': 'human', 'state': 'wet femoral compact bone',
        'method': 'tension and compression along and transverse to the osteon axis; transversely isotropic '
                  'elastic constants',
        'n': 'not verified here', 'ages': 'not verified here',
        'verification': 'secondary; the primary text was not retrievable, and the constants below were '
                        'read from a teaching table that cites Reilly and Burstein 1975 explicitly'},
    'bayraktar2004': {
        'citation': 'Bayraktar HH, Morgan EF, Niebur GL, Morris GE, Wong EK, Keaveny TM. Comparison of the '
                    'elastic and yield properties of human femoral trabecular and cortical bone tissue. '
                    'J Biomech 2004;37(1):27-35.',
        'doi': '10.1016/s0021-9290(03)00257-4', 'pmid': '14672565',
        'species': 'human', 'state': 'cadaveric femur',
        'method': 'micro-CT based finite-element back-calculation of trabecular TISSUE modulus, and '
                  'direct tensile test of cortical specimens',
        'n': '11 femoral neck donors (trabecular), 34 diaphyseal donors (cortical)',
        'ages': 'not read', 'verification': 'abstract'},
    'morgan2003': {
        'citation': 'Morgan EF, Bayraktar HH, Keaveny TM. Trabecular bone modulus-density relationships '
                    'depend on anatomic site. J Biomech 2003;36(7):897-904.',
        'doi': '10.1016/s0021-9290(03)00071-x', 'pmid': '12757797',
        'species': 'human', 'state': 'cadaveric',
        'method': 'apparent modulus and apparent density over 142 specimens from vertebra, proximal tibia, '
                  'greater trochanter and femoral neck',
        'n': '142 specimens', 'ages': 'not read',
        'verification': 'abstract; the site-specific power-law coefficients are in the paywalled tables '
                        'and were NOT retrieved. The abstract conclusion carried here is qualitative: '
                        'there is no universal modulus-density relationship for on-axis loading'},
    'vertebral2021': {
        'citation': 'Density and mechanical properties of vertebral trabecular bone - a review. '
                    'JOR Spine 2021;4(4):e1176.',
        'doi': '10.1002/jsp2.1176', 'pmid': '35005442',
        'species': 'human', 'state': 'review of cadaveric series',
        'method': 'literature review of vertebral trabecular apparent density and apparent modulus',
        'n': 'review', 'ages': 'review', 'verification': 'abstract'},
    'athanasiou1995': {
        'citation': 'Athanasiou KA, Niederauer GG, Schenck RC Jr. Biomechanical topography of human ankle '
                    'cartilage. Ann Biomed Eng 1995;23(5):697-704.',
        'doi': '10.1007/BF02584467', 'pmid': '7503470',
        'species': 'human', 'state': 'cadaveric',
        'method': 'automated biphasic creep indentation at 14 sites; aggregate modulus, Poisson ratio and '
                  'permeability solved simultaneously in the linear biphasic (KLM) model',
        'n': '14 ankles (7 pairs)', 'ages': 'not reported in abstract', 'verification': 'abstract'},
    'liu1997': {
        'citation': 'Liu GT, Lavery LA, Schenck RC Jr, Lanctot DR, Zhu CF, Athanasiou KA. Human articular '
                    'cartilage biomechanics of the second metatarsal intermediate cuneiform joint. '
                    'J Foot Ankle Surg 1997;36(5):367-374.',
        'doi': '10.1016/s1067-2516(97)80039-7', 'pmid': '9356916',
        'species': 'human', 'state': 'freshly frozen cadaveric',
        'method': 'automated biphasic creep indentation at 4 sites',
        'n': '14 specimens (7 pairs)', 'ages': 'not reported in abstract', 'verification': 'abstract'},
    'athanasiou1998': {
        'citation': 'Athanasiou KA, Liu GT, Lavery LA, Lanctot DR, Schenck RC Jr. Biomechanical topography '
                    'of human articular cartilage in the first metatarsophalangeal joint. Clin Orthop '
                    'Relat Res 1998;(348):269-281.',
        'doi': None, 'pmid': '9553561',
        'species': 'human', 'state': 'freshly frozen cadaveric',
        'method': 'automated biphasic creep indentation at 8 sites',
        'n': '7 pairs of metatarsophalangeal joints', 'ages': 'not reported in abstract',
        'verification': 'abstract'},
    'kanematsu2015': {
        'citation': 'Kanematsu N, Inaniwa T, Koba Y. Relationship between mass density, electron density '
                    'and elemental composition of body tissues for Monte Carlo simulation in radiation '
                    'treatment planning. arXiv:1508.00226 (2015). Table I reproduces the ICRP Publication '
                    '110 (2009) reference tissue densities.',
        'doi': '10.48550/arXiv.1508.00226', 'pmid': None,
        'species': 'human', 'state': 'reference composition, not a specimen',
        'method': 'reference tissue table', 'n': 'reference', 'ages': 'adult reference',
        'verification': 'fulltext (Table I)'},
})


# ---------------------------------------------------------------- values
# Every row: (value, unit, tier, source_key, note). `tier` is one of
#   measured     measured on THIS specimen                 -- unreachable here
#   transferred  a published value measured on other tissue/subjects
#   derived      arithmetic on a transferred value, stated in the note
#   assumed      an engineering number with no source
#   absent       no defensible value was found; the cell is deliberately empty
# An `absent` row keeps its note. Do not silently substitute a number for it.
A = 'absent'


def row(value, unit, tier, source, note):
    return {'value': value, 'unit': unit, 'tier': tier, 'source': source, 'note': note}


# Per canonical `role`. `hyperelastic` carries a published strain-energy fit when
# one exists; `linear` carries the isotropic reduction the current solver needs.
MATERIALS = {
    'muscle': {
        'tissue': 'passive skeletal muscle, transversely isotropic',
        'linear': {
            'shear_modulus': row(3910., 'Pa', 'transferred', 'chakouch2015',
                                 'rectus femoris at rest, the lowest of the nine thigh muscles; the band '
                                 'over the six bulk muscles is 3.91-4.23 kPa and the three strap muscles '
                                 '(gracilis 6.15, semitendinosus 5.32, sartorius 5.15 kPa) are given as '
                                 'entity overrides. This is a DYNAMIC MRE modulus at the driver frequency, '
                                 'not a quasi-static neo-Hookean mu; it will over-read a quasi-static fit'),
            'young_modulus': row(11730., 'Pa', 'derived', 'chakouch2015',
                                 'E = 3G under incompressibility from the 3.91 kPa MRE shear modulus; '
                                 'the source reports G, not E'),
            'poisson_ratio': row(None, '1', A, None,
                                 'no measurement retrieved; soft tissue is conventionally taken as nearly '
                                 'incompressible, which is an assumption, not a measured value'),
            'density': row(1050., 'kg/m3', 'transferred', 'icru44_nist',
                           'ICRU-44 skeletal muscle. Ward and Lieber measured 1112 +/- 6 kg/m3 in 4% and '
                           '1055 +/- 6 kg/m3 in 37% formaldehyde-FIXED human muscle and note that the '
                           'widely used 1059.7 originates in unfixed rabbit and canine tissue; no fresh '
                           'human value was retrieved'),
        },
        'anisotropy': {
            'status': 'UNRESOLVED, and the two available regimes disagree in sign',
            'tension': row(19.96, '1', 'transferred', 'morrow2010',
                           'RABBIT extensor digitorum longus, aponeurosis removed, fresh-frozen: '
                           'longitudinal 447 +/- 97.7 kPa against transverse 22.4 +/- 14.7 kPa, so along-'
                           'fibre is 20x stiffer in tension. Longitudinal shear 3.87 +/- 3.39 kPa'),
            'compression': row(0.447, '1', 'transferred', 'vanloocke2006',
                               'PORCINE, unconstrained compression at 0.3 strain: longitudinal 2.04 kPa '
                               'against transverse 4.56 kPa, so ACROSS-fibre is 2.2x stiffer, the opposite '
                               'ordering. These two numbers are quoted from the Morrow 2010 discussion, '
                               'not read in Van Loocke itself'),
            'human_value': row(None, '1', A, 'gennisson2010',
                               'no human along-versus-across fibre split was retrieved. Supersonic shear '
                               'imaging and multifrequency MRE both report human muscle anisotropy but '
                               'neither retrieved abstract carried the decomposed numbers'),
            'implication': 'a transversely isotropic law cannot be parameterised from human data today, '
                           'and the cross-species evidence says muscle is bimodular, so one anisotropy '
                           'ratio cannot serve both tension and compression',
        },
        'hyperelastic': row(None, None, A, 'takaza2013',
                            'no verified Ogden, Mooney-Rivlin or Yeoh coefficient set for passive skeletal '
                            'muscle of any species was retrieved. Takaza reports 77 kPa transverse stress '
                            'at stretch 1.1 in porcine longissimus dorsi, a stress point, not a fit'),
    },
    'tendon': {
        'tissue': 'tendon',
        'linear': {
            'young_modulus': row(1.16e9, 'Pa', 'transferred', 'maganaris2002',
                                 'human gastrocnemius tendon in vivo, 1.16 +/- 0.15 GPa, hysteresis '
                                 '18 +/- 3%. Independently corroborated: patellar 1.09 GPa in vivo, '
                                 'tibialis anterior 1.2 GPa in vivo. Three tendons, three cohorts, one '
                                 'band 1.09-1.20 GPa'),
            'poisson_ratio': row(None, '1', A, None, 'no measurement retrieved'),
            'density': row(None, 'kg/m3', A, None,
                           'no human whole-tendon density measurement was found. A bovine flexor value '
                           'and a collagen-microfibril value surfaced but neither is human whole tissue'),
        },
        'hyperelastic': row(None, None, A, None,
                            'no numeric HGO or Ogden coefficient set for human tendon was retrieved'),
    },
    'ligament': {
        'tissue': 'ligament, anisotropic',
        'linear': {
            'young_modulus_along_fibre': row(3.322e8, 'Pa', 'transferred', 'quapp1998',
                                             'human MCL tangent modulus 332.2 +/- 58.3 MPa; tensile '
                                             'strength 38.6 +/- 4.8 MPa; ultimate strain 17.1 +/- 1.5%'),
            'young_modulus_transverse': row(1.10e7, 'Pa', 'transferred', 'quapp1998',
                                            'human MCL transverse tangent modulus 11.0 +/- 3.6 MPa; the '
                                            'transverse curves are linear, the longitudinal ones are not. '
                                            'Along-fibre is 30x the transverse'),
            'poisson_ratio': row(None, '1', A, None, 'no measurement retrieved'),
            'density': row(None, 'kg/m3', A, None, 'no human ligament density measurement was found'),
        },
        'hyperelastic': row(None, None, A, 'quapp1998',
                            'Quapp and Weiss fitted three constitutive models by nonlinear regression but '
                            'the coefficients are not in the abstract and were not retrieved'),
    },
    'rigid_bone': {
        'tissue': 'bone; the atlas geometry is WHOLE bones, so it encloses marrow cavity as well as '
                  'mineralised tissue',
        'linear': {
            'young_modulus_cortical_longitudinal': row(1.70e10, 'Pa', 'transferred', 'reilly1975',
                                                       'human wet femoral compact bone, 17 GPa along the '
                                                       'osteon axis'),
            'young_modulus_cortical_transverse': row(1.15e10, 'Pa', 'transferred', 'reilly1975',
                                                     'human wet femoral compact bone, 11.5 GPa transverse; '
                                                     'cortical bone is transversely isotropic, not '
                                                     'isotropic'),
            'shear_modulus_cortical': row(3.6e9, 'Pa', 'transferred', 'reilly1975',
                                          'longitudinal shear 3.6 GPa, transverse shear 3.3 GPa'),
            'poisson_ratio_cortical': row(0.31, '1', 'transferred', 'reilly1975',
                                          'the teaching table quoting Reilly and Burstein gives 0.58 and '
                                          '0.31 for the two independent ratios of the transversely '
                                          'isotropic model; a second directly read compilation gives an '
                                          'isotropic 0.39. Nothing in the human bone literature retrieved '
                                          'supports 0.45'),
            'young_modulus_tissue_trabecular': row(1.80e10, 'Pa', 'transferred', 'bayraktar2004',
                                                   'human femoral neck trabecular TISSUE modulus '
                                                   '18.0 +/- 2.8 GPa against cortical tissue '
                                                   '19.9 +/- 1.8 GPa from the same laboratory; the '
                                                   'trabecular tissue is only 10% softer than cortical '
                                                   'tissue. This is tissue level, not apparent level'),
            'young_modulus_apparent_trabecular': row(None, 'Pa', A, 'vertebral2021',
                                                     'vertebral apparent modulus spans 0.1-976 MPa across '
                                                     'the literature, narrowing to about 10% spread only '
                                                     'under standardised testing. No single apparent '
                                                     'modulus is defensible; it is a function of apparent '
                                                     'density and site, and Morgan shows there is no '
                                                     'universal modulus-density law'),
            'density_whole_skeleton': row(1300., 'kg/m3', 'transferred', 'icrp89',
                                          'ICRP 89 Table 2.20, whole fresh adult skeleton. THIS is the '
                                          'density that matches whole-bone atlas geometry, because that '
                                          'geometry encloses the marrow cavity'),
            'density_cortical': row(1900., 'kg/m3', 'transferred', 'icrp89',
                                    'ICRP 89 Table 2.20, hydrated adult cortical bone; ICRU-44 gives 1920. '
                                    'Fresh bone free of marrow is 1900-2000 (ICRP 89 para 434)'),
            'density_trabecular_with_marrow': row(1100., 'kg/m3', 'transferred', 'icrp89',
                                                  'ICRP 89 para 434, range 800-1400'),
            'apparent_density_trabecular_vertebral': row(None, 'kg/m3', A, 'vertebral2021',
                                                         'wet apparent density 90-350 kg/m3; a range, not '
                                                         'a point value, and site dependent'),
        },
        'hyperelastic': row(None, None, A, None,
                            'bone is linear elastic over the physiological range; no hyperelastic fit is '
                            'needed or was sought'),
    },
    'cartilage': {
        'tissue': 'articular cartilage, biphasic',
        'linear': {
            'aggregate_modulus': row(1.0e6, 'Pa', 'transferred', 'liu1997',
                                     'human second metatarsal-intermediate cuneiform, H_A 0.99 MPa and '
                                     '1.05 MPa. Corroborated by human ankle 0.92-1.34 MPa and human first '
                                     'metatarsophalangeal 0.63-1.34 MPa, all by the same biphasic creep '
                                     'indentation method. Three human joints, one band around 1 MPa'),
            'poisson_ratio': row(0.08, '1', 'transferred', 'liu1997',
                                 'human, all test sites grouped. The human ankle study independently '
                                 'reports 0.08 as its LARGEST site value. This is the biphasic solid-phase '
                                 'ratio, and it is nowhere near 0.45'),
            'permeability': row(3.05e-15, 'm4/(N.s)', 'transferred', 'liu1997',
                                'human; the ankle study spans 0.80-1.79e-15 and the first '
                                'metatarsophalangeal 1.26-4.56e-15 m4/(N.s)'),
            'young_modulus': row(9.86e5, 'Pa', 'derived', 'liu1997',
                                 'E = H_A (1+nu)(1-2nu)/(1-nu) at H_A = 1.0 MPa and nu = 0.08. At nu = 0.08 '
                                 'the aggregate and Young moduli almost coincide; at nu = 0.45 the same E '
                                 'would give H_A = 3.79 MPa, so the Poisson ratio, not the modulus, is '
                                 'what the canonical table gets wrong'),
            'density': row(None, 'kg/m3', A, 'kanematsu2015',
                           'no cartilage-specific density was retrieved. ICRP 110 pools skin, cartilage '
                           'and spongiosa at 1090 kg/m3, which is a pooled category, not cartilage'),
        },
        'hyperelastic': row(None, None, A, None,
                            'no verified biphasic or hyperelastic parameter SET for human cartilage was '
                            'retrieved; the three sources above give H_A, nu and k, which is a complete '
                            'linear biphasic parameterisation but not a strain-energy function'),
    },
    'adipose': {
        'tissue': 'adipose / subcutaneous fat. NOT a canonical role: no adipose geometry exists in the '
                  'atlas. This row exists because 32.117 L of interior belongs to no entity and the fill '
                  'has to be given a material',
        'linear': {
            'young_modulus_initial': row(1600., 'Pa', 'transferred', 'alkhouli2013',
                                         'human SUBCUTANEOUS adipose, initial tangent modulus '
                                         '1.6 +/- 0.8 kPa; final tangent at 30% strain 11.7 +/- 6.4 kPa. '
                                         'Omental fat is stiffer: 2.9 +/- 1.5 initial, 32 +/- 15.6 final'),
            'young_modulus_final': row(11700., 'Pa', 'transferred', 'alkhouli2013',
                                       'human subcutaneous at 30% strain; the 7x rise from initial to '
                                       'final is exactly why a linear modulus is the wrong description'),
            'shear_modulus_mre': row(3040., 'Pa', 'transferred', 'chakouch2015',
                                     'human subcutaneous adipose in vivo, 3.04 +/- 0.12 kPa, the softest '
                                     'tissue in that MRE series; dynamic modulus caveat as for muscle'),
            'young_modulus_hypodermis_in_vivo': row(15000., 'Pa', 'transferred', 'feng2022',
                                                    'human hypodermis in vivo, about 15 kPa at 0.2-1 kHz'),
            'poisson_ratio': row(None, '1', A, None, 'no measurement retrieved'),
            'density': row(950., 'kg/m3', 'transferred', 'icru44_nist',
                           'ICRU-44 adipose tissue, and ICRP 110 adipose/marrow is also 950. Pure human '
                           'fat at 37 C is 900 kg/m3 (Fidanza, Keys and Anderson 1953, via Kanematsu). '
                           'Adipose TISSUE is about 80% fat in adults (ICRP 89 para 545), which is why '
                           'the tissue density exceeds the fat density'),
        },
        'hyperelastic': row(None, None, A, 'sommer2013',
                            'Sommer 2013 fitted a hyperelastic model to multiaxial human abdominal adipose '
                            'data but the coefficients are not in the abstract and were not retrieved. '
                            'Comley and Fleck is porcine and its numbers could not be verified'),
    },
    'skin_layer': {
        'tissue': 'skin, resolved by layer. The three canonical layers currently share ONE modulus',
        'linear': {
            'young_modulus_epidermis': row(4.0e6, 'Pa', 'transferred', 'feng2022',
                                           'human volar forearm in vivo, about 4 MPa at 4-10 kHz; the high '
                                           'band is what resolves the thin epidermis'),
            'young_modulus_dermis': row(40000., 'Pa', 'transferred', 'feng2022',
                                        'human volar forearm in vivo, about 40 kPa at 0.2-1 kHz'),
            'young_modulus_hypodermis': row(15000., 'Pa', 'transferred', 'feng2022',
                                            'human volar forearm in vivo, about 15 kPa at 0.2-1 kHz'),
            'young_modulus_excised_whole_skin': row(8.33e7, 'Pa', 'transferred', 'niannaidh2012',
                                                    'human excised back skin, 83.3 +/- 34.9 MPa; initial '
                                                    'slope 1.18 +/- 0.88 MPa; UTS 21.6 +/- 8.4 MPa; failure '
                                                    'strain 54 +/- 17%. This is the LARGE-STRAIN collagen-'
                                                    'engaged modulus of elderly excised tissue and is not '
                                                    'the same quantity as the in-vivo small-strain layer '
                                                    'moduli above; the paper own literature table spans '
                                                    '0.26-150 MPa across authors'),
            'poisson_ratio': row(None, '1', A, 'feng2022',
                                 'the OCE analysis ASSUMES 0.4999 and states so; that is an input to their '
                                 'model, not a measurement'),
            'density': row(1100., 'kg/m3', 'transferred', 'icrp89',
                           'ICRP 89 para 514, skin density approximately 1.1 g/cm3, the value that report '
                           'itself uses to convert mass thickness to linear thickness'),
        },
        'hyperelastic': row(None, None, A, 'niannaidh2012',
                            'no verified Ogden fit for human skin was retrieved. Ni Annaidh reports only '
                            'linear-region descriptors despite the constitutive framing'),
    },
    'soft_organ': {
        'tissue': 'soft organs; a single role covering 301 entities from cornea to lung',
        'linear': {
            'shear_modulus_liver': row(2000., 'Pa', 'transferred', 'rouviere2006',
                                       'healthy human liver MR elastography 2.0 +/- 0.3 kPa, '
                                       'independently 2.05 kPa (range 1.54-2.87) in 49 living donors'),
            'shear_modulus_kidney_cortex': row(5000., 'Pa', 'transferred', 'arda2011',
                                               'renal cortex 5.0 +/- 2.9 kPa; renal pelvis 23.6 +/- 5.4'),
            'shear_modulus_pancreas': row(4800., 'Pa', 'transferred', 'arda2011', '4.8 +/- 3 kPa'),
            'shear_modulus_thyroid': row(10970., 'Pa', 'transferred', 'arda2011', '10.97 +/- 3.1 kPa'),
            'shear_modulus_spleen': row(None, 'Pa', A, 'arda2011',
                                        'CONTESTED: 2.9 +/- 1.8 kPa in 127 volunteers by one shear-wave '
                                        'system against 16.6 +/- 2.5 kPa in 59 volunteers by another, a '
                                        'factor of 5.7. Elastography vendors do not agree on whether the '
                                        'printed kPa is a shear or a Young modulus. No point value is '
                                        'defensible until that is resolved'),
            'young_modulus_lung_parenchyma': row(86500., 'Pa', 'transferred', 'nelson2026',
                                                 'human donor lungs, non-smokers, final (large-strain) '
                                                 'tensile stiffness modulus 86.5 +/- 60.0 kPa; smokers '
                                                 '238.6 +/- 128.5 kPa. This is a large-strain modulus and '
                                                 'is not comparable to a small-strain one'),
            'density_lung': row(384., 'kg/m3', 'transferred', 'kanematsu2015',
                                'ICRP 110 reference lung, 0.384 g/cm3, 1.4% of body mass. Cross-check: '
                                'healthy whole-lung mass at functional residual capacity is 997 +/- 133 g '
                                'by CT in 16 subjects, and ICRP 89 gives lung with blood 1200 g'),
            'density_soft_tissue': row(1060., 'kg/m3', 'transferred', 'icru44_nist',
                                       'ICRU-44 soft tissue; brain grey and white matter 1040, whole blood '
                                       '1060. No ICRU-44 soft tissue exceeds 1060'),
            'poisson_ratio': row(None, '1', A, None, 'no measurement retrieved for any soft organ'),
        },
        'hyperelastic': {
            'brain': {'law': 'one-term modified Ogden, W = 2 mu / alpha^2 (l1^a + l2^a + l3^a - 3)',
                      'tier': 'transferred', 'source': 'budday2017',
                      'regions': {'corpus_callosum': {'mu_pa': 350., 'alpha': 25.3},
                                  'corona_radiata': {'mu_pa': 660., 'alpha': 24.3},
                                  'basal_ganglia': {'mu_pa': 700., 'alpha': 18.7},
                                  'cortex': {'mu_pa': 1430., 'alpha': 19.0}},
                      'note': 'calibrated SIMULTANEOUSLY against shear, compression and tension on the '
                              'same specimens, which is why this is the strongest entry in this table. '
                              'Shear-only fits give 330 +/- 180, 540 +/- 210, 560 +/- 200 and '
                              '1060 +/- 360 Pa for the same four regions. Note the very large alpha: this '
                              'is a strongly nonlinear fit and a neo-Hookean reduction of it is wrong'},
            'myocardium': {'law': 'Holzapfel-Ogden orthotropic', 'tier': 'transferred', 'source': 'gao2015',
                           'parameters': {'a_pa': 134.8, 'b': 3.243, 'af_pa': 3176.2, 'bf': 4.7435,
                                          'as_pa': 542.6, 'bs': 1.5998, 'afs_pa': 234.4, 'bfs': 3.39},
                           'note': 'IN VIVO inverse estimation from cine MRI and cavity pressure in one '
                                   'volunteer of three (ages 22, 28, 31); this is NOT the Sommer 2015 ex '
                                   'vivo fit, whose coefficients were not retrieved. Requires a '
                                   'reference-space fibre and sheet field, which the atlas does not carry'},
            'penis': {'law': 'see data/measurements/biomechanics/khorshidi_2024.json',
                      'tier': 'transferred', 'source': 'khorshidi2024',
                      'note': 'Ogden for corpus cavernosum and spongiosum, neo-Hookean split for fascia, '
                              'HGO for tunica albuginea. Already in the repo and already consumed by '
                              'ihm/calibration/penile.py; the only fitted hyperelastic material the model '
                              'has'},
        },
    },
    'nerve': {
        'tissue': 'peripheral nerve',
        'linear': {
            'young_modulus': row(None, 'Pa', A, 'tang2021',
                                 'no defensible human adult peripheral nerve modulus was found. The '
                                 'elastography literature reports shear wave VELOCITY and Tang explicitly '
                                 'argues velocity is the more meaningful quantity for these anisotropic '
                                 'nerves, declining to print a Young modulus'),
            'shear_wave_velocity': row(3.37, 'm/s', 'transferred', 'tang2021',
                                       'median nerve at the forearm, median 3.37 (IQR 0.76) m/s in 105 '
                                       'healthy volunteers; median nerve at the carpal tunnel 4.30 (0.80); '
                                       'tibial nerve at the tarsal tunnel 3.87 (0.76)'),
            'shear_modulus_if_converted': row(11810., 'Pa', 'derived', 'tang2021',
                                              'G = rho v^2 at rho = 1040 kg/m3 and v = 3.37 m/s. THE '
                                              'SOURCE AUTHORS DECLINE THIS CONVERSION for anisotropic '
                                              'nerve. Recorded so the number is not silently reinvented '
                                              'downstream, not endorsed as a material property'),
            'density': row(None, 'kg/m3', A, None, 'no human nerve density measurement was found'),
        },
        'hyperelastic': row(None, None, A, None, 'none retrieved'),
    },
    'vascular': {
        'tissue': 'artery and vein wall',
        'linear': {
            'young_modulus_aorta_media_circumferential': row(101200., 'Pa', 'transferred', 'teng2015',
                                                             'human abdominal aortic media, incremental '
                                                             'Young modulus at stretch 1.0, 101.2 '
                                                             '[87.7, 142.9] kPa'),
            'young_modulus_aorta_media_axial': row(107300., 'Pa', 'transferred', 'teng2015',
                                                   '107.3 [82.3, 114.8] kPa; the media is close to '
                                                   'isotropic at zero stretch and strongly anisotropic '
                                                   'above it'),
            'poisson_ratio': row(None, '1', A, None, 'no measurement retrieved'),
            'density': row(1060., 'kg/m3', 'transferred', 'icru44_nist',
                           'ICRU-44 soft tissue / whole blood; no vessel-wall-specific value was found'),
        },
        'hyperelastic': row(None, None, A, 'holzapfel2005',
                            'no verified HGO coefficient set for a human vessel was retrieved. Teng fitted '
                            'a modified Mooney-Rivlin, not HGO, and its coefficients were not read'),
    },
    'connective_tissue': {
        'tissue': 'fascia, aponeurosis, membrane, mesentery',
        'linear': {k: row(None, u, A, None, 'no source found for this class; the canonical 100 kPa is an '
                                            'engineering placeholder and this candidate table declines to '
                                            'replace it with another one')
                   for k, u in (('young_modulus', 'Pa'), ('poisson_ratio', '1'), ('density', 'kg/m3'))},
        'hyperelastic': row(None, None, A, None, 'none retrieved'),
    },
    'lymph_node_group': {
        'tissue': 'lymph node groups',
        'linear': {k: row(None, u, A, None, 'no source found; the canonical 3 kPa is transferred from an '
                                            'ex-vivo kidney and liver microindentation study')
                   for k, u in (('young_modulus', 'Pa'), ('poisson_ratio', '1'), ('density', 'kg/m3'))},
        'hyperelastic': row(None, None, A, None, 'none retrieved'),
    },
}

# Per-entity overrides, where a source names the structure. Keyed by canonical
# entity name substring so the table survives id churn; the emitter resolves them
# against the canonical entity list and reports any that match nothing.
ENTITY_OVERRIDES = {
    'gracilis': {'shear_modulus': row(6150., 'Pa', 'transferred', 'chakouch2015',
                                      'gracilis 6.15 +/- 0.45 kPa, the stiffest thigh muscle at rest')},
    'semitendinosus': {'shear_modulus': row(5320., 'Pa', 'transferred', 'chakouch2015',
                                            'semitendinosus 5.32 +/- 0.10 kPa')},
    'sartorius': {'shear_modulus': row(5150., 'Pa', 'transferred', 'chakouch2015',
                                       'sartorius 5.15 +/- 0.19 kPa')},
    'rectus femoris': {'shear_modulus': row(3910., 'Pa', 'transferred', 'chakouch2015',
                                            'rectus femoris 3.91 +/- 0.16 kPa')},
    'vastus intermedius': {'shear_modulus': row(4230., 'Pa', 'transferred', 'chakouch2015',
                                                'vastus intermedius 4.23 +/- 0.25 kPa')},
    'lobe of left lung': {'density': row(384., 'kg/m3', 'transferred', 'kanematsu2015', 'ICRP 110 lung')},
    'lobe of right lung': {'density': row(384., 'kg/m3', 'transferred', 'kanematsu2015', 'ICRP 110 lung')},
    'liver': {'shear_modulus': row(2000., 'Pa', 'transferred', 'rouviere2006', 'healthy human liver MRE')},
    'kidney': {'shear_modulus': row(5000., 'Pa', 'transferred', 'arda2011', 'renal cortex SWE')},
    'pancreas': {'shear_modulus': row(4800., 'Pa', 'transferred', 'arda2011', 'pancreas SWE')},
    'thyroid': {'shear_modulus': row(10970., 'Pa', 'transferred', 'arda2011', 'thyroid SWE')},
}

# Entities the canonical table gives a value that no retrieved source supports and
# that this table declines to replace. Recorded so the defect is not lost.
UNSOURCED_DEFECTS = {
    'cornea': 'canonical 3 kPa, transferred from an ex-vivo kidney and liver microindentation study. '
              'Cornea is a dense collagenous lamellar tissue and no candidate value was sourced here',
    'sclera': 'canonical 3 kPa, same transfer; sclera is dense collagen',
    'lens': 'canonical 3 kPa, same transfer; the crystalline lens has a graded modulus and no candidate '
            'value was sourced here',
    'vitreous body': 'canonical 3 kPa, same transfer; the vitreous is a dilute gel and no candidate value '
                     'was sourced here',
    'gingiva': 'canonical 3 kPa, same transfer',
}

# Why a source may sit in the registry without being cited by a value row: it
# corroborates a row that cites another paper, it supplies a cross-check used only
# in the ledger or the conflict list, or its numbers could not be retrieved and the
# record exists to say so. Every such source must appear here.
CORROBORATING = {
    'maganaris1999': 'corroborates the tendon row (tibialis anterior 1.2 GPa); it is also the only tendon '
                     'source the canonical table cites',
    'hansen2006': 'corroborates the tendon row (patellar 1.09 GPa)',
    'ward2005': 'corroborates the muscle density row and supplies its fixation caveat',
    'chakouch2016': 'corroborates the muscle MRE row at a second frequency set and cohort',
    'athanasiou1995': 'corroborates the cartilage rows on a second human joint (ankle)',
    'athanasiou1998': 'corroborates the cartilage rows on a third human joint (first MTP)',
    'lee2013': 'corroborates the liver MRE row in 49 living donors',
    'pawlus2016': 'the second side of the unresolved spleen dispute; no value is carried',
    'morgan2003': 'establishes that no universal trabecular modulus-density law exists, which is why the '
                  'apparent trabecular modulus cell is absent',
    'guenard1992': 'measured healthy whole-lung mass, used only as the ledger cross-check on lung density',
    'organs_pmid33176223': 'the single source the CANONICAL table leans on; carried so the audit can name '
                           'it, not because this table transfers a value from it',
    'sommer2015': 'the standard human ex-vivo myocardium reference; its coefficients were not retrieved',
    'sommer2013': 'the standard human adipose multiaxial reference; its coefficients were not retrieved',
    'pailler2008': 'a widely quoted in-vivo skin indentation study whose number could not be confirmed',
    'holzapfel2005': 'the standard layer-specific human artery reference; not retrieved',
    'icrp89': 'used by the ledger and the bone rows through SYSTEM_DENSITY and ICRP_ADULT_MALE',
    'icru44_nist': 'used by the ledger through SYSTEM_DENSITY',
    'kanematsu2015': 'used by the ledger through SYSTEM_DENSITY',
    'khorshidi2024': 'referenced by the soft_organ hyperelastic block, which is not a value row',
    'gao2015': 'referenced by the soft_organ hyperelastic block, which is not a value row',
    'budday2017': 'referenced by the soft_organ hyperelastic block, which is not a value row',
    'vanloocke2006': 'referenced by the muscle anisotropy block',
    'morrow2010': 'referenced by the muscle anisotropy block',
    'gennisson2010': 'referenced by the muscle anisotropy block as the human gap',
    'takaza2013': 'referenced by the muscle hyperelastic block as a stress point, not a fit',
    'niannaidh2012': 'referenced by the skin rows and by the skin literature dispute',
    'reilly1975': 'referenced by the bone rows',
    'quapp1998': 'referenced by the ligament rows',
    'bayraktar2004': 'referenced by the bone rows',
    'vertebral2021': 'referenced by the bone rows',
    'liu1997': 'referenced by the cartilage rows',
    'alkhouli2013': 'referenced by the adipose rows',
    'feng2022': 'referenced by the skin and adipose rows',
    'chakouch2015': 'referenced by the muscle and adipose rows',
    'maganaris2002': 'referenced by the tendon row',
    'rouviere2006': 'referenced by the liver row',
    'arda2011': 'referenced by the kidney, pancreas and thyroid rows',
    'nelson2026': 'referenced by the lung row',
    'tang2021': 'referenced by the nerve rows',
    'teng2015': 'referenced by the vascular rows',
}


# ---------------------------------------------------------------- ledger
# Densities used to weigh the geometry, one per canonical `system`. Every one is
# transferred; the source key is carried so the ledger is auditable line by line.
SYSTEM_DENSITY = {
    'muscular': (1050., 'icru44_nist', 'ICRU-44 skeletal muscle'),
    'skeletal': (1300., 'icrp89', 'ICRP 89 Table 2.20 whole fresh adult skeleton; the atlas bone geometry '
                                  'encloses the marrow cavity, so the whole-skeleton density is the one '
                                  'that matches it, not the 1900 of marrow-free cortical bone'),
    'digestive': (1050., 'icru44_nist', 'ICRU-44 soft tissue band; no organ-specific density retrieved'),
    'integumentary': (1100., 'icrp89', 'ICRP 89 para 514, skin about 1.1 g/cm3. CAVEAT: the owned volume '
                                       'in this row is the 10 mm voxel footprint of a zero-thickness '
                                       'double-sided skin slab, not a measured dermal volume. It lands '
                                       'near the ICRP reference skin mass by coincidence of grid spacing, '
                                       'and the three real skin-layer entities are absent from both voxel '
                                       'partitions entirely'),
    'respiratory': (384., 'kanematsu2015', 'ICRP 110 reference lung; the airway walls in this system are '
                                           'denser but the lobes dominate the volume'),
    'nervous': (1040., 'icru44_nist', 'ICRU-44 brain grey and white matter'),
    'arterial': (1060., 'icru44_nist', 'ICRU-44 whole blood; these entities are lumen-shaped'),
    'venous': (1060., 'icru44_nist', 'ICRU-44 whole blood'),
    'cardiac': (1050., 'icru44_nist', 'ICRU-44 soft tissue band'),
    'connective': (1050., 'icru44_nist', 'ICRU-44 soft tissue band'),
    'urinary': (1050., 'icru44_nist', 'ICRU-44 soft tissue band'),
    'lymphatic': (1030., 'icru44_nist', 'ICRU-44 soft tissue band, lower bound'),
    'reproductive': (1050., 'icru44_nist', 'ICRU-44 soft tissue band'),
    'sensory': (1050., 'icru44_nist', 'ICRU-44 soft tissue band'),
    'endocrine': (1050., 'icru44_nist', 'ICRU-44 soft tissue band'),
}
FILL_ADIPOSE = 950.       # ICRU-44 / ICRP 110 adipose tissue
FILL_INTERSTITIAL = 1010.  # interstitial fluid, taken as slightly above water; assumed, not sourced
ICRP_ADULT_MALE = {'total_body_kg': 73.0, 'height_m': 1.76, 'body_surface_m2': 1.90,
                   'skeletal_muscle_kg': 29.0, 'adipose_tissue_kg': 18.2,
                   'separable_adipose_excl_yellow_marrow_kg': 14.5, 'storage_fat_kg': 14.6,
                   'skin_kg': 3.3, 'total_skeleton_kg': 10.5, 'cortical_bone_kg': 4.4,
                   'trabecular_bone_kg': 1.1, 'liver_kg': 1.8, 'kidneys_kg': 0.31, 'brain_kg': 1.45,
                   'lung_with_blood_kg': 1.2, 'lung_tissue_only_kg': 0.5, 'heart_tissue_only_kg': 0.33,
                   'spleen_kg': 0.15, 'blood_kg': 5.6,
                   'source': 'icrp89', 'tables': 'Table 2.8 and Table 2.9'}


def ledger(mech, profile, unmod, part):
    env_m3 = unmod['envelope']['divergence_volume_m3']
    env_area = unmod['envelope']['surface_area_m2']
    pr = unmod['primary']
    occ, void = pr['occupied_volume_m3'], pr['void_volume_m3']
    claim = pr['overlap']['summed_entity_claim_volume_m3']
    diverg = unmod['summed_abs_divergence_volume_unique_meshes_m3']
    target = profile['mass_kg']

    rows, geom_v, geom_m = [], 0., 0.
    for name, s in sorted(part['systems'].items(), key=lambda kv: -kv[1]['owned_volume_m3']):
        rho, src, note = SYSTEM_DENSITY[name]
        v = s['owned_volume_m3']
        rows.append({'system': name, 'owned_volume_m3': v, 'density_kg_m3': rho, 'density_source': src,
                     'density_note': note, 'candidate_mass_kg': v * rho,
                     'canonical_mass_kg': s['allocated_mass_kg'],
                     'canonical_effective_density_kg_m3': s['owned_effective_density_kg_m3']})
        geom_v += v
        geom_m += v * rho
    void10 = env_m3 - geom_v

    # The three canonical skin layers occupy no voxel at all and their volume is
    # area times thickness against a slab area that is nearly twice the envelope.
    layers = {e['name']: e for e in mech['entities'] if e['role'] == 'skin_layer'}
    slab_area = layers['epidermis']['material_volume_m3'] / 1e-4  # 0.1 mm epidermis
    layer_v = sum(e['material_volume_m3'] for e in layers.values())
    layer_m = sum(e['mass_kg'] for e in layers.values())
    thickness = sum(e['material_volume_m3'] for e in layers.values()) / slab_area
    corrected_layer_v = env_area * thickness
    skin_only_v = env_area * 0.0016  # epidermis 0.1 mm + dermis 1.5 mm; ICRP skin excludes hypodermis

    fills = {}
    for label, rho in (('all_adipose', FILL_ADIPOSE), ('all_interstitial', FILL_INTERSTITIAL)):
        fills[label] = {'fill_density_kg_m3': rho, 'fill_mass_kg': void10 * rho,
                        'implied_total_mass_kg': geom_m + void10 * rho}
    required = (target - geom_m) / void10

    adipose_v = ICRP_ADULT_MALE['separable_adipose_excl_yellow_marrow_kg'] / FILL_ADIPOSE
    split = {'void_volume_m3': void10,
             'icrp_separable_adipose_kg': ICRP_ADULT_MALE['separable_adipose_excl_yellow_marrow_kg'],
             'that_adipose_volume_m3': adipose_v,
             'residual_void_after_adipose_m3': void10 - adipose_v,
             'residual_at_interstitial_kg': (void10 - adipose_v) * FILL_INTERSTITIAL,
             'implied_total_mass_kg': geom_m + ICRP_ADULT_MALE['separable_adipose_excl_yellow_marrow_kg']
                                      + (void10 - adipose_v) * FILL_INTERSTITIAL}

    # Two-compartment densitometry cross-check on the profile itself.
    fat = profile['body_fat_fraction']
    siri_density = 495. / (fat * 100. + 450.)
    siri_volume = target / (siri_density * 1000.)
    bmi_profile = target / profile['height_m'] ** 2

    return {
        'envelope': {'divergence_volume_m3': env_m3, 'voxel_interior_volume_m3': pr['interior_volume_m3'],
                     'surface_area_m2': env_area,
                     'agreement_note': 'divergence and 8 mm voxel interior agree to %.3f%%'
                                       % (100 * abs(env_m3 - pr['interior_volume_m3']) / env_m3)},
        'occupancy': {'exclusive_occupied_8mm_m3': occ, 'void_8mm_m3': void,
                      'void_fraction': pr['void_fraction_of_interior'],
                      'largest_connected_void_m3': pr['void_components']['largest_volume_m3'],
                      'overlapping_entity_claim_8mm_m3': claim,
                      'summed_abs_divergence_unique_meshes_m3': diverg,
                      'exclusive_owned_10mm_m3': geom_v,
                      'voxelisation_shortfall_note':
                          'the meshes integrate to %.4f L but only %.4f L of 8 mm cells are claimed, so '
                          'voxelisation loses %.4f L (%.1f%%) on thin structures; 1016 entities claim no '
                          '8 mm cell at all. Part of the "void" is discretisation, not anatomy'
                          % (diverg * 1000, claim * 1000, (diverg - claim) * 1000,
                             100 * (diverg - claim) / diverg)},
        'geometry_mass': {'rows': rows, 'total_owned_volume_m3': geom_v,
                          'candidate_mass_kg': geom_m,
                          'candidate_mean_density_kg_m3': geom_m / geom_v,
                          'canonical_mass_kg': sum(r['canonical_mass_kg'] for r in rows)},
        'skin_slab_defect': {
            'slab_area_m2': slab_area, 'envelope_area_m2': env_area,
            'area_ratio': slab_area / env_area,
            'dubois_area_m2_for_this_profile': 0.007184 * (profile['height_m'] * 100) ** 0.725
                                               * target ** 0.425,
            'icrp_adult_male_area_m2': ICRP_ADULT_MALE['body_surface_m2'],
            'layer_thickness_m': thickness,
            'canonical_layer_volume_m3': layer_v, 'canonical_layer_mass_kg': layer_m,
            'corrected_layer_volume_at_envelope_area_m3': corrected_layer_v,
            'phantom_volume_m3': layer_v - corrected_layer_v,
            'skin_epidermis_plus_dermis_volume_m3': skin_only_v,
            'skin_epidermis_plus_dermis_mass_kg': skin_only_v * 1100.,
            'icrp_skin_kg': ICRP_ADULT_MALE['skin_kg'],
            'note': 'the canonical skin slab is %.4f m2 against a %.5f m2 watertight envelope, a factor '
                    'of %.3f. The 1.6 mm epidermis-plus-dermis prior is itself well supported: at the '
                    'envelope area and ICRP skin density it weighs %.3f kg against the ICRP reference '
                    '%.1f kg. It is the AREA that is wrong, not the thickness. The three layers also '
                    'occupy no voxel in either partition, so their %.3f kg is mass with no place to be'
                    % (slab_area, env_area, slab_area / env_area, skin_only_v * 1100.,
                       ICRP_ADULT_MALE['skin_kg'], layer_m)},
        'fill': {'void_volume_m3': void10, 'uniform_fill_cases': fills,
                 'composed_fill_case': split,
                 'density_required_to_reach_profile_kg_m3': required,
                 'note': 'reaching the %.4f kg profile from this geometry needs the void to weigh '
                         '%.0f kg/m3. Adipose tissue is 950 and no ICRU-44 soft tissue exceeds 1060. '
                         'Only bone is denser than the required fill'
                         % (target, required)},
        'reconciliation': {
            'profile_mass_kg': target, 'profile_height_m': profile['height_m'],
            'profile_body_fat_fraction': fat, 'profile_bmi': bmi_profile,
            'profile_basis': 'inherited BioGears StandardMale prior; only height comes from the atlas',
            'siri_body_density_kg_m3': siri_density * 1000.,
            'siri_implied_body_volume_m3': siri_volume,
            'siri_note': 'Siri two-compartment densitometry, %%fat = 495/D - 450. At %.0f%% fat the '
                         'profile mass implies %.3f L of body, before adding residual lung and gut gas. '
                         'The measured envelope is %.3f L. The deficit is %.3f L' %
                         (fat * 100, siri_volume * 1000, env_m3 * 1000, (siri_volume - env_m3) * 1000),
            'envelope_deficit_m3': siri_volume - env_m3,
            'mean_density_required_kg_m3': target / env_m3,
            'siri_fat_at_that_density_percent': 495. / (target / env_m3 / 1000.) - 450.,
            'geometry_supported_mass_range_kg': [fills['all_adipose']['implied_total_mass_kg'],
                                                 fills['all_interstitial']['implied_total_mass_kg']],
            'geometry_supported_bmi_range': [fills['all_adipose']['implied_total_mass_kg']
                                             / profile['height_m'] ** 2,
                                             fills['all_interstitial']['implied_total_mass_kg']
                                             / profile['height_m'] ** 2],
            'icrp_reference': ICRP_ADULT_MALE,
            'verdict': 'NO. The %.4f kg profile and this geometry cannot both stand without inventing '
                       'volume. Two independent routes agree. Densitometry: at the profile 21%% fat the '
                       'profile mass needs %.2f L and the envelope measures %.2f L, short by %.2f L, and '
                       'forcing %.4f kg into %.2f L needs a mean density of %.0f kg/m3, which the Siri '
                       'relation maps to %.1f%% body fat. Composition: weighing every owned voxel at a '
                       'sourced tissue density and filling the void with anything between adipose and '
                       'interstitial fluid gives %.1f to %.1f kg. The geometry is a BMI %.1f to %.1f '
                       'body; the profile is BMI %.1f. The profile is a BioGears inheritance, the '
                       'geometry is BodyParts3D, and they were never the same person'
                       % (target, siri_volume * 1000, env_m3 * 1000, (siri_volume - env_m3) * 1000,
                          target, env_m3 * 1000, target / env_m3, 495. / (target / env_m3 / 1000.) - 450.,
                          fills['all_adipose']['implied_total_mass_kg'],
                          fills['all_interstitial']['implied_total_mass_kg'],
                          fills['all_adipose']['implied_total_mass_kg'] / profile['height_m'] ** 2,
                          fills['all_interstitial']['implied_total_mass_kg'] / profile['height_m'] ** 2,
                          bmi_profile)},
    }


# ---------------------------------------------------------------- conflicts

def conflicts(mech):
    """Canonical value against candidate, flagged where they differ by over 2x."""
    ents = {e['name']: e for e in mech['entities']}
    by_role = {}
    for e in mech['entities']:
        by_role.setdefault(e['role'], e)
    scale = mech['mass_allocation']['uniform_scale']
    lung_v = sum(e['material_volume_m3'] for e in mech['entities'] if 'lung' in e['name'].lower())
    lung_m = sum(e['mass_kg'] for e in mech['entities'] if 'lung' in e['name'].lower())
    bone_v = sum(e['material_volume_m3'] for e in mech['entities'] if e['role'] == 'rigid_bone')
    bone_m = sum(e['mass_kg'] for e in mech['entities'] if e['role'] == 'rigid_bone')

    def ratio(a, b):
        return None if not a or not b else max(a / b, b / a)

    out = []

    def add(field, canonical, candidate, unit, source, comment):
        r = ratio(canonical, candidate)
        out.append({'field': field, 'canonical': canonical, 'candidate': candidate, 'unit': unit,
                    'ratio': r, 'over_2x': bool(r and r > 2.0), 'candidate_source': source,
                    'comment': comment})

    add('skin_layer.dermis.young_modulus', by_role['skin_layer']['material']['young_modulus']['value'],
        40000., 'Pa', 'feng2022',
        'the canonical value is one number for epidermis, dermis and hypodermis alike, transferred from '
        'an ex-vivo kidney and liver microindentation study. The in-vivo layer measurement spans '
        '15 kPa hypodermis to 4 MPa epidermis, a factor of 267 across the three layers the canonical '
        'table gives one value to')
    add('skin_layer.epidermis.young_modulus', by_role['skin_layer']['material']['young_modulus']['value'],
        4.0e6, 'Pa', 'feng2022', 'same single canonical value against the in-vivo epidermis')
    add('soft_organ.lung.young_modulus', ents['inferior lobe of left lung']['material']['young_modulus']['value'],
        86500., 'Pa', 'nelson2026',
        'the candidate is a LARGE-STRAIN tensile modulus of excised parenchyma and the canonical is a '
        'small-strain prior, so the factor is partly a regime difference; but 2 kPa carries no source at '
        'all and lung is not the softest tissue in the body')
    add('soft_organ.lung.density_effective', 1000. * scale, 384., 'kg/m3', 'kanematsu2015',
        'the canonical gives lung the same 1000 kg/m3 as every other soft tissue and then scales it to '
        '1195.8. ICRP 110 reference lung is 384')
    add('soft_organ.lung.mass', lung_m, ICRP_ADULT_MALE['lung_with_blood_kg'], 'kg', 'icrp89',
        'the five canonical lobes hold %.4f m3 and are assigned %.3f kg. The ICRP reference lung with '
        'blood is 1.2 kg and healthy whole-lung mass by CT at functional residual capacity is '
        '0.997 +/- 0.133 kg' % (lung_v, lung_m))
    add('ligament.young_modulus', by_role['ligament']['material']['young_modulus']['value'], 3.322e8, 'Pa',
        'quapp1998',
        'the canonical 100 MPa is cited to a TENDON paper, which is the wrong tissue. Human MCL is '
        '332 MPa along the fibre and 11 MPa across it, so no isotropic number is right and the canonical '
        'one is 3.3x below the fibre value and 9x above the transverse value')
    add('cartilage.poisson_ratio', by_role['cartilage']['material']['poisson_ratio']['value'], 0.08, '1',
        'liu1997',
        'the canonical 0.45 is the blanket soft-tissue assumption. Two independent human joints give 0.08 '
        'by biphasic creep indentation. The consequence is quantitative: at E = 1 MPa the aggregate '
        'modulus is 1.01 MPa at nu = 0.08 and 3.79 MPa at nu = 0.45, so the canonical cartilage is 3.8x '
        'too stiff in confined compression')
    add('rigid_bone.young_modulus', None, 1.70e10, 'Pa', 'reilly1975',
        'the canonical rigid_bone entities carry NO elastic constant of any kind: no Young modulus, no '
        'Poisson ratio, no shear modulus, no Lame parameter. 257 entities and 5.585 L. For a soft-body '
        'materialization in which bone is a very stiff inclusion rather than a rigid body, this is the '
        'single largest hole in the table')
    add('rigid_bone.density_effective', bone_m / bone_v, 1300., 'kg/m3', 'icrp89',
        'the canonical nominal is 1900 with prior range 1500-2200, and the uniform scale pushes the '
        'effective density to %.0f, outside its own declared prior. 2272 exceeds fresh marrow-free bone '
        '(1900-2000, ICRP 89 para 434) and approaches dry mineralised matrix (2300). The atlas geometry '
        'is whole bones including the marrow cavity, for which ICRP 89 gives 1300' % (bone_m / bone_v))
    add('soft_organ.liver.young_modulus', ents['caudate lobe of liver']['material']['young_modulus']['value'],
        6000., 'Pa', 'rouviere2006',
        'candidate is E = 3G from the 2.0 kPa healthy-liver MRE shear modulus. Note also that the atlas '
        'has no whole-liver entity, only a caudate lobe and its duct')
    add('all_soft_tissue.density_effective', 1000. * scale, 1060., 'kg/m3', 'icru44_nist',
        'not a 2x conflict, but the canonical effective soft-tissue density of %.1f sits outside its own '
        'declared prior range of 900-1100 and above every soft tissue in ICRU-44, whose maximum is 1060 '
        'for whole blood' % (1000. * scale))
    add('muscle.young_modulus', by_role['muscle']['material']['young_modulus']['value'], 11730., 'Pa',
        'chakouch2015', 'NOT a conflict: the canonical 20 kPa assumed prior is within 1.7x of E = 3G from '
        'in-vivo MRE. It still carries no source')
    add('tendon.young_modulus', by_role['tendon']['material']['young_modulus']['value'], 1.16e9, 'Pa',
        'maganaris2002', 'NOT a conflict: three human in-vivo tendons agree at 1.09-1.20 GPa and the '
        'canonical 1.2 GPa sits in that band. This is the best-supported cell in the canonical table')
    add('vascular.young_modulus', by_role['vascular']['material']['young_modulus']['value'], 101200., 'Pa',
        'teng2015', 'NOT a conflict: the canonical 100 kPa assumed prior lands on the human abdominal '
        'aortic media incremental modulus almost exactly. It is still an assumed value that happens to '
        'be right, not a sourced one')

    literature_disputes = [
        {'quantity': 'spleen stiffness by shear-wave elastography',
         'values': ['2.9 +/- 1.8 kPa, n=127 (arda2011)', '16.6 +/- 2.5 kPa, n=59 (pawlus2016)'],
         'ratio': 5.7,
         'comment': 'both are human, in vivo, healthy, and published. Shear-wave systems differ in '
                    'whether the printed kPa is a shear or a Young modulus, and the two papers were not '
                    'reconciled here. No spleen value is carried into the candidate table'},
        {'quantity': 'passive skeletal muscle anisotropy, sign of it',
         'values': ['tension: along-fibre 20x stiffer, rabbit (morrow2010)',
                    'compression: across-fibre 2.2x stiffer, porcine (vanloocke2006)'],
         'ratio': None,
         'comment': 'the two regimes give opposite orderings, and neither is human. Muscle is bimodular. '
                    'A single transversely isotropic parameter set cannot serve both, and there is no '
                    'human number to choose between them'},
        {'quantity': 'human skin modulus',
         'values': ['in-vivo dermis about 40 kPa at 0.2-1 kHz (feng2022)',
                    'excised whole skin 83.3 +/- 34.9 MPa, large strain, donors aged 89 (niannaidh2012)',
                    'the Ni Annaidh literature table itself spans 0.26-150 MPa across authors'],
         'ratio': None,
         'comment': 'these are different quantities (small-strain in-vivo layer stiffness against '
                    'large-strain excised collagen-engaged stiffness) and the spread inside the excised '
                    'literature alone is a factor of 577. Skin cannot be given one modulus'},
    ]
    return {'canonical_vs_candidate': out,
            'over_2x_count': sum(1 for r in out if r['over_2x']),
            'literature_vs_literature': literature_disputes}


# ---------------------------------------------------------------- emit

def tier_census(mech):
    """Verified-versus-assumed split of the candidate table, by cell and by volume."""
    cells = {}

    def walk(node):
        if isinstance(node, dict):
            if 'tier' in node and 'value' in node:
                cells[node['tier']] = cells.get(node['tier'], 0) + 1
                return
            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)
    walk(MATERIALS)
    role_volume = {}
    for e in mech['entities']:
        role_volume[e['role']] = role_volume.get(e['role'], 0.) + (e.get('material_volume_m3') or 0.)
    covered, uncovered = {}, {}
    for role, v in role_volume.items():
        target = covered if role in MATERIALS else uncovered
        target[role] = v
    total = sum(role_volume.values())
    return {'cells_by_tier': cells,
            'roles_with_a_candidate_row': sorted(covered),
            'roles_without_a_candidate_row': sorted(uncovered),
            'volume_with_a_candidate_row_m3': sum(covered.values()),
            'volume_without_a_candidate_row_m3': sum(uncovered.values()),
            'volume_fraction_with_a_candidate_row': sum(covered.values()) / total,
            'note': 'no cell in this table is tier `measured`. Nothing here was measured on the '
                    'BodyParts3D reference specimen, so the verified-versus-assumed split is really '
                    'transferred-versus-absent, and every transferred row inherits its donor cohort'}


def resolve_overrides(mech):
    hits, misses = {}, []
    names = [(e['id'], e['name']) for e in mech['entities']]
    for key, block in ENTITY_OVERRIDES.items():
        matched = [i for i, n in names if key in n.lower()]
        if matched:
            hits[key] = {'entities': matched, 'count': len(matched), 'material': block}
        else:
            misses.append(key)
    defects = {}
    for key, note in UNSOURCED_DEFECTS.items():
        matched = [i for i, n in names if key in n.lower()]
        defects[key] = {'entities': matched, 'count': len(matched), 'note': note}
    return {'resolved': hits, 'unmatched_keys': misses, 'unsourced_defects': defects}


def build():
    mech = json.loads(CANONICAL.read_text())
    profile = json.loads(PROFILE.read_text())
    unmod = json.loads(UNMODELLED.read_text())
    part = json.loads(PARTITION.read_text())
    return {
        'schema': 'ihm.tissue-material-candidate.v1',
        'status': 'CANDIDATE. Not canonical, not promoted, not consumed by any runtime.',
        'inputs_sha256': {str(p.relative_to(ROOT)): sha(p) for p in INPUTS},
        'self_sha256': sha(Path(__file__)),
        'tier_vocabulary': {
            'measured': 'measured on the BodyParts3D reference specimen. No row here reaches this tier.',
            'transferred': 'a published value measured on other tissue, other subjects or other species.',
            'derived': 'arithmetic on a transferred value; the operation is stated in the note.',
            'assumed': 'an engineering number with no source.',
            'absent': 'no defensible value was found. The cell is deliberately empty and must not be '
                      'silently filled.'},
        'audit': audit(mech),
        'sources': SOURCES,
        'materials': MATERIALS,
        'entity_overrides': resolve_overrides(mech),
        'conflicts': conflicts(mech),
        'mass_ledger': ledger(mech, profile, unmod, part),
        'tier_census': tier_census(mech),
        'not_established': [
            'nothing here was measured on this specimen; every row is a cohort transfer',
            'no prior_range or +/- carried here is a probability interval',
            'elastography moduli are dynamic at the driver frequency and are not quasi-static '
            'neo-Hookean shear moduli; the two are conflated everywhere in the tissue literature',
            'shear-wave elastography vendors disagree on whether printed kPa is shear or Young modulus, '
            'which is the likely cause of the 5.7x spleen dispute recorded in conflicts',
            'the mass ledger uses exclusive voxel ownership at 10 mm, which loses thin structures; the '
            'meshes integrate to more volume than the voxels claim and that difference is recorded',
            'no fibre or sheet field exists in the atlas, so the transversely isotropic muscle, the '
            'anisotropic ligament and the orthotropic myocardium rows cannot yet be instantiated',
            'articular cartilage has no geometry in either atlas, so its row has nothing to attach to',
        ],
        'canonical_assets_modified': False,
    }


def emit(out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    doc = build()
    parts = {'audit.json': doc['audit'], 'sources.json': doc['sources'],
             'materials.json': {'materials': doc['materials'],
                                'entity_overrides': doc['entity_overrides'],
                                'tier_vocabulary': doc['tier_vocabulary']},
             'conflicts.json': doc['conflicts'], 'mass_ledger.json': doc['mass_ledger']}
    for name, payload in parts.items():
        (out_dir / name).write_text(json.dumps(payload, indent=2, sort_keys=False) + '\n')
    manifest = {k: v for k, v in doc.items()
                if k not in ('audit', 'sources', 'materials', 'entity_overrides', 'conflicts',
                             'mass_ledger')}
    manifest['artifacts_sha256'] = {n: sha(out_dir / n) for n in parts}
    (out_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return out_dir, doc, manifest


def self_test():
    doc = build()
    a, c, l, t = doc['audit'], doc['conflicts'], doc['mass_ledger'], doc['tier_census']

    # audit reproduces the shipped totals exactly
    assert abs(a['totals']['mass_kg'] - 77.1107029) < 1e-6, a['totals']
    assert a['roles']['rigid_bone']['fields']['young_modulus']['distinct_values'] == []
    assert a['roles']['rigid_bone']['fields']['young_modulus']['entities_without_field'] == 257
    assert set(a['roles']['cartilage']['fields']['poisson_ratio']['distinct_values']) == {0.45}
    assert a['young_modulus_basis_share']['absent']['volume_fraction'] > 0.09
    assert len(a['distinct_literature_sources_cited']) == 4

    # every source record is complete and every material row points at a real source
    for key, s in SOURCES.items():
        for f in ('citation', 'doi', 'pmid', 'species', 'state', 'method', 'n', 'ages', 'verification'):
            assert f in s, (key, f)
        assert s['doi'] or s['pmid'] or 'arXiv' in s['citation'] or 'NIST' in s['citation'], key
    seen = set()

    def check(node):
        if isinstance(node, dict):
            if 'tier' in node and 'value' in node:
                assert node['tier'] in TIERS + ('derived', 'absent'), node
                if node['tier'] == 'absent':
                    assert node['value'] is None, node
                else:
                    assert node['value'] is not None, node
                    assert node['source'] in SOURCES, node
                if node.get('source'):
                    seen.add(node['source'])
                    assert node['source'] in SOURCES, node
                assert node.get('note'), node
                return
            if node.get('source') in SOURCES:
                seen.add(node['source'])
            for v in node.values():
                check(v)
        elif isinstance(node, list):
            for v in node:
                check(v)
    check(MATERIALS)
    assert 'measured' not in t['cells_by_tier'], t
    assert t['cells_by_tier'].get('absent', 0) > 0, 'an honest table has empty cells'
    unused = sorted(set(SOURCES) - seen)
    assert set(unused) <= set(CORROBORATING), sorted(set(unused) - set(CORROBORATING))
    assert set(CORROBORATING) <= set(SOURCES), sorted(set(CORROBORATING) - set(SOURCES))

    # ledger arithmetic closes
    g = l['geometry_mass']
    assert abs(sum(r['candidate_mass_kg'] for r in g['rows']) - g['candidate_mass_kg']) < 1e-9
    assert abs(sum(r['owned_volume_m3'] for r in g['rows']) - g['total_owned_volume_m3']) < 1e-12
    assert abs(g['canonical_mass_kg'] - 49.468) < 0.01, g['canonical_mass_kg']
    f = l['fill']
    assert f['density_required_to_reach_profile_kg_m3'] > 1100., f
    lo, hi = l['reconciliation']['geometry_supported_mass_range_kg']
    assert 68. < lo < hi < 73., (lo, hi)
    assert l['reconciliation']['siri_fat_at_that_density_percent'] < 0., l['reconciliation']
    assert l['skin_slab_defect']['area_ratio'] > 1.9, l['skin_slab_defect']
    assert abs(l['skin_slab_defect']['layer_thickness_m'] - 0.0066) < 1e-6

    # conflicts: the four indefensible cells are all present and flagged
    fields = {r['field']: r for r in c['canonical_vs_candidate']}
    assert fields['rigid_bone.young_modulus']['canonical'] is None
    assert fields['cartilage.poisson_ratio']['over_2x']
    assert fields['soft_organ.lung.density_effective']['over_2x']
    assert fields['ligament.young_modulus']['over_2x']
    assert not fields['tendon.young_modulus']['over_2x']
    assert not fields['vascular.young_modulus']['over_2x']
    assert c['over_2x_count'] >= 6, c['over_2x_count']
    assert len(c['literature_vs_literature']) == 3

    # writing goes to a temporary root, never the repo
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        d, _, man = emit(Path(tmp) / 'tissue-material-candidate-v1')
        for name, digest in man['artifacts_sha256'].items():
            assert sha(d / name) == digest, name
        assert json.loads((d / 'manifest.json').read_text())['canonical_assets_modified'] is False
    assert sha(CANONICAL) == doc['inputs_sha256']['data/derived/canonical/mechanics.json']

    report = {'status': 'passed', 'assertions': 'audit totals, source completeness, tier vocabulary, '
                                                'ledger closure, conflict flags, manifest digests',
              'sources': len(SOURCES), 'material_roles': len(MATERIALS),
              'cells_by_tier': t['cells_by_tier'],
              'over_2x_conflicts': c['over_2x_count'],
              'geometry_supported_mass_kg': l['reconciliation']['geometry_supported_mass_range_kg'],
              'profile_mass_kg': l['reconciliation']['profile_mass_kg'],
              'canonical_assets_modified': False}
    print(json.dumps(report, indent=2))
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--output', type=Path)
    p.add_argument('--summary', action='store_true', help='print the ledger and conflict summary')
    a = p.parse_args()
    if not (a.self_test or a.output or a.summary):
        p.error('Select --self-test, --summary or --output')
    if a.self_test:
        self_test()
    if a.output:
        d, doc, man = emit(a.output)
        print(json.dumps({'output_dir': str(d), 'artifacts': man['artifacts_sha256'],
                          'over_2x_conflicts': doc['conflicts']['over_2x_count'],
                          'cells_by_tier': doc['tier_census']['cells_by_tier']}, indent=2))
    if a.summary:
        doc = build()
        print(json.dumps({'mass_ledger': doc['mass_ledger'], 'conflicts': doc['conflicts'],
                          'tier_census': doc['tier_census']}, indent=2))
