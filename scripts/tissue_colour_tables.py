"""Colour values, their sources, and the mucosal mapping. Data only; no logic.

Read by scripts/build_tissue_colour_palettes.py, which converts every CIE L*a*b* entry to sRGB
and hashes this file into the build manifest, so a changed number is a changed build.

Every entry declares a tier from ihm.structure-provenance.v1 applied to colour:
  measured      CIE L*a*b* colorimetry of this tissue, in humans, in the state claimed.
  transferred   colorimetry of a different tissue, species, site or state, mapped across.
  derived       computed from a measured value in this table by a stated operation.
  synthesized   built from a descriptive, non-colorimetric source. No measurement exists here.

An entry tiered 'measured' MUST cite a source whose identifier re-resolved. The build self-test
enforces that, and enforces that a source with resolution_verified False backs nothing measured.

Palette intent: FRESH, IN VIVO, PERFUSED TISSUE under neutral white light, as at operation.
Not a formalin-fixed cadaver (desaturated, browner, blood drained), not a plastinate, not an
atlas plate, not a histological stain.

The honest summary: three tissue classes in this table have real human in-vivo colorimetry of
that tissue -- lip vermilion, gingiva and tooth enamel. Everything visceral has none; the search
for it is recorded in UNMEASURED_SEARCHED. The evidence palette exists so that split is visible
on the model rather than buried in this file.

A note on synthesized values. They are given directly as sRGB hex, not as an invented L*a*b*,
because converting a made-up L*a*b* would dress an assertion up as colorimetry. The builder
records each one's L*a*b* under D65/2 degree so they sit in the same space for comparison.
"""

# ---------------------------------------------------------------------------- sources
#
# resolution_verified is True only where the identifier was re-resolved during this build and the
# returned title compared against the title recorded here. resolution_note records how.

_NCBI = ('NCBI E-utilities esummary, db=pubmed, batched request on 2026-09-07; the returned title, '
         'journal, year and DOI were compared against the record here')
_CROSSREF = ('api.crossref.org/works/<doi> fetched 2026-09-07; the returned title and container '
             'were compared against the record here')

SOURCES = {

    # ======================================================== oral and perioral, human, in vivo
    'vergnaud2024lip': {
        'citation': 'Vergnaud H, et al. Lip color diversity: an intricate study. Skin Res Technol. '
                    '2024;30(2):e13583.',
        'pmid': '38284291', 'doi': '10.1111/srt.13583',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo, lower (inferior) lip vermilion, 8.2 mm2 regions',
        'cohort': 'n=410 healthy women, 19-68 y (mean 42.5, SD 12.7), France and USA: 207 Caucasian '
                  '(French and American), 103 African American, 100 Hispanic American',
        'instrument': 'SpectraFace (Newtone Technologies) hyperspectral imaging, 30 bands 410-700 nm',
        'illuminant_observer_stated': True, 'illuminant_observer': 'D65_10',
        'measurement': 'Overall L*=41.9+/-6.6, C*=23.6+/-6.0, h=34.0+/-4.7 deg (range L* 24.6-54.4). '
                       'The paper reports C* and h, not a* and b*. By group: Caucasian L*46.0/C*26.7/h31.9; '
                       'African American L*33.2/C*16.2/h36.7; Hispanic American L*42.3/C*25.0/h35.6. '
                       'Darker lips are less saturated and more yellow; lighter lips more saturated and redder.'},

    'vergnaud2023lipdevice': {
        'citation': 'Vergnaud H, et al. Lip color measurement: a new hyperspectral imaging device. '
                    'Skin Res Technol. 2023;29(8):e13418.',
        'pmid': '37632193', 'doi': '10.1111/srt.13418',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo lip vermilion',
        'cohort': 'n=19 healthy French Caucasian women, 21-67 y (43.3+/-14.7), Fitzpatrick I-III',
        'instrument': 'three devices compared: SpectraFace hyperspectral; X-Rite VS3200 spectrophotometer '
                      '45:0 geometry, 400-700 nm at 10 nm; Canfield VISIA-CR cross-polarized camera',
        'illuminant_observer_stated': True, 'illuminant_observer': 'D65_10',
        'measurement': 'SpectraFace L*47.36+/-3.37, a*20.02+/-2.86, b*12.54+/-1.77. VS3200 L*48.29+/-3.17, '
                       'a*18.61+/-2.28, b*11.55+/-1.82. VISIA-CR L*42.31+/-5.78, a*29.03+/-3.85, b*11.93+/-2.73.',
        'caveat': 'On the same 19 lips, contact spectrophotometry and cross-polarized imaging differ by '
                  'about 10 units of a*. Instrument family matters more than cohort; do not mix spaces.'},

    'wang2025lipage': {
        'citation': 'Wang Y, et al. Biophysical characteristics of lip vermilion among healthy individuals '
                    'in southern China. Skin Res Technol. 2025;31(10):e70228.',
        'pmid': '41074489', 'doi': '10.1111/srt.70228',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo lower lip vermilion',
        'cohort': 'n=217 healthy volunteers, 18-55 y (37.9+/-8.2), 28 M / 189 F, southern China',
        'instrument': 'Canfield VISIA-CR (n=180) and Delfin SkinColorCatch contact colorimeter (n=151)',
        'illuminant_observer_stated': False,
        'measurement': 'VISIA-CR by age: 18-35 L*45.27/a*23.81/b*14.25; 36-45 43.51/22.55/13.60; '
                       '46-55 42.25/21.52/13.31. SkinColorCatch by age: 52/22/13, 51/21/12, 51/20/12. '
                       'L*, a* and b* all decline significantly with age (p<0.05). No upper-lip colour given.'},

    'ho2015gingiva': {
        'citation': 'Ho DK, et al. Color range and color distribution of healthy human gingiva: a '
                    'prospective clinical study. Sci Rep. 2015;5:18498.',
        'pmid': '26691598', 'doi': '10.1038/srep18498',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo keratinized gingiva, 2-3 mm apical to the mid-facial margin '
                                     'of tooth 8 or 9, 6 mm spot',
        'cohort': 'n=238 adults, 97 M / 141 F, USA: 42 African-American, 54 Asian, 82 Caucasian, 60 Hispanic',
        'instrument': 'Photo Research PR-670 spectroradiometer with MS-75 lens',
        'illuminant_observer_stated': False,
        'illuminant_observer_assumed': 'D65_2',
        'measurement': 'Overall L*=52.9+/-5.2, a*=23.3+/-3.4, b*=14.9+/-2.0 (ranges L* 37.2-64.0, '
                       'a* 13.4-31.7, b* 9.2-22.2). By group: African-American 50.6/20.4/14.3; Asian '
                       '50.8/24.8/15.7; Caucasian 54.7/23.3/14.4; Hispanic 53.8/24.1/15.1. Largest '
                       'inter-group dE* = 5.0. Ethnicity and age significant, sex not.',
        'caveat': 'The paper does not state an illuminant or standard observer. D65/2 degree is assumed '
                  'for the conversion here and that assumption is carried in the palette record.'},

    'gomezpolo2024gingiva': {
        'citation': 'Gomez-Polo C, et al. Explaining the colour of natural healthy gingiva. Odontology. '
                    '2024;112(4):1284-1295.',
        'pmid': '38403674', 'doi': '10.1007/s10266-024-00906-4',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo, free gingival margin / attached middle zone / mucogingival line',
        'cohort': 'n=360 (187 M / 173 F), 18-92 y (47.2+/-18.8), Caucasian, Spain',
        'instrument': 'SpectroShade Micro, non-contact, no tissue drying',
        'illuminant_observer_stated': False,
        'measurement': 'Free margin M 49.8/24.1/14.5, F 50.9/23.6/15.3; attached middle M 50.3/24.4/14.7, '
                       'F 50.9/24.0/15.8; mucogingival line M 49.8/23.7/14.5, F 49.5/23.5/15.2. The three '
                       'gingival zones differ by under 1.5 dE, so gingival colour is uniform along the '
                       'crown-apical axis. Only b* differs by sex, imperceptibly.'},

    'naranjo2023gingiva': {
        'citation': 'Naranjo MJ, et al. Study of attached gingiva space color according to gender and age '
                    'in Caucasian population. J Esthet Restor Dent. 2023;35(6):834-841.',
        'pmid': '36951233', 'doi': '10.1111/jerd.13038',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo attached gingiva, 2.5 mm apical to the upper central incisor zenith',
        'cohort': 'n=216 Caucasian (129 F / 87 M), three age bands, Spain',
        'instrument': 'SpectroShade Micro',
        'illuminant_observer_stated': False,
        'measurement': 'ENVELOPE ONLY (means are in paywalled tables): L* 40.4-61.2, a* 17.0-30.2, '
                       'b* 9.8-21.9. Age significantly affects b* (p=0.000); attached gingiva becomes '
                       'bluer with age.'},

    'huang2011gingiva': {
        'citation': 'Huang JW, et al. Using a spectrophotometric study of human gingival colour '
                    'distribution to develop a shade guide. J Dent. 2011;39 Suppl 3:e11-e16.',
        'pmid': '22005337', 'doi': '10.1016/j.jdent.2011.10.001',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo healthy gingiva',
        'cohort': 'n=362, grouped by sex and age',
        'instrument': 'reflectance spectrometer',
        'illuminant_observer_stated': False,
        'measurement': 'DIRECTION ONLY: ten gingival colour categories derived; significant sex difference '
                       'with dE > 3.7, female gingiva significantly lighter. Category coordinates paywalled.'},

    'hyun2017gingiva': {
        'citation': 'Hyun HK, et al. Colorimetric distribution of human attached gingiva and alveolar '
                    'mucosa. J Prosthet Dent. 2017;117(2):294-302.',
        'pmid': '27666499', 'doi': '10.1016/j.prosdent.2016.06.009',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo attached gingiva and alveolar (labial/buccal) mucosa, 23 sites',
        'cohort': 'n=40 periodontally healthy adults, 22 M / 18 F, 25-36 y',
        'instrument': 'colorimeter',
        'illuminant_observer_stated': False,
        'measurement': 'DIRECTION ONLY (tables paywalled): attached gingiva has higher L* and lower a* than '
                       'alveolar mucosa; attached gingiva is lighter and less red in the maxilla than the '
                       'mandible, and yellower in incisor than molar regions. This is the closest evidence '
                       'located for wet labial/buccal mucosa and it carries no absolute coordinates.'},

    'cho2025tongue': {
        'citation': 'Cho E, et al. Temporal changes in tongue color during immune checkpoint inhibitor '
                    'therapy in patients with non-small-cell lung cancer. Oncol Rev. 2025;19:1697252.',
        'pmid': '41445868', 'doi': '10.3389/or.2025.1697252',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo tongue, by region, calibrated imaging',
        'cohort': 'n=140, median 68.5 y (42-87), 115 M / 25 F, South Korea, 10 hospitals. NSCLC patients '
                  'at baseline, NOT healthy controls.',
        'instrument': 'KIOM computerized tongue image analysis system with in-frame ColorChecker calibration',
        'illuminant_observer_stated': False,
        'measurement': 'Baseline means (SE): body L*53.6 (0.278), a*25.3 (0.291); coating L*51.1, a*13.0; '
                       'centre L*56.6, a*23.0; tip a*28.6; side a*26.9; root a*18.6. b* not reported by '
                       'region. The tongue is reddest at the tip and desaturates posteriorly; coating drops '
                       'a* by about 12 units.',
        'caveat': 'Oncology cohort at baseline, not a healthy reference. No b* is published, so any b* used '
                  'with these numbers comes from elsewhere and the entry that does so is a transfer.'},

    'tian2024tongue': {
        'citation': 'Tian Z, et al. Association between color value of tongue and T2DM based on '
                    'dose-response analyses using restricted cubic splines in China. Medicine (Baltimore). '
                    '2024;103(25):e38575.',
        'pmid': '38905430', 'doi': '10.1097/MD.0000000000038575',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo tongue, mid-dorsum',
        'cohort': 'n=2439 total, healthy group n=1448 (717 M / 731 F), Tianjin, China',
        'instrument': 'TFDA-1 tongue diagnostic instrument, cool-white LED 5000 K at 2354 lux; RGB converted '
                      'to CIELAB by matrix, NOT a spectrophotometer',
        'illuminant_observer_stated': False,
        'measurement': 'Healthy L*=69.52+/-6.53, a*=14.19+/-4.28, b*=4.00 (median, IQR 2-5).',
        'caveat': 'NOT USED for a colour here. L* near 70 with b* near 4 is far lighter and far less yellow '
                  'than any contact-spectrophotometer mucosal reading and is very likely a white-balance '
                  'artefact of the RGB pipeline. Recorded so the rejection is on the record.'},

    'hosoki2007smoking': {
        'citation': 'Hosoki M. [Analysis of color changes of oral mucosa by smoking]. Kokubyo Gakkai Zasshi. '
                    '2007;74(2):108-118.',
        'pmid': '17682458', 'doi': '10.5357/koubyou.74.108',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo oral mucosa at multiple sites including buccal mucosa',
        'cohort': 'n=62 nonsmokers and 56 smokers, 30-83 y, Japan',
        'instrument': 'colorimeter, CIE L*a*b*',
        'illuminant_observer_stated': False,
        'measurement': 'DIRECTION ONLY (J-Stage returned HTTP 500; tables not retrieved): significant '
                       'differences in lip and gingival L* and a*, and buccal mucosa L* and b*. Smoking '
                       'lowers luminosity and shifts toward blue.'},

    'yamashiro1996munsell': {
        'citation': 'Yamashiro M. [A study on colorimetry of oral mucosal lesions]. Kokubyo Gakkai Zasshi. '
                    '1996;63(1):188-207.',
        'pmid': '8725366', 'doi': '10.5357/koubyou.63.188',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo oral mucosa',
        'cohort': 'not retrieved', 'instrument': 'Munsell colorimetry, not CIELAB',
        'illuminant_observer_stated': False,
        'measurement': 'Normal oral mucosa hue 5.0R-4.1YR, value 3.5-6.0, chroma 3.7-6.7. Munsell value '
                       '3.5-6.0 corresponds roughly to L* 35-60, consistent with the gingival data.',
        'caveat': 'Munsell, not CIELAB. Used only as a consistency check, never converted.'},

    'thibodeau1997lip': {
        'citation': 'Thibodeau EA, DAmbrosio JA. Measurement of lip and skin pigmentation using reflectance '
                    'spectrophotometry. Eur J Oral Sci. 1997;105(4):373-375.',
        'pmid': '9298371', 'doi': '10.1111/j.1600-0722.1997.tb00255.x',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo lip and adjacent facial skin',
        'cohort': 'white, olive and black skin types', 'instrument': 'reflectance spectrophotometry',
        'illuminant_observer_stated': False,
        'measurement': 'Melanin and haemoglobin indices, NOT L*a*b*. Establishes that both melanin and '
                       'haemoglobin are significantly higher in upper and lower lip than in adjacent facial '
                       'skin across skin types, i.e. that lip and skin must not share a colour.'},

    'takamoto2013oral': {
        'citation': 'Takamoto A, Sugahara K, Shibahara T, Katakura A, Matsuzaka K, Sugihara N. Screening for '
                    'oral mucosal diseases by a portable spectrophotometer. J Oral Maxillofac Surg Med '
                    'Pathol. 2013;25:314-327.',
        'pmid': None, 'doi': '10.1016/j.ajoms.2012.10.007',
        'resolution_verified': True, 'resolution_note': _CROSSREF+'; the DOI resolves to this exact article',
        'species': 'human', 'state': 'in vivo oral mucosa, multiple sites',
        'cohort': 'not retrieved', 'instrument': 'portable spectrophotometer',
        'illuminant_observer_stated': False,
        'measurement': 'NO NUMBERS TAKEN. Figures for buccal and lingual mucosa circulate in search snippets '
                       'and plausibly come from this paper, but the full text was not retrievable (publisher '
                       '403, no PMC copy, abstract elided) and the numbers were never seen in the source. '
                       'The citation is verified; the numbers are not, and none of them is used here.'},

    # ======================================================== skin pigmentation classification
    'delbino2018ita': {
        'citation': 'Del Bino S, Duval C, Bernerd F. Clinical and biological characterization of skin '
                    'pigmentation diversity and its consequences on UV impact. Int J Mol Sci. 2018;19(9):2668.',
        'pmid': '30205563', 'doi': '10.3390/ijms19092668',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'constitutive (sun-protected) skin colour',
        'cohort': 'review and cohort work establishing the six-group ITA classification',
        'instrument': 'reflectance colorimetry',
        'illuminant_observer_stated': False,
        'measurement': 'The six-group Individual Typology Angle classification used here. The formula and '
                       'the class boundaries were read verbatim from an independent open-access paper that '
                       'cites this one (see itoformula2026), not from this paper directly.'},

    'itoformula2026': {
        'citation': 'Ito S, et al. Unraveling UVA1-induced photomodifications of eumelanin and pheomelanin '
                    'in human skin. Int J Mol Sci. 2026;27(9):3973.',
        'pmid': '42123554', 'doi': '10.3390/ijms27093973',
        'resolution_verified': True,
        'resolution_note': 'Open-access full text fetched from PMC13163989 on 2026-09-07; the ITA formula '
                           'and the six class boundaries were read verbatim out of it, and it attributes the '
                           'classification to Del Bino 2018',
        'species': 'human', 'state': 'constitutive skin colour classification',
        'cohort': 'n/a for the classification itself',
        'instrument': 'n/a',
        'illuminant_observer_stated': False,
        'measurement': 'ITA = arctan((L* - 50)/b*) x 180/pi. Six groups: very light > 55 deg; light > 41; '
                       'intermediate > 28; tan > 10; brown > -30; dark <= -30.'},

    'chardon1991ita': {
        'citation': 'Chardon A, Cretois I, Hourseau C. Skin colour typology and suntanning pathways. Int J '
                    'Cosmet Sci. 1991;13(4):191-208.',
        'pmid': None, 'doi': '10.1111/j.1467-2494.1991.tb00561.x',
        'resolution_verified': True, 'resolution_note': _CROSSREF,
        'species': 'human', 'state': 'constitutive and facultative skin colour',
        'cohort': 'not retrieved', 'instrument': 'reflectance colorimetry',
        'illuminant_observer_stated': False,
        'measurement': 'The paper that introduced the Individual Typology Angle. Cited as the origin of the '
                       'construct; no numbers were taken from it.'},

    'delbino2013variation': {
        'citation': 'Del Bino S, Bernerd F. Variations in skin colour and the biological consequences of '
                    'ultraviolet radiation exposure. Br J Dermatol. 2013;169(s3):33-40.',
        'pmid': None, 'doi': '10.1111/bjd.12529',
        'resolution_verified': True, 'resolution_note': _CROSSREF,
        'species': 'human', 'state': 'constitutive skin colour across groups',
        'cohort': 'not retrieved (publisher returned 403)', 'instrument': 'reflectance colorimetry',
        'illuminant_observer_stated': False,
        'measurement': 'NO NUMBERS TAKEN. Cited as the standard reference for ITA-classified skin colour '
                       'variation; the full text was not retrievable in this build.'},

    'xiao2017skin': {
        'citation': 'Xiao K, Yates JM, Zardawi F, Sueeprasan S, Liao N, Gill L, Li C, Wuerger S. '
                    'Characterising the variations in ethnic skin colours: a new calibrated data base for '
                    'human skin. Skin Res Technol. 2017;23(1):21-29.',
        'pmid': '27273806', 'doi': '10.1111/srt.12295',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo skin at four sites: forehead, cheek, back of hand, inner forearm',
        'cohort': 'n=960: 187 Caucasian (UK), 202 Chinese, 145 Kurdish (Iraq), 426 Thai; 18-75 y, mostly 20-40',
        'instrument': 'Konica Minolta CM-2600d (3 mm aperture) and X-Rite SP62 (4 mm), d/8 geometry, '
                      'specular component included',
        'illuminant_observer_stated': True, 'illuminant_observer': 'D65_2',
        'measurement': 'Inner forearm means: Caucasian L*63.0 a*5.6 b*14.0; Chinese 60.9/7.0/15.0; Kurdish '
                       '60.6/6.5/16.4; Thai 61.9/7.1/17.4. Redness is nearly constant across groups but '
                       'strongly site-dependent (forehead highest, inner arm lowest); yellowness is strongly '
                       'group-dependent and weakly site-dependent; lightness is highest at the inner arm and '
                       'lowest at the forehead in all four groups.'},

    'lu2026issa': {
        'citation': 'Lu Y, Xiao K, Li C, Pointer M. Skin colour does not define ethnicity: quantifying '
                    'variation and overlap across diverse populations. Skin Res Technol. 2026;32(3):e70343.',
        'pmid': '41822993', 'doi': '10.1111/srt.70343',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo skin pooled across ten body sites including palm',
        'cohort': 'n=2113 subjects, 15256 spectra, International Skin Spectra Archive, 2012-2024, eight '
                  'population groups across the UK, Spain, China, Japan, Pakistan, Thailand, Iraq, Saudi Arabia',
        'instrument': 'Konica Minolta CM-2600d and CM-700d, X-Rite SP62; CIE di:8 degree, specular included',
        'illuminant_observer_stated': True, 'illuminant_observer': 'D65_2',
        'measurement': 'Group means (pooled across sites): Caucasian L*61.2 a*11.0 b*14.7; Chinese '
                       '59.8/10.1/16.7; Japanese 63.5/9.8/17.0; South Asian 52.7/10.4/17.9; African '
                       '39.6/10.2/14.4; Middle Eastern (Iraqi) 57.3/10.2/15.8; Southeast Asian (Thai) '
                       '56.4/10.0/19.0; Arabian 60.2/10.9/17.9. 89.4 % of individuals have a perceptually '
                       'indistinguishable counterpart in another group; median group gamut overlap 60.5 %.',
        'caveat': 'Pooled across body sites, so these are NOT constitutive (sun-protected) values and are '
                  'not used as a skin option here. Recorded because the overlap finding is the reason this '
                  'model offers skin colour as a declared parameter rather than a category.'},

    'lu2025issadata': {
        'citation': 'Lu Y, Xiao K, Pointer M, et al. The International Skin Spectra Archive (ISSA): a '
                    'multicultural human skin phenotype and colour spectra collection. Sci Data. '
                    '2025;12(1):487.',
        'pmid': '40122935', 'doi': '10.1038/s41597-025-04857-5',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo skin, ten body sites',
        'cohort': 'the archive behind lu2026issa',
        'instrument': 'Konica Minolta CM-2600d and CM-700d, X-Rite SP62; CIE di:8 degree, specular included',
        'illuminant_observer_stated': True, 'illuminant_observer': 'D65_2',
        'measurement': 'The data descriptor. Recorded for the instrument and geometry declaration, and '
                       'because it establishes that the palm-to-dorsum lightness contrast scales with '
                       'constitutive pigmentation rather than being a constant offset. No per-site numbers '
                       'are taken from it here.'},

    'sommers2019inguinal': {
        'citation': 'Sommers MS, Regueira Y, Tiller DA, et al. Understanding rates of genital-anal injury: '
                    'role of skin color and skin biomechanics. J Forensic Leg Med. 2019;66:120-128.',
        'pmid': '31299484', 'doi': '10.1016/j.jflm.2019.06.019',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human',
        'state': 'in vivo constitutive skin at a sun-protected site: the right inner upper thigh, two '
                 'inches below the groin',
        'cohort': 'n=341 women aged 21 and over, Philadelphia (USA) and San Juan (Puerto Rico); the study '
                  'reports its groups as non-Hispanic White (n=88), non-Hispanic Black (n=54), '
                  'Hispanic/Latina (n=190) and other (n=9)',
        'instrument': 'ColorTec PSM hand-held reflectance spectrophotometer',
        'illuminant_observer_stated': False, 'illuminant_observer_assumed': 'D65_2',
        'measurement': 'Non-Hispanic White L*64.39+/-3.39 a*7.71+/-1.42 b*18.12+/-2.66; Hispanic/Latina '
                       '55.93+/-6.92 / 9.15+/-1.67 / 20.26+/-2.55; non-Hispanic Black 41.05+/-6.23 / '
                       '10.16+/-0.94 / 19.40+/-2.97.',
        'caveat': 'One instrument, one sun-protected site, three groups measured in the same protocol, which '
                  'is why this study rather than a cross-study patchwork supplies the measured skin options. '
                  'The group labels are the study own descriptors and are reproduced verbatim; they are '
                  'social categories as recorded by that study, not colour classes. The authors state they '
                  'could not measure genital mucous membranes because moisture caused instrument error.'},

    'sommers2013genital': {
        'citation': 'Sommers M, Beacham B, Baker R, Fargo J. Intra- and inter-rater reliability of digital '
                    'image analysis for skin color measurement. Skin Res Technol. 2013;19(4):484-491.',
        'pmid': '23551208', 'doi': '10.1111/srt.12072',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo external genital mucosa (labia) and vaginal wall',
        'cohort': '210 colposcopic images, light- and temperature-controlled laboratory',
        'instrument': 'Leica DFC420 C on a colposcope, uncompressed TIFF, colour-corrected against a Munsell '
                      'ColorChecker Mini in Photoshop CS4. Image-derived, not a contact colorimeter.',
        'illuminant_observer_stated': False, 'illuminant_observer_assumed': 'D65_2',
        'measurement': 'Labia L*49.67+/-10.35 a*27.74+/-9.45 b*17.52+/-13.67 (range L* 35.4-67.2); vaginal '
                       'wall 56.39/24.05/9.58; adjacent skin 45.07/12.86/22.12.',
        'caveat': 'The best-controlled genital mucosal colorimetry located: chart-calibrated in a controlled '
                  'lab. It is still image-derived and states no illuminant or standard observer, and it does '
                  'not distinguish labia majora from minora. Every value in this domain is device-relative.'},

    'baker2010genital': {
        'citation': 'Baker RB, Fargo JD, Shambley-Ebron D, Sommers MS. A source of healthcare disparity: '
                    'race, skin color, and injuries after rape among adolescents and young adults. J '
                    'Forensic Nurs. 2010;6(3):144-150.',
        'pmid': '21175535', 'doi': '10.1111/j.1939-3938.2010.01070.x',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo external genital mucous membrane, vaginal wall and adjacent skin',
        'cohort': 'n=234 patients aged 14-29 (mean 20.9), USA; groups reported by the study as non-Hispanic '
                  'Black (n=100) and non-Hispanic White (n=131)',
        'instrument': 'Adobe Photoshop CS2 colorimetry on colposcope-mounted camera images; no colorimeter, '
                      'no stated calibration target',
        'illuminant_observer_stated': False,
        'measurement': 'External genital mucous membrane: non-Hispanic Black L*46.99+/-11.24 a*15.80+/-8.05 '
                       'b*17.94+/-13.00; non-Hispanic White 58.47+/-11.51 / 26.52+/-9.53 / 25.40+/-14.88. '
                       'Vaginal wall 61.81/25.23/21.72 and 59.43/27.10/22.26. Adjacent skin 53.82/14.79/20.91 '
                       'and 75.04/17.57/27.18.',
        'caveat': 'Uncalibrated image analysis with very large SDs. Recorded for the consistent direction it '
                  'shares with the better-controlled 2013 study -- genital mucosa is redder, less yellow and '
                  'darker than the adjacent keratinized skin, with the gap narrowing as constitutive '
                  'pigmentation rises -- not for its absolute numbers.'},

    'huang2018foreskin': {
        'citation': 'Huang WS, Wang YW, Hung KC, et al. High correlation between skin color based on CIELAB '
                    'color space, epidermal melanocyte ratio, and melanocyte melanin content. PeerJ. '
                    '2018;6:e4815.',
        'pmid': '29844968', 'doi': '10.7717/peerj.4815',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'ex vivo foreskin, OUTER cutaneous surface, measured immediately after '
                                     'circumcision with blood and fat removed',
        'cohort': 'n=15 Asian young adults, 24.5+/-1.0 y (range 21-38), Taipei, Taiwan',
        'instrument': 'Konica Minolta Color Reader CR-10, triplicate readings',
        'illuminant_observer_stated': False,
        'measurement': 'Per-donor L* 39.43 to 52.37; cohort L* about 47.3 (SD 3.6 across donors). a* and b* '
                       'were acquired and deliberately not reported.',
        'caveat': 'Outer cutaneous foreskin, not the inner preputial mucosa, and L* only, so it cannot give '
                  'a colour. The only instrumented measurement of any preputial tissue located.'},

    'horibata2025sites': {
        'citation': 'Horibata K, Kondo S, Hashimoto S, Takemura Y. An observational study to determine the '
                    'optimal physical evaluation site for detecting anemia. J Gen Fam Med. 2025;26(3):246-254.',
        'pmid': '40291064', 'doi': '10.1002/jgf2.776',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human',
        'state': 'in vivo palm (thenar eminence), inner upper arm, thumbnail, lower lip and lower palpebral '
                 'conjunctiva',
        'cohort': 'n=92 Japanese outpatients, mean age 68.5+/-15.5, 51 % men; 67 non-anaemic, 25 anaemic',
        'instrument': 'Konica Minolta CM-700d. The palpebral conjunctival value alone is indirect: '
                      'photographed against a CASMATCH chart, colour-corrected, inkjet-printed and the print '
                      'measured, because direct contact was avoided for infection control and retinal light '
                      'safety.',
        'illuminant_observer_stated': False, 'illuminant_observer_assumed': 'D65_2',
        'measurement': 'Non-anaemic group: palm (thenar) L*62.9+/-2.9 a*7.2+/-2.2 b*16.7+/-2.8; inner upper '
                       'arm 65.4+/-3.9 / 5.0+/-1.2 / 15.8+/-2.5; thumbnail 55.4+/-9.3 / 4.6+/-1.5 / '
                       '11.1+/-2.4; lower lip 59.1+/-12.6 / 34.6+/-6.2 / 12.7+/-8.9; lower palpebral '
                       'conjunctiva 50.0+/-7.4 / 36.7+/-5.7 / 21.5+/-7.5.',
        'caveat': 'A single Japanese cohort with a wide age range, and the paper is about anaemia detection, '
                  'so its sites were chosen for haemoglobin sensitivity. The conjunctival value went through '
                  'an inkjet round trip, which is a real gamut limitation.'},

    'schmalwieser2024sites': {
        'citation': 'Schmalwieser AW, Gotzinger S, Schwabel F. Exploratory study on the body distribution of '
                    'skin color, pigmentation and degree of tan in Central European Caucasian women. '
                    'Photochem Photobiol Sci. 2024;23(3):493-502.',
        'pmid': '38351275', 'doi': '10.1007/s43630-024-00533-6',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo skin at 18 body sites, measured in February and in May/June',
        'cohort': 'n=20 Austrian Caucasian women, 20-60 y (mean 39.7), Fitzpatrick I-III; sun-seekers and '
                  'solarium users excluded',
        'instrument': 'Minolta Chroma Meter CR-300, xenon flash',
        'illuminant_observer_stated': True, 'illuminant_observer': 'D65_2',
        'measurement': 'SITE VARIATION, values in figures only. The inner UPPER arm near the axilla is the '
                       'valid constitutive baseline: it did not change winter to summer and did not '
                       'correlate with age. The inner arm near the elbow and the inner forearm are '
                       'explicitly unsuitable as baselines. Winter site range L*68 a*6.7 (inner upper arm, '
                       'lightest) to L*62 a*12.6 (anterior thigh, nape, shoulder blade). Winter-to-summer '
                       'mean shift across all sites dL* -2.0, db* +0.7, da* negligible. Repeatability '
                       '+/-0.31 L*, +/-0.31 a*, +/-0.23 b*, +/-1.1 degrees ITA.',
        'caveat': 'The reason this palette uses one skin value for all 220 topographic regions and says so: '
                  'real site variation is roughly 6 L* units and 6 a* units, and no per-site table is '
                  'published in this paper to apply.'},

    'chien2016aging': {
        'citation': 'Chien AL, Suh J, Cesar SSA, et al. Pigmentation in African American skin decreases with '
                    'skin aging. J Am Acad Dermatol. 2016;75(4):782-787.',
        'pmid': '27318769', 'doi': '10.1016/j.jaad.2016.05.007',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'sun-protected buttock and sun-exposed dorsal forearm',
        'cohort': 'n=40 Caucasian and 43 African American, Baltimore USA',
        'instrument': 'tristimulus colorimetry',
        'illuminant_observer_stated': False,
        'measurement': 'DIRECTION ONLY (paywalled): in African Americans aged 18-30 the sun-protected '
                       'buttock was DARKER than the sun-exposed forearm (P<0.001), whereas in Caucasians the '
                       'buttock was LIGHTER than the forearm (P<0.001).',
        'caveat': 'Recorded because it inverts the usual constitutive-versus-facultative assumption in darkly '
                  'pigmented skin, so "sun-protected equals lightest" is not a safe global rule.'},

    'russell2014sclera': {
        'citation': 'Russell R, Sweda JR, Porcheron A, Mauger E. Sclera color changes with age and is a cue '
                    'for perceiving age, health, and beauty. Psychol Aging. 2014;29(3):626-635.',
        'pmid': '25244481', 'doi': '10.1037/a0036142',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo sclera, whole visible region, from studio portraits',
        'cohort': 'n=286 French Caucasian women, 20-70 y, no make-up, no coloured lenses',
        'instrument': 'studio photography and MATLAB image analysis; 8-bit Lab encoding, no absolute white '
                      'reference, no stated illuminant or observer',
        'illuminant_observer_stated': False,
        'measurement': 'DIRECTION ONLY. Sclera darkens, reddens and yellows with age: L* decreases '
                       '(F=79.0), a* increases (F=150.3), b* increases (F=167.2), all p<0.0001. No '
                       'coefficients or per-decade means are published.',
        'caveat': 'NOT USED as a value. The published numbers are an 8-bit MATLAB Lab encoding of '
                  'un-white-balanced portraits and include lid shadow and vessels, which is why they convert '
                  'to an implausibly dark sclera. Retained for the age direction and for the finding that '
                  'sclera is not neutral: a* and b* are both positive.'},

    'papas2000conjunctiva': {
        'citation': 'Papas EB. Key factors in the subjective and objective assessment of conjunctival '
                    'erythema. Invest Ophthalmol Vis Sci. 2000;41(3):687-691.',
        'pmid': '10711682', 'doi': None,
        'resolution_verified': False,
        'resolution_note': 'NOT re-resolved in this build. Recorded from the research pass without an '
                           'independent identifier check, so it backs nothing measured.',
        'species': 'human', 'state': 'in vivo bulbar conjunctiva',
        'cohort': 'n=21', 'instrument': 'image analysis',
        'illuminant_observer_stated': False,
        'measurement': 'Perceived conjunctival redness is driven by vessel coverage (vessel area R2=0.93, '
                       'vessel count R2=0.90) far more than by any colour variable (best R2=0.62).',
        'caveat': 'Recorded because it says something a colour palette cannot express: conjunctival redness '
                  'is vasculature, not a tissue tint. Even if a conjunctival entity existed, a flat colour '
                  'would be the wrong representation.'},

    # ======================================================== human colorimetry, other tissues
    'wee2023enamel': {
        'citation': 'Wee AG, Winkelmann DA, Gozalo DJ, Ito M, Johnston WM. Color and translucency of enamel '
                    'in vital maxillary central incisors. J Prosthet Dent. 2023;130(6):878-884.',
        'pmid': '35184886', 'doi': '10.1016/j.prosdent.2022.01.010',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'in vivo, vital unrestored maxillary central incisor, mid-incisal',
        'cohort': 'n=120 subjects, 60 M / 60 F, 4 ethnic groups x 5 age bands 18-85 y',
        'instrument': 'Photo Research PR-705 spectroradiometer, 0/45 geometry',
        'illuminant_observer_stated': True, 'illuminant_observer': 'D65_2',
        'measurement': 'Enamel at infinite thickness L*=73.5+/-7.6, a*=2.2+/-1.8, b*=11.9+/-8.4; '
                       'translucency parameter 10.1+/-3.6.'},

    'schafer2001skull': {
        'citation': 'Schafer AT. The colour of the human skull. Forensic Sci Int. 2001;117(1-2):53-56.',
        'pmid': '11230946', 'doi': '10.1016/s0379-0738(00)00448-5',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'fresh autopsy skull bone, defleshed',
        'cohort': 'n=124 skull samples at autopsy; no outer-vs-inner table difference; b* rises with age',
        'instrument': 'tristimulus colorimeter (Micro Color)',
        'illuminant_observer_stated': False,
        'measurement': 'L*=72.5, a*=-7.4, b*=16.4 (means; SDs not in the retrieved record).',
        'caveat': 'NOT USED as the bone colour. a*=-7.4 is green-shifted, which is what defleshed autopsy '
                  'skull looks like, not living cortical bone under periosteum in a bleeding field. '
                  'Converting it gives a green-tan (#b0b594) that would misrepresent the palette intent. '
                  'The rejection is recorded on the cortical_bone entry.'},

    'rubio2020bone': {
        'citation': 'Rubio L, et al. Spectrophotometric color measurement to assess temperature of exposure '
                    'in cortical and medullar heated human bones: a preliminary study. Diagnostics (Basel). '
                    '2020;10(11):979.',
        'pmid': '33233746', 'doi': '10.3390/diagnostics10110979',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'cadaveric long-bone sections; unheated 21 C control arm',
        'cohort': 'n=36 sections (femur, tibia, radius, ulna) from one male cadaver aged 67',
        'instrument': 'Spectro-color (Dr Lange), 8 mm tip',
        'illuminant_observer_stated': True, 'illuminant_observer': 'D65_10',
        'measurement': 'The unheated control L*a*b* means are presented only graphically and were not '
                       'extractable. Retained for the instrument and illuminant declaration, not a value.'},

    'popciutrila2016dentin': {
        'citation': 'Pop-Ciutrila IS, Ghinea R, Colosi HA, Dudea D. Dentin translucency and color evaluation '
                    'in human incisors, canines, and molars. J Prosthet Dent. 2016;115(4):475-481.',
        'pmid': '26548886', 'doi': '10.1016/j.prosdent.2015.07.015',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': '2 mm midcoronal dentin slabs, extracted teeth',
        'cohort': 'n=33 incisors, 7 canines, 33 molars',
        'instrument': 'VITA Easyshade Compact and SpectraScan PR-704',
        'illuminant_observer_stated': False,
        'measurement': 'DIRECTION ONLY: anterior dentin has higher L*, lower a* and lower b* than molar '
                       'dentin. Per-group means are in paywalled tables and were not retrieved.'},

    'ishimoto2009cartilage': {
        'citation': 'Ishimoto Y, Hattori K, Ohgushi H, Uematsu K, Tanikake Y, Tanaka Y, Takakura Y. '
                    'Spectrocolorimetric evaluation of human articular cartilage. Osteoarthritis Cartilage. '
                    '2009;17(9):1204-1208.',
        'pmid': '19328879', 'doi': '10.1016/j.joca.2009.02.014',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'articular cartilage graded by Outerbridge class',
        'cohort': 'not retrieved (paywalled tables)',
        'instrument': 'spectrocolorimeter; full reflectance spectrum plus L*a*b*',
        'illuminant_observer_stated': False,
        'measurement': 'DIRECTION ONLY: L*, a* and yellow/red spectral reflectance differ significantly '
                       'across the four macroscopic grades. Intact grade-1 cartilage reflectance rises '
                       'monotonically with wavelength; degenerate grades develop a dip near 580 nm and lose '
                       'reflectance overall.'},

    'shimbashi2024kidney': {
        'citation': 'Shimbashi S, Yoshimiya M, Tashiro A, Noriki S, Hyodoh H. "Shock kidney-like appearance": '
                    'objective evaluation of renal color changes in hemorrhagic shock deaths. Leg Med '
                    '(Tokyo). 2024;71:102521.',
        'pmid': '39191046', 'doi': '10.1016/j.legalmed.2024.102521',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'post-mortem bisected kidney, digital photograph analysed in L*a*b*',
        'cohort': 'n=122 autopsy cases, 83 M / 39 F, mean age 64.8 y',
        'instrument': 'digital camera plus ImageJ, not a contact colorimeter',
        'illuminant_observer_stated': False,
        'measurement': 'RELATIVE ONLY: cortico-medullary dL* and da*, not absolute organ colour. Normal '
                       'control kidneys have a cortex redder than medulla by more than the da* = -1.33 '
                       'shock cutoff (AUC 0.859).',
        'caveat': 'Establishes that cortex and medulla differ in a* and the sign of the difference, not '
                  'either absolute colour. It cannot be applied here anyway: the kidney is a single closed '
                  'surface in this geometry with no cortex/medulla partition.'},

    'gomezgavara2024liver': {
        'citation': 'Gomez-Gavara C, et al. Enhanced artificial intelligence methods for liver steatosis '
                    'assessment using machine learning and color image processing. Clin Transplant. '
                    '2024;38(10):e15465.',
        'pmid': '39382065', 'doi': '10.1111/ctr.15465',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'deceased-donor liver grafts, in situ and backbench, colour-calibrated '
                                     'photography',
        'cohort': 'n=192 livers, 362 photographs, 7240 patches, brain-death donors',
        'instrument': 'colour-calibrated smartphone photography; L*a*b* plus local binary pattern features',
        'illuminant_observer_stated': False,
        'measurement': 'DIRECTION ONLY: steatotic livers acquire a yellowish tone, i.e. steatosis raises L* '
                       'and b*. No mean L*a*b* per steatosis grade is published.'},

    'piella2024livercolor': {
        'citation': 'Piella G, et al. LiverColor: an artificial intelligence platform for liver graft '
                    'assessment. Diagnostics (Basel). 2024;14(15):1654.',
        'pmid': '39125531', 'doi': '10.3390/diagnostics14151654',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'donor liver graft photographs',
        'cohort': 'not retrieved (publisher returned 403)',
        'instrument': 'colour image analysis platform',
        'illuminant_observer_stated': False,
        'measurement': 'No numeric colour values were confirmed to exist in the retrievable record.'},

    'kneifel2022liverhsi': {
        'citation': 'Kneifel F, Wagner T, Flammang I, et al. Hyperspectral imaging for viability assessment '
                    'of human liver allografts during normothermic machine perfusion. Transplant Direct. '
                    '2022;8(12):e1420.',
        'pmid': '36406899', 'doi': '10.1097/TXD.0000000000001420',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'deceased-donor liver allograft under normothermic machine perfusion',
        'cohort': 'n=25 allografts imaged at 1, 2 and 4 h',
        'instrument': 'hyperspectral camera; oxygenation and haemoglobin indices, not colour coordinates',
        'illuminant_observer_stated': False,
        'measurement': 'StO2 median 49.8 [29.4-72.4] rising to 56.6 %; tissue haemoglobin index 65.5 '
                       '[44.5-86.8]. A perfused human liver sits near half-saturated haemoglobin at high '
                       'haemoglobin concentration, which is the physical reason it reads reddish-brown '
                       'rather than arterial red.'},

    'belasco2020urine': {
        'citation': 'Belasco R, Edwards T, Munoz AJ, Rayo V, Buono MJ. The effect of hydration on urine '
                    'color objectively evaluated in CIE L*a*b* color space. Front Nutr. 2020;7:576974.',
        'pmid': '33195369', 'doi': '10.3389/fnut.2020.576974',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human', 'state': 'voided urine across hydration states',
        'cohort': 'n=151 samples from 28 healthy adults (22 M, 6 F), 28.6+/-11.3 y',
        'instrument': 'HunterLab Vista spectrophotometer',
        'illuminant_observer_stated': False,
        'measurement': 'L*=95.8+/-3.4, a*=-2.2+/-1.5, b*=18.5+/-12.2 overall; b* rises near-linearly with '
                       'dehydration, a* is non-monotonic.',
        'caveat': 'Recorded but not applied: no urine-filled lumen is a separate entity in this body.'},

    # ======================================================== non-human, used as explicit transfers
    'papanikolopoulou2025beef': {
        'citation': 'Papanikolopoulou V, et al. Impact of breed and slaughter hygiene on beef carcass '
                    'quality traits in northern Greece. Foods. 2025;14(10):1776.',
        'pmid': '40428555', 'doi': '10.3390/foods14101776',
        'resolution_verified': True, 'resolution_note': _NCBI+'; DOI additionally confirmed via '+_CROSSREF,
        'species': 'BOVINE (Bos taurus)',
        'state': 'fresh longissimus dorsi, 24 h post-mortem, bloomed, retail day 1',
        'cohort': 'n=159 carcasses across 4 breeds',
        'instrument': 'Konica Minolta CR-410, 50 mm aperture',
        'illuminant_observer_stated': True, 'illuminant_observer': 'C_2',
        'measurement': 'Aberdeen Angus L*38.2+/-2.09 a*19.5+/-2.17 b*10.2+/-1.53; Limousin 40.1/21.3/7.0; '
                       'Holstein 37.7/17.6/7.8; crossbred 38.9/22.2/5.4. Four-breed mean L*38.7 a*20.2 b*7.6.',
        'caveat': 'Bloomed post-mortem beef is oxymyoglobin-rich at the cut surface, drained of blood and '
                  'pH-shifted. Living perfused human muscle is darker and less saturated than this, so the '
                  'transferred value is, if anything, too bright.'},

    'schelkopf2021beefmethod': {
        'citation': 'Schelkopf CS, et al. Nix Pro Color Sensor provides comparable color measurements to '
                    'HunterLab colorimeter for fresh beef. J Food Sci Technol. 2021;58(9):3661-3665.',
        'pmid': '34366483', 'doi': '10.1007/s13197-021-05077-6',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'BOVINE', 'state': 'fresh longissimus thoracis',
        'cohort': 'n=200 carcasses', 'instrument': 'HunterLab MiniScan vs Nix Pro',
        'illuminant_observer_stated': False,
        'measurement': 'Instrument agreement r=0.80-0.85 for L*, a*, b*. Correlations, not means. Retained '
                       'as method support for the transferred muscle value, not as a colour.'},

    'dunne2009fat': {
        'citation': "Dunne PG, Monahan FJ, O'Mara FP, Moloney AP. Colour of bovine subcutaneous adipose "
                    'tissue: a review of contributory factors, associations with carcass and meat quality. '
                    'Meat Sci. 2009;81(1):28-45.',
        'pmid': '22063959', 'doi': '10.1016/j.meatsci.2008.06.013',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'BOVINE', 'state': 'subcutaneous adipose tissue',
        'cohort': 'review', 'instrument': 'Minolta colorimeter across the reviewed studies',
        'illuminant_observer_stated': False,
        'measurement': 'MECHANISM ONLY: adipose b* (yellowness) tracks total carotenoid content, r=0.79 '
                       '(P<0.01); forage-fed carcass fat is markedly yellower than grain-fed.',
        'caveat': 'Bovine, and no usable absolute L*a*b*. The carotenoid mechanism carries over to humans; '
                  'the magnitude does not.'},

    'basson2021bloodcolour': {
        'citation': 'Basson EP, Zeiler GE, Kamerman PR, Meyer LCR. Use of blood colour for assessment of '
                    'arterial oxygen saturation in immobilized impala. Vet Anaesth Analg. 2021;48(5):725-733.',
        'pmid': '34362689', 'doi': '10.1016/j.vaa.2021.05.004',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'IMPALA (Aepyceros melampus)',
        'state': 'arterial blood, anaerobic draw, measured immediately',
        'cohort': 'n=10 adult females, 163 arterial samples',
        'instrument': 'spectrocolourimetry plus co-oximetry',
        'illuminant_observer_stated': False,
        'measurement': 'L*, a* and b* each fit SaO2 best as a quadratic over SaO2 15-95 %. A generated '
                       'colour palette was ordered correctly by observers over 45-95 % saturation. The '
                       'coefficients were not retrieved.',
        'caveat': 'Impala, not human, and no coefficients in hand. The only L*a*b*-versus-saturation mapping '
                  'located anywhere; recorded as the route a later build should take.'},

    # ======================================================== spectra and optics
    'zijlstra1991hb': {
        'citation': 'Zijlstra WG, Buursma A, Meeuwsen-van der Roest WP. Absorption spectra of human fetal '
                    'and adult oxyhemoglobin, de-oxyhemoglobin, carboxyhemoglobin, and methemoglobin. Clin '
                    'Chem. 1991;37(9):1633-1638.',
        'pmid': '1716537', 'doi': None,
        'resolution_verified': True, 'resolution_note': _NCBI+'; the PubMed record carries no DOI',
        'species': 'human', 'state': 'haemoglobin derivatives in solution',
        'cohort': 'n/a', 'instrument': 'spectrophotometry, 450-1000 nm',
        'illuminant_observer_stated': False,
        'measurement': 'Millimolar absorptivities of oxyHb, deoxyHb, COHb and metHb for fetal and adult human '
                       'haemoglobin. The canonical dataset from which arterial (bright crimson, 97-100 % '
                       'SaO2) and venous (dark burgundy, 70-75 % SaO2) colour can be computed by integrating '
                       'against the CIE colour-matching functions and an illuminant. Venous blood is dark '
                       'red, never blue; the blue-vein percept is subsurface scattering in skin.',
        'caveat': 'This build did not perform that integration, so the blood colours here stay synthesized '
                  'and only their ordering rests on this source.'},

    'jacques2013optics': {
        'citation': 'Jacques SL. Optical properties of biological tissues: a review. Phys Med Biol. '
                    '2013;58(11):R37-R61.',
        'pmid': '23666068', 'doi': '10.1088/0031-9155/58/11/R37',
        'resolution_verified': True, 'resolution_note': _NCBI,
        'species': 'human and animal, pooled review', 'state': 'various',
        'cohort': 'review', 'instrument': 'n/a',
        'illuminant_observer_stated': False,
        'measurement': 'Absorption and reduced-scattering formulae parameterised by blood, water, melanin, '
                       'fat and yellow-pigment volume fractions plus the small/large scatterer balance. '
                       'Colour follows from these by Monte Carlo or Kubelka-Munk to a reflectance spectrum '
                       'and then to CIE tristimulus values; no colour coordinates are published.',
        'caveat': 'The principled route to a procedural palette, not taken in this build. Recorded so the '
                  'next build knows where to start.'},

    'bashkatov2005optics': {
        'citation': 'Bashkatov AN, Genina EA, Kochubey VI, Tuchin VV. Optical properties of human skin, '
                    'subcutaneous and mucous tissues in the wavelength range from 400 to 2000 nm. J Phys D '
                    'Appl Phys. 2005;38:2543-2555.',
        'pmid': None, 'doi': '10.1088/0022-3727/38/15/004',
        'resolution_verified': True, 'resolution_note': _CROSSREF,
        'species': 'human', 'state': 'ex vivo skin, subcutaneous and mucous tissue',
        'cohort': 'not retrieved', 'instrument': 'integrating-sphere spectrophotometry',
        'illuminant_observer_stated': False,
        'measurement': 'Absorption and scattering coefficients 400-2000 nm, including mucous tissue, so the '
                       'visible band is in principle convertible to mucosal colour.',
        'caveat': 'Not converted in this build.'},

    # ======================================================== descriptive and relational
    'gross_description': {
        'citation': 'Standard gross-anatomical and surgical-pathology description of fresh tissue as given '
                    'in general reference works: Standring S (ed), Gray\'s Anatomy, 42nd ed., Elsevier, '
                    '2021; Kumar V, Abbas AK, Aster JC (eds), Robbins & Cotran Pathologic Basis of Disease, '
                    '10th ed., Elsevier, 2020.',
        'pmid': None, 'doi': None,
        'resolution_verified': False,
        'resolution_note': 'NOT an identifier and NOT verified at page level. These are printed books that '
                           'were not retrievable in this build, so no chapter or page locator is asserted. '
                           'The entry exists so a synthesized colour names the kind of source it rests on '
                           'instead of appearing from nowhere. It may not back a measured tier, and the '
                           'self-test enforces that.',
        'species': 'human', 'state': 'fresh, in vivo or at operation',
        'cohort': 'n/a', 'instrument': 'none; verbal description, no colorimetry',
        'illuminant_observer_stated': False,
        'measurement': 'Descriptions only. Liver reddish-brown under a glistening capsule; spleen dark '
                       'purple-red, softer and darker than liver; kidney cortex pale reddish-brown over a '
                       'darker, more red-purple medulla; pancreas pale pinkish-tan and lobulated; lung pink '
                       'in a child, darkening with age to grey-pink with black anthracotic streaking along '
                       'subpleural interlobular septa; tendon and aponeurosis pearly glistening white; '
                       'peripheral nerve pale cream-white and matte with visible fascicular striation; dura '
                       'mater off-white to pearly grey-white; fresh unfixed cortex grey-pink to pinkish-tan '
                       '(fixation is what greys it) over creamy-white white matter; bile golden-yellow to '
                       'olive-green; adipose buttery cream-yellow; normal articular cartilage smooth, '
                       'bluish-white and translucent; arterial adventitia pale pink-tan, not the red of the '
                       'blood inside; thin-walled veins dark red-purple from the desaturated blood showing '
                       'through them.'},

    'anatomical_relation': {
        'citation': 'Assignment by anatomical relation rather than by a source: a structure takes the colour '
                    'of the tissue that forms its visible surface, or of a neighbouring structure of the '
                    'same tissue type. The relation is stated on each entry that uses this.',
        'pmid': None, 'doi': None,
        'resolution_verified': False,
        'resolution_note': 'Not a citation. Recorded so no colour in this table lacks a stated basis.',
        'species': 'human', 'state': 'n/a', 'cohort': 'n/a', 'instrument': 'none',
        'illuminant_observer_stated': False, 'measurement': 'n/a'},
}


# ---------------------------------------------------------------------------- colours

_DESC = ['gross_description']
_REL = ['anatomical_relation']

COLOURS = {

    # ============================================ mucosa and mucocutaneous transitional zones
    # The specific gap this build closes. These surfaces are pink to red-brown and quite distinct
    # from skin; none of them inherits a skin or a system colour.
    'mucosa_lip_vermilion': {
        'label': 'Lip vermilion',
        'state': 'healthy adult lower lip vermilion, in vivo',
        'lab': [41.9, 19.565, 13.197], 'illuminant_observer': 'D65_10',
        'tier': 'measured',
        'tier_basis': 'Vergnaud 2024, n=410 women across four ethnic cohorts, hyperspectral imaging, '
                      'D65/10 degree. The paper reports L*C*h; a* and b* here are C*cos(h) and C*sin(h) '
                      'of the overall mean L*41.9 C*23.6 h=34.0 deg, computed in this build.',
        'sources': ['vergnaud2024lip', 'vergnaud2023lipdevice', 'wang2025lipage', 'thibodeau1997lip'],
        'note': 'The manifest entity is FMA "lip", the whole lip including its cutaneous part; the '
                'vermilion is not separately meshed, so a vermilion colour is applied to a region that is '
                'partly skin. Vermilion is thin non-keratinized epithelium over a dense capillary bed, '
                'which is why it reads red rather than as skin, and Thibodeau 1997 shows both melanin and '
                'haemoglobin are higher in lip than in adjacent facial skin. This is the multi-ethnic mean '
                'and it is deliberately NOT co-varied with the skin-tone option: the source cohorts are '
                'ethnic self-identification, not ITA class, and no mapping between the two is established. '
                'The measured group spread is L* 33.2 (African American) to 46.0 (Caucasian).'},

    'mucosa_gingival': {
        'label': 'Gingiva',
        'state': 'healthy keratinized attached and marginal gingiva, in vivo',
        'lab': [52.9, 23.3, 14.9], 'illuminant_observer': 'D65_2',
        'tier': 'measured',
        'tier_basis': 'Ho 2015, n=238 adults across four ethnic groups, PR-670 spectroradiometer. The paper '
                      'states no illuminant or observer; D65/2 degree is assumed for the conversion and the '
                      'assumption is carried in the palette record.',
        'sources': ['ho2015gingiva', 'gomezpolo2024gingiva', 'naranjo2023gingiva', 'huang2011gingiva'],
        'note': 'Three independent cohorts agree closely: Ho L*52.9/a*23.3/b*14.9, Gomez-Polo about '
                '50/24/15, Naranjo envelope L* 40.4-61.2 / a* 17.0-30.2 / b* 9.8-21.9. Gingival colour is '
                'uniform along the crown-apical axis (under 1.5 dE between free margin, attached zone and '
                'mucogingival line) but varies with pigmentation (dE 5.0 between the extreme groups), so '
                'pigmentation belongs on a global gingival tint, not a per-zone one. Gingiva is lighter '
                'than lip vermilion by about 11 L* units at similar a*.'},

    'mucosa_palatal': {
        'label': 'Palatal mucosa',
        'state': 'healthy soft palate and uvula, in vivo',
        'lab': [52.9, 23.3, 14.9], 'illuminant_observer': 'D65_2',
        'tier': 'transferred',
        'tier_basis': 'Gingival colorimetry (Ho 2015) transferred verbatim to palatal mucosa. Both are '
                      'masticatory mucosa: keratinized stratified squamous epithelium bound down to '
                      'periosteum, with the same vascular architecture. No palatal colorimetry exists -- a '
                      'full-text search of the open-access corpus for palatal mucosa with CIELAB returned '
                      'zero hits, and every relevant PubMed record measures graft or prosthesis colour '
                      'match rather than absolute mucosal coordinates.',
        'sources': ['ho2015gingiva', 'yamashiro1996munsell'],
        'note': 'No offset was invented. The soft palate is in fact lining mucosa rather than masticatory '
                'mucosa and is probably redder than this; nothing measured sets the magnitude.'},

    'mucosa_lingual': {
        'label': 'Tongue, dorsal mucosa',
        'state': 'healthy tongue dorsum, in vivo',
        'lab': [53.6, 25.3, 14.9], 'illuminant_observer': 'D65_2',
        'tier': 'transferred',
        'tier_basis': 'Composite. L*=53.6 and a*=25.3 are the tongue-body baseline of Cho 2025, measured in '
                      'an oncology cohort at baseline rather than in healthy controls. That paper publishes '
                      'no b*, so b*=14.9 is transferred from the gingival measurement. The illuminant is '
                      'not stated in either source; D65/2 degree is assumed.',
        'sources': ['cho2025tongue', 'ho2015gingiva', 'tian2024tongue'],
        'note': 'The tongue is reddest at the tip (a*28.6) and desaturates to the root (a*18.6), and the '
                'coating drops a* by about 12 units. One flat colour on one mesh cannot carry that '
                'gradient. The larger Tian 2024 healthy cohort (n=1448) was rejected: L*69.5 with b*4.0 is '
                'far lighter and less yellow than any contact-spectrophotometer mucosal reading and is '
                'almost certainly a white-balance artefact of its RGB pipeline.'},

    'mucosa_pharyngeal': {
        'label': 'Pharyngeal mucosa',
        'state': 'healthy naso-, oro- and laryngopharyngeal mucosa, in vivo',
        'srgb_hex': '#a2605c',
        'tier': 'synthesized',
        'tier_basis': 'No pharyngeal colorimetry exists. Built from the description of pharyngeal lining '
                      'mucosa as deeper and more saturated red than oral masticatory mucosa, placed near '
                      'L*48 a*26 b*14 -- inside the measured oral-mucosa envelope, on its red-dark side.',
        'sources': _DESC+['yamashiro1996munsell'],
        'note': 'Held inside the Munsell-derived normal oral mucosa envelope (value 3.5-6.0, i.e. roughly '
                'L* 35-60) so it is at least consistent with what was measured for neighbouring mucosa.'},

    'mucosa_nasal': {
        'label': 'Nasal respiratory mucosa',
        'state': 'healthy nasal cavity lining, in vivo',
        'srgb_hex': '#a2605c',
        'tier': 'synthesized',
        'tier_basis': 'No nasal mucosal colorimetry located. Built from the description of respiratory '
                      'mucosa as more vascular and darker red than oral mucosa, at the same point as the '
                      'pharyngeal value because nothing measured separates them.',
        'sources': _DESC,
        'note': 'The nasolacrimal duct carries the same lining and is given the same colour.'},

    'mucosa_gastric': {
        'label': 'Gastric mucosa',
        'state': 'healthy gastric body mucosa at endoscopy',
        'srgb_hex': '#aa6c65',
        'tier': 'synthesized',
        'tier_basis': 'No gastric mucosal colorimetry located. Built from the standard endoscopic '
                      'description of healthy gastric body mucosa as a uniform salmon pink-red, placed at '
                      'about L*52 a*24 b*14.',
        'sources': _DESC},

    'mucosa_glans': {
        'label': 'Glans penis and inner prepuce',
        'state': 'healthy glans, in vivo',
        'lab': [49.67, 27.74, 17.52], 'illuminant_observer': 'D65_2',
        'tier': 'transferred',
        'tier_basis': 'Female external genital mucosa (labia) transferred to the glans. Sommers 2013, 210 '
                      'colposcopic images colour-corrected against a Munsell ColorChecker Mini in a light- '
                      'and temperature-controlled laboratory. Both surfaces are mucocutaneous transitional '
                      'zones of the external genitalia: thin, minimally keratinized epithelium over a '
                      'vascular bed. No colorimetry of the glans exists in any colour system, in vivo or ex '
                      'vivo; the only preputial measurement is an ex-vivo OUTER foreskin L* with a* and b* '
                      'withheld. The illuminant is not stated in the source; D65/2 degree is assumed.',
        'sources': ['sommers2013genital', 'baker2010genital', 'huang2018foreskin', 'sommers2019inguinal'],
        'note': 'One of the three surfaces the owner named. Under the didactic palette it took the flat '
                'reproductive-system colour, which is what made it read wrong. The consistent signal across '
                'all three studies in this domain is that genital mucosa is redder (a* higher by roughly 10 '
                'to 15), less yellow and darker than the adjacent keratinized skin, with the gap narrowing '
                'as constitutive pigmentation rises. Every measurement in this domain is image-derived and '
                'none states an illuminant, so this value is device-relative. It is a female measurement '
                'applied to a male structure and does not co-vary with the skin-tone option.'},

    'mucosa_anal_margin': {
        'label': 'Anal margin',
        'state': 'anal verge and perianal skin, in vivo',
        'srgb_hex': '#7a564f',
        'tier': 'synthesized',
        'tier_basis': 'No anal-margin colorimetry located. Built from the description of the anal margin as '
                      'a pigmented transitional zone, browner and darker than adjacent perianal skin, at '
                      'about L*40 a*14 b*10.',
        'sources': _DESC,
        'note': 'The manifest entity is the topographic anal region, which contains both the pigmented '
                'margin and ordinary perianal skin; the anoderm itself is not meshed, so this colour is '
                'applied to a larger patch than the structure it describes.'},

    'lymphoid_tonsil': {
        'label': 'Palatine tonsil',
        'state': 'healthy tonsil, in vivo',
        'srgb_hex': '#a16864',
        'tier': 'synthesized',
        'tier_basis': 'No tonsillar colorimetry located. Built as mucosa-covered lymphoid tissue: slightly '
                      'deeper and duller than the surrounding oropharyngeal mucosa, about L*50 a*22 b*12.',
        'sources': _DESC},

    # ============================================ integument
    'skin': {
        'label': 'Skin',
        'state': 'constitutive (sun-protected) skin colour',
        'takes_skin_colour': True,
        'tier': 'synthesized',
        'tier_basis': 'Replaced at build time by the selected skin-tone option; the option carries its own '
                      'tier and basis. This placeholder is never emitted.',
        'sources': ['sommers2019inguinal', 'xiao2017skin'],
        'note': 'Constitutive skin colour is a declared parameter, not a default. The palette carries '
                'whichever option the caller selected and names it in the palette record. Site variation '
                'within one person is real and large -- roughly 6 L* and 6 a* units between the inner upper '
                'arm and the forehead -- and is not modelled: one skin value covers the whole integument.'},

    'skin_region': {
        'label': 'Skin, topographic region',
        'state': 'the same skin; these entities are named surface regions of it',
        'takes_skin_colour': True,
        'tier': 'synthesized',
        'tier_basis': 'Takes the selected skin colour by anatomical relation: these 220 entities are named '
                      'topographic regions of the integument, not separate tissues. The tier is that of the '
                      'selected option.',
        'sources': _REL,
        'note': 'Schmalwieser 2024 measured 18 sites and found the inner upper arm lightest and least red '
                'and the forehead, cheek, nape and shoulder darkest and reddest, but publishes its per-site '
                'values only as figures, so no per-region offset could be applied here.'},

    'skin_palmoplantar': {
        'label': 'Palmar and plantar skin',
        'state': 'in vivo palm, thenar eminence',
        'lab': [62.9, 7.2, 16.7], 'illuminant_observer': 'D65_2',
        'tier': 'measured',
        'tier_basis': 'Horibata 2025, thenar eminence, n=67 non-anaemic Japanese adults, Konica Minolta '
                      'CM-700d. The paper states no illuminant or observer; D65/2 degree is assumed.',
        'sources': ['horibata2025sites', 'lu2025issadata'],
        'note': 'A single Japanese cohort, and it deliberately does NOT co-vary with the selected skin '
                'option. That is a known limitation, not an oversight: the palm-to-dorsum lightness contrast '
                'scales with constitutive pigmentation rather than being a fixed offset, so a single global '
                '"palms are lighter" rule is wrong at both ends of the range. No colorimetry of the sole was '
                'found at all, so plantar skin is carrying a palmar value.'},

    'nail_plate': {
        'label': 'Nail plate',
        'state': 'in vivo thumbnail, plate over a perfused nail bed',
        'lab': [55.4, 4.6, 11.1], 'illuminant_observer': 'D65_2',
        'tier': 'measured',
        'tier_basis': 'Horibata 2025, thumbnail, n=67 non-anaemic Japanese adults, Konica Minolta CM-700d. '
                      'The paper states no illuminant or observer; D65/2 degree is assumed.',
        'sources': ['horibata2025sites'],
        'note': 'The only nail colorimetry located. It reads greyer and far less pink than intuition '
                'expects: a*=4.6 against a*=7.2 for the palm of the same subjects. Nail plate and nail bed '
                'have never been measured separately, and this single value covers both.'},

    'hair': {
        'label': 'Hair',
        'state': 'scalp and body hair',
        'srgb_hex': '#3e322c',
        'tier': 'synthesized',
        'tier_basis': 'A single dark-brown placeholder at about L*22 a*4 b*6. No colorimetry was sought.',
        'sources': _DESC,
        'note': 'Hair colour spans a range at least as wide as skin and is deliberately NOT parameterized '
                'here. Treat this as a placeholder, not a claim, and parameterize it the way skin is '
                'parameterized before anyone reads it as the body having brown hair.'},

    # ============================================ musculoskeletal
    'skeletal_muscle': {
        'label': 'Skeletal muscle',
        'state': 'living, perfused muscle at operation',
        'lab': [38.725, 20.15, 7.6], 'illuminant_observer': 'C_2',
        'tier': 'transferred',
        'tier_basis': 'Bovine longissimus dorsi, four-breed mean of Papanikolopoulou 2025 (n=159 carcasses, '
                      'Konica Minolta CR-410, illuminant C / 2 degree, chromatically adapted to D65 here). '
                      'No human skeletal-muscle colorimetry exists, in vivo or cadaveric: intraoperative '
                      'viability work uses ICG fluorescence, Doppler and tissue oximetry and never publishes '
                      'colour coordinates.',
        'sources': ['papanikolopoulou2025beef', 'schelkopf2021beefmethod'],
        'note': 'Bloomed post-mortem beef is oxymyoglobin-rich at the cut surface, drained and pH-shifted. '
                'Living perfused human muscle is darker and less saturated, so this transfer errs bright.'},

    'cardiac_muscle': {
        'label': 'Myocardium',
        'state': 'beating, perfused myocardium',
        'srgb_hex': '#783c39',
        'tier': 'synthesized',
        'tier_basis': 'No myocardial colorimetry located. Built at about L*33 a*26 b*14: darker and browner '
                      'than the transferred skeletal-muscle value, following the standard description of '
                      'myocardium as deeper red-brown than skeletal muscle.',
        'sources': _DESC},

    'valve_leaflet': {
        'label': 'Cardiac valve leaflet',
        'state': 'thin, translucent, pearly',
        'srgb_hex': '#dedacf',
        'tier': 'synthesized',
        'tier_basis': 'Built from the description of a normal valve cusp as thin, translucent and pearly '
                      'white, about L*87 a*-0.5 b*6.',
        'sources': _DESC},

    'tendon': {
        'label': 'Tendon and aponeurosis',
        'state': 'fresh tendon at operation',
        'srgb_hex': '#e2dccf',
        'tier': 'synthesized',
        'tier_basis': 'No tendon colorimetry exists; the spectrocolorimetric programme that covered '
                      'cartilage never covered tendon. Built from the description of tendon as pearly '
                      'glistening white with a fine longitudinal sheen, about L*88 a*-0.5 b*7.',
        'sources': _DESC},

    'ligament': {
        'label': 'Ligament',
        'state': 'fresh ligament at operation',
        'srgb_hex': '#dcd4c3',
        'tier': 'synthesized',
        'tier_basis': 'No ligament colorimetry exists. Built slightly duller and warmer than tendon, about '
                      'L*85 a*0 b*9, following the description of ligament as less specular than tendon.',
        'sources': _DESC},

    'fascia': {
        'label': 'Fascia, septum, fibrous capsule',
        'state': 'fresh fascia at operation',
        'srgb_hex': '#ddd3c2',
        'tier': 'synthesized',
        'tier_basis': 'No fascial or dural colorimetry exists. Built as off-white with a faint warm cast, '
                      'about L*85 a*0.5 b*10, following the description of dura and deep fascia as off-white '
                      'to pearly grey-white.',
        'sources': _DESC},

    'serous_membrane': {
        'label': 'Serous membrane and bursa',
        'state': 'pleura, pericardium, peritoneum, synovial bursa',
        'srgb_hex': '#dedacf',
        'tier': 'synthesized',
        'tier_basis': 'Built as a thin glistening near-transparent membrane, about L*87 a*-0.5 b*6.',
        'sources': _DESC},

    'cortical_bone': {
        'label': 'Cortical bone',
        'state': 'living bone surface at operation',
        'srgb_hex': '#e1cfb2',
        'tier': 'synthesized',
        'tier_basis': 'Built as ivory rather than white, about L*84 a*1.5 b*17.',
        'sources': _DESC,
        'rejected_measurement': {
            'source': 'schafer2001skull',
            'value': 'L*=72.5, a*=-7.4, b*=16.4 (n=124 human autopsy skulls)',
            'converted_would_be': '#b0b594',
            'reason': 'The only human bone colorimetry located, but it measures defleshed autopsy skull. '
                      'Its a*=-7.4 is green-shifted and converts to a green-tan that does not represent '
                      'living cortical bone under periosteum in a bleeding surgical field, which is what '
                      'this palette claims. Rejected on state, not on quality; recorded so the rejection '
                      'is auditable rather than silent.'}},

    'tooth': {
        'label': 'Tooth enamel',
        'state': 'vital unrestored maxillary central incisor, in vivo',
        'lab': [73.5, 2.2, 11.9], 'illuminant_observer': 'D65_2',
        'tier': 'measured',
        'tier_basis': 'Wee 2023, n=120 subjects across four ethnic groups and five age bands, PR-705 '
                      'spectroradiometer, 0/45 geometry, D65 and the CIE 2 degree observer, enamel at '
                      'infinite thickness.',
        'sources': ['wee2023enamel', 'popciutrila2016dentin'],
        'note': 'The flat colour of enamel at infinite thickness is a good deal darker than a tooth looks '
                'in a face, because a real tooth is translucent (translucency parameter 10.1) over lighter '
                'dentin and is seen against a dark oral cavity. A flat-shaded mesh cannot carry that.'},

    'cartilage_hyaline': {
        'label': 'Hyaline and elastic cartilage',
        'state': 'intact articular and laryngeal cartilage',
        'srgb_hex': '#d1d9da',
        'tier': 'synthesized',
        'tier_basis': 'Built as bluish-white and translucent, about L*86 a*-2.5 b*-1.5. Ishimoto 2009 '
                      'measured human cartilage spectrocolorimetrically but published only per-grade '
                      'differences; its intact-grade finding, that reflectance rises monotonically with '
                      'wavelength without the 580 nm dip of degenerate cartilage, is consistent with a '
                      'neutral-to-blue cast and no red component.',
        'sources': _DESC+['ishimoto2009cartilage']},

    'cartilage_fibro': {
        'label': 'Fibrocartilage',
        'state': 'meniscus, labrum, articular disc',
        'srgb_hex': '#d1cbc1',
        'tier': 'synthesized',
        'tier_basis': 'Built denser and warmer than hyaline cartilage, about L*82 a*0 b*6.',
        'sources': _DESC},

    'adipose': {
        'label': 'Adipose tissue',
        'state': 'fresh fat at operation',
        'srgb_hex': '#e9cf95',
        'tier': 'synthesized',
        'tier_basis': 'No human adipose colorimetry exists. Built as buttery cream-yellow, about L*84 a*1 '
                      'b*32. The high-b* direction is mechanistically supported: adipose yellowness tracks '
                      'carotenoid content (r=0.79) in the bovine review, and the same carotenoid reservoir '
                      'is what tints human fat. The magnitude is asserted.',
        'sources': _DESC+['dunne2009fat']},

    # ============================================ vessels, blood, lymph
    'artery_wall': {
        'label': 'Arterial wall',
        'state': 'external (adventitial) surface of an artery at operation',
        'srgb_hex': '#c79d8a',
        'tier': 'synthesized',
        'tier_basis': 'Built as pale pink-tan with a yellowish cast, about L*68 a*13 b*16.',
        'sources': _DESC,
        'note': 'What is visible on an artery is its adventitia, not the blood inside it. The didactic '
                'palette paints every artery a saturated red, which is the colour of arterial blood, not '
                'of an artery.'},

    'vein_wall': {
        'label': 'Venous wall',
        'state': 'external surface of a vein at operation',
        'srgb_hex': '#715768',
        'tier': 'synthesized',
        'tier_basis': 'Built as a dark red-purple, about L*40 a*14 b*-6: the colour of desaturated blood '
                      'seen through a thin wall.',
        'sources': _DESC+['zijlstra1991hb'],
        'note': 'Thin-walled veins take their colour from the blood inside them and read dark red-purple. '
                'Thick-walled veins look much like arteries; one class colour cannot hold both, and this '
                'value represents the thin-walled case, which is the common one in this model. Venous '
                'blood is never blue -- the blue-vein percept is subsurface scattering in skin.'},

    'blood_cavity': {
        'label': 'Blood-filled cavity',
        'state': 'cardiac chamber contents',
        'srgb_hex': '#672129',
        'tier': 'synthesized',
        'tier_basis': 'Built as desaturated venous blood, about L*24 a*32 b*12. The right-sided chambers '
                      'take this and the left-sided chambers take an arterial value; see ENTITY_OVERRIDES.',
        'sources': _DESC+['zijlstra1991hb'],
        'note': 'Zijlstra 1991 gives the human haemoglobin absorption spectra from which arterial and '
                'venous colour can be computed exactly. That integration was not performed here, so only '
                'the ordering -- arterial lighter and more crimson than venous -- rests on it.'},

    'lymph_node': {
        'label': 'Lymph node',
        'state': 'fresh node at operation',
        'srgb_hex': '#c4a58f',
        'tier': 'synthesized',
        'tier_basis': 'Built as tan-grey with a pink cast, about L*70 a*8 b*16.',
        'sources': _DESC},

    'lymph_vessel': {
        'label': 'Lymphatic vessel',
        'state': 'lymphatic carrying clear to straw-coloured lymph',
        'srgb_hex': '#e6ddb7',
        'tier': 'synthesized',
        'tier_basis': 'Built as pale straw, about L*88 a*-3 b*20.',
        'sources': _DESC},

    # ============================================ nervous
    'peripheral_nerve': {
        'label': 'Peripheral nerve',
        'state': 'fresh nerve at operation',
        'srgb_hex': '#e2d3b6',
        'tier': 'synthesized',
        'tier_basis': 'No peripheral-nerve colorimetry exists. Built as pale cream-white, about L*85 a*1 '
                      'b*16, matte rather than the specular white of tendon.',
        'sources': _DESC},

    'brain_grey': {
        'label': 'Grey matter',
        'state': 'fresh unfixed cortex and deep grey nuclei',
        'srgb_hex': '#b18e85',
        'tier': 'synthesized',
        'tier_basis': 'No colorimetry of fresh brain exists. Built as pinkish-tan, about L*62 a*12 b*10.',
        'sources': _DESC,
        'note': 'Fresh cortex is pinkish-tan. The grey of an anatomy specimen is what fixation does to it '
                'and is not what this palette represents.'},

    'brain_white': {
        'label': 'White matter',
        'state': 'fresh unfixed white matter',
        'srgb_hex': '#e8dbc6',
        'tier': 'synthesized',
        'tier_basis': 'Built as creamy white, about L*88 a*1 b*12.',
        'sources': _DESC},

    'choroid_plexus': {
        'label': 'Choroid plexus',
        'state': 'fresh, highly vascular fronds',
        'srgb_hex': '#bd7a78',
        'tier': 'synthesized',
        'tier_basis': 'Built as a vascular pink-red frond, about L*58 a*26 b*12.',
        'sources': _DESC},

    'csf_space': {
        'label': 'Cerebrospinal fluid space',
        'state': 'clear cerebrospinal fluid',
        'srgb_hex': '#d9e4e8',
        'tier': 'synthesized',
        'tier_basis': 'Built as clear with a faint blue cast, about L*90 a*-3 b*-3.',
        'sources': _DESC},

    # ============================================ viscera
    'liver': {
        'label': 'Liver',
        'state': 'perfused liver at operation',
        'srgb_hex': '#854538',
        'tier': 'synthesized',
        'tier_basis': 'No absolute human liver colour has been published, only the direction that steatosis '
                      'raises L* and b*. Built as reddish-brown, about L*37 a*26 b*20. The choice of a dark '
                      'brown-red rather than an arterial red is supported by the measured hyperspectral '
                      'state of a perfused human graft: roughly half-saturated haemoglobin at a high '
                      'haemoglobin index.',
        'sources': _DESC+['kneifel2022liverhsi', 'gomezgavara2024liver', 'piella2024livercolor']},

    'spleen': {
        'label': 'Spleen',
        'state': 'perfused spleen at operation',
        'srgb_hex': '#66323a',
        'tier': 'synthesized',
        'tier_basis': 'No spleen colorimetry of any species was located. Built as dark purple-red, about '
                      'L*28 a*24 b*6: darker and less yellow than liver.',
        'sources': _DESC},

    'kidney': {
        'label': 'Kidney',
        'state': 'capsular surface of a perfused kidney',
        'srgb_hex': '#8d5447',
        'tier': 'synthesized',
        'tier_basis': 'Built as reddish-brown, about L*42 a*22 b*18: lighter and yellower than liver.',
        'sources': _DESC+['shimbashi2024kidney'],
        'note': 'Cortex and medulla cannot be separated here. Each kidney is a single closed surface in '
                'this geometry with no internal partition, so only the capsular surface can be coloured. '
                'The measured cortico-medullary a* difference is in the source list and cannot be applied.'},

    'pancreas': {
        'label': 'Pancreas',
        'state': 'perfused pancreas at operation',
        'srgb_hex': '#cca380',
        'tier': 'synthesized',
        'tier_basis': 'No pancreatic colorimetry located. Built as pale pinkish-tan, about L*70 a*10 b*24.',
        'sources': _DESC},

    'gallbladder': {
        'label': 'Gallbladder',
        'state': 'bile-filled gallbladder seen through its wall',
        'srgb_hex': '#6b6d45',
        'tier': 'synthesized',
        'tier_basis': 'No colorimetry of bile or gallbladder exists; all bilirubin measurement is analytical '
                      'chemistry, not colour. Built as dark olive-green, about L*45 a*-8 b*22.',
        'sources': _DESC},

    'bile_duct': {
        'label': 'Biliary duct',
        'state': 'bile-containing duct',
        'srgb_hex': '#8a8556',
        'tier': 'synthesized',
        'tier_basis': 'Built lighter and yellower than the gallbladder, about L*55 a*-6 b*26, since duct '
                      'bile is less concentrated than gallbladder bile.',
        'sources': _DESC},

    'stomach_serosa': {
        'label': 'Stomach, serosal surface',
        'state': 'external surface at operation',
        'srgb_hex': '#bf9888',
        'tier': 'synthesized',
        'tier_basis': 'Built as pale pink-grey, about L*66 a*12 b*14.',
        'sources': _DESC,
        'note': 'The exposed surface of the stomach entity is serosa, not the red gastric mucosa, which is '
                'a separate entity in this body and is coloured from the mucosal set.'},

    'oesophagus': {
        'label': 'Oesophagus',
        'state': 'external muscular surface',
        'srgb_hex': '#c4a594',
        'tier': 'synthesized',
        'tier_basis': 'Built as pale pinkish-tan, about L*70 a*9 b*13.',
        'sources': _DESC},

    'small_intestine': {
        'label': 'Small intestine',
        'state': 'perfused small bowel at operation',
        'srgb_hex': '#bb8375',
        'tier': 'synthesized',
        'tier_basis': 'Built as well-perfused pink, about L*60 a*20 b*16: pinker than colon.',
        'sources': _DESC},

    'large_intestine': {
        'label': 'Large intestine',
        'state': 'perfused colon at operation',
        'srgb_hex': '#bc8f79',
        'tier': 'synthesized',
        'tier_basis': 'Built greyer and more tan than small bowel, about L*63 a*14 b*18.',
        'sources': _DESC,
        'note': 'Taeniae are paler and appendices epiploicae are yellow; a single class colour cannot carry '
                'that, and the taenia entities take this colour rather than a paler one.'},

    'mesentery_omentum': {
        'label': 'Mesentery and omentum',
        'state': 'peritoneal fold',
        'srgb_hex': '#e2c287',
        'tier': 'synthesized',
        'tier_basis': 'Built as adipose slightly darkened for the vessels running in it, about L*80 a*3 '
                      'b*34, by anatomical relation to the adipose entry.',
        'sources': _REL+_DESC},

    'lung': {
        'label': 'Lung',
        'state': 'ventilated adult lung at operation',
        'srgb_hex': '#b28d8c',
        'tier': 'synthesized',
        'tier_basis': 'No lung colorimetry exists; the anthracosis literature is entirely qualitative. '
                      'Built as adult grey-pink, about L*62 a*14 b*6.',
        'sources': _DESC,
        'note': 'Lung is pink in a child and darkens with age to grey-pink with black anthracotic streaking '
                'along subpleural interlobular septa, near-universal in urban adults and heavier in smokers '
                'and biomass-exposed populations. That is a reticulated pattern, not a uniform darkening; a '
                'flat colour gives the mean tone and not the look.'},

    'airway_wall': {
        'label': 'Airway wall',
        'state': 'trachea and bronchi, external surface',
        'srgb_hex': '#cbaa9b',
        'tier': 'synthesized',
        'tier_basis': 'Built as pale pink over cartilage rings, about L*72 a*10 b*12.',
        'sources': _DESC},

    'urothelium_wall': {
        'label': 'Urinary tract wall',
        'state': 'ureter, bladder and urethra, external surface',
        'srgb_hex': '#c5a592',
        'tier': 'synthesized',
        'tier_basis': 'Built as pale pinkish-tan, about L*70 a*9 b*14.',
        'sources': _DESC},

    'thymus': {
        'label': 'Thymus',
        'state': 'fresh thymus',
        'srgb_hex': '#cca998',
        'tier': 'synthesized',
        'tier_basis': 'Built as pale pinkish-grey and lobulated, about L*72 a*10 b*14.',
        'sources': _DESC},

    # ============================================ glands
    'gland_thyroid': {
        'label': 'Thyroid gland',
        'state': 'perfused thyroid at operation',
        'srgb_hex': '#8d4c42',
        'tier': 'synthesized',
        'tier_basis': 'Built as deep red-brown, about L*40 a*26 b*18, following the description of thyroid '
                      'as one of the most vascular organs in the body and distinctly darker than the strap '
                      'muscles over it.',
        'sources': _DESC},

    'gland_parathyroid': {
        'label': 'Parathyroid gland',
        'state': 'fresh parathyroid',
        'srgb_hex': '#bc8c69',
        'tier': 'synthesized',
        'tier_basis': 'Built as yellow-tan, about L*62 a*14 b*26 -- the tan-against-fat contrast surgeons '
                      'rely on to find it.',
        'sources': _DESC},

    'gland_adrenal': {
        'label': 'Adrenal gland',
        'state': 'fresh adrenal, cortical surface',
        'srgb_hex': '#ddb675',
        'tier': 'synthesized',
        'tier_basis': 'Built as chrome yellow, about L*76 a*6 b*38.',
        'sources': _DESC,
        'note': 'The cortex is strongly yellow and the medulla dark red-brown. This geometry does not '
                'separate them, so only the cortical surface is represented.'},

    'gland_pituitary': {
        'label': 'Pituitary gland',
        'state': 'fresh pituitary',
        'srgb_hex': '#b18e85',
        'tier': 'synthesized',
        'tier_basis': 'Built as pinkish-grey, about L*62 a*12 b*10.',
        'sources': _DESC},

    'gland_pineal': {
        'label': 'Pineal gland',
        'state': 'fresh pineal',
        'srgb_hex': '#a5847e',
        'tier': 'synthesized',
        'tier_basis': 'Built as grey-pink, about L*58 a*12 b*8.',
        'sources': _DESC},

    'gland_salivary': {
        'label': 'Salivary gland',
        'state': 'parotid, submandibular and sublingual glands and their ducts',
        'srgb_hex': '#cda996',
        'tier': 'synthesized',
        'tier_basis': 'Built as pale pinkish-tan and lobulated, about L*72 a*10 b*15.',
        'sources': _DESC},

    'gland_serous': {
        'label': 'Lacrimal gland',
        'state': 'fresh lacrimal gland and drainage apparatus',
        'srgb_hex': '#cda996',
        'tier': 'synthesized',
        'tier_basis': 'Takes the salivary-gland colour by anatomical relation: both are lobulated exocrine '
                      'glands of the same gross appearance.',
        'sources': _REL+_DESC},

    # ============================================ reproductive
    'testis': {
        'label': 'Testis',
        'state': 'tunica albuginea of a perfused testis',
        'srgb_hex': '#d3d2ca',
        'tier': 'synthesized',
        'tier_basis': 'Built as a dense white with a faint blue cast, about L*84 a*-1 b*4.',
        'sources': _DESC,
        'note': 'The visible surface is the tunica albuginea, not the tan seminiferous tissue beneath it.'},

    'epididymis': {
        'label': 'Epididymis',
        'state': 'fresh epididymis',
        'srgb_hex': '#c09884',
        'tier': 'synthesized',
        'tier_basis': 'Built as tan-pink, about L*66 a*12 b*16.',
        'sources': _DESC},

    'duct_muscular': {
        'label': 'Deferent and ejaculatory duct',
        'state': 'fresh, thick-walled white cord',
        'srgb_hex': '#d7d1c4',
        'tier': 'synthesized',
        'tier_basis': 'Built as an off-white muscular cord, about L*84 a*0 b*7.',
        'sources': _DESC},

    'seminal_vesicle': {
        'label': 'Seminal vesicle',
        'state': 'fresh seminal vesicle',
        'srgb_hex': '#be9981',
        'tier': 'synthesized',
        'tier_basis': 'Built as pale tan-pink, about L*66 a*10 b*18.',
        'sources': _DESC},

    'prostate': {
        'label': 'Prostate',
        'state': 'fresh prostate',
        'srgb_hex': '#c4ac9b',
        'tier': 'synthesized',
        'tier_basis': 'Built as firm pale grey-tan, about L*72 a*6 b*12.',
        'sources': _DESC},

    'erectile_tissue': {
        'label': 'Erectile body',
        'state': 'tunica albuginea over cavernous tissue',
        'srgb_hex': '#c5aba2',
        'tier': 'synthesized',
        'tier_basis': 'Built as the white tunica tinged by the dark cavernous tissue showing through it, '
                      'about L*72 a*8 b*8.',
        'sources': _DESC},

    # ============================================ eye
    'sclera': {
        'label': 'Sclera',
        'state': 'bulbar sclera under the conjunctiva',
        'srgb_hex': '#e5d9c6',
        'tier': 'synthesized',
        'tier_basis': 'No usable scleral colorimetry exists. Built at about L*87 a*3 b*13, off-white and '
                      'warm rather than neutral, because both sources that report scleral colour agree it is '
                      'never neutral: a* and b* are both positive, and both rise with age.',
        'sources': _DESC+['russell2014sclera'],
        'note': 'Russell 2014 publishes an 8-bit MATLAB Lab encoding of un-white-balanced portraits, over '
                'the whole visible sclera including lid shadow and vessels; converting it gives an '
                'implausibly dark sclera and it was rejected as a value. The bulbar conjunctiva that '
                'overlies the sclera is not a separate entity in this body, so the conjunctiva could not be '
                'coloured from the mucosal set even though a cohort measurement of palpebral conjunctiva '
                'exists. See UNMAPPED.'},

    'cornea': {
        'label': 'Cornea',
        'state': 'transparent cornea',
        'srgb_hex': '#d0d8dd',
        'tier': 'synthesized',
        'tier_basis': 'Built as a near-neutral pale blue-grey standing in for transparency, about L*86 '
                      'a*-2 b*-3. An opaque colour for a transparent tissue is a representational choice, '
                      'not a claim.',
        'sources': _DESC},

    'iris': {
        'label': 'Iris',
        'state': 'iris stroma',
        'srgb_hex': '#694f3b',
        'tier': 'synthesized',
        'tier_basis': 'A mid-brown placeholder at about L*36 a*8 b*16.',
        'sources': _DESC,
        'note': 'Iris colour is as variable between people as skin and is deliberately NOT parameterized '
                'here. Brown is the most common iris colour worldwide; it is not a claim about this body, '
                'and it should be turned into an option set alongside skin.'},

    'lens': {
        'label': 'Lens',
        'state': 'adult crystalline lens',
        'srgb_hex': '#e4dcc6',
        'tier': 'synthesized',
        'tier_basis': 'Built as near-clear with the faint yellow of an adult lens, about L*88 a*-1 b*12.',
        'sources': _DESC},

    'ocular_humour': {
        'label': 'Aqueous and vitreous',
        'state': 'clear ocular fluid',
        'srgb_hex': '#e3e9ea',
        'tier': 'synthesized',
        'tier_basis': 'Built as near-clear, about L*92 a*-2 b*-1.',
        'sources': _DESC},

    'ciliary_body': {
        'label': 'Ciliary body',
        'state': 'pigmented ciliary body',
        'srgb_hex': '#7b554b',
        'tier': 'synthesized',
        'tier_basis': 'Built as dark pigmented brown-pink, about L*40 a*14 b*12.',
        'sources': _DESC},

    'choroid': {
        'label': 'Choroid',
        'state': 'highly vascular pigmented choroid',
        'srgb_hex': '#6b3835',
        'tier': 'synthesized',
        'tier_basis': 'Built as dark red-brown, about L*30 a*22 b*12.',
        'sources': _DESC},

    'retina': {
        'label': 'Retina',
        'state': 'fundus appearance in vivo',
        'srgb_hex': '#ab644c',
        'tier': 'synthesized',
        'tier_basis': 'Built as the orange-red of the fundus, about L*50 a*26 b*26.',
        'sources': _DESC,
        'note': 'The neural retina is nearly transparent. What is seen through it is the retinal pigment '
                'epithelium and choroid, and that is what is coloured here.'},

    # ============================================ not tissue
    'field_overlay': {
        'label': 'Computed field overlay',
        'state': 'not a tissue surface',
        'srgb_hex': '#9aa2a6',
        'tier': 'synthesized',
        'tier_basis': 'Deliberately neutral, about L*66 a*-2 b*-3, so a computed overlay is not read as a '
                      'tissue colour in either palette.',
        'sources': _REL,
        'note': 'A computed field overlay, not anatomy.'},
}


# ---------------------------------------------------------------------------- skin as a parameter
#
# Constitutive skin colour varies enormously between people and this model does not hard-code one.
# The option set is declared, each option states its basis, and the build says which one it used.
#
# The axis is the Individual Typology Angle, ITA = arctan((L*-50)/b*) x 180/pi, whose six-group
# classification -- very light > 55, light > 41, intermediate > 28, tan > 10, brown > -30,
# dark <= -30 degrees -- was read verbatim out of an open-access source and traced to Del Bino 2018.
# Only L* and b* enter ITA; a* is unused.
#
# There are three kinds of option and they are tiered differently.
#
#   measured_*   real cohort means. The sommers_* options are the strongest: one instrument, one
#                sun-protected site (inner upper thigh below the groin), three groups measured under
#                the same protocol, n=332. The xiao_* options are inner forearm, a second instrument
#                family with a declared D65 / 2 degree observer, four cohorts, n=960.
#   ita_*        constructed points on the published ITA axis, for callers who want the full range
#                rather than the groups any one study happened to enrol. b*=17 and a*=10 are asserted
#                as typical of skin and L* is solved from L* = 50 + b* tan(ITA). No cohort table of
#                mean L*a*b* per ITA class has ever been published, so these are synthesized.
#   not_applied  the build default. ITA = 0 exactly, i.e. L* = 50. Not a claim about this body; it is
#                what the palette shows when no option has been chosen.
#
# On the group labels: where an option carries one, it is the descriptor the source study used for
# its own cohort, reproduced verbatim in cohort_as_reported. Those are social categories as recorded
# by that study, not colour classes, and nothing here maps one onto the other. The largest study in
# this literature reports that 89.4 % of individuals have a perceptually indistinguishable
# counterpart in another group and that the median gamut overlap between groups is 60.5 %, which is
# the reason skin colour is a parameter here and not a category.
#
# Two further facts that this option set does NOT model, recorded so they are not mistaken for
# oversights: constitutive colour is site-dependent by roughly 6 L* units within one person, and
# "sun-protected is lightest" is not universal -- in one cohort of African American participants
# aged 18-30 the sun-protected buttock was darker than the sun-exposed forearm, the reverse of the
# Caucasian pattern in the same study.

import math as _math


def _ita(L, b):
    return round(_math.degrees(_math.atan((L-50.0)/b)), 1)


def _measured_skin(key, label, lab, source, cohort, instrument, illuminant, assumed, note):
    return {'label': label,
            'description': note,
            'cohort_as_reported': cohort,
            'ita_degrees': _ita(lab[0], lab[2]),
            'ita_class_boundaries': 'very light >55, light >41, intermediate >28, tan >10, brown >-30, '
                                    'dark <=-30',
            'basis': 'Measured cohort mean. '+instrument+'. '+
                     ('Illuminant and standard observer stated as '+illuminant+'.' if not assumed else
                      'The source states no illuminant or observer; '+illuminant+' is assumed for the '
                      'conversion and the assumption is carried in the palette record.'),
            'colour': {'lab': lab, 'illuminant_observer': illuminant,
                       'tier': 'measured',
                       'tier_basis': 'cohort mean CIE L*a*b* of constitutive human skin measured by '
                                     'reflectance spectrophotometry; see basis',
                       'sources': [source]}}


def _ita_skin(ita, label, description):
    b, a = 17.0, 10.0
    return {'label': label, 'description': description, 'cohort_as_reported': None,
            'ita_degrees': ita,
            'ita_class_boundaries': 'very light >55, light >41, intermediate >28, tan >10, brown >-30, '
                                    'dark <=-30',
            'basis': 'Constructed, not measured. The ITA class boundaries are published and were verified; '
                     'the colour is built by asserting b*=17 and a*=10 as typical of skin and solving '
                     'L* = 50 + b* tan(ITA) at the representative ITA shown. No measured group-mean L*a*b* '
                     'per ITA class was retrievable in this build, so this is synthesized, not a '
                     'measurement, and a later build can replace it without changing the interface.',
            'colour': {'lab': [round(50+b*_math.tan(_math.radians(ita)), 3), a, b],
                       'illuminant_observer': 'D65_2',
                       'tier': 'synthesized',
                       'tier_basis': 'constructed from the verified ITA classification with two asserted '
                                     'chromatic coordinates; see basis',
                       'sources': ['itoformula2026', 'delbino2018ita', 'chardon1991ita',
                                   'delbino2013variation']}}


_SOMMERS_INSTR = ('ColorTec PSM hand-held reflectance spectrophotometer, right inner upper thigh two inches '
                  'below the groin, a sun-protected site, n=341 women, Philadelphia USA and San Juan')
_XIAO_INSTR = ('Konica Minolta CM-2600d or X-Rite SP62, d/8 geometry with the specular component included, '
               'inner forearm, n=960 across four countries')

SKIN_TONES = {
    'not_applied': {
        'label': 'Not applied (neutral)',
        'description': 'The build default. ITA = 0 degrees exactly, i.e. L* = 50, on the boundary between '
                       'the tan and brown classes. It is not a claim that this body has this skin colour; '
                       'it is what the palette shows when no option has been chosen.',
        'cohort_as_reported': None,
        'ita_degrees': 0.0,
        'ita_class_boundaries': 'very light >55, light >41, intermediate >28, tan >10, brown >-30, dark <=-30',
        'basis': 'ITA = 0 by construction, with the same asserted a* and b* as the constructed options. '
                 'Chosen so the shipped palette does not silently pick a pigmentation.',
        'colour': {'lab': [50.0, 10.0, 17.0], 'illuminant_observer': 'D65_2', 'tier': 'synthesized',
                   'tier_basis': 'neutral placeholder on the ITA axis; no measurement',
                   'sources': ['itoformula2026']}},

    'measured_inguinal_a': _measured_skin(
        'measured_inguinal_a', 'Measured, sun-protected inguinal skin, L* 64.4',
        [64.39, 7.71, 18.12], 'sommers2019inguinal', 'non-Hispanic White (n=88)', _SOMMERS_INSTR, 'D65_2',
        True, 'Lightest of the three groups measured in one protocol at one sun-protected site.'),
    'measured_inguinal_b': _measured_skin(
        'measured_inguinal_b', 'Measured, sun-protected inguinal skin, L* 55.9',
        [55.93, 9.15, 20.26], 'sommers2019inguinal', 'Hispanic/Latina (n=190)', _SOMMERS_INSTR, 'D65_2',
        True, 'Middle of the three groups measured in one protocol at one sun-protected site.'),
    'measured_inguinal_c': _measured_skin(
        'measured_inguinal_c', 'Measured, sun-protected inguinal skin, L* 41.1',
        [41.05, 10.16, 19.40], 'sommers2019inguinal', 'non-Hispanic Black (n=54)', _SOMMERS_INSTR, 'D65_2',
        True, 'Darkest of the three groups measured in one protocol at one sun-protected site.'),

    'measured_forearm_a': _measured_skin(
        'measured_forearm_a', 'Measured, inner forearm, L* 63.0',
        [63.0, 5.6, 14.0], 'xiao2017skin', 'Caucasian, UK (n=187)', _XIAO_INSTR, 'D65_2', False,
        'Inner forearm mean. Note that the same study measures the forehead about 4 L* darker and 6 a* '
        'redder in every cohort, so a whole-body single value is a simplification.'),
    'measured_forearm_b': _measured_skin(
        'measured_forearm_b', 'Measured, inner forearm, L* 60.9',
        [60.9, 7.0, 15.0], 'xiao2017skin', 'Chinese (n=202)', _XIAO_INSTR, 'D65_2', False,
        'Inner forearm mean.'),
    'measured_forearm_c': _measured_skin(
        'measured_forearm_c', 'Measured, inner forearm, L* 60.6',
        [60.6, 6.5, 16.4], 'xiao2017skin', 'Kurdish, Iraq (n=145)', _XIAO_INSTR, 'D65_2', False,
        'Inner forearm mean.'),
    'measured_forearm_d': _measured_skin(
        'measured_forearm_d', 'Measured, inner forearm, L* 61.9',
        [61.9, 7.1, 17.4], 'xiao2017skin', 'Thai (n=426)', _XIAO_INSTR, 'D65_2', False,
        'Inner forearm mean. Highest b* of the four cohorts; yellowness is the coordinate that separates '
        'groups in this study, while redness is nearly constant across them.'),

    'ita_very_light': _ita_skin(65.0, 'Constructed, very light (ITA > 55)',
                                'Representative point inside the published very-light class.'),
    'ita_light': _ita_skin(48.0, 'Constructed, light (ITA 41 to 55)',
                           'Representative point inside the published light class.'),
    'ita_intermediate': _ita_skin(34.5, 'Constructed, intermediate (ITA 28 to 41)',
                                  'Representative point inside the published intermediate class.'),
    'ita_tan': _ita_skin(19.0, 'Constructed, tan (ITA 10 to 28)',
                         'Representative point inside the published tan class.'),
    'ita_brown': _ita_skin(-10.0, 'Constructed, brown (ITA -30 to 10)',
                           'Representative point inside the published brown class.'),
    'ita_dark': _ita_skin(-45.0, 'Constructed, dark (ITA <= -30)',
                          'Representative point inside the published dark class. The class is open-ended '
                          'below, so -45 degrees is a stated choice, not a midpoint.'),
}


# ---------------------------------------------------------------------------- per-entity overrides
#
# Keyed by structure id or by normalized structure name. Used where the tissue class is too coarse
# to be honest.

ENTITY_OVERRIDES = {
    'cavity of left atrium': {
        'label': 'Left atrial blood', 'srgb_hex': '#8f2222', 'tier': 'synthesized',
        'tier_basis': 'Oxygenated blood, about L*32 a*45 b*28. Left-sided chambers carry blood at 97-100 % '
                      'saturation and are lighter and more crimson than the right; that ordering is what '
                      'the human haemoglobin spectra give, though the spectra were not integrated here.',
        'sources': ['zijlstra1991hb', 'gross_description']},
    'cavity of left ventricle': {
        'label': 'Left ventricular blood', 'srgb_hex': '#8f2222', 'tier': 'synthesized',
        'tier_basis': 'Oxygenated blood; see cavity of left atrium.',
        'sources': ['zijlstra1991hb', 'gross_description']},
    'cavity of right atrium': {
        'label': 'Right atrial blood', 'srgb_hex': '#672129', 'tier': 'synthesized',
        'tier_basis': 'Desaturated blood at about 70-75 % saturation, darker and more burgundy than the '
                      'left side. Never blue.',
        'sources': ['zijlstra1991hb', 'gross_description']},
    'cavity of right ventricle': {
        'label': 'Right ventricular blood', 'srgb_hex': '#672129', 'tier': 'synthesized',
        'tier_basis': 'Desaturated blood; see cavity of right atrium.',
        'sources': ['zijlstra1991hb', 'gross_description']},
}


# ---------------------------------------------------------------------------- mucosal targets
#
# The mucosal and transitional surfaces named in the request, and what happened to each. A target
# with structure_ids is coloured from the mucosal set; a target with none could not be coloured
# and says why.

UNMAPPED = [
    {'target': 'areola and nipple',
     'status': 'no entity',
     'reason': 'There is no areola, nipple, mammary gland or breast entity anywhere in the 7390 display '
               'structures. The nearest is z-anatomy "Mammary region" (left and right), a topographic skin '
               'region covering the whole breast surface; colouring that from the mucosal set would recolour '
               'the entire breast, which would be worse than leaving it as skin. The underlying atlas is an '
               'adult male reference specimen.',
     'nearest_entities': ['za-* Mammary region.l', 'za-* Mammary region.r'],
     'would_need': 'an areola/nipple surface patch on the outer envelope, or a segmented mammary entity'},

    {'target': 'inner prepuce (foreskin)',
     'status': 'no entity',
     'reason': 'No prepuce, foreskin or preputial entity exists. Only the glans, corpus cavernosum and '
               'corpus spongiosum are meshed. The glans itself IS coloured from the mucosal set.',
     'nearest_entities': ['bp3d-FJ3134 glans penis', 'body-bp3d-FJ3134', 'za-29295a3c5ede2af3 Glans penis'],
     'would_need': 'a preputial surface entity'},

    {'target': 'labia minora, vulvar vestibule, vaginal and cervical mucosa',
     'status': 'no entity',
     'reason': 'The reproductive set in this body is male only: testis, epididymis, deferent duct, seminal '
               'vesicle, prostate, ejaculatory duct and the three penile bodies. There is no vulva, vagina, '
               'cervix, uterus or ovary entity at all, so none of these mucosal surfaces can be coloured.',
     'nearest_entities': [],
     'would_need': 'a female reproductive tract in the display manifest'},

    {'target': 'bulbar and palpebral conjunctiva',
     'status': 'no entity',
     'reason': 'No conjunctival entity exists. The sclera is present and is the structure the bulbar '
               'conjunctiva overlies, but sclera is not mucosa and colouring it from the mucosal set would '
               'be wrong. This is the one unmapped target for which a cohort measurement DOES exist -- '
               'Horibata 2025 measured lower palpebral conjunctiva at L*50.0 a*36.7 b*21.5 in 67 '
               'non-anaemic adults -- so the gap here is geometry, not evidence. Even with an entity, a '
               'flat colour would be a poor representation: conjunctival redness is driven by vessel '
               'coverage (R2=0.93) far more than by any tissue colour variable (best R2=0.62).',
     'nearest_entities': ['left sclera', 'right sclera'],
     'would_need': 'a conjunctival surface entity over the sclera, and ideally a vessel layer rather than '
                   'a tint'},

    {'target': 'buccal mucosa',
     'status': 'no entity',
     'reason': 'No buccal mucosa entity exists. "Buccal region" (left and right) is the EXTERNAL cheek skin '
               'region and the buccinator is muscle; neither is the wet inner cheek. No colour class is '
               'defined for buccal mucosa, because a class with no structure to carry it would be dead data.',
     'nearest_entities': ['za-* Buccal region.l', 'za-* Buccal region.r'],
     'would_need': 'an intraoral cheek lining surface'},

    {'target': 'labial mucosa (wet inner lip)',
     'status': 'no entity',
     'reason': 'Not separately meshed. The single "lip" entity spans the cutaneous lip, the vermilion and '
               'the wet labial mucosa as one surface, and takes the vermilion colour.',
     'nearest_entities': ['bp3d-FJ2814 lip', 'body-bp3d-FJ2814'],
     'would_need': 'a lip split into cutaneous, vermilion and mucosal zones'},

    {'target': 'anal canal lining (anoderm and transitional zone)',
     'status': 'partially mapped',
     'reason': 'The anal margin is covered by the topographic "Anal region" entity, which is coloured from '
               'the mucosal set, but that region includes ordinary perianal skin as well. The anal canal '
               'lining itself is not meshed; the rectum entity is the tube seen from outside.',
     'nearest_entities': ['za-* Anal region.l', 'za-* Anal region.r', 'rectum'],
     'would_need': 'an anal canal luminal surface'},

    {'target': 'hard palate mucosa',
     'status': 'no entity',
     'reason': 'Only the soft palate and the uvula are meshed. Both are coloured from the mucosal set. No '
               'hard palate mucosal surface exists.',
     'nearest_entities': ['za-ad5c00cdd8042733 Soft palate', 'za-4b5c6f9187019cbf Uvula of palate'],
     'would_need': 'a hard palate lining surface'},

    {'target': 'renal cortex and medulla',
     'status': 'not a mucosal target, recorded because it was asked for',
     'reason': 'Each kidney is a single closed capsular surface with no internal partition, so the '
               'cortico-medullary colour difference -- which is the one internal contrast in this list that '
               'HAS been measured in humans -- cannot be shown. Only the capsular surface is coloured.',
     'nearest_entities': ['bp3d-* left kidney', 'bp3d-* right kidney'],
     'would_need': 'a kidney segmented into cortex and medulla'},

    {'target': 'adrenal cortex and medulla',
     'status': 'not a mucosal target, recorded because it was asked for',
     'reason': 'Single surface per gland. The cortex is strongly yellow and the medulla dark red-brown; '
               'only the cortical surface is represented.',
     'nearest_entities': ['left adrenal gland', 'right adrenal gland'],
     'would_need': 'an adrenal segmented into cortex and medulla'},
]


# Tissues for which a colorimetric measurement was searched for and NOT found. Recorded so the
# absence is a result rather than an omission.
UNMEASURED_SEARCHED = [
    # viscera and connective tissue
    'spleen (no colour measurement of any species located)',
    'lung parenchyma (the anthracosis literature is entirely qualitative)',
    'pancreas (graft assessment uses perfusion and histology, never colour)',
    'tendon and ligament (the spectrocolorimetric programme covered cartilage only)',
    'peripheral nerve',
    'dura mater and meninges',
    'brain grey and white matter, fresh (scattering data exists; colour does not)',
    'bile (all bilirubin measurement is analytical chemistry, not colour of the fluid)',
    'human skeletal muscle, in vivo or cadaveric (viability work uses ICG, Doppler and oximetry)',
    'human adipose tissue (only bovine)',
    'absolute human liver colour (only the relative steatosis direction is published)',
    # mucosa
    'palatal, buccal, labial, pharyngeal, nasal and gastric mucosa (absolute CIELAB)',
    'hard palate mucosa',
    # anogenital: searched by name, nothing found in any colour system
    'glans penis (nothing, in any colour system, in vivo or ex vivo)',
    'inner prepuce / preputial mucosa (the one prepuce measurement is the OUTER surface, L* only)',
    'labia minora reported separately (acquired and withheld in two studies)',
    'labia majora in CIELAB (only n=10 Mexameter arbitrary units)',
    'vulvar vestibule',
    'anal margin, perianal skin, anoderm and anal verge (a dedicated search returned zero records)',
    'scrotum',
    'Individual Typology Angle for any genital site',
    # integument and appendages
    'areola and nipple L*a*b* or ITA (a melanin-index ratio and Munsell notations exist, unretrieved)',
    'plantar / sole skin CIELAB',
    'palmar erythema CIELAB',
    'nail bed and nail plate measured separately',
    # eye
    'bulbar conjunctiva absolute CIELAB (only CIE u* and HSV indices, no absolutes)',
    'scleral b* in CIELAB, jaundiced or normal (the field publishes xy, XYZ and YCbCr only)',
    'limbus and caruncle',
    # classification
    'mean L*a*b* per Individual Typology Angle class (the classes are (L*,b*) regions only)',
]
