"""Literature behind every hair-field number, with the identifier check that verified it.

Tiers. measured: the cited cohort measured this quantity at this body site.
transferred: measured at another site, sex or population and carried here with the
donor named. assumed: an engineering value with no cohort behind it. A quantity
with no source stays None; no field is filled with a plausible number.

Every pmid was resolved through eutils esummary and every doi through
api.crossref.org on 2026-09-07. Two identifiers recalled from memory were WRONG and
are recorded as corrections below, which is why nothing here is unverified.
"""

SOURCES = {
    'jonsson2017': {
        'citation': 'Jonsson EH, Bendas J, Weidner K, Wessberg J, Olausson H, Backlund Wasling H, Croy I. The '
                    'relation between human hair follicle density and touch perception. Sci Rep 2017;7(1):2499.',
        'doi': '10.1038/s41598-017-02308-9', 'pmid': '28566678', 'pmcid': 'PMC5451466', 'species': 'human',
        'method': 'cyanoacrylate skin surface stripping over a marked shaved 2x2 cm area, follicles counted per '
                  'cm2 under a light microscope',
        'cohorts': {'pretest_nine_sites': 'n=15, 8 female, mean age 23.9 +/- 3.0 y, range 20-32 y',
                    'forearm_combined': 'n=138; study 1 n=58 (34 female) mean age 26.2 +/- 6.3 y, study 2 n=80 '
                                        '(51 female) mean age 24.9 +/- 4.1 y'},
        'ancestry': 'not reported by the authors; ethics approval Dresden University of Technology, recruitment '
                    'Dresden and Gothenburg',
        'verification': 'fulltext, Table 2 and Methods, read 2026-09-07 via Europe PMC PMC5451466',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'The relation between human hair follicle density and touch '
                                               'perception'},
        'values_read_here': {'unit': 'hair follicles per cm2, mean and SD, pretest n=15',
                             'forehead': [285.0, 84.1], 'neck': [47.3, 17.1], 'chest': [28.4, 6.4],
                             'upper arm': [45.9, 14.5], 'forearm': [40.3, 14.9], 'back': [24.7, 5.9],
                             'abdomen': [16.9, 6.2], 'thigh': [21.0, 8.0], 'calf': [15.6, 3.6],
                             'forearm_n138': [39.2, 12.1], 'forearm_n138_women': [41.4, 12.8],
                             'forearm_n138_men': [35.5, 10.0]},
        'note': 'the Methods name the trunk sites lower back and lower abdomen; Table 2 labels them back and '
                'abdomen. ihm/assembly/hair.py carries exactly these nine numbers as unattributed priors; this '
                'record is the source they came from.',
    },
    'otberg2004': {
        'citation': 'Otberg N, Richter H, Schaefer H, Blume-Peytavi U, Sterry W, Lademann J. Variations of hair '
                    'follicle size and distribution in different body sites. J Invest Dermatol 2004;122(1):14-19.',
        'doi': '10.1046/j.0022-202X.2003.22110.x', 'pmid': '14962084', 'species': 'human',
        'method': 'cyanoacrylate skin surface biopsy at lateral forehead, back, thorax, upper arm, forearm, thigh '
                  'and calf',
        'cohorts': {}, 'ancestry': 'Caucasian, per the same group in mangelsdorf2006',
        'verification': 'abstract only (eutils efetch 2026-09-07). The table was NOT obtained: jidonline.org and '
                        'sciencedirect.com both return HTTP 403 and no PMC copy exists.',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Variations of hair follicle size and distribution in different '
                                               'body sites',
                             'correction': 'PMID 14962082, recalled from memory for this paper, is in fact Leask A, '
                                           'Denton CP, Abraham DJ, Insights into the molecular mechanism of chronic '
                                           'fibrosis, J Invest Dermatol 2004;122(1):1-6. The adjacent article. '
                                           '14962084 is correct.'},
        'used_for': 'site ordering and method corroboration. Two of its numbers, forehead 292 and back 29 '
                    'follicles per cm2, were read second-hand in todo2017; they are recorded as comparanda and no '
                    'field takes its value from them.',
        'rejected_transcriptions': 'a 2026 burns paper (PMID 42494290) tabulates scalp, hand, foot, perineum and '
                                   'gluteal densities attributed to this paper. Otberg measured seven sites and '
                                   'none of those five, and its back value there (58) contradicts the 29 read in '
                                   'todo2017. Those numbers are not used.',
    },
    'mangelsdorf2006': {
        'citation': 'Mangelsdorf S, Otberg N, Maibach HI, Sinkgraven R, Sterry W, Lademann J. Ethnic variation in '
                    'vellus hair follicle size and distribution. Skin Pharmacol Physiol 2006;19(3):159-167.',
        'doi': '10.1159/000093050', 'pmid': '16679817', 'species': 'human',
        'method': 'cyanoacrylate skin surface biopsy at seven body sites',
        'cohorts': {'groups': 'Asians and African-Americans, compared with the Caucasian cohort of otberg2004'},
        'ancestry': 'the authors terms are Asians, African-Americans, Caucasians, Whites',
        'verification': 'abstract only (eutils efetch 2026-09-07); Karger full text returns HTTP 403',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Ethnic Variation in Vellus Hair Follicle Size and Distribution'},
        'read_here': 'follicular density on the forehead is significantly lower in Asians and African-Americans '
                     'than in Whites; smaller volume, surface, follicular orifice and hair shaft diameter on the '
                     'thigh and calf in Asians and African-Americans; follicular volume generally higher in '
                     'Caucasians.',
        'used_for': 'the existence and direction of ancestry variation, in the authors own terms. No number is '
                    'carried: the abstract gives significance and direction, not the table.',
    },
    'alsharif2022': {
        'citation': 'Alsharif SH, AlGhamdi KM. Evaluation of Scalp Hair Density and Diameter in the Arab '
                    'Population: Clinical Office-Based Phototrichogram Analysis. Clin Cosmet Investig Dermatol '
                    '2022;15:2737-2743.',
        'doi': '10.2147/CCID.S394045', 'pmid': '36545499', 'pmcid': 'PMC9762255', 'species': 'human',
        'method': 'clinical office-based phototrichogram (TrichoSciencePro, x25) at 12, 24 and 30 cm from the '
                  'glabella on the midline, defining frontal, vertex and occipital areas',
        'cohorts': {'scalp': 'n=120, 60 male and 60 female, healthy adults aged 18-60 y, Riyadh 2021-2022'},
        'ancestry': 'Arabian adults, in the authors own term',
        'verification': 'fulltext, Table 2 text read 2026-09-07 via Europe PMC PMC9762255; also carried in-repo at '
                        'data/measurements/hair/strand_reference.json scalp_geometry',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Evaluation of Scalp Hair Density and Diameter in the Arab '
                                               'Population: Clinical Office-Based Phototrichogram Analysis',
                             'correction': 'the in-repo card gives no page range; the record is 15:2737-2743, not '
                                           'the 2857-2865 first written here.'},
        'values_read_here': {'unit': 'hairs per cm2 and shaft diameter um, mean and SD',
                             'frontal': {'density_cm2': [143.9, 5.7], 'diameter_um': [83.5, 2.5]},
                             'vertex': {'density_cm2': [147.1, 7.8], 'diameter_um': [87.0, 4.9]},
                             'occipital': {'density_cm2': [153.6, 6.5], 'diameter_um': [90.7, 5.8]},
                             'sex': 'occipital male 152.3 +/- 5.0 against female 154.8 +/- 7.6, p=0.038; no other '
                                    'site differed by sex'},
        'note': 'a phototrichogram counts visible hairs, not follicles. Telogen follicles carrying no shaft are '
                'not in this number, so it is a lower bound on follicle density and the right number for mass.',
    },
    'leerunyakul2020': {
        'citation': 'Leerunyakul K, Suchonwanit P. Evaluation of Hair Density and Hair Diameter in the Adult Thai '
                    'Population Using Quantitative Trichoscopic Analysis. Biomed Res Int 2020;2020:2476890.',
        'doi': '10.1155/2020/2476890', 'pmid': '32104683', 'pmcid': 'PMC7035527', 'species': 'human',
        'method': 'quantitative trichoscopic analysis',
        'cohorts': {'scalp': 'n=239, 79 men and 160 women, mean age 37.9 y'},
        'ancestry': 'adult Thai population, in the authors own term',
        'verification': 'fulltext Table 1 via Europe PMC PMC7035527 (subagent read, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils', 'status': 'resolved',
                             'returned_title': 'Evaluation of Hair Density and Hair Diameter in the Adult Thai '
                                               'Population Using Quantitative Trichoscopic Analysis'},
        'values_read_here': {'unit': 'hairs per cm2 and shaft diameter um, mean and SD',
                             'frontal': {'density_cm2': [154.3, 12.8], 'diameter_um': [81.1, 6.9]},
                             'vertex': {'density_cm2': [162.9, 15.7], 'diameter_um': [80.8, 5.9]},
                             'temporoparietal': {'density_cm2': [133.7, 14.6], 'diameter_um': [80.1, 4.5]},
                             'occipital': {'density_cm2': [160.2, 15.3], 'diameter_um': [80.3, 5.1]}},
        'used_for': 'the temporal scalp, which alsharif2022 did not measure, and as an independent population '
                    'comparison for the scalp density this candidate carries.',
    },
    'birnbaum2018': {
        'citation': 'Birnbaum MR, McLellan BN, Shapiro J, Ye K, Reid SD. Evaluation of Hair Density in Different '
                    'Ethnicities in a Healthy American Population Using Quantitative Trichoscopic Analysis. Skin '
                    'Appendage Disord 2018;4(4):304-307.',
        'doi': '10.1159/000485522', 'pmid': '30410902', 'species': 'human',
        'method': 'quantitative trichoscopy; participants self-identified their ethnicity',
        'cohorts': {'scalp': 'n=166; 99 Americans of Hispanic descent, 44 of African descent, 23 Caucasians; '
                             'recruited at a New York City medical center'},
        'ancestry': 'Americans of Hispanic descent, individuals of African descent, Caucasians, self-identified',
        'verification': 'abstract (eutils efetch 2026-09-07)',
        'identifier_check': {'resolver': 'eutils', 'status': 'resolved',
                             'returned_title': 'Evaluation of Hair Density in Different Ethnicities in a Healthy '
                                               'American Population Using Quantitative Trichoscopic Analysis'},
        'values_read_here': {'unit': 'hairs per cm2, range of the per-site means',
                             'hispanic_descent': [169, 178], 'african_descent': [148, 160], 'caucasian': [214, 230]},
        'used_for': 'the ancestry spread of scalp hair density, stated in the authors own terms. It bounds how far '
                    'the single scalp density carried here can be from another population: 148 to 230 per cm2 '
                    'across the three groups, against the 148.2 this candidate uses.',
    },
    'jimenez1999': {
        'citation': 'Jimenez F, Ruifernandez JM. Distribution of human hair in follicular units. A mathematical '
                    'model for estimating the donor size in follicular unit transplantation. Dermatol Surg '
                    '1999;25(4):294-298.',
        'doi': '10.1046/j.1524-4725.1999.08114.x', 'pmid': '10417585', 'species': 'human',
        'method': 'digital photography of the occipital donor scalp',
        'cohorts': {'scalp': 'n=50 patients; sex, age and ancestry not stated in the abstract; clinic in Las '
                             'Palmas, Spain'},
        'ancestry': 'not stated',
        'verification': 'abstract (eutils efetch 2026-09-07)',
        'identifier_check': {'resolver': 'eutils', 'status': 'resolved',
                             'returned_title': 'Distribution of human hair in follicular units. A mathematical '
                                               'model for estimating the donor size in follicular unit '
                                               'transplantation'},
        'values_read_here': {'occipital_follicular_units_per_cm2': [65, 85],
                             'occipital_hair_density_per_cm2': [124, 200],
                             'inter_unit_distance_mm': [1.0, 1.4]},
        'used_for': 'the follicular-unit structure of the scalp and an independent bracket on occipital hair '
                    'density. The candidate carries hairs, not follicular units.',
    },
    'blume1991': {
        'citation': 'Blume U, Ferracin J, Verschoore M, Czernielewski JM, Schaefer H. Physiology of the vellus '
                    'hair follicle: hair growth and sebum excretion. Br J Dermatol 1991;124(1):21-28.',
        'doi': '10.1111/j.1365-2133.1991.tb03277.x', 'pmid': '1993141', 'species': 'human',
        'method': 'computerised image analysis of in vivo photographs, forehead, cheek, chest, shoulder and back',
        'cohorts': {'vellus': 'healthy men and women aged 15-30 y; n not stated in the abstract'},
        'ancestry': 'not stated; CIRD, Valbonne, France',
        'verification': 'abstract (eutils efetch 2026-09-07), first hand',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Physiology of the vellus hair follicle: hair growth and sebum '
                                               'excretion'},
        'values_read_here': {'forehead_hairs_per_cm2': 439, 'forehead_growing_fraction': .49,
                             'back_hairs_per_cm2': 85, 'back_growing_fraction': .315,
                             'growth_rate_mm_per_day': [.03, .13],
                             'growth_rate_note': 'from 0.03 mm/day on the forehead to 0.13 mm/day on the back'},
        'second_hand_values': {'source': 'todo2017 quotes this paper for per-sex site means',
                               'unit': 'hairs per cm2, mean and SD',
                               'forehead': {'female': [448, 30.8], 'male': [429, 37.7]},
                               'cheek': {'female': [426, 37.7]},
                               'chest': {'female': [53, 1.4], 'male': [61, 5.0]},
                               'back': {'female': [93, 6.1], 'male': [77, 5.6]},
                               'consistency_check': 'the mean of the quoted per-sex forehead values, 438.5, and of '
                                                    'the back values, 85.0, reproduce this papers own abstract '
                                                    'figures of 439 and 85, which supports the transcription',
                               'caveat': 'the cheek SD 37.7 equals the forehead-male SD, so the cheek dispersion '
                                         'may be a transcription artefact. The cheek mean is the only per-area '
                                         'number found anywhere for the beard region and it is second hand.'},
    },
    'todo2017': {
        'citation': 'Todo H. Transdermal Permeation of Drugs in Various Animal Species. Pharmaceutics '
                    '2017;9(3):33.',
        'doi': '10.3390/pharmaceutics9030033', 'pmid': '28878145', 'pmcid': 'PMC5620574', 'species': 'review',
        'method': 'review; quotes otberg2004 and blume1991 numerically',
        'cohorts': {}, 'ancestry': 'not applicable',
        'verification': 'fulltext via PMC5620574 (subagent read, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils', 'status': 'resolved',
                             'returned_title': 'Transdermal Permeation of Drugs in Various Animal Species'},
        'used_for': 'a secondary route to otberg2004 forehead 292 and back 29 per cm2 and to the blume1991 '
                    'per-sex site means. Everything taken through it is labelled second hand.',
    },
    'ali2023': {
        'citation': 'Ali RA, Fayek M, Noureldin M, El-Essawy NM. Eyebrow Restoration in Deep Facial Burn: '
                    'Follicular Unit Extraction Hair Transplantation after Nanofat Graft. Plast Reconstr Surg Glob '
                    'Open 2023;11(10):e5331.',
        'doi': '10.1097/GOX.0000000000005331', 'pmid': '37829100', 'pmcid': 'PMC10567046', 'species': 'human',
        'method': 'phototrichoscopy of the unburned contralateral eyebrow as the control side',
        'cohorts': {'eyebrow': 'n=17 patients, 15 women and 2 men, age 18-26 y, mean 21.7 +/- 2.7 y, Cairo'},
        'ancestry': 'not stated',
        'verification': 'fulltext Results read 2026-09-07 via Europe PMC PMC10567046',
        'identifier_check': {'resolver': 'eutils', 'status': 'resolved',
                             'returned_title': 'Eyebrow Restoration in Deep Facial Burn: Follicular Unit Extraction '
                                               'Hair Transplantation after Nanofat Graft'},
        'values_read_here': {'control_side_density_per_cm2': [133.95, 38.38],
                             'control_side_shaft_thickness_mm': [.06, .01],
                             'burned_recipient_side_density_per_cm2': [88.60, 29.96]},
        'caveat': 'a burn-reconstruction series, not a normative study. The control side is an unburned eyebrow of '
                  'a burn patient, mostly young women. It is the only per-cm2 eyebrow density found in a '
                  'retrievable source.',
    },
    'aumond2018': {
        'citation': 'Aumond S, Bitton E. The eyelash follicle features and anomalies: A review. J Optom '
                    '2018;11(4):211-222.',
        'doi': '10.1016/j.optom.2018.05.003', 'pmid': '30017866', 'pmcid': 'PMC6147748', 'species': 'review',
        'method': 'narrative review; the lash counts are a review statement citing its own refs 3, 6 and 7',
        'cohorts': {}, 'ancestry': 'not applicable',
        'verification': 'fulltext read 2026-09-07 via Europe PMC PMC6147748, verbatim: "The human lower lid '
                        'contains 75-80 lashes dispersed in three to four rows, whereas the upper lid has 90-160 '
                        'lashes scattered on five to six rows."',
        'identifier_check': {'resolver': 'eutils', 'status': 'resolved',
                             'returned_title': 'The eyelash follicle features and anomalies: A review'},
        'values_read_here': {'upper_lid_lashes': [90, 160], 'lower_lid_lashes': [75, 80]},
        'caveat': 'review level. The upstream primary sources were not read, so the tier is transferred and never '
                  'measured.',
    },
    'thibaut2010': {
        'citation': 'Thibaut S, De Becker E, Caisey L, Baras D, Karatas S, Jammayrac O, Pisella PJ, Bernard BA. '
                    'Human eyelash characterization. Br J Dermatol 2010;162(2):304-310.',
        'doi': '10.1111/j.1365-2133.2009.09487.x', 'pmid': '19804590', 'species': 'human',
        'method': 'high-resolution camera and image analysis, 4 volunteers followed weekly for 9 months for the '
                  'cycle; immunohistology on 17 ectropion-repair eyelid biopsies',
        'cohorts': {'eyelash': 'n=29 caucasian female volunteers aged 26-60 y'},
        'ancestry': 'caucasian, in the authors own term',
        'verification': 'abstract (eutils efetch 2026-09-07), first hand',
        'identifier_check': {'resolver': 'eutils', 'status': 'resolved',
                             'returned_title': 'Human eyelash characterization'},
        'values_read_here': {'anagen_days': [34, 9], 'full_cycle_days': [90, 5],
                             'growth_rate_mm_per_day': [.12, .05]},
        'note': 'no lash count and no lash length are reported. Length here is derived as anagen duration times '
                'growth rate, 34 d x 0.12 mm/d = 4.08 mm, and that derivation is marked as such.',
    },
    'vogt2007': {
        'citation': 'Vogt A, Hadam S, Heiderhoff M, Audring H, Lademann J, Sterry W, Blume-Peytavi U. Morphometry '
                    'of human terminal and vellus hair follicles. Exp Dermatol 2007;16(11):946-950.',
        'doi': '10.1111/j.1600-0625.2007.00602.x', 'pmid': '17927578', 'species': 'human',
        'method': 'histological sections, scalp and retroauricular region',
        'cohorts': {}, 'ancestry': 'not stated in the abstract',
        'verification': 'abstract (eutils efetch, subagent read 2026-09-07)',
        'identifier_check': {'resolver': 'eutils', 'status': 'resolved',
                             'returned_title': 'Morphometry of human terminal and vellus hair follicles',
                             'correction': 'PMID 17927577, recalled from memory, is not this paper; 17927578 is.'},
        'values_read_here': {'terminal_follicle_total_length_um': [3864, 605],
                             'terminal_infundibulum_length_um': [580, 84],
                             'vellus_follicle_total_length_um': [646, 140],
                             'vellus_infundibulum_length_um': [225, 34]},
        'used_for': 'the terminal against vellus follicle size separation. It reports no follicle density and no '
                    'shaft diameter, so no density or mass here rests on it.',
    },
    'muellner2020': {
        'citation': 'Muellner ARM, Pahl R, Brandhuber D, Peterlik H. Porosity at Different Structural Levels in '
                    'Human and Yak Belly Hair and Its Effect on Hair Dyeing. Molecules 2020;25(9):2143.',
        'doi': '10.3390/molecules25092143', 'pmid': '32375277', 'species': 'human and yak; the human value is used',
        'method': 'mass over volume regression on 10 human hairs, two groups of five from different persons',
        'cohorts': {'fibre': '10 human hairs from at least two persons'}, 'ancestry': 'not reported',
        'fibre_density_kg_m3': 1312.0, 'fibre_density_sd_kg_m3': 43.0,
        'verification': 'carried in-repo at data/measurements/hair/strand_reference.json density citing sections '
                        '2.5 and 4.6; identifiers resolved 2026-09-07',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Porosity at Different Structural Levels in Human and Yak Belly '
                                               'Hair and Its Effect on Hair Dyeing',
                             'correction': 'PMID 32370185, first written here for this paper, is Jarchi D et al, '
                                           'Recognition of Patient Groups with Sleep Related Disorders using '
                                           'Bio-signal Processing and Deep Learning, Sensors 2020;20(9):2594. '
                                           '32375277 is correct.'},
    },
    'szabo1967': {
        'citation': 'Szabo G. The regional anatomy of the human integument with special reference to the '
                    'distribution of hair follicles, sweat glands and melanocytes. Philos Trans R Soc Lond B Biol '
                    'Sci 1967;252(779):447-485.',
        'doi': '10.1098/rstb.1967.0029', 'pmid': None, 'species': 'human',
        'method': 'histological counts on cadaver integument', 'cohorts': {}, 'ancestry': 'not obtained',
        'verification': 'identifier only. royalsocietypublishing.org returns HTTP 403 and the paper has no PubMed '
                        'record, so pmid is null rather than omitted. No number is carried from it.',
        'identifier_check': {'resolver': 'crossref', 'status': 'resolved',
                             'returned_title': 'The regional anatomy of the human integument with special '
                                               'reference to the distribution of hair follicles, sweat glands and '
                                               'melanocytes',
                             'pubmed': 'empty id list for this title; no PubMed record exists'},
        'used_for': 'named as the classical regional source a later lane with library access should read. It '
                    'supplies nothing here.',
    },
}

SOURCES.update({
    'loussouarn2016': {
        'citation': 'Loussouarn G, Lozano I, Panhard S, Collaudin C, El Rawadi C, Genain G. Diversity in human '
                    'hair growth, diameter, colour and shape. An in vivo study on young adults from 24 different '
                    'ethnic groups observed in the five continents. Eur J Dermatol 2016;26(2):144-154.',
        'doi': '10.1684/ejd.2015.2726', 'pmid': '27019510', 'species': 'human',
        'method': 'phototrichogram at vertex, nape and temple; shaft diameter by LaserScan on 823 subjects',
        'cohorts': {'scalp': 'n=2249, 1065 men and 1184 women, age 18-35 y (mean 26 +/- 5), 24 groups on five '
                             'continents'},
        'ancestry': 'the authors 24 group labels, and their own caution, verbatim: "terms such as Danish or Thai '
                    'embraced in the present paper should be solely viewed as arbitrary shortcuts"; their Figure 7 '
                    'shows "a complete continuum of data points, where population clusters are not readily '
                    'apparent"',
        'verification': 'fulltext, Tables 1, 2 and 5 (subagent read the publisher open PDF, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Diversity in human hair growth, diameter, colour and shape. An in '
                                               'vivo study on young adults from 24 different ethnic groups '
                                               'observed in the five continents'},
        'values_read_here': {'total_density_range_per_cm2': [153, 233],
                             'clusters': {'african_type': {'density_cm2': 178, 'telogen_percent': 12,
                                                           'growth_um_per_day': 325, 'diameter_um': 71},
                                          'asian_type': {'density_cm2': 182, 'telogen_percent': 12,
                                                         'growth_um_per_day': 413, 'diameter_um': 84},
                                          'caucasian_type': {'density_cm2': 214, 'telogen_percent': 10,
                                                             'growth_um_per_day': 386, 'diameter_um': 73}},
                             'telogen_percent_range': [8, 14], 'telogen_men_women_percent': [12.2, 10.1],
                             'diameter_median_range_um': [69, 89],
                             'sex': 'men lower only at the vertex, by about 19 hairs per cm2',
                             'growth_range_um_per_day': [272, 426]},
        'used_for': 'the population spread of scalp density, diameter, growth rate and telogen fraction. The '
                    'candidate carries a different cohort for the scalp value itself; this source bounds how far '
                    'that value can be from another population.',
    },
    'freitag2014': {
        'citation': 'Freitag FM, Cahill K, Wu A, Nakra T, Wojno T, Woodward JA, Rubin PAD, Yoon MK. '
                    'Retrospective review of eyelash number in patients who have undergone full-thickness eyelid '
                    'resection. Ophthalmic Plast Reconstr Surg 2014;30(1):1-6.',
        'doi': '10.1097/IOP.0b013e3182a650bb', 'pmid': '24398479', 'species': 'human',
        'method': 'eyelash counting from postoperative photographs; the contralateral unoperated lid is the control',
        'cohorts': {'eyelash': 'n=38 patients, 10 men and 28 women, mean age 57.9 y, range 14-86 y'},
        'ancestry': 'not stated',
        'verification': 'abstract (subagent read via eutils efetch, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Retrospective review of eyelash number in patients who have '
                                               'undergone full-thickness eyelid resection'},
        'values_read_here': {'upper_lid_lashes': 72.1, 'lower_lid_lashes': 38.2,
                             'sex': 'no statistical significance between men and women'},
        'caveat': 'eyelid-lesion patients and an older cohort than the rest of this table. It is nevertheless the '
                  'only MEASURED lash count found; the 75-160 per lid in aumond2018 is a review statement whose '
                  'primaries were not read, and it disagrees with this measurement by roughly a factor of two.',
    },
    'elder1997': {
        'citation': 'Elder MJ. Anatomy and physiology of eyelash follicles: relevance to lash ablation '
                    'procedures. Ophthalmic Plast Reconstr Surg 1997;13(1):21-25.',
        'doi': '10.1097/00002341-199703000-00004', 'pmid': '9076779', 'species': 'human',
        'method': 'histological sections of upper and lower eyelids',
        'cohorts': {'eyelash': 'n=10 patients; the hair-cycle classification rests on a single patient'},
        'ancestry': 'not stated',
        'verification': 'abstract (subagent read via eutils efetch, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Anatomy and physiology of eyelash follicles: relevance to lash '
                                               'ablation procedures'},
        'values_read_here': {'upper_lid_follicle_depth_mm': [1.8, .3], 'upper_lid_bulb_width_um': [188, 44],
                             'upper_lid_shaft_width_um': [205, 28], 'lower_lid_follicle_depth_mm': [.9, .2],
                             'lower_lid_bulb_width_um': [132, 19], 'lower_lid_shaft_width_um': [158, 26],
                             'active_follicle_fraction': {'upper': .41, 'lower': .15, 'n': 1}},
        'used_for': 'follicle depth only. The 205 and 158 um shaft widths are NOT used as the fibre radius: they '
                    'are two to three times any reported hair-fibre diameter, so they are read as the follicular '
                    'shaft in section rather than the emergent fibre, and no mass rests on them.',
    },
    'floyd2018': {
        'citation': 'Floyd EL, Henry JB, Johnson DL. Influence of facial hair length, coarseness, and areal '
                    'density on seal leakage of a tight-fitting half-face respirator. J Occup Environ Hyg '
                    '2018;15(4):334-340.',
        'doi': '10.1080/15459624.2017.1416388', 'pmid': '29283316', 'species': 'human',
        'method': 'direct measurement of cheek and chin hair diameter and areal density in bearded men',
        'cohorts': {'beard': 'n=19 bearded male subjects'}, 'ancestry': 'not stated',
        'verification': 'abstract (subagent read via eutils efetch, 2026-09-07); the areal densities are in the '
                        'paywalled full text and were not obtained',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Influence of facial hair length, coarseness, and areal density on '
                                               'seal leakage of a tight-fitting half-face respirator'},
        'values_read_here': {'beard_shaft_diameter_um': [76, 7.4],
                             'note': 'normally distributed across subjects; chin and cheek areal densities differ '
                                     'significantly and correlate only weakly, but their values are not in the '
                                     'abstract'},
    },
    'seago1985': {
        'citation': 'Seago SV, Ebling FJ. The hair cycle on the human thigh and upper arm. Br J Dermatol '
                    '1985;113(1):9-16.',
        'doi': '10.1111/j.1365-2133.1985.tb02038.x', 'pmid': '4015973', 'species': 'human',
        'method': 'repeated observation of defined thigh and upper-arm areas',
        'cohorts': {'body_hair': 'n=20, 11 female and 9 male, aged 20-30 y'}, 'ancestry': 'not stated',
        'verification': 'abstract (subagent read via eutils efetch, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'The hair cycle on the human thigh and upper arm'},
        'values_read_here': {'thigh_anagen_days': {'male': 54, 'female': 22},
                             'thigh_cycle_days': {'male': 151, 'female': 84},
                             'upper_arm_anagen_days': {'male': 28, 'female': 22},
                             'upper_arm_cycle_days': {'male': 108, 'female': 106},
                             'sex_statement': 'follicle density did not differ between males and females; on the '
                                              'thigh the definitive length of hair was on average three times '
                                              'greater in males, mainly a longer anagen (x2.46) and partly a '
                                              'greater growth rate (x1.22)'},
        'used_for': 'the anagen fractions on thigh and upper arm, and the measured statement that the sex '
                    'difference in body hair is in anagen duration and shaft length, not in follicle number. '
                    'Absolute definitive lengths in mm are not in the abstract.',
    },
    'bouabbache2019': {
        'citation': 'Bouabbache S, Panhard S, Loussouarn G. Exploring some characteristics (density, anagen '
                    'ratio, growth rate) of human body hairs. Variations with skin sites, gender and ethnics. Int '
                    'J Cosmet Sci 2019;41(1):46-54.',
        'doi': '10.1111/ics.12510', 'pmid': '30580453', 'species': 'human',
        'method': 'photo-trichogram at axilla, cheek, chin, leg and upper lip, with a nape reference',
        'cohorts': {'body_hair': 'women and men in four groups the authors call African, Caucasian, Chinese and '
                                 'North African'},
        'ancestry': 'the authors terms: African, Caucasian, Chinese, North African',
        'verification': 'abstract (subagent read via eutils efetch, 2026-09-07); the per-site tables are '
                        'paywalled and were not obtained',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Exploring some characteristics (density, anagen ratio, growth '
                                               'rate) of human body hairs. Variations with skin sites, gender and '
                                               'ethnics'},
        'values_read_here': {'body_growth_rate_um_per_day': [180, 485],
                             'men_face_anagen_percent': 'greater than 85 on cheek and upper lip in Caucasians and '
                                                        'North Africans, lower in African and Chinese men',
                             'women': 'anagen phase percentage higher on leg and armpit; facial hairs in women '
                                      'were too thin for the photo-trichogram to resolve'},
        'used_for': 'the axilla, cheek and chin are measured here, but no per-cm2 density is in the abstract. '
                    'This is the nearest source to the axillary field and it supplies no number to it.',
    },
    'visessiri2020': {
        'citation': 'Visessiri Y, Pakornphadungsit K, Leerunyakul K, Rutnin S, Vachiramon V, Suchonwanit P. The '
                    'study of hair follicle counts from scalp histopathology in the Thai population. Int J '
                    'Dermatol 2020;59(8):978-981.',
        'doi': '10.1111/ijd.14989', 'pmid': '32501534', 'species': 'human',
        'method': 'occipital 4 mm punch biopsy, horizontal sections',
        'cohorts': {'scalp': 'n=90'}, 'ancestry': 'Thai population, in the authors own term',
        'verification': 'abstract (subagent read via eutils efetch, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'The study of hair follicle counts from scalp histopathology in the '
                                               'Thai population'},
        'values_read_here': {'total_hairs_per_4mm_punch': [20.5, 5.2], 'terminal_per_punch': [18.2, 4.1],
                             'vellus_median_per_punch': 2, 'follicular_units_per_punch': [9.1, 1.6],
                             'terminal_to_vellus': 8.9, 'anagen_percent': 91.9, 'telogen_percent': 7.9},
        'derived_here': 'a 4 mm punch is pi*(0.2 cm)^2 = 0.12566 cm2, so 20.5 hairs is 163.1 hairs per cm2, which '
                        'agrees with the 154-163 per cm2 that leerunyakul2020 measured by trichoscopy in the same '
                        'population.',
    },
    'patel2026': {
        'citation': 'Patel and Patel. Standardization and Validation of Scalp Surface Area Measurement and '
                    'Assessment of Inter-relationship Between Head Crown Area and Total Hair Count Measurement. '
                    'Cureus 2026;18(2):e104173.',
        'doi': '10.7759/cureus.104173', 'pmid': '41909352', 'species': 'human',
        'method': 'scalp surface area in n=25; hair count from n=50 archival global 90 degree photographs by '
                  'Image-Pro and by CASLite NOVA',
        'cohorts': {'scalp_area': 'n=25 adults aged 18-55, Ahmedabad', 'hair_count': 'n=50 archival photographs'},
        'ancestry': 'not stated',
        'verification': 'fulltext (subagent read, open access, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Standardization and Validation of Scalp Surface Area Measurement '
                                               'and Assessment of Inter-relationship Between Head Crown Area and '
                                               'Total Hair Count Measurement'},
        'values_read_here': {'scalp_surface_area_cm2': [492.36, 30.97], 'scalp_area_range_cm2': [447.29, 592.52],
                             'total_hair_count_image_pro': [63996, 12534.76],
                             'total_hair_count_caslite': [66426, 14643.13]},
        'used_for': 'an independent comparison for the measured scalp footprint area of this candidate. Their '
                    'counted region is the head crown including frontal and mid-scalp, not the whole hair-bearing '
                    'scalp, so it is a comparison and not a target.',
    },
    'barman1965': {
        'citation': 'Barman JM, Astore I, Pecoraro V. The normal trichogram of the adult. J Invest Dermatol '
                    '1965;44:233-236.',
        'doi': '10.1038/jid.1965.42', 'pmid': '14282465', 'species': 'human', 'method': 'trichogram',
        'cohorts': {}, 'ancestry': 'not obtained',
        'verification': 'identifier only, confirmed through PubMed, Europe PMC, Unpaywall and fatcat. PubMed and '
                        'Europe PMC carry no abstract, no open copy exists, and jidonline.org returns HTTP 403. '
                        'No number is carried from it.',
        'identifier_check': {'resolver': 'eutils+crossref+europepmc+fatcat', 'status': 'resolved',
                             'returned_title': 'THE NORMAL TRICHOGRAM OF THE ADULT'},
        'used_for': 'named as the classical adult trichogram a later lane should read. It supplies nothing here.',
    },
    'saitoh1970': {
        'citation': 'Saitoh M, Uzuka M, Sakamoto M. Human hair cycle. J Invest Dermatol 1970;54(1):65-81.',
        'doi': '10.1111/1523-1747.ep12551679', 'pmid': '5416680', 'species': 'human', 'method': 'hair cycle study',
        'cohorts': {}, 'ancestry': 'not obtained',
        'verification': 'identifier only; same four-route check and same result as barman1965. Any anagen or '
                        'telogen by-body-site table attributed to this paper is unverified here.',
        'identifier_check': {'resolver': 'eutils+crossref+europepmc+fatcat', 'status': 'resolved',
                             'returned_title': 'Human hair cycle'},
        'used_for': 'nothing. Recorded so that a later lane does not re-derive numbers from a citation this lane '
                    'could not read.',
    },
})

SOURCES.update({
    'coelho2026': {
        'citation': 'Efficacy and safety of intense pulsed light compared to diode Laser for hair removal: a '
                    'randomized controlled trial. Lasers Med Sci 2026;41(1):108.',
        'doi': '10.1007/s10103-026-04904-6', 'pmid': '42249955', 'pmcid': 'PMC13242373', 'species': 'human',
        'method': 'baseline of a randomised hair-removal trial; 75 cm2 template with a central 4 cm2 window, '
                  'magnified photographs, blinded manual count in ImageJ counting every emerging hair including '
                  'vellus',
        'cohorts': {'axilla': 'n=48, all female, Sao Paulo, median age 25.5 y (IQR 22-38.7), Fitzpatrick I-IV'},
        'ancestry': 'not stated beyond Fitzpatrick phototype',
        'verification': 'abstract and open text (subagent read, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Efficacy and safety of intense pulsed light compared to diode '
                                               'Laser for hair removal: a randomized controlled trial'},
        'values_read_here': {'axillary_visible_hairs_per_cm2': [[11.1, .50], [11.0, .50]],
                             'axillary_shaft_thickness_mm': [[.183, .005], [.178, .005]]},
        'caveat': 'a visible-hair count at trial baseline in women, not a follicle density, and the arms were '
                  'shaved into the protocol. It is the only MEASURED axillary number this sweep retrieved.',
    },
    'pecoraro1971': {
        'citation': 'Growth rate and hair density of the human axilla. A comparative study of normal males and '
                    'females and pregnant and post-partum females. J Invest Dermatol 1971;56(5):362-365.',
        'doi': '10.1111/1523-1747.ep12261236', 'pmid': '5556525', 'species': 'human', 'method': 'axillary trichogram',
        'cohorts': {}, 'ancestry': 'not obtained',
        'verification': 'identifier only. PubMed carries no abstract, Unpaywall reports not open access, and both '
                        'jidonline.org and sciencedirect return HTTP 403. No number is carried.',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Growth rate and hair density of the human axilla. A. Comparative '
                                               'study of normal males and females and pregnant and post-partum '
                                               'females'},
        'used_for': 'this is the single most on-target paper for the axillary field and it could not be read. '
                    'Named so a later lane with library access goes straight to it.',
    },
    'astore1979': {
        'citation': 'Astore IP, Pecoraro V, Pecoraro EG. The normal trichogram of pubic hair. Br J Dermatol '
                    '1979;101(4):441-445.',
        'doi': '10.1111/j.1365-2133.1979.tb00023.x', 'pmid': '508610', 'species': 'human', 'method': 'trichogram',
        'cohorts': {'pubic': 'males and non-pregnant, pregnant and post-partum females'}, 'ancestry': 'not stated',
        'verification': 'abstract, which carries direction only: density and growth rate decrease with age in '
                        'males and non-pregnant females, telogen percentage higher in women, hair thickness not '
                        'modified. No numbers. Full text closed.',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'The normal trichogram of pubic hair'},
        'used_for': 'the pubic field carries no density because this, its primary source, publishes none that '
                    'could be retrieved.',
    },
    'hwang1997': {
        'citation': 'Hwang K, Baik SH. Distribution of hairs and sweat glands on the bodies of Korean adults: a '
                    'morphometric study. Acta Anat (Basel) 1997;158(2):112-120.',
        'doi': '10.1159/000147920', 'pmid': '9311420', 'species': 'human',
        'method': 'punch biopsy at 30 body regions in 74 adult cadavers, hairs counted per cm2 in serial '
                  'transverse sections, means and SD per region, sex differences tested',
        'cohorts': {'regional': 'n=74 Korean adult cadavers, 30 body regions'},
        'ancestry': 'Korean adults, in the authors own term',
        'verification': 'abstract only; Karger full text closed. Abstract-level findings: hairs densest in the '
                        'buccal region in men and the mental region in women, sparsest on the dorsum of the '
                        'middle phalanx of the fingers in both, and only four of thirty regions showed '
                        'significant sexual dimorphism in hair density.',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Distribution of hairs and sweat glands on the bodies of Korean '
                                               'adults: a morphometric study'},
        'used_for': 'nothing numerically. It is the one retrievable-in-principle source that covers 30 regions '
                    'including the ones this candidate leaves absent, so it is the highest-value target for the '
                    'next lane.',
    },
    'xu2017': {
        'citation': 'Xu H, Yu S, Lin C, Dong D, Xiao J, Ye Y, Wang M. Reference values for skin microanatomy: A '
                    'systematic review and meta-analysis of ex vivo studies. J Am Acad Dermatol '
                    '2017;77(6):1133-1144.e4.',
        'doi': '10.1016/j.jaad.2017.06.009', 'pmid': '28716435', 'pmcid': 'PMC5685878', 'species': 'human',
        'method': 'systematic review and meta-analysis pooling 56 ex vivo studies',
        'cohorts': {'pooled': '56 articles, all anatomic locations combined'}, 'ancestry': 'pooled',
        'verification': 'abstract (subagent read via eutils efetch, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Reference values for skin microanatomy: A systematic review and '
                                               'meta-analysis of ex vivo studies'},
        'values_read_here': {'pooled_follicle_density_per_mm2': [1.40, [.91, 1.89]],
                             'as_per_cm2': [140, [91, 189]]},
        'used_for': 'a whole-body pooled comparator for the area-weighted mean density this candidate produces. '
                    'The authors own limitation is significant heterogeneity from differing histological '
                    'technique and absent standardised definitions.',
    },
    'berg2025': {
        'citation': 'Van den Berg C, Khumalo NP, Ngoepe MN. Quantifying whole human hair scalp fibres of varying '
                    'curl: A micro-computed tomographic study. J Microsc 2025;297(2):227-251.',
        'doi': '10.1111/jmi.13365', 'pmid': '39564786', 'pmcid': 'PMC11733847', 'species': 'human',
        'method': 'laser micrometry (FDAS) on 15-18 strands per type and micro-CT at 2 um on 3 strands per type',
        'cohorts': {'fibre': 'Caucasian types II and IV and African type VI, age 20-30 y; sex not stated'},
        'ancestry': 'Caucasian and African, in the authors own terms',
        'verification': 'fulltext tables (subagent read, open access, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Quantifying whole human hair scalp fibres of varying curl: A '
                                               'micro-computed tomographic study'},
        'values_read_here': {'type_II_low_curl': {'major_um': 68.2, 'minor_um': 53.6, 'ellipticity': 1.28,
                                                  'cross_section_um2': 2870},
                             'type_IV_medium_curl': {'major_um': 79.3, 'minor_um': 49.2, 'ellipticity': 1.6,
                                                     'cross_section_um2': 3070},
                             'type_VI_high_curl': {'major_um': 100.7, 'minor_um': 58.2, 'ellipticity': 1.7,
                                                   'cross_section_um2': 4680},
                             'micro_ct_cross_section_range_um2': [2800, 7420],
                             'micro_ct_ellipticity_range': [1.2, 1.9]},
        'used_for': 'the measured bracket on the circular-section idealisation this candidate uses for mass. A '
                    'hair is elliptical with an aspect ratio of 1.28 to 1.9, and its measured cross-section is '
                    '2870 to 4680 um2.',
    },
    'wikramanayake2012': {
        'citation': 'Wikramanayake TC, Mauro LM, Tabas IA, Chen AL, Llanes IC, Jimenez JJ. Cross-section '
                    'Trichometry: A Clinical Tool for Assessing the Progression and Treatment Response of '
                    'Alopecia. Int J Trichology 2012;4(4):259-264.',
        'doi': '10.4103/0974-7753.111221', 'pmid': '23766610', 'pmcid': 'PMC3681107', 'species': 'human',
        'method': 'cross-sectional area of an isolated hair bundle from a 2x2 cm scalp area',
        'cohorts': {'hair_mass_index': 'no n, sex, age or population given for the quoted normal range'},
        'ancestry': 'not stated',
        'verification': 'fulltext (subagent read, open access, 2026-09-07), verbatim: "In the absence of hair '
                        'loss, normal HMI ranged from 75 for fine hair to 100 for coarse hair with an average of '
                        '87", where HMI is mm2 of hair per cm2 of scalp times 100',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Cross-section Trichometry: A Clinical Tool for Assessing the '
                                               'Progression and Treatment Response of Alopecia'},
        'values_read_here': {'hair_cross_section_mm2_per_cm2_scalp': [.75, 1.00], 'average': .87},
        'used_for': 'an independent, count-free route to scalp hair mass: cross-section per unit scalp area times '
                    'length times fibre density. It checks the count-times-shaft route without sharing its inputs.',
        'caveat': 'the normal range is an uncharacterised clinical reference, not a described cohort.',
    },
    'gao2022': {
        'citation': 'Gao M, Wang Y, Xu H, Xu C, Yang X, Nie J, Zhang Z, Li Z, Hou W, Liu Y. Deep Learning-based '
                    'Trichoscopic Image Analysis and Quantitative Model for Predicting Basic and Specific '
                    'Classification in Male Androgenetic Alopecia. Acta Derm Venereol 2022;102:adv00635.',
        'doi': '10.2340/actadv.v101.564', 'pmid': '34935989', 'pmcid': 'PMC9631273', 'species': 'human',
        'method': 'annotation protocol over 2910 trichoscopic images',
        'cohorts': {}, 'ancestry': 'Chinese Academy of Medical Sciences, Nanjing',
        'verification': 'fulltext (subagent read, 2026-09-07), verbatim: "blue for vellus hairs < 0.03 mm '
                        'diameter; green for intermediate hairs 0.03-0.06 mm diameter, and red for terminal hairs '
                        '> 0.06 mm diameter"',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Deep Learning-based Trichoscopic Image Analysis and Quantitative '
                                               'Model for Predicting Basic and Specific Classification in Male '
                                               'Androgenetic Alopecia'},
        'values_read_here': {'vellus_diameter_um_max': 30, 'intermediate_diameter_um': [30, 60],
                             'terminal_diameter_um_min': 60},
        'caveat': 'the paper states these cutoffs as its own annotation protocol and cites no external source for '
                  'them. A length criterion exists too and the sources disagree: asfiya2025 uses under 2-3 mm for '
                  'vellus, yildiz2010 uses under 5 mm.',
    },
    'yildiz2010': {
        'citation': 'Yildiz BO, Bolour S, Woods K, Moore A, Azziz R. Visually scoring hirsutism. Hum Reprod '
                    'Update 2010;16(1):51-64.',
        'doi': '10.1093/humupd/dmp024', 'pmid': '19567450', 'pmcid': 'PMC2792145', 'species': 'human',
        'method': 'review; the terminal hair counts are computerised image analysis through a calibrated glass '
                  'plate, attributed by the authors to Hines et al. 2001',
        'cohorts': {'hirsute_women': 'n=20, 12 white and 8 Black'},
        'ancestry': 'white and Black, in the authors own terms',
        'verification': 'fulltext (subagent read, open access, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Visually scoring hirsutism'},
        'values_read_here': {'hirsute_facial_terminal_hairs_per_cm2': [15.6, 14.2],
                             'hirsute_abdominal_terminal_hairs_per_cm2': [5.4, 1.9],
                             'terminal_length_min_cm': .5, 'vellus_length_max_cm': .5},
        'used_for': 'the length-based terminal against vellus criterion, and as the only retrieved measurement of '
                    'terminal facial and abdominal hair per cm2 in any cohort. It is a hirsute female cohort, so '
                    'no field takes its value.',
    },
    'tsai2022': {
        'citation': 'Tsai J, Rostom M, Garza LA. Understanding and Harnessing Epithelial-Mesenchymal '
                    'Interactions in the Development of Palmoplantar Identity. J Invest Dermatol '
                    '2022;142(2):282-284.',
        'doi': '10.1016/j.jid.2021.06.016', 'pmid': '34366107', 'pmcid': 'PMC8792145', 'species': 'human',
        'method': 'commentary stating the defining characteristics of palmoplantar skin',
        'cohorts': {}, 'ancestry': 'not applicable',
        'verification': 'fulltext (subagent read, 2026-09-07), verbatim: "Palmoplantar skin has several unique '
                        'characteristics such as increased thickness, high resilience, hypopigmentation, and lack '
                        'of hair follicles"',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Understanding and Harnessing Epithelial-Mesenchymal Interactions '
                                               'in the Development of Palmoplantar Identity'},
        'used_for': 'the zero follicle density of the palmoplantar part of the glabrous field. Corroborated by '
                    'rittie2013 and fradette1995.',
    },
    'rittie2013': {
        'citation': 'Rittie L, Sachs DL, Orringer JS, Voorhees JJ, Fisher GJ. Eccrine sweat glands are major '
                    'contributors to reepithelialization of human wounds. Am J Pathol 2013;182(1):163-171.',
        'doi': '10.1016/j.ajpath.2012.09.019', 'pmid': '23159944', 'pmcid': 'PMC3538027', 'species': 'human',
        'method': 'histology; the statement used is definitional',
        'cohorts': {}, 'ancestry': 'not applicable',
        'verification': 'fulltext (subagent read, 2026-09-07), verbatim: "glabrous skin (devoid of hair '
                        'follicles) such as on palms and soles"',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Eccrine sweat glands are major contributors to '
                                               'reepithelialization of human wounds'},
    },
    'fradette1995': {
        'citation': 'Fradette J, Godbout MJ, Michel M, Germain L. Localization of Merkel cells at hairless and '
                    'hairy human skin sites using keratin 18. Biochem Cell Biol 1995;73(9-10):635-639.',
        'doi': '10.1139/o95-070', 'pmid': '8714683', 'species': 'human', 'method': 'immunohistochemistry',
        'cohorts': {}, 'ancestry': 'not stated',
        'verification': 'abstract (subagent read, 2026-09-07), verbatim: "their density is higher at hairless '
                        'anatomic sites such as palms and soles"',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Localization of Merkel cells at hairless and hairy human skin '
                                               'sites using keratin 18'},
    },
    'asfiya2025': {
        'citation': 'Asfiya A, et al. Trichoscopic Analysis of Healthy Indian Males for Standardization of the '
                    'Measurable Parameters: An Observational, Cross-sectional Study. Indian Dermatol Online J '
                    '2025;16(1):81-86.',
        'doi': '10.4103/idoj.idoj_201_24', 'pmid': '39850670', 'pmcid': 'PMC11753574', 'species': 'human',
        'method': 'videodermoscopy at frontal, occipital and temporal scalp',
        'cohorts': {'scalp': 'n=120 healthy males, age 18-50 y, mean 24.6 y'},
        'ancestry': 'Indian males, in the authors own term',
        'verification': 'fulltext (subagent read, open access, 2026-09-07)',
        'identifier_check': {'resolver': 'eutils+crossref', 'status': 'resolved',
                             'returned_title': 'Trichoscopic Analysis of Healthy Indian Males for Standardization '
                                               'of the Measurable Parameters: An Observational, Cross-sectional '
                                               'Study'},
        'values_read_here': {'frontal_density_cm2': [160.05, 86.89], 'occipital_density_cm2': [156.27, 97.82],
                             'temporal_density_cm2': [123.6, 64.8], 'shaft_diameter_mm': [.1, .04],
                             'anagen_percent': {'frontal': [73.84, 14.40], 'occipital': [74.65, 12.59],
                                                'temporal': [79.5, 10.94]},
                             'vellus_definition': 'nonmedullated, hypopigmented hairs which are under 0.03 mm '
                                                  'thick and under 2-3 mm long'},
        'caveat': 'its telogen fractions, 20-26 percent, are two to three times the phototrichogram values in '
                  'loussouarn2016 and the histology in visessiri2020. That is a method difference, not a '
                  'population fact, and the three must not be blended.',
    },
})

_VELLUS_R = {'value': 1.0e-5, 'band': [8.0e-6, 1.6e-5], 'tier': 'assumed', 'source': [],
             'basis': 'engineering scenario inside the conventional vellus ceiling of a 30 um shaft (gao2022); the '
                      'band is the radius range the repo already materialises in '
                      'data/derived/hair/elastic_v3. No measured regional vellus shaft diameter was retrievable: '
                      'the otberg2004 and mangelsdorf2006 tables that hold it are behind HTTP 403.'}
_VELLUS_L = {'value': 1.0e-3, 'band': [5.0e-4, 2.0e-3], 'tier': 'assumed', 'source': [],
             'basis': 'the 0.5 to 2 mm scenario the repo already materialises. yildiz2010 puts the vellus length '
                      'ceiling at 5 mm and asfiya2025 at 2 to 3 mm; the sources disagree and none measures a mean.'}
_FOLLICLE_FRACTION = {'value': 1.0, 'tier': 'assumed', 'source': ['jonsson2017'],
                      'basis': 'cyanoacrylate stripping counts follicular orifices, not shafts. Assuming every '
                               'counted follicle carries one shaft over-states mass by at most the telogen '
                               'fraction, which loussouarn2016 puts at 8 to 14 percent on the scalp and which no '
                               'source measures for these body sites.'}
_HAIR_FRACTION = {'value': 1.0, 'tier': 'measured', 'source': ['leerunyakul2020'],
                  'basis': 'trichoscopy counts visible hairs, so the density is already a shaft count and no '
                           'shaft-bearing correction applies.'}
_NO_LENGTH = {'value': None, 'band': None, 'tier': None, 'source': [],
              'basis': 'no measured length was retrieved for this field, so its mass is absent rather than '
                       'invented.'}


def _vellus(**extra):
    record = {'shaft_radius_m': dict(_VELLUS_R), 'shaft_length_m': dict(_VELLUS_L),
              'shaft_bearing_fraction': dict(_FOLLICLE_FRACTION), 'growth': {}, 'colour': None}
    record.update(extra)
    return record


MORPHOLOGY = {
    '_fibre_density_source': 'muellner2020',
    '_cross_section_model': 'circular section of the stated radius. berg2025 measures the real section as an '
                            'ellipse of aspect ratio 1.28 to 1.9 with an area of 2870 to 4680 um2; a circle of '
                            'the 80.6 um mean scalp diameter has 5102 um2, so this idealisation over-states scalp '
                            'fibre section by roughly 1.1 to 1.8 times.',
    'scalp': {'shaft_radius_m': {'value': 4.028750e-5, 'band': [4.005e-5, 4.055e-5], 'tier': 'measured',
                                 'source': ['leerunyakul2020'],
                                 'basis': 'half the mean of the four measured scalp-area diameters 81.1, 80.8, '
                                          '80.1 and 80.3 um; the band is half the min and max of those four'},
              'shaft_length_m': {'value': .03, 'band': [.01, .60], 'tier': 'assumed', 'source': [],
                                 'basis': 'the 30 mm haircut scenario the repo already materialises. No source in '
                                          'this sweep measures a mean scalp hair length, and the band spans a '
                                          'close crop to long hair, so scalp mass scales linearly with a choice '
                                          'that is not evidence.'},
              'shaft_bearing_fraction': dict(_HAIR_FRACTION),
              'growth': {'rate_m_per_day': {'value': None, 'band': [2.72e-4, 4.26e-4], 'tier': 'measured',
                                            'source': ['loussouarn2016'],
                                            'basis': 'group means from 272 to 426 um/day over 2249 subjects; the '
                                                     'single value is left absent because the specimen has no '
                                                     'declared population'},
                         'telogen_fraction': {'value': None, 'band': [.08, .14], 'tier': 'measured',
                                              'source': ['loussouarn2016', 'visessiri2020'],
                                              'basis': 'phototrichogram 8 to 14 percent (loussouarn2016, n=2249) '
                                                       'and histology 7.9 percent (visessiri2020). asfiya2025 '
                                                       'reports 20 to 26 percent by trichoscopy software; that is '
                                                       'a method difference and is not blended in.'},
                         'anagen_fraction': {'value': None, 'band': [.86, .92], 'tier': 'derived',
                                             'source': ['loussouarn2016'],
                                             'basis': '1 minus the telogen band'}},
              'colour': {'value': None, 'tier': None, 'source': ['loussouarn2016'],
                         'basis': 'no colour is assigned. loussouarn2016 measures that lighter tones are thinner, '
                                  '72 +/- 9 um against 81 +/- 10 um, so colour and diameter are not independent, '
                                  'but no source gives a colour for an unspecified reference subject.'}},
    'eyebrow': {'shaft_radius_m': {'value': 3.0e-5, 'band': [2.5e-5, 3.5e-5], 'tier': 'measured',
                                   'source': ['ali2023'],
                                   'basis': 'half the 0.06 +/- 0.01 mm control-side eyebrow shaft thickness'},
                'shaft_length_m': {'value': .01, 'band': [.005, .02], 'tier': 'assumed', 'source': [],
                                   'basis': 'no measured eyebrow length was retrieved; flament2019 measured '
                                            'eyebrow length in 3600 women but publishes no number in the abstract '
                                            'and is paywalled'},
                'shaft_bearing_fraction': {'value': 1.0, 'tier': 'measured', 'source': ['ali2023'],
                                           'basis': 'phototrichoscopy counts visible hairs'},
                'growth': {}, 'colour': None},
    'eyelash': {'shaft_radius_m': {'value': 4.028750e-5, 'band': [3.0e-5, 5.0e-5], 'tier': 'transferred',
                                   'source': ['leerunyakul2020', 'thibaut2010'],
                                   'basis': 'transferred from the scalp radius on thibaut2010, which found '
                                            'eyelash follicle morphology very close to the scalp hair follicle. '
                                            'The 205 and 158 um shaft widths in elder1997 are rejected as the '
                                            'fibre radius: they exceed every reported hair-fibre diameter and '
                                            'read as the follicular shaft in section.'},
                'shaft_length_m': {'value': 4.08e-3, 'band': [1.75e-3, 7.31e-3], 'tier': 'derived',
                                   'source': ['thibaut2010'],
                                   'basis': 'anagen duration times growth rate, 34 d x 0.12 mm/d = 4.08 mm; the '
                                            'band is (34-9) x (0.12-0.05) to (34+9) x (0.12+0.05). This is well '
                                            'below the 8 to 12 mm usually quoted for lashes, so either the '
                                            'calculated anagen is not the full growth window or the two '
                                            'parameters come from different subsets. Treat as unreliable.'},
                'shaft_bearing_fraction': {'value': 1.0, 'tier': 'measured', 'source': ['freitag2014'],
                                           'basis': 'lashes were counted from photographs, so the count is '
                                                    'already of visible shafts'},
                'growth': {'rate_m_per_day': {'value': 1.2e-4, 'band': [7e-5, 1.7e-4], 'tier': 'measured',
                                              'source': ['thibaut2010'], 'basis': '0.12 +/- 0.05 mm/day'},
                           'cycle_days': {'value': 90, 'band': [85, 95], 'tier': 'measured',
                                          'source': ['thibaut2010'], 'basis': '90 +/- 5 days'},
                           'anagen_fraction': {'value': .378, 'band': None, 'tier': 'derived',
                                               'source': ['thibaut2010'], 'basis': '34 d anagen over a 90 d cycle'}},
                'colour': None},
    'beard': {'shaft_radius_m': {'value': 3.8e-5, 'band': [3.43e-5, 4.17e-5], 'tier': 'measured',
                                 'source': ['floyd2018'],
                                 'basis': 'half the 76 +/- 7.4 um cheek and chin hair diameter measured in 19 '
                                          'bearded men'},
              'shaft_length_m': dict(_NO_LENGTH),
              'shaft_bearing_fraction': {'value': None, 'tier': None, 'source': [],
                                         'basis': 'the fraction of cheek follicles that produce a terminal beard '
                                                  'shaft in an adult male is not measured in any source this '
                                                  'sweep retrieved. Without it the beard mass is absent rather '
                                                  'than invented, even though its shaft diameter is measured.'},
              'growth': {'anagen_fraction': {'value': None, 'band': [.85, 1.0], 'tier': 'measured',
                                             'source': ['bouabbache2019'],
                                             'basis': 'greater than 85 percent on cheek and upper lip in '
                                                      'Caucasian and North African men, lower in African and '
                                                      'Chinese men, in the authors own terms'}},
              'colour': None},
    'axillary': _vellus(shaft_radius_m={'value': 9.15e-5, 'band': [8.9e-5, 9.4e-5], 'tier': 'measured',
                                        'source': ['coelho2026'],
                                        'basis': 'half the 0.183 +/- 0.005 mm axillary shaft thickness measured '
                                                 'at trial baseline in 48 women'},
                        growth={'rate_m_per_day': {'value': None, 'band': [1.8e-4, 4.85e-4], 'tier': 'measured',
                                                   'source': ['bouabbache2019'],
                                                   'basis': 'body hair growth 180 to 485 um/day over five sites; '
                                                            'on the axilla close to terminal hair'}}),
    'facial_vellus': _vellus(growth={'rate_m_per_day': {'value': 3.0e-5, 'band': [3.0e-5, 1.3e-4],
                                                        'tier': 'measured', 'source': ['blume1991'],
                                                        'basis': '0.03 mm/day on the forehead, the low end of '
                                                                 'this papers 0.03 to 0.13 mm/day range'},
                                     'anagen_fraction': {'value': .49, 'band': None, 'tier': 'measured',
                                                         'source': ['blume1991'],
                                                         'basis': '49 percent growing hairs on the forehead'}}),
    'back': _vellus(growth={'rate_m_per_day': {'value': 1.3e-4, 'band': None, 'tier': 'measured',
                                               'source': ['blume1991'], 'basis': '0.13 mm/day on the back'},
                            'anagen_fraction': {'value': .315, 'band': None, 'tier': 'measured',
                                                'source': ['blume1991'],
                                                'basis': '31.5 percent growing hairs on the back'}}),
    'thigh': _vellus(growth={'anagen_fraction': {'value': .358, 'band': [.262, .358], 'tier': 'derived',
                                                 'source': ['seago1985'],
                                                 'basis': 'anagen over cycle, 54/151 in men and 22/84 in women. '
                                                          'The male value is carried because the atlas specimen '
                                                          'is an adult male reference.'},
                             'cycle_days': {'value': 151, 'band': [84, 151], 'tier': 'measured',
                                            'source': ['seago1985'], 'basis': 'male 151 d, female 84 d'}}),
    'arm': _vellus(growth={'anagen_fraction': {'value': .259, 'band': [.208, .259], 'tier': 'derived',
                                               'source': ['seago1985'],
                                               'basis': 'anagen over cycle, 28/108 in men and 22/106 in women'},
                           'cycle_days': {'value': 108, 'band': [106, 108], 'tier': 'measured',
                                          'source': ['seago1985'], 'basis': 'male 108 d, female 106 d'}}),
    'glabrous': {'shaft_radius_m': {'value': None, 'band': None, 'tier': None, 'source': [],
                                    'basis': 'no shaft: the field carries no follicles'},
                 'shaft_length_m': {'value': None, 'band': None, 'tier': None, 'source': [],
                                    'basis': 'no shaft'},
                 'shaft_bearing_fraction': {'value': 0.0, 'tier': 'measured',
                                            'source': ['tsai2022', 'rittie2013', 'fradette1995'],
                                            'basis': 'palmoplantar skin lacks hair follicles'},
                 'growth': {}, 'colour': None},
}
for _field in ('auricular', 'neck', 'chest', 'abdomen', 'gluteal', 'forearm', 'hand_dorsum', 'leg', 'foot_dorsum',
               'pubic', 'perineal'):
    MORPHOLOGY[_field] = _vellus()

_JON = 'jonsson2017'


def _density(value, tier, source, basis, sd=None, **extra):
    record = {'value': value, 'sd': sd, 'tier': tier, 'source': list(source), 'basis': basis}
    record.update(extra)
    return record


FIELD_EVIDENCE = {
    'scalp': {'follicle_density_cm2': _density(
        None, 'measured', ['leerunyakul2020'],
        'area-weighted over the four scalp areas leerunyakul2020 measured, mapped one to one onto the four '
        'Terminologia Anatomica head regions inside the hair footprint',
        subregion_cm2={'Frontal region': 154.3, 'Parietal region': 162.9, 'Temporal region': 133.7,
                       'Occipital region': 160.2},
        subregion_mapping='frontal->Frontal region, vertex->Parietal region, temporoparietal->Temporal region, '
                          'occipital->Occipital region',
        population_spread='133.7 to 233 per cm2 across cohorts: leerunyakul2020 Thai 133.7-162.9, alsharif2022 '
                          'Arabian 143.9-153.6, asfiya2025 Indian males 123.6-160.1, loussouarn2016 clusters '
                          '178/182/214 with a 153-233 range over 2249 subjects, birnbaum2018 148-230 across three '
                          'self-identified American groups. The atlas specimen has no declared population, so no '
                          'cohort is the right one and the spread is the honest uncertainty.',
        comparanda={'jimenez1999_occipital_hairs_cm2': [124, 200],
                    'visessiri2020_derived_cm2': 163.1, 'xu2017_pooled_all_sites_cm2': 140}),
        'notes': ['ihm/assembly/hair.py and data/measurements/hair/strand_reference.json currently carry '
                  'alsharif2022 for the scalp. This candidate carries leerunyakul2020 instead for one reason: it '
                  'is the only cohort that measured all four of the surface regions the scalp footprint is made '
                  'of, including the temporal, so the density can be area-weighted on the geometry rather than '
                  'averaged over regions the source never measured.']},
    'eyebrow': {'follicle_density_cm2': _density(
        133.95, 'measured', ['ali2023'], 'unburned contralateral control eyebrows', sd=38.38),
        'notes': ['a burn-reconstruction series of mostly young women, not a normative study. It is the only '
                  'per-cm2 eyebrow density retrievable; flament2019 measured eyebrow density in 3600 women across '
                  'six groups and publishes no number in its abstract.']},
    'eyelash': {'follicle_density_cm2': _density(
        None, None, [], 'a lash line is a strip a few millimetres wide; a per-area density would be a worse '
                        'quantity than the count, which is measured'),
        'follicle_count': {'value': 220.6, 'range': [180.0, 260.0], 'tier': 'measured', 'source': ['freitag2014'],
                           'basis': 'two upper lids at 72.1 lashes and two lower lids at 38.2, so 2 x (72.1 + '
                                    '38.2) = 220.6 lashes. The range brackets the count over the reported '
                                    'variation of the cohort rather than a stated SD, which the abstract does '
                                    'not give.'},
        'notes': ['aumond2018 states 90-160 upper and 75-80 lower lashes per lid, which is roughly twice the only '
                  'measured count. That is a review statement whose primaries were not read, so the measured '
                  'count is carried and the review range is recorded as the disagreement it is.',
                  'the footprint area is measured from the authored eyelash shells, but the count does not use '
                  'it: an areal density on a lash line would be an artefact of the strip width.']},
    'beard': {'follicle_density_cm2': _density(
        426.0, 'transferred', ['blume1991', 'todo2017'],
        'cheek vellus hair density in women, read second-hand in todo2017 quoting blume1991, transferred to the '
        'beard field of a male reference atlas', sd=37.7,
        transfer='from female cheek vellus hair to the whole beard field. Follicle number is roughly sex-invariant '
                 '(seago1985 measured no sex difference in follicle density on thigh or upper arm); what differs '
                 'is how many follicles produce a terminal shaft, and that fraction is not measured anywhere '
                 'retrieved. So this is a follicle count, not a beard-hair count.',
        caveats=['second hand: the primary abstract gives only forehead 439 and back 85, not the cheek',
                 'the quoted cheek SD 37.7 equals the quoted forehead-male SD, so the dispersion may be a '
                 'transcription artefact',
                 'floyd2018 and bouabbache2019 both measured beard areal density directly and neither publishes '
                 'it outside a paywall']),
        'notes': ['beard mass is absent, not estimated: the shaft diameter is measured at 76 +/- 7.4 um but no '
                  'length and no terminal fraction exist.']},
    'facial_vellus': {'follicle_density_cm2': _density(
        285.0, 'measured', [_JON], 'forehead, cyanoacrylate stripping, n=15', sd=84.1,
        comparanda={'otberg2004_forehead_second_hand': 292, 'blume1991_forehead': 439,
                    'note': 'blume1991 photographic counting gives 1.5 times the cyanoacrylate figure. Method, '
                            'not population.'}),
        'notes': ['the forehead measurement is carried over the whole non-beard facial field, including the '
                  'orbital, nasal, zygomatic and perioral regions that no source measured separately.']},
    'auricular': {'follicle_density_cm2': _density(
        None, None, [], 'no source measures the auricle or the mastoid region'),
        'notes': ['the field is geometrically resolved and populated by nothing. hwang1997 measured 30 body '
                  'regions in 74 Korean cadavers and is the source that would close this; its table is paywalled.']},
    'neck': {'follicle_density_cm2': _density(47.3, 'measured', [_JON], 'neck, n=15', sd=17.1)},
    'axillary': {'follicle_density_cm2': _density(
        11.05, 'measured', ['coelho2026'],
        'mean of the two trial arms, 11.1 and 11.0 visible hairs per cm2 in a 4 cm2 window', sd=.50,
        caveats=['a visible-hair count in 48 women at the baseline of a hair-removal trial, not a follicle '
                 'density, and low against every impression of axillary hair',
                 'pecoraro1971 is exactly the axillary density and growth-rate study needed and could not be '
                 'read through any of four routes']),
        'unresolved_reason': 'the Z-Anatomy atlas carries 256 Terminologia Anatomica body surface regions and no '
                             'axillary region among them, so there is no authored skin label to anchor the field '
                             'to. A nearest-canonical-structure rule over the axillary vessels and node groups '
                             'was tried and rejected: the surrounding muscles win every envelope face but 13, '
                             'giving 4.18 cm2 in the supraclavicular fossae, which is not the axilla. No '
                             'literature source gives an axillary hair-field area either. The footprint is '
                             'therefore absent rather than boxed.',
        'notes': ['this is the one owner-named field with no geometry. Its follicle count and mass are null and '
                  'the area it would have taken stays inside the chest and arm fields.']},
    'chest': {'follicle_density_cm2': _density(
        28.4, 'measured', [_JON], 'chest, n=15', sd=6.4,
        comparanda={'blume1991_second_hand': {'female': 53, 'male': 61}})},
    'abdomen': {'follicle_density_cm2': _density(16.9, 'measured', [_JON], 'lower abdomen, n=15', sd=6.2)},
    'back': {'follicle_density_cm2': _density(
        24.7, 'measured', [_JON], 'lower back, n=15', sd=5.9,
        comparanda={'otberg2004_second_hand': 29, 'blume1991': 85,
                    'blume1991_second_hand': {'female': 93, 'male': 77}})},
    'gluteal': {'follicle_density_cm2': _density(
        24.7, 'transferred', [_JON], 'transferred from the lower back, the nearest measured trunk site', sd=5.9,
        transfer='lower back to gluteal region, gluteal fold and hip region. No source measures gluteal skin.')},
    'arm': {'follicle_density_cm2': _density(45.9, 'measured', [_JON], 'upper arm, n=15', sd=14.5)},
    'forearm': {'follicle_density_cm2': _density(
        39.2, 'measured', [_JON], 'forearm, the combined n=138 sample rather than the n=15 pretest', sd=12.1,
        sex={'women': [41.4, 12.8], 'men': [35.5, 10.0],
             'statement': 'women higher than men, which the authors explain by body size and weight'})},
    'hand_dorsum': {'follicle_density_cm2': _density(
        39.2, 'transferred', [_JON], 'transferred from the forearm, the nearest measured site', sd=12.1,
        transfer='forearm to dorsum of hand, dorsal surfaces of the digits and the radial foveola. hwang1997 '
                 'measured the dorsum of the middle phalanx of the fingers as the sparsest of its 30 regions, so '
                 'this transfer is very likely an over-estimate for the digits.')},
    'thigh': {'follicle_density_cm2': _density(21.0, 'measured', [_JON], 'thigh, n=15', sd=8.0)},
    'leg': {'follicle_density_cm2': _density(15.6, 'measured', [_JON], 'calf, n=15', sd=3.6)},
    'foot_dorsum': {'follicle_density_cm2': _density(
        15.6, 'transferred', [_JON], 'transferred from the calf, the nearest measured site', sd=3.6,
        transfer='calf to dorsum of foot, dorsal surfaces of the toes, the foot borders and the metatarsal region')},
    'pubic': {'follicle_density_cm2': _density(
        None, None, [], 'no retrievable source gives a pubic follicle density'),
        'notes': ['astore1979, the normal trichogram of pubic hair, publishes direction only in its abstract and '
                  'is closed. The field is geometrically resolved from the atlas pubic-hair shell and carries no '
                  'population.']},
    'perineal': {'follicle_density_cm2': _density(
        None, None, [], 'no perineal or perianal hair measurement of any kind was found')},
    'glabrous': {'follicle_density_cm2': _density(
        0.0, 'measured', ['tsai2022', 'rittie2013', 'fradette1995'],
        'palmoplantar skin lacks hair follicles', sd=0.0,
        caveats=['the three sources cover palms and soles. This field also holds the perionyx and the lip '
                 'vermilion regions (tubercle of upper lip, labial commissure), for which no clean follicle-free '
                 'statement was retrieved. Those are 13.1 cm2 of the fields 913.6 cm2, so the zero is sourced '
                 'over 98.6 percent of the field by area and asserted over the remaining 1.4 percent.',
                 'the heel region, 177.2 cm2, is placed here on its plantar surface. Its posterior aspect is not '
                 'glabrous and this is an assignment, not a measurement.'])},
}

NOTES = {
    'reconciliation': {
        'question': 'what happens to body-dynamic-scalp-hair and body-dynamic-body-hair, and to the hair-strands '
                    'experiment that renders them',
        'answer': 'keep them alongside, as the materialization of these anatomical fields. Do not supersede and '
                  'do not merge.',
        'reasoning': [
            'they are different objects. The two dynamic structures are elastic strand populations with '
            'centrelines, radii and a bending law, produced by scripts/build_hair_strands.py for rendering and '
            'for a 10 Hz dynamics loop. The 21 records here are anatomical entities: a named field, a footprint '
            'on the measured envelope, a follicle population and a provenance chain. One is a materialization of '
            'the other, in the sense docs/research/SOFT_BODY_MATERIALIZATION.md gives that word.',
            'nothing double-counts, because they count different things. The dynamic structures carry '
            'total_simulated_guide_mass_kg of 1.563e-05 and 2.507e-08 kg, which is the mass of 64 and 32 '
            'simulated guide strands, not of a hair population. The mass in this candidate is the population '
            'mass. Summing the two would be wrong; they are not disjoint parts of one quantity.',
            'both point at the same canonical parent. The dynamic structures already declare '
            'canonical_entity_id body-bp3d-FJ2810, the skin. These fields declare a field_on connection to the '
            'same entity. Adding a materializes field id to each dynamic structure would bind them without '
            'changing either.'],
        'supersedes': [],
        'retains': ['body-dynamic-scalp-hair', 'body-dynamic-body-hair'],
        'existing_canonical_hair_entities': {
            'body-bp3d-FJ2813': 'hair of head, an authored hair-shell surface of 0.15202 m2 with role soft_organ. '
                                'Not superseded: it is the source geometry the scalp footprint is projected from. '
                                'Its role is wrong and is recorded as such in role_vocabulary rather than '
                                'silently reused.',
            'body-bp3d-FJ2815': 'pubic hair, an authored hair-shell surface of 0.01234 m2 with role soft_organ. '
                                'Same relationship to the pubic field.',
            'body-bp3d-FJ2812': 'eyebrow, 0.00148 m2, role soft_organ. The eyebrow field is anchored to the '
                                'Z-Anatomy Hairs of eyebrow shells, which are a different mesh.'},
        'display': 'the manifest fragment gives every field system "hair", matching the two dynamic structures, '
                   'so the Layers list shows one expandable Hair group holding 21 anatomical fields and 2 '
                   'materializations. The canonical entity records carry system "integumentary", because that is '
                   'the canonical vocabulary and there is no "hair" system in it. The split is deliberate: '
                   'display grouping and anatomical system are different questions.',
        'changes_needed_elsewhere': [
            'ihm/assembly/hair.py PRIORS is the jonsson2017 nine-site table with no attribution and slices on raw '
            'x/y/z bounds. Once these fields are promoted it is superseded by them and should be retired, not '
            'edited. This lane does not touch it.',
            'scripts/build_hair_strands.py scalp_mask is an atlas-coordinate box. It should read the scalp field '
            'footprint from this candidate instead. Not done here: it would change a canonical artifact.'],
    },
    'limits': [
        'The footprint is measured on one generic authored body. The follicle densities come from cohorts that '
        'are not that body and, in seven of eleven cases, not even the same continent.',
        'Nine of the density values come from one n=15 pretest. It is a real measurement with a real method, but '
        'it is fifteen people.',
        'Every non-scalp shaft radius and length is an assumed engineering scenario. The mass of those fields is '
        'therefore an assumption with a geometry attached, and its band is reported alongside it.',
        'In an adult male a substantial fraction of limb and trunk follicles bear terminal rather than vellus '
        'shafts. seago1985 measures that the sex difference acts through anagen duration and shaft length, not '
        'follicle number, but no retrievable source gives terminal shaft dimensions by site. The non-scalp mass '
        'here is therefore a lower bound.',
        'Follicle counts assume one shaft per follicle. Telogen follicles carry none; loussouarn2016 puts that at '
        '8 to 14 percent on the scalp and nothing measures it for the body sites.',
        'The circular section over-states the real elliptical section: berg2025 measures 2870 to 4680 um2 against '
        'the 5102 um2 a circle of the mean scalp diameter would have.',
        'Region assignment is nearest-sample over registered surface meshes, so a face on a boundary between two '
        'Terminologia regions goes to one of them by distance. The median assignment distance is reported.',
        'The registration that places the Z-Anatomy regions is the existing z_anatomy landmark fit, whose own '
        'held-out RMS is 4.67 mm. Region boundaries inherit that.',
        'Three owner-named fields carry no population: axillary has no footprint and no follicle density, pubic '
        'and perineal have footprints and no density.',
    ],
    'unsourced_after_search': [
        'axillary hair-field area, follicle density, count and length',
        'pubic follicle density, escutcheon area and hair length',
        'any perineal or perianal hair measurement',
        'beard follicle density by facial zone, and total beard hair count',
        'eyebrow hair count, and eyebrow hair length',
        'measured eyelash length in mm from a primary study',
        'auricular, gluteal, hand dorsum and foot dorsum follicle density',
        'terminal shaft diameter and length by body site',
        'total scalp follicle count in an adult from a primary measurement',
        'any study that has weighed the hair on a head or a body',
        'regional hair colour',
    ],
    'sources_identified_but_unread': {
        'hwang1997': '30 body regions in 74 Korean cadavers. The single highest-value target: it covers auricular, '
                     'gluteal, hand, foot, axillary and perineal regions that are absent here.',
        'otberg2004': 'seven-site density, orifice and shaft diameter table.',
        'mangelsdorf2006': 'the same seven sites in Asian and African-American cohorts.',
        'pecoraro1971': 'axillary density and growth rate by sex.',
        'astore1979': 'pubic trichogram.',
        'szabo1967': 'the classical regional integument survey.',
        'barman1965, saitoh1970': 'the classical adult trichogram and per-site cycle durations.',
        'floyd2018, bouabbache2019': 'measured beard and axillary areal densities, both paywalled.',
        'flament2019': 'eyebrow morphometry in 3600 women.',
    },
}

ASSUMPTION_LEDGER = [
    {'id': 'HAIR-FIELD-LEXICON', 'kind': 'spatial partition prior',
     'statement': 'Hair fields are named groups of Terminologia Anatomica body surface regions carried verbatim '
                  'from the Z-Anatomy atlas and registered by the existing z_anatomy landmark fit. Every outer '
                  'envelope face is assigned to the nearest region sample, so the fields partition the measured '
                  '1.78125 m2 envelope exactly, with no gap and no overlap. Grouping regions into fields is a '
                  'declared lexicon, not a measured boundary: a face between two regions goes to one of them by '
                  'distance, and region boundaries inherit the 4.67 mm held-out RMS of the registration. Scalp, '
                  'eyebrow, eyelash and pubic footprints override the lexicon where the atlas authors the hair '
                  'itself, within a 6 mm standoff whose sensitivity is tabulated.'},
    {'id': 'HAIR-DENSITY-TRANSFER', 'kind': 'population transfer',
     'statement': 'Follicle density, shaft diameter, shaft length, growth rate and phase fractions come from '
                  'published cohorts that are not this specimen, and in most cases not the same continent. Each '
                  'value carries a measured, transferred or assumed tier, the cohort that produced it, and the '
                  'identifier check that verified it. Values with no retrievable source are absent, not filled. '
                  'Every non-scalp shaft radius and length is an assumed engineering scenario, so the mass of '
                  'those fields is an assumption with geometry attached and is reported with its band.'},
]
NOTES['assumption_ledger'] = ASSUMPTION_LEDGER
NOTES['observed_repo_state'] = {
    'hair_strands_gate': 'ihm.app.experiments.read_experiment(root, "hair-strands") currently raises "Hair source '
                         'changed; rebuild materialization". The differing hash is app/src/hair_dynamics.js, '
                         'which a concurrent front-end lane is rewriting. Every other source in that fragment '
                         'still matches. This lane did not touch app/src and does not fix it; it is recorded so '
                         'the failure is not attributed to this candidate.',
}
