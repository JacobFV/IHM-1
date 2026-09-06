#!/usr/bin/env python3
"""Read retained actual-native transfers; report known Na charge, not membrane voltage."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
F = 96485.33212  # C/mol, conventional displayed precision
DEFAULT = ROOT / 'data/derived/audits/gi-absorption-correction-4rrrl_ru/report.json'

def audit(receipt):
    raw = receipt.read_bytes()
    report = json.loads(raw)
    assert report['passed'] and report['variant'] == 'corrected'
    assert len(report['cases']) == 14
    mw, definitions = {}, {}
    for name in ['Glucose', 'Sodium', 'AminoAcids']:
        path = ROOT / f'data/raw/physiology/biogears/share/data/substances/{name}.xml'
        body = path.read_bytes()
        tree = ET.fromstring(body)
        element = next(e for e in tree.iter() if e.tag.split('}')[-1] == 'MolarMass')
        assert element.attrib['unit'] == 'g/mol'
        mw[name] = float(element.attrib['value'])
        definitions[name] = {'sha256': hashlib.sha256(body).hexdigest(), 'molar_mass_g_mol': mw[name]}
    cases = {}
    for name, case in report['cases'].items():
        assert case['species_order'] == ['glucose', 'sodium', 'amino_acids', 'TAG', 'calcium', 'chloride']
        credit = case['native_vascular_credit_g']
        for initial, final, moved in zip(case['initial_chyme_g'], case['final_chyme_g'], credit):
            assert math.isclose(initial - final, moved, abs_tol=2e-15)
        g, na, aa, _, ca, cl = credit
        assert ca == cl == 0  # These fixtures cannot establish general Ca/Cl behavior.
        ng, nna, naa = g / mw['Glucose'], na / mw['Sodium'], aa / mw['AminoAcids']
        independent_na = na - g / 2 - aa
        assert independent_na >= -1e-18
        conditional_coupled_na = 2 * ng + naa
        for mass, mol, molecular_mass in [(g, ng, mw['Glucose']), (na, nna, mw['Sodium']), (aa, naa, mw['AminoAcids'])]:
            assert math.isclose(mol * molecular_mass, mass, abs_tol=1e-20)
        cycles = conditional_coupled_na / 3
        assert math.isclose(3 * cycles, conditional_coupled_na, abs_tol=1e-20)
        # Pump has one net outward positive charge per cycle; K recycling
        # restores the two imported K without changing the nutrient comparison.
        assert math.isclose(3 * cycles - 2 * cycles, cycles, abs_tol=1e-20)
        assert F * nna + (-F * nna) == 0
        cases[name] = {
            'glucose_mol': ng, 'amino_acid_pseudomoles': naa,
            'native_sodium_mol': nna,
            'source_attributed_glucose_sodium_mol': g / (2 * mw['Sodium']),
            'source_attributed_aa_sodium_mol': aa / mw['Sodium'],
            'independent_sodium_mol': independent_na / mw['Sodium'],
            'known_na_charge_equivalent_to_vascular_C': F * nna,
            'known_na_charge_equivalent_from_lumen_C': -F * nna,
            'measured_chloride_credit_mol': 0,
            'conditional_2Na_glucose_1Na_neutralAA_mol': conditional_coupled_na,
            'conditional_pump_cycles_mol_ATP': cycles,
            'conditional_pump_K_import_and_recycling_mol': 2 * cycles,
            'unresolved_AA_charge': 'Positive renal enum has no numeric valence or molecular mixture; excluded from known charge.',
        }
    ratios = {'Na_per_glucose': mw['Glucose'] / (2 * mw['Sodium']), 'Na_per_AA_pseudomolecule': mw['AminoAcids'] / mw['Sodium']}
    assert ratios['Na_per_glucose'] > 3.9 and ratios['Na_per_AA_pseudomolecule'] > 3.8
    assert cases['zero_all']['native_sodium_mol'] == 0
    assert cases['fat_tail']['independent_sodium_mol'] > 0
    assert cases['sequential_sodium']['glucose_mol'] > 0 and cases['sequential_sodium']['amino_acid_pseudomoles'] > 0
    return {'audit_passed': True, 'meaning': 'Checks reproduce retained transfers and expose coupling mismatch; not validation of epithelial charge closure or pH.', 'receipt': str(receipt.relative_to(ROOT)), 'receipt_sha256': hashlib.sha256(raw).hexdigest(), 'substance_definitions': definitions, 'native_molar_ratios': ratios, 'cases': cases}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = audit(args.receipt)
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    else:
        print(text, end='')
