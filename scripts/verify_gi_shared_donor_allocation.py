#!/usr/bin/env python3
"""Isolated allocation design regression; Python algebra, not corrected native code."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import tempfile

ROOT = Path(__file__).resolve().parents[1]
NATIVE_RECEIPT = ROOT/'data/derived/audits/gi-shared-donor-gepygiov/report.json'


def allocate(masses, ve, ei, corrected):
    if len(masses) != 3 or not all(math.isfinite(v) for v in [*masses, ve, ei]):
        raise ValueError('Allocation requires three finite donor masses and finite transfers')
    if any(m < 0 for m in masses):
        raise ValueError('Starting donor mass must be nonnegative')
    ve = min(ve, masses[0]) if ve > 0 else -min(-ve, masses[1])
    ei = min(ei, masses[1]) if ei > 0 else -min(-ei, masses[2])
    if corrected:
        out_to_v = max(-ve, 0.)
        out_to_i = max(ei, 0.)
        shared_requested = out_to_v + out_to_i
        if shared_requested > masses[1]:
            scale = masses[1] / shared_requested
            if ve < 0: ve *= scale
            if ei > 0: ei *= scale
    return ve, ei


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['original', 'corrected'], default='corrected')
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (30, 30))
    out = Path(tempfile.mkdtemp(prefix='gi-donor-allocation-', dir=ROOT/'data/derived/audits'))
    native = json.loads(NATIVE_RECEIPT.read_text())
    checks = {}; rows = {}
    # The same admissibility gate is deliberately false for the retained native calls.
    for layout in ['leaf', 'parent']:
        row = native['cases']['shared_scarce_'+layout]
        checks['retained_native_'+layout+'_fails_admissibility'] = not (
            row['native_nonnegative'] and abs(row['native_mass_residual_ug']) <= 1e-9)
    cases = [
        ('scarce_equal', [0.,1000.,0.], -1000., 1000., [500.,0.,500.]),
        ('scarce_unequal', [0.,1000.,0.], -800., 400., [2000/3,0.,1000/3]),
        ('ample', [0.,1000.,0.], -25., 20., [25.,955.,20.]),
        ('zero', [1.,2.,3.], 0., 0., [1.,2.,3.]),
        ('empty_E', [1.,0.,3.], -5., 8., [1.,0.,3.]),
        ('forward_chain', [100.,50.,0.], 100., 100., [0.,100.,50.]),
        ('reverse_chain', [0.,50.,100.], -100., -100., [50.,100.,0.]),
        ('two_inflows', [100.,0.,100.], 100., -100., [0.,200.,0.]),
    ]
    for name, initial, requested_ve, requested_ei, expected in cases:
        ve, ei = allocate(initial, requested_ve, requested_ei, args.mode == 'corrected')
        delta = [-ve, ve-ei, ei]
        for parent in [False, True]:
            final = [m+d for m,d in zip(initial, delta)]
            if parent: final[1] = max(final[1], 0.)  # held parent debit cap
            key = name+'/'+('parent' if parent else 'leaf')
            checks[key+'/nonnegative'] = all(math.isfinite(m) and m >= -1e-9 for m in final)
            checks[key+'/conservation'] = abs(sum(final)-sum(initial)) <= 1e-9
            checks[key+'/owner_incidence'] = all(abs(b-a-d) <= 1e-9 for a,b,d in zip(initial,final,delta))
            checks[key+'/expected_allocation'] = all(math.isclose(a,b,abs_tol=1e-9,rel_tol=1e-12) for a,b in zip(final,expected))
            rows[key] = {'initial_ug':initial,'requested_ve_ei_ug':[requested_ve,requested_ei],
                         'allocated_ve_ei_ug':[ve,ei],'final_ug':final,'mass_residual_ug':sum(final)-sum(initial)}
    report = {'passed':all(checks.values()),'mode':args.mode,'checks':checks,'rows':rows,
              'scope':'Independent Python allocation fixture; corrected native source was not compiled or executed.',
              'native_receipt':str(NATIVE_RECEIPT),'native_receipt_sha256':hashlib.sha256(NATIVE_RECEIPT.read_bytes()).hexdigest(),
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'receipt':str(out/'report.json'),'passed':report['passed'],
                      'check_count':len(checks),'failed_checks':[k for k,v in checks.items() if not v]},indent=2))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
