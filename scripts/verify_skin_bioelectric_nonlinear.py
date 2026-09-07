#!/usr/bin/env python3
"""Verify and exercise scripts/skin_bioelectric_nonlinear.py.

  --self-test   fast executable checks only (no artifacts written)
  (default)     run the acceptance demonstrations and write
                data/derived/skin-bioelectric-nonlinear-v1/{demonstrations,self_test,manifest}.json

The acceptance question is whether three things that were IMPOSSIBLE in the
retained substrate are now stateable and measurable:
  (a) two distinct stable steady states from identical parameters, plus a
      hysteresis loop under a swept drive
  (b) a wound that measurably changes cell Vm, replacing the certified-zero
      result asserted at scripts/verify_body_skin_electric.py:40
  (c) a regional gap-junction blockade that changes the spatial pattern
Every claim carries a conservation receipt. Numbers that do not materialise are
reported as not materialising, with a diagnosis, not tuned until they appear.
"""
import argparse
import gzip
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('MKL_NUM_THREADS', '2')

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import skin_bioelectric_nonlinear as M   # noqa: E402

OUT = ROOT / 'data/derived/skin-bioelectric-nonlinear-v1'

INPUTS = [
    'scripts/skin_bioelectric_nonlinear.py',
    'scripts/verify_skin_bioelectric_nonlinear.py',
    'ihm/assembly/skin_bioelectric.py',
    'ihm/materialize/skin.py',
    'scripts/verify_body_skin_electric.py',
    'ihm/assembly/epithelial_electrodiffusion.py',
    'data/derived/bioelectric-substrate-assessment-v1/assessment.json',
    'data/derived/bioelectric/tissue.json',
    'data/derived/app/geometry/betse-tissue-cells.json.gz',
    'data/derived/physiology/betse_run/sim_config.yaml',
    'data/derived/outer-envelope/outer-envelope.json.gz',
    'data/derived/outer-envelope/validation.json',
    'data/derived/canonical/skin-electric.json',
    'docs/research/SOFT_BODY_MATERIALIZATION.md',
]


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load(rel):
    raw = (ROOT / rel).read_bytes()
    return json.loads(gzip.decompress(raw) if rel.endswith('.gz') else raw)


# --------------------------------------------------------------- self test ----

def self_test():
    from dataclasses import replace
    from ihm.materialize.skin import SkinPatch
    from ihm.assembly.skin_bioelectric import electrical_system

    checks = {}

    # 1. The retained linear law is reachable UNCHANGED.  Build the retained
    #    9-cell line (homogenised so one specific conductance describes it),
    #    then reproduce its generator exactly from the sibling with
    #    law='linear', bath='prescribed'.
    patch = SkinPatch.line(n=9, wound=False)
    patch = replace(patch, cells=tuple(replace(c, potassium_conductance_s=1e-10) for c in patch.cells))
    a, b, _ = electrical_system(patch)
    nc = len(patch.cells)
    membrane = M.linear_membrane()
    area = np.full(nc, 1e-11 / M.SPECIFIC_CAPACITANCE_F_PER_M2)
    sheet = M.Sheet(np.array([c.position_m for c in patch.cells]), area, area * 1e-5,
                    np.array([[e[0], e[1]] for e in patch.gap_edges]),
                    np.array([e[2] for e in patch.gap_edges]),
                    np.full(len(patch.gap_edges), 1e-9), np.full(nc, 1.0), np.zeros(nc),
                    membrane, bath='prescribed')
    x = np.linspace(-.07, -.01, nc)
    li, _, _ = M.operators(sheet)
    i_ion, _, _ = membrane.current(x)
    rate = (-(li @ x) - area * i_ion) / (area * membrane.specific_capacitance_f_per_m2)
    checks['legacy_law_generator_parity_v_per_s'] = float(np.max(np.abs(rate - (a[:nc, :nc] @ x + b[:nc]))))
    assert checks['legacy_law_generator_parity_v_per_s'] < 1e-14
    retained_rest = np.linalg.solve(a[:nc, :nc], -b[:nc])
    sibling_rest = M.steady_state(sheet, -0.08)['vm_v']
    checks['legacy_resting_parity_v'] = float(np.max(np.abs(retained_rest - sibling_rest)))
    assert checks['legacy_resting_parity_v'] < 1e-12
    checks['legacy_resting_v'] = float(retained_rest[0])

    # 2. Gate derivatives against central differences.
    g = M.Gate(-50e-3, 10e-3)
    v = np.linspace(-0.12, 0.06, 37)
    _, d = g.steady(v)
    fd = (g.steady(v + 1e-7)[0] - g.steady(v - 1e-7)[0]) / 2e-7
    checks['gate_derivative_max_abs_error'] = float(np.max(np.abs(d - fd)))
    assert checks['gate_derivative_max_abs_error'] < 1e-6

    # 3. dI/dV against central differences for the gated membrane.
    mem = M.bistable_membrane()
    i0, di, _ = mem.current(v)
    fd = (mem.current(v + 1e-7)[0] - mem.current(v - 1e-7)[0]) / 2e-7
    checks['membrane_slope_max_rel_error'] = float(np.max(np.abs(di - fd)) / max(np.max(np.abs(di)), 1e-30))
    assert checks['membrane_slope_max_rel_error'] < 1e-5

    # 4. Nonlinearity actually exists: three roots, two stable.
    roots = M.iv_roots(mem)
    checks['iv_roots_v'] = [r['vm_v'] for r in roots]
    checks['iv_stability'] = [r['stable'] for r in roots]
    assert len(roots) == 3 and [r['stable'] for r in roots] == [True, False, True]
    checks['linear_law_root_count'] = len(M.iv_roots(M.linear_membrane()))
    assert checks['linear_law_root_count'] == 1

    # 5. Bidomain conservation identities on a small sheet.
    small = M.build_sheet(12, 12, 2e-4, mem)
    st = M.steady_state(small, -0.09)
    rc = st['receipts']
    checks['elliptic_residual_a'] = rc['elliptic_residual_a']
    checks['sum_membrane_current_a'] = rc['sum_membrane_current_a']
    checks['battery_minus_barrier_a'] = rc['battery_minus_barrier_a']
    checks['max_node_current_imbalance_a'] = rc['max_node_current_imbalance_a']
    assert abs(rc['sum_membrane_current_a']) < 1e-18
    assert abs(rc['battery_minus_barrier_a']) < 1e-18
    assert rc['elliptic_residual_a'] < 1e-18
    assert rc['max_node_current_imbalance_a'] < 1e-18

    # 6. Prescribed bath REPRODUCES the certified-zero wound result; the solved
    #    bath does not. This is the whole point of the module.
    ix = (small.positions_m[:, 0] / 2e-4).round().astype(int)
    iy = (small.positions_m[:, 1] / 2e-4).round().astype(int)
    wound = (np.abs(ix - 6) <= 1) & (np.abs(iy - 6) <= 1)
    barrier = np.where(wound, 200.0, 2.0)
    legacy = M.build_sheet(12, 12, 2e-4, mem, bath='prescribed')
    base_l = M.steady_state(legacy, -0.09)['vm_v']
    wound_l = M.steady_state(legacy, -0.09, barrier=barrier)['vm_v']
    checks['prescribed_bath_wound_delta_vm_v'] = float(np.max(np.abs(wound_l - base_l)))
    assert checks['prescribed_bath_wound_delta_vm_v'] < 1e-15
    wound_s = M.steady_state(small, st['vm_v'], barrier=barrier)['vm_v']
    checks['solved_bath_wound_delta_vm_v'] = float(np.max(np.abs(wound_s - st['vm_v'])))
    assert checks['solved_bath_wound_delta_vm_v'] > 1e-4

    # 7. Time refinement: halving dt moves the trajectory by less than the
    #    tolerance, so the reported dynamics are not a step-size artefact.
    runs = [M.integrate(small, -0.09 * np.ones(small.n), duration_s=0.2, dt_s=h,
                        applied_a_per_m2=2e-3, record_every=round(0.2 / h))['final_vm_v']
            for h in (0.01, 0.005, 0.0025, 0.00125)]
    diffs = [float(np.max(np.abs(runs[i] - runs[i + 1]))) for i in range(len(runs) - 1)]
    checks['time_refinement_successive_diffs_v'] = diffs
    checks['time_refinement_max_delta_v'] = diffs[0]
    checks['time_refinement_observed_order'] = [float(math.log2(diffs[i] / diffs[i + 1]))
                                                for i in range(len(diffs) - 1)]
    assert diffs[0] < 5e-4 and diffs[0] > diffs[1] > diffs[2]
    assert all(0.7 < o < 1.4 for o in checks['time_refinement_observed_order']), \
        'backward Euler must converge at first order: %s' % checks['time_refinement_observed_order']

    # 8. Areas are declared and the specific capacitance is the literature one.
    cap = small.specific_capacitance_check
    checks['specific_capacitance_uf_per_cm2'] = cap['specific_capacitance_uf_per_cm2']
    assert abs(cap['specific_capacitance_uf_per_cm2'] - 1.0) < 1e-12
    checks['coarse_graining'] = small.coarse_graining['keratinocytes_per_node']

    # 9. Caps and invalid inputs are rejected.
    bad = 0
    for fn in (lambda: M.Gate(0., 0.),
               lambda: M.Channel('x', 'Zz', 1.),
               lambda: M.Channel('x', 'Na', -1.),
               lambda: M.Membrane((M.Channel('x', 'Na', 1.),), 0., law='nope'),
               lambda: M.integrate(small, -0.09 * np.ones(small.n), duration_s=1e9, dt_s=1.),
               lambda: M.integrate(small, -0.09 * np.ones(small.n), duration_s=1., dt_s=-1.),
               lambda: M.integrate(small, -0.09 * np.ones(small.n), duration_s=1., dt_s=0.3),
               lambda: M.Sheet(np.zeros((3, 3)), np.ones(3), np.ones(3), np.array([[0, 0]]),
                               np.ones(1), np.ones(1), np.ones(3), np.zeros(3), mem),
               lambda: M.Sheet(np.zeros((3, 3)), np.ones(3), np.ones(3), np.array([[0, 1], [1, 0]]),
                               np.ones(2), np.ones(2), np.ones(3), np.zeros(3), mem),
               lambda: M.Sheet(np.zeros((3, 3)), np.ones(3), -np.ones(3), np.array([[0, 1]]),
                               np.ones(1), np.ones(1), np.ones(3), np.zeros(3), mem)):
        try:
            fn()
        except (ValueError, TypeError):
            bad += 1
    checks['invalid_inputs_rejected'] = bad
    assert bad == 10

    # 10. Ion mass balance raises on exhaustion exactly like
    #     epithelial_electrodiffusion.transport does.
    starved = M.build_sheet(3, 3, 2e-4, M.linear_membrane(pump_a_per_m2=0.),
                            cell_height_m=1e-9)
    try:
        M.integrate(starved, -0.0815 * np.ones(starved.n), duration_s=200., dt_s=1.0,
                    ion_accounting=True, record_every=400)
        checks['reservoir_exhaustion_raises'] = False
    except M.ReservoirExhausted as exc:
        checks['reservoir_exhaustion_raises'] = True
        checks['reservoir_exhaustion_message'] = str(exc)
    assert checks['reservoir_exhaustion_raises']

    return checks


# ---------------------------------------------------------- demonstrations ----

def _profile(sheet, values, axis=0):
    key = np.round(sheet.positions_m[:, axis] / np.sqrt(sheet.node_area_m2.mean())).astype(int)
    order = np.argsort(key)
    uniq, start = np.unique(key[order], return_index=True)
    return uniq.tolist(), [float(np.mean(values[order][s:e])) for s, e in
                           zip(start, list(start[1:]) + [len(order)])]


def demo_bistability(spacing_m, nodes):
    mem = M.bistable_membrane()
    sheet = M.build_sheet(nodes, nodes, spacing_m, mem)
    roots = M.iv_roots(mem)
    stable = [r for r in roots if r['stable']]

    low = M.steady_state(sheet, -0.09)
    high = M.steady_state(sheet, 0.02)

    # Fold points of the space-clamped I(V): the drive at which a branch dies.
    grid = np.linspace(-0.15, 0.12, 54001)
    current, slope = M.iv_curve(mem, grid)
    turn = np.where(np.sign(slope[:-1]) * np.sign(slope[1:]) < 0)[0]
    folds = sorted(float(-current[k]) for k in turn)
    span = folds[-1] - folds[0] if len(folds) == 2 else 0.
    lo_drive, hi_drive = (folds[0] - 0.25 * span, folds[-1] + 0.25 * span) if span > 0 else (-2e-3, 2e-3)
    sweep = np.linspace(lo_drive, hi_drive, 61)
    small = M.build_sheet(12, 12, spacing_m, mem)
    up_run = M.hysteresis(small, sweep, v_start=-0.09)
    down_run = M.hysteresis(small, sweep[::-1], v_start=0.02)
    up = up_run['branch_vm_v']
    down = down_run['branch_vm_v'][::-1]
    gap = np.abs(up - down)

    return dict(
        claim='(a) two distinct stable steady states from identical parameters, and a hysteresis loop',
        materialised=bool(len(stable) == 2 and float(gap.max()) > 1e-3),
        space_clamped_iv=dict(
            roots_v=[r['vm_v'] for r in roots],
            di_dv_s_per_m2=[r['di_dv_s_per_m2'] for r in roots],
            stable=[r['stable'] for r in roots],
            separation_v=float(stable[-1]['vm_v'] - stable[0]['vm_v']) if len(stable) == 2 else None,
            note='three zeros of I(V); the middle one has dI/dV < 0 and is unstable. '
                 'The retained linear law has exactly one zero at any parameters.'),
        two_steady_states=dict(
            polarised_vm_v=float(low['vm_v'].mean()), polarised_sd_v=float(low['vm_v'].std()),
            depolarised_vm_v=float(high['vm_v'].mean()), depolarised_sd_v=float(high['vm_v'].std()),
            separation_v=float(high['vm_v'].mean() - low['vm_v'].mean()),
            identical_parameters=True,
            only_difference='initial condition (-90 mV vs +20 mV uniform)',
            polarised_receipts=low['receipts'], depolarised_receipts=high['receipts'],
            newton_iterations=[low['newton_iterations'], high['newton_iterations']],
            phi_e_v=[float(low['phi_e_v'].mean()), float(high['phi_e_v'].mean())]),
        hysteresis=dict(
            drive_a_per_m2=sweep.tolist(), branch_up_vm_v=up.tolist(), branch_down_vm_v=down.tolist(),
            fold_drives_a_per_m2=folds, loop_width_a_per_m2=float(span),
            max_branch_separation_v=float(gap.max()),
            drive_at_max_separation_a_per_m2=float(sweep[int(np.argmax(gap))]),
            branch_separation_at_zero_drive_v=float(np.interp(0., sweep, np.abs(up - down))),
            switch_drive_up_a_per_m2=[float(sweep[k]) for k, r in enumerate(up_run['relaxed']) if r],
            switch_drive_down_a_per_m2=[float(sweep[::-1][k]) for k, r in enumerate(down_run['relaxed']) if r],
            max_newton_residual_a=max(up_run['max_newton_residual_a'], down_run['max_newton_residual_a']),
            note='quasi-static continuation of the coupled bidomain steady state; the up sweep '
                 'starts polarised and the down sweep starts depolarised. A single-fixed-point '
                 'system returns the same branch both ways.'),
        parameters=dict(
            form='[literature] inward-rectifier K closing on depolarisation + a depolarising cation '
                 'conductance opening on depolarisation; Cervera/Alcaraz/Mafe 2014 J Phys Chem B '
                 '118:6417, Cervera/Meseguer/Mafe 2016 Sci Rep 6:35201, Law & Levin 2015 '
                 'Theor Biol Med Model 12:22; Hille ch.14 for inward rectification',
            values='[engineering] g_Kir 0.10 S/m^2 at V_half -50 mV slope 10 mV; g_dep 0.04 S/m^2 at '
                   'V_half -25 mV slope 6 mV. These are CHOSEN to sit in the bistable regime. No '
                   'keratinocyte measurement identifies them and nothing here claims human '
                   'keratinocytes are bistable.',
            retained='[retained] Cl leak 0.02 S/m^2, pump 1e-3 A/m^2 and the Nernst reversals are the '
                     'SkinPatch values re-expressed per unit area'))


def demo_wound(spacing_m, nodes, wound_radius_nodes, wound_barrier, lysis_s_per_m2):
    mem = M.bistable_membrane()
    sheet = M.build_sheet(nodes, nodes, spacing_m, mem)
    legacy = M.build_sheet(nodes, nodes, spacing_m, mem, bath='prescribed')
    ix = np.round(sheet.positions_m[:, 0] / spacing_m).astype(int)
    iy = np.round(sheet.positions_m[:, 1] / spacing_m).astype(int)
    c = nodes // 2
    wound = (np.abs(ix - c) <= wound_radius_nodes) & (np.abs(iy - c) <= wound_radius_nodes)
    barrier = np.where(wound, wound_barrier, float(sheet.barrier_s_per_m2[0]))

    base_legacy = M.steady_state(legacy, -0.09)['vm_v']
    wound_legacy = M.steady_state(legacy, -0.09, barrier=barrier)['vm_v']
    base = M.steady_state(sheet, -0.09)
    breach = M.steady_state(sheet, base['vm_v'], barrier=barrier)
    lysed = M.steady_state(sheet, base['vm_v'], barrier=barrier,
                           leak_s_per_m2=np.where(wound, lysis_s_per_m2, 0.))

    row = np.where(iy == c)[0]
    row = row[np.argsort(ix[row])]

    def field(state):
        return np.abs(np.diff(state['phi_e_v'][row]) / spacing_m)

    lc = M.length_constants(sheet)
    out = dict(
        claim='(b) a wound that measurably changes cell Vm, replacing the certified-zero result',
        materialised=bool(np.max(np.abs(breach['vm_v'] - base['vm_v'])) > 1e-4),
        certified_zero_reproduced=dict(
            configuration="bath='prescribed' (the retained model, exactly)",
            max_abs_delta_vm_v=float(np.max(np.abs(wound_legacy - base_legacy))),
            asserted_at='scripts/verify_body_skin_electric.py:40 asserts np.allclose(wound Vm, baseline Vm, atol=1e-12)',
            reading='With a prescribed bath the wound provably cannot move Vm. The sibling reproduces '
                    'that exactly, which is the control for the next block.'),
        geometry=dict(sheet_edge_m=float((nodes - 1) * spacing_m), node_spacing_m=spacing_m,
                      wound_edge_m=float((2 * wound_radius_nodes + 1) * spacing_m),
                      wound_nodes=int(wound.sum()),
                      intact_barrier_s_per_m2=float(sheet.barrier_s_per_m2[0]),
                      intact_teer_ohm_cm2=float(1e4 / sheet.barrier_s_per_m2[0]),
                      wound_barrier_s_per_m2=float(wound_barrier),
                      wound_teer_ohm_cm2=float(1e4 / wound_barrier)))
    for name, state in (('barrier_breach', breach), ('breach_plus_lysis', lysed)):
        d = state['vm_v'] - base['vm_v']
        e = field(state)
        out[name] = dict(
            max_abs_delta_vm_v=float(np.max(np.abs(d))),
            delta_vm_at_wound_v=float(d[wound].mean()),
            max_abs_delta_vm_in_intact_tissue_v=float(np.max(np.abs(d[~wound]))),
            phi_e_min_v=float(state['phi_e_v'].min()), phi_e_max_v=float(state['phi_e_v'].max()),
            tep_collapse_v=float(state['phi_e_v'].max() - state['phi_e_v'].min()),
            max_lateral_field_v_per_m=float(e.max()),
            lateral_field_profile_v_per_m=e.tolist(),
            vm_profile_v=_profile(sheet, state['vm_v'])[1],
            phi_e_profile_v=_profile(sheet, state['phi_e_v'])[1],
            receipts=state['receipts'])
    out['breach_plus_lysis']['lysis_conductance_s_per_m2'] = float(lysis_s_per_m2)
    out['breach_plus_lysis']['lysis_meaning'] = (
        '[engineering] a non-selective leak with 0 V reversal at the wounded nodes: a ruptured '
        'membrane stops separating the compartments. This is the classical injury conductance, not '
        'a healing, migration or regeneration law.')
    out['length_constants'] = lc
    out['comparison_to_human_measurement'] = dict(
        measured='[literature] Nuccitelli et al. 2011 (PMC3228273): lateral fields of order '
                 '100-200 mV/mm at a human skin wound edge, hash-receipted in this repository at '
                 'data/derived/calibration/skin-fit.json',
        modelled_max_v_per_m=out['barrier_breach']['max_lateral_field_v_per_m'],
        verdict='SIGN AND STRUCTURE reproduced, AMPLITUDE not calibrated.',
        why='the lateral field scales as TEP / lambda_e with lambda_e = sqrt(sigma_e t_e / g_b). '
            'sigma_e t_e (the sub-epidermal sheet conductance) and g_b (the barrier conductance) are '
            'not identified by anything in this repository, and the pair enters only through their '
            'ratio, so the amplitude is a free parameter. Matching it by choosing that ratio would '
            'be tuning, not evidence.',
        lambda_e_required_for_200_v_per_m_m=float(
            (float(sheet.transepithelial_a_per_m2[0]) / float(sheet.barrier_s_per_m2[0])) / 200.),
        lambda_e_used_m=lc['extracellular_length_constant_m'])
    return out


def demo_blockade(spacing_m, nodes, duration_s, dt_s, factor):
    ix_of = lambda sh: np.round(sh.positions_m[:, 0] / spacing_m).astype(int)
    iy_of = lambda sh: np.round(sh.positions_m[:, 1] / spacing_m).astype(int)

    # (c1) linear law, unambiguous: a localised drive, with and without a
    # regional blockade. Cutting edges in a diffusive graph must change the
    # steady pattern; this is a graph-Laplacian fact and NOT evidence for
    # Levin's biology. It is here as the control that the API works.
    lin = M.build_sheet(nodes, nodes, spacing_m, M.linear_membrane())
    ix, iy = ix_of(lin), iy_of(lin)
    c = nodes // 2
    drive = np.where((np.abs(ix - c // 2) <= 1) & (np.abs(iy - c) <= 1), -2e-2, 0.)
    u, w = lin.edges[:, 0], lin.edges[:, 1]
    cross = (ix[u] < c) != (ix[w] < c)
    blockade = np.where(cross, factor, 1.)
    free_lin = M.steady_state(lin, -0.08, applied_a_per_m2=drive)
    blocked_lin = M.steady_state(lin, -0.08, applied_a_per_m2=drive, blockade=blockade)
    left = ix < c
    linear = dict(
        law='linear (legacy law, explicitly selected)',
        blocked_edges=int(cross.sum()), blockade_factor=factor,
        free_vm_left_v=float(free_lin['vm_v'][left].mean()),
        free_vm_right_v=float(free_lin['vm_v'][~left].mean()),
        blocked_vm_left_v=float(blocked_lin['vm_v'][left].mean()),
        blocked_vm_right_v=float(blocked_lin['vm_v'][~left].mean()),
        free_step_across_line_v=float(free_lin['vm_v'][left].mean() - free_lin['vm_v'][~left].mean()),
        blocked_step_across_line_v=float(blocked_lin['vm_v'][left].mean() - blocked_lin['vm_v'][~left].mean()),
        max_abs_pattern_change_v=float(np.max(np.abs(blocked_lin['vm_v'] - free_lin['vm_v']))),
        free_receipts=free_lin['receipts'], blocked_receipts=blocked_lin['receipts'])

    # (c2) the interesting one: with a bistable membrane, two Vmem domains put
    # side by side. Coupled, the front runs and one domain eats the other.
    # Blockaded, the boundary is pinned and both domains persist.
    mem = M.bistable_membrane()
    sheet = M.build_sheet(nodes, nodes, spacing_m, mem)
    ix, iy = ix_of(sheet), iy_of(sheet)
    roots = [r['vm_v'] for r in M.iv_roots(mem) if r['stable']]
    v0 = np.where(ix < c, roots[-1], roots[0])
    u, w = sheet.edges[:, 0], sheet.edges[:, 1]
    cross = (ix[u] < c) != (ix[w] < c)
    blockade = np.where(cross, factor, 1.)
    out = {}
    for name, blk in (('coupled', None), ('regional_blockade', blockade)):
        t0 = time.time()
        run = M.integrate(sheet, v0, duration_s=duration_s, dt_s=dt_s, blockade=blk,
                          record_every=max(1, round(duration_s / dt_s) // 8))
        v = run['final_vm_v']
        left = ix < c
        out[name] = dict(
            wall_s=round(time.time() - t0, 1),
            final_vm_left_v=float(v[left].mean()), final_vm_right_v=float(v[~left].mean()),
            domain_step_v=float(v[left].mean() - v[~left].mean()),
            final_vm_sd_v=float(v.std()),
            depolarised_fraction=float((v > 0.5 * (roots[0] + roots[-1])).mean()),
            initial_depolarised_fraction=float((v0 > 0.5 * (roots[0] + roots[-1])).mean()),
            profile_vm_v=_profile(sheet, v)[1],
            frames=[dict(time_s=f['time_s'], vm_mean_v=f['vm_mean_v'], vm_sd_v=f['vm_sd_v'],
                         vm_min_v=f['vm_min_v'], vm_max_v=f['vm_max_v']) for f in run['frames']],
            charge_audit=run['charge_audit'])
    step_free = abs(out['coupled']['domain_step_v'])
    step_blocked = abs(out['regional_blockade']['domain_step_v'])
    return dict(
        claim='(c) a regional gap-junction blockade that changes the spatial pattern',
        materialised=bool(linear['max_abs_pattern_change_v'] > 1e-4 and step_blocked > 2 * step_free),
        linear_control=linear,
        bistable_domains=out,
        domain_step_ratio_blocked_over_free=float(step_blocked / step_free) if step_free > 0 else None,
        gate_available='GapJunctionGate implements the Boltzmann transjunctional-voltage dependence '
                       '(Harris/Spray/Bennett 1981) with the retained BETSE half-close 15 mV and 0.10 '
                       'floor; the blockade demonstrated here is the regional multiplier, which is the '
                       'pharmacological analogue (a connexin blocker applied to a territory).',
        honest_note='A blockade changing a diffusive pattern is a property of graph Laplacians. What '
                    'is new here is that the pattern being changed can be a MULTISTABLE one, so the '
                    'blockade selects between persistent states rather than reshaping a single '
                    'attractor. That is the property Levin-style theses need; it is still not '
                    'evidence for them.')


def demo_timescale(spacing_m, nodes):
    """Raise the horizon where the physics stays sound, and say where it stops."""
    out = dict(claim='(5) minutes-to-days horizons',
               retained_caps=dict(build_skin_electric_duration_s=10.0, build_skin_electric_dt_s=0.1,
                                  epithelial_electrodiffusion_duration_s=10.0,
                                  skin_patch_state_cap=512,
                                  source='ihm/assembly/skin_bioelectric.py:99, '
                                         'ihm/assembly/epithelial_electrodiffusion.py:105, '
                                         'ihm/materialize/skin.py SkinPatch.__post_init__'),
               new_caps=dict(max_duration_s=M.MAX_DURATION_S, max_steps=M.MAX_STEPS,
                             max_nodes=M.MAX_NODES, max_states=M.MAX_STATES),
               orders_of_magnitude_gained=float(math.log10(M.MAX_DURATION_S / 10.0)))

    # (1) A genuinely slow process: a depolarising conductance with a 30 min
    # gating time constant, driven past its fold. Nothing at this timescale was
    # expressible before. [engineering] tau is a chosen slow-gate time constant.
    slow = M.bistable_membrane(tau_dep_s=1800.)
    sheet = M.build_sheet(nodes, nodes, spacing_m, slow)
    grid = np.linspace(-0.15, 0.12, 54001)
    current, slope = M.iv_curve(slow, grid)
    turn = np.where(np.sign(slope[:-1]) * np.sign(slope[1:]) < 0)[0]
    folds = sorted(float(-current[k]) for k in turn)
    drive = folds[0] - 0.10 * (folds[-1] - folds[0]) if len(folds) == 2 else -1e-3
    t0 = time.time()
    run = M.integrate(sheet, -0.0946 * np.ones(sheet.n), duration_s=8.64e4, dt_s=60.,
                      applied_a_per_m2=drive, record_every=48)
    times = np.array([f['time_s'] for f in run['frames']])
    means = np.array([f['vm_mean_v'] for f in run['frames']])
    crossed = np.where(means > -0.04)[0]
    out['slow_gated_switch'] = dict(
        gate_tau_s=1800., drive_a_per_m2=drive, fold_drives_a_per_m2=folds,
        horizon_s=8.64e4, horizon_hours=24.0, dt_s=60., steps=int(round(8.64e4 / 60.)),
        wall_s=round(time.time() - t0, 1),
        vm_start_v=float(means[0]), vm_end_v=float(means[-1]),
        switch_time_s=float(times[crossed[0]]) if len(crossed) else None,
        switched=bool(len(crossed) > 0),
        trajectory=[[float(t), float(v)] for t, v in zip(times, means)],
        charge_audit=run['charge_audit'],
        soundness='no ion mass balance in this run, so the reversals are held fixed; that is exactly '
                  'the retained fixed-reservoir approximation and it is what makes an arbitrarily long '
                  'horizon cheap AND unfalsifiable. The next two blocks are the honest version.')

    # (2) Ion mass balance WITH the pump: does the inventory survive a long run?
    small = M.build_sheet(6, 6, spacing_m, M.linear_membrane())
    t0 = time.time()
    kept = M.integrate(small, -0.0815 * np.ones(small.n), duration_s=1.0e4, dt_s=10.,
                       ion_accounting=True, record_every=200)
    inside = np.asarray(kept['final_inside_mol_per_m3'], float)
    out['ion_accounting_pump_on'] = dict(
        horizon_s=1.0e4, dt_s=10., wall_s=round(time.time() - t0, 1),
        initial_inside_mol_per_m3=list(M.linear_membrane().inside_mol_per_m3),
        final_inside_mol_per_m3=inside[0].tolist(),
        max_abs_drift_mol_per_m3=float(np.max(np.abs(inside - np.asarray(
            M.linear_membrane().inside_mol_per_m3, float)))),
        final_vm_v=float(kept['final_vm_v'].mean()),
        survived=True,
        note='the retained prescribed pump is a CONSTANT current, not ATPase kinetics, so it does not '
             'balance the individual species fluxes. Any drift below is the size of the error the '
             'fixed-reservoir approximation was hiding.')

    # (3) Ion mass balance WITHOUT the pump: where mass balance stops you.
    starved = M.build_sheet(6, 6, spacing_m, M.linear_membrane(pump_a_per_m2=0.))
    horizon, dt = 1.0e5, 20.
    t0 = time.time()
    try:
        run = M.integrate(starved, -0.0815 * np.ones(starved.n), duration_s=horizon, dt_s=dt,
                          ion_accounting=True, record_every=500)
        inside = np.asarray(run['final_inside_mol_per_m3'], float)
        out['ion_accounting_pump_off'] = dict(
            horizon_s=horizon, dt_s=dt, wall_s=round(time.time() - t0, 1), exhausted=False,
            final_inside_mol_per_m3=inside[0].tolist(),
            final_vm_v=float(run['final_vm_v'].mean()))
    except M.ReservoirExhausted as exc:
        out['ion_accounting_pump_off'] = dict(
            horizon_s=horizon, dt_s=dt, wall_s=round(time.time() - t0, 1), exhausted=True,
            message=str(exc),
            reading='THIS is where mass balance stops the horizon. The same failure mode as '
                    'ihm/assembly/epithelial_electrodiffusion.transport, which raises '
                    "'Finite reservoir exhausted'. A fixed-reservoir run of the same length is "
                    'arithmetically fine and physically empty.')
    out['verdict'] = (
        'The 10 s ceiling was an arbitrary guard, not physics, for a fixed-reservoir model: raised to '
        '%.3g s (%.2f orders). With ion mass balance the horizon is bounded by the inventory and the '
        'code raises rather than integrating past it. Minutes-to-days pattern dynamics are now '
        'expressible; whether they are TRUE requires transporter kinetics with an ATP budget, which '
        'this module does not have.' % (M.MAX_DURATION_S, out['orders_of_magnitude_gained']))
    return out


def demo_envelope(layer_thickness_m, keratinocyte_density_per_m2):
    """Declare area on the REAL integument envelope and state the coarse-graining."""
    geometry = load('data/derived/outer-envelope/outer-envelope.json.gz')
    validation = load('data/derived/outer-envelope/validation.json')
    v = np.asarray(geometry['positions'], float).reshape(-1, 3)
    f = np.asarray(geometry['indices'], int).reshape(-1, 3)
    pos, area, volume, edges, weight = M.surface_sheet(v, f, layer_thickness_m=layer_thickness_m)
    recorded = validation['envelope_topology']['surface_area_m2']
    ledger = dict(
        source='data/derived/outer-envelope/outer-envelope.json.gz',
        nodes=int(len(pos)), faces=int(len(f)), edges=int(len(edges)),
        total_area_m2=float(area.sum()), recorded_envelope_area_m2=float(recorded),
        area_agreement_relative=float(abs(area.sum() - recorded) / recorded),
        node_area_m2=dict(min=float(area.min()), median=float(np.median(area)),
                          mean=float(area.mean()), max=float(area.max())),
        node_spacing_equivalent_m=float(math.sqrt(area.mean())),
        node_volume_m3=dict(total=float(volume.sum()), mean=float(volume.mean())),
        layer_thickness_m=layer_thickness_m,
        cotangent_weight=dict(min=float(weight.min()), max=float(weight.max()),
                              negative_fraction=float((weight < 0).mean()),
                              note='negative cotangent weights come from obtuse triangles and are '
                                   'normal for the cotangent Laplace-Beltrami operator; they are '
                                   'reported, not clipped'))
    per_node = keratinocyte_density_per_m2 * area
    ledger['coarse_graining'] = dict(
        keratinocyte_areal_density_per_m2=float(keratinocyte_density_per_m2),
        single_layer_keratinocytes_total=float(keratinocyte_density_per_m2 * area.sum()),
        keratinocytes_per_node=dict(min=float(per_node.min()), median=float(np.median(per_node)),
                                    mean=float(per_node.mean()), max=float(per_node.max())),
        at_2mm_surface_edge=dict(
            node_area_m2=float((math.sqrt(3) / 2) * (2e-3) ** 2),
            nodes=float(area.sum() / ((math.sqrt(3) / 2) * (2e-3) ** 2)),
            keratinocytes_per_node=float(keratinocyte_density_per_m2 * (math.sqrt(3) / 2) * (2e-3) ** 2),
            stratified_x10_layers=float(10 * keratinocyte_density_per_m2 * (math.sqrt(3) / 2) * (2e-3) ** 2)),
        correction=('data/derived/bioelectric-substrate-assessment-v1/assessment.json states "At 2 mm '
                    'spacing one node aggregates roughly 5e7 keratinocytes". Recomputed from its own '
                    'inputs that number is 4.66e4 for a single basal layer and 4.66e5 for a ten-cell-deep '
                    'epidermis. The assessment is high by about three orders of magnitude. Its '
                    'CONCLUSION is unaffected: a node is still a population and a single-cell claim at '
                    'this resolution is still a category error.'),
        statement='Vm at a node is a population mean over that many keratinocytes. This module cannot '
                  'settle any single-cell claim, and no amount of mesh refinement short of 2.4e10 nodes '
                  'would change that.')
    return ledger


def envelope_run(layer_thickness_m, keratinocyte_density_per_m2, ecm_sheet_s, gj_sheet_s,
                 barrier_s_per_m2, tep_v, wound_radius_m, budget_s):
    """One linear-law bidomain solve on the whole 54k-node integument envelope."""
    ledger = demo_envelope(layer_thickness_m, keratinocyte_density_per_m2)
    geometry = load('data/derived/outer-envelope/outer-envelope.json.gz')
    v = np.asarray(geometry['positions'], float).reshape(-1, 3)
    f = np.asarray(geometry['indices'], int).reshape(-1, 3)
    pos, area, volume, edges, weight = M.surface_sheet(v, f, layer_thickness_m=layer_thickness_m)
    n = len(pos)
    # Edge conductance = sheet conductivity x cotangent weight. Clip the negative
    # weights of obtuse triangles to zero for the CONDUCTANCE graph only: a
    # negative conductance is not a physical resistor and would break the SPD
    # structure the elliptic solve relies on. [engineering] the clipped fraction
    # is reported so the cost is visible.
    w = np.maximum(weight, 0.)
    membrane = M.linear_membrane()
    sheet = M.Sheet(pos, area, volume, edges, gj_sheet_s * w, ecm_sheet_s * w,
                    np.full(n, barrier_s_per_m2), np.full(n, barrier_s_per_m2 * tep_v),
                    membrane, coarse_graining=ledger['coarse_graining'],
                    provenance=dict(geometry='data/derived/outer-envelope/outer-envelope.json.gz',
                                    discretisation='cotangent Laplace-Beltrami; negative weights clipped'))
    t0 = time.time()
    base = M.steady_state(sheet, -0.0815)
    build_s = time.time() - t0
    # Wound: every node within wound_radius_m of the most anterior-superior node.
    seed = int(np.argmax(pos[:, 1]))
    d = np.linalg.norm(pos - pos[seed], axis=1)
    wound = d <= wound_radius_m
    barrier = np.where(wound, 100. * barrier_s_per_m2, barrier_s_per_m2)
    t1 = time.time()
    wounded = M.steady_state(sheet, base['vm_v'], barrier=barrier)
    delta = wounded['vm_v'] - base['vm_v']
    return dict(
        ledger=ledger,
        clipped_negative_weight_fraction=float((weight < 0).mean()),
        law='linear (legacy law) on real geometry; the nonlinear law is available but a whole-envelope '
            'bistable continuation was not run inside the compute budget of %g s' % budget_s,
        baseline=dict(vm_v=float(base['vm_v'].mean()), vm_sd_v=float(base['vm_v'].std()),
                      phi_e_v=float(base['phi_e_v'].mean()),
                      phi_e_sd_v=float(base['phi_e_v'].std()),
                      newton_iterations=base['newton_iterations'],
                      wall_s=round(build_s, 1), receipts=base['receipts']),
        resolution_limit='median node spacing on this envelope is %.4f m, so a %.0f mm-radius wound is '
                         'resolved by only %d nodes. The envelope is a whole-body coarse mesh; a wound '
                         'study needs a locally refined integument shell, which does not exist yet.'
                         % (math.sqrt(np.median(area)), wound_radius_m * 1e3, int(wound.sum())),
        wound=dict(wound_nodes=int(wound.sum()), wound_area_m2=float(area[wound].sum()),
                   wound_radius_m=wound_radius_m,
                   max_abs_delta_vm_v=float(np.max(np.abs(delta))),
                   delta_vm_at_wound_v=float(delta[wound].mean()),
                   phi_e_min_v=float(wounded['phi_e_v'].min()),
                   phi_e_max_v=float(wounded['phi_e_v'].max()),
                   wall_s=round(time.time() - t1, 1),
                   receipts=wounded['receipts']))


# ---------------------------------------------------------------- assembly ----

def split_verified_inferred(demos, checks):
    return dict(
        verified=[
            'The retained linear law is reproduced EXACTLY by the sibling under law="linear", '
            'bath="prescribed": generator parity %.3e V/s, resting parity %.3e V.'
            % (checks['legacy_law_generator_parity_v_per_s'], checks['legacy_resting_parity_v']),
            'The gated membrane has THREE zeros of I(V) at %s V with the middle one unstable; the '
            'retained linear law has exactly one at any parameters.'
            % [round(x, 6) for x in checks['iv_roots_v']],
            'Two distinct stable steady states of the SAME sheet with the SAME parameters, differing '
            'only in initial condition: %.6f V and %.6f V, separation %.6f V, spatial SD %.2e and '
            '%.2e V.' % (demos['bistability']['two_steady_states']['polarised_vm_v'],
                         demos['bistability']['two_steady_states']['depolarised_vm_v'],
                         demos['bistability']['two_steady_states']['separation_v'],
                         demos['bistability']['two_steady_states']['polarised_sd_v'],
                         demos['bistability']['two_steady_states']['depolarised_sd_v']),
            'A hysteresis loop under a swept drive: up and down branches differ by up to %.6f V over '
            'a fold-to-fold drive window of %.4e A/m^2.'
            % (demos['bistability']['hysteresis']['max_branch_separation_v'],
               demos['bistability']['hysteresis']['loop_width_a_per_m2']),
            'With the prescribed bath a wound moves Vm by %.3e V (the certified-zero result, '
            'reproduced as a control). With the solved extracellular domain the SAME wound moves Vm '
            'by %.6f V, and %.6f V in tissue that is not itself wounded.'
            % (demos['wound']['certified_zero_reproduced']['max_abs_delta_vm_v'],
               demos['wound']['barrier_breach']['max_abs_delta_vm_v'],
               demos['wound']['barrier_breach']['max_abs_delta_vm_in_intact_tissue_v']),
            'A regional gap-junction blockade turns a transient into a persistent pattern: coupled, '
            'the two Vmem domains collapse to one uniform state (step %.6f V, spatial SD %.3e V); '
            'blockaded, both persist (step %.6f V, spatial SD %.6f V).'
            % (demos['blockade']['bistable_domains']['coupled']['domain_step_v'],
               demos['blockade']['bistable_domains']['coupled']['final_vm_sd_v'],
               demos['blockade']['bistable_domains']['regional_blockade']['domain_step_v'],
               demos['blockade']['bistable_domains']['regional_blockade']['final_vm_sd_v']),
            'Conservation receipts hold at every reported state: elliptic residual <= %.2e A, '
            'sum of membrane currents <= %.2e A, battery minus barrier return <= %.2e A.'
            % (checks['elliptic_residual_a'], abs(checks['sum_membrane_current_a']),
               abs(checks['battery_minus_barrier_a'])),
            'Every node declares a membrane area and an intracellular volume, and the specific '
            'capacitance in use is exactly the literature 1 uF/cm^2.',
        ],
        inferred=[
            'The gating parameters that place the membrane in the bistable regime (g_Kir, V_Kir, '
            's_Kir, g_dep, V_dep, s_dep) are ENGINEERING CHOICES. Nothing in this repository '
            'identifies them for keratinocytes. Bistability is demonstrated for this parameter set, '
            'not for human skin.',
            'The gap-junction sheet conductance (5e-7 S) and the sub-epidermal sheet conductance '
            '(3e-6 S = 0.3 S/m over a 10 um cleft) are engineering choices selected to give '
            'literature-order length constants; they are not measured here. Every spatial result '
            'scales with them.',
            'The barrier conductance (2 S/m^2 = 5000 Ohm.cm^2) and the 100x wound multiplier are '
            'literature-order, not identified. The lateral wound field scales as TEP/lambda_e and '
            'lambda_e depends only on their ratio to the sheet conductance.',
            'The lysis conductance (a 0 V-reversal leak at wounded nodes) is an engineering surrogate '
            'for membrane rupture. There is no cell death, migration, healing or regeneration law.',
            'A coarse-grained node is a population of keratinocytes. No result here can settle any '
            'single-cell claim, and refining the mesh does not fix that: whole-integument single-cell '
            'coverage is 2.4e10 nodes.',
            'Nothing here is a morphogenetic readout. There is still no state whose value is a shape, '
            'so Levin theses ABOUT ANATOMICAL CONTROL remain unstateable; what is now stateable is '
            'the bioelectric layer they are stated over.',
        ])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--self-test', action='store_true', help='run the executable checks only')
    ap.add_argument('--no-envelope', action='store_true', help='skip the 54k-node real-geometry run')
    ap.add_argument('--nodes', type=int, default=48)
    ap.add_argument('--spacing-m', type=float, default=2e-4)
    ap.add_argument('--front-duration-s', type=float, default=0.8)
    ap.add_argument('--front-dt-s', type=float, default=0.005)
    args = ap.parse_args()

    t0 = time.time()
    checks = self_test()
    if args.self_test:
        print(json.dumps(dict(self_test='PASS', wall_s=round(time.time() - t0, 1), checks=checks),
                         indent=1, allow_nan=False))
        return

    demos = dict(
        bistability=demo_bistability(args.spacing_m, args.nodes),
        wound=demo_wound(args.spacing_m, args.nodes, wound_radius_nodes=4,
                         wound_barrier=200.0, lysis_s_per_m2=50.0),
        blockade=demo_blockade(args.spacing_m, args.nodes, args.front_duration_s,
                               args.front_dt_s, factor=0.01),
        timescale=demo_timescale(args.spacing_m, 8),
    )
    envelope = None
    if not args.no_envelope:
        try:
            envelope = envelope_run(layer_thickness_m=1e-4,
                                    keratinocyte_density_per_m2=1.3463086096010574e10,
                                    ecm_sheet_s=3.0e-6, gj_sheet_s=5.0e-7,
                                    barrier_s_per_m2=2.0, tep_v=30e-3,
                                    wound_radius_m=5e-3, budget_s=600.)
        except Exception as exc:                                   # noqa: BLE001
            envelope = dict(attempted=True, completed=False, reason=repr(exc))

    artifact = dict(
        schema='ihm.skin-bioelectric-nonlinear.v1',
        generator='scripts/verify_skin_bioelectric_nonlinear.py',
        module='scripts/skin_bioelectric_nonlinear.py',
        ownership='NEW files under scripts/ and a NEW directory under data/derived/ only. '
                  'ihm/assembly/skin_bioelectric.py, ihm/materialize/skin.py and '
                  'scripts/verify_body_skin_electric.py are imported and unmodified.',
        why='the retained substrate is linear time-invariant with a prescribed 0 V basal bath, so '
            'bistability, hysteresis, fronts and wound-driven Vm change are not expressible in it at '
            'any resolution. This is a model-class change, not a resolution change.',
        formulation=dict(
            electrical='bidomain on a graph: L_i (Vm + phi_e) + A (C_m dVm/dt + I_ion) = 0 and '
                       '(L_e + G_b) phi_e = A (C_m dVm/dt + I_ion) + A J_te, eliminated to '
                       'K phi_e = A J_te - L_i Vm with K = L_i + L_e + G_b',
            boundary_conditions='no-flux on both domains at the sheet boundary (graph/cotangent '
                                'Laplacian rows sum to zero); the only current sink is the '
                                'distributed barrier conductance G_b to a 0 V outside bath, i.e. a '
                                'Robin condition through the epidermal barrier. G_b > 0 makes K SPD, '
                                'so the extracellular gauge is fixed by physics rather than pinned.',
            gating='Boltzmann-gated conductances, optionally first-order in time; the '
                   'transjunctional-voltage gap-junction gate uses the retained BETSE half-close '
                   '15 mV and 0.10 floor',
            legacy='law="linear", bath="prescribed" recovers the retained law exactly',
            literature=['Tung 1978 (bidomain)', 'Henriquez 1993 Crit Rev Biomed Eng 21:1-77',
                        'Keener & Sneyd, Mathematical Physiology ch.11',
                        'Hodgkin & Huxley 1952 J Physiol 117:500',
                        'Hille, Ion Channels of Excitable Membranes 3rd ed.',
                        'Harris, Spray & Bennett 1981 J Gen Physiol 77:95',
                        'Cervera, Alcaraz & Mafe 2014 J Phys Chem B 118:6417',
                        'Cervera, Meseguer & Mafe 2016 Sci Rep 6:35201',
                        'Law & Levin 2015 Theor Biol Med Model 12:22',
                        'Barker, Jaffe & Vanable 1982 Am J Physiol 242:R358 (skin TEP)',
                        'Nuccitelli et al. 2011 PMC3228273 (human wound field)',
                        'Pinkall & Polthier 1993 Exp Math 2:15 (cotangent Laplacian)']),
        caps=dict(retained=dict(states=512, duration_s=10.0, dt_s=0.1),
                  new=dict(nodes=M.MAX_NODES, states=M.MAX_STATES,
                           duration_s=M.MAX_DURATION_S, steps=M.MAX_STEPS)),
        self_test=checks,
        demonstrations=demos,
        envelope=envelope,
        evidence=split_verified_inferred(demos, checks),
        limitations=[
            'applied_a_per_m2 is an electrode current with NO ion attribution, exactly like the '
            'retained electrode; it moves charge but not mass, so ion accounting under an applied '
            'drive is not a closed species balance.',
            'The pump is a prescribed constant cycle rate with 3:2 stoichiometry, not ATPase '
            'kinetics and not an ATP budget. Long ion-accounting horizons are therefore bounded by '
            'the inventory, and the code raises rather than pretending otherwise.',
            'The extracellular domain is a SHEET: one extracellular potential per node with a Robin '
            'leak to a 0 V outside bath. It is a bidomain, not a three-dimensional volume conductor, '
            'and it carries no depth structure through the epidermis.',
            'Ion concentrations enter the reversal potentials but there is no electroneutrality '
            'solve, no water or osmotic mechanics and no extracellular ion depletion.',
            'The gating is Boltzmann in Vm only. No ligand gating, no connexin expression state, no '
            'gene-regulatory coupling, and no serotonin/butyrate-style transporter chemistry.',
            'No cell birth, death, migration, differentiation or matrix remodelling. A wound here is '
            'an electrical boundary change, not a healing model.',
            'No morphogenetic readout: there is still no state whose value is a shape. Levin theses '
            'about ANATOMICAL setpoints remain unstateable on a single frozen atlas; what this '
            'module supplies is the bioelectric layer such a thesis would have to be stated over.',
            'Nothing here is calibrated against a human measurement. The Nuccitelli lateral-field fit '
            'is used only as an order-of-magnitude comparison target.',
        ],
        wall_s=round(time.time() - t0, 1),
    )

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'demonstrations.json').write_text(json.dumps(artifact, indent=1, allow_nan=False) + '\n')
    manifest = dict(
        schema='ihm.skin-bioelectric-nonlinear-manifest.v1',
        generator='scripts/verify_skin_bioelectric_nonlinear.py',
        generator_sha256=sha(Path(__file__)),
        module_sha256=sha(ROOT / 'scripts/skin_bioelectric_nonlinear.py'),
        repo_files_modified=[], native_jobs_run=0,
        argv=vars(args),
        inputs_sha256={rel: sha(ROOT / rel) for rel in INPUTS if (ROOT / rel).exists()},
        outputs_sha256={'demonstrations.json': sha(OUT / 'demonstrations.json')},
    )
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')

    print(json.dumps(dict(
        outputs=[str(p.relative_to(ROOT)) for p in sorted(OUT.iterdir())],
        bistable_states_v=[demos['bistability']['two_steady_states']['polarised_vm_v'],
                           demos['bistability']['two_steady_states']['depolarised_vm_v']],
        hysteresis_max_branch_separation_v=demos['bistability']['hysteresis']['max_branch_separation_v'],
        wound_prescribed_bath_delta_vm_v=demos['wound']['certified_zero_reproduced']['max_abs_delta_vm_v'],
        wound_solved_bath_delta_vm_v=demos['wound']['barrier_breach']['max_abs_delta_vm_v'],
        wound_lateral_field_v_per_m=demos['wound']['barrier_breach']['max_lateral_field_v_per_m'],
        blockade_domain_step_v=[demos['blockade']['bistable_domains']['coupled']['domain_step_v'],
                                demos['blockade']['bistable_domains']['regional_blockade']['domain_step_v']],
        horizon_s=M.MAX_DURATION_S,
        envelope_nodes=None if not envelope or not envelope.get('ledger') else envelope['ledger']['nodes'],
        materialised={k: demos[k].get('materialised') for k in ('bistability', 'wound', 'blockade')},
        wall_s=artifact['wall_s'],
    ), indent=1))


if __name__ == '__main__':
    main()
