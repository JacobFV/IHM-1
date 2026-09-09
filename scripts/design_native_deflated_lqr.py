#!/usr/bin/env python3
"""Design muscle excitation feedback after deflating the ground-plane symmetry.

At a static equilibrium on flat ground the native tangent map has an exact
three dimensional invariant subspace with unit discrete eigenvalue that no
muscle can reach: translating or yawing the whole body in the ground plane
changes neither gravity nor the contact geometry. Penalising those directions
in the LQR stage cost asks for a stabilising Riccati solution that provably
does not exist, which is what makes solve_discrete_are return a large residual
on the loaded-landing linearisations.

This designer measures that subspace from the retained A/B, removes it from the
design coordinates and the stage cost, solves the reduced algebraic Riccati
equation and lifts the gain back. The resulting feedback ignores absolute
ground-plane position and yaw by construction, and the closed loop keeps those
three unit eigenvalues untouched. No physical model parameter is changed; only
the feedback gain K in the artifact differs from the source design.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
from scipy.linalg import expm, solve_discrete_are

ROOT = Path(__file__).resolve().parents[1]


def deflation_basis(Ad, Bd, *, tolerance, input_tolerance, gap, maximum):
    """Right/left bases of the unit-eigenvalue subspace unreachable from Bd.

    The subspace is taken from the smallest singular directions of Ad - I, and
    is accepted only when it is separated from the rest of the spectrum by a
    clear gap and every direction in it is numerically unreachable from Bd.
    """
    U, s, Vt = np.linalg.svd(Ad - np.eye(Ad.shape[0]))
    order = np.argsort(s)
    relative = s[order] / s[order[-1]]
    reach_all = np.linalg.norm(U[:, order].T @ Bd, axis=1) / np.linalg.norm(U[:, order].T, axis=1)
    k = 0
    while k < min(maximum, len(s) - 1) and relative[k] <= tolerance and reach_all[k] <= input_tolerance:
        k += 1
    if not k:
        raise ValueError('No unreachable unit-eigenvalue subspace found; nothing to deflate')
    if relative[k] < gap * max(relative[k - 1], np.finfo(float).eps):
        raise ValueError('Unit-eigenvalue subspace is not separated from the rest of the spectrum')
    keep = order[:k]
    N = Vt[keep].T.copy()           # (A-I) N ~ 0
    W = U[:, keep].T.copy()         # W (A-I) ~ 0
    return N, W, relative[:k], reach_all[:k], float(relative[k]), float(reach_all[k])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', help='linearization.npz carrying A, B, Q, R, u0, minimum_activation')
    p.add_argument('--output-name', default='deflated_margin')
    p.add_argument('--dt', type=float, default=.01)
    p.add_argument('--singular-tolerance', type=float, default=1e-8)
    p.add_argument('--input-tolerance', type=float, default=1e-5)
    p.add_argument('--spectral-gap', type=float, default=100.)
    p.add_argument('--maximum-deflation', type=int, default=8)
    a = p.parse_args()
    if not np.isfinite(a.dt) or a.dt <= 0:
        raise ValueError('Positive finite sampling interval required')

    source = Path(a.source)
    with np.load(source, allow_pickle=False) as loaded:
        data = {k: loaded[k] for k in loaded.files}
    A, B = data['A'], data['B']
    n, m = B.shape
    if 'minimum_activation' not in data:
        raise ValueError('Margin-aware design requires native minimum activation values')
    transition = expm(np.block([[A, B], [np.zeros((m, n + m))]]) * a.dt)
    Ad, Bd = transition[:n, :n], transition[:n, n:]
    Q = data['Q'] * a.dt
    R = np.diag(np.diag(data['R'] * a.dt) * (.05 / np.maximum(data['u0'] - data['minimum_activation'], .005)) ** 2)

    N, W, singular, reach, next_singular, next_reach = deflation_basis(
        Ad, Bd, tolerance=a.singular_tolerance, input_tolerance=a.input_tolerance,
        gap=a.spectral_gap, maximum=a.maximum_deflation)
    k = N.shape[1]
    # Complete N to a full basis and read off the complementary coordinates.
    Qr, _ = np.linalg.qr(np.hstack([N, np.eye(n)]))
    T = Qr[:, :n]
    T[:, :k] = N
    Tinv = np.linalg.inv(T)
    S, C = Tinv[k:, :], T[:, k:]
    coupling = float(np.abs(S @ Ad @ N).max())
    A22, B2, Q22 = S @ Ad @ C, S @ Bd, C.T @ Q @ C
    P = solve_discrete_are(A22, B2, Q22, R)
    K2 = np.linalg.solve(R + B2.T @ P @ B2, B2.T @ P @ A22)
    K = K2 @ S
    residual = float(np.linalg.norm(A22.T @ P @ A22 - P - A22.T @ P @ B2 @ K2 + Q22) / np.linalg.norm(Q22))
    reduced = np.abs(np.linalg.eigvals(A22 - B2 @ K2))
    full = np.abs(np.linalg.eigvals(Ad - Bd @ K))

    names = [str(x) for x in data['state_names']]
    weight = np.abs(N).sum(axis=1)
    dominant = [names[i] for i in np.argsort(-weight)[:8]]

    report_path = source.with_name('report.json')
    report = json.loads(report_path.read_text())
    registration = json.loads((ROOT / report['registration']).read_text())
    retained_model = hashlib.sha256((source.parent / 'inputs/subject_walk_scaled.osim').read_bytes()).hexdigest()
    if retained_model != registration['model_sha256']:
        raise ValueError('Retained linearization model differs from current registration')
    report.update(
        schema='ihm.native-stance-deflated-lqr.v1',
        design='ground-plane-symmetry-deflated discrete LQR; margin-aware stage cost',
        dt_s=a.dt, margin_aware=True, model_sha256=retained_model,
        source_artifact_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        deflated_dimension=int(k),
        deflated_relative_singular_values=[float(x) for x in singular],
        deflated_input_reachability=[float(x) for x in reach],
        first_retained_relative_singular_value=next_singular,
        first_retained_input_reachability=next_reach,
        deflated_dominant_states=dominant,
        deflated_block_coupling_max=coupling,
        reduced_dare_residual_relative=residual,
        reduced_closed_loop_spectral_radius=float(reduced.max()),
        reduced_closed_loop_second_radius=float(np.sort(reduced)[-2]),
        full_closed_loop_spectral_radius=float(full.max()),
        full_closed_loop_modes_outside_unit_disc=int((full > 1 + 1e-9).sum()),
        full_closed_loop_radius_excluding_deflated=float(np.sort(full)[-(k + 1)]),
        discrete_gain_max=float(abs(K).max()),
        lqr_solved=True,
    )
    out = source.parent / a.output_name
    out.mkdir(exist_ok=False)
    data.update(K=K, Ad=Ad, Bd=Bd, Q_discrete=Q, R_discrete=R, dt_s=np.array(a.dt),
                model_sha256=np.array(retained_model), deflated_dimension=np.array(k))
    np.savez_compressed(out / 'linearization.npz', **data)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'output': str(out), **{k2: v for k2, v in report.items()
                                             if k2.startswith(('deflat', 'reduced', 'full', 'discrete_gain'))}}, indent=2))


if __name__ == '__main__':
    main()
