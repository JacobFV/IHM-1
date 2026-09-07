#!/usr/bin/env python3
"""Nonlinear bidomain epithelial bioelectric sheet.

A SIBLING to ihm/assembly/skin_bioelectric.py + ihm/materialize/skin.py. It imports
them and mutates nothing. It exists because the retained substrate is linear
time-invariant with a prescribed 0 V basal bath, so bistability, hysteresis,
fronts and wound-driven Vm change are not expressible in it at any resolution.

WHAT IS DIFFERENT
  1. Boltzmann-gated conductances (voltage-dependent), so the current-voltage
     relation can be N-shaped and more than one steady state can exist.
  2. A SOLVED extracellular potential (bidomain) instead of a prescribed bath,
     so a barrier breach shunts real current and moves cell Vm.
  3. Every node declares a membrane area and an intracellular volume, so every
     parameter is a specific quantity comparable with literature.
  4. The retained linear law is reachable unchanged as an explicit selection
     (law='linear', bath='prescribed'), the way ihm/brain keeps its legacy
     regional law behind source_pin=None.

FORMULATION (bidomain on a graph; Tung 1978, Henriquez 1993 Crit Rev Biomed Eng
21:1-77, Keener & Sneyd Mathematical Physiology ch.11) [literature]

  nodes i with membrane area A_i; L_i = graph Laplacian of gap-junction edge
  conductances (S); L_e = graph Laplacian of interstitial edge conductances (S);
  G_b = diag(A_i g_b,i) the barrier leak from the sub-epidermal extracellular
  space to the outside bath held at 0 V; J_te = active transepithelial current
  density driven into the extracellular space (the epidermal battery).

  Vm = phi_i - phi_e
  intracellular node balance   L_i (Vm + phi_e) + A (C_m dVm/dt + I_ion) = 0
  extracellular node balance   (L_e + G_b) phi_e = A (C_m dVm/dt + I_ion) + A J_te
  eliminating the membrane current
                               K phi_e = A J_te - L_i Vm,   K = L_i + L_e + G_b

BOUNDARY CONDITIONS. Both domains are no-flux at the sheet boundary (natural for
a graph/cotangent Laplacian: rows sum to zero). The ONLY current sink is the
distributed barrier conductance G_b to the outside bath at 0 V, i.e. a Robin
condition applied through the epidermal barrier rather than a Dirichlet patch.
G_b > 0 everywhere is what makes K symmetric positive definite, so the
extracellular gauge is fixed by physics and not by an arbitrary pin. A wound is
g_b -> g_b_wound at the wounded nodes. With bath='prescribed', phi_e == 0 is
imposed instead and the retained model is recovered exactly.

EXACT CONSERVATION IDENTITIES USED AS RECEIPTS
  sum_i A_i I_m,i          == 0                  (1^T L_i == 0)
  sum_i A_i g_b,i phi_e,i  == sum_i A_i J_te,i   (the battery returns through the barrier)
  A_i C_m dVm_i            == -integral (L_i phi_i)_i dt - A_i integral I_ion,i dt

EVIDENCE TAGS  [literature] cited form or value; [retained] recomputed from a
hashed repository input; [engineering] a choice made here and defensible only as
a choice.
"""
from dataclasses import dataclass, replace
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

R_GAS = 8.314462618
FARADAY = 96485.33212

# [literature] Hodgkin & Huxley 1952 J Physiol 117:500; Cole 1968; Gentet, Stuart
# & Clements 2000 Biophys J 79:314 measure 0.9 uF/cm^2. 1 uF/cm^2 = 1e-2 F/m^2.
SPECIFIC_CAPACITANCE_F_PER_M2 = 1.0e-2

# [engineering] hard caps. Raised from the retained 512 states / 10 s ceiling.
MAX_NODES = 200_000
MAX_STATES = 1_500_000
MAX_DURATION_S = 8.64e5          # 10 days; Levin-style pattern work is 1e2-1e5 s
MAX_STEPS = 400_000


class ReservoirExhausted(RuntimeError):
    """Ion mass balance ran an intracellular species to the floor."""


@dataclass(frozen=True)
class Gate:
    """Boltzmann-gated open fraction, optionally first-order in time.

    m_inf(V) = 1 / (1 + exp((V - v_half_v) / slope_v)) ** power
    slope_v < 0 opens on depolarisation; slope_v > 0 closes on depolarisation
    (inward-rectifier sense). tau_s = 0 is instantaneous gating.

    [literature] Boltzmann conductance gating: Hodgkin & Huxley 1952; Hille,
    Ion Channels of Excitable Membranes 3rd ed. ch.2-3. Two-channel Boltzmann
    membranes producing bistable Vmem and multicellular Vmem domains:
    Cervera, Alcaraz & Mafe 2014 J Phys Chem B 118:6417; Cervera, Meseguer &
    Mafe 2016 Sci Rep 6:35201; Law & Levin 2015 Theor Biol Med Model 12:22.
    """
    v_half_v: float
    slope_v: float
    tau_s: float = 0.0
    power: int = 1

    def __post_init__(self):
        if not all(math.isfinite(v) for v in (self.v_half_v, self.slope_v, self.tau_s)):
            raise ValueError('finite gate parameters required')
        if self.slope_v == 0 or self.tau_s < 0 or not isinstance(self.power, int) or self.power < 1:
            raise ValueError('nonzero slope, nonnegative tau, positive integer power required')

    def steady(self, v):
        z = np.clip((np.asarray(v, float) - self.v_half_v) / self.slope_v, -500., 500.)
        s = 1. / (1. + np.exp(z))
        d = -s * (1. - s) / self.slope_v
        if self.power == 1:
            return s, d
        return s ** self.power, self.power * s ** (self.power - 1) * d

    def advance(self, previous, v, dt):
        """Backward-Euler gate update; returns (m, dm/dV) with the chain rule."""
        s, ds = self.steady(v)
        if self.tau_s == 0.:
            return s, ds
        w = (dt / self.tau_s) / (1. + dt / self.tau_s)
        return (previous + (dt / self.tau_s) * s) / (1. + dt / self.tau_s), w * ds


@dataclass(frozen=True)
class Channel:
    ident: str
    species: str                      # 'Na' | 'K' | 'Cl'
    g_max_s_per_m2: float
    gate: object = None               # None -> ohmic leak (the retained law)
    evidence: str = 'engineering'

    def __post_init__(self):
        if self.species not in ('Na', 'K', 'Cl'):
            raise ValueError('species must be Na, K or Cl')
        if not math.isfinite(self.g_max_s_per_m2) or self.g_max_s_per_m2 <= 0:
            raise ValueError('positive maximal specific conductance required')
        if self.gate is not None and not isinstance(self.gate, Gate):
            raise ValueError('gate must be a Gate or None')


VALENCE = {'Na': 1, 'K': 1, 'Cl': -1}


@dataclass(frozen=True)
class Membrane:
    """Per-unit-area membrane. Every quantity is specific; nothing is lumped."""
    channels: tuple
    pump_cycles_mol_per_m2_s: float          # Na/K-ATPase; 3 Na out, 2 K in, +1 e out
    specific_capacitance_f_per_m2: float = SPECIFIC_CAPACITANCE_F_PER_M2
    temperature_k: float = 310.15            # [literature] core temperature
    inside_mol_per_m3: tuple = (15., 140., 30.)      # Na, K, Cl  [retained] SkinPatch
    outside_mol_per_m3: tuple = (140., 4., 110.)     # Na, K, Cl  [retained] SkinPatch
    law: str = 'gated'                       # 'gated' | 'linear' (legacy selection)

    def __post_init__(self):
        if self.law not in ('gated', 'linear'):
            raise ValueError("law must be 'gated' or 'linear'")
        if not self.channels or any(not isinstance(c, Channel) for c in self.channels):
            raise ValueError('at least one Channel required')
        for v in (self.pump_cycles_mol_per_m2_s, self.specific_capacitance_f_per_m2, self.temperature_k):
            if not math.isfinite(v):
                raise ValueError('finite membrane scalars required')
        if self.specific_capacitance_f_per_m2 <= 0 or self.temperature_k <= 0:
            raise ValueError('positive capacitance and temperature required')
        for name in ('inside_mol_per_m3', 'outside_mol_per_m3'):
            c = getattr(self, name)
            if len(c) != 3 or any(not math.isfinite(x) or x <= 0 for x in c):
                raise ValueError('three positive concentrations required in ' + name)

    @property
    def thermal_v(self):
        return R_GAS * self.temperature_k / FARADAY

    def reversals(self, inside=None):
        """Nernst reversal per species from concentrations (mol/m^3)."""
        inside = np.asarray(self.inside_mol_per_m3 if inside is None else inside, float)
        outside = np.asarray(self.outside_mol_per_m3, float)
        z = np.array([VALENCE['Na'], VALENCE['K'], VALENCE['Cl']], float)
        return self.thermal_v / z * np.log(outside / np.maximum(inside, 1e-12))

    @property
    def pump_current_a_per_m2(self):
        return FARADAY * self.pump_cycles_mol_per_m2_s

    def current(self, v, gates=None, reversals=None):
        """Return (I_ion A/m^2, dI/dV S/m^2, per-species charge flux A/m^2).

        Positive current is outward. gates is a dict channel-ident -> open
        fraction, ignored when law == 'linear'.
        """
        v = np.asarray(v, float)
        e = self.reversals() if reversals is None else np.asarray(reversals, float)
        e_by = {'Na': e[..., 0], 'K': e[..., 1], 'Cl': e[..., 2]}
        total = np.zeros_like(v)
        slope = np.zeros_like(v)
        per = {'Na': np.zeros_like(v), 'K': np.zeros_like(v), 'Cl': np.zeros_like(v)}
        for c in self.channels:
            drive = v - e_by[c.species]
            if c.gate is None or self.law == 'linear':
                g = np.full_like(v, c.g_max_s_per_m2)
                dg = np.zeros_like(v)
            else:
                m = gates[c.ident] if gates is not None and c.ident in gates else c.gate.steady(v)[0]
                dm = gates.get(c.ident + '/dV') if gates is not None else None
                if dm is None:
                    dm = c.gate.steady(v)[1]
                g = c.g_max_s_per_m2 * np.asarray(m, float)
                dg = c.g_max_s_per_m2 * np.asarray(dm, float)
            total = total + g * drive
            slope = slope + g + dg * drive
            per[c.species] = per[c.species] + g * drive
        pump = self.pump_current_a_per_m2
        total = total + pump
        per['Na'] = per['Na'] + 3. * FARADAY * self.pump_cycles_mol_per_m2_s
        per['K'] = per['K'] - 2. * FARADAY * self.pump_cycles_mol_per_m2_s
        return total, slope, per


@dataclass(frozen=True)
class GapJunctionGate:
    """Transjunctional-voltage-gated gap junction conductance.

    g(Vj) = g_min + (g_max - g_min) / (1 + exp((|Vj| - v_half_v) / slope_v))

    [literature] Boltzmann Vj dependence of junctional conductance: Harris,
    Spray & Bennett 1981 J Gen Physiol 77:95; Spray, Harris & Bennett 1981
    Science 211:712. [retained] the half-closing voltage 15 mV and the 0.10
    residual open fraction are the values in the retained BETSE run
    (data/derived/physiology/betse_run/sim_config.yaml, gj_vthresh / gj_min).
    """
    v_half_v: float = 15e-3
    slope_v: float = 3e-3
    min_fraction: float = 0.10

    def __post_init__(self):
        if not all(math.isfinite(v) for v in (self.v_half_v, self.slope_v, self.min_fraction)):
            raise ValueError('finite gap-junction gate parameters required')
        if self.slope_v <= 0 or not 0. <= self.min_fraction <= 1.:
            raise ValueError('positive slope and open fraction in [0,1] required')

    def fraction(self, vj):
        z = np.clip((np.abs(np.asarray(vj, float)) - self.v_half_v) / self.slope_v, -500., 500.)
        return self.min_fraction + (1. - self.min_fraction) / (1. + np.exp(z))


@dataclass(frozen=True)
class Sheet:
    """A coarse-grained epithelial sheet with DECLARED per-node area and volume.

    node_area_m2      membrane area represented by the node (m^2)
    node_volume_m3    intracellular volume represented by the node (m^3)
    edges             (m,2) int node pairs
    gj_s              (m,) gap-junction edge conductance (S) at full opening
    ecm_s             (m,) interstitial edge conductance (S)
    barrier_s_per_m2  (n,) barrier conductance from extracellular space to the
                      0 V outside bath (S/m^2); a wound raises it locally
    transepithelial_a_per_m2 (n,) active current density into the extracellular
                      space (the epidermal battery)
    coarse_graining   explicit statement of how many cells one node aggregates
    """
    positions_m: object
    node_area_m2: object
    node_volume_m3: object
    edges: object
    gj_s: object
    ecm_s: object
    barrier_s_per_m2: object
    transepithelial_a_per_m2: object
    membrane: Membrane
    gj_gate: object = None
    bath: str = 'solved'                # 'solved' (bidomain) | 'prescribed' (legacy 0 V)
    coarse_graining: object = None
    provenance: object = None

    def __post_init__(self):
        n = len(self.positions_m)
        if not 1 <= n <= MAX_NODES:
            raise ValueError('node count must be in [1, %d]' % MAX_NODES)
        if self.bath not in ('solved', 'prescribed'):
            raise ValueError("bath must be 'solved' or 'prescribed'")
        for name, shape, positive in (('positions_m', (n, 3), False), ('node_area_m2', (n,), True),
                                      ('node_volume_m3', (n,), True), ('barrier_s_per_m2', (n,), True),
                                      ('transepithelial_a_per_m2', (n,), False)):
            a = np.asarray(getattr(self, name), float)
            if a.shape != shape or not np.isfinite(a).all() or (positive and (a <= 0).any()):
                raise ValueError('invalid ' + name)
            a = np.ascontiguousarray(a)
            a.setflags(write=False)
            object.__setattr__(self, name, a)
        e = np.asarray(self.edges, int)
        if e.ndim != 2 or e.shape[1] != 2 or (e < 0).any() or (e >= n).any() or (e[:, 0] == e[:, 1]).any():
            raise ValueError('invalid edge list')
        key = np.sort(e, axis=1)
        if len(np.unique(key[:, 0] * n + key[:, 1])) != len(e):
            raise ValueError('duplicate edges')
        e.setflags(write=False)
        object.__setattr__(self, 'edges', e)
        for name in ('gj_s', 'ecm_s'):
            a = np.asarray(getattr(self, name), float)
            if a.shape != (len(e),) or not np.isfinite(a).all() or (a < 0).any():
                raise ValueError('invalid ' + name)
            a.setflags(write=False)
            object.__setattr__(self, name, a)
        if self.gj_gate is not None and not isinstance(self.gj_gate, GapJunctionGate):
            raise ValueError('gj_gate must be a GapJunctionGate or None')
        gated = sum(1 for c in self.membrane.channels if c.gate is not None and c.gate.tau_s > 0
                    and self.membrane.law == 'gated')
        if n * (1 + gated) > MAX_STATES:
            raise ValueError('state count exceeds %d' % MAX_STATES)
        object.__setattr__(self, 'coarse_graining', dict(self.coarse_graining or {}))
        object.__setattr__(self, 'provenance', dict(self.provenance or {}))

    @property
    def n(self):
        return len(self.positions_m)

    @property
    def specific_capacitance_check(self):
        """Falsifiability receipt: the specific capacitance actually in use."""
        m = self.membrane
        return dict(specific_capacitance_f_per_m2=m.specific_capacitance_f_per_m2,
                    specific_capacitance_uf_per_cm2=m.specific_capacitance_f_per_m2 * 1e2,
                    literature_uf_per_cm2=1.0,
                    literature='Hodgkin & Huxley 1952; Gentet et al. 2000 Biophys J 79:314 (0.9 uF/cm2)',
                    lumped_capacitance_range_f=[float(self.node_area_m2.min() * m.specific_capacitance_f_per_m2),
                                                float(self.node_area_m2.max() * m.specific_capacitance_f_per_m2)],
                    total_membrane_area_m2=float(self.node_area_m2.sum()),
                    total_intracellular_volume_m3=float(self.node_volume_m3.sum()))


def laplacian(edges, weights, n):
    u, v = edges[:, 0], edges[:, 1]
    w = np.asarray(weights, float)
    rows = np.concatenate([u, v, u, v])
    cols = np.concatenate([v, u, u, v])
    data = np.concatenate([-w, -w, w, w])
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n))


def gj_weights(sheet, phi_i=None, blockade=None):
    w = sheet.gj_s if blockade is None else sheet.gj_s * np.asarray(blockade, float)
    if sheet.gj_gate is not None and phi_i is not None:
        vj = np.asarray(phi_i, float)[sheet.edges[:, 0]] - np.asarray(phi_i, float)[sheet.edges[:, 1]]
        w = w * sheet.gj_gate.fraction(vj)
    return w


def operators(sheet, phi_i=None, blockade=None, barrier=None):
    n = sheet.n
    li = laplacian(sheet.edges, gj_weights(sheet, phi_i, blockade), n)
    le = laplacian(sheet.edges, sheet.ecm_s, n)
    gb = sp.diags(sheet.node_area_m2 * (sheet.barrier_s_per_m2 if barrier is None else np.asarray(barrier, float)))
    return li, le, gb


def extracellular(sheet, vm, li=None, le=None, gb=None, blockade=None, barrier=None):
    """Solve K phi_e = A J_te - L_i Vm.  Returns (phi_e, elliptic residual)."""
    if sheet.bath == 'prescribed':
        return np.zeros(sheet.n), 0.0
    if li is None:
        li, le, gb = operators(sheet, None, blockade, barrier)
    k = (li + le + gb).tocsc()
    rhs = sheet.node_area_m2 * sheet.transepithelial_a_per_m2 - li @ vm
    phi = spla.spsolve(k, rhs)
    return phi, float(np.max(np.abs(k @ phi - rhs)))


def _gate_state(sheet, v, previous, dt):
    gates = {}
    if sheet.membrane.law != 'gated':
        return gates
    for c in sheet.membrane.channels:
        if c.gate is None:
            continue
        if c.gate.tau_s == 0. or previous is None or dt is None:
            m, dm = c.gate.steady(v)
        else:
            m, dm = c.gate.advance(previous[c.ident], v, dt)
        gates[c.ident] = m
        gates[c.ident + '/dV'] = dm
    return gates


def _residual_and_jacobian(sheet, v, phi_e, gates, li, le, gb, reversals, applied, leak, mass_term):
    """mass_term is (A*C_m/dt, A*C_m/dt*v_prev) or (0, 0) for the steady state."""
    m = sheet.membrane
    i_ion, di_dv, _ = m.current(v, gates, reversals)
    i_ion = i_ion + applied + leak * v
    di_dv = di_dv + leak
    a = sheet.node_area_m2
    scale, rhs_prev = mass_term
    f1 = scale * v - rhs_prev + li @ (v + phi_e) + a * i_ion
    if sheet.bath == 'prescribed':
        j11 = sp.diags(scale + a * di_dv) + li
        return f1, np.zeros(sheet.n), j11, None, None, None
    k = li + le + gb
    f2 = k @ phi_e + li @ v - a * sheet.transepithelial_a_per_m2
    j11 = sp.diags(scale + a * di_dv) + li
    return f1, f2, j11, li, li, k


def _solve_newton(sheet, v0, phi0, gates_prev, dt, reversals, applied, leak, mass_scale, mass_rhs,
                  blockade, barrier, tol=1e-12, max_iter=60):
    v = np.array(v0, float)
    phi = np.array(phi0, float)
    n = sheet.n
    a = sheet.node_area_m2
    # A fixed current scale in amperes: the current a 100 mV excursion drives
    # through one node's membrane, plus its capacitive current over one step.
    g_tot = sum(c.g_max_s_per_m2 for c in sheet.membrane.channels)
    scale = float(np.max(a)) * (g_tot * 0.1 + float(np.max(mass_scale)) * 0.1 / max(float(np.max(a)), 1e-300))
    scale = max(scale, 1e-30)
    last = float('inf')
    it = 0
    for it in range(max_iter):
        gates = _gate_state(sheet, v, gates_prev, dt)
        li, le, gb = operators(sheet, v + phi, blockade, barrier)
        f1, f2, j11, j12, j21, j22 = _residual_and_jacobian(
            sheet, v, phi, gates, li, le, gb, reversals, applied, leak, (mass_scale, mass_rhs))
        last = max(float(np.max(np.abs(f1))), float(np.max(np.abs(f2))) if sheet.bath == 'solved' else 0.)
        if last <= tol * scale:
            break
        if sheet.bath == 'prescribed':
            step_v = spla.spsolve(j11.tocsc(), -f1)
            step_p = np.zeros(n)
        else:
            block = sp.bmat([[j11, j12], [j21, j22]], format='csc')
            step = spla.spsolve(block, -np.concatenate([f1, f2]))
            step_v, step_p = step[:n], step[n:]
        if not np.isfinite(step_v).all() or not np.isfinite(step_p).all():
            raise RuntimeError('Newton step is not finite')
        # Backtracking line search on the infinity-norm residual: the N-shaped
        # I(V) makes an undamped Newton oscillate across the unstable branch.
        lam = 1.0
        for _ in range(30):
            vt, pt = v + lam * step_v, phi + lam * step_p
            gt = _gate_state(sheet, vt, gates_prev, dt)
            lt, let_, gbt = operators(sheet, vt + pt, blockade, barrier)
            g1, g2, *_ = _residual_and_jacobian(sheet, vt, pt, gt, lt, let_, gbt, reversals,
                                                applied, leak, (mass_scale, mass_rhs))
            trial = max(float(np.max(np.abs(g1))),
                        float(np.max(np.abs(g2))) if sheet.bath == 'solved' else 0.)
            if np.isfinite(trial) and trial < last:
                break
            lam *= 0.5
        v, phi = v + lam * step_v, phi + lam * step_p
        move = lam * max(float(np.max(np.abs(step_v))), float(np.max(np.abs(step_p))))
        if not np.isfinite(v).all() or not np.isfinite(phi).all():
            raise RuntimeError('Newton diverged')
        if move < 1e-14:
            break
    else:
        raise RuntimeError('Newton did not converge; residual %.3e against scale %.3e' % (last, scale))
    gates = _gate_state(sheet, v, gates_prev, dt)
    return v, phi, gates, last, it


def steady_state(sheet, v0, *, applied_a_per_m2=0., leak_s_per_m2=0., blockade=None, barrier=None,
                 reversals=None, tol=1e-12):
    """Solve the coupled nonlinear steady state (no capacitive term).

    leak_s_per_m2 is a non-selective membrane leak with reversal 0 V, the
    lysed-cell (injury) conductance: a ruptured membrane no longer separates
    the compartments, so Vm at that node is dragged to zero. [engineering]
    """
    v0 = np.broadcast_to(np.asarray(v0, float), (sheet.n,)).copy()
    applied = np.broadcast_to(np.asarray(applied_a_per_m2, float), (sheet.n,)).astype(float)
    leak = np.broadcast_to(np.asarray(leak_s_per_m2, float), (sheet.n,)).astype(float)
    rev = sheet.membrane.reversals() if reversals is None else np.asarray(reversals, float)
    phi0, _ = extracellular(sheet, v0, blockade=blockade, barrier=barrier)
    v, phi, gates, residual, iters = _solve_newton(
        sheet, v0, phi0, None, None, rev, applied, leak, np.zeros(sheet.n), np.zeros(sheet.n),
        blockade, barrier, tol=tol)
    return dict(vm_v=v, phi_e_v=phi, phi_i_v=v + phi, gates={k: g for k, g in gates.items() if '/dV' not in k},
                newton_residual_a=residual, newton_iterations=int(iters) + 1,
                receipts=_receipts(sheet, v, phi, gates, rev, applied, leak, blockade, barrier))


def _receipts(sheet, v, phi, gates, reversals, applied, leak, blockade, barrier, capacitive=None):
    """Conservation receipts. capacitive is A*C_m*dV/dt (zero at a steady state)."""
    a = sheet.node_area_m2
    li, le, gb = operators(sheet, v + phi, blockade, barrier)
    i_ion, _, per = sheet.membrane.current(v, gates, reversals)
    i_ion = i_ion + applied + leak * v
    membrane_total = -(li @ (v + phi))
    cap = np.zeros_like(v) if capacitive is None else np.asarray(capacitive, float)
    out = dict(
        sum_membrane_current_a=float(membrane_total.sum()),
        max_node_current_imbalance_a=float(np.max(np.abs(membrane_total - a * i_ion - cap))),
        battery_a=float((a * sheet.transepithelial_a_per_m2).sum()),
        barrier_return_a=float((gb @ phi).sum()),
        elliptic_residual_a=float(np.max(np.abs((li + le + gb) @ phi + li @ v
                                                - a * sheet.transepithelial_a_per_m2)))
        if sheet.bath == 'solved' else 0.0)
    out['battery_minus_barrier_a'] = out['battery_a'] - out['barrier_return_a']
    out['species_current_a_per_m2'] = {k: float(np.abs(x).max()) for k, x in per.items()}
    return out


def integrate(sheet, v0, *, duration_s, dt_s, applied_a_per_m2=0., leak_s_per_m2=0., blockade=None,
              barrier=None, ion_accounting=False, record_every=1, tol=1e-12, concentration_floor=1e-3):
    """Backward-Euler integration of the coupled bidomain + gated membrane.

    Timescale: duration_s up to MAX_DURATION_S (10 days), against the retained
    10 s ceiling. Without ion accounting the reversals are held fixed and the
    horizon is limited only by the honesty of a fixed-reservoir approximation.
    With ion_accounting=True the intracellular inventory is integrated from the
    channel and pump fluxes, the reversals follow it, and the run raises
    ReservoirExhausted at concentration_floor mol/m^3 -- the same failure mode as
    ihm/assembly/epithelial_electrodiffusion.transport.

    A bistable membrane can make one implicit step non-unique at large dt; the
    step is then subdivided by 4 (up to 5 levels) rather than the tolerance
    loosened. Every subdivision is counted in the audit.
    """
    if not math.isfinite(duration_s) or not 0 < duration_s <= MAX_DURATION_S:
        raise ValueError('duration must be in (0, %g] s' % MAX_DURATION_S)
    if not math.isfinite(dt_s) or dt_s <= 0:
        raise ValueError('positive dt required')
    steps = round(duration_s / dt_s)
    if abs(duration_s / dt_s - steps) > 1e-8 or not 1 <= steps <= MAX_STEPS:
        raise ValueError('dt must divide duration into at most %d steps' % MAX_STEPS)
    if not isinstance(record_every, int) or record_every < 1:
        raise ValueError('record_every must be a positive integer')
    n = sheet.n
    m = sheet.membrane
    a = sheet.node_area_m2
    cm = m.specific_capacitance_f_per_m2
    v = np.broadcast_to(np.asarray(v0, float), (n,)).astype(float).copy()
    applied = np.broadcast_to(np.asarray(applied_a_per_m2, float), (n,)).astype(float)
    leak = np.broadcast_to(np.asarray(leak_s_per_m2, float), (n,)).astype(float)
    inside = np.tile(np.asarray(m.inside_mol_per_m3, float), (n, 1))
    rev = np.tile(m.reversals(), (n, 1)) if ion_accounting else m.reversals()
    phi, _ = extracellular(sheet, v, blockade=blockade, barrier=barrier)
    gates_prev = {c.ident: c.gate.steady(v)[0] for c in m.channels
                  if c.gate is not None and m.law == 'gated'}
    frames = []
    worst = dict(newton_residual_a=0., elliptic_residual_a=0., node_charge_residual_c=0.,
                 battery_minus_barrier_a=0., newton_iterations=0, substep_refinements=0)
    z = np.array([1., 1., -1.])

    def snap(t, vv, pp, gg, cap):
        return dict(time_s=float(t), vm_v=vv.tolist(), phi_e_v=pp.tolist(),
                    vm_mean_v=float(vv.mean()), vm_sd_v=float(vv.std()),
                    vm_min_v=float(vv.min()), vm_max_v=float(vv.max()),
                    phi_e_min_v=float(pp.min()), phi_e_max_v=float(pp.max()),
                    receipts=_receipts(sheet, vv, pp, gg, rev, applied, leak, blockade, barrier, cap))

    frames.append(snap(0., v, phi, _gate_state(sheet, v, gates_prev, dt_s), None))
    charge0 = a * cm * v
    network_integral = np.zeros(n)
    ionic_integral = np.zeros(n)

    def substep(vs, ps, gp, h, ins, reversals, depth):
        """One backward-Euler substep. Returns the new state and its charge flow."""
        scale = a * cm / h
        vn, pn, gates, r, k = _solve_newton(sheet, vs, ps, gp, h, reversals, applied, leak,
                                            scale, scale * vs, blockade, barrier, tol=tol)
        li, _, _ = operators(sheet, vn + pn, blockade, barrier)
        i_ion, _, per = m.current(vn, gates, reversals)
        i_ion = i_ion + applied + leak * vn
        net = h * (-(li @ (vn + pn)))
        ion = h * (a * i_ion)
        gp2 = {kk: g for kk, g in gates.items() if '/dV' not in kk}
        ins2 = ins
        rev2 = reversals
        if ion_accounting:
            dn = np.stack([-a * per['Na'] / FARADAY, -a * per['K'] / FARADAY,
                           a * per['Cl'] / FARADAY], axis=1)
            ins2 = ins + h * dn / sheet.node_volume_m3[:, None]
            if (ins2 <= concentration_floor).any():
                bad = int(np.argmin(ins2.min(axis=1)))
                sp = ('Na', 'K', 'Cl')[int(np.argmin(ins2[bad]))]
                raise ReservoirExhausted(
                    'intracellular %s reached the %g mol/m^3 floor at node %d; step-end '
                    'concentrations Na/K/Cl = %s mol/m^3'
                    % (sp, concentration_floor, bad, np.round(ins2[bad], 6).tolist()))
            rev2 = m.thermal_v / z * np.log(np.asarray(m.outside_mol_per_m3, float) / ins2)
        return vn, pn, gp2, ins2, rev2, net, ion, r, k

    for step in range(steps):
        v_start = v
        net_sum = np.zeros(n)
        ion_sum = np.zeros(n)
        residual = 0.
        iters = 0
        queue = [(dt_s, 1)]
        while queue:
            h, depth = queue.pop(0)
            try:
                v, phi, gates_prev, inside, rev, net, ion, r, k = substep(
                    v, phi, gates_prev, h, inside, rev, depth)
            except ReservoirExhausted as exc:
                raise ReservoirExhausted('%s (t <= %.6g s)' % (exc, (step + 1) * dt_s)) from None
            except RuntimeError:
                if depth >= 5:
                    raise
                worst['substep_refinements'] += 1
                queue = [(h / 4., depth + 1)] * 4 + queue
                continue
            net_sum += net
            ion_sum += ion
            residual = max(residual, r)
            iters = max(iters, k + 1)
        network_integral += net_sum
        ionic_integral += ion_sum
        capacitive = a * cm * (v - v_start) / dt_s
        worst['newton_residual_a'] = max(worst['newton_residual_a'], residual)
        worst['newton_iterations'] = max(worst['newton_iterations'], iters)
        worst['node_charge_residual_c'] = max(
            worst['node_charge_residual_c'],
            float(np.max(np.abs(a * cm * (v - v_start) - (net_sum - ion_sum)))))
        if (step + 1) % record_every == 0 or step == steps - 1:
            f = snap((step + 1) * dt_s, v, phi, _gate_state(sheet, v, gates_prev, dt_s), capacitive)
            worst['elliptic_residual_a'] = max(worst['elliptic_residual_a'],
                                               f['receipts']['elliptic_residual_a'])
            worst['battery_minus_barrier_a'] = max(worst['battery_minus_barrier_a'],
                                                   abs(f['receipts']['battery_minus_barrier_a']))
            frames.append(f)
    charge = a * cm * v
    audit = dict(
        method='backward Euler; membrane current, gates and the elliptic constraint all evaluated at '
               'the new state; charge flows accumulated per substep',
        identity='A*C_m*dVm == integral of (intracellular network current - ionic current) dt, per node',
        capacitor_delta_c=float((charge - charge0).sum()),
        integrated_network_c=float(network_integral.sum()),
        integrated_ionic_c=float(ionic_integral.sum()),
        global_residual_c=float((charge - charge0).sum() - (network_integral - ionic_integral).sum()),
        max_node_charge_residual_c=worst['node_charge_residual_c'],
        network_sums_to_zero_c=float(network_integral.sum()),
        worst=worst)
    return dict(frames=frames, final_vm_v=v, final_phi_e_v=phi,
                final_inside_mol_per_m3=inside.tolist() if ion_accounting else None,
                reversal_v=np.asarray(rev).tolist(),
                charge_audit=audit, ion_accounting=bool(ion_accounting),
                clock=dict(dt_s=dt_s, duration_s=duration_s, steps=steps,
                           max_duration_s=MAX_DURATION_S, max_steps=MAX_STEPS))


def iv_curve(membrane, v_grid, applied=0.):
    v = np.asarray(v_grid, float)
    i, di, _ = membrane.current(v)
    return i + applied, di


def iv_roots(membrane, lo=-0.15, hi=0.12, samples=54001, applied=0.):
    """Bracket and bisect every zero of the space-clamped I(V); classify stability."""
    v = np.linspace(lo, hi, samples)
    i, _ = iv_curve(membrane, v, applied)
    idx = np.where(np.sign(i[:-1]) * np.sign(i[1:]) < 0)[0]
    out = []
    for k in idx:
        a, b = v[k], v[k + 1]
        fa = iv_curve(membrane, np.array([a]), applied)[0][0]
        for _ in range(200):
            mid = 0.5 * (a + b)
            fm = iv_curve(membrane, np.array([mid]), applied)[0][0]
            if fm == 0:
                a = b = mid
                break
            if np.sign(fm) == np.sign(fa):
                a, fa = mid, fm
            else:
                b = mid
        root = 0.5 * (a + b)
        slope = iv_curve(membrane, np.array([root]), applied)[1][0]
        out.append(dict(vm_v=float(root), di_dv_s_per_m2=float(slope),
                        stable=bool(slope > 0)))
    return out


def hysteresis(sheet, drive_a_per_m2, *, v_start=None, relax_s=5.0, relax_steps=25, tol=1e-12):
    """Quasi-static continuation of the coupled steady state along a drive sweep.

    Past a fold the branch being followed ceases to exist and a pure steady
    Newton has nowhere near to land, so the step is completed by relaxing the
    dynamics (backward Euler for relax_s seconds) and then polishing. Points
    that needed relaxation are reported, not hidden: they are exactly the
    switching events.
    """
    drive = np.asarray(drive_a_per_m2, float)
    v = np.full(sheet.n, -0.09 if v_start is None else float(v_start))
    branch, relaxed, residual = [], [], []
    for d in drive:
        used = False
        try:
            state = steady_state(sheet, v, applied_a_per_m2=float(d), tol=tol)
        except RuntimeError:
            used = True
            run = integrate(sheet, v, duration_s=relax_s, dt_s=relax_s / relax_steps,
                            applied_a_per_m2=float(d), record_every=relax_steps)
            v = run['final_vm_v']
            try:
                state = steady_state(sheet, v, applied_a_per_m2=float(d), tol=tol)
            except RuntimeError:
                state = dict(vm_v=v, newton_residual_a=float('nan'))
        v = state['vm_v']
        branch.append(float(v.mean()))
        relaxed.append(used)
        residual.append(float(state['newton_residual_a']))
    return dict(branch_vm_v=np.asarray(branch), relaxed=relaxed,
                max_newton_residual_a=float(np.nanmax(residual)))


# ---------------------------------------------------------------- geometry ----

def square_sheet(nx, ny, spacing_m, *, cell_height_m=1e-5):
    """Regular square lattice. Node membrane area = spacing^2 exactly, and the
    5-point Laplacian edge weight for a sheet conductance sigma_sheet (S) is
    exactly sigma_sheet on a square lattice (width/length = 1)."""
    ix, iy = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
    pos = np.stack([ix.ravel() * spacing_m, iy.ravel() * spacing_m, np.zeros(nx * ny)], axis=1)
    idx = np.arange(nx * ny).reshape(nx, ny)
    edges = np.concatenate([np.stack([idx[:-1, :].ravel(), idx[1:, :].ravel()], axis=1),
                            np.stack([idx[:, :-1].ravel(), idx[:, 1:].ravel()], axis=1)])
    area = np.full(nx * ny, spacing_m ** 2)
    return pos, area, area * cell_height_m, edges


def surface_sheet(vertices, faces, *, layer_thickness_m=1e-4):
    """Cotangent-weighted Laplacian on a triangulated surface, with a barycentric
    per-node area and a per-node volume = area x layer thickness.

    [literature] cotangent discretisation of the Laplace-Beltrami operator:
    MacNeal 1949; Pinkall & Polthier 1993 Exp Math 2:15; Meyer et al. 2003.
    Edge weight for a sheet conductivity sigma_sheet (S) is
    sigma_sheet * (cot alpha + cot beta) / 2.
    """
    v = np.asarray(vertices, float)
    f = np.asarray(faces, int)
    n = len(v)
    e0 = v[f[:, 2]] - v[f[:, 1]]
    e1 = v[f[:, 0]] - v[f[:, 2]]
    e2 = v[f[:, 1]] - v[f[:, 0]]
    cross = np.cross(e2, -e1)
    twice_area = np.linalg.norm(cross, axis=1)
    tri_area = 0.5 * twice_area
    cot = np.stack([(-np.einsum('ij,ij->i', e1, e2)),
                    (-np.einsum('ij,ij->i', e2, e0)),
                    (-np.einsum('ij,ij->i', e0, e1))], axis=1) / np.maximum(twice_area, 1e-300)[:, None]
    area = np.zeros(n)
    np.add.at(area, f.ravel(), np.repeat(tri_area / 3., 3))
    pairs = np.concatenate([f[:, [1, 2]], f[:, [2, 0]], f[:, [0, 1]]])
    w = np.concatenate([cot[:, 0], cot[:, 1], cot[:, 2]]) / 2.
    key = np.sort(pairs, axis=1)
    order = np.lexsort((key[:, 1], key[:, 0]))
    key, w = key[order], w[order]
    uniq, start = np.unique(key, axis=0, return_index=True)
    weight = np.add.reduceat(w, np.sort(start))
    return v, area, area * layer_thickness_m, uniq, weight


# ------------------------------------------------------ standard membranes ----

def linear_membrane(*, g_na=1e-2, g_k=1e-1, g_cl=2e-2, pump_a_per_m2=1e-3):
    """The RETAINED law, re-expressed per unit area and nothing else.

    [retained] SkinPatch cell defaults are C = 10 pF, g_Na/g_K/g_Cl =
    10/100/20 pS, pump 1 pA. Dividing by the membrane area those lumped values
    imply at the literature 1 uF/cm^2 (10 pF -> 1e-9 m^2) gives exactly the
    specific values above. The dynamics are unchanged; the parameters become
    falsifiable because an area is now declared.
    """
    return Membrane(channels=(Channel('Na_leak', 'Na', g_na, None, 'retained'),
                              Channel('K_leak', 'K', g_k, None, 'retained'),
                              Channel('Cl_leak', 'Cl', g_cl, None, 'retained')),
                    pump_cycles_mol_per_m2_s=pump_a_per_m2 / FARADAY, law='linear')


def bistable_membrane(*, g_kir=0.10, v_kir=-50e-3, s_kir=10e-3,
                      g_dep=4e-2, v_dep=-25e-3, s_dep=6e-3,
                      g_cl=2e-2, pump_a_per_m2=1e-3, tau_dep_s=0.0):
    """Two Boltzmann-gated conductances plus the retained Cl leak and pump.

    [literature] the FORM: an inward-rectifier K conductance that closes on
    depolarisation together with a depolarising conductance that opens on
    depolarisation gives an N-shaped I(V) with two stable zeros. This is the
    standard minimal bistable non-excitable membrane (Cervera, Alcaraz & Mafe
    2014 J Phys Chem B 118:6417; Cervera, Meseguer & Mafe 2016 Sci Rep 6:35201;
    Law & Levin 2015 Theor Biol Med Model 12:22), and the inward-rectifier
    voltage dependence is Hille ch.14.
    [retained] g_cl and the pump current density are the SkinPatch values
    re-expressed per area; the reversals are the retained concentrations.
    [engineering] g_kir, v_kir, s_kir, g_dep, v_dep, s_dep are CHOSEN. They are
    of literature order for keratinocyte Kir and for a depolarising cation
    conductance, but no keratinocyte measurement identifies them. They place the
    membrane in the bistable regime; nothing here claims human keratinocytes are
    bistable at these values.
    """
    return Membrane(channels=(
        Channel('Kir', 'K', g_kir, Gate(v_kir, s_kir), 'engineering'),
        Channel('Cat_dep', 'Na', g_dep, Gate(v_dep, -s_dep, tau_dep_s), 'engineering'),
        Channel('Cl_leak', 'Cl', g_cl, None, 'retained')),
        pump_cycles_mol_per_m2_s=pump_a_per_m2 / FARADAY, law='gated')


def build_sheet(nx, ny, spacing_m, membrane, *, gj_sheet_s=5.0e-7, ecm_sheet_s=3.0e-6,
                barrier_s_per_m2=2.0, tep_v=30e-3, cell_height_m=1e-5,
                keratinocyte_density_per_m2=1.3463086096010574e10, gj_gate=None,
                bath='solved'):
    """A coarse-grained square patch with a fully declared area/volume ledger."""
    pos, area, volume, edges = square_sheet(nx, ny, spacing_m, cell_height_m=cell_height_m)
    n = len(pos)
    j_te = barrier_s_per_m2 * tep_v
    cg = dict(
        node_membrane_area_m2=float(spacing_m ** 2),
        node_intracellular_volume_m3=float(spacing_m ** 2 * cell_height_m),
        keratinocyte_areal_density_per_m2=float(keratinocyte_density_per_m2),
        keratinocytes_per_node=float(keratinocyte_density_per_m2 * spacing_m ** 2),
        density_evidence='[retained] recomputed from the BETSE polygon areas in '
                         'data/derived/app/geometry/betse-tissue-cells.json.gz',
        statement='One node is a POPULATION of keratinocytes sharing one membrane potential. '
                  'Vm at a node is a population mean. Any single-cell claim at this resolution '
                  'is a category error.')
    return Sheet(pos, area, volume, edges,
                 np.full(len(edges), gj_sheet_s), np.full(len(edges), ecm_sheet_s),
                 np.full(n, barrier_s_per_m2), np.full(n, j_te), membrane,
                 gj_gate=gj_gate, bath=bath, coarse_graining=cg,
                 provenance=dict(geometry='synthetic square lattice %dx%d at %g m' % (nx, ny, spacing_m),
                                 gj_sheet_s='[engineering] %g S' % gj_sheet_s,
                                 ecm_sheet_s='[engineering] %g S' % ecm_sheet_s,
                                 barrier_s_per_m2='[literature order] %g S/m^2 = %g Ohm.cm^2 TEER'
                                                  % (barrier_s_per_m2, 1e4 / barrier_s_per_m2),
                                 tep_v='[literature] %g V; Barker, Jaffe & Vanable 1982 '
                                       'Am J Physiol 242:R358 report 15-50 mV skin TEP' % tep_v))


def length_constants(sheet):
    """Analytic receipts that say WHY the numbers come out as they do."""
    m = sheet.membrane
    g_m = sum(c.g_max_s_per_m2 for c in m.channels)
    h = float(np.sqrt(sheet.node_area_m2.mean()))
    gj = float(np.mean(sheet.gj_s))
    ecm = float(np.mean(sheet.ecm_s))
    gb = float(np.mean(sheet.barrier_s_per_m2))
    lam_i = math.sqrt(gj / g_m) if g_m > 0 else float('inf')
    lam_e = math.sqrt(ecm / gb) if gb > 0 else float('inf')
    return dict(node_spacing_m=h, membrane_conductance_s_per_m2=g_m,
                intracellular_length_constant_m=lam_i,
                extracellular_length_constant_m=lam_e,
                lam_i_over_spacing=lam_i / h, lam_e_over_spacing=lam_e / h,
                vm_response_ratio_lam_i_over_lam_e_squared=(lam_i / lam_e) ** 2 if lam_e > 0 else float('inf'),
                derivation='steady cable: -sigma grad^2 Vm + g_m Vm = g_m V_rest + sigma grad^2 phi_e, so a '
                           'change in phi_e reaches Vm only through its CURVATURE and the transfer is '
                           'O((lam_i/lam_e)^2). This is the analytic reason a wound moves Vm by little '
                           'when the intracellular sheet is far less coupled than the extracellular one.')
