"""Pinned regional non-neural skin RC experiment with explicit reservoir ledger.

Coefficients are the existing illustrative SkinPatch coefficients, never inferred
from the human lateral-field regression. Membrane and apical networks are separate.
"""
from dataclasses import asdict, replace
import hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from ihm.materialize.skin import SkinPatch
from .regional_touch import forearm_anchor


def electrical_system(patch, electrode_a=None):
    """Re-express source Ohmic/Nernst laws with fixed ionic reservoirs in SI."""
    nc = len(patch.cells); n = nc + len(patch.surface)
    c = np.array([x.capacitance_f for x in (*patch.cells, *patch.surface)])
    a = np.zeros((n, n)); b = np.zeros(n)
    thermal = 8.314462618 * patch.temperature_k / 96485.33212
    for i, cell in enumerate(patch.cells):
        ions = ((1, cell.sodium_inside_mM, patch.sodium_outside_mM, cell.sodium_conductance_s),
                (1, cell.potassium_inside_mM, patch.potassium_outside_mM, cell.potassium_conductance_s),
                (-1, cell.chloride_inside_mM, patch.chloride_outside_mM, cell.chloride_conductance_s))
        a[i, i] = -sum(g for _, _, _, g in ions)
        b[i] = sum(g * thermal / z * np.log(outside / inside) for z, inside, outside, g in ions) - cell.pump_current_a
    for j, site in enumerate(patch.surface):
        a[nc+j, nc+j] = -site.barrier_conductance_s
        b[nc+j] = site.pump_current_a
    if electrode_a is not None:
        electrode_a = np.asarray(electrode_a, float)
        if electrode_a.shape != (len(patch.surface),) or not np.isfinite(electrode_a).all():
            raise ValueError('finite electrode currents required at all surface sites')
        b[nc:] += electrode_a
    for edges, offset in ((patch.gap_edges, 0), (patch.surface_edges, nc)):
        for u, v, g in edges:
            u += offset; v += offset
            a[u, v] += g; a[v, u] += g; a[u, u] -= g; a[v, v] -= g
    return a / c[:, None], b / c, c


def _experiment(patch, initial, dt, duration, electrode):
    a, b, capacitance = electrical_system(patch, electrode)
    n = len(b); nc = len(patch.cells)
    generator = np.zeros((n+1, n+1)); generator[:n, :n] = a; generator[:n, n] = b
    step = expm(generator * dt)
    nodes, weights = np.polynomial.legendre.leggauss(12)
    quadrature = [expm(generator * dt * (t+1)/2) for t in nodes]
    x = np.r_[initial, 1.]; integrated = {k: np.zeros(n) for k in ('ionic_leak', 'pump', 'barrier', 'electrode', 'internal_edges')}
    max_pair = 0.
    thermal = 8.314462618 * patch.temperature_k / 96485.33212

    def currents(v):
        channels = {k: np.zeros(n) for k in integrated}
        for i, cell in enumerate(patch.cells):
            for z, inside, outside, g in ((1, cell.sodium_inside_mM, patch.sodium_outside_mM, cell.sodium_conductance_s),
                    (1, cell.potassium_inside_mM, patch.potassium_outside_mM, cell.potassium_conductance_s),
                    (-1, cell.chloride_inside_mM, patch.chloride_outside_mM, cell.chloride_conductance_s)):
                channels['ionic_leak'][i] += g * (thermal/z*np.log(outside/inside)-v[i])
            channels['pump'][i] = -cell.pump_current_a
        for j, site in enumerate(patch.surface):
            channels['barrier'][nc+j] = -site.barrier_conductance_s*v[nc+j]
            channels['pump'][nc+j] = site.pump_current_a
            channels['electrode'][nc+j] = electrode[j]
        for edges, offset in ((patch.gap_edges, 0), (patch.surface_edges, nc)):
            for u, w, g in edges:
                u += offset; w += offset
                current = g*(v[w]-v[u])
                channels['internal_edges'][u] += current
                channels['internal_edges'][w] -= current
        return channels

    def frame(t):
        surface = x[nc:n]
        field = [(surface[u]-surface[v])/np.linalg.norm(np.subtract(patch.surface[v].position_m, patch.surface[u].position_m))
                 for u, v, _ in patch.surface_edges]
        return dict(time_s=t, membrane_voltage_V=x[:nc].tolist(), apical_voltage_V=surface.tolist(), edge_field_V_m=field)

    frames = [frame(0.)]
    for i in range(round(duration/dt)):
        for w, transition in zip(weights, quadrature):
            channel = currents((transition @ x)[:n])
            max_pair = max(max_pair, abs(float(channel['internal_edges'].sum())))
            for k in integrated: integrated[k] += w*dt/2*channel[k]
        x = step @ x
        frames.append(frame((i+1)*dt))
    delta = capacitance*(x[:n]-initial)
    balance = delta-sum(integrated.values())
    external = sum(v for k, v in integrated.items() if k != 'internal_edges')
    return dict(frames=frames, electrode_current_a=list(electrode), source_patch=asdict(patch),
                charge_audit=dict(capacitor_delta_c=delta.tolist(), integrated_currents_c={k:v.tolist() for k,v in integrated.items()},
                    node_residual_c=balance.tolist(), max_node_residual_c=float(max(abs(balance))),
                    global_residual_c=float(delta.sum()-external.sum()), max_internal_pair_residual_a=max_pair,
                    method='Exact matrix exponential states; independent 12-point Gauss current quadrature per interval; positive current enters capacitor',
                    reservoirs='Fixed intracellular Na+, K+, Cl- and extracellular baths; prescribed pumps and electrode supply exchange charge. No ion mass-balance or finite reservoir depletion.'))


def build_skin_electric(root, *, dt_s=.05, duration_s=2.):
    if not np.isfinite([dt_s, duration_s]).all() or not 0 < dt_s <= .1 or not 0 < duration_s <= 10 or abs(duration_s/dt_s-round(duration_s/dt_s)) > 1e-8:
        raise ValueError('finite bounded clock aligned to duration required')
    root = Path(root); anchor = forearm_anchor(root)
    patch = SkinPatch.line(n=9, wound=False)
    a, b, _ = electrical_system(patch)
    initial = np.linalg.solve(a, -b)
    shunt = replace(patch, surface=tuple(replace(s, barrier_conductance_s=1e-4, pump_current_a=0., wounded=True) if i==4 else s for i,s in enumerate(patch.surface)))
    electrode = np.zeros(9); electrode[0] = 3e-8; electrode[-1] = -3e-8
    perturbed = initial.copy(); perturbed[0] += .01
    experiments = {name: _experiment(p, state, dt_s, duration_s, current) for name, p, state, current in (
        ('baseline', patch, initial, np.zeros(9)), ('wound_shunt', shunt, initial, np.zeros(9)),
        ('electrode', patch, initial, electrode), ('membrane_perturbation', patch, perturbed, np.zeros(9)))}
    world = lambda points: ((np.asarray(points)-[.0008,0,0])@np.asarray(anchor['local_axes'])+anchor['origin_m']).tolist()
    paths = ['ihm/assembly/skin_bioelectric.py', 'ihm/assembly/regional_touch.py', 'ihm/materialize/skin.py',
             'ihm/calibration/skin.py', 'scripts/build_body_skin_electric.py', 'scripts/verify_body_skin_electric.py',
             'docs/research/REGIONAL_SKIN_ELECTRIC.md']
    sources = ['data/sources/human-skin-wound-fields.json', 'data/sources/human-skin-field-age.json',
               'data/derived/integumentary/human-wound-observations.jsonl', 'data/derived/calibration/skin-fit.json']
    return dict(schema_version=1, id='forearm-skin-electric', anchor=anchor,
        geometry=dict(cell_positions_m=world([c.position_m for c in patch.cells]), surface_positions_m=world([s.position_m for s in patch.surface]),
            gap_edges=[list(e) for e in patch.gap_edges], surface_edges=[list(e) for e in patch.surface_edges],
            domain_kind='synthetic 1.6 mm tangent line, membrane layer depth 0.1 mm; fixed pinned face reference attachment'),
        clock=dict(dt_s=dt_s, duration_s=duration_s, intervention_onset_s=0., semantics='Independent experiments start at equilibrium; interventions held from t=0; membrane perturbation is an initial condition'),
        experiments=experiments, recorded_baseline=experiments['baseline']['frames'][0],
        voltage_reference='Vm = intracellular minus basal extracellular bath (0 V); apical = apical bath minus basal bath; TEP = -apical',
        parameter_status='All RC, pump, concentration, geometry and intervention values are explicit illustrative priors inherited from SkinPatch.line or specified here; no human Vm calibration',
        intervention_priors=dict(electrode_current_a=electrode.tolist(), membrane_initial_delta_V=.01, shunt_conductance_s=1e-4),
        runtime_sources={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths},
        source_hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sources},
        observation_binding=dict(expression='E_ab = (apical[a] - apical[b]) / distance_m', unit='V/m',
            human_fit_path='data/derived/calibration/skin-fit.json', calibration_scope='Human lateral wound-field phenotype only; no cell Vm or RC parameter identification'),
        coupling=dict(native_blood_storage_connected=False, neural_connection=False, whole_body_skin_activation=False),
        limitations=['Fixed intracellular and extracellular ion reservoirs; no ion mass-balance or electroneutrality solution.',
            'Prescribed basal bath decouples intracellular gap-junction and apical extracellular graphs.',
            'Barrier shunt marks an electrical wound surrogate; cells retained, no ablation, migration, healing or regeneration law.',
            'Deterministic mean trajectory; source Gaussian covariance is not propagated or claimed as parameter uncertainty.',
            'Pinned planar line is a regional alternative materialization, not curved skin coverage or whole-body skin activation.',
            'Human field measurements support the observable only and do not calibrate this patch or its membrane voltage.'])


def native_skin_ionic_boundary(exchange):
    """Observed bulk native skin concentrations and their Nernst potentials.

    This does not replace keratinocyte concentrations or infer membrane voltage.
    No current or ion mass is returned to native state. Any epithelial transfer
    requires an explicit additional constitutive assumption.
    """
    import math
    tissue = exchange['native_compartments']
    inside = tissue['Skin.intracellular']['ionic_molarity_mmol_per_l']
    outside = tissue['Skin.extracellular']['ionic_molarity_mmol_per_l']
    celsius = exchange['thermal_boundary_c']['skin']
    if isinstance(celsius, bool) or not isinstance(celsius, (int, float)) or not math.isfinite(celsius) or celsius <= -273.15:
        raise ValueError('Native skin temperature is missing or invalid')
    kelvin = celsius + 273.15
    nernst = {}
    for ion, charge in [('Sodium', 1), ('Potassium', 1), ('Chloride', -1)]:
        a, b = inside.get(ion), outside.get(ion)
        if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0 for v in (a, b)):
            raise ValueError('Positive native intracellular and extracellular '+ion+' are required')
        nernst[ion] = 8.314462618 * kelvin / (charge * 96485.33212) * math.log(b / a)
    return {'time_s': exchange['time_s'], 'temperature_k': kelvin,
            'inside_mmol_per_l': dict(inside), 'outside_mmol_per_l': dict(outside),
            'nernst_potential_v': nernst, 'owner': 'native bulk Skin tissue',
            'membrane_voltage_predicted': False, 'native_mass_mutated': False,
            'limitation': 'Bulk tissue gradients are not measured keratinocyte or nerve membrane states.'}
