"""Energy, force and penetration audit of the production supine blanket.

Nothing here accepts motion as evidence. The body is held stationary and the
supports are fixed, so the blanket's total mechanical energy may not rise, the
contact models must return equal and opposite impulses, the self-collision
passes must carry no net force or torque, and the sheet must end every step
without passing through itself or through the sampled skin envelope.

This audits the PRODUCTION explicit environment path. It is not a passivity
proof: the alternating position projections still do net positive work, which
the damped integrator then removes, and the numbers below report both. The
separately verified implicit path (``verify_cloth_passive_skin.py``) is the one
with a per-step energy and momentum acceptance test. Neither is volumetric
soft-body FEM or a calibrated textile model.
"""
import argparse, hashlib, json, tempfile, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCOPE = ('Held canonical body, fixed supports, production explicit cloth path with hinge bending, '
         'triangle/edge self collision and minimum-norm envelope repair. Engineering coefficients, '
         'sampled sphere skin envelope, 437-node sheet. Not volumetric soft-body FEM, not a calibrated '
         'textile model, and not a per-step passivity proof of the projection stage.')

HISTORICAL_BASELINE = {
    'source': 'data/derived/unified-world-b_xo1l9j/cloth_jitter_diagnosis.json',
    'note': 'Same held-body fixture and 500 one-millisecond substeps, measured before the hinge '
            'bending term, the minimum-norm envelope repair and the strain-sweep increase.',
    'mechanical_energy_change_j': 0.9354982300645869,
    'mechanical_energy_trend_w': 0.976089,
    'mechanical_energy_trend_standard_error_w': 0.107320,
    'peak_transient_energy_rise_j': 1.43721,
    'projection_elastic_work_j': 17.94228183235141,
    'projection_unresolved_work_j': 11.06396325967851,
    'mean_peak_speed_last_100_substeps_m_s': 3.422690097909144}

# Reported, deliberately NOT asserted, because it is still not true. The
# explicit production path is not a passive solver: its energy trend is still
# significantly positive and individual substeps still gain energy from the
# alternating position projections. The required check below is an explicit
# regression bound against the retained baseline, not a passivity claim. The
# implicit path verified by verify_cloth_passive_skin.py is the one with a
# per-step energy and momentum acceptance test.
NOT_YET_MET = ('per_step_passivity_of_the_projection_stage',)
# One percent of the retained baseline energy growth rate.
ENERGY_TREND_REGRESSION_BOUND_W = .01*HISTORICAL_BASELINE['mechanical_energy_trend_w']

TOLERANCE = {'mechanical_energy_trend_w': ENERGY_TREND_REGRESSION_BOUND_W,
             'cover_reciprocity_impulse_ns': 1e-12,
             'self_collision_momentum_ns': 1e-12,
             'self_collision_torque_nms': 1e-9,
             'skin_penetration_m': 1e-9,
             'extension_ratio': 1.121,
             'max_fold_angle_deg': 90.}


def fixture():
    import ihm.assembly.environment_dynamics as env
    from ihm.assembly.articulated import CanonicalRegistration
    payload = json.loads((ROOT/'data/derived/canonical/mechanics.json').read_bytes())
    specs = {e['id']: e for e in payload['entities']}
    groups = {}
    for name, row in payload['registration'].items():
        bones = [specs[i] for i in row.get('canonical_bones', []) if i in specs]
        if bones:
            groups[name] = {'canonical_bones': [e['id'] for e in bones],
                            'bounds_min_m': np.min([e['bounds_m']['min'] for e in bones], axis=0),
                            'bounds_max_m': np.max([e['bounds_m']['max'] for e in bones], axis=0)}
    registration = SimpleNamespace(specs=specs, groups=groups)
    registration._ranking = lambda point: CanonicalRegistration._ranking(registration, point)
    entities = {i: {'centroid_m': e['centroid_m'], 'rotation_matrix': np.eye(3).tolist()}
                for i, e in specs.items()}
    owner = env.EnvironmentDynamics(ROOT, 'supine', {'scene': 'bedroom', 'objects': []}, registration)
    return env, owner, entities


def instrument(mesh):
    """Wrap the cover projection to measure cloth/body impulse reciprocity."""
    from ihm.assembly import cloth_cover
    original = cloth_cover.SupineClothCover.project
    record = {'max_reciprocity_ns': 0., 'contact_impulse_ns': np.zeros(3)}

    def traced(self, target, h, friction=.45):
        before = target.v.copy()
        impulses = original(self, target, h, friction=friction)
        mass = np.broadcast_to(np.asarray(target.mass, float), (len(target.x),))
        cloth = np.sum(mass[:, None]*(target.v-before), axis=0)
        skin = np.sum(impulses, axis=0)
        record['max_reciprocity_ns'] = max(record['max_reciprocity_ns'],
                                           float(np.abs(cloth+skin).max()))
        record['contact_impulse_ns'] = record['contact_impulse_ns']+skin
        return impulses

    cloth_cover.SupineClothCover.project = traced
    return record, original


def energy(mesh, gravity):
    return float(.5*mesh.mass*np.sum(mesh.v**2)+mesh.elastic_energy()
                 - mesh.mass*float(np.sum(mesh.x@gravity)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=.5)
    parser.add_argument('--output', default=None,
                        help='Retained receipt directory. A fresh temporary directory under '
                             'data/derived is created when omitted; never pass a source root.')
    arguments = parser.parse_args()
    start = time.monotonic()
    output = Path(arguments.output) if arguments.output else Path(
        tempfile.mkdtemp(prefix='cloth-energy-audit-', dir=ROOT/'data/derived'))
    output.mkdir(parents=True, exist_ok=True)
    env, owner, entities = fixture()
    mesh = next(o for o in owner.soft if o.kind == 'cloth')
    from ihm.assembly.cloth_cover import SupineClothCover
    cover = SupineClothCover(owner.skin_points(entities))
    reciprocity, original = instrument(mesh)
    sources = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in (
        'ihm/assembly/environment_dynamics.py', 'ihm/assembly/cloth_bending.py',
        'ihm/assembly/cloth_self_collision.py', 'ihm/assembly/cloth_cover.py',
        'ihm/assembly/cloth_stretch.py', 'ihm/assembly/cloth_contact.py',
        'scripts/verify_cloth_energy_audit.py')}
    report = {'passed': False, 'scope': SCOPE, 'tolerance': TOLERANCE, 'source_sha256': sources,
              'preparation': dict(mesh.preparation),
              'parameters': {'bending_rigidity_n_m': env.CLOTH_BENDING_RIGIDITY_N_M,
                             'sheet_thickness_m': env.CLOTH_THICKNESS_M,
                             'stretch_sweeps': env.CLOTH_STRETCH_SWEEPS,
                             'cover_passes': env.CLOTH_COVER_PASSES,
                             'vertex_count': int(len(mesh.x)), 'substep_s': .001}}
    rows = []
    history = []
    try:
        initial = energy(mesh, owner.gravity)
        totals = {'self_collision_momentum_ns': np.zeros(3), 'self_collision_torque_nms': 0.,
                  'self_collision_kinetic_j': 0., 'projection_unresolved_j': 0.,
                  'projection_spring_j': 0., 'projection_kinetic_j': 0.}
        speeds = []
        unconverged = 0
        worst = {'energy_j': initial, 'extension_ratio': 0., 'fold_angle_deg': 0.,
                 'skin_penetration_m': 0., 'speed_m_s': 0., 'self_collision_kinetic_j': 0.,
                 'min_primitive_distance_m': float('inf')}
        for index in range(round(arguments.seconds/.001)):
            owner.advance(.001, entities)
            state = mesh.constraint_state
            triangle = mesh.triangle_contact_state
            totals['projection_unresolved_j'] += state['unresolved_energy_change_j']
            totals['projection_spring_j'] += state['spring_projection_energy_change_j']
            totals['projection_kinetic_j'] += state['kinetic_energy_change_j']
            if triangle.get('evaluated'):
                totals['self_collision_momentum_ns'] += np.asarray(triangle['linear_momentum_change_ns'])
                totals['self_collision_torque_nms'] += triangle['impulse_torque_residual_nms']
                totals['self_collision_kinetic_j'] += triangle['kinetic_energy_change_j']
                worst['self_collision_kinetic_j'] = max(worst['self_collision_kinetic_j'],
                                                        triangle['kinetic_energy_change_j'])
                worst['min_primitive_distance_m'] = min(worst['min_primitive_distance_m'],
                                                        triangle['residual_min_distance_m'])
            depth, _, _ = cover.query(mesh.x)
            current = energy(mesh, owner.gravity)
            history.append(current)
            speeds.append(float(np.linalg.norm(mesh.v, axis=1).max()))
            unconverged += 0 if mesh.stretch_state['converged'] else 1
            worst['energy_j'] = max(worst['energy_j'], current)
            worst['extension_ratio'] = max(worst['extension_ratio'], mesh.stretch_state['max_extension_ratio'])
            worst['fold_angle_deg'] = max(worst['fold_angle_deg'], mesh.bending_state['max_fold_angle_deg'])
            worst['skin_penetration_m'] = max(worst['skin_penetration_m'], float(depth.max(initial=0.)))
            worst['speed_m_s'] = max(worst['speed_m_s'], float(np.linalg.norm(mesh.v, axis=1).max()))
            if (index+1) % 50 == 0:
                rows.append({'time_s': (index+1)*.001, 'mechanical_energy_j': current,
                             'max_speed_m_s': float(np.linalg.norm(mesh.v, axis=1).max()),
                             'max_extension_ratio': mesh.stretch_state['max_extension_ratio'],
                             'stretch_converged': mesh.stretch_state['converged'],
                             'bending': mesh.bending_state, 'triangle_self_contact': triangle,
                             'constraint_ledger': state})
                print(json.dumps({k: rows[-1][k] for k in ('time_s', 'mechanical_energy_j',
                      'max_speed_m_s', 'max_extension_ratio')}), flush=True)
        final = energy(mesh, owner.gravity)
        # The energy oscillates, so a quarter-mean difference is dominated by
        # where the window cuts the oscillation. Fit a trend instead and ask
        # whether it is significantly positive against its own residual scatter.
        # That criterion calibrates itself; it is not a chosen threshold.
        clock = np.arange(len(history))*.001
        slope, intercept = np.polyfit(clock, history, 1)
        residual = np.asarray(history)-(slope*clock+intercept)
        spread = np.sum((clock-clock.mean())**2)
        error = float(np.sqrt(np.sum(residual**2)/max(1, len(history)-2)/max(spread, 1e-300)))
        drift = float(slope)
        measured = {'initial_mechanical_energy_j': initial, 'final_mechanical_energy_j': final,
                    'mechanical_energy_change_j': final-initial,
                    'mechanical_energy_trend_w': drift,
                    'mechanical_energy_trend_standard_error_w': error,
                    'peak_transient_energy_rise_j': worst['energy_j']-initial,
                    'substeps': len(history),
                    'substeps_reporting_unconverged_strain': unconverged,
                    'max_cover_reciprocity_residual_ns': reciprocity['max_reciprocity_ns'],
                    'total_skin_contact_impulse_ns': reciprocity['contact_impulse_ns'].tolist(),
                    'self_collision_net_momentum_ns': float(np.abs(totals['self_collision_momentum_ns']).max()),
                    'self_collision_torque_residual_nms': totals['self_collision_torque_nms'],
                    'self_collision_kinetic_energy_j': totals['self_collision_kinetic_j'],
                    'worst_single_self_collision_kinetic_j': worst['self_collision_kinetic_j'],
                    'min_primitive_distance_m': worst['min_primitive_distance_m'],
                    'max_skin_penetration_m': worst['skin_penetration_m'],
                    'max_extension_ratio': worst['extension_ratio'],
                    'max_fold_angle_deg': worst['fold_angle_deg'],
                    'max_speed_m_s': worst['speed_m_s'],
                    'mean_peak_speed_last_100_substeps_m_s': float(np.mean(speeds[-100:])),
                    'projection_unresolved_work_j': totals['projection_unresolved_j'],
                    'projection_elastic_work_j': totals['projection_spring_j'],
                    'projection_kinetic_work_j': totals['projection_kinetic_j'],
                    'final_bending': mesh.bending_state,
                    'final_triangle_self_contact': mesh.triangle_contact_state}
        report['measured'] = measured
        checks = {
            'energy_growth_rate_is_at_most_one_percent_of_the_retained_baseline':
                measured['mechanical_energy_trend_w'] <= TOLERANCE['mechanical_energy_trend_w'],
            'cloth_body_contact_impulses_are_reciprocal':
                measured['max_cover_reciprocity_residual_ns'] <= TOLERANCE['cover_reciprocity_impulse_ns'],
            'self_collision_carries_no_net_force':
                measured['self_collision_net_momentum_ns'] <= TOLERANCE['self_collision_momentum_ns'],
            'self_collision_carries_no_net_torque':
                measured['self_collision_torque_residual_nms'] <= TOLERANCE['self_collision_torque_nms'],
            'self_collision_never_adds_kinetic_energy':
                measured['worst_single_self_collision_kinetic_j'] <= 0.,
            'sheet_does_not_pass_through_itself':
                measured['min_primitive_distance_m'] >= env.CLOTH_THICKNESS_M-1e-9,
            'sheet_does_not_penetrate_the_skin_envelope':
                measured['max_skin_penetration_m'] <= TOLERANCE['skin_penetration_m'],
            'material_strain_limit_holds':
                measured['max_extension_ratio'] <= TOLERANCE['extension_ratio'],
            'strain_solver_reports_convergence_on_every_substep':
                measured['substeps_reporting_unconverged_strain'] == 0,
            'no_sharp_creases_remain':
                measured['max_fold_angle_deg'] <= TOLERANCE['max_fold_angle_deg'],
        }
        report['checks'] = checks
        report['not_yet_met'] = {
            'per_step_passivity_of_the_projection_stage': {
                'mechanical_energy_trend_w': measured['mechanical_energy_trend_w'],
                'mechanical_energy_trend_standard_error_w': measured['mechanical_energy_trend_standard_error_w'],
                'upward_trend_is_statistically_significant':
                    bool(measured['mechanical_energy_trend_w'] > 2*measured['mechanical_energy_trend_standard_error_w']),
                'peak_transient_energy_rise_j': measured['peak_transient_energy_rise_j'],
                'projection_unresolved_work_j': measured['projection_unresolved_work_j'],
                'statement': 'The energy trend is still significantly positive and individual substeps '
                             'still gain energy from the alternating position projections, so this path '
                             'is NOT a passive cloth solver. The required check above is a regression '
                             'bound against the retained baseline, not a passivity claim. The implicit '
                             'path in verify_cloth_passive_skin.py is the one that audits every step.'}}
        report['historical_baseline'] = HISTORICAL_BASELINE
        failed = sorted(name for name, ok in checks.items() if not ok)
        report['failed_checks'] = failed
        report['passed'] = not failed
    except Exception as error:
        report['error'] = f'{type(error).__name__}: {error}'
    finally:
        from ihm.assembly import cloth_cover
        cloth_cover.SupineClothCover.project = original
        report['wall_s'] = time.monotonic()-start
        report['rows'] = rows
        report['mechanical_energy_series_j'] = [float(v) for v in history]
        report['final_positions'] = mesh.x.tolist()
        (output/'receipt.json').write_text(json.dumps(report, indent=2, default=float))
        print(json.dumps({k: v for k, v in report.items() if k not in ('rows', 'final_positions')},
                         indent=2, default=float), flush=True)
        print(str(output), flush=True)
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
