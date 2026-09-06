"""Observe installed native Skin leaves; never integrate or allocate their state."""
from copy import deepcopy
import math

PREFIX = 'tissue.regional_skin.'
REGIONS = ('region_a', 'region_b', 'residual')
FRACTIONS = (.2, .3, .5)  # Engineering installation allocation, not anatomy.
PATHS = ('SkinVToSkinE1', 'SkinE1ToSkinE2', 'SkinE2ToSkinE3',
         'SkinE3ToGround', 'SkinE3ToSkinI', 'SkinE3ToSkinL1',
         'SkinL1ToSkinL2', 'SkinToLymphValve', 'SkinSweating')


def _number(value, key, nonnegative=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or (nonnegative and value < 0)):
        raise ValueError('Invalid regional native quantity: ' + key)
    return float(value)


def observe_regional_skin(snapshot, *, native_identity=None, source_hashes=None):
    """Return mergeable native owners and a non-owning aggregate, or unavailable.

    A partially exported installation fails closed. Missing species observations
    remain unknown. Conservation checks compare observations, not inferred flows
    of solutes. Native last-step residuals are not interval integrations.
    """
    values = snapshot.get('values', {})
    raw = {k: v for k, v in values.items() if k.startswith(PREFIX)}
    if not raw:
        return {'schema': 'native_regional_skin_exchange_v1', 'available': False}
    for key in raw:
        scope = key[len(PREFIX):].split('.')[0]
        if scope not in (*REGIONS, 'aggregate', 'sweat', 'completed_native_steps'):
            raise ValueError('Unknown regional owner namespace: ' + scope)
    time = _number(snapshot.get('time_s'), 'time_s', True)
    unknown = []

    def read(suffix, required=True, nonnegative=False):
        key = PREFIX + suffix
        value = values.get(key)
        if value is None and not required:
            unknown.append(key)
            return None
        return _number(value, key, nonnegative)

    checks = []

    def check(label, residual, scale=0., unit='mL'):
        tolerance = max(1e-8, abs(scale) * 1e-10)
        checks.append({'check': label, 'residual': residual, 'unit': unit,
                       'absolute_tolerance': tolerance,
                       'status': 'unobserved' if residual is None else 'passed'})
        if residual is not None and abs(residual) > tolerance:
            raise ValueError('Regional ownership/conservation failed: ' + label)

    steps = read('completed_native_steps', nonnegative=True)
    if steps != int(steps):
        raise ValueError('Nonintegral completed_native_steps')
    species = sorted({'Albumin', 'Glucose', 'Oxygen', 'CarbonDioxide',
                      'Sodium', 'Potassium', 'Chloride'} | {key.split('.species.', 1)[1].split('.')[0]
                      for key in raw if '.species.' in key})
    owners, views, paths, transfers, rates = {}, {}, [], [], {}
    masses_ug = {}
    for region, fraction in zip(REGIONS, FRACTIONS):
        observed_fraction = read(region + '.fraction', nonnegative=True)
        if abs(observed_fraction - fraction) > 1e-12:
            raise ValueError('Unsupported native installation fraction: ' + region)
        owner = 'Skin.' + region + '.extracellular'
        volume = read(region + '.volume_ml', nonnegative=True)
        masses_ug[region] = {s: read(region + '.species.' + s + '.mass_ug',
                                   required=False, nonnegative=True) for s in species}
        mass_g = {s: None if m is None else m * 1e-6
                  for s, m in masses_ug[region].items()}
        owners[owner] = {
            'native_owner': owner,
            'native_compartment': 'IHM_' + region + '_SkinTissueExtracellular',
            'accounting_owner': True, 'independent_store': True,
            'independent_chemical_law': False,
            'installation_fraction': fraction,
            'fraction_evidence': 'engineering_initial_volume_and_compliance_allocation',
            'volume_ml': volume, 'pressure_mmhg': read(region + '.pressure_mmhg'),
            'mass_g': mass_g,
            'concentration_g_per_l': dict.fromkeys(species),
            'mass_concentration_residual_g': dict.fromkeys(species),
            'ionic_molarity_mmol_per_l': dict.fromkeys(('Sodium', 'Potassium', 'Chloride')),
            'gas_partial_pressure_mmhg': dict.fromkeys(('Oxygen', 'CarbonDioxide')),
            'quantity_sources': {
                'volume_ml': {'source_key': PREFIX + region + '.volume_ml', 'source_unit': 'mL'},
                'pressure_mmhg': {'source_key': PREFIX + region + '.pressure_mmhg', 'source_unit': 'mmHg'},
                'mass_g': {s: {'source_key': PREFIX + region + '.species.' + s + '.mass_ug',
                              'source_unit': 'ug', 'unit': 'g', 'conversion_factor': 1e-6}
                           for s in species}},
        }
        flows = {}
        for path in PATHS:
            key = region + '.path.' + path + '.flow_ml_per_s'
            flows[path] = read(key)
            paths.append({'native_path': 'IHM_' + region + '_' + path,
                          'source_key': PREFIX + key, 'flow_ml_per_s': flows[path],
                          'role': 'compliance_boundary' if path == 'SkinE3ToGround'
                          else 'external_loss' if path == 'SkinSweating'
                          else 'native_circuit_path'})
        for path, source, target in (
                ('SkinE2ToSkinE3', 'Skin.vascular', owner),
                ('SkinE3ToSkinI', owner, 'Skin.intracellular'),
                ('SkinE3ToSkinL1', owner, 'Lymph'),
                ('SkinSweating', owner, 'external.sweat')):
            flow = flows[path]
            transfers.append({'native_path': 'IHM_' + region + '_' + path,
                              'source_key': PREFIX + region + '.path.' + path + '.flow_ml_per_s',
                              'source': source, 'target': target, 'flow_ml_per_s': flow,
                              'kind': 'observed_native_fluid_flow'})
            rates[source] = rates.get(source, 0.) - flow
            rates[target] = rates.get(target, 0.) + flow
        check(region + '.native_last_step_fluid',
              read(region + '.fluid_step_residual_ml', required=steps > 0))
    aggregate_volume = read('aggregate.volume_ml', nonnegative=True)
    check('children_minus_aggregate_volume',
          math.fsum(o['volume_ml'] for o in owners.values()) - aggregate_volume,
          aggregate_volume)
    for field in ('volume_ownership_residual_ml', 'install_volume_residual_ml'):
        check('aggregate.' + field, read('aggregate.' + field), aggregate_volume)
    check('aggregate.native_last_step_fluid',
          read('aggregate.fluid_step_residual_ml', required=steps > 0))
    aggregate_mass = {}
    for sub in species:
        total = read('aggregate.species.' + sub + '.mass_ug', False, True)
        children = [masses_ug[r][sub] for r in REGIONS]
        residual = None if total is None or None in children else math.fsum(children) - total
        check('children_minus_aggregate_mass.' + sub, residual, total or 0., 'ug')
        for field in ('ownership_residual_ug', 'install_residual_ug'):
            check('aggregate.' + sub + '.' + field,
                  read('aggregate.species.' + sub + '.' + field, False), total or 0., 'ug')
        aggregate_mass[sub] = None if total is None else total * 1e-6
        legacy = values.get('tissue.Skin.extracellular.' + sub + '.mass_g')
        if legacy is not None and total is not None:
            check('legacy_aggregate_mass.' + sub,
                  _number(legacy, sub, True) * 1e6 - total, total, 'ug')
    legacy_volume = values.get('tissue.Skin.extracellular.volume_ml')
    if legacy_volume is not None:
        check('legacy_aggregate_volume', _number(legacy_volume, 'legacy volume', True)
              - aggregate_volume, aggregate_volume)
    views['Skin.extracellular'] = {
        'native_owner': 'Skin.extracellular', 'native_compartment': 'SkinTissueExtracellular',
        'accounting_owner': False, 'independent_store': False,
        'classification': 'native_aggregate_view', 'children': list(owners),
        'volume_ml': aggregate_volume, 'mass_g': aggregate_mass}
    return {
        'schema': 'native_regional_skin_exchange_v1', 'available': True, 'time_s': time,
        'native_identity': deepcopy(native_identity or {}), 'source_hashes': dict(source_hashes or {}),
        'completed_native_steps': int(steps), 'native_compartments': owners,
        'aggregate_views': views, 'native_path_observations': paths,
        'native_circuit': deepcopy(raw), 'fluid_transfers': transfers,
        'internal_volume_rate_ml_per_s': {k: v for k, v in rates.items() if k != 'external.sweat'},
        'external_volume_rate_ml_per_s': {'external.sweat': rates['external.sweat']},
        'ownership_checks': checks, 'unobserved_source_keys': sorted(set(unknown)),
        'solute_fluxes': [], 'whole_body_mass_closure_claimed': False,
        'merge_contract': {'replace_accounting_owners': ['Skin.extracellular'],
                          'shared_owners': ['Skin.vascular', 'Skin.intracellular', 'Lymph'],
                          'omit_legacy_paths': ['SkinE1ToSkinE2', 'SkinE3ToSkinI', 'SkinE3ToSkinL1'],
                          'replace_parent_extracellular_partitions': True},
        'limitations': ['Regional extracellular stores share native aggregate Tissue/Diffusion chemical laws.',
                        'Installation fractions do not prescribe current regional quantities or measured anatomy.',
                        'Selected incidence uses E2ToE3 inflow and E3 outflows; serial paths are not extra transfers.',
                        'Compliance flow is not an additional drain. Sweat is an external boundary.',
                        'Concentrations, molarities, gas partial pressures and regional solute fluxes are unobserved.',
                        'Native last-step residuals do not establish whole-body or arbitrary-interval closure.']}
