"""Strict server-owned controller choices; clients cannot supply model paths."""

# Kinds that own the identified upright native model and its fixed initial pose.
STANCE_KINDS = ('engineering_stance', 'implicit_cortical_stance', 'implicit_curriculum16_stance')
# Kinds whose cortex is a persistent trained IBM E/I stance policy at 10 ms.
CORTICAL_STANCE_KINDS = ('implicit_cortical_stance', 'implicit_curriculum16_stance')
# Kinds that run the raw kernel through the segmental cord with no trained motor
# policy. `None` keeps the historical donor default; a path is a retained kernel.
RAW_IMPLICIT_KERNELS = {'implicit': None,
                        'implicit_curriculum16': 'data/models/ibm_curriculum16_kernel_v1/kernel.pt'}
RAW_IMPLICIT_KINDS = tuple(RAW_IMPLICIT_KERNELS)
KINDS = ('regional', 'implicit', 'implicit_curriculum16', 'implicit_ankle_primitive',
         'implicit_cortical_ankle', 'engineering_stance', 'implicit_cortical_stance',
         'implicit_curriculum16_stance')


def resolve_controller(value=None):
    if value is None:
        return {'kind': 'regional', 'sever': False, 'no_cord': False}
    if not isinstance(value, dict) or set(value) - {'kind', 'sever', 'no_cord', 'target_rad'}:
        raise ValueError('Unknown controller configuration')
    kind = value.get('kind', 'regional')
    if kind not in KINDS:
        raise ValueError('Unknown controller kind')
    result = {'kind': kind}
    for name in ('sever', 'no_cord'):
        flag = value.get(name, False)
        if type(flag) is not bool:
            raise ValueError(f'Controller {name} must be boolean')
        result[name] = flag
    if kind in ('regional','engineering_stance') and (result['sever'] or result['no_cord']):
        raise ValueError('Kernel/cord ablations require the implicit controller')
    if kind in CORTICAL_STANCE_KINDS and result['no_cord']:
        raise ValueError('Cortical stance has no cord motor path to ablate')
    if kind in ('implicit_ankle_primitive', 'implicit_cortical_ankle'):
        import math
        target = value.get('target_rad', .12)
        if type(target) not in (float, int) or not math.isfinite(target) or not -.25 <= target <= .25:
            raise ValueError('Ankle target must be finite within [-.25,.25] radians')
        if kind == 'implicit_ankle_primitive' and result['no_cord']:
            raise ValueError('Ankle primitive uses its trained motor output; no-cord is not an ablation of that policy')
        result['target_rad'] = float(target)
    elif 'target_rad' in value:
        raise ValueError('Ankle target requires the trained ankle primitive')
    return result
