"""State ownership for the executable reduced body, separate from display meshes."""
from dataclasses import asdict
from .contracts import StateSpec
from .interfaces import ContractRegistry


def body_state_registry(assets,body_version):
    entities={e['id']:e for e in assets['anatomy']['entities']}
    registry=ContractRegistry(body_version,entities)
    thoracic={b['entity_id'] for b in assets['respiration']['bindings']}
    affine={b['entity_id'] for b in assets['respiration']['bindings'] if b['kind'] in ['lung','diaphragm']}
    def state(id,owner,support,unit,space,domain):
        registry.register_state(StateSpec(id,owner,support,unit,None,space,'time-domain samples with explicit integrator clock',domain))
    for id,e in entities.items():
        state(id+'/translation','thoracic' if id in thoracic else 'mechanics',id,'m','centroid vector','Reduced reference-centroid displacement; constrained orientation')
        state(id+'/deformation','thoracic' if id in affine else 'mechanics',id,'1','affine 3x3 tensor','One affine tensor per entity; no volumetric contact')
    for node in assets['brain']['nodes']:
        for variable,unit in [('potential','V'),('activity','Hz'),('adaptation','1')]:
            state(node['id']+'/'+variable,'brain',node['body_entity_id'],unit,'population','Preserved IBM population reduction; export potential mV converts to V')
    skin=next(e['id'] for e in entities.values() if e['role']=='skin')
    for patch in assets['peripheral']['receptor_patches']:
        for modality in ['pressure','stretch','warm','cold']:
            state(patch['id']+'/'+modality,'peripheral',patch['body_entity_id'],'Hz','regional receptor patch','Excess-rate generic prior; not individual afferent population fit')
    for binding in assets['peripheral']['muscle_bindings']:
        state(binding['muscle_id']+'/activation','peripheral',binding['canonical_entity_id'],'1','actuator','Delayed explicit command; activation relaxation prior')
    state('peripheral/events','peripheral',skin,'1','timestamped event queue','Queue is checkpointed; signals have conduction delay')
    state('native/physiology','native',skin,'1','heterogeneous external source record','Native CSV units retained per named column; sole gas/blood/thermal storage owner; opaque replay boundary')
    # Replay transfers are observations/Dirichlet boundaries, not Exchanges.
    # Registering a fictitious energy exchange would conceal absent feedback.
    registry.validate_scenario({'thoracic','mechanics','brain','peripheral','native'},set(),set(),set())
    return registry


def describe_states(registry):
    return {'body_version':registry.body_version,'states':[asdict(s) for s in registry.states.values()],
        'transfer_semantics':'native recorded values and prescribed mechanical boundaries; no conservative native return exchange',
        'unresolved':['native per-column typed state expansion','empirical laws for cortical motor recruitment','reference contact and articulated joints','two-way native pressure/work exchange']}
