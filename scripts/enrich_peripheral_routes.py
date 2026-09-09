"""Versioned geometry contract and a preserved IBM fibre-velocity snapshot."""
import ast
import hashlib
import json
import math
from pathlib import Path
from peripheral_route_catalog import SOMATIC, VISCERAL, SPECIAL

LIMITATIONS = [
 'Authored schematic geometry, not dissected nerves or measured axon lengths.',
 'One representative endpoint-to-relay branch; not every branch of a nerve or plexus.',
 'Length excludes downstream cortical transport and synaptic delays.',
]

def enrich(data, lines, root):
    source = root.parent/'IBM-1/ibm/topologies/nerve.py'
    raw = source.read_bytes()
    tree = ast.parse(raw)
    velocities = next(ast.literal_eval(n.value) for n in tree.body
                      if isinstance(n, ast.AnnAssign) and n.target.id == 'FIBRE_VELOCITY_M_S')
    snapshot = root/'data/derived/canonical/peripheral-sources/ibm-nerve.py'
    snapshot.write_bytes(raw)
    data['source_receipts'].append(dict(source_path=str(source),preserved_path=str(snapshot.relative_to(root)),
        sha256=hashlib.sha256(raw).hexdigest(),role='fibre_velocity_prior_snapshot_not_executed'))
    data['fibre_velocity_m_s'] = {k:dict(low=v[0],typical=v[1],high=v[2]) for k,v in velocities.items()}
    data['fibre_velocity_scope'] = 'IBM nerve.py snapshot; unvalidated physiology priors, not independent IHM measurements. Composition is owned by IBM; local channels do not assert whole-trunk composition.'
    relays = {r['id']:r for r in data['relays']}
    nerves = {n['id']:n for n in data['nerves']}
    records = data['muscle_bindings']+data['receptor_patches']
    for r in records:
        r.update(geometry_kind='schematic_route',measured_axon_geometry=False,
            path_length_scope='representative_endpoint_to_relay',
            length_method='max(0.03, 1.15 * Euclidean distance(endpoint, relay))',
            endpoint_id=r.get('canonical_entity_id',r.get('id')),limitations=list(LIMITATIONS))
        if 'muscle_id' in r:
            r['delays_s']={c:r['path_length_m']/velocities[c][1] for c in ('ia','ib','ii','alpha','gamma')}
            r['proprioceptor_scope']='Single illustrative Ia-timed stretch/force proxy; separate Ib, II and gamma dynamics are not implemented.'
            r['central_delay_s']=.012
            r['afferent_delay_s']=r['delays_s']['ia']+.012
            r['motor_delay_s']=r['delays_s']['alpha']+.012
            r['scalar_delay_scope']='Compatibility aliases: Ia afferent and alpha efferent, each plus central_delay_s; never all fibre classes.'
    for n in nerves.values():
        candidates=sorted((r for r in records if r['nerve_id']==n['id']),key=lambda r:(r['path_length_m'],r.get('muscle_id',r.get('id'))))
        representative=candidates[len(candidates)//2]
        n.update(path_length_m=representative['path_length_m'],path_length_scope='representative_endpoint_to_relay',
            representative_binding_id=representative.get('muscle_id',representative.get('id')),
            endpoint_id=representative['endpoint_id'],length_method='upper median of existing endpoint-to-relay binding estimates',
            limitations=list(LIMITATIONS),runtime_support='existing_somatic_bindings')
    def add(side, name, level, endpoint, pos, relay_pos=None, cortex=None, kind='somatic'):
        if f'peripheral-nerve-{side}-{name}' in nerves:
            return  # Existing bound routes retain their identity and length.
        sign=1 if side=='left' else -1
        pos=[sign*pos[0],*pos[1:]]
        rid=f'peripheral-relay-{side}-{level}'
        if rid not in relays:
            xyz=relay_pos or [.012,.655,-.035]
            relays[rid]=dict(id=rid,name=f'{side} {level.replace("_"," ")} relay',side=side,
                position_m=[sign*xyz[0],*xyz[1:]],kind='named_relay_group',evidence_kind='authored_anatomical_prior')
        end=relays[rid]['position_m']
        # Somatic/visceral via an authored proximal waypoint; special senses direct.
        points=[pos,end] if kind=='special_sense' else [pos,[sign*.04,end[1],end[2]],end]
        length=sum(math.dist(a,b) for a,b in zip(points,points[1:]))
        nid=f'peripheral-nerve-{side}-{name}'
        n=dict(id=nid,name=f'{side} {name.replace("_"," ")} nerve',side=side,relay_id=rid,
            evidence_kind='authored_anatomical_prior',geometry_kind='schematic_route',measured_axon_geometry=False,
            path_length_m=length,path_length_scope='representative_endpoint_to_relay',
            endpoint_label=endpoint,points_m=points,length_method='sum of authored polyline segment lengths',
            route_kind=kind,runtime_support='topology_only',limitations=list(LIMITATIONS))
        if kind=='visceral':
            n['limitations'].append('Representative visceral afferent path only; efferent ganglionic stages and organ-specific branches are unresolved. Do not apply one velocity across a multi-neuron autonomic chain.')
        if cortex:
            n['downstream_brain_target_ids']=[f'brain-{h}-{cortex}' for h in ('lh','rh')]
            n['limitations'].append('Cortical targets are downstream annotations only; bilateral projections, decussation and central conduction are not resolved.')
        nerves[nid]=n
        lines.append(dict(id=f'peripheral-route-{side}-{name}',nerve_id=nid,entity_ids=[],points_m=points,
                          kind=kind,evidence_kind='schematic_anatomical_prior'))
    for side in ('left','right'):
        for name,level,endpoint,pos in SOMATIC: add(side,name,level,endpoint,pos)
        for name,level,endpoint,pos in VISCERAL: add(side,name,level,endpoint,pos,kind='visceral')
        for name,level,endpoint,pos,rpos,cortex in SPECIAL: add(side,name,level,endpoint,pos,rpos,cortex,'special_sense')
    data['nerves']=list(nerves.values());data['relays']=list(relays.values())
    data['counts'].update(nerves=len(nerves),relays=len(relays))
    data['units'].update(path_length='m',delay='s',velocity='m/s')
    data['route_contract']=dict(version=2,join_key='strip peripheral-nerve-{side}- from id',
        aliases={'sciatic_tibial':'sciatic'},unmatched_ihm_names=['sciatic_fibular'],
        length_field='nerves[].path_length_m',length_scope='representative_endpoint_to_relay',
        cortical_transport_included=False,synaptic_delays_included=False,
        missing_length_policy='error; never silently substitute a trunk length',
        new_routes_runtime_support='topology_only; no special-sense transduction or autonomic controller')
    data['limitations'] += LIMITATIONS + ['New routes are topology-only declarations; canonical cortical targets are downstream annotations, not executable special-sense inputs.']
    data['sources']['route_topology']=dict(url='https://www.ncbi.nlm.nih.gov/books/NBK537359/',finding='CN VIII brainstem relays and downstream auditory cortex; no coordinates or lengths transferred.')
    data['sources']['visceral_topology']=dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC9720663/',finding='Pelvic dissections support distinction of splanchnic pathways; no measured geometry transferred.')
