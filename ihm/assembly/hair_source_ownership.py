"""Explicit engineered hair ownership and non-applying native inertia partition.

Exact source nodes/faces never move. Nearest named bone-envelope ownership is
the existing neutral contact prior, with a disclosed retained-component rigid
clamp coherence rule. It is not anatomical compartment calibration.
"""
import copy,math,re
import numpy as np
from .hair_source_factory import guide_mass_properties,transform_mass_properties,partition_native_inertia


def classify_root_eligibility(mapping,face_component_ids,contact_eligible_ids,*,selected_component_id):
    labels=np.asarray(face_component_ids);faces=np.asarray(mapping['source_face_indices']);samples=mapping['sample_ids']
    if labels.ndim!=1 or labels.dtype.kind not in 'iu' or faces.ndim!=1 or faces.dtype.kind not in 'iu' or len(faces)!=len(samples) or not 0<len(faces)<=512 or np.any(faces<0) or faces.max()>=len(labels):raise ValueError('Exact bounded retained root/source component mapping required')
    permitted=set(map(int,contact_eligible_ids));eligible=[];rejected=[]
    for i,(face,sample) in enumerate(zip(faces,samples)):
        component=int(labels[face])
        if component==selected_component_id and int(face) in permitted:eligible.append(i)
        else:rejected.append({'guide_index':i,'sample_id':sample,'source_face_index':int(face),'component_id':component,'reason':'Root is on an excluded inner/seam source component or ineligible exterior face; no reattachment'})
    return {'original_guide_count':len(faces),'eligible_guide_indices':eligible,'ineligible_guides':rejected,'selection_status':'Engineering exterior proxy eligibility, not closed-volume or calibrated anatomical proof'}


def select_guide_geometry(geometry,indices):
    selected=np.asarray(indices);strands=geometry['strands'];offsets=np.asarray(strands['strand_offsets']);x=np.asarray(strands['centerlines_m']).reshape(-1,3);attachment=geometry['attachment'];count=len(strands['radius_m'])
    if selected.dtype.kind not in 'iu' or selected.ndim!=1 or not len(selected) or np.any(selected<0) or selected.max()>=count or len(np.unique(selected))!=len(selected):raise ValueError('Explicit nonempty retained guide selection required')
    centers=[];next_offsets=[0]
    for i in selected:centers.extend(x[offsets[i]:offsets[i+1]].tolist());next_offsets.append(len(centers))
    triangles=np.asarray(attachment['reference_triangles_m']).reshape(count,9)
    return {'strands':{**copy.deepcopy(strands),'centerlines_m':centers,'strand_offsets':next_offsets,'radius_m':np.asarray(strands['radius_m'])[selected].tolist()},
            'attachment':{**copy.deepcopy(attachment),'sample_ids':np.asarray(attachment['sample_ids'])[selected].tolist(),'barycentric':np.asarray(attachment['barycentric'])[selected].tolist(),'reference_triangles_m':triangles[selected].ravel().tolist()}}


def select_source_faces(patch,source_faces):
    wanted=set(map(int,source_faces));face_indices=np.array([i for i,f in enumerate(patch['source_face_indices']) if f in wanted]);tri=np.asarray(patch['triangles'])[face_indices]
    if not len(face_indices) or len(face_indices)!=len(wanted):raise ValueError('Eligible exact face is absent from retained patch')
    nodes=np.unique(tri)
    return {**copy.deepcopy(patch),'source_face_indices':np.asarray(patch['source_face_indices'])[face_indices].tolist(),'source_node_indices':np.asarray(patch['source_node_indices'])[nodes].tolist(),
            'positions_m':np.asarray(patch['positions_m'])[nodes].tolist(),'triangles':np.searchsorted(nodes,tri).tolist(),'coverage':'Exact exterior-eligible retained follicle faces only; excluded original roots remain separately recorded'}


def infer_source_owners(patches,registration,*,registration_sha256,ambiguity_margin_m=.001):
    names=list(registration['groups'])
    if not names or len(names)>128 or not re.fullmatch('[0-9a-f]{64}',registration_sha256) or not np.isfinite(ambiguity_margin_m) or ambiguity_margin_m<0:raise ValueError('Bounded registered owner set and ambiguity margin required')
    low=np.array([registration['groups'][name]['bounds_min_m'] for name in names],float);high=np.array([registration['groups'][name]['bounds_max_m'] for name in names],float)
    if low.shape!=(len(names),3) or high.shape!=low.shape or not np.isfinite(low).all() or not np.isfinite(high).all() or np.any(low>high):raise ValueError('Finite named bone envelopes required')
    coordinates={};faces={};source=None
    for patch in patches.values():
        identity=(patch['source_id'],patch['source_sha256'])
        if source is not None and source!=identity:raise ValueError('Patches have different source identities')
        source=identity;nodes=np.asarray(patch['source_node_indices']);x=np.asarray(patch['positions_m'],float);tri=np.asarray(patch['triangles']);face_ids=np.asarray(patch['source_face_indices'])
        if nodes.dtype.kind not in 'iu' or face_ids.dtype.kind not in 'iu' or len(np.unique(nodes))!=len(nodes) or len(np.unique(face_ids))!=len(face_ids) or x.shape!=(len(nodes),3) or not np.isfinite(x).all() or tri.dtype.kind not in 'iu' or tri.shape!=(len(face_ids),3) or not len(tri) or tri.min()<0 or tri.max()>=len(nodes):raise ValueError('Invalid exact source patch')
        for node,point in zip(nodes,x):
            key=int(node)
            if key in coordinates and not np.array_equal(coordinates[key],point):raise ValueError('Shared source node has incompatible coordinates')
            coordinates[key]=point.copy()
        for face,indices in zip(face_ids,tri):
            key=int(face);global_nodes=tuple(map(int,nodes[indices]))
            if key in faces and faces[key]!=global_nodes:raise ValueError('Shared source face has incompatible topology')
            faces[key]=global_nodes
    if not coordinates or len(coordinates)>12288 or len(faces)>4096:raise ValueError('Retained ownership patch exceeds budget')
    ids=sorted(coordinates);index={node:i for i,node in enumerate(ids)};x=np.array([coordinates[node] for node in ids]);parent={node:node for node in ids}
    def find(node):
        while parent[node]!=node:parent[node]=parent[parent[node]];node=parent[node]
        return node
    for nodes in faces.values():
        for a,b in zip(nodes,nodes[1:]):
            ra,rb=find(a),find(b);parent[max(ra,rb)]=min(ra,rb)
    components={}
    for node in ids:components.setdefault(find(node),[]).append(node)
    distance=np.linalg.norm(np.maximum(np.maximum(low[None]-x[:,None],x[:,None]-high[None]),0),axis=2);rank=np.argsort(distance,axis=1,kind='stable');nearest=rank[:,0];assigned=nearest.copy();component_records=[]
    for component,nodes in components.items():
        indices=np.array([index[node] for node in nodes]);split=len(np.unique(nearest[indices]))>1;cost=np.sum(distance[indices]**2,axis=0);order=np.argsort(cost,kind='stable')
        owner=int(order[0]) if split else int(nearest[indices[0]]);assigned[indices]=owner
        component_records.append({'component_source_node':component,'source_nodes':nodes,'owner':names[owner],'split_nearest_owners':split,
                                  'candidate_costs_m2':[{'owner':names[i],'sum_squared_distance_m2':float(cost[i])} for i in order[:2]],
                                  'rule':'Minimum summed squared envelope distance across retained connected source nodes' if split else 'Unmodified common nearest-envelope owner'})
    rows=[]
    for i,node in enumerate(ids):
        second=int(rank[i,1]) if len(names)>1 else int(rank[i,0]);margin=float(distance[i,second]-distance[i,nearest[i]])
        rows.append({'source_node_index':node,'owner':names[assigned[i]],'nearest_owner':names[nearest[i]],'nearest_distance_m':float(distance[i,nearest[i]]),
                     'second_nearest_owner':names[second],'second_nearest_distance_m':float(distance[i,second]),'distance_margin_m':margin,
                     'ambiguous_within_margin':len(names)>1 and margin<=ambiguity_margin_m,'coherence_override':bool(assigned[i]!=nearest[i]),'component_source_node':find(node)})
    by_node={r['source_node_index']:r for r in rows};basis='Engineered fixed nearest named bone-envelope prior; only split retained-face components receive explicit rigid-clamp coherence by minimum summed squared distance. Not anatomical calibration.'
    bindings={name:{'source_id':patch['source_id'],'source_sha256':patch['source_sha256'],'registration_sha256':registration_sha256,'basis':basis,
                    'nodes':[copy.deepcopy(by_node[int(node)]) for node in patch['source_node_indices']]} for name,patch in patches.items()}
    return {'native_enabled':False,'bindings':bindings,'node_evidence':rows,'components':component_records,'ambiguity_margin_m':ambiguity_margin_m,
            'ambiguity_status':'Engineering distance margin, not a confidence probability','summary':{'nodes':len(rows),'faces':len(faces),'coherent_components':len(components),'resolved_split_components':sum(c['split_nearest_owners'] for c in component_records),
            'overridden_nodes':sum(r['coherence_override'] for r in rows),'ambiguous_nodes':sum(r['ambiguous_within_margin'] for r in rows)},
            'limits':['Ownership domain is the current reduced named-body topology only','Unselected adjacent source faces and whole-skin boundary continuity are not validated','Source positions and exact source face/node correspondence remain unchanged']}


def _moments(mass,first,second,*,frame,owner=None):
    center=first/mass;central=second-mass*np.outer(center,center)
    return {'mass_kg':float(mass),'first_moment_kg_m':first.tolist(),'second_moment_origin_kg_m2':second.tolist(),'centroid_m':center.tolist(),
            'inertia_com_kg_m2':(np.trace(central)*np.eye(3)-central).tolist(),'inertia_origin_kg_m2':(np.trace(second)*np.eye(3)-second).tolist(),'frame':frame,'owner':owner}


def _combine(items,*,frame,owner=None):
    return _moments(math.fsum(p['mass_kg'] for p in items),np.sum([p['first_moment_kg_m'] for p in items],axis=0),np.sum([p['second_moment_origin_kg_m2'] for p in items],axis=0),frame=frame,owner=owner)


def partition_represented_guides(prepared,native_bodies):
    names=[b['owner'] for b in native_bodies]
    if not names or len(set(names))!=len(names):raise ValueError('Unique actual native owners required')
    by_owner={name:[] for name in names};references={name:[] for name in names};all_guides=[];seen=set()
    for group,data in prepared.items():
        strands=data['strands'];x=np.asarray(strands['centerlines_m'],float).reshape(-1,3);offsets=np.asarray(strands['strand_offsets']);radii=strands['radius_m'];owners=data['guide_owner_names'];samples=data['sample_ids']
        if len(owners)!=len(radii) or len(samples)!=len(radii) or len(offsets)!=len(radii)+1:raise ValueError('Complete retained guide ownership required')
        for i,(a,b) in enumerate(zip(offsets[:-1],offsets[1:])):
            owner=owners[i];key=(group,samples[i])
            if owner not in by_owner or key in seen:raise ValueError('Unknown owner or duplicate represented guide')
            seen.add(key);single={**strands,'centerlines_m':x[a:b].tolist(),'strand_offsets':[0,b-a],'radius_m':[radii[i]]};properties=guide_mass_properties(single,owner=owner)
            by_owner[owner].append(properties);all_guides.append(properties);references[owner].append({'group':group,'guide_index':i,'sample_id':samples[i],'mass_kg':properties['mass_kg']})
    if not all_guides:raise ValueError('At least one represented independent guide required')
    before_global=[];after_global=[];partitions=[]
    for body in native_bodies:
        owner=body['owner'];mass=body['mass_kg'];center=np.asarray(body['centroid_m']);inertia=np.asarray(body['inertia_com_kg_m2']);matrix=np.asarray(body['canonical_reference_to_body_local'])
        if not np.isfinite(mass) or mass<=0 or center.shape!=(3,) or inertia.shape!=(3,3) or not np.isfinite(center).all() or not np.isfinite(inertia).all() or not np.allclose(inertia,inertia.T,rtol=0,atol=1e-14):raise ValueError('Finite actual native spatial inertia required')
        central=.5*np.trace(inertia)*np.eye(3)-inertia
        if np.linalg.eigvalsh(central).min()<=0:raise ValueError('Actual native second moment is nonpositive')
        before=_moments(mass,mass*center,central+mass*np.outer(center,center),frame=body['frame'],owner=owner)
        if by_owner[owner]:
            canonical=_combine(by_owner[owner],frame='canonical-reference-m',owner=owner)
            selected=transform_mass_properties(canonical,matrix,frame=body['frame'],owner=owner)
            receipt=partition_native_inertia(body,selected,partition_basis='Explicit inferred partition of this actual reduced native body for retained independent guides only; not canonical proxy subtraction')
            after=receipt['remaining'];partitions.append({'owner':owner,'native_before':before,'represented_guides':selected,'native_after':after,
                'guide_references':references[owner],'admissibility':receipt['residual_second_moment_eigenvalues_kg_m2'],'native_enabled':False})
        else:after=before
        inverse=np.linalg.inv(matrix)
        before_global.append(transform_mass_properties(before,inverse,frame='canonical-reference-m',owner=owner))
        after_global.append(transform_mass_properties(after,inverse,frame='canonical-reference-m',owner=owner))
    before=_combine(before_global,frame='canonical-reference-m');after=_combine([*after_global,*all_guides],frame='canonical-reference-m')
    mass_residual=after['mass_kg']-before['mass_kg'];first=np.array(after['first_moment_kg_m'])-before['first_moment_kg_m'];second=np.array(after['second_moment_origin_kg_m2'])-before['second_moment_origin_kg_m2']
    if abs(mass_residual)>1e-12 or np.linalg.norm(first)>1e-11 or np.linalg.norm(second)>1e-11:raise ValueError('Represented-guide/native residual moment conservation failed')
    return {'native_enabled':False,'represented_guide_count':len(all_guides),'represented_mass_kg':math.fsum(g['mass_kg'] for g in all_guides),'modified_owner_count':len(partitions),'partitions':partitions,
            'conservation':{'before':before,'residual_native_plus_guides':after,'mass_residual_kg':mass_residual,'first_moment_residual_kg_m':first.tolist(),'second_moment_residual_kg_m2':second.tolist(),
                            'inertia_origin_residual_kg_m2':(np.trace(second)*np.eye(3)-second).tolist()},
            'limits':['This is an inferred preparation, not an installed native mass change','Unrepresented population remains in the lumped native bodies; guides have no population weights',
                      'Moving-state activation requires matched rigid initial guide velocities and momentum/energy receipts','A changed cervical/head topology invalidates this reduced-body partition']}
