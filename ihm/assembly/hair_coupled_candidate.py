"""One opt-in guide owner for the frozen reduced native hair candidate."""
import copy
import numpy as np
from .hair_dynamics import ElasticHairState
from .hair_feedback import HairFeedback


class HairReferenceRegistration:
    def __init__(self,base,registration):
        self.records={b['owner']:b for b in base['actual_frozen_native_inertia']['bodies']};self.bodies=list(self.records)
        self.global_map=np.asarray(registration['global_rigid_fit']['source_to_canonical_ground']);self.basis=self.global_map[:3,:3]
    def transforms(self,state):
        return {name:self.global_map@np.asarray(state['bodies'][name]['transform_ground'])@np.asarray(b['canonical_reference_to_body_local']) for name,b in self.records.items()}


def coupled_candidate(partition,base,registration,snapshot,*,additional_contact_faces=()):
    reg=HairReferenceRegistration(base,registration);nodes={};faces={};owners={};positions=[];offsets=[0];radii=[];follicles=[];beta=[];guide_owners=[];samples=[];material=None
    evidence={r['source_node_index']:r['owner'] for r in partition['ownership']['node_evidence']}
    for name,group in partition['groups'].items():
        patch=group['source_patch'];prepared=group['prepared'];strands=prepared['strands'];x=np.asarray(strands['centerlines_m']).reshape(-1,3)
        current={k:strands[k] for k in ('density_kg_m3','tensile_modulus_pa','bending_modulus_pa')}
        if material is not None and material!=current:raise ValueError('Combined owner requires identical declared materials')
        material=current
        for source,point in zip(patch['source_node_indices'],patch['positions_m']):
            if source in nodes and nodes[source]!=point:raise ValueError('Shared exact source node mismatch')
            nodes[source]=point;owners[source]=reg.bodies.index(evidence[source])
        for source,tri in zip(patch['source_face_indices'],patch['triangles']):
            mapped=[patch['source_node_indices'][i] for i in tri]
            if source in faces and faces[source]!=mapped:raise ValueError('Shared exact face mismatch')
            faces[source]=mapped
        for i,(a,b) in enumerate(zip(strands['strand_offsets'][:-1],strands['strand_offsets'][1:])):
            positions.extend(x[a:b].tolist());offsets.append(len(positions));radii.append(strands['radius_m'][i]);follicles.append(prepared['source_face_indices'][i]);beta.append(prepared['barycentric'][i]);guide_owners.append(prepared['guide_owner_names'][i]);samples.append((name,prepared['sample_ids'][i]))
    if len(additional_contact_faces)>8:raise ValueError('At most8explicit contact faces')
    for face in additional_contact_faces:
        body=face['owner']
        if body not in reg.bodies or len(face['source_node_indices'])!=3:raise ValueError('Exact registered contact triangle required')
        for source,point in zip(face['source_node_indices'],face['positions_m']):
            if source in nodes and (nodes[source]!=point or owners[source]!=reg.bodies.index(body)):raise ValueError('Additional face changes exact source node or owner')
            nodes[source]=point;owners[source]=reg.bodies.index(body)
        source=face['source_face_index']
        if source in faces and faces[source]!=face['source_node_indices']:raise ValueError('Additional face changes exact topology')
        faces[source]=face['source_node_indices']
    ids=sorted(nodes);index={n:i for i,n in enumerate(ids)};faceids=sorted(faces);tri=np.array([[index[n] for n in faces[f]] for f in faceids]);surface=np.array([nodes[n] for n in ids])
    hair=ElasticHairState({**material,'centerlines_m':positions,'strand_offsets':offsets,'radius_m':radii},max_substep_s=1e-4)
    feedback=HairFeedback(hair,surface,tri,np.array([owners[n] for n in ids]),reg,np.array([faceids.index(f) for f in follicles]),beta,
        identity={'source_id':next(iter(partition['groups'].values()))['source_patch']['source_id'],'sample_references':samples,'status':'Inferred reduced-body exterior guide candidate; no live enablement'},friction_static=0.,friction_kinetic=0.)
    transforms=reg.transforms(snapshot)
    for i,(a,b) in enumerate(zip(offsets[:-1],offsets[1:])):
        body=guide_owners[i];transform=transforms[body];hair.position_m[a:b]=hair.position_m[a:b]@transform[:3,:3].T+transform[:3,3]
        state=snapshot['bodies'][body];native_origin=np.asarray(state['transform_ground'])[:3,3];canonical_origin=reg.basis@native_origin+reg.global_map[:3,3]
        hair.velocity_m_s[a:b]=reg.basis@np.asarray(state['origin_velocity_m_s'])+np.cross(reg.basis@np.asarray(state['angular_velocity_rad_s']),hair.position_m[a:b]-canonical_origin)
    hair.root_m,hair.tangent=feedback.roots(snapshot);hair.time_s=snapshot['time_s']
    if feedback.clamp_residual(snapshot)>1e-10:raise ValueError('Initial exact material clamp mismatch')
    return feedback
