"""Exact retained hair material mapping and non-applying inertia preparation.

Canonical proxy masses are metadata, not separately integrated native masses.
Any debit must explicitly partition actual native segment inertia in one frame.
This module never changes a native body or enables a coupled session.
"""
import copy
from types import SimpleNamespace
import numpy as np


def map_retained_guides(attachment,population,*,source_id,source_sha256):
    if attachment.get('skin_entity_id')!=source_id or attachment.get('source_geometry_sha256')!=source_sha256:raise ValueError('Retained hair/source identity mismatch')
    ids=np.asarray(population['ids']);wanted=np.asarray(attachment['sample_ids']);faces=np.asarray(population['face_index'])
    if ids.ndim!=1 or ids.dtype.kind not in 'iu' or not 0<len(ids)<=2000000 or len(np.unique(ids))!=len(ids):raise ValueError('Population IDs are absent, ambiguous or outside budget')
    if wanted.ndim!=1 or wanted.dtype.kind not in 'iu' or not 0<len(wanted)<=512 or len(np.unique(wanted))!=len(wanted):raise ValueError('Retained guide IDs must be unique and bounded')
    if faces.shape!=ids.shape or faces.dtype.kind not in 'iu' or np.any(faces<0):raise ValueError('Retained population face indices required')
    order=np.argsort(ids);locations=np.searchsorted(ids[order],wanted)
    if np.any(locations>=len(ids)) or not np.array_equal(ids[order][locations],wanted):raise ValueError('Retained sample ID is absent from its original population')
    selected=order[locations]
    if 'source_geometry_sha256' in population and str(np.asarray(population['source_geometry_sha256']).item())!=source_sha256:raise ValueError('Population source hash mismatch')
    if 'source_entity_id' in population and str(np.asarray(population['source_entity_id']).item())!=source_id:raise ValueError('Population source entity mismatch')
    bary=np.asarray(attachment['barycentric'],float);population_bary=np.asarray(population['barycentric'],float)
    if population_bary.shape!=(len(ids),3) or bary.shape!=(len(wanted),3) or not np.isfinite(bary).all() or np.any(bary<0) or not np.allclose(bary.sum(1),1,rtol=0,atol=1e-12) or not np.allclose(bary,population_bary[selected],rtol=0,atol=1e-12):raise ValueError('Retained barycentric material weights differ from population')
    return {'sample_ids':wanted.tolist(),'population_indices':selected.tolist(),'source_face_indices':faces[selected].tolist()}


def validate_retained_source(geometry,population,source_patch):
    attachment=geometry['attachment'];strands=geometry['strands'];patch=source_patch
    if isinstance(patch,dict):patch=SimpleNamespace(**{k:np.asarray(v) if k in ('positions_m','triangles','source_node_indices','source_face_indices') else v for k,v in patch.items()})
    mapping=map_retained_guides(attachment,population,source_id=patch.source_id,source_sha256=patch.source_sha256)
    local_faces={int(face):i for i,face in enumerate(patch.source_face_indices)}
    if any(face not in local_faces for face in mapping['source_face_indices']):raise ValueError('Retained follicle face is absent from exact source patch')
    follicles=np.array([local_faces[face] for face in mapping['source_face_indices']]);triangles=patch.positions_m[patch.triangles[follicles]]
    expected=np.asarray(attachment['reference_triangles_m'],float)
    if expected.size!=triangles.size or not np.allclose(expected.reshape(triangles.shape),triangles,rtol=0,atol=1e-12):raise ValueError('Retained follicle triangle differs from exact source face')
    roots=np.sum(np.asarray(attachment['barycentric'])[:,:,None]*triangles,axis=1)
    normal=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]);length=np.linalg.norm(normal,axis=1)
    if np.any(length<=1e-16):raise ValueError('Degenerate retained follicle face')
    normal/=length[:,None];selected=np.array(mapping['population_indices']);count=len(selected)
    for key,expected in [('roots_m',roots),('normals',normal)]:
        values=np.asarray(population[key],float)
        if values.shape!=(len(population['ids']),3) or not np.allclose(values[selected],expected,rtol=0,atol=1e-10):raise ValueError('Population root or face-normal direction mismatch')
    x=np.asarray(strands['centerlines_m'],float).reshape(-1,3);offsets=np.asarray(strands['strand_offsets']);radii=np.asarray(strands['radius_m'],float)
    if offsets.dtype.kind not in 'iu' or offsets.shape!=(count+1,) or offsets[0]!=0 or offsets[-1]!=len(x) or np.any(np.diff(offsets)<3) or np.any(np.diff(offsets)>32) or radii.shape!=(count,):raise ValueError('Retained beam population/topology mismatch')
    lengths=np.asarray(population['length_m'],float)[selected];expected_radii=np.asarray(population['radius_m'],float)[selected]
    if not np.isfinite(radii).all() or np.any(radii<=0) or not np.allclose(radii,expected_radii,rtol=0,atol=1e-15) or not np.isfinite(lengths).all() or np.any(lengths<=0):raise ValueError('Retained physical radius/length mismatch')
    material=np.asarray([strands[k] for k in ('density_kg_m3','tensile_modulus_pa','bending_modulus_pa')],float)
    if not np.isfinite(material).all() or np.any(material<=0):raise ValueError('Positive finite beam material required')
    for s,(a,b) in enumerate(zip(offsets[:-1],offsets[1:])):
        expected=roots[s]+np.linspace(0,lengths[s],b-a)[:,None]*normal[s]
        if not np.allclose(x[a:b],expected,rtol=0,atol=1e-10):raise ValueError('Beam centerline differs from retained physical root/normal/length')
    return {**mapping,'follicle_triangles':follicles.tolist(),'source_node_indices':patch.source_node_indices.tolist(),
            'source_id':patch.source_id,'source_sha256':patch.source_sha256,'source_correspondence_validated':True}


def prepare_hair_factory(geometry,population,source_patch):
    validated=validate_retained_source(geometry,population,source_patch);patch=source_patch
    owners=patch.node_owner_indices[patch.triangles[validated['follicle_triangles']]]
    if np.any(owners!=owners[:,0,None]):raise ValueError('Cross-owner follicle clamp derivative is unresolved')
    return {**validated,'guide_owner_names':[patch.owner_names[i] for i in owners[:,0]],
            'registration_sha256':patch.registration_sha256,'owner_basis':patch.owner_basis,'strands':copy.deepcopy(geometry['strands']),'barycentric':copy.deepcopy(geometry['attachment']['barycentric']),
            'native_enabled':False,'blockers':['Actual native segment inertia partition/debit is not applied','Single-interval common-clock composition and authoritative viewer bridge remain required']}


def _properties(mass,first,second,*,frame,owner=None):
    if not np.isfinite(mass) or mass<=0 or first.shape!=(3,) or second.shape!=(3,3) or not np.isfinite(first).all() or not np.isfinite(second).all():raise ValueError('Finite positive mass and compatible moments required')
    center=first/mass;central=second-mass*np.outer(center,center)
    return {'mass_kg':float(mass),'first_moment_kg_m':first.tolist(),'second_moment_origin_kg_m2':second.tolist(),
            'centroid_m':center.tolist(),'inertia_origin_kg_m2':(np.trace(second)*np.eye(3)-second).tolist(),
            'inertia_com_kg_m2':(np.trace(central)*np.eye(3)-central).tolist(),'frame':frame,'owner':owner}


def guide_mass_properties(strands,*,owner=None,frame='canonical-reference-m'):
    x=np.asarray(strands['centerlines_m'],float).reshape(-1,3);offsets=np.asarray(strands['strand_offsets']);radii=np.asarray(strands['radius_m'],float);density=strands['density_kg_m3']
    if not 0<len(radii)<=512 or offsets.dtype.kind not in 'iu' or offsets.shape!=(len(radii)+1,) or offsets[0]!=0 or offsets[-1]!=len(x) or not np.isfinite(x).all() or not np.isfinite(radii).all() or np.any(radii<=0) or not np.isfinite(density) or density<=0:raise ValueError('Finite retained strand mass inputs required')
    masses=np.zeros(len(x))
    for s,(a,b) in enumerate(zip(offsets[:-1],offsets[1:])):
        if not 3<=b-a<=32:raise ValueError('Bounded guide nodes required')
        segments=np.linalg.norm(np.diff(x[a:b],axis=0),axis=1)
        if np.any(segments<=0) or not np.allclose(segments,segments[0],rtol=1e-6,atol=0):raise ValueError('Uniform guide segments required')
        masses[a:b]=density*np.pi*radii[s]**2*segments[0];masses[[a,b-1]]*=.5
    result=_properties(masses.sum(),np.sum(masses[:,None]*x,axis=0),x.T@(masses[:,None]*x),frame=frame,owner=owner)
    result['representation']='Exact lumped nodal inertia of retained independent beam guides; no bundle or render-population weighting'
    result['guide_count']=len(radii);return result


def population_mass_properties(population,density_kg_m3):
    roots=np.asarray(population['roots_m']);normal=np.asarray(population['normals']);radii=np.asarray(population['radius_m']);lengths=np.asarray(population['length_m']);mass=0.;first=np.zeros(3);second=np.zeros((3,3))
    if roots.ndim!=2 or roots.shape[1]!=3 or not 0<len(roots)<=2000000 or normal.shape!=roots.shape or radii.shape!=(len(roots),) or lengths.shape!=radii.shape or not np.isfinite(density_kg_m3) or density_kg_m3<=0:raise ValueError('Bounded compatible physical population inputs required')
    for a in range(0,len(roots),2048):
        b=min(a+2048,len(roots));r=radii[a:b];length=lengths[a:b];direction=normal[a:b];norm=np.linalg.norm(direction,axis=1)
        if np.any(norm<=0) or np.any(r<=0) or np.any(length<=0) or not all(np.isfinite(v).all() for v in (roots[a:b],direction,r,length)):raise ValueError('Finite physical population strands required')
        direction=direction/norm[:,None];m=density_kg_m3*np.pi*r*r*length;center=roots[a:b]+.5*length[:,None]*direction
        mass+=m.sum();first+=np.sum(m[:,None]*center,axis=0);second+=center.T@(m[:,None]*center)
        # Exact second moment of each straight circular cylinder, accumulated
        # in bounded batches. This population is a geometric/material prior.
        axial=m*length*length/12;radial=m*r*r/4
        second+=np.eye(3)*radial.sum()+direction.T@((axial-radial)[:,None]*direction)
    result=_properties(mass,first,second,frame='canonical-reference-m');result['strand_count']=len(roots)
    result['representation']='Full retained independent straight-cylinder population under stated density/length/radius priors; not guide weights or measured subject mass'
    return result


def transform_mass_properties(properties,transform,*,frame,owner):
    matrix=np.asarray(transform,float)
    if matrix.shape!=(4,4) or not np.isfinite(matrix).all() or not np.allclose(matrix[3],[0,0,0,1],rtol=0,atol=1e-12):raise ValueError('Finite rigid common-frame transform required')
    rotation=matrix[:3,:3];translation=matrix[:3,3]
    if not np.allclose(rotation.T@rotation,np.eye(3),rtol=0,atol=1e-10) or np.linalg.det(rotation)<0:raise ValueError('Inertia registration must be proper rigid, not a scale fit')
    mass=properties['mass_kg'];rotated_first=rotation@properties['first_moment_kg_m'];first=rotated_first+mass*translation
    second=rotation@properties['second_moment_origin_kg_m2']@rotation.T+np.outer(rotated_first,translation)+np.outer(translation,rotated_first)+mass*np.outer(translation,translation)
    result=_properties(mass,first,second,frame=frame,owner=owner);result['representation']=properties.get('representation');return result


def partition_native_inertia(native,selected,*,partition_basis):
    if not isinstance(partition_basis,str) or not partition_basis.strip() or not isinstance(native.get('frame'),str) or not native['frame'] or native.get('frame')!=selected.get('frame') or not isinstance(native.get('owner'),str) or not native['owner'] or native['owner']!=selected.get('owner'):raise ValueError('Explicit same-owner/common-frame inferred native inertia partition required')
    mass=float(native['mass_kg']);center=np.asarray(native['centroid_m'],float);inertia=np.asarray(native['inertia_com_kg_m2'],float)
    if center.shape!=(3,) or inertia.shape!=(3,3) or not np.isfinite(mass) or mass<=0 or not np.isfinite(center).all() or not np.isfinite(inertia).all() or not np.allclose(inertia,inertia.T,rtol=0,atol=1e-14):raise ValueError('Finite native segment spatial inertia required')
    native_second=.5*np.trace(inertia)*np.eye(3)-inertia+mass*np.outer(center,center)
    selected_mass=float(selected['mass_kg']);selected_first=np.asarray(selected['first_moment_kg_m'],float);selected_second=np.asarray(selected['second_moment_origin_kg_m2'],float)
    if not np.isfinite(selected_mass) or selected_mass<=0 or selected_first.shape!=(3,) or selected_second.shape!=(3,3) or not np.isfinite(selected_first).all() or not np.isfinite(selected_second).all() or not np.allclose(selected_second,selected_second.T,rtol=0,atol=1e-14):raise ValueError('Finite positive selected guide spatial inertia required')
    if np.linalg.eigvalsh(selected_second-np.outer(selected_first,selected_first)/selected_mass).min()<-1e-12*max(np.linalg.norm(selected_second),1e-30):raise ValueError('Selected guide inertia is nonphysical')
    remaining_mass=mass-selected_mass;remaining_first=mass*center-selected_first;remaining_second=native_second-selected_second
    if not np.isfinite(remaining_mass) or remaining_mass<=0:raise ValueError('Guide debit exhausts native segment mass')
    center_remaining=remaining_first/remaining_mass;central=remaining_second-remaining_mass*np.outer(center_remaining,center_remaining)
    tolerance=1e-12*max(float(np.linalg.norm(inertia)),1e-15)
    if not np.isfinite(central).all() or np.linalg.eigvalsh(central).min()<=tolerance:raise ValueError('Guide debit leaves a nonpositive residual native inertia tensor')
    remaining=_properties(remaining_mass,remaining_first,remaining_second,frame=native['frame'],owner=native['owner'])
    return {'native_enabled':False,'partition_basis':partition_basis,'status':'Algebraically admissible inferred partition only; no anatomical calibration or native mutation',
            'remaining':remaining,'selected':copy.deepcopy(selected),'residual_second_moment_eigenvalues_kg_m2':np.linalg.eigvalsh(central).tolist()}
