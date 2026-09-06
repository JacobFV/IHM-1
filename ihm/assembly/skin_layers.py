"""Explicit exterior-proxy support for synthesized skin layer volumes.

The source asset remains intact. A reviewed component mask supplies engineering
quadrature, not measured body surface area or a proof of an exclusive surface.
"""
import gzip
import hashlib
import json
import math
import numpy as np


def physical_skin_support(reference, geometry_bytes, evidence_bytes):
    digest=hashlib.sha256(geometry_bytes).hexdigest()
    if reference.get('sha256')!=digest or reference.get('units')!='m' or reference.get('frame')!='bodyparts3d-display-m':
        raise ValueError('Physical skin support requires exact source geometry in canonical metres')
    evidence=json.loads(evidence_bytes)
    receipts=[r for r in evidence.get('source_receipts',[]) if r.get('path')==reference.get('path') and r.get('sha256')==digest]
    if evidence.get('schema')!='engineered_skin_territories_v1' or len(receipts)!=1:
        raise ValueError('Matching reviewed skin component evidence required')
    geometry=json.loads(gzip.decompress(geometry_bytes))
    vertices=np.asarray(geometry['positions'],float).reshape(-1,3)
    source_faces=np.asarray(geometry['indices'])
    if source_faces.dtype.kind not in 'iu':raise ValueError('Integer source faces required')
    faces=source_faces.reshape(-1,3)
    if not len(faces) or not np.isfinite(vertices).all() or faces.min()<0 or faces.max()>=len(vertices):
        raise ValueError('Invalid skin source geometry')
    diagnostic=evidence['surface_diagnostic'];ids=evidence['contact_eligible_triangle_ids']
    if not isinstance(ids,list) or not ids or any(type(i) is not int or i<0 or i>=len(faces) for i in ids) or ids!=sorted(set(ids)):
        raise ValueError('Unique ordered physical skin face indices required')
    components=diagnostic['face_component_ids'];selected=diagnostic['selected_component_id']
    if type(selected) is not int or len(components)!=len(faces) or any(type(c) is not int or c<0 for c in components):
        raise ValueError('Source face component identity required')
    if ids!=[i for i,c in enumerate(components) if c==selected]:
        raise ValueError('Physical skin mask must match declared whole component')
    if diagnostic.get('selection_evidence')!='inferred_largest_component_positive_signed_volume_integral':
        raise ValueError('Unsupported physical skin selection basis')
    with np.errstate(over='ignore',invalid='ignore'):
        area=np.linalg.norm(np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]]),axis=1)/2
        total=float(area.sum())
    if not np.isfinite(area).all() or not math.isfinite(total):
        raise ValueError('Nonfinite source triangle area')
    chosen=float(area[ids].sum());declared=diagnostic.get('selected_area_m2')
    if isinstance(declared,bool) or not isinstance(declared,(int,float)) or not math.isfinite(declared) or chosen<=0 or not math.isclose(chosen,declared,rel_tol=1e-12,abs_tol=1e-12):
        raise ValueError('Physical skin support area disagrees with exact selected faces')
    # These source diagnostics are unresolved. Never promote the prior to an
    # independently validated physical surface merely by editing a flag.
    for key in ('physical_surface_exclusivity_validated','self_intersections_tested','signed_integral_is_closed_volume'):
        if diagnostic.get(key) is not False:raise ValueError('Unsupported physical skin validation claim')
    return {'schema':'ihm.physical-skin-support.v1','area_m2':chosen,
            'raw_source_area_m2':total,'selected_triangle_count':len(ids),
            'source_triangle_count':len(faces),'selected_component_id':selected,
            'source_geometry_sha256':digest,'component_evidence_sha256':hashlib.sha256(evidence_bytes).hexdigest(),
            'area_basis':'Explicit inferred exterior component; uniform inward layer depth is a separate prior',
            'clinical_body_surface_area':False,'physical_surface_exclusivity_validated':False,
            'self_intersections_tested':False,'signed_integral_is_closed_volume':False}
