"""Explicit engineered three-store configuration; no native mutation or anatomy claim."""
import gzip
import numpy as np
import hashlib
import json
import math
from pathlib import Path
import re

MATERIALIZATION_SHA = 'de669f96f6db5b150dd467f149464677e19f22ec2e7be5a77d895d873b3fef1a'


def digest_document(document):
    return hashlib.sha256(json.dumps(document,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def build_configuration(materialization_path, *, selected_unions, region_names, priors,
                        expected_sha256=MATERIALIZATION_SHA):
    raw=Path(materialization_path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected_sha256:
        raise ValueError('Changed engineered materialization')
    m=json.loads(raw)
    if len(region_names)!=3 or len(set(region_names))!=3 or any(not isinstance(x,str) or not re.fullmatch('[a-z][a-z0-9_]{0,79}',x) or x in ('region_a','region_b','residual') for x in region_names):
        raise ValueError('Three explicit unique semantic names required')
    if len(selected_unions)!=2 or any(not union or len(union)!=len(set(union)) for union in selected_unions) or set(selected_unions[0])&set(selected_unions[1]):
        raise ValueError('Two exclusive nonempty territory unions required')
    if len(priors)!=3:raise ValueError('Three explicit volume priors required')
    for prior in priors:
        if set(prior)!= {'thickness_m','extracellular_fraction','evidence'} or not isinstance(prior['evidence'],str) or not prior['evidence'].strip():
            raise ValueError('Incomplete engineering volume prior')
        for key in ['thickness_m','extracellular_fraction']:
            x=prior[key]
            if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<=0:
                raise ValueError('Invalid positive engineering volume prior')
        if prior['extracellular_fraction']>1:raise ValueError('EC fraction exceeds unity')
    regions={r['id']:r for r in m['regions']}
    if any(name not in regions for union in selected_unions for name in union):
        raise ValueError('Unknown selected geometric territory')
    eligible=m['contact_eligible_triangle_ids']; count=m['inventory']['triangle_count']
    if len(eligible)!=len(set(eligible)) or any(type(i) is not int or not 0<=i<count for i in eligible):
        raise ValueError('Invalid exterior inventory')
    mapping=[-1]*count
    for i in eligible:mapping[i]=2
    for index,union in enumerate(selected_unions):
        for rid in union:
            for face in regions[rid]['triangle_ids']:
                if mapping[face]!=2:raise ValueError('Overlapping or ineligible territory')
                mapping[face]=index
    areas={r['id']:r['area_m2'] for r in m['inventory']['regions']}
    a=math.fsum(areas[rid] for rid in selected_unions[0]); b=math.fsum(areas[rid] for rid in selected_unions[1])
    exterior_area=m['surface_diagnostic']['selected_area_m2']
    proxy_areas=[a,b,exterior_area-a-b]
    weights=[area*p['thickness_m']*p['extracellular_fraction'] for area,p in zip(proxy_areas,priors)]
    if any(not math.isfinite(w) or w<=0 for w in weights):raise ValueError('Nonpositive or nonfinite region volume proxy')
    total=math.fsum(weights); fractions=[weights[0]/total,weights[1]/total,0.]
    fractions[2]=1.-fractions[0]-fractions[1]
    if any(not math.isfinite(f) or f<=0 for f in fractions):raise ValueError('Degenerate normalized fraction')
    config={'schema':'engineered_native_skin_configuration_v1','materialization_sha256':expected_sha256,
            'skin_geometry_sha256':m['source_receipts'][0]['sha256'],
            'source_receipts':m['source_receipts'],'regions':[
                {'native_index':i,'name':name,'selected_territories':list(selected_unions[i]) if i<2 else None,
                 'selection':'explicit_union' if i<2 else 'entire_exterior_complement',
                 'contact_area_m2':proxy_areas[i],'volume_proxy_m3':weights[i],
                 'volume_prior':dict(priors[i]),'initial_fraction':fractions[i]}
                for i,name in enumerate(region_names)],
            'face_to_native_index':mapping,'contact_ineligible_index':-1,
            'source_face_count':count,'contact_eligible_face_count':len(eligible),
            'surface_evidence':'largest_connected_component_exterior_proxy_with_open_boundaries',
            'allocation_evidence':'normalize_engineered_area_times_thickness_times_EC_fraction',
            'full_native_inventory_preserved':True,
            'inventory_interpretation':'All existing Skin volume/species allocated once across three stores. Inner/seam exclusion removes contact representations only; no native tissue/fluid mass is discarded.',
            'measured_extracellular_volumes':False,'validated_lymph_watersheds':False,
            'native_activation_allowed':False,'native_commands':[]}
    config['configuration_sha256']=digest_document(config)
    return config


def validate_configuration(config):
    body={k:v for k,v in config.items() if k!='configuration_sha256'}
    if digest_document(body)!=config.get('configuration_sha256'):
        raise ValueError('Changed configuration receipt')
    fractions=[r['initial_fraction'] for r in config['regions']]
    if len(fractions)!=3 or any(not math.isfinite(f) or f<=0 for f in fractions) or abs(math.fsum(fractions)-1.)>1e-15:
        raise ValueError('Invalid configured fractions')
    return fractions


def map_contacts(config, contacts, *, root):
    """Return named pressure proposals with conserved supplied area/force, no commands."""
    validate_configuration(config)
    source=config['source_receipts'][0]
    raw=(Path(root)/source['path']).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=config['skin_geometry_sha256']:raise ValueError('Changed contact geometry')
    geometry=json.loads(gzip.decompress(raw))
    xyz=np.asarray(geometry['positions'],dtype=float).reshape(-1,3)
    faces=np.asarray(geometry['indices'],dtype=int).reshape(-1,3)
    rows={r['name']:{'contact_area_m2':0.,'normal_force_n':0.,'force_n':[0.,0.,0.]} for r in config['regions']}
    seen=set()
    for c in contacts:
        i=c['triangle_id']
        if type(i) is not int or not 0<=i<len(config['face_to_native_index']) or i in seen:
            raise ValueError('Invalid or repeated contact face')
        seen.add(i); region=config['face_to_native_index'][i]
        if region<0:raise ValueError('Contact on excluded inner/seam surface')
        area=c['area_m2']; pressure=c['normal_pressure_pa']; force=c['force_n']
        if len(force)!=3 or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in [area,pressure,*force]) or area<=0 or pressure<0:
            raise ValueError('Invalid contact observations')
        triangle=xyz[faces[i]]
        source_area=float(np.linalg.norm(np.cross(triangle[1]-triangle[0],triangle[2]-triangle[0]))/2)
        if area>source_area*(1+1e-12):raise ValueError('Contact exceeds pinned source face')
        row=rows[config['regions'][region]['name']]
        row['contact_area_m2']+=area; row['normal_force_n']+=area*pressure
        row['force_n']=[x+y for x,y in zip(row['force_n'],force)]
    for r in config['regions']:
        row=rows[r['name']]
        if row['contact_area_m2']>r['contact_area_m2']*(1+1e-12):raise ValueError('Contact exceeds regional surface')
        row['surface_pressure_proposal_pa']=row['normal_force_n']/r['contact_area_m2']
    return {'configuration_sha256':config['configuration_sha256'],'regions':rows,
            'interstitial_pressure_identified':False,'native_commands':[]}
