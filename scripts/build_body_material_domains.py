"""Materialize source-shaped pelvic solids and a whole-body ownership registry."""
from pathlib import Path
import argparse,hashlib,json,sys
import numpy as np
BASE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(BASE))
from ihm.assembly.material_domains import source_surface,voxel_partition,MaterialOwnership
from ihm.assembly.contact_dynamics import DynamicTetrahedra

SOURCES=['body-bp3d-FJ3132','body-bp3d-FJ3133','body-bp3d-FJ3134']
LITERATURE=[
 {'doi':'10.1016/j.actbio.2024.06.035','url':'https://pubmed.ncbi.nlm.nih.gov/38945188/','title':'Characterisation of human penile tissue properties using experimental testing combined with multi-target inverse finite element modelling','species':'human','preparation':'fresh-frozen; individual tissue and whole segment tests','status':'primary abstract inspected; numerical tables and constitutive parameters not acquired','use':'priority for calibration; no coefficients imported'},
 {'doi':'10.1016/j.actbio.2024.03.013','url':'https://pubmed.ncbi.nlm.nih.gov/38494081/','title':'Experimental testing combined with inverse-FE for mechanical characterisation of penile tissues','species':'horse','status':'primary abstract inspected; numerical tables not acquired','use':'testing-method evidence; no human coefficient substitution'},
 {'doi':'10.1016/j.compbiomed.2023.107524','url':'https://research.tudelft.nl/en/publications/development-of-in-silico-models-to-guide-the-experimental-charact/','title':'Development of in silico models to guide the experimental characterisation of penile tissue and inform surgical treatment of erectile dysfunction','species':'idealized human geometry','status':'primary abstract inspected; publisher/repository fulltext fetch not available in this run','use':'tunica/material identification gap; not quantitative calibration'},
 {'url':'https://pubmed.ncbi.nlm.nih.gov/33945173/','title':'Assessment of the Rigidity Changes of Corpus Cavernosum Penis in Vasculary Erectile Dysfunction (ED) Subtypes by Shear Wave Elastography (SWE)','species':'human','status':'primary abstract inspected','use':'measurement context; SWE apparent modulus not copied into a quasi-static large-strain constitutive law'}]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def build(spacing_m=.004,output_dir=None):
    out=Path(output_dir or BASE/f'data/derived/material-domains/pelvis-{spacing_m:g}m')
    if out.exists():raise ValueError('Choose a fresh material-domain directory')
    anatomy_path=BASE/'data/derived/canonical/anatomy.json';mechanics_path=BASE/'data/derived/canonical/mechanics.json'
    anatomy=json.loads(anatomy_path.read_text());mechanics=json.loads(mechanics_path.read_text())
    by_id={e['id']:e for e in anatomy['entities']};mech={e['id']:e for e in mechanics['entities']}
    surfaces=[];surface_metadata=[]
    for source_id in SOURCES:
        e=by_id[source_id];path=BASE/e['reference_geometry']['path']
        if sha(path)!=e['reference_geometry']['sha256']:raise ValueError('Stale canonical source surface')
        x,t,metadata=source_surface(path);surfaces.append((source_id,x,t));surface_metadata.append(dict(metadata,id=source_id,name=e['name'],path=str(path.relative_to(BASE))))
    partition=voxel_partition(surfaces,spacing_m=spacing_m)
    ownership=MaterialOwnership({s:mech[s]['mass_kg'] for s in SOURCES});ownership.claim('pelvic-erectile-tissue-volume',SOURCES)
    owner=partition['material_index'];volumes=np.bincount(owner,minlength=len(SOURCES))*spacing_m**3/6
    materials=[]
    for i,s in enumerate(SOURCES):
        material=mech[s]['material']
        materials.append({'source_id':s,'name':by_id[s]['name'],'volume_m3':float(volumes[i]),
            'allocated_mass_kg':mech[s]['mass_kg'],'effective_inertial_density_kg_m3':float(mech[s]['mass_kg']/volumes[i]),
            'constitutive':'compressible neo-Hookean','mu_pa':material['shear_modulus']['value'],
            'lambda_pa':material['lame_lambda']['value'],'source_material':material,
            'calibration_status':'existing generic canonical prior; penile tissue tables not acquired',
            'density_basis':'existing exclusive canonical mass allocation / synthesized voxel volume; not independently measured density'})
    mu=np.array([m['mu_pa'] for m in materials])[owner];lam=np.array([m['lambda_pa'] for m in materials])[owner]
    density=np.array([m['effective_inertial_density_kg_m3'] for m in materials])[owner]
    x=partition['vertices_m'];t=partition['tetrahedra']
    # Posterior/proximal surface nodes form an explicit support fixture. The
    # threshold is not an inferred dissected tendon or a claimed insertion site.
    body=DynamicTetrahedra(x,t,mu_pa=mu,lambda_pa=lam,density_kg_m3=density)
    fixed=body.surface_nodes[x[body.surface_nodes,2]<x[:,2].min()+1.01*spacing_m]
    registry=[]
    for e in anatomy['entities']:
        m=mech[e['id']]
        registry.append({'entity_id':e['id'],'name':e['name'],'system':e['system'],'role':e['role'],
            'material_owner':e['id'],'mass_kg':m['mass_kg'],'mass_role':m['mass_role'],
            'source_geometry':e['reference_geometry'],'evidence_kind':e['evidence_kind'],
            'material_prior':m['material'],'resolved_volume_domain':'pelvic-erectile-tissue-volume' if e['id'] in SOURCES else None,
            'runtime_status':'detached dynamic replacement available; existing owner must be deactivated before coupling' if e['id'] in SOURCES else 'existing reduced owner; volumetric dynamics and collision interfaces unresolved'})
    out.mkdir(parents=True)
    np.savez_compressed(out/'pelvic-domain.npz',vertices_m=x,tetrahedra=t,material_index=owner,mu_pa=mu,lambda_pa=lam,density_kg_m3=density,fixed_nodes=fixed,boundary_triangles=body.triangles)
    (out/'registry.json').write_text(json.dumps({'schema_version':1,'frame':anatomy['frame'],'entities':registry,'single_body_mass_kg':sum(e['mass_kg'] for e in registry)},separators=(',',':'))+'\n')
    receipt={'schema_version':1,'domain_id':'pelvic-erectile-tissue-volume','frame':anatomy['frame'],'source_surfaces':surface_metadata,
        'sources':{str(p.relative_to(BASE)):sha(p) for p in [anatomy_path,mechanics_path,BASE/'ihm/assembly/material_domains.py',BASE/'ihm/assembly/contact_dynamics.py',Path(__file__)]},
        'geometry_basis':'cell-center occupancy of exact-welded source surfaces; synthesized voxel/tetra volume, not measured voxel segmentation',
        'spacing_m':spacing_m,'boundary_discretization_diagonal_m':partition['boundary_discretization_diagonal_m'],
        'overlap_rule':partition['overlap_rule'],'source_priority':SOURCES,'overlapping_source_cells':partition['overlapping_source_cells'],
        'vertices':len(x),'tetrahedra':len(t),'boundary_triangles':len(body.triangles),'material_regions':materials,
        'mass_kg':float(body.mass_kg.sum()),'mass_claim_kg':ownership.mass_kg('pelvic-erectile-tissue-volume'),
        'mass_activation_contract':{'replaces_source_owners':SOURCES,'additive_mass_allowed':False,'canonical_runtime_handoff_applied':False},
        'minimum_reference_tetrahedron_volume_m3':float(body.region.volumes.min()),'max_explicit_dt_s':body.max_explicit_dt_s,
        'fixed_node_count':len(fixed),'support_basis':'engineering proximal posterior fixture within one cell of posterior extent; not measured anatomical insertion',
        'interfaces':['shared nodes bond internal material labels; sliding interfaces unresolved','tunica albuginea, fascial/skin shell, urethral lumen and perfused poroelastic state not resolved','no global organ or pelvis/garment collision activated'],
        'literature_review_date':'2026-09-05','literature':LITERATURE,'parameter_confidence_percent':None,
        'artifacts':{p.name:sha(p) for p in [out/'pelvic-domain.npz',out/'registry.json']}}
    (out/'manifest.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:receipt[k] for k in ['vertices','tetrahedra','mass_kg','max_explicit_dt_s','overlapping_source_cells']},indent=2))
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--spacing-m',type=float,default=.004);p.add_argument('--output-dir');a=p.parse_args();build(a.spacing_m,a.output_dir)
