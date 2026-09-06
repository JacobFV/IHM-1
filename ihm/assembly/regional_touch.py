"""Body-attached reference contact and causal IBM tactile materialization.

The regional solid is a separate explicit materialization. It does not add a
second mass or force to the whole-body affine mechanics. No guessed response
unit conversion into neural firing or muscle force is made.
"""
from pathlib import Path
import gzip
import hashlib
import json
import numpy as np
from .mechanics_backend import DeformableRegion,tetra_box


def forearm_anchor(root):
    root=Path(root)
    anatomy_path=root/'data/derived/canonical/anatomy.json'
    anatomy=json.loads(anatomy_path.read_bytes())
    micro_path=root/'data/derived/canonical/microvascular.json'
    micro=json.loads(micro_path.read_bytes())
    attachment=micro['units'][0]['material_attachment']
    skin=next(e for e in anatomy['entities'] if e['id']==attachment['entity_id'])
    source=skin['reference_geometry'];raw=(root/source['path']).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=source['sha256']:raise ValueError('Canonical skin geometry changed')
    g=json.loads(gzip.decompress(raw));vertices=np.array(g['positions']).reshape(-1,3);faces=np.array(g['indices']).reshape(-1,3)
    tri=vertices[faces[attachment['face_index']]]
    origin=np.array(attachment['barycentric'])@tri
    tangent=tri[1]-tri[0];tangent/=np.linalg.norm(tangent)
    normal=np.cross(tangent,tri[2]-tri[0]);normal/=np.linalg.norm(normal)
    frame=np.array([tangent,np.cross(normal,tangent),normal])
    return {'entity_id':skin['id'],'face_index':attachment['face_index'],'barycentric':attachment['barycentric'],
        'origin_m':origin.tolist(),'local_axes':frame.tolist(),'frame':anatomy['frame'],
        'geometry_sha256':source['sha256'],'geometry_path':source['path'],
        'region':'left forearm','source_hashes':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [anatomy_path,micro_path]}}


def material_config():
    # Apparent shear stiffness under light indentation is a transfer prior for
    # this homogeneous solid, not a fitted hyperelastic dermal parameter.
    return {'mu_pa':2800.,'lambda_pa':25200.,'density_kg_m3':1000.,
        'parameter_evidence':{
            'mu_pa':{'value':2800.,'reported_sd_pa':800.,'kind':'transferred_human_measurement_prior',
                'source':'https://pubmed.ncbi.nlm.nih.gov/19152581/','doi':'10.1111/j.1600-0846.2008.00329.x',
                'cohort':'20 subjects aged55–70; in vivo total skin, light indentation/static friction',
                'transfer':'Apparent G* is used as a homogeneous neo-Hookean small-strain prior; not this generic subject or a layer fit.'},
            'lambda_pa':{'value':25200.,'kind':'constitutive_prior','assumption':'Poisson ratio0.45; lambda=2*mu*nu/(1-2*nu). No bulk-modulus measurement.'},
            'density_kg_m3':{'value':1000.,'kind':'waterlike_density_prior','used_for':'Mass only; experiment is quasistatic with zero gravity.'}}}


def run_touch(root,*,divisions=(32,32,16),indentation_m=.0001,dt_s=.002,duration_s=.8,onset_s=.1,release_s=.4,delay_s=.025):
    from ihm.brain.ibm_backend import IBMBackend
    from ihm.brain.causal import CausalIBM
    values=[indentation_m,dt_s,duration_s,onset_s,release_s,delay_s]
    if not np.isfinite(values).all() or not 0<indentation_m<=.0003 or not 0<dt_s<=.01 or not 0<=onset_s<release_s<duration_s<=10 or delay_s<0:
        raise ValueError('Invalid bounded regional pulse')
    if any(abs(x/dt_s-round(x/dt_s))>1e-8 for x in [onset_s,release_s,duration_s]):raise ValueError('Pulse boundaries must align to sample clock')
    root=Path(root).resolve()
    backend=IBMBackend(root=root)
    anchor=forearm_anchor(root);config=material_config()
    x,t=tetra_box((.01,.01,.003),divisions)
    solid=DeformableRegion(x,t,**{k:config[k] for k in ['mu_pa','lambda_pa','density_kg_m3']})
    bottom=np.flatnonzero(x[:,2]==0)
    reference=solid.solve(fixed_nodes=bottom)
    loaded=solid.solve(fixed_nodes=bottom,indenter={'center_m':[.005,.005],'radius_m':.003,'height_m':.003-indentation_m})
    released=solid.solve(fixed_nodes=bottom)
    site=int(np.argmin(np.linalg.norm(x-[.005,.005,.003],axis=1)))
    world=lambda p:(np.asarray(p)-[.005,.005,.003])@np.array(anchor['local_axes'])+anchor['origin_m']
    sites=world(x[[site]])
    receptors={name:CausalIBM(backend,sites_m=sites,kind=name,delay_s=delay_s) for name in ['rapid','slow']}
    frames=[]
    def sample(time,state,responses):
        depth=float(x[site,2]-state['positions_m'][site][2])
        return {'time_s':time,'mechanical_state':'loaded' if state is loaded else 'released' if state is released else 'reference',
            'indentation_m':depth,'reaction_n':state['indenter_reaction_n'],'elastic_energy_j':state['elastic_energy_j'],
            'rapid_response':float(responses['rapid']),'slow_response':float(responses['slow'])}
    frames.append(sample(0,reference,{'rapid':0.,'slow':0.}))
    for i in range(round(duration_s/dt_s)):
        time=i*dt_s
        state=loaded if onset_s<=time<release_s else released if time>=release_s else reference
        indentation=max(0.,x[site,2]-state['positions_m'][site][2])*1e6
        response={kind:model.advance([indentation],dt_s)[0] for kind,model in receptors.items()}
        frames.append(sample((i+1)*dt_s,state,response))
    states={name:{**state,'positions_m':world(state['positions_m']).tolist()} for name,state in [('reference',reference),('loaded',loaded),('released',released)]}
    paths=['ihm/assembly/regional_touch.py','ihm/assembly/mechanics_backend.py','ihm/brain/causal.py','ihm/brain/ibm_backend.py','ihm/brain/source_loader.py','ihm/brain/active_source.py','ihm/brain/candidate.py']
    _manifest=backend.artifact_dir/'manifest.json'
    return {'schema_version':1,'id':'forearm-contact-ibm','anchor':anchor,'material':config,
        'geometry':{'reference_positions_m':world(x).tolist(),'tetrahedra':t.tolist(),'local_dimensions_m':[.01,.01,.003],
            'domain_kind':'homogeneous reference volume tangent to pinned skin face; inferred thickness and planar footprint'},
        'clock':{'dt_s':dt_s,'duration_s':duration_s,'onset_s':onset_s,'release_s':release_s,'semantics':'zero-order held mechanical pulse; exact causal endpoint receptor response'},
        'states':states,'frames':frames,'receptors':{k:v.audit for k,v in receptors.items()},
        'checkpoints':{k:v.checkpoint() for k,v in receptors.items()},
        'runtime_sources':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths},
        'source_hashes':{str(_manifest.relative_to(root)):hashlib.sha256(_manifest.read_bytes()).hexdigest()},
        'executed_source':backend.identity_summary(),
        'coupling':{'solid_owner':'regional neo-Hookean reference solve','input':'computed normal surface displacement in m, converted exactly to um',
            'output':'signed donor response','native_blood_storage_connected':False,'whole_body_force_feedback':False,'calibrated_afferent_firing':False},
        'limitations':['Local planar homogeneous patch is an explicitly synthesized alternative materialization; no claim of exact layer partition or curved-surface conformity.',
            'Quasistatic endpoints omit dynamic wave propagation, friction, contact work during ramps and viscoelastic relaxation.',
            'Receptor class dynamics are source matched; gains, receptor census and nerve delay are priors.',
            'Regional response is not converted to brain Hz, neural mV or muscle force without an identified observation law.']}
