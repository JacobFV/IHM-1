#!/usr/bin/env python3
"""Build mechanical bindings on the single canonical anatomy; preserve sources."""
from pathlib import Path
import gzip
from copy import deepcopy
import hashlib
import json
import sys
import numpy as np
from scipy.spatial import cKDTree
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))

SOURCES={
 'organs':{'url':'https://pubmed.ncbi.nlm.nih.gov/33176223/','finding':'Ex-vivo microindentation: kidney/liver E 0.5–3 kPa, heart spatially 1–30 kPa. Scale/species/measurement transfer is uncertain.'},
 'muscle_tension':{'url':'https://pubmed.ncbi.nlm.nih.gov/11181594/','finding':'In-vivo human soleus/tibialis anterior specific tension study. Generic 0.3 MPa is a synthesis prior, not every-muscle measurement.'},
 'tendon':{'url':'https://pubmed.ncbi.nlm.nih.gov/10562354/','finding':'In-vivo tibialis anterior tendon tangent Young modulus at maximal load 1.2 GPa. Transfer to other tendons is an assumption.'},
 'opensim':{'url':'https://github.com/opensim-org/opensim-models/tree/master/Models/Rajagopal','finding':'Rajagopal2016 source model muscle force, fiber, slack length and pennation parameters; native corrected wrapping export.'}}


def param(value,unit,prior,basis,group,sources=()):
    return {'value':float(value),'unit':unit,'prior_range':prior,'basis':basis,'dependency_group':group,'sources':list(sources),'confidence_percent':None}


def material(e):
    role=e['role'];name=e['name'].lower()
    if role=='rigid_bone':return {'density':param(1900,'kg/m3',[1500,2200],'assumed_prior','mass-allocation')}
    young=3000.;ran=[500.,30000.];source=['organs']
    if role=='muscle':young=20000.;ran=[5000.,100000.];source=[]
    elif role in ('tendon','ligament'):young=1.2e9 if role=='tendon' else 1e8;ran=[1e8,1.5e9] if role=='tendon' else [1e7,1e9];source=['tendon']
    elif role=='cartilage':young=1e6;ran=[1e5,1e7];source=[]
    elif role=='skin':young=50000.;ran=[5000.,1e6];source=[]
    elif role in ('vascular','connective_tissue'):young=100000.;ran=[10000.,1e6];source=[]
    elif role=='nerve':young=2000.;ran=[100.,100000.];source=[]
    elif 'liver' in name or 'kidney' in name:young=1500.;ran=[500.,3000.]
    elif 'heart' in name or 'ventricle' in name:young=10000.;ran=[1000.,30000.]
    elif 'lung' in name:young=2000.;ran=[500.,10000.];source=[]
    # nu is a regularization prior; homogeneous tissue bulk compliance is not
    # the gas/fluid cavity pressure-volume law, which is driven separately.
    nu=.45;mu=young/(2*(1+nu));lam=young*nu/((1+nu)*(1-2*nu))
    group='material-'+role
    out={'young_modulus':param(young,'Pa',ran,'source_informed_prior' if source else 'assumed_prior',group,source),
         'poisson_ratio':param(nu,'1',[.3,.49],'assumed_prior',group),
         'shear_modulus':param(mu,'Pa',None,'derived_from_E_nu',group,source),
         'lame_lambda':param(lam,'Pa',None,'derived_from_E_nu',group,source),
         'density':param(1900 if role=='rigid_bone' else 1000,'kg/m3',[1500,2200] if role=='rigid_bone' else [900,1100],'assumed_prior','mass-allocation')}
    if role=='muscle':out['active_specific_tension']=param(300000,'Pa',[100000,600000],'source_informed_prior','muscle-specific-tension',['muscle_tension'])
    return out


def readgeom(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:return json.load(f)


def points(e):
    g=readgeom(BASE/e['reference_geometry']['path'])
    return np.asarray(g.get('positions',g.get('vertices')),float).reshape(-1,3)


def main():
    folder=BASE/'data/derived/canonical';anatomy_path=folder/'anatomy.json'
    anatomy=json.loads(anatomy_path.read_text());entities=anatomy['entities'];byid={e['id']:e for e in entities}
    bones=[e for e in entities if e['role']=='rigid_bone'];bone_names={e['name'].lower():e for e in bones}
    specs=[]
    for e in entities:
        ext=np.maximum(np.array(e['bounds_m']['max'])-e['bounds_m']['min'],.0001)
        volume=e.get('volume_m3');basis=e.get('volume_method')
        if volume is None or volume<=1e-12:
            g=readgeom(BASE/e['reference_geometry']['path']);vertices=np.asarray(g.get('positions',[])).reshape(-1,3);faces=np.asarray(g.get('indices',[]),int).reshape(-1,3)
            signed=0.
            if len(faces):
                triangles=vertices[faces]-np.asarray(e['centroid_m'])
                signed=abs(float(np.einsum('ij,ij->i',triangles[:,0],np.cross(triangles[:,1],triangles[:,2])).sum()/6))
            volume=signed if 1e-12<signed<float(np.prod(ext)) else float(np.prod(ext)*np.pi/6)
            basis='unclosed-surface signed integral volume prior; open boundaries and orientation are unresolved' if volume==signed else 'ellipsoid bounds volume prior; no usable oriented surface integral'
        if e['role']=='skin':volume=float(e.get('surface_area_m2',0)*.000001 or volume);basis='1 micrometer numerical carrier for skin boundary; physical layer masses belong to skin-layer entities'
        elif e['role']=='vascular':volume=float(e.get('surface_area_m2',0)*.0003 or volume);basis='surface area times assumed 0.3 mm wall thickness; lumen excluded'
        mat=material(e);carrier=e['role'] in ('skin','lymphatic_network') or 'cavity of' in e['name'].lower()
        mass=0. if carrier else volume*mat['density']['value']
        # Isotropic moment with mean box principal inertia; rigid SO(3) dynamics
        # avoids an unmeasured body inertia tensor and its false precision.
        inertia=float(mass*np.sum(ext*ext)/18)
        axis=np.asarray(e.get('principal_axis',[0,1,0]));axis=axis/max(np.linalg.norm(axis),1e-12)
        specs.append({'id':e['id'],'name':e['name'],'role':e['role'],'constitutive':'rigid' if e['role']=='rigid_bone' else 'affine_boundary_carrier' if carrier else 'affine_neo_hookean',
                      'centroid_m':e['centroid_m'],'bounds_m':e['bounds_m'],'volume_m3':volume,'volume_basis':basis,'material':mat,
                      'mass_kg':mass,'mass_role':'numerical_boundary_carrier' if carrier else 'material_partition_proxy','inertia_diagonal_kg_m2':[inertia]*3,'fiber_axis':axis.tolist(),
                      'uncertainty':{'biological':'unquantified generic prior; no empirical calibration','discretization':'one affine solid per surface, rigid isotropic-inertia reduction'},
                      'reference_geometry':e['reference_geometry'],
                      **({'shell':deepcopy(e['shell'])} if 'shell' in e else {}),
                      **({'physical_surface_support':deepcopy(e['physical_surface_support'])} if 'physical_surface_support' in e else {})})
    total=sum(e['mass_kg'] for e in specs);carriers=sum(e['mass_role']=='numerical_boundary_carrier' for e in specs);factor=(77.1107029-carriers*1e-6)/total
    for e in specs:
        e['mass_kg']*=factor;e['inertia_diagonal_kg_m2']=[v*factor for v in e['inertia_diagonal_kg_m2']]
        if e['mass_role']=='numerical_boundary_carrier':
            e['mass_kg']=1e-6;e['inertia_diagonal_kg_m2']=[1e-8]*3
            e['material_volume_m3']=0.
        else:e['material_volume_m3']=e['volume_m3']
    specs_byid={e['id']:e for e in specs}
    # Full bone vertices, not bounding centers, determine attachment proximity.
    bv=[];bi=[];bone_arrays={}
    for e in bones:
        p=points(e);bone_arrays[e['id']]=p;bv.extend(p);bi.extend([e['id']]*len(p))
    bv=np.asarray(bv);tree=cKDTree(bv)
    def nearest(point,candidates=None):
        if candidates:
            best=min(((float(np.linalg.norm(p-point,axis=1).min()),e['id'],p[np.linalg.norm(p-point,axis=1).argmin()]) for e in candidates for p in [bone_arrays[e['id']]]),key=lambda z:z[0])
            return best[1],best[2],best[0]
        distance,index=tree.query(point);return bi[index],bv[index],float(distance)
    def link(a,b,point,kind,young=30000.,area=1e-4,length=.01):
        k=young*area/length
        ma=specs_byid[a]['mass_kg'];mb=specs_byid[b]['mass_kg'];reduced=ma*mb/(ma+mb)
        return {'a':a,'b':b,'point_m':np.asarray(point).tolist(),'kind':kind,'stiffness_n_m':k,'damping_ns_m':.5*np.sqrt(k*reduced),
                'parameter_basis':{'young_modulus_pa':young,'effective_area_m2':area,'effective_length_m':length,'stiffness':'EA/L','damping_ratio':.25,'basis':'assumed effective support prior','prior_range_multiplier':[.1,10.],'dependency_group':'connective-supports'},
                'attachment_basis':'geometric support prior; not a dissected anatomical insertion'}
    # A minimum-distance bone tree gives connected whole-body structural support.
    # Its edges are explicit generic connective constraints, not named joint claims.
    centers=np.array([e['centroid_m'] for e in bones]);n=len(bones);visited={0};edges=[]
    distance=np.linalg.norm(centers[:,None]-centers[None,:],axis=2);np.fill_diagonal(distance,np.inf)
    while len(visited)<n:
        best=min((distance[a,b],a,b) for a in visited for b in range(n) if b not in visited)
        _,a,b=best;visited.add(b);pa=centers[a];pb=centers[b]
        _,surface_a,_=nearest(pb,[bones[a]]);_,surface_b,_=nearest(pa,[bones[b]])
        edges.append(link(bones[a]['id'],bones[b]['id'],(surface_a+surface_b)/2,'inferred_skeletal_support',young=1e6,area=1e-4,length=max(np.linalg.norm(surface_a-surface_b),.005)))
    for e in entities:
        if e['role']=='rigid_bone':continue
        c=np.asarray(e['centroid_m']);bone_id,_,_=nearest(c)
        # Supports attach at tissue center to keep prescribed inflation free of
        # spurious net translation; they transmit externally applied forces.
        edges.append(link(e['id'],bone_id,c,'soft_tissue_support'))
    # Transfer source muscle path points using each native bone's geometry fit.
    app=json.loads((BASE/'data/derived/app/manifest.json').read_text());native=json.loads((BASE/'data/derived/opensim/native_corrected/baseline/mechanics.json').read_text())
    source=json.loads((BASE/'data/derived/anatomy/opensim__Rajagopal__Rajagopal2016.json').read_text());params={m['name']:m for m in source['muscles']}
    model=next(m for m in app['models'] if m['id']=='opensim-rajagopal');rot=np.array(model['display_transform']['rotation']);translation=np.array(model['display_transform']['translation'])
    bodygroups={}
    for body in native['body_transforms_ground']:
        side='right' if body.endswith('_r') else 'left';base=body.rsplit('_',1)[0]
        if body=='pelvis':candidates=[e for e in bones if 'hip bone' in e['name'] or e['name']=='sacrum']
        elif body=='torso':candidates=[e for e in bones if 'vertebra' in e['name'] or 'rib' in e['name'] or 'sternum' in e['name']]
        elif base=='toes':candidates=[e for e in bones if side in e['name'] and ('toe' in e['name'] or 'metatarsal' in e['name'])]
        elif base=='calcn':candidates=[bone_names[side+' calcaneus']]
        elif base=='hand':candidates=[e for e in bones if side in e['name'] and ('finger' in e['name'] or 'thumb' in e['name'] or 'metacarpal' in e['name'])]
        else:candidates=[e for e in bones if e['name'].lower()==side+' '+base]
        if not candidates:continue
        structures=[s for s in app['structures'] if s.get('body_frame')==body and s['id'].startswith('opensim-rajagopal-bone')]
        if not structures:continue
        p=np.concatenate([np.asarray(readgeom(BASE/'data/derived/app/geometry'/f"{s['id']}.json.gz")['positions']).reshape(-1,3) for s in structures]);q=np.concatenate([bone_arrays[e['id']] for e in candidates])
        low,high=p.min(0),p.max(0);ql,qh=q.min(0),q.max(0);scale=(qh-ql)/np.maximum(high-low,1e-5);offset=(ql+qh)/2-scale*(low+high)/2
        bodygroups[body]={'scale':scale,'offset':offset,'candidates':candidates,'source_bounds':{'min':low.tolist(),'max':high.tolist()},'target_bounds':{'min':ql.tolist(),'max':qh.tolist()}}
    names={'addbrev':'adductor brevis','addlong':'adductor longus','addmagDist':'adductor magnus','addmagIsch':'adductor magnus','addmagMid':'adductor magnus','addmagProx':'adductor magnus','bflh':'long head of {side} biceps femoris','bfsh':'short head of {side} biceps femoris','edl':'extensor digitorum longus','ehl':'extensor hallucis longus','fdl':'flexor digitorum longus','fhl':'flexor hallucis longus','gaslat':'lateral head of {side} gastrocnemius','gasmed':'medial head of {side} gastrocnemius','glmax1':'gluteus maximus','glmax2':'gluteus maximus','glmax3':'gluteus maximus','glmed1':'gluteus medius','glmed2':'gluteus medius','glmed3':'gluteus medius','glmin1':'gluteus minimus','glmin2':'gluteus minimus','glmin3':'gluteus minimus','grac':'gracilis','iliacus':'iliacus','perbrev':'fibularis brevis','perlong':'fibularis longus','piri':'piriformis','psoas':'psoas major','recfem':'rectus femoris','sart':'sartorius','semimem':'semimembranosus','semiten':'semitendinosus','soleus':'soleus','tfl':'tensor fasciae latae','tibant':'tibialis anterior','tibpost':'tibialis posterior','vasint':'vastus intermedius','vaslat':'vastus lateralis','vasmed':'vastus medialis'}
    muscle_entities=[e for e in entities if e['role']=='muscle'];muscle_names={e['name'].lower():e for e in muscle_entities};muscles=[];native_records=[];mapped=set()
    for m in native['muscles']:
        source_name=m['name'];base,sidecode=source_name.rsplit('_',1);side='right' if sidecode=='r' else 'left';label=names[base];label=label.format(side=side) if '{side}' in label else side+' '+label
        entity=muscle_names.get(label)
        if entity is None:raise ValueError(f'Missing canonical muscle {source_name} -> {label}')
        mapped.add(entity['id']);anchors=[];residuals=[]
        path=np.array(m['path_ground_m']);nodes=m['path_nodes'];nodepoints=np.array([p['point_ground_m'] for p in nodes])
        for j,p in enumerate(path):
            closest=int(np.argmin(np.linalg.norm(nodepoints-p,axis=1)));body=nodes[closest]['body'];group=bodygroups[body]
            transferred=(rot@p+translation)*group['scale']+group['offset'];id,surface,dist=nearest(transferred,group['candidates'])
            endpoint=j in (0,len(path)-1);position=surface if endpoint else transferred
            anchors.append({'entity_id':id,'point_m':position.tolist(),'source_body':body,'source_point_ground_m':p.tolist(),'surface_projection':endpoint,'projection_distance_m':dist if endpoint else None})
            if endpoint:residuals.append(dist)
        par=params[source_name]['parameters'];value=lambda key:float(par[key]['value'])
        length=sum(np.linalg.norm(np.array(b['point_m'])-a['point_m']) for a,b in zip(anchors,anchors[1:]));scale=length/m['length_m']
        fmax=value('max_isometric_force');fiber=value('optimal_fiber_length')*scale;slack=value('tendon_slack_length')*scale
        item={'id':'body-muscle-opensim-'+source_name,'source_name':source_name,'canonical_entity_id':entity['id'],'anchors':anchors,'rest_path_length_m':float(length),'max_isometric_force_n':fmax,'optimal_fiber_length_m':fiber,'tendon_slack_length_m':slack,'pennation_angle_rad':value('pennation_angle_at_optimal'),'passive_stiffness_n_m':fmax/max(fiber,1e-5)*.05,
              'parameter_basis':'source Fmax/pennation; source fiber/slack length scaled by registered path ratio; passive 5% Fmax per fiber length is assumed linearized prior',
              'prior_range_multiplier':[.5,2.],'dependency_group':'rajagopal-transfer','source_parameters':par,'source_sha256':source['sha256'],
              'registration':{'method':'per-bone axis-aligned bounds fit followed by endpoint projection to canonical bone surface','path_length_scale':scale,'endpoint_projection_distances_m':residuals,'uncertainty':'no anatomical insertion accuracy calibration; AABB fit is a synthesis assumption; wrapping frozen from native rest configuration'}}
        muscles.append(item);native_records.append({'source_name':source_name,'canonical_entity_id':entity['id'],'source_native_record':m,'registered_muscle_id':item['id']})
    for e in muscle_entities:
        if e['id'] in mapped:continue
        spec=specs_byid[e['id']];p=points(e);axis=np.asarray(spec['fiber_axis']);projection=(p-np.array(e['centroid_m']))@axis;ends=[p[projection.argmin()],p[projection.argmax()]]
        side='right' if 'right' in e['name'] else 'left' if 'left' in e['name'] else None
        candidates=[b for b in bones if side is None or ('left' not in b['name'] and 'right' not in b['name']) or side in b['name']]
        anchors=[]
        for endpoint in ends:
            id,point,dist=nearest(endpoint,candidates)
            anchors.append({'entity_id':id,'point_m':point.tolist(),'projection_distance_m':dist,'surface_projection':True})
        if anchors[0]['entity_id']==anchors[1]['entity_id']:
            # Distinct nearest insertion avoids a force pair entirely internal
            # to one bone; this is intentionally flagged as a weak prior.
            candidates=[b for b in candidates if b['id']!=anchors[0]['entity_id']]
            id,point,dist=nearest(ends[1],candidates);anchors[1]={'entity_id':id,'point_m':point.tolist(),'projection_distance_m':dist,'surface_projection':True}
        length=max(float(np.linalg.norm(np.array(anchors[1]['point_m'])-anchors[0]['point_m'])),1e-5);fiber=max(float(projection.max()-projection.min())*.5,1e-4)
        fmax=300000*spec['volume_m3']/fiber
        muscles.append({'id':'body-muscle-prior-'+e['id'],'canonical_entity_id':e['id'],'anchors':anchors,'rest_path_length_m':length,'max_isometric_force_n':fmax,'optimal_fiber_length_m':fiber,'tendon_slack_length_m':length*.1,'pennation_angle_rad':0.,'passive_stiffness_n_m':fmax/max(fiber,1e-5)*.05,
                        'parameter_basis':'PCSA=estimated volume / assumed fiber length (half surface principal extent); Fmax=PCSA*0.3 MPa; zero pennation, slack=10% path; passive=5% Fmax/fiber length',
                        'attachment_basis':'principal surface extrema projected onto nearest two same-side/axial bones; anatomical insertion unverified',
                        'prior_range_multiplier':[.1,10.],'dependency_group':'geometry-muscle-synthesis','sources':['muscle_tension']})
    # Native tendons are integrated into actuator force paths; each anatomical
    # tendon/ligament also has an explicit tension-only load path.
    connective=[]
    for e in entities:
        if e['role'] not in ('tendon','ligament'):continue
        spec=specs_byid[e['id']];p=points(e);axis=np.asarray(spec['fiber_axis']);proj=p@axis;ends=[p[proj.argmin()],p[proj.argmax()]]
        anchors=[]
        for endpoint in ends:
            id,point,dist=nearest(endpoint);anchors.append({'entity_id':id,'point_m':point.tolist(),'projection_distance_m':dist})
        if anchors[0]['entity_id']==anchors[1]['entity_id']:
            candidates=[b for b in bones if b['id']!=anchors[0]['entity_id']]
            id,point,dist=nearest(ends[1],candidates);anchors[1]={'entity_id':id,'point_m':point.tolist(),'projection_distance_m':dist}
        length=max(float(np.linalg.norm(np.array(anchors[1]['point_m'])-anchors[0]['point_m'])),1e-4)
        connective.append({'id':'body-connective-'+e['id'],'canonical_entity_id':e['id'],'anchors':anchors,'rest_path_length_m':length,'max_isometric_force_n':0.,'optimal_fiber_length_m':length,'tendon_slack_length_m':length,'pennation_angle_rad':0.,'passive_stiffness_n_m':spec['material']['young_modulus']['value']*spec['volume_m3']/length**2,'parameter_basis':'tension-only EA/L with area=estimated volume/length; stress-free reference; anatomical insertions inferred','prior_range_multiplier':[.1,10.],'dependency_group':'connective-geometry-synthesis'})
    muscles.extend(connective)
    # Ensure muscle-shape centers receive the motion of origin/insertion rather
    # than only an unrelated nearest static support.
    for m in muscles:
        id=m['canonical_entity_id'];c=np.array(byid[id]['centroid_m'])
        for anchor in [m['anchors'][0],m['anchors'][-1]]:
            edges.append(link(id,anchor['entity_id'],c,'muscle_or_connective_attachment_support',young=20000.,area=1e-5,length=.02))
    payload={'schema_version':1,'model_id':'ihm-body','frame':'bodyparts3d-display-m','units':{'length':'m','mass':'kg','time':'s','force':'N','stress':'Pa','energy':'J'},'entities':specs,'links':edges,'muscles':muscles,'native_muscles':native_records,
             'source_files':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [anatomy_path,BASE/'data/derived/opensim/native_corrected/baseline/mechanics.json',BASE/'data/derived/anatomy/opensim__Rajagopal__Rajagopal2016.json']},'sources':SOURCES,
             'registration':{key:{'scale':g['scale'].tolist(),'translation_m':g['offset'].tolist(),'canonical_bones':[e['id'] for e in g['candidates']],'source_bounds':g['source_bounds'],'target_bounds':g['target_bounds']} for key,g in bodygroups.items()},
             'mass_allocation':{'target_mass_kg':77.1107029,'numerical_carrier_mass_kg':carriers*1e-6,'numerical_carrier_count':carriers,'exclusions':'skin parent, lymphatic network structural graph and cardiac cavities carry only 1 milligram numerical inertia each; material volume zero; this inertia is debited from the material mass allocation','unscaled_proxy_mass_kg':total,'uniform_scale':factor,'basis':'generic physiology mass constraint, allocated by estimated tissue volume and density; overlapping atlas representations and bounding-volume approximations are not a measured compartment partition'},
             'rest_state':'stress-free anatomy, zero activation, no gravity; no native pretension transplanted; unsupported gait/postural predictions',
             'reduction_priors':{'active_force_length_width':param(.45,'normalized fiber length',[.25,.75],'assumed_prior','reduced-muscle-law'),'max_log_shortening':param(.35,'1',[.1,.5],'assumed_prior','reduced-muscle-law'),'passive_fmax_fraction_per_fiber_length':param(.05,'1',[.01,.2],'assumed_prior','reduced-muscle-law')},
             'solver':{'rigid':'rigid translation; reference orientations constrained with audited holding moments; linearly implicit attachment stiffness/damping, 2 ms substeps','soft':'affine compressible neo-Hookean; quasi-static hydrostatic boundary solve and volume-preserving muscle shape','muscle':'registered polyline force gradient; active Gaussian force-length reduction, passive tension-only linear spring; NOT reimplementation of Millard equilibrium','coupling':'equal/opposite forces on every path segment and support, moments about body centers','gravity':'disabled; whole body is in unloaded reference configuration'},
             'limitations':['No collision/contact, sliding fascia, joint rotation or physiological postural equilibrium; bone orientations are constrained.','Geometric support tree and synthesized insertions require anatomical refinement; connected does not mean every connection anatomically measured.','One affine solid per surface cannot resolve local strain, folds, organ contact, wall cavities or tissue heterogeneity.','Native wrapped paths frozen at source rest pose and transferred; generic force-length response is an explicitly different reduced model.','Inherited nonlinear/coupled uncertainty is unquantified; prior ranges are sensitivity envelopes, not probability intervals.'],
             'counts':{'entities':len(specs),'rigid_bones':len(bones),'soft_solids':sum(e['constitutive']=='affine_neo_hookean' for e in specs),'boundary_carriers':sum(e['constitutive']=='affine_boundary_carrier' for e in specs),'muscle_actuators':len(muscles)-len(connective),'native_muscle_paths':len(native_records),'synthesized_muscles':len(muscles)-len(connective)-len(native_records),'tendon_ligament_paths':len(connective),'support_links':len(edges)}}
    folder.mkdir(parents=True,exist_ok=True);(folder/'mechanics.json').write_text(json.dumps(payload,separators=(',',':'),allow_nan=False)+'\n')
    print(json.dumps(payload['counts']))

if __name__=='__main__':main()
