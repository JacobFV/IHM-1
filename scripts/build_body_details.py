#!/usr/bin/env python3
"""Build reproducible skin detail; --append explicitly adds display to app manifest.

Canonical anatomy.json and manifest_fragment.json are always read-only inputs.
"""
from pathlib import Path
import argparse,gzip,hashlib,io,json,platform,sys,zipfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.details import physical_mesh,sample_hair,microvascular_unit,solve_network,shaft_mesh
from ihm.assembly.hair import regional_density,display_indices,PRIORS,outward_surface_mask
OUT=ROOT/'data/derived/canonical'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write_json(path,data):path.write_text(json.dumps(data,separators=(',',':'),allow_nan=False)+'\n')
def write_npz(path,arrays):
    # Fixed ZIP timestamps and sorted members yield byte reproducibility.
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for key,value in sorted(arrays.items()):
            buf=io.BytesIO();np.lib.format.write_array(buf,np.asarray(value),allow_pickle=False)
            info=zipfile.ZipInfo(key+'.npy',date_time=(1980,1,1,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;z.writestr(info,buf.getvalue())
def geometry(sid,g):
    path=OUT/'geometry'/f'{sid}.json.gz';path.write_bytes(gzip.compress(json.dumps(g,separators=(',',':'),allow_nan=False).encode(),mtime=0))
    return dict(geometry_url='/api/geometry/'+sid,geometry_path=str(path.relative_to(ROOT)),geometry_sha256=digest(path))
def serial(state):return {k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in state.items()}

def build(display_count=30000,append=False):
    import scipy,trimesh
    anatomy=json.loads((OUT/'anatomy.json').read_text());skin=next(e for e in anatomy['entities'] if e['id']=='body-bp3d-FJ2810')
    path=ROOT/skin['reference_geometry']['path'];g=json.loads(gzip.decompress(path.read_bytes()))
    v=np.asarray(g['positions'],float).reshape(-1,3);f=np.asarray(g['indices'],int).reshape(-1,3)
    mesh=physical_mesh(v,f)
    centers=v[f].mean(axis=1);normal=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]);normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-30)
    density,region=regional_density(v[f].mean(axis=1))
    # Require every corner to agree, preventing straddling into exclusion masks.
    _,corner=regional_density(v[f].reshape(-1,3));valid=np.all(corner.reshape(-1,3)==region[:,None],axis=1)
    valid &= outward_surface_mask(centers,normal,region)
    region[~valid]='boundary_or_inward_excluded';density[~valid]=0
    h=sample_hair(v,f,density,seed=20260905)
    h['region']=region[h['face_index']];h['face_density_per_cm2']=density
    h['source_entity_id']=np.array(skin['id']);h['source_geometry_sha256']=np.array(digest(path))
    full=OUT/'hair_samples.npz';write_npz(full,h)
    sample_ids=display_indices(h,display_count);display=OUT/'hair_display_attachment.npz'
    write_npz(display,{k:h[k][sample_ids] for k in ['ids','face_index','barycentric','radius_m','length_m']})
    sid='body-detail-hair';hair_g=shaft_mesh(h['roots_m'][sample_ids],h['normals'][sample_ids],h['length_m'][sample_ids],h['radius_m'][sample_ids])
    hair_g['attachment']={'kind':'MaterialPoint','skin_entity_id':skin['id'],'samples_path':str(display.relative_to(ROOT)),'source_geometry_sha256':digest(path),'vertices_per_shaft':8,'deformation_rule':'root=sum(barycentric_i * deformed_skin_vertex_i); rotate shaft local frame with deformed face tangent/normal; physical radius and length unchanged','rigid_follow_only':False}
    structures=[dict(id=sid,name='Regional vellus hair — sampled prior',model_id='ihm-body',system='hair',kind='mesh',color='#65513e',default_visible=True,evidence_kind='regional_human_density_prior',canonical_entity_id=skin['id'],**geometry(sid,hair_g))]
    units=[];allp=[];alle=[];allr=[]
    centers=v[f].mean(axis=1);normal=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]);normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-30)
    # Four bounded demonstration units per forearm. Their count is not capillary density.
    for side,sign in [('left',1),('right',-1)]:
        for j in range(4):
            target=np.array([sign*.28,.13+j*.035,.04]);eligible=np.flatnonzero((region=='forearm')&(centers[:,0]*sign>0))
            face=int(eligible[np.argmin(np.linalg.norm(centers[eligible]-target,axis=1))]);c=centers[face]-.001*normal[face]
            p,e,r,k=microvascular_unit(np.zeros(3),.003,32,seed=101+j)
            # Align local z to skin normal, putting the schematic unit below surface.
            n=normal[face];a=np.cross(n,[0.,1.,0.]);a/=np.linalg.norm(a);b=np.cross(n,a);p=p@np.array([a,b,n])+c
            s=solve_network(p,e,r,{0:4500.,1:1500.})
            uid=f'body-detail-microvascular-{side}-{j}';offset=sum(len(x) for x in allp);allp.append(p);alle.append(e+offset);allr.append(r)
            units.append(dict(id=uid,territory=f'{side} forearm skin',positions_m=p.tolist(),edges=e.tolist(),radius_m=r.tolist(),edge_kind=k.tolist(),edge_kind_legend={'0':'arteriolar_supply','1':'venular_return','2':'capillary'},state=serial(s),
                material_attachment={'entity_id':skin['id'],'face_index':face,'barycentric':[1/3]*3,'normal_depth_m':.001},
                boundary_association={'kind':'named_territory_boundary','supply':f'{side} forearm cutaneous arterial territory','return':f'{side} forearm cutaneous venous territory','pressure_nodes':{'supply':0,'return':1},'evidence_kind':'inferred_territory_prior','native_blood_storage_connected':False,'native_parent_vessel_id':None,'uncertainty':'No measured terminal vessel correspondence; prescribed isolated pressure boundaries'},
                assumptions=['Synthetic paired tree; capillary count/length/radius not measured human skin morphology','Pressure boundary 4500/1500 Pa is a coefficient prior','Steady passive resistors; no storage or exchange coupling']))
    p=np.concatenate(allp);e=np.concatenate(alle);r=np.concatenate(allr);delta=p[e[:,1]]-p[e[:,0]]
    sid='body-detail-microvascular';structures.append(dict(id=sid,name='Forearm microvascular units — inferred paired networks',model_id='ihm-body',system='microvascular',kind='mesh',color='#b75d79',default_visible=True,evidence_kind='inferred_territory_network',**geometry(sid,shaft_mesh(p[e[:,0]],delta,np.linalg.norm(delta,axis=1),r))))
    vascular={'schema':'ihm.microvascular.v1','units':units,'native_blood_storage_connected':False,'rheology':{'kind':'heterogeneous_apparent_viscosity_coefficient_prior','diameter_bins_um':[12,30],'viscosities_pa_s':[.004,.0024,.0018],'audited_pries_law':False,'hematocrit_model':None,'limitations':'Coefficient priors, not calibrated human capillary rheology; no RBC phase separation, shear dependence or hematocrit transport'},'density_per_mm2':None}
    write_json(OUT/'microvascular.json',vascular)
    fragment={'structures':structures};write_json(OUT/'details_manifest_fragment.json',fragment)
    inputs=[OUT/'anatomy.json',path,ROOT/'ihm/assembly/details.py',ROOT/'ihm/assembly/hair.py',Path(__file__),ROOT/'docs/research/REGIONAL_TRANSPORT_AND_SKIN.md']
    report={'schema':'ihm.details.v1','frame':'bodyparts3d-display-m','seed':20260905,'source_hashes':{str(x.relative_to(ROOT)):digest(x) for x in inputs},'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'trimesh':trimesh.__version__,'python_executable_sha256':digest(Path(sys.executable).resolve())},
        'skin_physical_mesh':{'source_vertices':len(v),'welded_vertices':len(mesh.vertices),'faces':len(f),'watertight':bool(mesh.is_watertight),'volume_m3':float(abs(mesh.volume)) if mesh.is_watertight else None,'volume_interpretation':'enclosed skin tissue shell only, NOT whole-body occupancy','welding':'exact coordinate equality within single canonical skin entity; no proximity weld'},
        'hair':{'full_count':len(h['ids']),'display_count':len(sample_ids),'samples_path':str(full.relative_to(ROOT)),'samples_sha256':digest(full),'display_attachment_path':str(display.relative_to(ROOT)),'display_attachment_sha256':digest(display),'identity':'(source entity ID, pinned geometry hash, seed, integer sample ID); independent of display count','density_prior':{k:{'mean_per_cm2':x[0],'sd_per_cm2':x[1],'n':15} for k,x in PRIORS.items()},'region_counts':{k:int(np.sum(h['region']==k)) for k in PRIORS},'display_reduction_applied':len(sample_ids)<len(h['ids']),'physical_radius_multiplier':1,'morphology_prior':'Vellus shaft radius uniform 8–16 um; length uniform .5–2 mm; neither measured cohort distribution nor complete terminal/scalp hair model','attachment':'Source triangle index and barycentric MaterialPoint; apply deformed source skin vertices with hair.attached_roots, transport shaft frame with face tangents. Renderer integration required for nonrigid animation.'},
        'microvascular':{'unit_count':len(units),'capillary_count':sum(sum(np.asarray(u['edge_kind'])==2) for u in units),'edge_count':len(e),'node_count':len(p),'maximum_internal_residual_m3_s':max(u['state']['maximum_internal_residual_m3_s'] for u in units),'native_blood_storage_connected':False,'territory':'bilateral forearm only; eight illustrative units, no density claim'},
        'limitations':['Atlas-coordinate masks are inferred regions, not measured skin classification','Outward anatomical-axis face filter excludes inward skin shell; inferred and conservative, not measured segmentation','Hands/feet completely excluded to protect palms/soles; face except forehead, scalp, groin, axilla and mask boundaries unresolved/excluded','No invented canonical native artery connection; boundary-only transport','Material samples are attached; app must consume deformation metadata to move shafts nonrigidly']}
    report['microvascular']['capillary_count']=int(report['microvascular']['capillary_count'])
    write_json(OUT/'details.json',report)
    if append:
        manifest_path=ROOT/'data/derived/app/manifest.json';manifest=json.loads(manifest_path.read_text());ids={s['id'] for s in structures}
        manifest['structures']=[s for s in manifest['structures'] if s['id'] not in ids]+structures;write_json(manifest_path,manifest)
    print(json.dumps({'hair_samples':len(h['ids']),'display_hairs':len(sample_ids),'vascular_edges':len(e),'appended':append}))
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--display-count',type=int,default=30000);ap.add_argument('--append',action='store_true');args=ap.parse_args();build(args.display_count,args.append)
