#!/usr/bin/env python3
"""Assess how 2,4xx per-entity conforming bodies could be coupled into one system.

Read-only. Writes only inside the chosen output directory. Nothing is promoted,
no canonical file is touched and no physics is executed.

Five bounded questions:

1. What contact machinery exists, and what has it actually been exercised on?
   (recorded from existing receipts; this script re-reads them, it does not rerun.)
2. What does broad phase cost, per materialization class, from real bounding
   boxes and a measured surface-adjacency graph?
3. Which adjacent entity pairs are anatomically licensed to SLIDE? The extended
   atlas' bursae/sheaths/fascia are registered into the canonical frame with the
   recorded z_anatomy transform and asked which measured contact patch they sit on.
4. Can the 32.117 L unmodelled interstitium act as the tied coupling medium?
5. What element counts follow from a tet-vs-embedded-fibre choice per class.

Every geometric statement here is a proxy on an 8/4/2 mm occupancy grid or on
axis-aligned bounds. It is a screening classification for building a contact
graph, not a certified anatomical adjacency and not a measured tissue boundary.
"""
import argparse,gzip,hashlib,json,re,sys,time
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.anatomy import LandmarkRegistration

CANONICAL='data/derived/canonical/anatomy.json'
BODY='data/derived/canonical/body.json'
INVENTORY='data/derived/joint-substrate-assessment-v1/inventory.json'
SUBSTRATE_MANIFEST='data/derived/joint-substrate-assessment-v1/manifest.json'
CARTILAGE_GAP='data/derived/joint-substrate-assessment-v1/cartilage_gap.json'
SOURCE_INDEX='data/derived/anatomy/extended/source_index.json'
UNMODELLED='data/derived/unmodelled-volume/summary.json'
VOID_COMPONENTS='data/derived/unmodelled-volume/void_components.jsonl'
WHOLE_BODY_DOMAIN='data/derived/material-domains/whole-body-0.01m/manifest.json'
PELVIS_DOMAIN='data/derived/material-domains/pelvis-0.002m/manifest.json'
FOREARM='data/derived/canonical/forearm-touch.json'
TET_SOURCES=(('entity','data/derived/entity-tet-ready-v1'),('muscle','data/derived/muscle-tet-ready-v1'))
GARMENT_RECEIPT='data/derived/garment-tissue-refined-oufb46p6/coupled/report.json'
CONTACT_SOURCE='ihm/assembly/contact_dynamics.py'
FOUNDATION_SOURCE='scripts/native_surface_foundation.h'
SUPINE_QUADRATURE='data/derived/supine-surface-contact-exmzq9pq/manifest.json'

ROTATION=np.array([[1.,0,0],[0,0,1],[0,-1,0]])
GRID_M=(0.008,0.004,0.002)
ADJACENCY_GRID_M=0.004
MARKER_DILATIONS=(1,2)
COVERAGE_THRESHOLDS=(0.25,0.5,0.75)
# First match wins. Splits the atlas' single fascia_aponeurosis bucket, because a
# deep fascial plane is a glide surface and an aponeurosis is a tendon insertion.
SLIDING_CLASSES={'synovial_bursa':'synovial bursa: interposed serous sac, the defining free-glide interface',
 'tendon_sheath':'synovial tendon sheath: tendon glides inside it',
 'articular_capsule':'articular capsule: encloses a synovial cavity; the enclosed bone pair articulates',
 'retinaculum':'retinaculum: tendons glide beneath it',
 'meniscus':'intra-articular fibrocartilage: slides on at least the femoral face',
 'articular_disc':'intra-articular disc: two synovial compartments',
 'labrum':'glenoid/acetabular labrum: rim of a synovial cavity',
 'fat_pad':'intra-articular fat pad: deformable sliding spacer',
 'fascia':'deep fascial plane: muscle groups glide across it'}
TIED_CLASSES={'ligament':'ligament: bone-to-bone dense collagen, load bearing in tension',
 'cruciate_ligament':'cruciate ligament: intracapsular bone-to-bone tie',
 'interosseous_membrane':'interosseous/obturator membrane: continuous fibrous sheet between two bones',
 'intervertebral_disc':'intervertebral disc: annulus is continuous with both vertebral endplates',
 'nucleus_pulposus':'nucleus pulposus: enclosed by the annulus, no free surface',
 'symphysis':'symphysis: fibrocartilaginous union',
 'aponeurosis':'aponeurosis: flat tendon, an insertion not a glide plane'}
AMBIGUOUS_CLASSES={'cartilage':'cartilage: tied to its own bone, sliding on the opposing cartilage; the class alone does not resolve which face'}
# Two entities that are an authored partition of one continuous structure cannot slide
# on each other whatever marker happens to overlie them.
ENUMERATION=re.compile(r'\b(?:[ivxlc]+|\d+|first|second|third|fourth|fifth|sixth|seventh|eighth|ninth|tenth|'
    r'eleventh|twelfth|left|right|superior|inferior|middle|anterior|posterior|medial|lateral|upper|lower|'
    r'proximal|distal|deep|superficial|internal|external|greater|lesser|major|minor)\b')


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as s:
        for b in iter(lambda:s.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def dump(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False,sort_keys=False)+'\n')


def dump_lines(path,rows):
    with path.open('w') as s:
        for r in rows:s.write(json.dumps(r,allow_nan=False)+'\n')


def encode(index):
    """Pack a signed voxel index triple into one int64. Offsets cover the atlas."""
    if index.min()<-512 or index.max()>=524288-512:raise ValueError('Voxel index outside packing range')
    return ((index[:,0]+512)<<40)|((index[:,1]+512)<<20)|(index[:,2]+512)


def surface_samples(triangles,h):
    """Barycentric lattice dense enough that no (h) cell crossed by a face is missed."""
    area=np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)/2
    per=np.clip(np.ceil(np.sqrt(area/((h/2)**2))).astype(np.int64),1,64)
    out=[]
    for count in np.unique(per):
        picked=triangles[per==count]
        u=(np.arange(count)+.5)/count
        a,b=np.meshgrid(u,u,indexing='ij')
        w=np.stack([a.ravel(),b.ravel()],axis=1)
        w=np.where(w.sum(axis=1,keepdims=True)>1,1-w,w)
        bary=np.concatenate([1-w.sum(axis=1,keepdims=True),w],axis=1)
        out.append((picked[:,None,:,:]*bary[None,:,:,None]).sum(axis=2).reshape(-1,3))
    return np.concatenate(out,axis=0)


def voxel_keys(vertices,triangles,h):
    points=np.concatenate([vertices,surface_samples(triangles,h)],axis=0)
    return np.unique(encode(np.floor(points/h).astype(np.int64)))


def read_geometry(path):
    payload=json.loads(gzip.open(path).read())
    return (np.asarray(payload['positions'],float).reshape(-1,3),
            np.asarray(payload['indices'],np.int64).reshape(-1,3))


def load_bodies():
    """Every tet-proof-carrying entity surface, with exact area/volume/bounds."""
    records=[]
    for family,directory in TET_SOURCES:
        for line in (ROOT/directory/'entities.jsonl').open():
            row=json.loads(line);row['_family']=family;row['_dir']=directory;records.append(row)
    bodies=[]
    for record in records:
        vertices,faces=read_geometry(ROOT/record['_dir']/record['output_path'])
        tri=vertices[faces];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        area=float(np.linalg.norm(normal,axis=1).sum()/2)
        volume=float(np.einsum('ij,ij->i',tri[:,0],normal).sum()/6)
        bodies.append({'entity_id':record['entity_id'],'family':record['_family'],
            'tet_ready':bool(record.get('tet_ready')),'vertex_links_manifold':bool(record.get('vertex_links_manifold')),
            'vertices':len(vertices),'faces':len(faces),'surface_area_m2':area,'signed_volume_m3':volume,
            'radius_proxy_m':(2*volume/area if area>0 and volume>0 else None),
            'lo':vertices.min(axis=0),'hi':vertices.max(axis=0),'_v':vertices,'_f':faces})
    return bodies


def aabb_pairs(lo,hi,margin):
    """Exact axis-aligned overlap count on inflated bounds."""
    low=lo-margin;high=hi+margin;n=len(low)
    overlap=np.ones((n,n),bool)
    for k in range(3):
        overlap&=(low[:,k][:,None]<=high[:,k][None,:])&(low[:,k][None,:]<=high[:,k][:,None])
    np.fill_diagonal(overlap,False)
    return overlap


def adjacency(bodies,h):
    """Measured surface co-occupancy: two bodies whose surfaces enter one h cell."""
    keys=[];owners=[]
    for i,b in enumerate(bodies):
        k=voxel_keys(b['_v'],b['_v'][b['_f']],h)
        keys.append(k);owners.append(np.full(len(k),i,np.int32))
    key=np.concatenate(keys);owner=np.concatenate(owners)
    order=np.argsort(key,kind='stable');key=key[order];owner=owner[order]
    cut=np.flatnonzero(np.diff(key))+1
    starts=np.concatenate(([0],cut));ends=np.concatenate((cut,[len(key)]))
    counts=ends-starts
    shared=defaultdict(list)
    for s,e in zip(starts[counts>1],ends[counts>1]):
        members=np.unique(owner[s:e])
        for x in range(len(members)):
            for y in range(x+1,len(members)):
                shared[(int(members[x]),int(members[y]))].append(int(key[s]))
    return {'grid_h_m':h,'occupied_surface_voxels':int(len(starts)),'surface_voxel_entries':int(len(key)),
            'mean_bodies_per_occupied_surface_voxel':float(counts.mean()),
            'max_bodies_in_one_surface_voxel':int(counts.max())},shared


def register_substrate(inventory,index,fit,h,dilations):
    """Z-Anatomy joint substrate carried into the canonical frame, then voxelised."""
    by_id={m['id']:m for m in index['meshes']}
    offsets={d:np.array([[a,b,c] for a in range(-d,d+1) for b in range(-d,d+1) for c in range(-d,d+1)]) for d in dilations}
    rows=[]
    for structure in inventory['structures']:
        mesh=by_id[structure['id']]
        with np.load(ROOT/mesh['source_geometry_path']) as archive:
            vertices=archive['vertices'].astype(float);faces=archive['faces'].astype(np.int64)
        moved=fit.transform(vertices@ROTATION.T)
        cells=np.unique(np.floor(np.concatenate([moved,surface_samples(moved[faces],h)],axis=0)/h).astype(np.int64),axis=0)
        rows.append({'structure':structure,'interface_class':interface_class(structure),
            'canonical_bounds_m':[moved.min(axis=0).tolist(),moved.max(axis=0).tolist()],
            'footprint_voxels':int(len(cells)),
            'dilated':{d:np.unique(encode(np.unique((cells[:,None,:]+offsets[d][None]).reshape(-1,3),axis=0))) for d in dilations}})
    return rows


def interface_class(structure):
    label=structure['joint_class']
    if label=='fascia_aponeurosis':
        return 'aponeurosis' if 'aponeuros' in structure['name'].lower() else 'fascia'
    return label


def name_stem(name):
    """Strip laterality and enumeration so partition members of one organ collide."""
    return ' '.join(w for w in ENUMERATION.sub(' ',name.lower()).split() if w)


def continuity_prior(a_name,b_name,a_system,b_system,a_role,b_role):
    """Name/system priors for interfaces that cannot slide however they are marked."""
    if a_system==b_system:
        sa,sb=name_stem(a_name),name_stem(b_name)
        if sa and sa==sb:return 'tied_by_partition'
        if a_role=='vascular' and b_role=='vascular':return 'tied_by_tree_continuity'
    return None


def classify(shared,markers,dilation,thresholds):
    """A marker claims a contact patch when it covers enough of that patch's cells."""
    pairs=sorted(shared)
    index={p:i for i,p in enumerate(pairs)}
    flat_pair=[];flat_key=[]
    for p in pairs:
        flat_pair.append(np.full(len(shared[p]),index[p],np.int32));flat_key.append(np.asarray(shared[p],np.int64))
    flat_pair=np.concatenate(flat_pair);flat_key=np.concatenate(flat_key)
    patch=np.array([len(shared[p]) for p in pairs],np.int64)
    slide_cover=np.zeros(len(pairs));tie_cover=np.zeros(len(pairs))
    slide_by=[None]*len(pairs);tie_by=[None]*len(pairs)
    ambiguous_cover=np.zeros(len(pairs));ambiguous_by=[None]*len(pairs)
    for marker in markers:
        keys=marker['dilated'][dilation]
        hit=flat_pair[np.isin(flat_key,keys,assume_unique=False)]
        if not len(hit):continue
        counted=np.bincount(hit,minlength=len(pairs))
        cover=counted/patch
        label=marker['interface_class']
        if label in SLIDING_CLASSES:target,by=slide_cover,slide_by
        elif label in TIED_CLASSES:target,by=tie_cover,tie_by
        elif label in AMBIGUOUS_CLASSES:target,by=ambiguous_cover,ambiguous_by
        else:continue
        better=cover>target
        for i in np.flatnonzero(better):by[i]=(marker['structure']['name'],label,float(cover[i]),int(marker['footprint_voxels']))
        marker.setdefault('_claims',{})[dilation]=int((cover>=0.5).sum())
        target[better]=cover[better]
    out=[]
    for i,p in enumerate(pairs):
        out.append({'pair':p,'contact_voxels':int(patch[i]),'slide_cover':float(slide_cover[i]),
            'tie_cover':float(tie_cover[i]),'ambiguous_cover':float(ambiguous_cover[i]),
            'slide_marker':slide_by[i],'tie_marker':tie_by[i],'ambiguous_marker':ambiguous_by[i]})
    return out


def coverage_table(rows,thresholds):
    table={}
    for t in thresholds:
        counter=Counter()
        for r in rows:
            slide=r['slide_cover']>=t;tie=r['tie_cover']>=t;amb=r['ambiguous_cover']>=t
            if slide and tie:counter['conflicted']+=1
            elif slide:counter['slide']+=1
            elif tie:counter['tied']+=1
            elif amb:counter['articular_cartilage_ambiguous']+=1
            else:counter['unknown']+=1
        total=sum(counter.values())
        table[str(t)]={**dict(counter),'total':total,
            'classified_fraction':float((total-counter['unknown'])/total) if total else 0.}
    return table


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default='data/derived/soft-body-coupling-assessment-v1')
    args=parser.parse_args()
    out=ROOT/args.out
    if out.exists() and any(out.iterdir()):raise SystemExit('Choose a fresh output directory: '+str(out))
    out.mkdir(parents=True,exist_ok=True)
    started=time.time()

    inputs={p:sha256(ROOT/p) for p in (CANONICAL,BODY,INVENTORY,SUBSTRATE_MANIFEST,CARTILAGE_GAP,SOURCE_INDEX,
        UNMODELLED,VOID_COMPONENTS,WHOLE_BODY_DOMAIN,PELVIS_DOMAIN,FOREARM,GARMENT_RECEIPT,SUPINE_QUADRATURE,
        CONTACT_SOURCE,FOUNDATION_SOURCE,'data/derived/entity-tet-ready-v1/manifest.json',
        'data/derived/muscle-tet-ready-v1/manifest.json')}
    canonical=json.loads((ROOT/CANONICAL).read_bytes())
    inventory=json.loads((ROOT/INVENTORY).read_bytes())
    index=json.loads((ROOT/SOURCE_INDEX).read_bytes())
    unmodelled=json.loads((ROOT/UNMODELLED).read_bytes())
    entities={e['id']:e for e in canonical['entities']}

    recorded=canonical['registrations']['z_anatomy']
    fit=LandmarkRegistration(np.array([l['source_rotated_m'] for l in recorded['landmarks']]),
                             np.array([l['target_m'] for l in recorded['landmarks']]))
    report=fit.report()
    registration={'reproduced_exactly':{k:bool(abs(report[k]-recorded[k])<=1e-12) for k in
                    ('affine_rms_m','fit_rms_m','fit_max_m','affine_determinant')},
        'fit_rms_m':report['fit_rms_m'],'fit_max_m':report['fit_max_m'],
        'held_out_rms_m':json.loads((ROOT/'data/derived/joint-substrate-assessment-v1/registration.json').read_bytes())['recorded']['held_out_rms_m'],
        'held_out_max_m':json.loads((ROOT/'data/derived/joint-substrate-assessment-v1/registration.json').read_bytes())['recorded']['held_out_max_m'],
        'meaning':'In-sample landmark residual is 1.05 mm; the deterministic five-fold held-out residual is 4.67 mm RMS and 18.6 mm max. Marker placement uncertainty, not tissue boundary uncertainty.'}
    if not all(registration['reproduced_exactly'].values()):raise ValueError('Recorded z_anatomy transform did not reproduce')

    bodies=load_bodies()
    lo=np.array([b['lo'] for b in bodies]);hi=np.array([b['hi'] for b in bodies])
    n=len(bodies)
    names=[entities[b['entity_id']]['name'] if b['entity_id'] in entities else b['entity_id'] for b in bodies]
    systems=[entities[b['entity_id']]['system'] if b['entity_id'] in entities else 'unknown' for b in bodies]
    roles=[entities[b['entity_id']]['role'] if b['entity_id'] in entities else 'unknown' for b in bodies]

    # --- broad phase ---------------------------------------------------------
    broad={'bodies':n,'unordered_pairs':n*(n-1)//2,'aabb':{}}
    overlap=None
    for margin in (0.,0.001,0.003,0.005):
        ov=aabb_pairs(lo,hi,margin)
        broad['aabb'][f'{margin:g}']={'margin_m':margin,'candidate_pairs':int(ov.sum()//2),
            'fraction_of_all_pairs':float(ov.sum()/(n*(n-1)))}
        if margin==0.:overlap=ov
    # empirical scaling of AABB candidate pairs with body count
    rng=np.random.default_rng(20260906);curve=[]
    for size in (100,200,400,800,1200,1800,n):
        counts=[]
        for _ in range(3 if size<n else 1):
            pick=rng.choice(n,size=size,replace=False) if size<n else np.arange(n)
            counts.append(int(aabb_pairs(lo[pick],hi[pick],0.).sum()//2))
        curve.append({'bodies':int(size),'mean_candidate_pairs':float(np.mean(counts)),
                      'pairs_per_body':float(np.mean(counts)/size)})
    broad['aabb_scaling_curve']=curve
    broad['aabb_scaling_note']=('Candidate pairs per body grows with body count: the atlas is not a uniform '
        'point cloud, so AABB overlap is superlinear in N over this range. Fit the two endpoints for the exponent.')
    a,b=curve[0],curve[-1]
    broad['aabb_empirical_exponent']=float(np.log(b['mean_candidate_pairs']/a['mean_candidate_pairs'])/np.log(b['bodies']/a['bodies']))

    # AABB tightness: how much of its own box does a body fill?
    extent=hi-lo;box=np.prod(np.maximum(extent,1e-9),axis=1)
    volume=np.array([b['signed_volume_m3'] if b['signed_volume_m3'] and b['signed_volume_m3']>0 else np.nan for b in bodies])
    fill=volume/box
    by_system=defaultdict(list)
    for f,s in zip(fill,systems):
        if np.isfinite(f):by_system[s].append(float(f))
    broad['aabb_fill_fraction']={'definition':'signed mesh volume divided by axis-aligned bounding box volume',
        'all_bodies':{'count':int(np.isfinite(fill).sum()),'median':float(np.nanmedian(fill)),
                      'p05':float(np.nanpercentile(fill,5)),'p95':float(np.nanpercentile(fill,95))},
        'by_system':{k:{'count':len(v),'median':float(np.median(v))} for k,v in sorted(by_system.items())}}

    # --- measured adjacency --------------------------------------------------
    grids={}
    for h in GRID_M:
        stats,shared=adjacency(bodies,h)
        stats['adjacent_pairs']=len(shared)
        grids[f'{h:g}']=stats
        if abs(h-ADJACENCY_GRID_M)<1e-12:contact=shared
    broad['surface_adjacency_by_grid']=grids
    adjacent=set(contact)
    candidate=set()
    ii,jj=np.nonzero(np.triu(overlap,1))
    for x,y in zip(ii.tolist(),jj.tolist()):candidate.add((x,y))
    broad['narrow_phase']={'grid_h_m':ADJACENCY_GRID_M,
        'aabb_candidate_pairs':len(candidate),'surface_adjacent_pairs':len(adjacent),
        'aabb_precision':float(len(adjacent&candidate)/len(candidate)),
        'adjacent_pairs_missed_by_zero_margin_aabb':len(adjacent-candidate),
        'meaning':'Precision is the fraction of zero-margin AABB candidates whose surfaces actually share a 4 mm cell.'}

    total_area=sum(b['surface_area_m2'] for b in bodies)
    total_vertices=sum(b['vertices'] for b in bodies)
    node_density=total_vertices/total_area
    contact_area=sum(len(v) for v in contact.values())*ADJACENCY_GRID_M**2
    broad['contact_extent']={'total_body_surface_area_m2':total_area,
        'total_surface_vertices':int(total_vertices),'native_surface_node_density_per_m2':node_density,
        'contact_patch_area_proxy_m2':contact_area,
        'contact_fraction_of_all_surface':contact_area/total_area,
        'contact_nodes_at_native_resolution':int(round(contact_area*node_density)),
        'meaning':'Contact patch area is 4 mm cells shared by two surfaces; node count is that area at the atlas native surface node density.'}

    # --- sliding classification ---------------------------------------------
    markers=register_substrate(inventory,index,fit,ADJACENCY_GRID_M,MARKER_DILATIONS)
    marker_classes=Counter(m['interface_class'] for m in markers)
    classification={}
    rows_at={}
    for dilation in MARKER_DILATIONS:
        rows=classify(contact,markers,dilation,COVERAGE_THRESHOLDS)
        rows_at[dilation]=rows
        classification[str(dilation)]={'dilation_cells':dilation,
            'tolerance_m':dilation*ADJACENCY_GRID_M,'coverage':coverage_table(rows,COVERAGE_THRESHOLDS)}

    chosen=rows_at[MARKER_DILATIONS[0]];threshold=0.5
    detail=[]
    for r in chosen:
        i,j=r['pair']
        slide=r['slide_cover']>=threshold;tie=r['tie_cover']>=threshold;amb=r['ambiguous_cover']>=threshold
        raw='conflicted' if slide and tie else 'slide' if slide else 'tied' if tie else 'articular_cartilage_ambiguous' if amb else 'unknown'
        prior=continuity_prior(names[i],names[j],systems[i],systems[j],roles[i],roles[j])
        # Precedence, stated once: a continuity prior beats any overlying marker; a
        # sliding and a tied marker on one patch are resolved to the larger coverage,
        # ties to slide, because a ligament beside a capsule is its own body and will
        # be coupled as one, whereas the enclosed articulation must not be welded.
        if prior is not None:label=prior
        elif raw=='conflicted':label='slide' if r['slide_cover']>=r['tie_cover'] else 'tied'
        else:label=raw
        detail.append({'a':bodies[i]['entity_id'],'b':bodies[j]['entity_id'],'a_name':names[i],'b_name':names[j],
            'a_system':systems[i],'b_system':systems[j],'a_role':roles[i],'b_role':roles[j],
            'contact_voxels':r['contact_voxels'],
            'contact_area_proxy_m2':r['contact_voxels']*ADJACENCY_GRID_M**2,
            'interface':label,'marker_label':raw,'continuity_prior':prior,
            'slide_cover':r['slide_cover'],'tie_cover':r['tie_cover'],
            'ambiguous_cover':r['ambiguous_cover'],
            'slide_marker':r['slide_marker'],'tie_marker':r['tie_marker'],'ambiguous_marker':r['ambiguous_marker']})
    dump_lines(out/'interface_classification.jsonl',detail)

    label_counts=Counter(d['interface'] for d in detail)
    raw_counts=Counter(d['marker_label'] for d in detail)
    slide_pairs=[d for d in detail if d['interface']=='slide']
    unknown_by_role=Counter('|'.join(sorted((d['a_role'],d['b_role']))) for d in detail if d['interface']=='unknown')
    slide_area=sum(d['contact_area_proxy_m2'] for d in slide_pairs)
    by_role=Counter(tuple(sorted((d['a_role'],d['b_role']))) for d in slide_pairs)
    sliding={'threshold_coverage':threshold,'dilation_cells':MARKER_DILATIONS[0],
        'tolerance_m':MARKER_DILATIONS[0]*ADJACENCY_GRID_M,
        'marker_structures':dict(marker_classes),
        'marker_class_meaning':{**{k:('slide: '+v) for k,v in SLIDING_CLASSES.items()},
            **{k:('tied: '+v) for k,v in TIED_CLASSES.items()},
            **{k:('ambiguous: '+v) for k,v in AMBIGUOUS_CLASSES.items()}},
        'labels_after_precedence':dict(label_counts),
        'labels_marker_evidence_only':dict(raw_counts),
        'precedence_rule':['continuity prior (authored partition of one structure, or two vessels of one tree) overrides any marker',
            'a patch carrying both a sliding and a tied marker takes the larger coverage, ties to slide',
            'otherwise the marker label stands'],
        'classified_fraction':float(1-label_counts['unknown']/sum(label_counts.values())),
        'unclassified_by_role':dict(unknown_by_role.most_common(15)),
        'marker_promiscuity':sorted([{'name':m['structure']['name'],'interface_class':m['interface_class'],
            'footprint_voxels':m['footprint_voxels'],
            'patches_claimed_at_half_coverage':int(m.get('_claims',{}).get(MARKER_DILATIONS[0],0))}
            for m in markers],key=lambda r:-r['patches_claimed_at_half_coverage'])[:25],
        'sliding_contact_area_proxy_m2':slide_area,
        'sliding_area_fraction_of_all_contact':slide_area/contact_area if contact_area else 0.,
        'sliding_pairs_by_role':{'|'.join(k):v for k,v in by_role.most_common(20)},
        'sensitivity':classification,
        'top_sliding_pairs':[{k:d[k] for k in ('a_name','b_name','contact_area_proxy_m2','slide_cover','slide_marker')}
            for d in sorted(slide_pairs,key=lambda d:-d['contact_area_proxy_m2'])[:30]]}
    dump(out/'sliding_interfaces.json',sliding)
    dump_lines(out/'substrate_landing.jsonl',[{'id':m['structure']['id'],'name':m['structure']['name'],
        'joint_class':m['structure']['joint_class'],'interface_class':m['interface_class'],
        'canonical_bounds_m':m['canonical_bounds_m'],'footprint_voxels':m['footprint_voxels'],
        'source_triangles':m['structure']['source_triangles'],
        'closed_oriented_manifold':m['structure']['closed_oriented_manifold']} for m in markers])

    # --- representation crossover -------------------------------------------
    radius=np.array([b['radius_proxy_m'] if b['radius_proxy_m'] else np.nan for b in bodies])
    finite=np.isfinite(radius)
    bins=(0.000125,0.00025,0.0005,0.001,0.002,0.004,0.008,0.016,0.032)
    tabulate=[]
    for t in bins:
        pick=finite&(radius<t)
        tabulate.append({'radius_proxy_below_m':t,'bodies':int(pick.sum()),
            'summed_volume_m3':float(np.nansum(volume[pick])),
            'summed_surface_area_m2':float(sum(b['surface_area_m2'] for b,p in zip(bodies,pick) if p)),
            'systems':dict(Counter(s for s,p in zip(systems,pick) if p).most_common(6))})
    def tet_count(vol,edge):
        return float(vol/(edge**3/6))  # 6 tets per cube of side edge

    def per_body_cost(mask,across):
        """Tets when each body is resolved with `across` elements over its own radius."""
        edge=np.where(mask,radius/max(across,1),np.nan)
        return float(np.nansum(np.where(mask,volume/(edge**3/6),0.)))

    def beam_cost(mask,per_radius):
        """Centreline elements for a circular tube of the same volume and radius."""
        length=np.where(mask,volume/(np.pi*np.maximum(radius,1e-9)**2),0.)
        return float(np.nansum(np.where(mask,np.ceil(length/(per_radius*np.maximum(radius,1e-9))),0.)))
    crossover={'definition':'radius_proxy_m = 2V/A of the repaired tet-ready surface; the same quantity the tet-ready builds record.',
        'measured_distribution':{'bodies_with_positive_volume':int(finite.sum()),
            'min_m':float(np.nanmin(radius)),'median_m':float(np.nanmedian(radius)),'max_m':float(np.nanmax(radius)),
            'cumulative_below':tabulate},
        'element_cost_below_threshold':[{'radius_proxy_below_m':t,
            'bodies':int((finite&(radius<t)).sum()),
            'tets_one_element_across_own_radius':per_body_cost(finite&(radius<t),1),
            'tets_two_elements_across_own_radius':per_body_cost(finite&(radius<t),2),
            'tets_three_elements_across_own_radius':per_body_cost(finite&(radius<t),3),
            'beam_elements_five_radii_long':beam_cost(finite&(radius<t),5),
            'volumetric_to_beam_element_ratio':(per_body_cost(finite&(radius<t),2)/beam_cost(finite&(radius<t),5)
                if beam_cost(finite&(radius<t),5)>0 else None)} for t in bins],
        'note':('A tube resolved with three elements across its radius needs edge <= r/3; the tet count then scales '
                'as L/r^2 per unit length, while a beam needs one element per L/r_axial. The crossover is set by the '
                'question, not by the geometry: a structure whose transverse stress field is not being asked for is '
                'a fibre.')}
    dump(out/'representation_crossover.json',crossover)

    # --- interstitium ---------------------------------------------------------
    components=[json.loads(l) for l in (ROOT/VOID_COMPONENTS).open()]
    largest=max(components,key=lambda c:c['volume_m3'])
    void=unmodelled['primary']
    interstitial={'source':UNMODELLED,'grid_h_m':unmodelled['primary_grid_h_m'],
        'interior_volume_m3':void['interior_volume_m3'],'void_volume_m3':void['void_volume_m3'],
        'void_fraction_of_interior':void['void_fraction_of_interior'],
        'void_components':void['void_components'],
        'largest_component':{'volume_L':largest['volume_L'],'voxels':largest['voxels'],
            'depth_below_skin_mm':largest['depth_below_skin_mm'],
            'nearest_entity_share':largest['nearest_entity_share']},
        'mesh_cost_if_filled':[{'edge_m':e,'tetrahedra':void['void_volume_m3']/(e**3/6),
            'nodes_approx':void['void_volume_m3']/(e**3)} for e in (0.008,0.004,0.002,0.001)],
        'entity_claims_outside_envelope_m3':void['entity_claims_outside_envelope_m3'],
        'overlap':void['overlap']}
    dump(out/'interstitium.json',interstitial)

    # --- materialization classes ---------------------------------------------
    whole=json.loads((ROOT/WHOLE_BODY_DOMAIN).read_bytes())
    pelvis=json.loads((ROOT/PELVIS_DOMAIN).read_bytes())
    forearm=json.loads((ROOT/FOREARM).read_bytes())

    def region_box(entity_names,pad):
        low=[];high=[]
        for name in entity_names:
            match=[e for e in canonical['entities'] if e['name']==name]
            if not match:raise ValueError('Unknown region anchor: '+name)
            low.append(np.asarray(match[0]['bounds_m']['min']));high.append(np.asarray(match[0]['bounds_m']['max']))
        return np.min(low,axis=0)-pad,np.max(high,axis=0)+pad

    def select(low,high):
        inside=np.flatnonzero(np.all(hi>=low,axis=1)&np.all(lo<=high,axis=1))
        return inside

    def inside_box(keys,low,high):
        k=np.asarray(keys,np.int64)
        idx=np.stack([((k>>40)&0xFFFFF)-512,((k>>20)&0xFFFFF)-512,(k&0xFFFFF)-512],axis=1)
        centre=(idx+.5)*ADJACENCY_GRID_M
        return int(np.sum(np.all(centre>=low,axis=1)&np.all(centre<=high,axis=1)))

    def region_report(label,low,high,description):
        pick=select(low,high)
        sub=set(int(i) for i in pick)
        pairs=[(p,v) for p,v in contact.items() if p[0] in sub and p[1] in sub]
        area=sum(inside_box(v,low,high) for _,v in pairs)*ADJACENCY_GRID_M**2
        ov=aabb_pairs(lo[pick],hi[pick],0.)
        return {'label':label,'description':description,
            'box_min_m':low.tolist(),'box_max_m':high.tolist(),
            'box_volume_L':float(np.prod(high-low)*1000),
            'bodies_intersecting':int(len(pick)),
            'aabb_candidate_pairs':int(ov.sum()//2),
            'measured_adjacent_pairs':len(pairs),
            'contact_patch_area_proxy_m2_clipped_to_box':area,
            'contact_nodes_at_native_resolution':int(round(area*node_density)),
            'summed_volume_of_intersecting_bodies_m3':float(np.nansum(volume[pick])),
            'bodies_by_system':dict(Counter(systems[i] for i in pick).most_common(8)),
            'bodies_with_radius_proxy_below_1mm':int(np.sum(finite[pick]&(radius[pick]<0.001)))}

    forearm_low,forearm_high=region_box(['left radius','left ulna'],0.02)
    implant_centre=(forearm_low+forearm_high)/2
    implant_low=implant_centre-0.015;implant_high=implant_centre+0.015

    coarse=np.flatnonzero(finite&(radius>=0.005))
    coarse_pairs=[(p,v) for p,v in contact.items() if p[0] in set(int(i) for i in coarse) and p[1] in set(int(i) for i in coarse)]
    coarse_area=sum(len(v) for _,v in coarse_pairs)*ADJACENCY_GRID_M**2
    coarse_ov=aabb_pairs(lo[coarse],hi[coarse],0.)
    classes={'a_whole_body_soft':{'label':'whole-body soft materialization, major structures only',
        'selection':'radius_proxy_m >= 0.005 (five millimetre characteristic half-thickness or thicker)',
        'bodies':int(len(coarse)),'aabb_candidate_pairs':int(coarse_ov.sum()//2),
        'measured_adjacent_pairs':len(coarse_pairs),
        'contact_patch_area_proxy_m2':coarse_area,
        'summed_body_volume_m3':float(np.nansum(volume[coarse])),
        'tets_at_edge_0.005m':float(np.nansum(volume[coarse])/(0.005**3/6)),
        'tets_at_edge_0.0025m':float(np.nansum(volume[coarse])/(0.0025**3/6)),
        'contact_nodes_at_5mm_surface_resolution':int(round(coarse_area/(0.005**2))),
        'existing_comparable_receipt':{'path':WHOLE_BODY_DOMAIN,'nodes':whole['mesh']['nodes'],
            'tetrahedra':whole['mesh']['tetrahedra'],'spacing_m':whole['spacing_m'],
            'entities_owning_no_cell':whole['entities_owning_no_cell'],
            'note':'That build is one shared-node voxel partition: every material label is welded to its neighbour and it silently drops 1,210 entities at 10 mm.'}},
      'b_regional_high_fidelity':{**region_report('left forearm segment',forearm_low,forearm_high,
            'bounds of left radius and left ulna padded by 20 mm'),
        'tets_at_edge_0.001m':None,'existing_comparable_receipt':{
            'path':FOREARM,'tetrahedra':len(forearm['geometry']['tetrahedra']),
            'nodes':len(forearm['geometry']['reference_positions_m']),
            'local_dimensions_m':forearm['geometry']['local_dimensions_m'],
            'note':'A single homogeneous 10x10x3 mm block, one body, no contact at all; whole_body_force_feedback is false in its own coupling record.'}},
      'c_local_implant':{**region_report('30 mm implant box in the forearm',implant_low,implant_high,
            '30 mm cube at the centre of the forearm segment box'),
        'existing_comparable_receipt':{'path':PELVIS_DOMAIN,'tetrahedra':pelvis['tetrahedra'],
            'vertices':pelvis['vertices'],'spacing_m':pelvis['spacing_m'],
            'interfaces':pelvis['interfaces']}}}
    for label,edge in (('b_regional_high_fidelity',0.001),):
        classes[label]['tets_at_edge_0.001m']=float(classes[label]['summed_volume_of_intersecting_bodies_m3']/(edge**3/6))
    classes['c_local_implant']['tets_at_edge_graded_50um_to_1mm']={
        'basis':'20 percent of the box volume at 50 um and the remainder at 1 mm; 6 tets per cube',
        'tetrahedra':float(0.2*classes['c_local_implant']['box_volume_L']*1e-3/(50e-6**3/6)
                           +0.8*classes['c_local_implant']['box_volume_L']*1e-3/(0.001**3/6))}
    composition={'what_exists':[
        {'evidence':BODY+' coupling_contract','finding':'Eight owners are declared with a direction and an explicit feedback flag. Three carry feedback_applied: false. The contract records who owns what state; it does not transport state.'},
        {'evidence':FOREARM+' coupling','finding':'whole_body_force_feedback is false and native_blood_storage_connected is false. The 98,304-tet forearm block is driven by a prescribed indenter, not by a body.'},
        {'evidence':PELVIS_DOMAIN+' interfaces','finding':'"shared nodes bond internal material labels; sliding interfaces unresolved" and "no global organ or pelvis/garment collision activated".'},
        {'evidence':PELVIS_DOMAIN+' mass_activation_contract','finding':'replaces_source_owners lists the canonical entities the fine domain supersedes, additive_mass_allowed is false and canonical_runtime_handoff_applied is false. This is the only existing machinery that names a coarse owner a fine domain replaces, and it is a mass ledger, not a state transfer.'},
        {'evidence':CONTACT_SOURCE+' step_coupled','finding':'Refuses to advance owners whose clocks differ by more than 1e-12 s and rejects a duplicated material owner. A common-clock multi-owner contract already exists; it just has no spatial hierarchy.'},
        {'evidence':WHOLE_BODY_DOMAIN,'finding':'The one existing whole-body volume mesh is a shared-node voxel partition at 10 mm: 71,861 nodes, 229,692 tets, 1,210 entities owning no cell. Every interface in it is welded.'}],
      'recommended_mechanism':{'name':'prescribed-boundary handoff with a returned residual receipt',
        'steps':['Run the coarse materialization and checkpoint it (DynamicTetrahedra.checkpoint already returns position/velocity/time).',
            'Cut the fine domain box out of the coarse mesh. Every coarse node on the cut surface becomes a Dirichlet driver for the fine domain: interpolate coarse displacement and velocity onto the fine boundary nodes each coarse step, hold them zero-order between coarse steps, and record the interpolation residual.',
            'Extend the existing mass_activation_contract so the fine domain declares replaces_source_owners over the coarse cells it covers, additive_mass_allowed false, and now also declares the boundary node set and the coarse checkpoint digest it was driven from.',
            'Return the fine domain traction integral over the cut surface as a receipt, and compare it against the coarse internal traction the cut removed. That difference is the one-way coupling error and it must be published, not absorbed.',
            'Only if that residual is small relative to the coarse load may the fine result be interpreted as the body response; otherwise it is a probe under prescribed boundary motion, which is what forearm-touch already honestly says it is.'],
        'why_this_one':'It reuses three things that already exist: the checkpoint/restore contract, the common-clock guard in step_coupled, and the mass_activation_contract ownership ledger. It adds no new solver.'},
      'what_would_make_it_unsound':[
        'Two-way feedback without a shared clock: the coarse and fine domains have different explicit stability limits (whole-body 10 mm vs pelvic 2 mm, whose recorded max_explicit_dt_s is 4.596e-05 s). Subcycling the fine domain inside one coarse step and pushing force back is an unconserved staggered scheme unless the exchanged impulse is audited the way resolve_node_triangle_contact audits its pair impulses.',
        'Double counting mass. The fine domain must replace, never add. additive_mass_allowed is already false in the pelvic contract; a composition mechanism that forgets it inflates body mass.',
        'A cut surface that crosses a sliding interface. Prescribing displacement across a bursa or fascial plane welds it by construction and destroys exactly the behaviour the per-entity representation was chosen for. Cut surfaces must be checked against the sliding classification in this artifact.',
        'Driving a fine domain from a coarse mesh that never resolved the fine structures: 1,210 of 2,403 entities own no cell at 10 mm, so a boundary condition read from that mesh carries no information about them.',
        'Reading a coarse displacement field as a measured boundary condition. It is a model output; the residual receipt is the only thing that makes the composition falsifiable.']}
    dump(out/'materialization_classes.json',{'classes':classes,
        'coupling_contract':json.loads((ROOT/BODY).read_bytes())['coupling_contract'],
        'composition':composition,
        'note':'Every existing regional domain records its own detachment. Nothing in the repo passes mechanical state between two domains.'})

    recommendation={'ranked_architectures':[
      {'rank':1,'name':'Per-entity bodies with a classified interface graph: tied interfaces as bilateral node-pair constraints, licensed sliding interfaces as declared contact pairs',
       'decisive_tradeoff':('Buys the sliding behaviour the representation was chosen for and costs one declared pair per sliding interface. At the measured 4 mm adjacency graph that is '
            +str(label_counts['slide'])+' sliding pairs out of '+str(len(detail))+' adjacent pairs, so roughly '
            +str(round(100*label_counts['slide']/len(detail),1))+' percent of interfaces need a contact solve and the rest can be welded or tied. '
            'The cost is that '+str(label_counts['unknown'])+' pairs have no marker evidence at all and must be given a default.'),
       'default_for_unknown':'Tie. A wrongly tied interface is stiff and visibly wrong; a wrongly freed interface lets tissue pass through tissue and destroys the conservation receipts.'},
      {'rank':2,'name':'Interstitial matrix as the tied medium: mesh the 32.117 L void conformally against every structure surface, share nodes with tied neighbours, keep licensed sliding interfaces as contact pairs',
       'decisive_tradeoff':('Removes every tied constraint by construction and gives the body a physically real filler that currently carries 46.1 percent of the interior volume. '
            'Costs 3.0e6 tets at 4 mm or 2.4e7 at 2 mm for the void alone, and requires a conforming tetrahedralization against 2,403 surfaces at once - a single global mesh generation whose failure mode is all-or-nothing. '
            'The 1,001 disconnected void components and the 0.492 L of entity claims outside the envelope are pre-existing defects it would have to absorb.')},
      {'rank':3,'name':'Embedded-fibre coupling for thin structures inside a coarser continuum',
       'decisive_tradeoff':('Measured: 876 of 2,403 bodies have radius_proxy below 1 mm and 1,455 below 2 mm. Resolving all of them volumetrically at two elements across their own radius costs 3.5e7 tets; '
            'as embedded beams at five radii per element they cost 4.7e4 elements, a ratio of about 743 to 1. The tradeoff is that a fibre has no transverse stress field and no surface, so it can neither be indented nor slide.')},
      {'rank':4,'name':'One monolithic shared-node mesh over the whole body',
       'decisive_tradeoff':('Cheapest and already built once at 10 mm (71,861 nodes, 229,692 tets), but it welds every interface. '
            'It is a clearly argued negative: it cannot express any of the '+str(label_counts['slide'])+' sliding interfaces this artifact identifies, and at 10 mm it drops 1,210 of 2,403 entities entirely.')},
      {'rank':5,'name':'Penalty contact on every adjacent pair with no classification',
       'decisive_tradeoff':('Needs no anatomy, but at '+str(len(detail))+' adjacent pairs and a measured '
            +str(round(contact_area,3))+' m2 of contact patch - 64.9 percent of all body surface area - it makes almost the entire body surface a contact problem, '
            'and penalty stiffness on a tied interface is a spring where anatomy has collagen. A clearly argued negative.')}],
      'biggest_technical_risk':{'name':'Interface classification error at scale, not solver cost',
        'statement':('The solver side is bounded and measurable: '+str(len(detail))+' adjacent pairs, '+str(round(contact_area,3))+' m2 of patch, '
            +str(int(round(contact_area*node_density)))+' contact nodes at native surface resolution. What is not bounded is being wrong about which interfaces slide. '
            'This artifact classifies '+str(round(100*(1-label_counts['unknown']/len(detail)),1))+' percent of adjacent pairs, and it does so through a registration whose held-out landmark residual is '
            '4.67 mm RMS and 18.6 mm maximum - larger than most bursae. A 2 mm bursa placed 5 mm from its true interface frees the wrong pair. '
            'Every wrongly freed interface is a place where tissue passes through tissue under load, and the existing kernel will report it only as unresolved edge contacts, not as an error.'),
        'mitigation':'Default unknown to tied, publish the per-marker promiscuity table with the classification, and require a per-interface receipt before any interface is freed.'}}
    dump(out/'recommendation.json',recommendation)

    # --- what exists ----------------------------------------------------------
    garment=json.loads((ROOT/GARMENT_RECEIPT).read_bytes())
    supine=json.loads((ROOT/SUPINE_QUADRATURE).read_bytes())
    existing={'nodal_contact_kernel':{'path':CONTACT_SOURCE,
        'formulation':'node against target boundary triangle; sequential zero-restitution impacts with finite target inertia; position projection is explicit and reported, not absorbed',
        'broad_phase':'per contacting node, one cKDTree ball query over the target body triangle centroids, radius = enclosing sphere + search distance + accumulated target motion',
        'friction':'Coulomb, static and kinetic, applied as a tangential impulse capped by the normal impulse',
        'integration':'symplectic Euler under a common clock with an explicit elastic stability limit per body',
        'conserved_receipts':['paired_impulse_residual_ns','angular_impulse_residual_nms','momentum_residual_ns',
            'numerical_energy_defect_j','contact_dissipation_j','kinetic_transfer_a_j + kinetic_transfer_b_j + dissipation_j'],
        'explicit_gaps':['no edge/edge contact (reported as unresolved_edge_contacts, never silently resolved)',
            'no self collision','no swept/continuous collision','no body/body broad phase: every pair must be declared',
            'contact is one-directional per NodeTriangleContact, so each interface needs two declared pairs']},
      'largest_exercised_contact':{'receipt':GARMENT_RECEIPT,'bodies':2,
        'owners':['pelvic tissue tet body','generated-shorts front panel cloth'],
        'panel_nodes':135,'contact_resolutions':garment['contact_resolutions'],
        'unresolved_edge_contact_candidates':garment['unresolved_edge_contact_candidates'],
        'wall_seconds':garment['wall_seconds'],'simulated_seconds':garment['seconds'],
        'maximum_pair_impulse_residual_ns':garment['maximum_pair_impulse_residual_ns'],
        'energy_audit_closure_j':garment['energy_audit_closure_j'],
        'statement':'Two deformable bodies is the largest body-body contact ever exercised in this repo.'},
      'surface_foundation':{'path':FOUNDATION_SOURCE,
        'formulation':'per-quadrature-point compressible neo-Hookean layer against a single rigid half plane, with tanh Coulomb plus viscous tangential drag and a dissipation guard that throws if the tangential law creates energy',
        'declared_point_cap':50000,'largest_exercised_points':21381,
        'largest_exercised_bodies':len(supine['bodies']),
        'receipt':SUPINE_QUADRATURE,
        'statement':'This is body-against-ground, not body-against-body: 22 RIGID OpenSim bodies each carrying skin quadrature against one plane. It contains no pair search and cannot express tissue-tissue sliding.'}}
    dump(out/'existing_machinery.json',existing)
    embree_lib=ROOT/'data/runtime/geometry/libigl-2.6.2/venv/lib/python3.10/site-packages/lib64/libembree4.a'
    try:
        sys.path.insert(0,str(ROOT/'data/runtime/geometry/libigl-2.6.2/venv/lib/python3.10/site-packages'))
        embree_symbols=None
        import subprocess
        embree_symbols=('rtcCollide' in subprocess.run(['nm','-g','--defined-only',str(embree_lib)],
            capture_output=True,text=True).stdout)
    except Exception:
        embree_symbols=None
    broad['broad_phase_options']=[
      {'option':'declare every pair by hand (today)','cost':'O(declared)','verdict':
       'What contact_dynamics.py requires. Fine for 2 bodies, impossible for '+str(len(adjacent))+'.'},
      {'option':'dense AABB overlap matrix (this script)','cost':'O(N^2) bytes and time; 2403 bodies is a 5.8 M entry boolean, 1.4 s',
       'measured_exponent_over_100_to_2403_bodies':broad['aabb_empirical_exponent'],
       'verdict':'Adequate to build a static interface graph once, at any atlas size the repo will plausibly reach. Not a per-step structure.'},
      {'option':'uniform spatial hash on surface samples (this script, 4 mm)','cost':
       str(int(grids[f'{ADJACENCY_GRID_M:g}']['surface_voxel_entries']))+' int64 entries, one sort, '
       +str(int(grids[f'{ADJACENCY_GRID_M:g}']['occupied_surface_voxels']))+' occupied cells, mean '
       +str(round(grids[f'{ADJACENCY_GRID_M:g}']['mean_bodies_per_occupied_surface_voxel'],3))+' bodies per cell',
       'verdict':('Directly gives the narrow-phase pair set with '+str(round(100*(1-len(adjacent&candidate)/len(candidate)),1))
        +' percent fewer false positives than AABB, and it is the structure a per-step broad phase should use. '
        'Mean occupancy near one means the hash is not degenerate at body scale.')},
      {'option':'Embree from the vendored libigl venv','embree_static_library':str(embree_lib.relative_to(ROOT)),
       'library_present':embree_lib.is_file(),'library_bytes':embree_lib.stat().st_size if embree_lib.is_file() else 0,
       'python_bindings_exposed':['EmbreeIntersector','ambient_occlusion','reorient_facets_raycast','shape_diameter_function'],
       'rtcCollide_defined_in_static_library':embree_symbols,
       'verdict':('A clearly argued negative for broad phase. The igl Python binding exposes a ray-cast intersector only: '
        'no scene-scene collision, no pair query, no BVH refit. Embree 4 itself has rtcCollide, but reaching it means '
        'linking libembree4.a from C++, and Embree BVHs are built for static or rebuilt geometry, which is the wrong '
        'shape for '+str(len(bodies))+' deforming bodies refit every step. Use it for the one thing it is unmatched at - '
        'signed distance and closest-point queries against a frozen reference surface - not for finding contact pairs.')}]
    dump(out/'broadphase.json',broad)

    summary={'schema':'ihm.soft-body-coupling-assessment.v1',
        'scope':'Read-only screening assessment. No promotion, no canonical write, no physics executed.',
        'wall_seconds':time.time()-started,
        'bodies':{'tet_proof_surfaces':n,'tet_ready':int(sum(b['tet_ready'] for b in bodies)),
            'tet_ready_with_manifold_vertex_links':int(sum(b['tet_ready'] and b['vertex_links_manifold'] for b in bodies)),
            'total_surface_faces':int(sum(b['faces'] for b in bodies)),
            'total_surface_vertices':int(total_vertices),'total_surface_area_m2':total_area},
        'registration':registration,
        'headline':{
            'unordered_pairs':n*(n-1)//2,
            'zero_margin_aabb_candidate_pairs':broad['aabb']['0']['candidate_pairs'],
            'measured_surface_adjacent_pairs_4mm':len(adjacent),
            'aabb_precision':broad['narrow_phase']['aabb_precision'],
            'sliding_pairs':label_counts['slide'],'tied_pairs':label_counts['tied'],
            'tied_by_partition_pairs':label_counts['tied_by_partition'],
            'tied_by_tree_continuity_pairs':label_counts['tied_by_tree_continuity'],
            'ambiguous_pairs':label_counts['articular_cartilage_ambiguous'],
            'unknown_pairs':label_counts['unknown'],
            'marker_conflicted_pairs_before_precedence':raw_counts['conflicted'],
            'classified_fraction':sliding['classified_fraction'],
            'contact_patch_area_proxy_m2':contact_area,
            'interstitial_void_L':void['void_volume_m3']*1000,
            'largest_void_component_L':largest['volume_L']},
        'verified':[ 'Every body count, face/vertex count, surface area and bounding box is read from the repaired tet-ready surfaces on disk.',
            'AABB candidate pair counts and the 8/4/2 mm surface co-occupancy adjacency are exact over those surfaces.',
            'The recorded z_anatomy affine + thin-plate-spline transform reproduces its recorded affine_rms_m, fit_rms_m, fit_max_m and affine_determinant to 1e-12.',
            'Joint-substrate class counts are read from the existing verified inventory.',
            'Contact machinery behaviour, its conserved receipts and the two-body exercise limit are read from source and from existing run receipts.',
            'Interstitial volumes and component statistics are read from the existing unmodelled-volume build.'],
        'inferred':[ 'Surface co-occupancy in a 4 mm cell is a proxy for anatomical adjacency; it neither proves contact nor excludes it.',
            'A joint-substrate marker covering a contact patch is evidence that the interface is that structure. Marker placement carries the 4.67 mm held-out registration RMS, so small bursae near several interfaces can be assigned to the wrong one.',
            'The slide/tied class map is an anatomical reading of each structure class, not a measured mechanical boundary condition.',
            'The continuity priors are name-stem and system rules, not geometry: they assert that an authored partition of one organ and two vessels of one tree cannot slide on each other.',
            'Conflict precedence is a stated modelling choice; 3,032 patches at 4 mm carried both a sliding and a tied marker.',
            'All tetrahedron and node counts for hypothetical materializations are volume/edge arithmetic at 6 tets per cube, not meshed results.',
            'Contact node counts scale a measured contact area by the atlas native surface node density; a real contact mesh will differ.'],
        'artifacts':['broadphase.json','sliding_interfaces.json','interface_classification.jsonl',
            'substrate_landing.jsonl','representation_crossover.json','interstitium.json',
            'materialization_classes.json','existing_machinery.json','recommendation.json']}
    dump(out/'summary.json',summary)
    dump(out/'manifest.json',{'schema':'ihm.soft-body-coupling-assessment-manifest.v1',
        'builder':'scripts/assess_soft_body_coupling.py',
        'builder_sha256':sha256(Path(__file__)),
        'inputs_sha256':inputs,
        'grid_m':list(GRID_M),'adjacency_grid_m':ADJACENCY_GRID_M,
        'marker_dilations_cells':list(MARKER_DILATIONS),'coverage_thresholds':list(COVERAGE_THRESHOLDS),
        'canonical_assets_modified':False,
        'artifacts_sha256':{p.name:sha256(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='manifest.json'}})
    print(json.dumps(summary['headline'],indent=2))


if __name__=='__main__':
    main()
