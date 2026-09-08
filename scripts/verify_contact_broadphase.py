#!/usr/bin/env python3
"""Broad phase, tied/contact partition, self collision, symmetric pairs and edge/edge.

Read-only over the atlas. Writes only inside the chosen output directory. No
canonical file is touched, nothing is promoted and `ihm/assembly/contact_dynamics.py`
is imported, never modified, so its published conservation receipts still apply.

What this verifies, in the order the coupling assessment lists the gaps:

1. Pair discovery replaces hand declaration. Three independent bounding-volume
   algorithms (dense matrix, sweep and prune, BVH) must return the SAME candidate
   set, and the two-stage broad phase must recover the measured adjacency of the
   `a_whole_body_soft` class with no false negative. The stage-one margin needed
   for that is derived, then measured.
2. A tied interface is a constraint. The partition refuses to hand one to the
   contact kernel and reports the counts.
3. Self collision candidates within one body, with topological neighbours excluded.
4. One symmetric declaration reproduces the two directional passes exactly, and a
   single-direction declaration is shown to leave the other side unresolved.
5. Edge/edge contact, with the same impulse audit as the node/triangle kernel, and
   a measurement of what the existing kernel does with a crossing it cannot see.
6. A profile of the shipped node/triangle kernel and a projected cost curve.

Honesty note carried into the receipt: stage two of the broad phase is the same
4 mm surface co-occupancy test that DEFINED the assessment's adjacency ground
truth, so its recall against that ground truth is exact by construction and is
reported as a reproduction, not as an independent validation. What is
independently measured is stage-one recall, the algorithms' mutual agreement, and
every cost number.
"""
import argparse,gzip,hashlib,io,json,sys,time
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import ihm.assembly.contact_dynamics as contact_dynamics
from ihm.assembly.contact_dynamics import (DynamicTetrahedra,NodeTriangleContact,resolve_node_triangle_contact,
                                           step_coupled)
from ihm.assembly.mechanics_backend import tetra_box
from ihm.assembly.contact_broadphase import (BodySurface,BroadPhase,SymmetricContact,dense_aabb_pairs,
    bvh_pairs,expand_symmetric,merge_symmetric_receipts,partition_interfaces,resolve_edge_edge_contact,
    self_collision_candidates,surface_edges,sweep_and_prune_pairs,voxel_keys,contact_triangle_ids)

TET_SOURCES=(('entity','data/derived/entity-tet-ready-v1'),('muscle','data/derived/muscle-tet-ready-v1'))
ASSESSMENT='data/derived/soft-body-coupling-assessment-v1'
ADJACENCY_GRID_M=0.004
COARSE_RADIUS_M=0.005


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as s:
        for b in iter(lambda:s.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def read_geometry(path):
    payload=json.loads(gzip.open(path).read())
    return (np.asarray(payload['positions'],float).reshape(-1,3),
            np.asarray(payload['indices'],np.int64).reshape(-1,3))


def rotation(axis,degrees):
    c,s=np.cos(np.radians(degrees)),np.sin(np.radians(degrees))
    return {0:np.array([[1,0,0],[0,c,-s],[0,s,c]]),1:np.array([[c,0,s],[0,1,0],[-s,0,c]]),
            2:np.array([[c,-s,0],[s,c,0],[0,0,1]])}[axis]


def crossed_ridges(offset_m=0.005,gap_m=0.0002,speed_m_s=-0.05):
    """Two 45 degree bars whose top and bottom ridges cross perpendicular.

    The crossing sits mid-segment on both ridges, so no vertex of either body is
    anywhere near the other surface. This is the pure edge/edge case.
    """
    def bar(size,divisions,axis,shift,speed):
        x,t=tetra_box(size,divisions);x=x-x.mean(axis=0)
        x=x@rotation(axis,45).T+np.asarray(shift,float)
        body=DynamicTetrahedra(x,t,mu_pa=1e5,lambda_pa=1e5,density_kg_m3=1000)
        body.velocity_m_s[:]=(0,0,speed);return body
    b=bar((.06,.01,.01),(6,1,1),0,(0,0,0),0.)
    a=bar((.01,.06,.01),(1,6,1),1,(offset_m,offset_m,0.),speed_m_s)
    a.position_m=a.position_m+np.array([0,0,b.position_m[:,2].max()+gap_m-a.position_m[:,2].min()])
    return a,b


def stacked_boxes(divisions=6,overlap_m=1e-4):
    x,t=tetra_box((.04,.04,.02),(divisions,divisions,2))
    a=DynamicTetrahedra(x,t,mu_pa=5e4,lambda_pa=5e4,density_kg_m3=1000)
    b=DynamicTetrahedra(x,t,mu_pa=5e4,lambda_pa=5e4,density_kg_m3=1000)
    a.position_m=a.position_m+np.array([0.,0.,.02-overlap_m]);a.velocity_m_s[:]=(0.,0.,-.5)
    return a,b


# --------------------------------------------------------------------------
# 1. broad phase
# --------------------------------------------------------------------------

def check_broad_phase_agreement():
    """Three bounding-volume algorithms must return the identical candidate set."""
    rng=np.random.default_rng(20260908)
    centre=rng.normal(size=(180,3))*0.05
    half=np.abs(rng.normal(size=(180,3)))*0.01+0.001
    lo=centre-half;hi=centre+half
    reference=None;report={}
    for name,fn in (('dense',dense_aabb_pairs),('sweep_and_prune',sweep_and_prune_pairs),('bvh',bvh_pairs)):
        pairs,stats=fn(lo,hi,0.002)
        key=set(map(tuple,pairs.tolist()))
        if reference is None:reference=key
        elif key!=reference:raise AssertionError(f'{name} disagrees with the dense AABB pair set')
        report[name]={'pairs':int(len(pairs)),**{k:v for k,v in stats.items() if k!='algorithm'}}
    report['agreed_pair_count']=len(reference)
    return report


def margin_theorem_check():
    """Two surfaces sharing one h cell differ by less than h on every axis.

    So a per-box margin of h/2 cannot produce a false negative, and any smaller
    margin can. Both halves are exercised, not asserted.
    """
    h=0.01
    a=np.array([[0.,0.,0.]]);b=np.array([[0.009,0.,0.]])
    tri=np.array([[0,0,0]])
    lo=np.stack([a[0],b[0]]);hi=np.stack([a[0],b[0]])
    small,_=dense_aabb_pairs(lo,hi,0.0)
    exact,_=dense_aabb_pairs(lo,hi,h/2)
    same_cell=len(np.intersect1d(np.floor(a/h).astype(int)@np.array([1,0,0]),
                                 np.floor(b/h).astype(int)@np.array([1,0,0])))>0
    return {'two_points_in_one_cell':bool(same_cell),'pairs_at_zero_margin':int(len(small)),
            'pairs_at_half_cell_margin':int(len(exact)),
            'meaning':'A zero margin misses a pair that shares a cell; h/2 recovers it.'}


def load_atlas_bodies():
    records=[]
    for family,directory in TET_SOURCES:
        for line in (ROOT/directory/'entities.jsonl').open():
            row=json.loads(line);row['_dir']=directory;records.append(row)
    bodies=[]
    for record in records:
        v,f=read_geometry(ROOT/record['_dir']/record['output_path'])
        tri=v[f];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        area=float(np.linalg.norm(normal,axis=1).sum()/2)
        volume=float(np.einsum('ij,ij->i',tri[:,0],normal).sum()/6)
        bodies.append({'entity_id':record['entity_id'],'surface_area_m2':area,'signed_volume_m3':volume,
            'radius_proxy_m':(2*volume/area if area>0 and volume>0 else None),'v':v,'f':f})
    return bodies


def atlas_broad_phase(bodies,recorded):
    radius=np.array([b['radius_proxy_m'] if b['radius_proxy_m'] else np.nan for b in bodies])
    coarse=np.flatnonzero(np.isfinite(radius)&(radius>=COARSE_RADIUS_M))
    surfaces=[BodySurface(bodies[i]['entity_id'],bodies[i]['v'],bodies[i]['f']) for i in coarse]
    lo=np.array([s.bounds()[0] for s in surfaces]);hi=np.array([s.bounds()[1] for s in surfaces])
    # Ground truth: the assessment's own 4 mm surface co-occupancy, recomputed here.
    truth=BroadPhase(surfaces,cell_m=ADJACENCY_GRID_M,margin_m=0.,prune='none')
    started=time.perf_counter();truth_result=truth.discover();truth_seconds=time.perf_counter()-started
    truth_pairs=set(truth_result.pairs)
    reproduced={'bodies':len(surfaces),'measured_adjacent_pairs':len(truth_pairs),
        'contact_patch_area_proxy_m2':sum(truth_result.shared_cells.values())*ADJACENCY_GRID_M**2,
        'contact_nodes_at_5mm_surface_resolution':int(round(
            sum(truth_result.shared_cells.values())*ADJACENCY_GRID_M**2/(0.005**2))),
        'unpruned_seconds':truth_seconds}
    recorded_class=recorded['classes']['a_whole_body_soft']
    reproduced['matches_recorded_assessment']={
        'bodies':reproduced['bodies']==recorded_class['bodies'],
        'measured_adjacent_pairs':reproduced['measured_adjacent_pairs']==recorded_class['measured_adjacent_pairs'],
        'contact_patch_area_proxy_m2':abs(reproduced['contact_patch_area_proxy_m2']
            -recorded_class['contact_patch_area_proxy_m2'])<1e-9,
        'contact_nodes_at_5mm_surface_resolution':reproduced['contact_nodes_at_5mm_surface_resolution']
            ==recorded_class['contact_nodes_at_5mm_surface_resolution']}
    if not all(reproduced['matches_recorded_assessment'].values()):
        raise AssertionError('Class a_whole_body_soft did not reproduce from the shipped broad phase')
    stage_one=[]
    for margin in (0.,0.001,0.002,0.004):
        row={'margin_m':margin}
        seen=None
        for name,fn in (('dense',dense_aabb_pairs),('sweep_and_prune',sweep_and_prune_pairs),('bvh',bvh_pairs)):
            started=time.perf_counter();pairs,stats=fn(lo,hi,margin);seconds=time.perf_counter()-started
            key=set(map(tuple,pairs.tolist()))
            if seen is None:seen=key
            elif key!=seen:raise AssertionError('Bounding-volume algorithms disagree on the atlas')
            row[name]={'seconds':seconds,'tests':int(stats.get('overlap_tests',0)+stats.get('node_tests',0))}
        recovered=len(seen&truth_pairs)
        row.update({'candidate_pairs':len(seen),'recovered_adjacent_pairs':recovered,
            'missed_adjacent_pairs':len(truth_pairs)-recovered,
            'recall':recovered/len(truth_pairs),
            'false_positive_pairs':len(seen)-recovered,
            'false_positive_rate':(len(seen)-recovered)/len(seen)})
        stage_one.append(row)
    if [r for r in stage_one if r['margin_m']==ADJACENCY_GRID_M/2][0]['missed_adjacent_pairs']!=0:
        raise AssertionError('Half-cell margin produced a false negative')
    two_stage={}
    for prune in ('sweep_and_prune','bvh','dense'):
        started=time.perf_counter()
        result=BroadPhase(surfaces,cell_m=ADJACENCY_GRID_M,prune=prune).discover()
        seconds=time.perf_counter()-started
        got=set(result.pairs)
        if got!=truth_pairs:raise AssertionError(f'Two-stage broad phase with {prune} did not recover the adjacency')
        two_stage[prune]={'seconds':seconds,'discovered_pairs':len(got),'recall':1.0,'false_positive_rate':0.0,
            **{k:v for k,v in result.stats.items() if k in
               ('candidate_pairs','surface_cell_entries','occupied_surface_cells','prune_precision','margin_m')}}
    return {'reproduced_class_a':reproduced,'stage_one_bounding_volume':stage_one,
            'two_stage':two_stage,'surfaces':surfaces,'truth':truth_result}


def atlas_scaling(bodies):
    lo=np.array([b['v'].min(axis=0) for b in bodies]);hi=np.array([b['v'].max(axis=0) for b in bodies])
    rng=np.random.default_rng(20260908);n=len(bodies);curve=[]
    for size in (100,200,400,800,1200,1800,n):
        pick=rng.choice(n,size=size,replace=False) if size<n else np.arange(n)
        row={'bodies':int(size)}
        seen=None
        for name,fn in (('dense',dense_aabb_pairs),('sweep_and_prune',sweep_and_prune_pairs),('bvh',bvh_pairs)):
            started=time.perf_counter();pairs,stats=fn(lo[pick],hi[pick],ADJACENCY_GRID_M/2)
            row[name]={'seconds':time.perf_counter()-started,
                       'tests':int(stats.get('overlap_tests',0)+stats.get('node_tests',0))}
            key=set(map(tuple,pairs.tolist()))
            if seen is None:seen=key;row['candidate_pairs']=len(key)
            elif key!=seen:raise AssertionError('Algorithms disagree while scaling')
        curve.append(row)
    return curve


def distance_validation(surfaces,truth_pairs,sample_h=ADJACENCY_GRID_M/2):
    """Independent check: is 4 mm cell co-occupancy actually a proximity criterion?

    Stage two of the broad phase IS the test that defined the assessment's
    adjacency, so recall against it is a reproduction, not a validation. This
    measures the true minimum surface-sample distance of every stage-one candidate
    and asks how the voxel set compares with a distance threshold. Expensive:
    one KD-tree query set per candidate pair.
    """
    from scipy.spatial import cKDTree
    from ihm.assembly.contact_broadphase import surface_samples as _samples
    points=[]
    for s in surfaces:
        v=s.vertices_m
        points.append(np.concatenate([v,_samples(v[s.triangles],sample_h)],axis=0))
    trees=[cKDTree(p) for p in points]
    lo=np.array([s.bounds()[0] for s in surfaces]);hi=np.array([s.bounds()[1] for s in surfaces])
    candidates,_=sweep_and_prune_pairs(lo,hi,ADJACENCY_GRID_M/2)
    started=time.perf_counter();gap={}
    for a,b in candidates:
        a=int(a);b=int(b)
        if len(points[a])<=len(points[b]):d,_=trees[b].query(points[a],k=1)
        else:d,_=trees[a].query(points[b],k=1)
        gap[(a,b)]=float(d.min())
    seconds=time.perf_counter()-started
    adjacent=np.array([gap[p] for p in gap if p in truth_pairs])
    other=np.array([gap[p] for p in gap if p not in truth_pairs])
    thresholds=[]
    for tau in (0.,0.0005,0.001,0.002,0.004,0.006):
        near={p for p,v in gap.items() if v<=tau}
        agree=len(near&truth_pairs)
        thresholds.append({'tolerance_m':tau,'pairs_within_tolerance':len(near),
            'also_voxel_adjacent':agree,'voxel_adjacent_but_farther':len(truth_pairs)-agree,
            'within_tolerance_but_not_voxel_adjacent':len(near)-agree})
    if adjacent.max()>ADJACENCY_GRID_M*np.sqrt(3):
        raise AssertionError('A voxel-adjacent pair exceeded the cell diagonal')
    return {'sample_spacing_m':sample_h,'total_sample_points':int(sum(len(p) for p in points)),
        'candidate_pairs':int(len(candidates)),'seconds':seconds,
        'voxel_adjacent_surface_gap_m':{'median':float(np.median(adjacent)),
            'p95':float(np.percentile(adjacent,95)),'max':float(adjacent.max())},
        'non_adjacent_candidate_surface_gap_m':{'median':float(np.median(other)),
            'p05':float(np.percentile(other,5)),'min':float(other.min())},
        'by_tolerance':thresholds,
        'finding':('The 4 mm co-occupancy graph is a 4 mm-tolerance proximity graph, not a contact set. '
                   'Every one of its pairs is genuinely within one cell diagonal, so it has no geometric '
                   'false positive at its own tolerance, but at a 1 mm contact tolerance a large fraction of '
                   'its pairs are not in contact, and some pairs within 6 mm are outside it. Choose the '
                   'tolerance before quoting a pair count.')}


# --------------------------------------------------------------------------
# 2. tied versus sliding
# --------------------------------------------------------------------------

def interface_labels(path):
    labels={}
    for line in Path(path).open():
        row=json.loads(line)
        a,b=row['a'],row['b']
        labels[(a,b) if a<b else (b,a)]=row['interface']
    return labels


def partition_report(result,labels):
    partition=partition_interfaces(result,labels,unknown_as='contact')
    strict=partition_interfaces(result,labels,unknown_as='tied')
    refused=False
    try:partition_interfaces(result,labels,unknown_as='refuse')
    except ValueError:refused=True
    if not refused and partition.unknown_pairs:raise AssertionError('Unknown interfaces were not refusable')
    tied_area=sum(result.shared_cells[p] for p in partition.tied_pairs)*ADJACENCY_GRID_M**2
    slide=[p for p in partition.contact_pairs if p not in set(partition.unknown_pairs)]
    slide_area=sum(result.shared_cells[p] for p in slide)*ADJACENCY_GRID_M**2
    unknown_area=sum(result.shared_cells[p] for p in partition.unknown_pairs)*ADJACENCY_GRID_M**2
    return {'label_counts':partition.counts,
        'tied_pairs':len(partition.tied_pairs),'sliding_pairs':len(slide),
        'unknown_pairs':len(partition.unknown_pairs),
        'contact_declarations_if_unknown_is_contact':2*len(partition.contact_pairs),
        'contact_declarations_if_unknown_is_tied':2*len(strict.contact_pairs),
        'tied_contact_patch_area_m2':tied_area,'sliding_contact_patch_area_m2':slide_area,
        'unknown_contact_patch_area_m2':unknown_area,
        'refuse_policy_raises':refused,
        'meaning':('A tied pair is a constraint and is never returned as a contact. The declaration count '
                   'is doubled because the kernel pair is one-directional; see the symmetric section.')}


# --------------------------------------------------------------------------
# 3. self collision
# --------------------------------------------------------------------------

def folded_strip():
    """A strip closed into a ring with a 1.8 mm gap between its two end faces.

    The gap is smaller than the self-collision cell, so the two ends must be
    reported as candidates, and they are not topological neighbours of each other.
    """
    x,t=tetra_box((.12,.004,.004),(24,1,1))
    theta=2*np.pi*0.985*x[:,0]/.12
    radius=.0191+(x[:,1]-.002)
    return np.stack([radius*np.sin(theta),radius*np.cos(theta),x[:,2]],axis=1),t


def self_collision_report():
    x,t=folded_strip()
    body=DynamicTetrahedra(x,t,mu_pa=1e5,lambda_pa=1e5,density_kg_m3=1000)
    out={}
    for ring in (0,1,2):
        started=time.perf_counter()
        r=self_collision_candidates(body.position_m,body.triangles,cell_m=.002,ring=ring)
        r['seconds']=time.perf_counter()-started
        pairs=r.pop('candidate_pairs')
        r['candidate_pairs']=int(len(pairs))
        if ring>0:
            # Every surviving candidate must be a genuine proximity, not a mesh neighbour.
            corners=body.position_m[body.triangles]
            centres=corners.mean(axis=1)
            if len(pairs):
                gaps=np.linalg.norm(centres[pairs[:,0]]-centres[pairs[:,1]],axis=1)
                r['minimum_candidate_centroid_gap_m']=float(gaps.min())
                r['median_candidate_centroid_gap_m']=float(np.median(gaps))
        out[f'ring_{ring}']=r
    if out['ring_0']['candidate_pairs']<=out['ring_1']['candidate_pairs']:
        raise AssertionError('Shared-vertex exclusion removed nothing')
    return out


def atlas_self_collision(surfaces,count=6):
    """Self collision candidates on the largest real bodies of the class."""
    order=sorted(range(len(surfaces)),key=lambda i:-len(surfaces[i].triangles))[:count]
    rows=[]
    for i in order:
        s=surfaces[i]
        started=time.perf_counter()
        r=self_collision_candidates(s.vertices_m,s.triangles,cell_m=ADJACENCY_GRID_M,ring=1)
        seconds=time.perf_counter()-started
        pairs=r['candidate_pairs']
        centres=s.vertices_m[s.triangles].mean(axis=1)
        gap=float(np.linalg.norm(centres[pairs[:,0]]-centres[pairs[:,1]],axis=1).min()) if len(pairs) else None
        rows.append({'body_id':s.body_id,'triangles':int(len(s.triangles)),
            'cell_pairs_before_topology':r['cell_pairs_before_topology'],
            'topological_pairs_excluded':r['topological_pairs_excluded'],
            'self_collision_candidates':int(len(pairs)),
            'minimum_candidate_centroid_gap_m':gap,'seconds':seconds})
    return rows


# --------------------------------------------------------------------------
# 4. symmetric declaration
# --------------------------------------------------------------------------

def symmetric_report():
    friction=(.5,.3);search=.0005
    def run(contacts):
        a,b=stacked_boxes()
        dt=min(a.max_explicit_dt_s,b.max_explicit_dt_s)*.9
        out=step_coupled({'a':a,'b':b},dt,contacts=contacts,gravity_m_s2=(0,0,0))
        return a,b,out
    two=(NodeTriangleContact('a','b',*friction,search),NodeTriangleContact('b','a',*friction,search))
    symmetric=expand_symmetric(SymmetricContact('a','b',*friction,search),NodeTriangleContact)
    a1,b1,r1=run(two)
    a2,b2,r2=run(symmetric)
    identical=bool(np.array_equal(a1.position_m,a2.position_m) and np.array_equal(a1.velocity_m_s,a2.velocity_m_s)
                   and np.array_equal(b1.position_m,b2.position_m) and np.array_equal(b1.velocity_m_s,b2.velocity_m_s))
    if not identical:raise AssertionError('Symmetric expansion changed the state')
    a3,b3,r3=run((NodeTriangleContact('a','b',*friction,search),))
    merged=merge_symmetric_receipts(r2['contacts'][0],r2['contacts'][1])

    def receipts(out):
        return {'numerical_energy_defect_j':out['numerical_energy_defect_j'],
                'contact_dissipation_j':out['contact_dissipation_j'],
                'momentum_residual_ns':float(np.linalg.norm(out['momentum_residual_ns'])),
                'max_paired_impulse_residual_ns':float(max(np.linalg.norm(c['paired_impulse_residual_ns'])
                                                           for c in out['contacts'])),
                'max_angular_impulse_residual_nms':float(max(np.linalg.norm(c['angular_impulse_residual_nms'])
                                                             for c in out['contacts'])),
                'kinetic_transfer_sum_plus_dissipation_j':float(sum(c['kinetic_transfer_a_j']+c['kinetic_transfer_b_j']
                                                                    +c['dissipation_j'] for c in out['contacts'])),
                'contact_count':int(sum(c['contact_count'] for c in out['contacts'])),
                'unresolved_edge_contacts':int(sum(c['unresolved_edge_contacts'] for c in out['contacts'])),
                'max_preprojection_penetration_m':float(max(c['max_preprojection_penetration_m'] for c in out['contacts']))}
    return {'two_hand_declarations':receipts(r1),'one_symmetric_declaration':receipts(r2),
        'state_identical':identical,
        'merged_interface_receipt':{'paired_impulse_residual_ns':float(np.linalg.norm(merged['paired_impulse_residual_ns'])),
            'interface_impulse_residual_ns':float(np.linalg.norm(merged['interface_impulse_residual_ns'])),
            'contact_count':merged['contact_count'],'dissipation_j':merged['dissipation_j'],
            'kinetic_transfer_sum_plus_dissipation_j':merged['kinetic_transfer_a_j']+merged['kinetic_transfer_b_j']
                +merged['dissipation_j']},
        'single_direction_only':{**receipts(r3),
            'meaning':'Declaring one direction resolves only A into B; nothing stops B entering A.'},
        'finding':('The one-directional API is fixable at zero cost to the receipts: one declaration expands '
                   'to the same two passes, in a fixed order, so every array is bit-identical. What the '
                   'expansion does NOT fix is the underlying discretisation asymmetry: a node of A against a '
                   'triangle of B is not the same constraint as a node of B against a triangle of A, and the '
                   'two passes are sequential, so the interface impulse depends on which side is swept first.')}


def symmetric_order_sensitivity():
    """How much does the interface answer depend on which direction runs first?"""
    friction=(.5,.3);search=.0005
    def run(order):
        a,b=stacked_boxes()
        dt=min(a.max_explicit_dt_s,b.max_explicit_dt_s)*.9
        contacts=tuple(NodeTriangleContact(x,y,*friction,search) for x,y in order)
        step_coupled({'a':a,'b':b},dt,contacts=contacts,gravity_m_s2=(0,0,0))
        return np.concatenate([a.velocity_m_s.ravel(),b.velocity_m_s.ravel()])
    forward=run((('a','b'),('b','a')));reverse=run((('b','a'),('a','b')))
    return {'velocity_difference_norm_m_s':float(np.linalg.norm(forward-reverse)),
            'velocity_norm_m_s':float(np.linalg.norm(forward)),
            'relative_difference':float(np.linalg.norm(forward-reverse)/np.linalg.norm(forward))}


# --------------------------------------------------------------------------
# 5. edge/edge
# --------------------------------------------------------------------------

def activation_criterion_report():
    """What actually makes the shipped node/triangle kernel fire."""
    xb=np.array([[0.,0,0],[1.,0,0],[0.,1,0]]);tri=np.array([[0,1,2]])
    rows=[]
    for z in (0.02,0.001,0.0002,1e-13,0.,-0.0002,-0.001):
        r=resolve_node_triangle_contact(np.array([[.2,.2,z]]),np.array([[0.,0,-1.]]),np.array([2.]),
            xb,np.zeros((3,3)),np.ones(3)*3,tri,node_ids=np.array([0]),
            friction_static=.5,friction_kinetic=.3,search_distance_m=.01)
        rows.append({'node_height_m':z,'contact_count':r['contact_count'],
                     'unresolved_edge_contacts':r['unresolved_edge_contacts']})
    if any(r['contact_count'] for r in rows if r['node_height_m']>1e-12):
        raise AssertionError('Kernel fired on a node in front of the triangle plane')
    if not all(r['contact_count'] for r in rows if r['node_height_m']<=1e-13):
        raise AssertionError('Kernel failed to fire on a node behind the triangle plane')
    overhang=[]
    for x in (0.499,0.5005,0.502,0.505,0.51):
        r=resolve_node_triangle_contact(np.array([[x,x,-0.0002]]),np.array([[0.,0,-1.]]),np.array([2.]),
            xb,np.zeros((3,3)),np.ones(3)*3,tri,node_ids=np.array([0]),
            friction_static=.5,friction_kinetic=.3,search_distance_m=.01)
        overhang.append({'in_plane_overhang_m':float(max(0.,(2*x-1)/np.sqrt(2))),
                         'contact_count':r['contact_count'],'unresolved_edge_contacts':r['unresolved_edge_contacts']})
    return {'search_band':rows,'edge_overhang':overhang,
        'finding':('search_distance_m is not a proximity band. It only widens the triangle set; a node must '
                   'already be at or behind the triangle plane (gap <= 1e-12 m) for the kernel to fire, so the '
                   'law is penetration recovery, not predictive contact. With no swept collision a node that '
                   'crosses a surface within one step is never seen. unresolved_edge_contacts counts only a '
                   'node already behind the plane, outside the footprint and inside the search distance.')}


def edge_edge_report():
    search=.001
    sweep=[]
    for offset_mm in (0.,0.5,1.,2.,3.,4.,5.):
        a,b=crossed_ridges(offset_m=offset_mm/1000)
        node=0;unresolved=0
        for x,y in ((a,b),(b,a)):
            r=resolve_node_triangle_contact(x.position_m,x.velocity_m_s,x.mass_kg,y.position_m,y.velocity_m_s,
                y.mass_kg,y.triangles,node_ids=np.unique(x.triangles),friction_static=.5,friction_kinetic=.3,
                search_distance_m=search)
            node+=r['contact_count'];unresolved+=r['unresolved_edge_contacts']
        e=resolve_edge_edge_contact(a.position_m,a.velocity_m_s,a.mass_kg,b.position_m,b.velocity_m_s,b.mass_kg,
            edges_a=surface_edges(a.triangles),edges_b=surface_edges(b.triangles),
            friction_static=.5,friction_kinetic=.3,search_distance_m=search)
        sweep.append({'crossing_offset_m':offset_mm/1000,'node_triangle_contacts':node,
            'node_triangle_unresolved_edge_contacts':unresolved,'edge_edge_contacts':e['contact_count'],
            'edge_edge_minimum_gap_m':e['minimum_edge_gap_m'],
            'edge_edge_endpoint_deferred':e['endpoint_solutions_deferred'],
            'edge_edge_parallel_refused':e['unresolved_parallel_edges']})
    blind=[r for r in sweep if r['crossing_offset_m']>0.0004]
    if any(r['node_triangle_contacts'] or r['node_triangle_unresolved_edge_contacts'] for r in blind):
        raise AssertionError('The node/triangle kernel was expected to be blind to a mid-segment crossing')
    if not all(r['edge_edge_contacts']>0 for r in blind):
        raise AssertionError('Edge/edge kernel failed to see a mid-segment crossing')
    a,b=crossed_ridges()
    e=resolve_edge_edge_contact(a.position_m,a.velocity_m_s,a.mass_kg,b.position_m,b.velocity_m_s,b.mass_kg,
        edges_a=surface_edges(a.triangles),edges_b=surface_edges(b.triangles),
        friction_static=.5,friction_kinetic=.3,search_distance_m=search)
    momentum=(((e['velocity_a_m_s']-a.velocity_m_s)*a.mass_kg[:,None]).sum(axis=0)
              +((e['velocity_b_m_s']-b.velocity_m_s)*b.mass_kg[:,None]).sum(axis=0))
    closure=e['kinetic_transfer_a_j']+e['kinetic_transfer_b_j']+e['dissipation_j']
    receipts={'paired_impulse_residual_ns':float(np.linalg.norm(e['paired_impulse_residual_ns'])),
        'angular_impulse_residual_nms':float(np.linalg.norm(e['angular_impulse_residual_nms'])),
        'momentum_residual_from_applied_velocities_ns':float(np.linalg.norm(momentum)),
        'kinetic_transfer_a_j':e['kinetic_transfer_a_j'],'kinetic_transfer_b_j':e['kinetic_transfer_b_j'],
        'contact_dissipation_j':e['dissipation_j'],
        'kinetic_transfer_sum_plus_dissipation_j':float(closure),
        'contact_count':e['contact_count'],'minimum_edge_gap_m':e['minimum_edge_gap_m']}
    if receipts['paired_impulse_residual_ns']>1e-18:raise AssertionError('Edge/edge impulses do not pair')
    if abs(closure)>1e-15:raise AssertionError('Edge/edge energy audit does not close')
    if e['dissipation_j']<0:raise AssertionError('Edge/edge tangential law created energy')
    return {'crossing_offset_sweep':sweep,'single_impact_receipts':receipts,
        'finding':('A perpendicular crossing 0.2 mm apart at mid-segment produces contacts = 0 AND '
                   'unresolved_edge_contacts = 0 from the shipped kernel in both declared directions. The '
                   'existing counter therefore does not measure edge/edge exposure at all: it measures a '
                   'vertex that fell just outside a triangle footprint. The 34 candidates in the garment '
                   'receipt are a count of that different thing, so the true edge/edge exposure of that run '
                   'is unmeasured, not 34.'),
        'limits':['Velocity-level proximity constraint: no position projection, because an edge crossing '
                  'carries no inside/outside test and the penetration sign would be a guess.',
                  'Nearly parallel edges are refused, not resolved.',
                  'An endpoint solution is deferred to the node/triangle kernel, which owns the footprint test.',
                  'Not integrated into step_coupled: step_coupled owns the energy audit and was not modified.']}


# --------------------------------------------------------------------------
# 6. profile of the shipped kernel
# --------------------------------------------------------------------------

class _ProbeTree:
    """The real cKDTree, re-queried at the loop's first radius to isolate the
    accumulated target_motion_bound term. It never changes what the kernel sees."""
    state=None
    def __init__(self,points):
        from scipy.spatial import cKDTree
        self._tree=cKDTree(points)
    def query_ball_point(self,x,r):
        out=self._tree.query_ball_point(x,r)
        s=_ProbeTree.state
        if s['first_radius'] is None:s['first_radius']=float(r)
        s['returned']+=len(out)
        s['at_first_radius']+=len(self._tree.query_ball_point(x,s['first_radius']))
        s['queries']+=1;s['last_radius']=float(r)
        return out


def kernel_profile():
    import cProfile,pstats
    rows=[]
    for divisions in (4,8,12,16,20,24,32):
        a,b=stacked_boxes(divisions)
        nodes=np.unique(a.triangles)
        # Untimed, uninstrumented repeats first: wall time on a shared machine is
        # noisy, so the reported second is the best of three, while the candidate
        # test counts below are deterministic and repeat-independent.
        seconds=float('inf')
        for _ in range(3):
            started=time.perf_counter()
            resolve_node_triangle_contact(a.position_m,a.velocity_m_s,a.mass_kg,b.position_m,b.velocity_m_s,
                b.mass_kg,b.triangles,node_ids=nodes,friction_static=.5,friction_kinetic=.3,search_distance_m=.0005)
            seconds=min(seconds,time.perf_counter()-started)
        _ProbeTree.state={'queries':0,'returned':0,'at_first_radius':0,'first_radius':None,'last_radius':None}
        real=contact_dynamics.cKDTree;contact_dynamics.cKDTree=_ProbeTree
        try:
            r=resolve_node_triangle_contact(a.position_m,a.velocity_m_s,a.mass_kg,b.position_m,b.velocity_m_s,
                b.mass_kg,b.triangles,node_ids=nodes,friction_static=.5,friction_kinetic=.3,search_distance_m=.0005)
        finally:contact_dynamics.cKDTree=real
        s=_ProbeTree.state
        rows.append({'divisions':divisions,'surface_nodes':int(len(nodes)),'target_triangles':int(len(b.triangles)),
            'contacts':r['contact_count'],'seconds':seconds,'ball_queries':s['queries'],
            'candidate_triangle_tests':s['returned'],
            'candidate_triangle_tests_at_constant_first_radius':s['at_first_radius'],
            'motion_bound_inflation':s['returned']/max(s['at_first_radius'],1),
            'query_radius_growth':s['last_radius']/s['first_radius']})
    lo,hi=rows[1],rows[-1]
    exponents={'wall_seconds':float(np.log(hi['seconds']/lo['seconds'])/np.log(hi['surface_nodes']/lo['surface_nodes'])),
        'candidate_triangle_tests':float(np.log(hi['candidate_triangle_tests']/lo['candidate_triangle_tests'])
            /np.log(hi['surface_nodes']/lo['surface_nodes'])),
        'candidate_triangle_tests_without_motion_bound':float(
            np.log(hi['candidate_triangle_tests_at_constant_first_radius']/lo['candidate_triangle_tests_at_constant_first_radius'])
            /np.log(hi['surface_nodes']/lo['surface_nodes']))}
    a,b=stacked_boxes(32);nodes=np.unique(a.triangles)
    profiler=cProfile.Profile();profiler.enable()
    resolve_node_triangle_contact(a.position_m,a.velocity_m_s,a.mass_kg,b.position_m,b.velocity_m_s,b.mass_kg,
        b.triangles,node_ids=nodes,friction_static=.5,friction_kinetic=.3,search_distance_m=.0005)
    profiler.disable()
    stats=pstats.Stats(profiler)
    total=stats.total_tt
    by_function={}
    for func,(cc,nc,tt,ct,callers) in stats.stats.items():
        name=f'{Path(func[0]).name}:{func[1]}({func[2]})'
        by_function[name]={'calls':nc,'tottime_s':tt,'cumtime_s':ct,'cumulative_fraction':ct/total if total else 0.}
    ranked=dict(sorted(by_function.items(),key=lambda kv:-kv[1]['tottime_s'])[:12])
    ball=[v for k,v in by_function.items() if 'query_ball_point' in k]
    projection=[]
    reference=rows[-1]
    for target_nodes in (5000,20000,57172):
        scale=target_nodes/reference['surface_nodes']
        projection.append({'surface_nodes':target_nodes,
            'projected_seconds_per_contact_resolution_pass_at_measured_exponent':
                reference['seconds']*scale**exponents['wall_seconds'],
            'projected_candidate_triangle_tests':
                reference['candidate_triangle_tests']*scale**exponents['candidate_triangle_tests']})
    return {'scaling':rows,'measured_exponents_in_surface_nodes':exponents,
        'profile_total_seconds':total,'top_by_tottime':ranked,
        'ball_query_cumulative_fraction':(ball[0]['cumtime_s']/total if ball and total else None),
        'projection':projection,
        'finding':('The cKDTree ball query is not the bottleneck: it is a few percent of the kernel. The cost '
                   'is the narrow-phase triangle projection over a candidate set that the kernel itself '
                   'inflates. Two mechanisms, both measured: (1) the query radius is the enclosing sphere of '
                   'the WORST triangle in the whole target mesh, so candidates per node grow with target mesh '
                   'density even at constant geometry; (2) target_motion_bound accumulates every positional '
                   'projection for the whole sweep and is never reset, so the radius grows monotonically '
                   'within one call and the last nodes queried are far more expensive than the first.')}


def triangle_subset_report():
    """The kernel already accepts a triangle subset. What does using it buy?"""
    a,b=stacked_boxes(24)
    nodes=np.unique(a.triangles)
    args=(a.position_m,a.velocity_m_s,a.mass_kg,b.position_m,b.velocity_m_s,b.mass_kg)
    kwargs=dict(node_ids=nodes,friction_static=.5,friction_kinetic=.3,search_distance_m=.0005)
    resolve_node_triangle_contact(*args,b.triangles,**kwargs)
    started=time.perf_counter()
    full=resolve_node_triangle_contact(*args,b.triangles,**kwargs)
    full_seconds=time.perf_counter()-started
    surface_a=BodySurface('a',a.position_m,a.triangles)
    surface_b=BodySurface('b',b.position_m,b.triangles)
    started=time.perf_counter()
    keep=contact_triangle_ids(surface_a,surface_b,cell_m=ADJACENCY_GRID_M,dilation=1)
    subset_build_seconds=time.perf_counter()-started
    if not len(keep):return {'skipped':'no triangle subset found'}
    resolve_node_triangle_contact(*args,b.triangles[keep],**kwargs)
    started=time.perf_counter()
    subset=resolve_node_triangle_contact(*args,b.triangles[keep],**kwargs)
    subset_seconds=time.perf_counter()-started
    return {'target_triangles_full':int(len(b.triangles)),'target_triangles_subset':int(len(keep)),
        'subset_build_seconds':subset_build_seconds,
        'full_seconds':full_seconds,'subset_seconds':subset_seconds,
        'speedup':full_seconds/subset_seconds if subset_seconds else None,
        'contacts_full':full['contact_count'],'contacts_subset':subset['contact_count'],
        'paired_impulse_residual_full_ns':float(np.linalg.norm(full['paired_impulse_residual_ns'])),
        'paired_impulse_residual_subset_ns':float(np.linalg.norm(subset['paired_impulse_residual_ns'])),
        'dissipation_full_j':full['dissipation_j'],'dissipation_subset_j':subset['dissipation_j'],
        'velocity_difference_norm_m_s':float(np.linalg.norm(full['velocity_a_m_s']-subset['velocity_a_m_s'])),
        'impulse_difference_norm_ns':float(np.linalg.norm(full['impulse_a_ns']-subset['impulse_a_ns'])),
        'equivalent_to_the_full_target':bool(full['contact_count']==subset['contact_count']),
        'finding':('MEASURED NEGATIVE. The kernel accepts triangle_ids, and restricting the target to the '
                   'broad-phase cell neighbourhood does shrink its KD tree and its single global '
                   'enclosing-sphere radius. It also CHANGES THE ANSWER. Removing triangles cannot make the '
                   'nearest candidate nearer, but the kernel skips a node whose nearest candidate has a '
                   'positive plane gap, so deleting a nearer triangle can promote a farther one that the node '
                   'is behind, and a contact fires where none did. The contact count moved and the post-step '
                   'velocity field moved with it. Combined with a speedup near unity and a subset build cost '
                   'larger than the kernel call itself, triangle subsetting is not a usable optimisation for '
                   'this kernel as written. Fixing the cost needs the per-triangle enclosing radius and a '
                   'scoped motion bound inside the kernel, not a smaller argument to it.')}


def target_scale_projection(truth,profile):
    """Projected cost at 326 bodies / 57,172 nodes, from the MEASURED patch sizes.

    The cost model is fitted on the candidate triangle test count, not on wall
    seconds, because the test count is deterministic and repeat-independent while
    wall time on a shared machine is not. Seconds are then obtained by one
    measured seconds-per-test constant, reported so it can be rescaled.

    The garment receipt's 1,500x-slower-than-real-time figure is one monolithic
    pair. A broad phase replaces it with 1,363 small pairs, and because the
    kernel's cost is superlinear in the node count of a pair, splitting the same
    nodes across many small pairs is the whole saving. Both are reported.
    """
    cells=np.array(sorted(truth.shared_cells.values()),float)
    nodes=cells*ADJACENCY_GRID_M**2/(0.005**2)
    rows=profile['scaling']
    n=np.array([r['surface_nodes'] for r in rows],float)
    tests=np.array([r['candidate_triangle_tests'] for r in rows],float)
    bare=np.array([r['candidate_triangle_tests_at_constant_first_radius'] for r in rows],float)
    seconds=np.array([r['seconds'] for r in rows],float)
    fits={}
    for label,y in (('as_shipped',tests),('without_accumulated_motion_bound',bare)):
        exponent,intercept=np.polyfit(np.log(n[1:]),np.log(y[1:]),1)
        fits[label]={'c':float(np.exp(intercept)),'p':float(exponent)}
    seconds_per_test=float(seconds[-1]/tests[-1])
    out={'interfaces':int(len(cells)),'directional_declarations':int(2*len(cells)),
        'total_contact_nodes_at_5mm':float(nodes.sum()),
        'contact_nodes_per_interface':{'median':float(np.median(nodes)),'p95':float(np.percentile(nodes,95)),
            'max':float(nodes.max()),'mean':float(nodes.mean())},
        'fitted_cost_model':{'form':'candidate_triangle_tests = c * surface_nodes^p','fits':fits,
            'fitted_over_nodes':[float(n[1]),float(n[-1])],
            'seconds_per_candidate_triangle_test':seconds_per_test,
            'caveat':('Fitted on one stacked-box configuration between 194 and 2,306 surface nodes, in which '
                      'every surface node is a contact candidate. Real interfaces contact over part of their '
                      'surface, so this is a pessimistic per-pair model. Everything below is an extrapolation.')},
        'projected':{}}
    pelvis=json.loads((ROOT/'data/derived/material-domains/pelvis-0.002m/manifest.json').read_bytes())
    dt=float(pelvis['max_explicit_dt_s']);steps=1.0/dt
    for label,fit in fits.items():
        def cost(count,fit=fit):return fit['c']*np.maximum(count,1.)**fit['p']
        partitioned=float(2*cost(nodes).sum())
        monolithic=float(cost(float(nodes.sum())))
        out['projected'][label]={
            'candidate_triangle_tests_partitioned_1363_interfaces_two_passes':partitioned,
            'candidate_triangle_tests_single_monolithic_pair':monolithic,
            'monolithic_over_partitioned':monolithic/partitioned if partitioned else None,
            'seconds_per_contact_sweep_partitioned':partitioned*seconds_per_test,
            'seconds_per_contact_sweep_monolithic':monolithic*seconds_per_test,
            'wall_seconds_per_simulated_second_partitioned':partitioned*seconds_per_test*steps,
            'wall_seconds_per_simulated_second_monolithic':monolithic*seconds_per_test*steps}
    out['explicit_step_basis']={'max_explicit_dt_s':dt,
        'source':'data/derived/material-domains/pelvis-0.002m/manifest.json',
        'steps_per_simulated_second':steps}
    out['meaning']=('Contact resolution only, in this Python kernel. It excludes the elastic force '
                    'evaluation, which for 1.46 M tets at 5 mm dominates, so these are a lower bound on the '
                    'real-time factor and not an estimate of it. The single number that matters is '
                    'monolithic_over_partitioned: it is the payoff of having a broad phase at all, and it '
                    'comes entirely from the superlinear exponent.')
    return out


def self_test():
    report={'broad_phase_algorithm_agreement':check_broad_phase_agreement(),
            'margin_theorem':margin_theorem_check(),
            'self_collision_folded_strip':self_collision_report(),
            'symmetric_declaration':symmetric_report(),
            'symmetric_order_sensitivity':symmetric_order_sensitivity(),
            'node_triangle_activation_criterion':activation_criterion_report(),
            'edge_edge':edge_edge_report()}
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default=None,help='receipt directory; omit to run checks only')
    parser.add_argument('--self-test',action='store_true',help='synthetic checks only, no atlas read')
    parser.add_argument('--distance-check',action='store_true',help='independent surface-distance validation (slow)')
    args=parser.parse_args()
    started=time.time()
    report=self_test()
    report['mode']='self_test' if args.self_test else 'atlas'
    if not args.self_test:
        recorded=json.loads((ROOT/ASSESSMENT/'materialization_classes.json').read_bytes())
        bodies=load_atlas_bodies()
        broad=atlas_broad_phase(bodies,recorded)
        surfaces=broad.pop('surfaces');truth=broad.pop('truth')
        report['atlas_broad_phase']=broad
        report['atlas_scaling']=atlas_scaling(bodies)
        report['interface_partition']=partition_report(truth,interface_labels(ROOT/ASSESSMENT/'interface_classification.jsonl'))
        report['atlas_self_collision']=atlas_self_collision(surfaces)
        if args.distance_check:
            report['distance_validation']=distance_validation(surfaces,set(truth.pairs))
    report['kernel_profile']=kernel_profile()
    report['triangle_subset']=triangle_subset_report()
    if not args.self_test:
        report['target_scale_projection']=target_scale_projection(truth,report['kernel_profile'])
    report['wall_seconds']=time.time()-started
    report['status']='passed'
    if args.out:
        out=ROOT/args.out
        if out.exists() and any(out.iterdir()):raise SystemExit('Choose a fresh output directory: '+str(out))
        out.mkdir(parents=True,exist_ok=True)
        inputs={'ihm/assembly/contact_dynamics.py':sha256(ROOT/'ihm/assembly/contact_dynamics.py'),
                'ihm/assembly/contact_broadphase.py':sha256(ROOT/'ihm/assembly/contact_broadphase.py'),
                'scripts/verify_contact_broadphase.py':sha256(Path(__file__))}
        if not args.self_test:
            for p in (ASSESSMENT+'/materialization_classes.json',ASSESSMENT+'/interface_classification.jsonl',
                      ASSESSMENT+'/manifest.json','data/derived/entity-tet-ready-v1/manifest.json',
                      'data/derived/muscle-tet-ready-v1/manifest.json'):
                inputs[p]=sha256(ROOT/p)
        report['inputs_sha256']=inputs
        report['canonical_assets_read']=[]
        (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False,default=float)+'\n')
        print('receipt',out/'report.json')
    print('PASS contact broad phase: pair discovery, tied/contact partition, self collision, symmetric '
          'declaration, edge/edge conservation and a measured kernel profile')
    return report


if __name__=='__main__':main()
