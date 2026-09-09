"""Turn a segment's own real SURFACES into contact geometry.

The upright plant contacts the floor through `fall_proxy_<body>` ContactSpheres
whose radius is INSCRIBED IN THE SEGMENT'S INERTIA ELLIPSOID and whose centre is
the segment's mass centre.  That is a simplification of exactly the kind the
programme direction forbids: a femur is not a ball at its centre of mass.

The replacement does not need a registration.  The scaled source model already
carries, per Body, the bone meshes that segment IS -- 81 of them over the 22
bodies, expressed in the segment's own frame, with the subject's per-segment
`scale_factors` beside them.  `femur_r` carries `r_femur.vtp`; `tibia_r` carries
the tibia AND the fibula; `hand_r` carries 27 carpals, metacarpals and phalanges.
Taking contact geometry from there costs zero registration error, because these
are the surfaces the segment frames were defined against.  (The BodyParts3D
anatomy binding is a DIFFERENT mapping, for display, and carries 24.7 mm RMS of
registration error; nothing here goes through it.)

Two things have to be done to those meshes before OpenSim will accept them as
`ContactMesh`:

1. **Bake the scale.**  `ContactMesh` has a filename and a transform and NO
   scale factors, while the source geometry is generic and every Body carries
   its own subject scale.  Pointing `ContactMesh` at `r_femur.vtp` would put an
   UNSCALED femur on a scaled subject -- silently, since nothing checks.  So the
   scale is baked into a written mesh file, one per (body, mesh).
2. **State what the surface is.**  Per mesh this records vertices, faces,
   watertightness, surface area, volume, and the ratio of that volume to its own
   convex hull's.  That last number is the concavity of the bone, and it is here
   because the question "does the solver silently take the convex hull?" is only
   answerable against a mesh that is measurably NOT convex.

The control that says the concavity number means something: a sphere and a cube
are convex, so both must print a volume/hull ratio of 1.0.  They are measured by
the same code path and written into the report as `controls`.
"""
from pathlib import Path
import argparse,hashlib,json,sys,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def read_vtp(path):
    """Vertices and triangles of an ASCII VTK PolyData file.

    Deliberately narrow: it accepts exactly the shape the OpenSim bone library
    ships (single Piece, ASCII, Polys with connectivity/offsets) and raises on
    anything else rather than guessing.  A geometry file this fails on is a file
    whose contact behaviour would not have been the file's.
    """
    root=ET.parse(path).getroot()
    pieces=root.findall('./PolyData/Piece')
    if len(pieces)!=1:raise ValueError('expected exactly one PolyData Piece: '+str(path))
    piece=pieces[0]
    for name in ('NumberOfVerts','NumberOfLines','NumberOfStrips'):
        if int(piece.get(name,'0')):raise ValueError('only polygons are supported: '+str(path))
    points=piece.find('./Points/DataArray')
    if points.get('format')!='ascii':raise ValueError('only ASCII VTP is supported: '+str(path))
    vertices=np.fromstring(points.text,sep=' ').reshape(-1,3)
    if len(vertices)!=int(piece.get('NumberOfPoints')):raise ValueError('point count disagrees with header: '+str(path))
    arrays={a.get('Name'):a for a in piece.findall('./Polys/DataArray')}
    if set(arrays)!={'connectivity','offsets'}:raise ValueError('unexpected Polys arrays: '+str(path))
    if any(a.get('format')!='ascii' for a in arrays.values()):raise ValueError('only ASCII VTP is supported: '+str(path))
    connectivity=np.fromstring(arrays['connectivity'].text,sep=' ').astype(np.int64)
    offsets=np.fromstring(arrays['offsets'].text,sep=' ').astype(np.int64)
    if len(offsets)!=int(piece.get('NumberOfPolys')):raise ValueError('polygon count disagrees with header: '+str(path))
    starts=np.concatenate(([0],offsets[:-1]))
    sizes=offsets-starts
    faces=[]
    for start,size in zip(starts,sizes):
        polygon=connectivity[start:start+size]
        if size<3:raise ValueError('degenerate polygon: '+str(path))
        # Fan triangulation, which is what SimTK's own PolygonalMesh ->
        # TriangleMesh conversion does; recorded so the face count in the report
        # is the face count the solver will put springs on.
        for k in range(1,size-1):faces.append([polygon[0],polygon[k],polygon[k+1]])
    faces=np.asarray(faces,dtype=np.int64)
    if faces.size and (faces.min()<0 or faces.max()>=len(vertices)):raise ValueError('face index out of range: '+str(path))
    return vertices,faces

def measure(vertices,faces):
    """Shape facts, including the one the convex-hull question turns on."""
    import trimesh
    mesh=trimesh.Trimesh(vertices=vertices,faces=faces,process=False)
    hull=mesh.convex_hull
    volume=float(abs(mesh.volume));hull_volume=float(abs(hull.volume))
    # Deepest point of the mesh surface below its own hull's surface: how far a
    # hull would have to move a bone's surface to get rid of its concavity.
    # The hull is convex, so its interior distance-to-boundary is exactly the
    # smallest of the supporting-plane distances -- no proximity tree needed,
    # and no approximation.  A convex mesh must give 0 here, which is the point
    # of the icosphere and box controls.
    normals=np.asarray(hull.face_normals)
    offsets=(normals*np.asarray(hull.triangles)[:,0,:]).sum(axis=1)
    slack=offsets[None,:]-np.asarray(mesh.vertices)@normals.T
    depth=float(slack.min(axis=1).max()) if len(mesh.vertices) else 0.
    edges=mesh.edges_sorted.reshape(-1,2)
    _,counts=np.unique(edges,axis=0,return_counts=True)
    return dict(vertices=int(len(mesh.vertices)),faces=int(len(mesh.faces)),
                surface_area_m2=float(mesh.area),volume_m3=volume,
                convex_hull_volume_m3=hull_volume,convex_hull_faces=int(len(hull.faces)),
                volume_over_hull_volume=volume/hull_volume if hull_volume>0 else float('nan'),
                maximum_depth_below_hull_m=depth,
                watertight=bool(mesh.is_watertight),
                euler_number=int(mesh.euler_number),
                boundary_edges=int((counts==1).sum()),
                nonmanifold_edges=int((counts>2).sum()),
                bounding_box_m=[float(v) for v in mesh.extents],
                centroid_local_m=[float(v) for v in mesh.centroid])

def simtk_precondition(vertices,faces):
    """Exactly what SimTK::ContactGeometry::TriangleMesh's constructor demands.

    Replicated here rather than inferred, from
    SimTKmath/Geometry/src/ContactGeometry_TriangleMesh.cpp: no repeated vertex
    within a face, no zero-area face, no two faces sharing the same DIRECTED
    edge, and equal forward/backward edge counts -- i.e. a closed, consistently
    oriented, edge-2-manifold surface.  A mesh that fails ANY of these does not
    become a contact geometry; the engine throws while building the model.

    This is a gate against a case whose answer is known: the icosphere and box
    controls must pass it, and they do.
    """
    forward=set();backward=set()
    for a,b,c in faces:
        if a==b or b==c or c==a:return 'repeated vertex within a face'
        if np.linalg.norm(np.cross(vertices[b]-vertices[a],vertices[c]-vertices[a]))<=0:return 'degenerate (zero-area) face'
        for u,v in ((a,b),(b,c),(c,a)):
            key=(u,v) if u<v else (v,u)
            side=forward if u<v else backward
            if key in side:return 'two faces share the same directed edge'
            side.add(key)
    # SimTK checks the COUNTS and then, per edge, that the opposite half-edge
    # actually exists.  Equal counts are not enough -- a mesh can have as many
    # unmatched forward as backward half-edges -- and testing only the counts is
    # how the first version of this gate passed meshes the engine then threw on.
    if forward!=backward:
        return 'not closed: %d forward and %d backward half-edges, %d unmatched'%(
            len(forward),len(backward),len(forward^backward))
    return None

def repair(vertices,faces):
    """Weld, drop degenerate/duplicate faces, unify winding, fill holes.

    Reported, never silent: the caller records the surface-area and volume the
    repair moved, because a repaired bone is a DIFFERENT surface from the one
    the model ships and the difference is contact error.
    """
    import trimesh
    mesh=trimesh.Trimesh(vertices=vertices,faces=faces,process=True,validate=True)
    mesh.update_faces(mesh.unique_faces());mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    trimesh.repair.fix_winding(mesh);trimesh.repair.fill_holes(mesh);trimesh.repair.fix_normals(mesh)
    return np.asarray(mesh.vertices),np.asarray(mesh.faces)


def cap_boundaries(vertices,faces,maximum_loops=64):
    """Close every open boundary loop with a fan to its own centroid.

    Cutting a surface into per-segment pieces leaves each piece open, and SimTK
    refuses an open mesh.  trimesh's `fill_holes` only closes small holes -- on
    these cuts it adds zero to sixteen triangles and leaves the loop open -- so
    the cap is built explicitly: chain the unmatched directed edges into loops
    and fan each loop to a new vertex at its centroid, with the winding that
    makes each cap face carry the REVERSE of its boundary half-edge.

    The caps are invented surface.  They sit at the joint, between two segments,
    which is where nothing outside the body can reach; the area they add is
    returned so it is reported rather than absorbed into the skin's own.
    """
    faces=[tuple(int(v) for v in f) for f in faces]
    directed={}
    for a,b,c in faces:
        for u,v in ((a,b),(b,c),(c,a)):directed[(u,v)]=directed.get((u,v),0)+1
    boundary=[e for e in directed if (e[1],e[0]) not in directed]
    if not boundary:return np.asarray(vertices),np.asarray(faces,dtype=np.int64),0.
    successors={}
    for u,v in boundary:successors.setdefault(u,[]).append(v)
    remaining=set(boundary);loops=[]
    while remaining:
        if len(loops)>maximum_loops:raise ValueError('boundary is not a small set of loops')
        start=next(iter(remaining))[0];loop=[];node=start
        while True:
            options=successors.get(node)
            if not options:raise ValueError('open boundary chain does not close')
            following=options.pop()
            if (node,following) not in remaining:raise ValueError('boundary edge visited twice')
            remaining.discard((node,following));loop.append((node,following));node=following
            if node==start:break
        loops.append(loop)
    vertices=list(np.asarray(vertices));added=0.
    for loop in loops:
        ring=np.asarray([vertices[u] for u,_ in loop])
        centre=len(vertices);vertices.append(ring.mean(axis=0))
        for u,v in loop:
            faces.append((v,u,centre))
            added+=float(np.linalg.norm(np.cross(vertices[u]-vertices[centre],vertices[v]-vertices[centre])))/2
    return np.asarray(vertices),np.asarray(faces,dtype=np.int64),added

def controls():
    """Cases whose answer is known: convex bodies must print a hull ratio of 1."""
    import trimesh
    out={}
    for name,mesh in (('icosphere_r0.1',trimesh.creation.icosphere(subdivisions=3,radius=.1)),
                      ('box_0.1',trimesh.creation.box(extents=(.1,.1,.1)))):
        out[name]=measure(np.asarray(mesh.vertices),np.asarray(mesh.faces))
    return out

def build(model_path,geometry_dirs,out_dir,scope,layer='bone'):
    model_path=Path(model_path).resolve();out_dir=Path(out_dir).resolve()
    tree=ET.parse(model_path)
    bodies=tree.getroot().find('.//BodySet/objects')
    if bodies is None:raise ValueError('model has no BodySet')
    meshes=out_dir/'meshes';meshes.mkdir(parents=True,exist_ok=True)
    # The source foot contacts are 12 anatomically placed spheres on calcn/toes
    # and are NOT inertia proxies; `proxy` replaces only what the inertia
    # ellipsoid invented, which is exactly the bodies the fall_proxy loop covers.
    foot_bodies={'talus_r','talus_l','calcn_r','calcn_l','toes_r','toes_l'}
    records=[];rejected=[]
    for body in bodies:
        if body.tag!='Body':continue
        name=body.get('name')
        if scope=='proxy' and name in foot_bodies:continue
        attached=body.find('attached_geometry')
        if attached is None:raise ValueError('body carries no attached geometry: '+name)
        for geometry in attached:
            if geometry.tag!='Mesh':continue
            frame=geometry.findtext('socket_frame')
            if frame!='..':raise ValueError('mesh is not attached to the body frame: '+name+'/'+geometry.get('name'))
            filename=geometry.findtext('mesh_file')
            scale=np.fromstring(geometry.findtext('scale_factors'),sep=' ')
            if scale.shape!=(3,) or not np.isfinite(scale).all() or (scale<=0).any():raise ValueError('invalid scale factors: '+filename)
            source=None
            for directory in geometry_dirs:
                candidate=Path(directory)/filename
                if candidate.is_file():source=candidate;break
            if source is None:raise ValueError('mesh file not found: '+filename)
            vertices,faces=read_vtp(source)
            unscaled=measure(vertices,faces)
            # SimTK will not accept a surface that is not a closed, consistently
            # oriented 2-manifold, so the mesh is tested BEFORE it is written and
            # repaired only if it has to be.  61 of the 81 bones OpenSim ships
            # with this model pass untouched; repair recovers 7 more; 13 cannot
            # be made admissible by welding, rewinding and hole filling, and
            # those are refused rather than quietly reshaped further.
            reason=simtk_precondition(vertices,faces);repaired=False
            if reason is not None:
                fixed_vertices,fixed_faces=repair(vertices,faces)
                fixed_reason=simtk_precondition(fixed_vertices,fixed_faces)
                if fixed_reason is None:vertices,faces,repaired=fixed_vertices,fixed_faces,True
                rejected.append(dict(body=name,mesh_file=filename,
                                     raw_reason=reason,repaired_reason=fixed_reason,
                                     admitted=fixed_reason is None,
                                     raw=unscaled,repaired=measure(fixed_vertices,fixed_faces)))
                if fixed_reason is not None:continue
            scaled_vertices=vertices*scale
            scaled=measure(scaled_vertices,faces)
            stem=name+'__'+Path(filename).stem
            target=meshes/(stem+'.obj')
            with target.open('w') as handle:
                handle.write('# '+stem+' from '+filename+' scaled by '+' '.join(map(repr,scale))+'\n')
                for v in scaled_vertices:handle.write('v %.9g %.9g %.9g\n'%tuple(v))
                for f in faces:handle.write('f %d %d %d\n'%(f[0]+1,f[1]+1,f[2]+1))
            records.append(dict(body=name,element=layer+'_'+stem,mesh_file=target.name,
                                source_geometry=str(source.relative_to(ROOT)),
                                source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                                scale_factors=[float(v) for v in scale],
                                written_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                                simtk_precondition_repaired=repaired,
                                unscaled=unscaled,scaled=scaled))
    if not records:raise ValueError('no bone meshes selected')
    report=dict(schema='ihm.segment-contact-meshes.v1',layer=layer,model=str(model_path),scope=scope,
                model_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(),
                geometry_dirs=[str(Path(d)) for d in geometry_dirs],
                bodies=sorted({r['body'] for r in records}),
                meshes=len(records),
                total_faces=sum(r['scaled']['faces'] for r in records),
                total_vertices=sum(r['scaled']['vertices'] for r in records),
                watertight_meshes=sum(1 for r in records if r['scaled']['watertight']),
                concave_meshes=sum(1 for r in records if r['scaled']['volume_over_hull_volume']<.999),
                faces_by_body={b:sum(r['scaled']['faces'] for r in records if r['body']==b) for b in sorted({r['body'] for r in records})},
                basis='Bone surfaces taken from the scaled source model\'s own attached_geometry, in the segment frame, with the subject scale_factors baked into the written file because ContactMesh has no scale property. No atlas registration is involved.',
                repaired_meshes=sum(1 for r in records if r['simtk_precondition_repaired']),
                refused_meshes=sum(1 for r in rejected if not r['admitted']),
                refused=[dict(body=r['body'],mesh_file=r['mesh_file'],reason=r['repaired_reason'],
                              faces=r['raw']['faces'],surface_area_m2=r['raw']['surface_area_m2'])
                         for r in rejected if not r['admitted']],
                simtk_precondition='SimTK::ContactGeometry::TriangleMesh requires a closed, consistently oriented, non-degenerate edge-2-manifold. A surface that fails it is not contact geometry at all -- the engine throws while building the model -- so a refused mesh means that bone is absent from contact, and the segment falls back to whatever else it carries.',
                controls=controls(),records=records,rejected=rejected)
    (out_dir/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--model',required=True)
    parser.add_argument('--geometry',action='append',default=[])
    parser.add_argument('--out',required=True)
    parser.add_argument('--scope',choices=('proxy','all'),default='all')
    parser.add_argument('--layer',choices=('bone','skin'),default='bone',
        help='What these surfaces ARE. Contact with the world is never bone against world: '
             'it is skin, over fat and muscle, over bone. A `bone` bundle is a collider against '
             'other bones and its own soft tissue, and a run that stands on one is standing on '
             'its skeleton -- which is a measurement, not the shipping contact path.')
    args=parser.parse_args()
    directories=args.geometry or [ROOT/'data/raw/anatomy/opensim-models/source/Models/Rajagopal/Geometry',
                                  ROOT/'data/raw/anatomy/opensim-models/source/Geometry']
    report=build(args.model,directories,args.out,args.scope,args.layer)
    print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
