"""Bounded exact exterior-face candidate scan; no root reassignment or native run."""
from pathlib import Path
import json,sys
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.source_skin import file_sha256,read_small_member

def scan():
    partition=json.loads((ROOT/'data/research/hair_native_partition.json').read_bytes());topology_path=ROOT/'data/research/engineered_skin_territories/materialization.json';raw=topology_path.read_bytes()
    if file_sha256(raw)!=partition['topology_evidence']['sha256']:raise ValueError('Topology identity changed')
    topology=json.loads(raw);source=next(r for r in topology['source_receipts'] if r['entity_id']=='body-bp3d-FJ2810');raw=(ROOT/source['path']).read_bytes()
    if len(raw)>6*1024*1024 or file_sha256(raw)!=source['sha256']:raise ValueError('Bounded source identity mismatch')
    x=np.asarray(read_small_member(raw,'positions',max_chars=12*1024*1024)).reshape(-1,3);faces=np.asarray(read_small_member(raw,'indices',max_chars=8*1024*1024)).reshape(-1,3)
    eligible=np.asarray(topology['contact_eligible_triangle_ids']);tri=x[faces[eligible]];centers=tri.mean(1);tree=cKDTree(centers);rows=[]
    for group_name,group in partition['groups'].items():
        strands=group['prepared']['strands'];points=np.asarray(strands['centerlines_m']).reshape(-1,3);offsets=strands['strand_offsets']
        for guide,(a,b) in enumerate(zip(offsets[:-1],offsets[1:])):
            tangent=points[a+1]-points[a];tangent/=np.linalg.norm(tangent)
            for node in range(a+2,b):
                _,near=tree.query(points[node],k=16)
                for candidate in near:
                    face=int(eligible[candidate])
                    if face==group['prepared']['source_face_indices'][guide]:continue
                    t=tri[candidate];e0=t[1]-t[0];e1=t[2]-t[0];n=np.cross(e0,e1);n/=np.linalg.norm(n);r=points[node]-t[0];gap=float(r@n)
                    beta12=np.linalg.solve(np.array([[e0@e0,e0@e1],[e0@e1,e1@e1]]),[r@e0,r@e1]);beta=np.r_[1-beta12.sum(),beta12]
                    lateral=abs(float(n@tangent))
                    if beta.min()>.02 and 0<gap<.004 and lateral<.8:
                        rows.append({'group':group_name,'guide_index':guide,'sample_id':group['prepared']['sample_ids'][guide],'node_offset':node-a,'source_face_index':face,'gap_m':gap,'normal_tangent_cosine_abs':lateral,'face_barycentric':beta.tolist(),'normal':n.tolist(),'positions_m':t.tolist(),'source_node_indices':faces[face].tolist()})
    registration=json.loads((ROOT/partition['input_registration_path']).read_bytes()) if 'input_registration_path' in partition else json.loads((ROOT/'data/derived/audits/cutaneous-factory-ylro2d66/body/mechanics/registration.json').read_bytes())
    group=registration['groups']['hand_r'];center=(np.asarray(group['bounds_min_m'])+group['bounds_max_m'])/2;hand=[]
    names=list(registration['groups']);low=np.array([registration['groups'][n]['bounds_min_m'] for n in names]);high=np.array([registration['groups'][n]['bounds_max_m'] for n in names])
    for candidate in np.argsort(np.linalg.norm(centers-center,axis=1))[:128]:
        vertices=tri[candidate];distance=np.linalg.norm(np.maximum(np.maximum(low[None]-vertices[:,None],vertices[:,None]-high[None]),0),axis=2)
        if all(names[i]=='hand_r' for i in np.argmin(distance,axis=1)):
            normal=np.cross(vertices[1]-vertices[0],vertices[2]-vertices[0]);normal/=np.linalg.norm(normal)
            hand.append({'source_face_index':int(eligible[candidate]),'source_node_indices':faces[eligible[candidate]].tolist(),'positions_m':vertices.tolist(),'normal':normal.tolist(),'owner':'hand_r','owner_basis':'Same fixed nearest bone-envelope prior on all3exact exterior nodes'})
        if len(hand)==8:break
    rows.sort(key=lambda r:r['gap_m']);return {'source_sha256':source['sha256'],'topology_sha256':file_sha256(topology_path),'engineered_pose_hand_faces':hand,'candidates':rows[:20],'candidate_count':len(rows),'native_run':False,'source_roots_unchanged':True}

if __name__=='__main__':
    report=scan();print(json.dumps(report,indent=2));(ROOT/'data/research/hair_collision_face_scan.json').write_text(json.dumps(report,indent=2)+'\n')
