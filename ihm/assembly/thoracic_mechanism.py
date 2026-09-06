"""Executable source-only thoracic material mechanism, without recoil or drive.

Coordinates: torso linear/angular velocity plus26 internal coordinates. Two
ambiguous rib axes are constrained. Positive triangle quadrature represents
all debited material; no independent effective-mass point is introduced.
"""
import gzip,hashlib,json
from pathlib import Path
import numpy as np
from .thoracic_anatomy import hinge_motion
from .cervical_inertia import _topology


def closed_volume(vertices,faces):
    v=np.asarray(vertices,float);f=np.asarray(faces,int)
    if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all() or f.ndim!=2 or f.shape[1]!=3 or not len(f) or f.min()<0 or f.max()>=len(v):raise ValueError('Invalid cavity geometry')
    if not _topology(f)['closed_oriented_edge_manifold']:raise ValueError('Cavity must have closed oriented indexed topology')
    x=v-v.mean(0);a,b,c=np.moveaxis(x[f],1,0)
    volume=float(np.einsum('ij,ij->i',a,np.cross(b,c)).sum()/6)
    if volume<=0:raise ValueError('Nonpositive oriented cavity volume')
    gradient=np.zeros_like(v)
    for slot,value in enumerate((np.cross(b,c),np.cross(c,a),np.cross(a,b))):np.add.at(gradient,f[:,slot],value/6)
    return volume,gradient


def spatial_jacobian(position,internal):
    n=len(position);result=np.zeros((n,3,32));result[:,:,:3]=np.eye(3)
    result[:,:,3:6]=np.cross(np.eye(3)[None,:,:],position[:,None,:]).transpose(0,2,1)
    result[:,:,6:]=internal
    return result


class ThoracicMechanism:
    def __init__(self,manifest):
        self.path=Path(manifest);self.root=Path(__file__).resolve().parents[2]
        self.recipe=json.loads(self.path.read_bytes());r=self.recipe
        if r['native_activation_allowed'] is not False:raise ValueError('Research mechanism cannot authorize native activation')
        self.anatomy_path=self.root/r['anatomy_manifest']['path'];raw=self.anatomy_path.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=r['anatomy_manifest']['sha256']:raise ValueError('Anatomy recipe changed')
        self.anatomy=json.loads(raw);self.material={}
        for ident,entry in r['materials'].items():
            source=self.anatomy['entities'][ident];raw=(self.anatomy_path.parent/source['source_geometry']['retained_copy']).read_bytes()
            if hashlib.sha256(raw).hexdigest()!=source['source_geometry']['sha256']:raise ValueError('Material geometry changed')
            payload=json.loads(gzip.decompress(raw));v=np.asarray(payload['positions'],float).reshape(-1,3);faces=np.asarray(payload['indices'],int).reshape(-1,3)
            t=np.asarray(source['source_geometry_to_torso']);v=v@t[:3,:3].T+t[:3,3]
            raw=(self.path.parent/entry['map_file']).read_bytes()
            if hashlib.sha256(raw).hexdigest()!=entry['map_sha256']:raise ValueError('Material map changed')
            import io
            with np.load(io.BytesIO(raw),allow_pickle=False) as z:bindings={k:z[k] for k in z.files}
            tri=v[faces];area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
            self.material[ident]={'reference':v,'faces':faces,'triangle_mass':source['source_proxy_mass_kg']*area/area.sum(),'entry':entry,'bindings':bindings}
        self.hinges={h['coordinate_index']:h for h in self.anatomy['rib_hinges'].values()}
        self.locked=set(r['locked_internal_coordinates']);self.active=[i for i in range(32) if i<6 or i-6 not in self.locked]
        self.sternum_axis=np.array(self.anatomy['sternum_mode']['axis_in_torso']);self.descent_axis=np.array(r['diaphragm_axis_in_torso'])
        raw=(self.path.parent/r['cavity']['file']).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=r['cavity']['sha256']:raise ValueError('Cavity identity changed')
        import io
        with np.load(io.BytesIO(raw),allow_pickle=False) as z:self.cavity_data={k:z[k] for k in z.files}
        self.cavity_ids=r['cavity']['material_entity_order']
        self.reference_volume=closed_volume(self.cavity_data['reference'],self.cavity_data['faces'])[0]

    def coordinates(self,q):
        q=np.asarray(q,float)
        if q.shape!=(26,) or not np.isfinite(q).all():raise ValueError('Expected26finite coordinates')
        if any(abs(q[k])>1e-14 for k in self.locked):raise ValueError('Unresolved rib hinge is constrained')
        if np.max(np.abs(q[:24]))>.05 or abs(q[24])>.005 or abs(q[25])>.02:raise ValueError('Outside declared engineering geometry domain')
        return q

    def driver(self,points,indices,q):
        points=np.asarray(points,float);indices=np.asarray(indices,int)
        out=points.copy();jac=np.zeros((len(points),3,26))
        for index in np.unique(indices):
            mask=indices==index
            if index==-1 or index in self.locked:continue
            if index==24:out[mask]+=q[24]*self.sternum_axis;jac[mask,:,24]=self.sternum_axis
            elif index in self.hinges:
                h=self.hinges[index];out[mask],jac[mask,:,index]=hinge_motion(points[mask],h['origin_m'],h['axis'],q[index])
            else:raise ValueError('Unknown rigid material driver')
        return out,jac

    def material_state(self,ident,q,vertex_indices=None):
        q=self.coordinates(q);m=self.material[ident];b=m['bindings'];select=slice(None) if vertex_indices is None else np.asarray(vertex_indices,int)
        reference=m['reference'][select];kind=m['entry']['map_kind']
        if kind=='rigid':return self.driver(reference,np.full(len(reference),m['entry']['driver']),q)
        if kind!='moving_anchors':raise ValueError('Unknown material map')
        anchor,aj=self.driver(b['anchor_reference'],b['anchor_driver'],q)
        delta=anchor-b['anchor_reference'];indices=b['indices'][select];weights=b['weights'][select]
        out=reference.copy();jac=np.zeros((len(reference),3,26))
        for k in range(indices.shape[1]):
            out+=weights[:,k,None]*delta[indices[:,k]]
            jac+=weights[:,k,None,None]*aj[indices[:,k]]
        if 'descent_weights' in b:
            j=b['descent_weights'][select,None]*self.descent_axis
            out+=q[25]*j;jac[:,:,25]+=j
        return out,jac

    def kinetic(self,q,velocity=None):
        q=self.coordinates(q);speed=np.zeros(32) if velocity is None else np.asarray(velocity,float)
        if speed.shape!=(32,) or not np.isfinite(speed).all() or any(abs(speed[6+k])>1e-14 for k in self.locked):raise ValueError('Invalid generalized velocity')
        core=self.anatomy['replace_reduced_torso_with'];mass=core['mass_kg'];center=np.asarray(core['center_m']);inertia=np.asarray(core['inertia_kg_m2'])
        a=spatial_jacobian(center[None,:],np.zeros((1,3,26)))[0]
        matrix=mass*a.T@a;matrix[3:6,3:6]+=inertia
        totalmass=mass;first=mass*center;second=np.trace(inertia)/2*np.eye(3)-inertia+mass*np.outer(center,center)
        direct=.5*mass*np.sum((a@speed)**2)+.5*speed[3:6]@inertia@speed[3:6]
        bary=np.full((3,3),1/6);np.fill_diagonal(bary,2/3)
        for ident,m in self.material.items():
            x,j=self.material_state(ident,q)
            for start in range(0,len(m['faces']),1500):
                f=m['faces'][start:start+1500];weights=np.repeat(m['triangle_mass'][start:start+1500]/3,3)
                p=np.einsum('ab,tbc->tac',bary,x[f]).reshape(-1,3)
                pj=np.einsum('ab,tbcd->tacd',bary,j[f]).reshape(-1,3,26)
                jac=spatial_jacobian(p,pj);flat=jac.reshape(-1,32)
                matrix+=flat.T@(flat*np.repeat(weights,3)[:,None])
                totalmass+=weights.sum();first+=weights@p;second+=p.T@(weights[:,None]*p)
                velocities=np.einsum('ijk,k->ij',jac,speed);direct+=.5*float(weights@np.einsum('ij,ij->i',velocities,velocities))
        com=first/totalmass;second_com=second-totalmass*np.outer(com,com);combined_inertia=np.trace(second_com)*np.eye(3)-second_com
        active=matrix[np.ix_(self.active,self.active)];eigen=np.linalg.eigvalsh(active)
        if eigen[0]<=0:raise ValueError('Independent generalized mass is not positive definite; constrain geometry instead of adding stiffness')
        return {'mass_matrix':matrix,'active_indices':self.active,'active_eigenvalues':eigen,'mass_kg':float(totalmass),'center_m':com,
            'inertia_kg_m2':combined_inertia,'kinetic_energy_J':float(.5*speed@matrix@speed),'direct_material_energy_J':float(direct)}

    def cavity(self,q,*,pressure_pa=0.,gas_reference_volume_m3=None):
        q=self.coordinates(q);c=self.cavity_data;out=c['reference'].copy();jac=np.zeros((len(out),3,26))
        for code,ident in enumerate(self.cavity_ids):
            for slot in range(c['weights'].shape[1]):
                mask=c['entity_codes'][:,slot]==code
                if not mask.any():continue
                ids=c['vertex_indices'][mask,slot];x,j=self.material_state(ident,q,ids)
                delta=x-self.material[ident]['reference'][ids];w=c['weights'][mask,slot]
                out[mask]+=w[:,None]*delta;jac[mask]+=w[:,None,None]*j
        f=c['faces'];original=c['reference'][f];current=out[f]
        normal0=np.cross(original[:,1]-original[:,0],original[:,2]-original[:,0]);normal=np.cross(current[:,1]-current[:,0],current[:,2]-current[:,0])
        if (np.einsum('ij,ij->i',normal0,normal)<=0).any():raise ValueError('Cavity facet folds outside geometric domain')
        volume,gradient=closed_volume(out,f);jv=np.einsum('ij,ijk->k',gradient,jac)
        if not np.isfinite(pressure_pa):raise ValueError('Finite pressure required')
        forces=pressure_pa*gradient
        result={'geometric_volume_m3':volume,'volume_jacobian_m3_per_coordinate':jv,'pressure_generalized_force':pressure_pa*jv,
            'pressure_resultant_force_N':forces.sum(0),'pressure_resultant_moment_Nm':np.cross(out,forces).sum(0),
            'positions':out,'material_jacobian':jac,'pressure_nodal_forces':forces}
        if gas_reference_volume_m3 is not None:
            if not np.isfinite(gas_reference_volume_m3) or gas_reference_volume_m3<=0:raise ValueError('Explicit positive native reference required')
            adjusted=gas_reference_volume_m3+volume-self.reference_volume
            if adjusted<=0:raise ValueError('Offset cavity gives nonpositive gas-facing volume')
            result['gas_facing_volume_candidate_m3']=adjusted
            result['gas_offset_m3']=gas_reference_volume_m3-self.reference_volume
        return result

    def point_load(self,ident,vertex,q,force_N):
        x,j=self.material_state(ident,q,[vertex]);f=np.asarray(force_N,float)
        if f.shape!=(3,) or not np.isfinite(f).all():raise ValueError('Finite point force required')
        generalized=spatial_jacobian(x,j)[0].T@f
        return {'generalized_force':generalized,'fixed_parent_reaction_force_N':-f,
            'fixed_parent_reaction_moment_Nm':-np.cross(x[0],f),
            'reaction_scope':'Resultant required to hold the parent fixed under this applied point force; not an additional force on a free parent'}
