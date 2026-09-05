"""Conservative atlas region priors; regions are not measured skin labels."""
import numpy as np
PRIORS={'forehead':(285.,84.1),'neck':(47.3,17.1),'chest':(28.4,6.4),'upper_arm':(45.9,14.5),'forearm':(40.3,14.9),'back':(24.7,5.9),'abdomen':(16.9,6.2),'thigh':(21.,8.),'calf':(15.6,3.6)}


def regional_density(points):
    """Source-frame x left/y superior/z anterior. Zero means excluded OR unresolved.

    Entire hands, feet, central groin, face except forehead, scalp and axilla
    are conservatively omitted: no claim these excluded surfaces are all glabrous.
    """
    p=np.asarray(points);x=np.abs(p[:,0]);y=p[:,1];z=p[:,2]
    region=np.full(len(p),'unresolved_excluded',dtype='<U32')
    for name,mask in [
        ('back',(x<.19)&(y>.04)&(y<.54)&(z<-.025)),
        ('abdomen',(x<.16)&(y>.015)&(y<.27)&(z>.025)),
        ('chest',(x<.18)&(y>=.27)&(y<.52)&(z>.025)),
        ('neck',(x<.07)&(y>.54)&(y<.62)),
        ('upper_arm',(x>.20)&(y>.32)&(y<.50)),
        ('forearm',(x>.235)&(y>.065)&(y<=.32)),
        ('thigh',(x>.045)&(y<-.10)&(y>-.40)),
        ('calf',(x>.04)&(y<=-.40)&(y>-.73)),
        ('forehead',(x<.068)&(y>.745)&(y<.79)&(z>.075))]:region[mask]=name
    density=np.zeros(len(p))
    for name,(mean,sd) in PRIORS.items():density[region==name]=mean
    return density,region


def attached_roots(deformed_vertices,source_faces,samples):
    """Evaluate MaterialPoint root coordinates after skin vertex deformation."""
    return np.einsum('ni,nij->nj',samples['barycentric'],np.asarray(deformed_vertices)[np.asarray(source_faces)[samples['face_index']]])


def display_indices(samples,count):
    """Stable subset order; display count never changes full samples or their IDs."""
    ids=samples['ids'].astype(np.uint64)
    score=(ids+np.uint64(0x9e3779b97f4a7c15))*np.uint64(0xbf58476d1ce4e5b9)
    return np.argsort(score,kind='stable')[:max(0,count)]


def outward_surface_mask(points,normals,regions):
    """Conservative outward-face prior for the thin shell, using regional body axes.

    This is not a measured outer epidermis segmentation. Ambiguous tangential
    faces are excluded; axis approximation is a documented source of omissions.
    """
    p=np.asarray(points);axis=p.copy();axis[:,0]=0;axis[:,2]=0
    for name,x in [('upper_arm',.22),('forearm',.28),('thigh',.10),('calf',.10)]:
        mask=regions==name;axis[mask,0]=np.sign(p[mask,0])*x
    delta=p-axis;delta/=np.maximum(np.linalg.norm(delta,axis=1)[:,None],1e-30)
    return np.einsum('ij,ij->i',delta,normals)>.35
