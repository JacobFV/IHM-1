"""Native-volume display constraint and its articulated force Jacobian.

No passive recoil, mass, gas store or independently integrated respiratory DOF.
"""
from copy import deepcopy
import math
import numpy as np
from .respiratory_feedback import RespiratoryLoadPort


def array(value,shape,label):
 a=np.asarray(value,float)
 if a.shape!=shape or not np.isfinite(a).all():raise ValueError('Invalid '+label)
 return a

def scalar(value,label,positive=False):
 if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or (positive and value<=0):raise ValueError('Invalid '+label)
 return float(value)

class EmbodiedRespiration:
 def __init__(self,payload,initial_lung_volume_ml):
  self.payload=deepcopy(payload);self.reference_m3=scalar(initial_lung_volume_ml,'initial lung volume',True)*1e-6
  port=RespiratoryLoadPort(payload);self.b=port.b.copy();self.bindings={x['entity_id']:deepcopy(x) for x in payload['bindings']}
  self.dimensions=np.array([payload['reference_dimensions_m'][x] for x in ('transverse','anteroposterior','height')],float)
  if not np.isfinite(self.dimensions).all() or (self.dimensions<=0).any():raise ValueError('Invalid reference thorax dimensions')
  for row in self.bindings.values():
   array(row['centroid_m'],(3,),'reference centroid')
   if row['kind'] not in ('rib','sternum','diaphragm','lung','posterior_support'):raise ValueError('Unsupported respiratory binding kind')
  self.skin=deepcopy(payload['skin_field']);self.skin_ids=set(self.skin['entity_ids'])
  self._skin_basis(array(self.skin['center_m'],(3,),'skin center'))

 def _coordinates(self,volume_ml):
  volume=scalar(volume_ml,'lung volume',True)*1e-6;delta=volume-self.reference_m3;q=self.b*delta
  if np.any(np.abs(q)>.25*np.array([self.dimensions[0]/2,self.dimensions[1],self.dimensions[2]])):raise ValueError('Respiratory display exceeds source small-deformation domain')
  return volume,q

 def _pose(self,ident,entities):
  if ident not in entities:raise ValueError('Missing articulated respiratory entity '+ident)
  state=entities[ident];r=array(state.get('rotation_matrix',np.eye(3)),(3,3),'rotation');f=array(state.get('deformation_gradient',np.eye(3)),(3,3),'mechanical deformation')
  if not np.allclose(r.T@r,np.eye(3),rtol=0,atol=1e-8) or np.linalg.det(r)<=0 or np.linalg.det(f)<=0:raise ValueError('Nonproper rotation or collapsed mechanical shape')
  current=array(state['centroid_m'],(3,),'current centroid')
  if ident in self.bindings:reference=array(self.bindings[ident]['centroid_m'],(3,),'reference centroid')
  else:reference=current-array(state['translation_m'],(3,),'skin translation')
  if 'translation_m' in state and not np.allclose(current,reference+array(state['translation_m'],(3,),'translation'),rtol=0,atol=1e-8):raise ValueError('Inconsistent articulated centroid/translation')
  return reference,current,r,f

 def _shape(self,ident,volume,q,f):
  binding=self.bindings[ident];kind=binding['kind'];translation=array(binding['translation_basis'],(3,3),'translation basis')@q;jtranslation=np.asarray(binding['translation_basis'])@self.b
  d=np.eye(3);derivative=np.zeros((3,3));width,depth,height=self.dimensions
  if kind=='lung':
   raw=np.array([1+2*q[0]/width,1+q[2]/height,1+q[1]/depth]);dr=np.array([2*self.b[0]/width,self.b[2]/height,self.b[1]/depth])
   if (raw<=0).any():raise ValueError('Collapsed respiratory lung shape')
   # Native volume controls total determinant, including mechanical shape.
   shape=raw*(volume/self.reference_m3/(np.linalg.det(f)*np.prod(raw)))**(1/3)
   ds=shape*(dr/raw+(1/volume-np.sum(dr/raw))/3)
   d=np.diag(shape);derivative=np.diag(ds)
  elif kind=='diaphragm':
   exponents=np.array([-1.,2.,-1.])/height;shape=np.exp(exponents*q[2]);d=np.diag(shape);derivative=np.diag(shape*exponents*self.b[2])
  return translation,jtranslation,d,derivative

 def _skin_basis(self,point):
  c=array(self.skin['center_m'],(3,),'skin center');low=array(self.skin['bounds_m']['min'],(3,),'skin lower bound');high=array(self.skin['bounds_m']['max'],(3,),'skin upper bound');radii=array(self.skin['reference_radii_m'],(2,),'skin radii');a,b=radii
  thorax=array(self.skin['thorax_y_m'],(2,),'thorax vertical band') if 'thorax_y_offsets_m' not in self.skin else c[1]+array(self.skin['thorax_y_offsets_m'],(2,),'thorax offsets')
  span=max(high[0]-c[0],c[0]-low[0])-a
  if (low>=high).any() or min(a,b,span)<=0 or not low[1]<thorax[0]<thorax[1]<high[1]:raise ValueError('Invalid skin field support')
  def smooth(x):
   x=max(0.,min(1.,x));return x*x*(3-2*x)
  x,y,z=point;weight=smooth((y-low[1])/(thorax[0]-low[1]))*smooth((high[1]-y)/(high[1]-thorax[1]))*(1-smooth((abs(x-c[0])-a)/span))
  basis=np.zeros((3,3));basis[0,0]=weight*max(-1.,min(1.,(x-c[0])/a));basis[2,1]=weight*max(0.,min(1.,(z-c[2]+b)/(2*b)))
  return basis

 def point_position(self,ident,reference_point_m,volume_ml,mechanical_entities):
  volume,q=self._coordinates(volume_ml);point=array(reference_point_m,(3,),'material point');reference,current,r,f=self._pose(ident,mechanical_entities)
  if ident in self.skin_ids:return current+r@f@(point-reference+self._skin_basis(point)@q)
  if ident not in self.bindings:raise ValueError('Unknown respiratory material point')
  translation,_,d,_=self._shape(ident,volume,q,f)
  return current+r@f@(translation+d@(point-reference))

 def geometry(self,current_lung_volume_ml,mechanical_entities,time_s=0.):
  time=scalar(time_s,'time')
  if time<0:raise ValueError('Negative respiratory clock')
  volume,q=self._coordinates(current_lung_volume_ml);result=deepcopy(mechanical_entities)
  for ident in self.bindings:
   reference,current,r,f=self._pose(ident,mechanical_entities);translation,_,d,_=self._shape(ident,volume,q,f);new_center=current+r@f@translation
   result[ident].update(centroid_m=new_center.tolist(),translation_m=(new_center-reference).tolist(),deformation_gradient=(f@d).tolist())
  for ident in self.skin_ids:self._pose(ident,mechanical_entities)
  skin=deepcopy(self.skin);skin.update(displacement_m=q.tolist(),lateral_expansion_m=float(q[0]),anterior_expansion_m=float(q[1]),diaphragm_descent_m=float(q[2]),coordinate_frame='canonical_reference_before_entity_transform')
  skin['entity_transforms']={ident:deepcopy(mechanical_entities[ident]) for ident in self.skin_ids if ident in mechanical_entities}
  return {'schema':'ihm.embodied-respiration.v1','time_s':time,'entities':result,'skin_field':skin,'displacement_m':q.tolist(),'native_volume_change_m3':volume-self.reference_m3,'lung_volume_ratio':volume/self.reference_m3,'independent_mass_or_recoil':False,'owner':'BioGears lung gas volume and passive chest/lung recoil','geometry_scope':'Articulated native-volume constraint; generic static modal shape, not a pleural/FEM solve'}

 def project_load(self,forces,mechanical_entities,current_lung_volume_ml):
  volume,q=self._coordinates(current_lung_volume_ml);ports=[];ignored=[];generalized=0.
  for force in forces:
   ident=force['id'];vector=array(force['force_n'],(3,),'force');point=array(force['point_m'],(3,),'world force point')
   moment=array(force.get('moment_nm',[0.,0.,0.]),(3,),'moment')
   if ident not in self.bindings and ident not in self.skin_ids:ignored.append(ident);continue
   reference,current,r,f=self._pose(ident,mechanical_entities);a=r@f
   if ident in self.skin_ids:
    if 'moment_nm' in force:raise ValueError('Skin strain work requires resolved point tractions, not a resultant wrench')
    local=reference+np.linalg.solve(a,point-current);material=local.copy()
    for _ in range(20):
     residual=material+self._skin_basis(material)@q-local
     if np.linalg.norm(residual)<1e-11:break
     h=1e-6;jac=np.eye(3)
     for k in range(3):
      step=np.zeros(3);step[k]=h;jac[:,k]+=(self._skin_basis(material+step)-self._skin_basis(material-step))@q/(2*h)
     if np.linalg.det(jac)<=0:raise ValueError('Folded local skin display map')
     material-=np.linalg.solve(jac,residual)
    else:raise ValueError('Skin force material-point inversion did not converge')
    j=a@self._skin_basis(material)@self.b
   else:
    translation,jt,d,dd=self._shape(ident,volume,q,f)
    if self.bindings[ident]['kind'] in ('lung','diaphragm') and 'moment_nm' in force:raise ValueError('Deforming tissue strain work requires resolved point tractions')
    relative=np.linalg.solve(d,np.linalg.solve(a,point-current)-translation)
    j=a@(jt+dd@relative)
   # Existing rib/sternum modes have no virtual angular motion. A resultant
   # moment is retained, but contributes zero for those translation-only modes.
   contribution=float(vector@j);generalized+=contribution
   ports.append({'id':ident,'point_m':point.tolist(),'force_n':vector.tolist(),'moment_nm':moment.tolist(),'point_jacobian_m_per_m3':j.tolist(),'angular_jacobian_rad_per_m3':[0.,0.,0.],'generalized_force_pa':contribution})
  return {'external_pressure_pa':-generalized,'force_ports':ports,'ignored_nonrespiratory_ids':ignored,'work_convention':'At fixed articulation, external virtual work = -pressure * d(native lung volume)','passive_recoil_added':False,'limitations':['Resultant wrenches do not determine distributed deforming-tissue strain work','Rigid respiratory bindings have translation-only modes, not rib rotation','Force points must lie in the current composed display geometry; material inversion follows that map','Contact-body to respiratory-entity assignment must be supplied explicitly']}
