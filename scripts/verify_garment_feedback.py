"""Tiny finite-body coupling fixture: no native launch, no full body trajectory."""
from pathlib import Path
import copy,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.garment_feedback import GarmentFeedback,wrench_point_loads
from ihm.assembly.garment_body_ports import BarycentricAttachment
from ihm.assembly.clothing import Cloth

class Registration:
    bodies=['body'];basis=np.eye(3);global_map=np.eye(4)
    def transforms(self,state):return {'body':np.array(state['bodies']['body']['transform_ground'])}
class Native:
    def __init__(self):self.x=np.zeros(3);self.v=np.zeros(3);self.time=0.;self.mass=1.;self.tokens=0
    def snapshot(self):
        t=np.eye(4);t[:3,3]=self.x
        return {'time_s':self.time,'gravity_m_s2':[0,0,0],'bodies':{'body':{'transform_ground':t.tolist(),'origin_velocity_m_s':self.v.tolist(),'angular_velocity_rad_s':[0,0,0]}}}
    def checkpoint(self):self.tokens+=1;return self.snapshot()
    def restore(self,s):self.x=np.array(s['bodies']['body']['transform_ground'])[:3,3];self.v=np.array(s['bodies']['body']['origin_velocity_m_s']);self.time=s['time_s']
    def release(self,s):self.tokens-=1
    def advance(self,dt,forces,actuation):
        f=np.sum([p['force_n'] for p in forces],axis=0) if forces else np.zeros(3);a=f/self.mass
        self.x=self.x+self.v*dt+.5*a*dt*dt;self.v=self.v+a*dt;self.time+=dt;return self.snapshot()
def fixture():
    b=np.array([[-1.,-1,0],[1,-1,0],[0,1,0]]);tri=np.array([[0,1,2]])
    x=np.array([[-.1,0,.03],[.1,0,.03],[0,.1,.03]]);cloth=Cloth(x,tri,areal_density_kg_m2=1.,edge_stiffness_n_m=.1)
    attachment=BarycentricAttachment(0,(0,1,2),(.3,.2,.5),0.,1.)
    return Native(),GarmentFeedback({'patch':cloth},b,tri,np.zeros(3,int),Registration(),{'patch':[attachment]},friction_static=.4,friction_kinetic=.3)
def main():
    loads=wrench_point_loads('body',[.3,-.2,.7],[2,3,4],[.3,-.8,.2]);origin=np.array([.3,-.2,.7])
    assert np.allclose(np.sum([p['force_n'] for p in loads],axis=0),[2,3,4])
    assert np.allclose(np.sum([np.cross(np.array(p['point_m'])-origin,p['force_n']) for p in loads],axis=0),[.3,-.8,.2])
    native,garment=fixture();result,report=garment.advance(native,.001,[],{},tolerance_m=1e-11)
    cloth=garment.garments['patch'];momentum=native.mass*native.v+(cloth.mass_kg[:,None]*cloth.velocity_m_s).sum(0)
    assert np.linalg.norm(momentum)<1e-14 and np.linalg.norm(native.v)>1e-6
    assert report['boundary_position_residual_m']<=1e-11 and abs(cloth.time_s-native.time)<1e-14 and native.tokens==0
    # Numerical power mismatch must be exposed, not silently treated as heat.
    assert np.isfinite(report['interface_work_quadrature_defect_j'])
    native,garment=fixture();before=garment.checkpoint();initial=native.snapshot()
    try:garment.advance(native,.001,[],{},tolerance_m=1e-20,max_iterations=1)
    except RuntimeError:pass
    else:raise AssertionError('Nonconverged coupling accepted')
    assert native.snapshot()==initial and native.tokens==0
    assert np.array_equal(garment.garments['patch'].position_m,before['garments']['patch']['position_m'])
    assert garment.garments['patch'].time_s==0
    print('PASS exact wrench reduction, finite-body garment reaction and total momentum, common clock, convergence receipt, atomic rollback')
if __name__=='__main__':main()
