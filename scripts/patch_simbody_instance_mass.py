"""Exact retained-source patch recipe; never modifies the source installation."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
RECIPE=ROOT/'scripts/simbody_instance_mass_source_identity.json'

def patch(overlay):
    expected=json.loads(RECIPE.read_text())['source_sha256']
    for name,digest in expected.items():
        if hashlib.sha256((overlay/name).read_bytes()).hexdigest()!=digest:raise ValueError('Unmatched pinned Simbody source: '+name)
    def edit(name,old,new,count=1):
        p=overlay/name;s=p.read_text()
        if s.count(old)!=count:raise ValueError('Patch anchor mismatch: '+name+' '+old)
        p.write_text(s.replace(old,new))
    edit('MobilizedBodyImpl.h','return getMyRigidBodyNode().getMassProperties_OB_B();','return iv.bodyMassProperties[getMyMobilizedBodyIndex()];')
    edit('RigidBodyNode.h','void calcJointIndependentKinematicsPos(\n    SBTreePositionCache& pc) const;','void calcJointIndependentKinematicsPos(\n    const MassProperties& mass, SBTreePositionCache& pc) const;')
    edit('RigidBodyNode.cpp','void RigidBodyNode::calcJointIndependentKinematicsPos(\n    SBTreePositionCache&    pc) const','void RigidBodyNode::calcJointIndependentKinematicsPos(\n    const MassProperties& mass, SBTreePositionCache& pc) const')
    edit('RigidBodyNode.cpp','getUnitInertia_OB_B().reexpress(~R_GB)','mass.getUnitInertia().reexpress(~R_GB)')
    edit('RigidBodyNode.cpp','R_GB*getCOM_B()','R_GB*mass.getMassCenter()')
    edit('RigidBodyNode.cpp','SpatialInertia(getMass(), p_BBc_G, G_Bo_G)','SpatialInertia(mass.getMass(), p_BBc_G, G_Bo_G)')
    edit('RigidBodyNode.cpp','const SpatialVec b  = getMass() *','const SpatialVec b  = getMk_G(pc).getMass() *')
    edit('RigidBodyNodeSpec.h','calcJointIndependentKinematicsPos(pc);','calcJointIndependentKinematicsPos(iv.bodyMassProperties[nodeNum],pc);')
    edit('RigidBodyNode_Weld.cpp','getUnitInertia_OB_B().reexpress(~R_GB)','iv.bodyMassProperties[nodeNum].getUnitInertia().reexpress(~R_GB)')
    edit('RigidBodyNode_Weld.cpp','R_GB*getCOM_B()','R_GB*iv.bodyMassProperties[nodeNum].getMassCenter()')
    edit('RigidBodyNode_Weld.cpp','SpatialInertia(getMass(), p_BBc_G, G_Bo_G)','SpatialInertia(iv.bodyMassProperties[nodeNum].getMass(), p_BBc_G, G_Bo_G)')
    edit('RigidBodyNode_LoneParticle.cpp','Mat33(1/getMass())','Mat33(1/sbs.getInstanceVars().bodyMassProperties[nodeNum].getMass())')
    edit('RigidBodyNode_LoneParticle.cpp','q + getCOM_B()','q + sbs.getInstanceVars().bodyMassProperties[nodeNum].getMassCenter()')
    edit('RigidBodyNode_LoneParticle.cpp','SpatialInertia(getMass(), getCOM_B(), getUnitInertia_OB_B())','SpatialInertia(sbs.getInstanceVars().bodyMassProperties[nodeNum].getMass(), sbs.getInstanceVars().bodyMassProperties[nodeNum].getMassCenter(), sbs.getInstanceVars().bodyMassProperties[nodeNum].getUnitInertia())')
    edit('RigidBodyNode_LoneParticle.cpp','const SBTreePositionCache&,\n        const SBArticulatedBodyInertiaCache&,\n        const SBArticulatedBodyVelocityCache&,','const SBTreePositionCache& pc,\n        const SBArticulatedBodyInertiaCache&,\n        const SBArticulatedBodyVelocityCache&,')
    edit('RigidBodyNode_LoneParticle.cpp','getMass()*udot','getMk_G(pc).getMass()*udot')
    edit('RigidBodyNode_LoneParticle.cpp','eps/getMass()','eps/getMk_G(pc).getMass()',2)
    # Append an exported additive ABI, leaving the installed public header alone.
    p=overlay/'MobilizedBody.cpp';p.write_text(p.read_text()+r'''
// IHM diagnostic State-instance mass ABI, separately built source variant.
extern "C" void ihm_simbody_set_instance_mass(SimTK::SimbodyMatterSubsystem* matter,
 SimTK::State* state,int index,const SimTK::MassProperties* proposed) {
 using namespace SimTK;
 if(!matter||!state||!proposed||index<=0||index>=matter->getNumBodies())
  throw std::invalid_argument("invalid instance mass target");
 const auto valid=[](const MassProperties& mp) {
  if(!mp.isFinite()||mp.getMass()<=0)return false;
  const auto central=mp.calcCentralInertia().asSymMat33();
  // Conservative first ABI: nondegenerate physical rigid-body inertia only.
  const auto covariance=Real(.5)*(central(0,0)+central(1,1)+central(2,2))*Mat33(1)-Mat33(central);
  return Inertia::isValidInertiaMatrix(central) && covariance(0,0)>0
   && covariance(0,0)*covariance(1,1)-covariance(0,1)*covariance(1,0)>0
   && det(covariance)>0;
 };
 const auto& body=matter->getMobilizedBody(MobilizedBodyIndex(index));
 if(!valid(body.getDefaultMassProperties())||!valid(*proposed))
  throw std::invalid_argument("instance mass requires nondegenerate physical rigid-body inertia");
 matter->updRep().updBodyMassProperties(*state,MobilizedBodyIndex(index))=*proposed;
 state->invalidateAllCacheAtOrAbove(Stage::Instance);
}
''')
    return sorted(expected)
